# These are routines that were stripped from other python files.
#
# Known contributors include:
#
# Charles Chen
# Matthew Holman
# Michael Lackner
# Ed Lin
# Matthew Payne
# Pavlos Protopapas
#

import sys
import warnings
import math
import random
from collections import defaultdict
from collections import Counter

import scipy as sp
from astropy.time import Time
#import ephem
import numpy as np
from scipy import spatial

from astropy.io import fits


#def sunradec(jd):
#        day = jd -2415020
#        sun = ephem.Sun(day)
#        return math.degrees(sun.a_ra), math.degrees(sun.a_dec)

def compute_snr(mag_sig):
    small_sig = 0.0001
    nsr = 1.0 - np.power(10.0, -0.4*(np.abs(mag_sig) + small_sig))
    snr = 1.0/nsr
    return snr

# =============================
# JD Convertor (use pyephem)
# This piece came from one of Ed Lin's codes.
# =============================a
def date2JD(Cinput):
    #2010-07-22T10:17:00.187239
    '''
    Cinput=str(Cinput)
    Cdate = Cinput.split('T')[0]
    Ctime = Cinput.split('T')[1]
    Cyear = Cdate.split('-')[0]
    Cmonth= Cdate.split('-')[1]
    Cday  = Cdate.split('-')[2]
    d = ephem.Date('%s/%s/%s %s' %(Cyear,Cmonth,Cday,Ctime))
    print (float(d)+2415020)
    '''
    # Alternative approach by Payne that uses astropy instead
    t = Time(Cinput, format='isot')
    return t.jd

 #  Below this point is mostly stuff added by MJH.


def get_chip_names(xs, ys, xw=5045., x0=-20180, yw=5205., y0=-20820.):
        def corners(x, y):
                xc = 7-int(math.floor((x-x0)/xw))
                yc = 7-int(math.floor((y-y0)/yw))
                if xc<0 or xc>7 or yc<0 or yc>7 or (xc==0 and yc==0) or (xc==7 and yc==7) or (xc==0 and yc==7) or (xc==7 and yc==0):
                        return 'Corner'
                else:
                        return 'XY' + '%d%d' %(xc, yc)
        xys = zip(xs, ys)
        names = [(i, corners(x, y)) for i, (x, y) in enumerate(xys)]
        return names

def inside_chip(x, y, xw=4846, yw=4868):
 return np.logical_and(np.logical_and(x>0.0, x<xw), np.logical_and(y>0.0, y<yw))

def inside_cell_interior(x, y, xmin=13., xmax=580., ymin=6., ymax=587):
        return foldx(x)>xmin and foldx(x)<xmax and foldy(y)>ymin and foldy(y)<ymax

# Note the direction is flipped.  There
# might be a better way of dealing with this.
def foldx(x, xw=590, xgap=18):
    #xp = xw+xgap - np.mod(x, (xw+xgap))
    xp = np.mod(x, (xw+xgap))
    return xp

def foldy(y, yw=598, ygap=12):
    yp = np.mod(y, (yw+ygap))
    return yp

# There should be a way to do this in numpy
# Can include a trimmed region here too.
def inside_cell(x, y, xw=590, xgap=18, yw=598, ygap=12):
    return np.logical_and(foldx(x, xw=xw, xgap=xgap) < xw, foldy(y, yw=yw, ygap=ygap) < yw)

# There should be a way to do this in numpy
# Can include a trimmed region here too.
def inside_cell_x(x, xw=590, xgap=18):
    return foldx(x, xw=xw, xgap=xgap) < xw

# There should be a way to do this in numpy
# Can include a trimmed region here too.
def inside_cell_y(y, yw=598, ygap=12):
    return foldy(y, yw=yw, ygap=ygap) < yw


# There should be a way to do this in numpy
def get_cell_name(x, y, xgap=18, ygap=12):
    nx = (x/(590+xgap)).astype(int)
    ny = (y/(598+ygap)).astype(int)
    return [('xy' + '%d%d' %(nx1, ny1)) for nx1, ny1 in zip(nx, ny)]

def get_cell_number(x, y, xgap=18, ygap=12):
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
        
    nx = (x/(590+xgap)).astype(int)
    ny = (y/(598+ygap)).astype(int)
    return nx, ny

def mask_good(chip, cell, xpix, ypix, mask, xw=590, yw=598):
        cell = np.atleast_1d(cell)
        xpix = np.atleast_1d(xpix)
        ypix = np.atleast_1d(ypix)                
        pix_ok = np.logical_and(np.logical_and(xw-xpix-1>=0, xw-xpix-1<xw),
                                np.logical_and(ypix>=0, ypix<yw))
        # If the pixel is not ok but the cell number is ok, then the detection is in
        # intra-chip gaps.
        cell_nums = [(int(cn[2:3]), int(cn[3:4])) for cn in cell]
        # If the cell is not ok, then the detection is off the chip.
        cell_ok =[(nx>=0 and nx<8 and ny>=0 and ny<8) for nx, ny in cell_nums]
        return np.array([pok and cok and\
                (mask[chip][cn].data[ypx][xw-xpx-1] ==0) for pok, cok, xpx, ypx, cn in zip(pix_ok, cell_ok, xpix, ypix, cell)])

def mask_value(chip, cell, xpix, ypix, mask, xw=590):
        cell = np.atleast_1d(cell)
        xpix = np.atleast_1d(xpix)
        ypix = np.atleast_1d(ypix)                
        return mask[chip][cell].data[ypix][xw-xpix-1]

def cell_good(chip, cell, cell_summary_dict):
        cell = np.atleast_1d(cell)
        cell_nums = [(int(cn[2:3]), int(cn[3:4])) for cn in cell]        
        cell_ok =[(nx>=0 and nx<8 and ny>=0 and ny<8) for nx, ny in cell_nums]        
        return np.array([cok and (cell_summary_dict[chip][cll]=='S') for cok, cll in zip(cell_ok, cell)])

## This function grabs the summary about cells from the
## headers of each chip in the set of HDUs.
def get_cell_summaries(hdus):
    chip_cell_dict=defaultdict(dict)
    for hdu in hdus:
        if '.hdr' in hdu.name:
            chip = hdu.name.rstrip('.hdr')
            cell_strings = [c for c in list(hdu.header['COMMENT']) if 'xy' in c]
            cell_dict={}
            for line in cell_strings:
                key_string = line[0:7]
                value_string = line[7:25]
                y_string = key_string[5]
                for i, res in enumerate(value_string.strip().split()):
                    cell_name = 'xy%s%s' % (i, y_string)
                    cell_dict[cell_name] = res
            chip_cell_dict[chip] = cell_dict
    return chip_cell_dict
    
def extract_header_info(filename):
    info_dict={}
    with open(filename) as file:
        for fn in file:
            fn = fn.rstrip()
            with fits.open(fn) as hdus:
                phdr = hdus['PRIMARY'].header
                xhdr = hdus['XY33.hdr'].header
                psf_hdr = hdus['XY33.psf'].header

                try:
                    shutter_open = date2JD(hdus['XY11.hdr'].header['SHUTOUTC'].replace(' ', 'T'))
                    exptime = phdr['exptime']
                    jdObs = shutter_open + 0.5*exptime/(24*3600.)
                    filt = phdr['FILTERID']
                    exttype = psf_hdr['EXTTYPE']

                    ra = phdr['RA']
                    dec = phdr['Dec']
                    
                    info_dict[fn] = (fn, exttype, jdObs, float(ra), float(dec), filt, exptime)
                    
                except:
                    info_dict[fn]=None
                    pass
                
    return info_dict


def ra_dec2pix(ra, dec, nside=32, nested=True):
    phi   = ra * np.pi / 180.
    theta = np.pi / 2.0 - dec * np.pi / 180.
    pix = hp.ang2pix(nside, theta, phi, nested)
    return pix

def ra_dec2vec(ra, dec):
    radeg=180./np.pi
    x = np.cos(ra/radeg)*np.cos(dec/radeg)
    y = np.sin(ra/radeg)*np.cos(dec/radeg)
    z = np.sin(dec/radeg)
    return np.array((x,y,z)).T

def ra_dec2tuple(ra, dec):
    radeg=180./np.pi
    x = np.cos(ra/radeg)*np.cos(dec/radeg)
    y = np.sin(ra/radeg)*np.cos(dec/radeg)
    z = np.sin(dec/radeg)
    return (x,y,z)
