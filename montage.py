# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 22:33:29 2022

@author: Hina Goto
"""
import numpy as np
import matplotlib.pyplot as py
import MontagePy as Montage
from MontagePy.main import *
from MontagePy.archive import *
from IPython.display import Image
from astroquery.mast import Catalogs
from astroquery.mast import Observations
import pandas as pd
from astropy.table import QTable
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve, convolve_fft
from astropy.utils.data import get_pkg_data_filename
from scipy.signal import convolve as scipy_convolve
import os
import shutil
import matplotlib.pyplot as plt
from matplotlib import colors

### This code assumes that all the files are in "raw" directory in order for this to run

path_init = 'C:/Users/Amiti/Documents/Spring2022/GALEX'
os.chdir( path_init )
os.getcwd()
home = os.getcwd()
# Enter the same coordinate used in query below
c = SkyCoord(ra=195*u.degree, dec=21*u.degree, frame='icrs') # Copy and paste from the query

# These are the parameters defining the mosaic we want to make
location = c.to_string('hmsdms')
size     = 5.0
dataset  = "2MASS J"
workdir  = 'raw_files/Test_py/RA195/nan_fits' # where raw directory is located

# We create and move into a subdirectory in this notebook
# but we want to come back to the original startup directory
# whenever we restart the processing
try:
    home
except:
    home = os.getcwd()

os.chdir(home)
print("Start up folder: " + home)

# Clean out any old copy of the work tree, then remake it 
# and the set of the subdirectories we will need.

print("Work directory: " + workdir, flush=True)

os.chdir(workdir)

#os.makedirs("raw")
os.makedirs("projected")
os.makedirs("diffs")
os.makedirs("corrected")

# Create the FITS header for the mosaic
rtn = mHdr(location,size,size, "region.hdr")
print("mHdr: " + str(rtn), flush=True)

# Retreive archive images covering the region from the raw 
# directory and then scan the images for their coverage metadata. 
rtn = mImgtbl("raw", "rimages.tbl")
print("mImgtbl (raw): " + str(rtn), flush=True)

# Reproject the original images to the frame of the 
# output FITS header we created
rtn = mProjExec("raw", "rimages.tbl", "region.hdr", projdir="projected")
print("mProjExec: " + str(rtn), flush=True)

mImgtbl("projected", "pimages.tbl")
print("mImgtbl (projected): " + str(rtn), flush=True)

# Coadd the projected images without background correction.
# This step is just to illustrate the need for background correction
# and can be omitted
rtn = mAdd("projected", "pimages.tbl", "region.hdr", "uncorrected.fits", coadd = 3)
print("mAdd: " + str(rtn), flush=True)

# Make a PNG rendering of the data and display it.
rtn = mViewer("-ct 1 -gray uncorrected.fits -2s max gaussian-log -out uncorrected.png", "", mode=2)
print("mViewer: " + str(rtn), flush=True)
Image(filename='uncorrected.png')