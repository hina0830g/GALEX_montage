# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 00:20:28 2022

@author: Hina Goto
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import colors
from datetime import date
import argparse
import ast
import os
import sys

# Description: This code combines the two mosaic images 
# and creates one mosaic image.
# The new image will be mosaic(cnt)/mosaic(rrhr); it shows
# the combined image and saves the image if wanted

today = date.today()
print("Today's date:", today)

# Swicth to the raw file directory
path_init = os.getcwd() + '/raw_files'
os.chdir( path_init )
print("Path changed to: ", os.getcwd())

# Create an ArgumentParser object
parser = argparse.ArgumentParser()

# Add an argument for the list of directory and filename pairs
# nargs='+' indicates that the argument can take one or more values
parser.add_argument('--pair-list', nargs=2, type=str, help='Pair of integers (int1 int2)')

# Parse the command-line arguments
args = parser.parse_args()

# Use the list of directory and filename pairs in your script
if args.pair_list:
    # Split a pair into a directory and a filename
    RA, DEC = map(int, args.pair_list)
    print("RA:", RA)
    print("DEC:", DEC)
    try:
        # Get the path from the coordinate information
        dir_name = str(today) + "-RA" + str(RA) + "-DEC" + str(DEC)
        
        # Switch to the new directory 
        path_new = str(os.getcwd()) + '/' + dir_name
        os.chdir( path_new )
        
        print("New location: ", path_new)
    except FileNotFoundError:
        print("Error")

# Shows the mosaic of cnt 
cnt_fits = 'cnt/uncorrected.fits'
hdu = fits.open(cnt_fits)
cnt_data = (hdu[0].data).astype('float64')

# Shows the mosaic of rrhr
rrhr_fits = 'rrhr/uncorrected.fits'
hdu = fits.open(rrhr_fits)
rrhr_data = (hdu[0].data).astype('float64')

# Shows the mosaic of the two combined (cnt divided by rrhr)
divided = cnt_data/rrhr_data

file_name = dir_name + '_combined.fits'
hdu[0].data = divided
hdu.writeto(file_name, overwrite=True)

