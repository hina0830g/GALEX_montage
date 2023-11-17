# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 16:41:06 2022

@author: Hina
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 23:50:15 2022

@author: Hina
"""
import os
import sys
import shutil
from IPython.display import Image
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import glob
# Description: This program finds cnt and rrhr files in your directory
# and replaces pixes outside of the circle with nans 
# It splits an image into 4 quadrants and applies the distance formula
# Recommended radius is 1450 or smaller

# Specify the path to start the program
path_init = 'C:/Users/Amiti/Documents/Spring2022/GALEX/raw_files/Test_py/RA195'
os.chdir( path_init )
os.getcwd()

radius = int(input('Give a radius; anything outside of the radius will be replaced with nan. Recommended radius is 1450 or smaller: '))

# Defines the distance formula
def dist(x, y):
    d = np.sqrt((1919-x)**2 + (1919-y)**2)
    return d
 
# Determines how many cnt/rrhr files are in the directory
i = 0
for file in glob.glob("*cnt.fits"):
    i += 1
    
print('\n')
print(i , 'cnt/rrhr files available')

# Iterates through all the files ending with _cnt.fits
for file in glob.glob("*cnt.fits"):
    # Opens cnt file and the data
    cnt_fn = file.removesuffix('.fits')
    hdu = fits.open(file)
    cnt_data = (hdu[0].data).astype('float64')
    # Below iterates through each row and column and replaces pixels with nans
    # Top left corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            cnt_data[i][j] = np.where(dist(i,j) > radius, np.nan, cnt_data[i][j])
    # bottom left corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            cnt_data[1919+i][j] = np.where(dist(1919+i,j) > radius, np.nan, cnt_data[1919+i][j])
    # Top right corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            cnt_data[i][1919+j] = np.where(dist(i,1919+j) > radius, np.nan, cnt_data[i][1919+j])
    # bottom right corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            cnt_data[1919+i][1919+j] = np.where(dist(1919+i,1919+j) > radius, np.nan, cnt_data[1919+i][1919+j])
  
    # Display the final image; outside of the radius should look white
    fig, ax = plt.subplots()
    plt.imshow(cnt_data, cmap='rainbow', norm=colors.LogNorm(vmin=1, vmax=50))
    plt.colorbar()
    # Saves the new data as a new file called XXX_cnt_nan.fits
    hdu[0].data = cnt_data
    print(cnt_fn + '_nan.fits')
    hdu.writeto(cnt_fn + '_nan.fits') 

# Iterates through all the files ending with _rrhr.fits
for file in glob.glob("*rrhr.fits"):
    # Opens cnt file and the data
    rrhr_fn = file.removesuffix('.fits')
    hdu = fits.open(file)
    rrhr_data = (hdu[0].data).astype('float64')
    # Below iterates through each row and column and replaces pixels with nans
    # Top left corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            rrhr_data[i][j] = np.where(dist(i,j) > radius, np.nan, rrhr_data[i][j])
    # bottom left corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            rrhr_data[1919+i][j] = np.where(dist(1919+i,j) > radius, np.nan, rrhr_data[1919+i][j])
    # Top right corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            rrhr_data[i][1919+j] = np.where(dist(i,1919+j) > radius, np.nan, rrhr_data[i][1919+j])
    # bottom right corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            rrhr_data[1919+i][1919+j] = np.where(dist(1919+i,1919+j) > radius, np.nan, rrhr_data[1919+i][1919+j])
  
    # Display the final image; outside of the radius should look white
    fig, ax = plt.subplots()
    plt.imshow(rrhr_data, cmap='rainbow', norm=colors.LogNorm(vmin=1, vmax=50))
    plt.colorbar()
    
    # Saves the new data as a new file called XXX_rrhr_nan.fits
    hdu[0].data = rrhr_data
    print(rrhr_fn + '_nan.fits')
    hdu.writeto(rrhr_fn + '_nan.fits') 
