# -*- coding: utf-8 -*-
"""
Created on Fri May 27 21:47:21 2022

@author: Hina Goto
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
from glob import glob
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve, convolve_fft
from astropy.utils.data import get_pkg_data_filename

# Specify the initial path that contains edited cnt and rrhr files
# as well as raw skybg and objmask files 
path_init = 'C:/Users/Amiti/Documents/Spring2022/GALEX/raw_files/Test_py/test_gz/test_mask'
os.chdir( path_init )
os.getcwd()

# Retrieves each type of files in the directory
fn_cnt, fn_rrhr, fn_objmask, fn_skybg = glob('*cnt_nan.fits'), glob('*rrhr_nan.fits'), glob('*objmask.fits'),glob('*skybg.fits')

# The radius should be consistent with the radius in prior codes
sigma = float(input('Enter the value of sigma for the Gaussian filter; 5~10 is recommended '))
radius = int(input('Give a radius; anything outside of the radius will be replaced with nan. Recommended radius is 1450 or smaller: '))

# Defines the distance formula
def dist(x, y):
    d = np.sqrt((1919-x)**2 + (1919-y)**2)
    return d

# Below iterates through sets of each type of files to apply masking;
# Each type of files should be corresponding to each other (same file name before the suffix)
for index in range(len(fn_cnt)): 
    # cnt file
    hdu = fits.open(fn_cnt[index])
    cnt_data = (hdu[0].data).astype('float64') # REtr
  
    # Shows cnt files before masking and filtering 
    fig, ax = plt.subplots()
    plt.imshow(cnt_data, cmap='rainbow', norm=colors.LogNorm(vmin=.1, vmax=1300))
    plt.colorbar()
  
    # rrhr file
    hdu = fits.open(fn_rrhr[index])
    rrhr_data = (hdu[0].data).astype('float64')
    
    #objmask file
    hdu = fits.open(fn_objmask[index])
    objmask_data = (hdu[0].data).astype('float64')
  
    # skybg file
    hdu = fits.open(fn_skybg[index])
    skybg_data = (hdu[0].data).astype('float64')

    # Now apply the masking by quadrant 
    # Top left corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            if np.isnan(cnt_data[i][j]) == False: 
                cnt_data[i][j] = np.where(objmask_data[i][j] !=0, skybg_data[i][j]*rrhr_data[i][j], cnt_data[i][j])
              
    # bottom left corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            if np.isnan(cnt_data[1919+i][j]) == False: 
                cnt_data[1919+i][j] = np.where(objmask_data[1919+i][j] !=0, skybg_data[1919+i][j]*rrhr_data[1919+i][j], cnt_data[1919+i][j])
              
    # top right corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            if np.isnan(cnt_data[i][1919+j]) == False: 
                cnt_data[i][1919+j] = np.where(objmask_data[i][1919+j] !=0, skybg_data[i][1919+j]*rrhr_data[i][1919+j], cnt_data[i][1919+j])
              
    # bottom right corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            if np.isnan(cnt_data[1919+i][1919+j]) == False: 
                cnt_data[1919+i][1919+j] = np.where(objmask_data[1919+i][1919+j] !=0, skybg_data[1919+i][1919+j]*rrhr_data[1919+i][1919+j], cnt_data[1919+i][1919+j])

    fig, ax = plt.subplots()
    plt.imshow(cnt_data, cmap='rainbow', norm=colors.LogNorm(vmin=.1, vmax=1300))
    plt.colorbar()
    
    # Applies gaussian filter 
    gauss_kernel = Gaussian2DKernel(sigma)
    smoothed_data_gauss = convolve(cnt_data, gauss_kernel)
    
    # Replaces outer pixels with nans again to clean it up 
    # Top left corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            smoothed_data_gauss[i][j] = np.where(dist(i,j) > radius, np.nan, smoothed_data_gauss[i][j])
    # bottom right corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            smoothed_data_gauss[1919+i][j] = np.where(dist(1919+i,j) > radius, np.nan, smoothed_data_gauss[1919+i][j])
    # Top right corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            smoothed_data_gauss[i][1919+j] = np.where(dist(i,1919+j) > radius, np.nan, smoothed_data_gauss[i][1919+j])
    # bottom left corner
    for i in range(0, 1919):
        for j in range(0, 1919):
            smoothed_data_gauss[1919+i][1919+j] = np.where(dist(1919+i,1919+j) > radius, np.nan, smoothed_data_gauss[1919+i][1919+j])
    
    # Displays the filtered image
    fig, ax = plt.subplots()
    plt.imshow(smoothed_data_gauss, cmap='rainbow', norm=colors.LogNorm(vmin=1, vmax=50))
    plt.colorbar()
    
    # Saves the new data (masked and filtered) as a new file called 
    # XXX_cnt_mask_filtered.fits; we run Montage on this
    hdu[0].data = smoothed_data_gauss
    hdu.writeto(fn_cnt[index].removesuffix('.fits')+'_mask_filtered.fits')
 
os.chdir( path_init )
os.getcwd()

# Below creates new directories so we would have cnt/raw and rrhr/raw
# and necessary files will be placed in each "raw" directory;
# we will run MontagePy on files in raw directories.

# Create new directories    
os.mkdir('cnt')
os.mkdir('rrhr')

os.chdir(path_init + '/cnt')
os.mkdir('raw') # cnt/raw

os.chdir(path_init + '/rrhr')
os.mkdir('raw') # rrhr/raw

os.chdir( path_init )

# Transfers cnt and rrhr files into raw directories

# masked and filtered cnt files will be placed in /cnt/raw
for file in glob('*mask_filtered.fits'): # Filtered cnt files 
    src_path = path_init + '/' + file  
    dst_path = path_init + '/cnt/raw/' + file 
    shutil.move(src_path, dst_path)
    
os.chdir( path_init ) 
# rrhr files whose outer pixels have been replaced with nans will
# be placed in rrhr/raw   
for file in glob('*rrhr_nan.fits'): # Edited rrhr files 
    src_path = path_init + '/' + file  
    dst_path = path_init + '/rrhr/raw/' + file 
    shutil.move(src_path, dst_path)
