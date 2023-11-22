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

radius = int(1400)

# Create a grid of coordinates
x, y = np.arange(0, len(data)), np.arange(0, len(data))
x_grid, y_grid = np.meshgrid(x, y)
x_cent, y_cent = int(round(len(data) / 2)), int(round(len(data) / 2))

# Calculate distances for all pixels at once
distances = np.sqrt((x_cent- x_grid)**2 + (y_cent- y_grid)**2)

# Create a mask for pixels outside the circle
mask = distances > radius

# Iterates through all the files ending with _cnt.fits
#for file in glob.glob("*cnt.fits"):

def process_file(filename):
    # Opens cnt file and the data
    hdu = fits.open(filename)
    data = (hdu[0].data).astype('float64')

    if "rrhr" in filename:
        print("rrhr file: ", filename)
        # Replace pixels below 0 with nan
        data[data < 0] = np.nan
    
    # Replace pixels outside the radius with nans
    data[mask] = np.nan

    # Display the final image; outside of the radius should look white
    fig, ax = plt.subplots()
    plt.imshow(cnt_data, cmap='rainbow', norm=colors.LogNorm(vmin=1, vmax=50))
    plt.colorbar()
 
    # Save the file    
    hdu[0].data = data
    hdu.writeto(filename.removesuffix('.fits')+ '_nan.fits', overwrite=True)
    print(f'File saved for {filename}')
    return data

if __name__ == "__main__":
    fn = sorted(glob(*cnt.fits)) + sorted(glob(*rrhr.fits))
    output = process_file(fn) 