import os, time
from multiprocessing import Pool
from datetime import date
import argparse
import ast
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

#### Input
# The radius should be consistent with the radius in prior codes
radius = int(1400)

num_cpus=int(os.environ["SLURM_CPUS_ON_NODE"])
print("Running a pool with %s workers"%num_cpus)

today = date.today()
print("Today's date:", today)

# Swicth to the raw file directory
path_init = os.getcwd() + '/raw_files'
os.chdir( path_init )

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

# Create a grid of coordinates
x, y = np.arange(0, 3840), np.arange(0, 3840)
x_grid, y_grid = np.meshgrid(x, y)

# Calculate distances for all pixels at once
distances = np.sqrt((1919 - x_grid)**2 + (1919 - y_grid)**2)

# Create a mask for pixels outside the circle
mask = distances > radius

def process_file(filename):

    # Open the file and convert the data type
    hdu = fits.open(filename)
    data = (hdu[0].data).astype('float64') 
    
    if "rrhr" in filename:
        print("rrhr file: ", filename)
        data[data < 0] = np.nan

    # Replace the outer pixels with nans
    data[mask] = np.nan
    
    # Save the file
    hdu[0].data = data
    hdu.writeto(filename.removesuffix('.fits')+ '_nan.fits', overwrite=True)    
    print(f'File saved for {filename}')

    return data

if __name__=="__main__":
    with Pool(num_cpus) as p:
        # Retrieves each type of files in the directory
        fn = sorted(glob('*cnt.fits')) + sorted(glob('*rrhr.fits'))
        output = p.map(process_file, fn)


