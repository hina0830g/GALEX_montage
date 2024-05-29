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
from photutils.segmentation import make_2dgaussian_kernel

from reproject import reproject_interp

home_dir = "/xdisk/hamden/hina0830/venv39"

#### Input
# The radius should be consistent with the radius in prior codes
radius = int(1400)

num_cpus = int(os.environ["SLURM_CPUS_ON_NODE"])
print("Running a pool with %s workers" % num_cpus)

today = date.today()
print("Today's date:", today)

# Swicth to the raw file directory
path_init = home_dir + "/raw_files"
os.chdir(path_init)

# Create an ArgumentParser object
parser = argparse.ArgumentParser()

# Add an argument for the list of directory and filename pairs
# nargs='+' indicates that the argument can take one or more values
parser.add_argument(
    "--pair-list", nargs=2, type=str, help="Pair of integers (int1 int2)"
)

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
        path_new = str(os.getcwd()) + "/" + dir_name
        os.chdir(path_new)

        print("New location: ", path_new)
    except FileNotFoundError:
        print("Error")

#os.chdir("2024-04-21-RA195-DEC21")
print("New location: ", os.getcwd())

def make_mask(data):
    radius = 1400 # Default. Should not be changed
    
    # Create a grid of coordinates
    x, y = np.arange(0, int(len(data))), np.arange(0, int(len(data[0])))
    x_grid, y_grid = np.meshgrid(x, y)
    x_cent, y_cent = int(round(len(data) / 2)), int(round(len(data[0]) / 2))

    # Calculate distances for all pixels at once
    distances = np.sqrt((x_cent - x_grid) ** 2 + (x_cent - y_grid) ** 2)

    # Create a mask for pixels outside the circle
    mask = distances > radius
    
    return mask


def process_file(filename):
    print("Processing: ", filename)
    # Open the file and convert the data type
    hdu = fits.open(filename)
    data = (hdu[0].data).astype("float64")

    if "rrhr" in filename:
        print("rrhr file: ", filename)
        data[data < 0] = np.nan

    elif "int" in filename:
        print("int file: ", filename)
        data[data < 0] = np.nan
        # Apply Gaussian filter to intensity file
        kernel = make_2dgaussian_kernel(7.0, size=21)  # FWHM = 7
        data = convolve(data, kernel)
    
    # Replace the outer pixels with nans
    mask = make_mask(data)
    data[mask] = np.nan

    # Save the file
    hdu[0].data = data
    hdu.writeto(filename.removesuffix(".fits") + "_preprocessed.fits", overwrite=True)
    print(f"File saved for {filename}")

    return data

def reproject(file_tuple):
    fn_cnt, fn_flags = file_tuple
    print("Processing: ", fn_cnt, fn_flags)
    hdu1 = fits.open(fn_cnt)[0]
    hdu2 = fits.open(fn_flags)[0] # flag, 480 by 480
    
    array, footprint = reproject_interp(hdu2, hdu1.header, order = "nearest-neighbor")  # reproject flag file to cnt (3840 by 3840) # order = " reproject_interp "
    print(np.shape(array))
    
    print(fn_flags.removesuffix('.fits') + '_wcs.fits', " reprojected and saved. ")
    fits.writeto(fn_flags.removesuffix('.fits') + '_wcs.fits', array, hdu1.header, overwrite=True)
    
    return array


if __name__ == "__main__":
    with Pool(num_cpus) as p:
        # Retrieves each type of files in the directory
        fn = sorted(glob("*cnt.fits")) + sorted(glob("*rrhr.fits")) + sorted(glob("*int.fits"))
        print("nans filename: ", fn)
        # Nan the outside
        output = p.map(process_file, fn)
        
        fn_reproject =  sorted(glob("*cnt.fits")) + sorted(glob("*flags.fits"))
        
        mapped_argument_reproject = list(zip(sorted(glob("*cnt.fits")), sorted(glob("*flags.fits"))))
        
        print("mapped_argument_reproject: ", mapped_argument_reproject)
        # Reproject flag files
        output_reproject = p.map(reproject, mapped_argument_reproject)
        
        
