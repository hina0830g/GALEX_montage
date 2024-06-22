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
area_circle = np.pi * (radius ** 2)
artifacts_cutoff = 20
edge_thickness = 500

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

#os.chdir("2024-05-16-RA10-DEC41")
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
    # drop the negative pixels in rrhr and int.fits bty replacing them with nans
    # and apply a Gaussian smoothing filter to inte
    
    print("Process_file at: ", filename)
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
    
    # Trim 500 blank pixels from each side
    data = data[edge_thickness:len(data)-edge_thickness, edge_thickness:len(data[0])-edge_thickness]

    # Save the file
    hdu[0].data = data
    hdu.writeto(filename.removesuffix(".fits") + "_preprocessed.fits", overwrite=True)
    print(f"File saved for {filename}")

    return data

def preprocess(file_tuple):
    print("Processing tuple ", file_tuple)
    cnt_file, rrhr_file, int_file, skybg_file, flags_file = file_tuple

    cnt_data, rrhr_data, int_data, skybg_data, flags_data = (
        fits.getdata(cnt_file).astype("float64"),
        fits.getdata(rrhr_file).astype("float64"),
        fits.getdata(int_file).astype("float64"),
        fits.getdata(skybg_file).astype("float64"),
        fits.getdata(flags_file).astype("float64"))
    
    hdu1 = fits.open(cnt_file)[0]
    hdu2 = fits.open(flags_file)[0] # flag, 480 by 480
    
    flags_data, footprint = reproject_interp(hdu2, hdu1.header, order = "nearest-neighbor") 
    flag_largeval = np.where(flags_data < 128, 0 , flags_data)
    nonzero_pix = np.count_nonzero(flag_largeval)
    artifacts_frac = (nonzero_pix / area_circle) * 100
    
    print("frac: " + (f"{artifacts_frac:.2f}") + "%")
    
    # Preprocess the files only if fraction of the flagged pixels is above our cutoff value 
    if artifacts_frac < artifacts_cutoff:
        prerocessed_cnt = process_file(cnt_file)
        prerocessed_rrhr = process_file(rrhr_file)
        prerocessed_int = process_file(int_file)
        
        flags_data = flags_data[edge_thickness:len(flags_data)-edge_thickness, edge_thickness:len(flags_data[0])-edge_thickness]
        
        # Save reprojected flags
        fits.writeto(flags_file.removesuffix('.fits') + '_wcs.fits', flags_data, hdu1.header, overwrite=True)
        print(flags_file.removesuffix('.fits') + '_wcs.fits', " reprojected and saved. ")
        
        hdu = fits.open(skybg_file)
        
        skybg_data = skybg_data[edge_thickness:len(skybg_data)-edge_thickness, edge_thickness:len(skybg_data[0])-edge_thickness]
        hdu[0].data = skybg_data
        hdu.writeto(skybg_file, overwrite=True)
        
    else:
        print(cnt_file.removesuffix('.fits') + " discarded at " + (f"{artifacts_frac:.2f}") + "%")
        os.remove(skybg_file)

if __name__ == "__main__":
    with Pool(num_cpus) as p:
        fn_cnt, fn_rrhr, fn_int, fn_skybg, fn_flags = (
            sorted(glob("*-cnt.fits")),
            sorted(glob("*rrhr.fits")),
            sorted(glob("*int.fits")),
            sorted(glob("*skybg.fits")),
            sorted(glob("*flags.fits")),
        )
        
        mapped_new = list(zip(fn_cnt, fn_rrhr, fn_int, fn_skybg, fn_flags)) 
        if len(fn_cnt) == len(fn_rrhr) == len(fn_int) == len(fn_skybg) == len(fn_flags):
        
            print("mapped_new: ", mapped_new)
            print(len(mapped_new), "files total.")

            output = p.map(preprocess, mapped_new)

        
