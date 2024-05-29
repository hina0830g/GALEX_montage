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
from photutils.datasets import make_100gaussians_image
from photutils.background import Background2D, MedianBackground
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import simple_norm
from photutils.detection import DAOStarFinder
from photutils.segmentation import SourceFinder
from photutils.segmentation import SourceCatalog
from astropy.stats import sigma_clipped_stats
from photutils.segmentation import detect_sources, make_2dgaussian_kernel
from photutils.aperture import CircularAperture
import pickle
from astropy.table import QTable

from photutils.aperture import EllipticalAperture
from matplotlib import colors
from astropy.coordinates import Angle

from matplotlib.patches import Ellipse
import scipy
from scipy import signal


sys.path.append("/xdisk/hamden/hina0830/venv39/")
# from cleaning import combined
import poisson_segment as ps
import cleaning as cl

print("Running object_detection.py")

### Adjust the parameters here ###
radius = 1400
sigma = 3.0
edge_thickness = 500

home_dir = "/xdisk/hamden/hina0830/venv39"

num_cpus = int(os.environ["SLURM_CPUS_ON_NODE"])
print("Running a pool with %s workers" % num_cpus)

os.chdir(home_dir + "/raw_files")

today = date.today()
print("Today's date:", today)

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
        path_init = os.getcwd() + "/" + dir_name
        os.chdir(path_init)
        print("Path changed to: ", os.getcwd())

    except FileNotFoundError:
        print("Error")

# os.chdir("2024-05-14-RA180-DEC12")
print("New location: ", os.getcwd())


def segmtantion(file_tuple):
    # Unpack the tuple and get lists of file names
    cnt_file, rrhr_file, skybg_file, int_file, flags_file = file_tuple
    print("segmentation running at: ", cnt_file)

    # Read all the files and get data
    cnt_data, rrhr_data, skybg_data, int_data, flags_data = (
        fits.getdata(cnt_file).astype("float64"),
        fits.getdata(rrhr_file).astype("float64"),
        fits.getdata(skybg_file).astype("float64"),
        fits.getdata(int_file).astype("float64"),
        fits.getdata(flags_file).astype("float64"),
    )

    print("processing: ", file_tuple)

    # Read cnt file here to get hdu
    hdu = fits.open(cnt_file)
    hdr = hdu[0].header

    fn_current = cnt_file.removesuffix("-cnt_nan.fits")  # fn as filename

    exptime, project_name = hdr["EXPTIME"], hdr["MPSTYPE"]
    print(fn_current, exptime, project_name)

    ### The default parameters ###
    edge_thickness = 500  # Trim 500 (blank) pixels on each side
    radius = 1400
    sigma = 3.0
    npixels_value = 500

    # If statements to determine the coefficient based on
    # exposure time
    if exptime >= 30000:  # DIS
        th_coeff = 15
    elif exptime >= 1500:  # MIS & NGS
        th_coeff = 10
    elif exptime >= 100:  # "AIS
        th_coeff = 6
    else:
        print("exception; exposure time is: ", exptime)
        th_coeff = 6

    # int_data[flags_data  >= 256] = np.nans

    mask = ps.segmentation(fn_current, int_data, npixels_value, th_coeff)  # Perform segmentations on large sources
    print("segmentation finished for ", cnt_file)

    cnt_noise = ps.poisson_noise(cnt_data, rrhr_data, mask, skybg_data)  # Perform Poisson infill
    print("poisson infill finished for ", cnt_file)

    # Trim each side (empty pixels) of the array by the edge_thickness (=500pix)
    cnt_infilled_trimmed = cnt_noise[edge_thickness:len(cnt_noise)-edge_thickness, edge_thickness:len(cnt_noise[0])-edge_thickness]
    rrhr_trimmed = rrhr_data[edge_thickness:len(rrhr_data)-edge_thickness, edge_thickness:len(rrhr_data[0])-edge_thickness]
    flags_trimmed = flags_data[edge_thickness:len(flags_data)-edge_thickness, edge_thickness:len(flags_data[0])-edge_thickness]
    
    # Apply Gaussian filter and clean the edges (cleaning)
    final_cnt = cl.combined(sigma=sigma, radius=radius, data=cnt_infilled_trimmed)
    print("cl.combined finished for ", cnt_file)

    # Flag the infilled image
    final_cnt[flags_trimmed >= 256] = np.nan

    # Save the new files
    hdu[0].data = final_cnt
    hdu.writeto(cnt_file.removesuffix("_preprocessed.fits") + "_Pinfilled_trimmed.fits", overwrite=True)
    print(f"new file saved for {cnt_file}")

    hdu[0].data = rrhr_trimmed
    hdu.writeto(rrhr_file.removesuffix(".fits") + "_trimmed.fits", overwrite=True)
    return cnt_data


def starfinder_new(file_tuple):
    cnt_file, int_file, flags_file, rrhr_file = file_tuple
    print("star finder running at: ", cnt_file)

    edge_thickness = 500  # Trim 500 (blank) pixels on each side
    mask_size = 16  # size of mask for each star in pixels

    th_coeff = 25.0  # for AIS
    DAO_fwhm = 10.385

    hdu = fits.open(cnt_file)

    cnt_data, rrhr_data, flags_data = (
        fits.getdata(cnt_file).astype("float64"),
        fits.getdata(rrhr_file).astype("float64"),
        fits.getdata(flags_file).astype("float64"))

    flags_data = flags_data[edge_thickness : len(flags_data) - edge_thickness, edge_thickness : len(flags_data[0]) - edge_thickness]
    
    # Manually create intensity file (infilled cnt/rrhr)
    divided_data = cnt_data / rrhr_data

    # Flag the bad pixels (non-zeros for now, but the threshold value tbd; 128 or 256)
    divided_data[flags_data != 0] = np.nan
    divided_data[np.isnan(divided_data)] = 0  # nan = 0

    # Run the point source finder and extract masked image, mask, and a list of coordinates
    masked_data, mask_data, coord_lis = ps.psfinder(divided_data, flags_data, th_coeff, DAO_fwhm, mask_size)
    print("point source finder completed for ", cnt_file)

    # --- Save files here ---
    
    cnt_file_int = cnt_file.replace('cnt', 'int')
    
    with open( cnt_file_int.removesuffix('Pinfilled_trimmed.fits') + str(len(coord_lis)) + '_coord.pkl', "wb") as f:
        pickle.dump(coord_lis,f)
        
    # Save masked data (out_image)
    hdu[0].data = divided_data
    hdu.writeto( cnt_file_int, overwrite=True)
    print( cnt_file.replace('cnt', 'int'), " saved.")

    # Save masked data (out_image)
    hdu[0].data = masked_data
    hdu.writeto( cnt_file_int.removesuffix('.fits') + '_' + str(len(coord_lis)) + '_masked.fits', overwrite=True)
    print(cnt_file_int.removesuffix('.fits') + '_' + str(len(coord_lis)) + '_masked.fits', " saved.")

    # Save mask file (bimage)
    hdu[0].data = mask_data
    hdu.writeto( cnt_file_int.removesuffix('Pinfilled_trimmed.fits') + str(len(coord_lis)) + '_mask.fits', overwrite=True)
    print(cnt_file_int.removesuffix('Pinfilled_trimmed.fits') + str(len(coord_lis)) + '_mask.fits', " saved.")
    
    return mask_data


if __name__ == "__main__":
    print("__name__ == __main__; currently at: ", os.getcwd())
    with Pool(num_cpus) as p:
        # Retrieves each type of files in the directory
        fn_cnt, fn_rrhr, fn_skybg, fn_int, fn_flags = (
            sorted(glob("*cnt_preprocessed.fits")),
            sorted(glob("*rrhr_preprocessed.fits")),
            sorted(glob("*skybg.fits")),
            sorted(glob("*int_preprocessed.fits")),
            sorted(glob("*flags_wcs.fits")),
        )
        mapped_argument = list(zip(fn_cnt, fn_rrhr, fn_skybg, fn_int, fn_flags))
        print("mapped_argument for for process_file : ", mapped_argument)

        ### --- run the segmentation code here --- ###
        if len(fn_cnt) == len(fn_rrhr)  == len(fn_skybg) == len(fn_skybg) == len(fn_flags):
            print(len(fn_cnt), "files per each")
            output = p.map(segmtantion, mapped_argument)
        else:
            print("Error: some files are missing!")

        fn_cnt_infilled, fn_int, fn_flags, fn_rrhr_trimmed = (
            sorted(glob("*-cnt_Pinfilled_trimmed.fits")),
            sorted(glob("*int_preprocessed.fits")),
            sorted(glob("*flags_wcs.fits")),
            sorted(glob("*rrhr_preprocessed_trimmed.fits")),
        )

        mapped_argument_star = list(zip(fn_cnt_infilled, fn_int, fn_flags, fn_rrhr_trimmed))
        print("mapped_argument for star finder: ", mapped_argument_star)

        ### --- run the coordinate finder for point sources and flaggd pixels --- ###
        if len(fn_cnt_infilled) == len(fn_int) == len(fn_flags) == len(fn_rrhr_trimmed):
            print(len(fn_cnt_infilled), "files per each")
            output = p.map(starfinder_new, mapped_argument_star)
        else:
            print("Error: some files are missing!")
           
