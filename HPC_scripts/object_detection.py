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

#os.chdir("2024-05-14-RA180-DEC12")
print("New location: ", os.getcwd())


def segmtantion(file_tuple):
    # Unpack the tuple and get lists of file names
    cnt_file, rrhr_file, skybg_file, int_file, flags_file = file_tuple
    
    # Read all the files and get data
    cnt_data, rrhr_data, skybg_data, int_data, flags_data = fits.getdata(cnt_file).astype('float64'), fits.getdata(rrhr_file).astype('float64'), fits.getdata(skybg_file).astype('float64'), fits.getdata(int_file).astype('float64'), fits.getdata(flags_file).astype('float64')

    print("processing: ", file_tuple)

    # Read cnt file here to get hdu
    hdu = fits.open(cnt_file)
    hdr = hdu[0].header
    
    fn_current = cnt_file.removesuffix('-cnt_nan.fits') # fn as filename
    
    exptime, project_name = hdr['EXPTIME'], hdr['MPSTYPE'] 
    print(fn_current, exptime, project_name)
    
    ### The default parameters ###
    edge_thickness = 500 # Trim 500 (blank) pixels on each side
    radius = 1400
    sigma = 3.0
    npixels_value = 500
    
    # If statements to determine the coefficient based on 
    # exposure time
    if exptime >= 30000: # DIS
        th_coeff = 15
    elif exptime >= 1500: # MIS & NGS
        th_coeff = 10
    elif exptime >= 100: # "AIS
        th_coeff = 6
    else:
        print("exception; exposure time is: ", exptime)
        th_coeff = 6
        
    #int_data[flags_data  >= 256] = np.nans
     
    mask = ps.segmentation(fn_current, int_data, npixels_value, th_coeff) # Perform segmentations on large sources
    print("ps.segmentation finished")

    cnt_noise = ps.poisson_noise(cnt_data, rrhr_data, mask, skybg_data) # Perform Poisson infill
    print("ps.poisson_noise finished ")
    
    # Nan out bad pixels?
    #cnt_noise[flags_data  >= 256] = np.nan
    
    # Trim each side (empty pixels) of the array by the edge_thickness (=500pix)
    cnt_infilled_trimmed = cnt_noise[edge_thickness:len(cnt_noise)-edge_thickness, edge_thickness:len(cnt_noise[0])-edge_thickness]
    rrhr_trimmed = rrhr_data[edge_thickness:len(rrhr_data)-edge_thickness, edge_thickness:len(rrhr_data[0])-edge_thickness]
    
    #flags_trimmed = flags_data[edge_thickness:len(flags_data)-edge_thickness, edge_thickness:len(flags_data[0])-edge_thickness]
    
    # Apply Gaussian filter and clean the edges
    final_cnt = cl.combined(sigma=sigma, radius=radius, data=cnt_infilled_trimmed)
    print("cl.combined finished")
    
    #final_cnt[flags_trimmed  >= 256] = np.nan

    # Save the new files
    hdu[0].data = final_cnt
    hdu.writeto(
        cnt_file.removesuffix("_preprocessed.fits") + "_Pinfilled_trimmed.fits",
        overwrite=True,
    )
    print(f"Processing completed for {cnt_file}")
    
    hdu[0].data = rrhr_trimmed
    hdu.writeto(
        rrhr_file.removesuffix(".fits") + "_trimmed.fits",
        overwrite=True,
    )
    return cnt_data

def starfinder_new(file_tuple):
    cnt_file, int_file, flags_file, rrhr_file = file_tuple
    
    edge_thickness = 500 # Trim 500 (blank) pixels on each side
    mask_size = 16 # size of mask for each star in pixels

    th_coeff = 25.0 # for AIS
    DAO_fwhm = 10.385
    
    hdu = fits.open(cnt_file)
    
    cnt_data, rrhr_data, flags_data  = (hdu[0].data).astype("float64"), fits.getdata(rrhr_file).astype('float64'), fits.getdata(flags_file).astype('float64')
    
    flags_data = flags_data[edge_thickness:len(flags_data)-edge_thickness, edge_thickness:len(flags_data[0])-edge_thickness]
    
    divided_data = cnt_data/rrhr_data 

    divided_data[flags_data  != 0] = np.nan # flag the data 
    # Drops nans
    divided_data[np.isnan(divided_data)] = 0 # nan = 0

    masked_data, mask_data, coord_lis = ps.psfinder(divided_data, flags_data, th_coeff, DAO_fwhm, mask_size)

    # --- Save files here ---
    
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

def starfinder(file_tuple):
    cnt_file, int_file, flags_file, rrhr_file = file_tuple
    
    edge_thickness = 500 # Trim 500 (blank) pixels on each side
    mask_size = 16 # size of mask for each star in pixels

    th_coeff = 25.0 # for AIS
    DAO_fwhm = 10.385
    inner_r = 1400
    
    hdu = fits.open(cnt_file)
    data = (hdu[0].data).astype("float64") # Get poisson infilled cnt 
    flags_data = fits.getdata(flags_file).astype('float64') # Get flags file
    flags_data = flags_data[edge_thickness:len(flags_data)-edge_thickness, edge_thickness:len(flags_data[0])-edge_thickness]
    
    rrhr_data = fits.getdata(rrhr_file).astype('float64') 
    
    divided_data = data/rrhr_data 

    divided_data[flags_data  != 0] = np.nan # flag the data 
    # Drops nans
    divided_data[np.isnan(divided_data)] = 0 # nan = 0
    
    mean, median, std = sigma_clipped_stats(divided_data)  
    print((mean, median, std)) 

    daofind = DAOStarFinder(fwhm=DAO_fwhm, threshold=th_coeff*std) # fwhm=3.0, 3.0*std 
    sources = daofind(divided_data - median)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=4.0)
    norm = ImageNormalize(stretch=SqrtStretch())

    # Define custom radii for the stars 
    star_radii = np.ones(len(sources)) * mask_size

    # Mask stars by replacing pixels with 0s
    mask_data = np.copy(divided_data)

    # create an epmty array, same size as the image
    canvas = np.zeros(shape=(divided_data.shape[0], divided_data.shape[1]))

    y_indices, x_indices = np.indices(data.shape)
    x_cent, y_cent = int(round(len(divided_data) / 2)), int(round(len(divided_data) / 2))
    coord_lis = []

    for star, radius in zip(sources, star_radii):
        x, y = int(star['xcentroid']), int(star['ycentroid'])

        # Distances between the star and every pixel in the image
        dist = np.sqrt((x_indices - x) ** 2 + (y_indices - y) ** 2) 

        # Find the distance between the center of the image and the star
        star_dist = np.sqrt((x_cent - x) ** 2 + (y_cent - y) ** 2)

        # Ignore the stars near the edge (outside the inner radius)
        if star_dist <= inner_r: 
            # fill in the pixel around the star and mask the star in the image and the mask
            # the size of masking for each star is mask_size defined at the top
            # masking the image; masked region = 0
            mask_data[dist <= radius] = 0

            # mask file; masked = 1 and nonmasked = 0
            canvas[dist <= radius] = 1

            coord = [x, y]
            coord_lis.append(coord)

            
    #mask_data[flags_data  != 0] = np.nan # flag the data 
   
    #mask_data[np.isnan(mask_data)] = 0 # nan = 0

    masked_data = mask_data

    mask_data = np.where(mask_data==0, 1, canvas)
    print(len(coord_lis), " stars detected")
    
    
    ######### Flags segmentation

    npixels_value_wcs = 25
    threshold_wcs = 50
    
    segm = detect_sources(
    flags_data, threshold=threshold_wcs, npixels=npixels_value_wcs, connectivity=8
    )
    
    print("segm: ", segm)
    print(type(segm))
    print("detect_sources done.")
    finder = SourceFinder(npixels=npixels_value_wcs, progress_bar=False, deblend=False)
    segment_map = finder(flags_data, threshold=threshold_wcs)

    canvas = np.zeros(shape=(len(flags_data), len(flags_data[0])))

    if segm is not None:
        cat = SourceCatalog(flags_data, segment_map)

        # Extract positions(x/ycentroid, flux, and eliptical aperture information),
        columns = ["label", "xcentroid", "ycentroid", "segment_flux", "kron_flux", "kron_aperture", "semimajor_sigma", "semiminor_sigma", "orientation"]

        tbl = cat.to_table(columns=columns)

        # Assuming tbl is your QTable, you can convert it to a DataFrame if needed
        df = tbl.to_pandas()

        # Find the index of the row with the largest value of "semimajor_sigma"
        indices_to_drop = df[df['semimajor_sigma'] > 50].index

        # Drop the row using the index
        df.drop(indices_to_drop, inplace=True)

        # Convert back to QTable if needed
        tbl = QTable.from_pandas(df)

        mask_shape = flags_data.shape
        mask = np.zeros(mask_shape, dtype=bool)

        # Loop over ellipses and draw them onto the mask
        for i in range(len(tbl)):
            x, y = int(tbl['xcentroid'][i]), int(tbl['ycentroid'][i])
            # Find the distance between the center of the image and the star
            star_dist = np.sqrt((x_cent - x) ** 2 + (y_cent - y) ** 2)

            if star_dist < 1400: # Save the coordinates if the coordinates are inside the circle
                coord = [int(tbl['xcentroid'][i]), int(tbl['ycentroid'][i])]
                coord_lis.append(coord)

    # Save a list of coordinates of stars
    
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
            sorted(glob('*int_preprocessed.fits')),
            sorted(glob('*flags_wcs.fits')),
        )
        mapped_argument = list(zip(fn_cnt, fn_rrhr, fn_skybg, fn_int, fn_flags))
        print("mapped_argument for for process_file : ", mapped_argument )
        

        if len(fn_cnt) == len(fn_rrhr) == len(fn_skybg) == len(fn_skybg) == len(fn_flags):
            print(len(fn_cnt), "files per each")
            output = p.map(segmtantion, mapped_argument)
                       
        fn_cnt_infilled, fn_int, fn_flags, fn_rrhr_trimmed = (
            sorted(glob("*-cnt_Pinfilled_trimmed.fits")),
            sorted(glob('*int_preprocessed.fits')),
            sorted(glob('*flags_wcs.fits')),
            sorted(glob('*rrhr_preprocessed_trimmed.fits')),
        )

        
        mapped_argument_star = list(zip(fn_cnt_infilled, fn_int, fn_flags, fn_rrhr_trimmed))
        print("mapped_argument for star finder: ", mapped_argument_star )

        if len(fn_cnt_infilled) == len(fn_int) == len(fn_flags) == len(fn_rrhr_trimmed):
            print(len(fn_cnt_infilled), "files per each")
            output = p.map(starfinder_new, mapped_argument_star)

