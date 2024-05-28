from datetime import date
import argparse
import ast
import sys
import MontagePy as Montage
from MontagePy.main import *
from MontagePy.archive import *
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import os
import shutil
from multiprocessing import Pool
from glob import glob 

import sys
sys.path.append("/xdisk/hamden/hina0830/venv39/")
import cleaning as cl

home_dir = "/xdisk/hamden/hina0830/venv39"
num_cpus=int(os.environ["SLURM_CPUS_ON_NODE"])
print("Running a pool with %s workers"%num_cpus)

os.chdir(home_dir + "/raw_files")
print("current path: ", os.getcwd())
today = date.today()
print("Today's date:", today)

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
    c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='icrs')
    
    try:
        # Switch to the new directory
        dir_name = str(today) + "-RA" + str(RA) + "-DEC" + str(DEC)
        os.chdir(dir_name)
        path_init = os.getcwd()
        print("Path changed to: ", os.getcwd())    
        
    except FileNotFoundError:
        print("Error")
        
else:
    # Switch to the new directory
    dir_name = "2024-05-15-RA180-DEC12" # Manually enter a path

    # Get coordinates from the directory name
    i = dir_name.find("RA") # Find the starting position of "RA"
    RA_string, DEC_string = astronomy_string[i:].split('-')
    RA, DEC = int(RA_string.removeprefix("RA")), int(DEC_string.removeprefix("DEC"))
    
    c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='icrs')
    
    os.chdir(dir_name)
    path_init = os.getcwd()
    print("Path changed to: ", os.getcwd())
    
home = os.getcwd()

# Clean the edge of draw.fits and save them using headers from intensity files
draw_files, int_files = sorted(glob("*draw.fits")), sorted(glob("*-int_Pinfilled_trimmed.fits"))
print(len(draw_files), len(int_files))

os.makedirs("raw", exist_ok=True)

for i in range(len(draw_files)):
    hdu = fits.open(int_files[i])
    data = fits.getdata(draw_files[i])
    
    result_cleaned = cl.nan_outside(1400, data)
    hdu[0].data = result_cleaned 
    
    filename = draw_files[i].removesuffix(".fits") + "_cleaned.fits"
    
    hdu.writeto(
    path_init + "/raw/" + filename,
    overwrite=True,
)


os.chdir( path_init )

# These are the parameters defining the mosaic we want to make
location = c.to_string('hmsdms')
size     = 8.0
dataset  = "GALEX"
workdir  = dir_name # where raw directory is locat$

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

os.makedirs("projected", exist_ok=True)
os.makedirs("diffs", exist_ok=True)
os.makedirs("corrected", exist_ok=True)

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
rtn = mAdd("projected", "pimages.tbl", "region.hdr", "uncorrected.fits", coadd = 0)
print("mAdd: " + str(rtn), flush=True)
