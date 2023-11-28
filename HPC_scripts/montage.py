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

num_cpus=int(os.environ["SLURM_CPUS_ON_NODE"])
print("Running a pool with %s workers"%num_cpus)

# Swicth to the raw file directory
path_init = os.getcwd() + '/raw_files'
os.chdir( path_init )
print('initial path (path_init): ', path_init)
home = os.getcwd()

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
    try:
        # Specify the coordinates using the input information
        c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='icrs')
        dir_name = str(today) + "-RA" + str(RA) + "-DEC" + str(DEC)
        
    except FileNotFoundError:
        print("Error")


# These are the parameters defining the mosaic we want to make
location = c.to_string('hmsdms')
size     = 10.0
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

def montage(type):
    if type == "cnt":
        os.chdir( home + '/' + workdir + '/' + "cnt")
        print("cnt cwd: ", os.getcwd() )
    elif type == "rrhr":
        os.chdir( home + '/' + workdir + '/' + "rrhr" )
        print("rrhr cwd: ", os.getcwd() )

    os.makedirs("projected")
    os.makedirs("diffs")
    os.makedirs("corrected")

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
    rtn = mAdd("projected", "pimages.tbl", "region.hdr", "uncorrected.fits", coadd = 3)
    print("mAdd: " + str(rtn), flush=True)

    # Read and return the result (mosaic)
    hdu = fits.open("uncorrected.fits")
    data = (hdu[0].data).astype('float64')

    return data

if __name__=="__main__":
    with Pool(num_cpus) as p:
        type_list = ["cnt", "rrhr"]
        data = p.map(montage, type_list)
