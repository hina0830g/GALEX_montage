<h1 align="center"> 🌌 Space Image Processing with High-Performance Computing 🌌 </h1>

This is a Python & Julia based astronomical image processing project, parallelized and designed to run in an HPC environment. We developed a multi-step algorithm that aims to identify all bright objects (galaxies, stars, and other extended/saturated sources) from individual images and cleanly fill them in with predicted background pixels. Processed images are montaged together to create an all-sky map of FUV galactic dust. The final map will be stored as a hierarchy of FITS files called HiPS (Hierarchical Progressive Survey). 

<h2 align="center"> Project Goals </h2>

- Updating and improving the previously published work (Hamden et al 2013); optimized for HPC, written in a more modern language, multi-step approach of object identification, and usage of machine learning
- Criteria-based query for FUV data extraction
- Identification of large celestial objects (extended sources and galaxies) via segmentation
- Identification of remaining point sources via 2D Gaussian fit 
- Clean removal and infill of the identified sources 
- Creation of Far UV all-sky map for public access and future science
- Parallelizing the scripts for faster and more efficient image processing to process large data
- Usage argparse to pass a list of coordinates to HPC server and automate the slurm job schedulings

<h2 align="center"> Dataset </h2>

Images are queried from [*GALEX sky survey archive*](https://archive.stsci.edu/missions-and-data/galex) using astropy MAST query in Python. The types of files queried are the following: cnt.fits, rrhr.fits, skybg.fits, int.fits, and flags.fits.


<h2 align="center"> Installation </h2>

**Requirements**

- Python 3.8 or newer
- Julia Environment 
- HPC Environment (I request 470 GB of RAM, however, this can be reduced)

First, create a Python virtual environment for Python 3.8 or newer in an HPC system:  

```
module load python/<version>
virtualenv --system-site-packages </path/to/virtual/env>
```
Then activate the virtual environment (source venv/bin/activate).


To run the program, clone this repository by running the following command:

```
git clone dev https://github.com/hina0830g/GALEX_montage.git
```

<h2 align="center"> Getting Started </h2>

- **Create a new new directory "raw_files" in your virtual environment**  

query.py will try to locate your raw_files directory and make a new directory inside of it where all the files will be downlaoded. The format of the folder name is DATE-RAvalue-DECvalue

- **Create/Modify the input file**  

An input file should include the target coordinates.  Each line should include 2 integers with a single space in-between, Right Ascension followed by Declination in degrees. 


Example:
```
195 21  
195 26
```
This would submit an array job at RA=195, DEC=21 & RA=195, DEC=21. 

- **Open argparse.slurm and change the path to your input file path**     

CurrentCoordinates="$( sed "${SLURM_ARRAY_TASK_ID}q;d" **/path/to/input** )"

- **Open query.py and preprocessing.py and change the home directory path**    

home_dir = "path/to/your/virtualenvironment"

- **import packages**

To import poisson_segments.py and cleaning.py as packages, first locate the files in Python:  
```
sys.path.append("/path/to/packages")
```

Then import the packages:  

```
import poisson_segment as ps  
import cleaning as cl
```

<h2 align="center"> Codes </h2>

- [**argparse.slurm**](https://github.com/hina0830g/GALEX_montage/blob/dev/HPC_scripts/argparse.slurm)  
SLURM script for submitting an array job in the following order: query.py -> preprocessing.py -> mask_step1_ver3.py

- [**query.py**](https://github.com/hina0830g/GALEX_montage/blob/dev/HPC_scripts/query.py)  
This code reads an input file (a list of coordinates) and performs a criteria-based query at each coordinate. A new directory is created for each coordinate, and all the downloaded files get transferred to the designated directory.

- [**preprocessing.py**](https://github.com/hina0830g/GALEX_montage/blob/dev/HPC_scripts/preprocessing.py)
This code preprocesses 3 types of files, cnt.fits, rrhr.fits, and int.fits. 0s in all files are replaced with Nans to speed up the future calculation. These pixels are found outside the r~1400 [pix] of the images. Negative pixels in rrhr and int files also get replaced with Nans. Finally, a Gaussian filter is applied to int.fits files to smooth out the images before running segmentation. The default parameters for the Gaussian filter is fwhm=7 and kernel size=21 pixels.

- [**object_detection.py**] file name tbd  
The primary purposes of this code are to:
1. Perform segmentation on the images to identify large, bright objects (extended sources + galaxies) and create masks that cover them (function segmentation)
2. Generate Poisson noise from the background values and fill in the masked region to eliminate the bright sources (function segmentation)
3. Run a star finder to identify the location of the remaining bright sources; these are point sources aka stars. It stores the list of coordinates as a .pkl file and masked image and star mask as FITS. (function starfinder)  

In every step, flags.fits files are used in order to flag artifacts (=bad pixels).

- [**pointsource_infill.jl**](https://github.com/hina0830g/GALEX_montage/blob/dev/HPC_scripts/pointsource_infill.jl)
This code uses a Julia infill package CloudClean to fill in point sources and artifact pixels that were flagged in the previous step.

- [**post processing.py**](https://github.com/hina0830g/GALEX_montage/blob/dev/HPC_scripts/)
This code cleans the output images of point_source_infill.jl. It eliminates any undesired pixels outside of the actual images and transfers the files.

- [**montage.py**](https://github.com/hina0830g/GALEX_montage/blob/dev/HPC_scripts/montage.py)  
This code combines all the processed images to create a large mosaic. The primary function it uses is the coadd function from MontagePy and takes the mean for the regions where images overlap. The mosaic is saved as uncorrect.fits.

- [**Cleaning.py**](https://github.com/hina0830g/GALEX_montage/blob/dev/cleaning.py)  
A module for cleaning and smoothing the data. nan_outside replaces all the pixels outside a radius with nan in an image (2d array). gaussian_filter_2d applies a 2d Gaussian filter to an image to smooth it out. Combined performs both; gaussian filter first and cleans the edges using nan_outside.

