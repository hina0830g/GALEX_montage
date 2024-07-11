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


**Data Products used**


| Filename | Units | size | Description |
| :---: | --- | ------- | --- |
| fd-cnt.fits | counts/pixel | 3840 x 3840 | The raw number of counts per pixel, not corrected for the exposure time or flat field |
| fd-rrhr.fits| seconds | 3840 x 3840 | The high resolution relative response. This is the rr image linearly interpolated to the same pixel scale as the cnt map. |
|fd-int.fits| counts/sec/pixel |  3840 x 3840 | Intensity map (cnt / rrhr) |
| fd-skybg.fits | counts/sec/pixel |  3840 x 3840 | The sky background map subtracted from the data before identifying sources. |
| fd-flags.fits | flag value |  480 x 480 | Flag map indicating regions of the map likely contaminated by artifacts or regions where various types of artifacts have been removed. This file type needs to be reprojected in order to match the dimension to other files (preprocessing.fits) |

Image resolution for GALEX: 4" FWHM

Reference: GALEX Chapter 4 - [Imamging Data Products](http://www.galex.caltech.edu/researcher/techdoc-ch4.html)

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
git clone --branch dev https://github.com/hina0830g/GALEX_montage.git

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

- **Specify the path to your input file path in the SLURM file (input_path="/xdisk/hamden/hina0830/venv39/input/coordinates1")**     

input_path="/xdisk/hamden/hina0830/venv39/input/coordinates1"

- **Specify the home directory path in .py files**    

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
SLURM script for submitting an array job in the following order: query.py -> preprocessing.py -> object_detection.py -> pointsource_infill.jl -> query.py

- [**query.py**](https://github.com/hina0830g/GALEX_montage/blob/dev/HPC_scripts/query.py)  
This code reads an input file (a list of coordinates) and performs a criteria-based query at each coordinate. A new directory is created for each coordinate, and all the downloaded files get transferred to the designated directory.

- [**preprocessing.py**](https://github.com/hina0830g/GALEX_montage/blob/dev/HPC_scripts/preprocessing.py)  
This code uses the flags files to determine the quality of the observation. If more than 20% of the image is contaminated ( above 127 in flags.fits ) preprocesses the following types of files, cnt.fits, rrhr.fits, int.fits, skybg.fits, and flags.fits.

| File type | execution |
| :---: | --- |
| cnt.fits |  Replace the pixels outside of the circle with Numpy Nans. Cut 500 blank pixels from each side so 3840 x 3840 -> 2840 x 2840 |
| rrhr.fits |  Replace negative values with Numpy Nans. Cut 500 blank pixels from each side. |
| int.fits | Replace negative values with Numpy Nans. Apply a Gaussian smoothing filter (parameters: fwhm=7 & kernel size=21 in pixels). Cut 500 blank pixels from each side. |
| skybg.fits | Cut 500 blank pixels from each side of the image. |
| flags.fits | Reproject and resize the file to be 480 x 480 -> 3840 x 3840. Cut 500 blank pixels from each side.  |

- [**object_detection.py**](https://github.com/hina0830g/GALEX_montage/blob/dev/HPC_scripts/object_detection.py) 
The primary purposes of this code are to:
1. Perform segmentation on the images to identify large, bright objects (extended sources + galaxies) and create masks that cover them (function segmentation)
2. Generate Poisson noise from the background values and fill in the masked region to eliminate the bright sources (function segmentation)
3. Run a star finder to identify the location of the remaining bright sources; these are point sources aka stars. It stores the list of coordinates as a .pkl file and masked image and star mask as FITS. (function starfinder)  

In every step, flags.fits files are used in order to flag artifacts (=bad pixels).

- [**pointsource_infill.jl**](https://github.com/hina0830g/GALEX_montage/blob/dev/HPC_scripts/pointsource_infill.jl)
This code uses a Julia infill package CloudClean to fill in point sources and artifact pixels that were flagged in the previous step.

- [**montage.py**](https://github.com/hina0830g/GALEX_montage/blob/dev/HPC_scripts/montage.py)  
This code combines all the processed images to create a large mosaic. The primary function it uses is the coadd function from MontagePy and takes the mean for the regions where images overlap. The mosaic is saved as uncorrected.fits.

