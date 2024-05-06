<h1 align="center"> 🌌 Space Image Processing with High-Performance Computing 🌌 </h1>

This is a Python & Julia based astronomical image processing project, parallelized and designed to run in an HPC environment. We aim to remove all bright sources from individual images and create an all-sky map of FUV galactic dust. The final map will be stored as  a hierarchy of FITS files called HiPS (Hierarchical Progressive Survey). 

<h2 align="center"> Dataset </h2>

Images are queried from [*GALEX sky survey archive*](https://archive.stsci.edu/missions-and-data/galex) using astropy MAST query in Python. The types of files queried are the following: cnt.fits, rrhr.fits, skybg.fits, int.fits, and flags.fits.


<h2 align="center"> Installation </h2>

**Requirements**

- Python 3.8 or newer
- Julia Environment 
- HPC Environment ( I request 470 GB of RAM, however, this can be reduced)

To run the program, clone this repository by running the following command:

``
git clone thhps://github.com 
``

Then run the following line:

``
sudo python xx.py install
``


<h2 align="center"> Getting Started </h2>


- **Create a Python virtual environment for Python 3.8 or newer**  

- **Create/Modify the input file**  

An input file should include the target coordinates.  Each line should include 2 integers with a single space in-between, Right Ascension followed by Declination in degrees. 


Example:

195 21  
195 26

This would submit an array job at RA=195, DEC=21 & RA=195, DEC=21. 

- **argparse.slurm**

Change the path to your input file path:   

CurrentCoordinates="$( sed "${SLURM_ARRAY_TASK_ID}q;d" **/path/to/input** )"



<h2 align="center"> Codes </h2>

- [**query.py**](query.py)  
This code reads an input file (a list of coordinates) and performs a criteria-based query at each coordinate. A new directory is created for each coordinate, and all the downloaded files get transferred to the designated directory.

- [**preprocessing.py**]  
This code preprocesses 3 types of files, cnt.fits, rrhr.fits, and int.fits. Blank pixels (0s) in all files are replaced with np.nan in order to speed up the future calculation. These pixels are found outside the r~1400 [pix] of the images. Negative pixels found in rrhr and int files also get replaced with np.nan. Finally, a Gaussian filter is applied to int.fits files to smooth out the images. The default parameters for the Gaussian filter is fwhm=7 and kernel size=21 pixels.

- [**mask_step1_ver3.py**] file name tbd  
The primary purposes of this code are to:
1. Perform segmentation on the images to identify large, bright objects (extended sources + galaxies) and create masks that cover them (function segmentation)
2. Generate Poisson noise from the background values and fill in the masked region to eliminate the bright sources (function segmentation)
3. Run a star finder to identify the location of the remaining bright sources; these are point sources aka stars. It stores the list of coordinates as a .pkl file and masked image and star mask as FITS. (function starfinder)  
4. 
In every step, flags.fits files are used in order to flag artifacts (=bad pixels).

- [**Cleaning.py**](https://github.com/hina0830g/GALEX_montage/blob/dev/cleaning.py) <be>
A module for cleaning and smoothing the data. nan_outside replaces all the pixels outside a radius with nan in an image (2d array). gaussian_filter_2d applies a 2d Gaussian filter to an image to smooth it out. Combined performs both; gaussian filter first and cleans the edges using nan_outside.
