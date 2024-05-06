<h1 align="center"> 🌌 Space Image Processing with High-Performance Computing 🌌 </h1>

This is a Python & Julia based astronomical image processing project, parallelized and designed to run in an HPC environment. We aim to remove all bright sources from individual images and create an all-sky map of FUV galactic dust. The final map will be stored as  a hierarchy of FITS files called HiPS (Hierarchical Progressive Survey). 

<h2 align="center"> Dataset </h2>

Images are queried from [*GALEX sky survey archive*](https://archive.stsci.edu/missions-and-data/galex) using astropy MAST query in Python. The types of files queried are the following: cnt.fits, rrhr.fits, int.fits, skybg.fits, flags.fits.


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

- **Create/Modify the input file**  

An input file should include the target coordinates.  Each line should include Right Ascension followed by Declination in degrees. Must be integers.
The format is RA in degree DEC in degree per line. 

Example:

195 21  
195 26

This would submit an array job at RA=195, DEC=21 & RA=195, DEC=21. 

- **argparse.slurm**

Change the path to your input file path:   

CurrentCoordinates="$( sed "${SLURM_ARRAY_TASK_ID}q;d" **/path/to/input** )"



<h2 align="center"> Codes </h2>

- [**Cleaning.py**](query.py) <be>

- [**Cleaning.py**](https://github.com/hina0830g/GALEX_montage/blob/dev/cleaning.py) <be>
A module for cleaning and smoothing the data. nan_outside replaces all the pixels outside a radius with nan in an image (2d array). gaussian_filter_2d applies a 2d Gaussian filter to an image to smooth it out. Combined performs both; gaussian filter first and cleans the edges using nan_outside.
