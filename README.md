<h1 align="center"> 🌌 Space Image Processing with High-Performance Computing 🌌 </h1>
I am developing Python & HPC algorithms for processing space images taken by the GALEX space telescope. I am aiming to remove all the bright sources from individual images and create a smooth all-sky map that covers the entire sky so we can closely study the galactic dust and gas.

<h2 align="center"> Dataset </h2>

[*GALEX sky survey archive*](https://archive.stsci.edu/missions-and-data/galex). I queried images from this archive using Python.

<h2 align="center"> Code </h2>

- [**Cleaning.py**](https://github.com/hina0830g/GALEX_montage/blob/dev/cleaning.py) <be>
A module for cleaning and smoothing the data. nan_outside replaces all the pixels outside a radius with nan in an image (2d array). gaussian_filter_2d applies a 2d Gaussian filter to an image to smooth it out. Combined performs both; gaussian filter first and cleans the edges using nan_outside.
