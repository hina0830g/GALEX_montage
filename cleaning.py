from astropy.io import fits
import numpy as np
import scipy
from scipy import stats
from matplotlib import colors
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve, convolve_fft
import pandas as pd

# radius is around 1250 for GALEX data
def nan_outside(radius, data):
    """
    Replace all pixels outside of the bounds of the circle
    with numpy Nans.
    
    Parameters
    ----------
    radius : int
        The radius of the bounds of the circe 
        
    data : numpy array (float 64)
        The image to be edited
        
    Returns
    -------
    data : numpy array (float 64)
        The final image (edited)

    Description
    -----------
    This function cleans an image by defining a radius 
    and replacing outside the boundary with numpy nans.
    
    """   

    # Create a grid of coordinates that match the dimension of the image
    x, y = np.arange(0, len(data)), np.arange(0, len(data))
    x_grid, y_grid = np.meshgrid(x, y)

    ## Calculate distances for all pixels at once

    # Find the center of the image in x and y
    x_cent, y_cent = int(round(len(data) / 2)), int(round(len(data) / 2))
    # Distance formula
    distances = np.sqrt((x_cent - x_grid) ** 2 + (y_cent - y_grid) ** 2)
    # Create a mask for pixels outside the circle
    mask = distances > radius
    # Apply the mask ( masked pixels = nan )
    print("r = ", radius, ", cleaned.")
    data[mask] = np.nan

    return data

def zeros_outside(radius, data):
    """
    Replace all pixels outside of the bounds of the circle
    with 0s.
    
    Parameters
    ----------
    radius : int
        The radius of the bounds of the circe 
        
    data : numpy array (float 64)
        The image to be edited
        
    Returns
    -------
    data : numpy array (float 64)
        The final image (edited)

    Description
    -----------
    This function cleans an image by defining a radius 
    and replacing outside the boundary with 0s.
    
    """   
    
    # Create a grid of coordinates that match the dimension of the image
    x, y = np.arange(0, len(data)), np.arange(0, len(data))
    x_grid, y_grid = np.meshgrid(x, y)

    ## Calculate distances for all pixels at once

    # Find the center of the image in x and y
    x_cent, y_cent = int(round(len(data) / 2)), int(round(len(data) / 2))
    # Distance formula
    distances = np.sqrt((x_cent - x_grid) ** 2 + (y_cent - y_grid) ** 2)
    # Create a mask for pixels outside the circle
    mask = distances > radius
    # Apply the mask ( masked pixels = nan )
    print("r = ", radius, ", cleaned.")
    data[mask] = 0

    return data


def gaussian_filter_2d(sigma, data):
    """
    Replace all pixels outside of the bounds of the circle
    with 0s.
    
    Parameters
    ----------
    sigma : float
         kernel size for the gaussian filter
        
    data : numpy array (float 64)
        The image to be edited
        
    Returns
    -------
    data : numpy array (float 64)
        The final image (edited)

    Description
    -----------
    This function cleans an image by defining a radius 
    and replacing outside the boundary with 0s.
    
    """   
    gauss_kernel = Gaussian2DKernel(sigma, x_size=21, y_size=21)
    data = convolve(data, gauss_kernel)
    print("2D Gaussian filter completed.")
    return data


def combined(sigma, radius, data):

    
    """
    ### input(s)
    radius: radius of the circle (round 1250 for GALEX data).
    pixels outside will be cleaned (=np.nan)

    sigma: kernel size for the gaussian filter
    data: raw image (2d array) to be smoothed

    ### output(s)
    data_final : image smoothed/filtered and edges
    cleaned
    """
    data_gf = gaussian_filter_2d(sigma, data)
    data_final = nan_outside(radius, data_gf)
    return data_final
