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

    x_cent, y_cent = int(round(len(data) / 2)), int(round(len(data) / 2)) # Center of the image in x and y axis
    distances = np.sqrt((x_cent - x_grid) ** 2 + (y_cent - y_grid) ** 2) # Calculate distances for all pixels at once

    mask = distances > radius
    data[mask] = np.nan # masked regions are now Nans

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

    x_cent, y_cent = int(round(len(data) / 2)), int(round(len(data) / 2))  # Center of the image in x and y axis
    distances = np.sqrt((x_cent - x_grid) ** 2 + (y_cent - y_grid) ** 2) # Calculate distances for all pixels at once
   
    mask = distances > radius # Create a bool mask that indicates which pixels are outside the circle
    data[mask] = 0 # masked regions are now 0s

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
    gauss_kernel = Gaussian2DKernel(sigma)#, x_size=21, y_size=21)
    data = convolve(data, gauss_kernel)
    print("2D Gaussian filter completed.")
    return data


def combined(sigma, radius, data):
    """
    Combination of gaussian_filter_2d and nan_outside.
    
    Parameters
    ----------
    sigma : float
         kernel size for the gaussian filter
         
    radius : int
        The radius of the bounds of the circe 
        
    data : numpy array (float 64)
        The image to be edited
        
    Returns
    -------
    data_final : numpy array (float 64)
        The final image (edited)

    Description
    -----------
    This function Gaussian blurs an image
    and cleans it afterward.
    """  
    
    data_gf = gaussian_filter_2d(sigma, data)
    data_final = nan_outside(radius, data_gf)
    return data_final
