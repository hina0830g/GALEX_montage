from astropy.io import fits
import numpy as np
import scipy
from scipy import stats
from matplotlib import colors
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve, convolve_fft


# radius is around 1250 for GALEX data
def nan_outside(radius, data):
    """
    ### input(s)
    radius: radius of the circle (round 1250 for GALEX data).
    pixels outside of this circle will be cleaned and
    replaced with nan for faster processing

    data: raw image (2d array) to be cleaned

    ### output(s)
    data: smoothed/filtered image
    """

    # Create a grid of coordinates that match the dimension of the image
    x, y = np.arange(0, len(data)), np.arange(0, len(data[0]))
    x_grid, y_grid = np.meshgrid(x, y)

    ## Calculate distances for all pixels at once

    # Find the center of the image in x and y
    x_cent, y_cent = int(round(len(data) / 2)), int(round(len(data[0]) / 2))
    # Distance formula
    distances = np.sqrt((x_cent - x_grid) ** 2 + (y_cent - y_grid) ** 2)
    # Create a mask for pixels outside the circle
    mask = distances > radius
    # Apply the mask ( masked pixels = nan )
    print("r = ", radius, ", cleaned.")
    data[mask] = np.nan

    return data


def gaussian_filter_2d(sigma, data):
    """
    ### input(s)
    sigma: kernel size for the gaussian filter
    data: raw image (2d array) to be smoothed

    ### output(s)
    data: smoothed/filtered image
    """

    gauss_kernel = Gaussian2DKernel(sigma, x_size=75, y_size=75)
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
