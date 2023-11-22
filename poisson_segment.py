from astropy.io import fits
from photutils.datasets import make_100gaussians_image
from photutils.background import Background2D, MedianBackground
from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from photutils.segmentation import detect_sources
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import simple_norm
from photutils.segmentation import deblend_sources
from photutils.segmentation import SourceFinder
from photutils.segmentation import SourceCatalog
from astropy.stats import sigma_clipped_stats
from photutils.segmentation import detect_sources, make_2dgaussian_kernel

from photutils.aperture import EllipticalAperture
from matplotlib import colors
from astropy.coordinates import Angle
from astropy.coordinates import Angle
from photutils.aperture import EllipticalAperture

# radius is around 1250 for GALEX data
def zeros_outside(radius, data):
    """
    ### input(s)
    radius: radius of the circle (round 1400 for GALEX data).
    pixels outside of this circle will be 
    replaced with 0s

    data: raw image (2d array) to be cleaned

    ### output(s)
    data: cleaned data
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




def segmentation(cnt_data, npixels_value=10000, th_coeff=0.50):
    """
    ### input(s)
    cnt(2d array, float64): Nanned cnt files

    npixel_value(int64): number of pixels connected. Default = 10,000 for large objects

    th_coeff(float64): coefficient for determing the threshold of pixel value for detection.
    threshold is determined by this vlaue * standard diviation of the image

    ### output(s)
    bool_mask(0s and 1s): segmented mask. 0 = no mask, 1 = mask

    """
    data_orig = cnt_data

    # use sigma-clipped statistics to (roughly) estimate the background
    # background noise levels
    mean, _, std = sigma_clipped_stats(data_orig)

    # subtract the background
    data = data_orig - mean

    ### ---Change parameters here--- ###
    #npixels_value = 10000
    #th_coeff = 0.50  # 0.50

    # detect the sources
    threshold = th_coeff * std #1 * std

    kernel = make_2dgaussian_kernel(7.0, size=21)  # FWHM = 7
    convolved_data = convolve(data, kernel)

    segm = detect_sources(
        convolved_data, threshold, npixels=npixels_value, connectivity=8
    )

    finder = SourceFinder(npixels=npixels_value, progress_bar=False)
    segment_map = finder(convolved_data, threshold= threshold)

    cat = SourceCatalog(data_orig, segment_map, convolved_data=convolved_data)

    # Extract positions(x/ycentroid, flux, and eliptical aperture information),
    columns = ["label", "xcentroid", "ycentroid", "segment_flux", "kron_flux", "kron_aperture"]
    tbl = cat.to_table(columns=columns)
    # Convert to pandas dataframe 

    #df = tbl.to_pandas()
    #df_dropped = df[(df["kron_aperture"].a < 500) & (df["kron_aperture"].a < 500)]

    norm = simple_norm(data, "sqrt")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 12.5))
    ax1.imshow(data_orig, origin="lower", norm=colors.LogNorm(vmin=0.1, vmax=1300))
    ax1.set_title("Data", fontsize=16)
    ax2.imshow(np.array(segm), origin="lower", cmap=segm.cmap)  #
    ax2.set_title("Segmentation Image", fontsize=16)

    final_mask = segm
    
    mask_only = np.zeros(segm.shape)
    canvas = np.zeros(shape=(3840, 3840))
    counter = 0
    for i in range(len(tbl)):
        if tbl['segment_flux'][i] > 900:
            a, b = tbl["kron_aperture"][i].a, tbl["kron_aperture"][i].b

            #if tbl["kron_aperture"][i].a < 500 or tbl["kron_aperture"][i].b < 500:
            if a < 200 and b < 200:
                if tbl['segment_flux'][i] < 100000:
                    ap_size = 1.5
                    print("Flux ", tbl['segment_flux'][i], "aperture size ", ap_size)
                elif tbl['segment_flux'][i] > 100000:
                    ap_size = 2.5
                    print("Flux ", tbl['segment_flux'][i], "aperture size ", ap_size)

                aperture = EllipticalAperture(
                    tbl["kron_aperture"][i].positions,
                    a=ap_size* tbl["kron_aperture"][i].a,
                    b=ap_size*1.5* tbl["kron_aperture"][i].b,
                    theta=tbl["kron_aperture"][i].theta,
                )

                print(
                "original a & b: ", tbl["kron_aperture"][i].a, tbl["kron_aperture"][i].b
            )
                print(
                "enlarged a & b: ", ap_size *tbl["kron_aperture"][i].a, ap_size * tbl["kron_aperture"][i].b
            )

                mask = aperture.to_mask()
                ap_patches = aperture.plot(color="white", lw=1, ax=ax1)
                ap_patches = aperture.plot(color="white", lw=1, ax=ax2)

                subimage = mask

                x, y = (
                    tbl["kron_aperture"][i].positions[1],
                    tbl["kron_aperture"][i].positions[0],
                )
                width, height = np.shape(mask)[0], np.shape(mask)[1]

                top_left_x = int(x - (width / 2))
                top_left_y = int(y - (height / 2))
                bottom_right_x = top_left_x + width
                bottom_right_y = top_left_y + height

                canvas[top_left_x:bottom_right_x, top_left_y:bottom_right_y] = mask

                final_mask += canvas
                counter += 1

    if counter == 0:
        final_mask = canvas

    # Trim all the pixels that went over the circle
    final_mask = zeros_outside(radius=1400, data=final_mask)

    bool_mask = np.where(final_mask != 0, 1, final_mask)
    fig, ax = plt.subplots()
    plt.imshow(bool_mask, origin="lower")
    plt.title("final mask, aperture + segment")
    plt.clim(0, 1)
    plt.colorbar()

    return bool_mask


def poisson_noise(cnt, rrhr, mask, skybg):
    """
    ### input(s)
    cnt(2d array, float64): Nanned cnt files

    mask(2d array, float64): mask generated using segmentation

    skybg(2d array, float64): .objmask.fits obtained from GALEX

    ### output(s)
    cnt_noise: cnt data, masked and poisson noise added
    """

    data = skybg * rrhr  # New pixels for the segmented regions

    segmented_bg = np.where(mask == 0, 0, data)  # remove outside the segments

    # Feed background values from skybg and generate poisson distribution
    
    noise_added = np.random.poisson(lam=segmented_bg)

    # Replace the pixels in cnt with the generated noise
    cnt_noise = np.where(segmented_bg != 0, noise_added, cnt)

    fig, ax = plt.subplots()
    plt.imshow(cnt_noise, origin="lower", norm=colors.LogNorm(vmin=0.1, vmax=1300))
    plt.title("final image, noise added")
    plt.colorbar()

    # Return the final image
    return cnt_noise
