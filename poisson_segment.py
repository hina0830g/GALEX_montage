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
from photutils.aperture import EllipticalAperture

from matplotlib.patches import Ellipse

import gc

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

#def segmentation_ver2(fn_current, data_orig, npixels_value=10000, th_coeff=0.50):
def segmentation_ver2(fn_current, data_orig, npixels_value, th_coeff):
    """
    ### input(s)
    cnt(2d array, float64): Nanned cnt files

    npixel_value(int64): number of pixels connected. Default = 10,000 for large objects

    th_coeff(float64): coefficient for determing the threshold of pixel value for detection.
    threshold is determined by this vlaue * standard diviation of the image

    ### output(s)
    bool_mask(0s and 1s): segmented mask. 0 = no mask, 1 = mask

    """
    
    fig, ax = plt.subplots()
    plt.imshow(data_orig, origin="lower", norm=colors.LogNorm(vmin=1e-4, vmax=1e-2))
    plt.title(fn_current + " (flagged image, original)") 
    plt.colorbar()

    #kernel = make_2dgaussian_kernel(7.0, size=21)  # FWHM = 7
    #convolved_data = convolve(data_orig, kernel)

    # use sigma-clipped statistics to (roughly) estimate the background
    # background noise levels
    mean, _, std = sigma_clipped_stats(convolved_data)

    print("numpy std: ", np.nanstd(convolved_data))
    
    # subtract the background
    data = convolved_data - mean

    # detect the sources
    threshold = th_coeff * std 
    
    print("standard dev: ", std)
    print("threshold: ", threshold)
    print("n_pixels: ", npixels_value)

    segm = detect_sources(
        convolved_data, threshold, npixels=npixels_value, connectivity=8
    )
    
    print("segm: ", segm)
    print(type(segm))
    print("detect_sources done.")
    finder = SourceFinder(npixels=npixels_value, progress_bar=False, deblend=False)
    segment_map = finder(convolved_data, threshold= threshold)

    canvas = np.zeros(shape=(len(data), len(data[0])))
    
    if segm is not None:
        cat = SourceCatalog(data_orig, segment_map, convolved_data=convolved_data)

        # Extract positions(x/ycentroid, flux, and eliptical aperture information),
        columns = ["label", "xcentroid", "ycentroid", "segment_flux", "kron_flux", "kron_aperture", "semimajor_sigma", "semiminor_sigma", "orientation"]

        tbl = cat.to_table(columns=columns)

        mask_shape = data_orig.shape
        mask = np.zeros(mask_shape, dtype=bool)

        # Loop over ellipses and draw them onto the mask
        for i in range(len(tbl)):
            ellipse = Ellipse((tbl['xcentroid'][i], tbl['ycentroid'][i]), 
                              width=10 * tbl['semimajor_sigma'][i].value,  
                              height=10 * tbl['semiminor_sigma'][i].value, 
                              angle=tbl['orientation'][i].value, edgecolor = 'red', facecolor= 'none')
            
            # Convert the ellipse into a binary mask
            x, y = ellipse.get_verts().T
            x_min, x_max = int(x.min()), int(x.max())
            y_min, y_max = int(y.min()), int(y.max())
            xx, yy = np.meshgrid(np.arange(x_min, x_max), np.arange(y_min, y_max))
            ellipse_mask = ellipse.contains_points(np.vstack((xx.flatten(), yy.flatten())).T).reshape((y_max-y_min, x_max-x_min))

            # Get the bounding box of the ellipse in the original image coordinates
            x0, y0 = max(0, int(tbl['xcentroid'][i] - 0.5 * ellipse.width)), max(0, int(tbl['ycentroid'][i] - 0.5 * ellipse.height))
            x1, y1 = min(mask_shape[1], x0 + ellipse_mask.shape[1]), min(mask_shape[0], y0 + ellipse_mask.shape[0])

            # Add the ellipse mask to the main mask
            mask[y0:y1, x0:x1] += ellipse_mask[:y1-y0, :x1-x0]
            ax.add_artist(ellipse)
            
    elif segm is None:
        mask = canvas
        print("No extended sourced detected.")

    plt.close()
    gc.collect()
    return mask



def segmentation(fn, cnt_data, npixels_value=10000, th_coeff=0.50):
    """
    ### input(s)
    cnt(2d array, float64): Nanned cnt files

    npixel_value(int64): number of pixels connected. Default = 10,000 for large objects

    th_coeff(float64): coefficient for determing the threshold of pixel value for detection.
    threshold is determined by this vlaue * standard diviation of the image

    ### output(s)
    bool_mask(0s and 1s): segmented mask. 0 = no mask, 1 = mask

    """
    hdu = fits.open(fn)
    fn_current = fn.removesuffix('-cnt_nan.fits')
    print(fn_current)
    
    data_orig = cnt_data
    
    kernel = make_2dgaussian_kernel(7.0, size=21)  # FWHM = 7
    convolved_data = convolve(data_orig, kernel)

    # use sigma-clipped statistics to (roughly) estimate the background
    # background noise levels
    #mean, _, std = sigma_clipped_stats(data_orig, sigma=10)
    mean, _, std = sigma_clipped_stats(convolved_data)

    #print("numpy std: ", np.nanstd(data_orig))
    print("numpy std: ", np.nanstd(convolved_data))
    
    # subtract the background
    #data = data_orig - mean
    data = convolved_data - mean

    # detect the sources
    threshold = th_coeff * std 
    
    print("standard dev: ", std)
    print("threshold: ", threshold)
    print("n_pixels: ", npixels_value)

    #kernel = make_2dgaussian_kernel(7.0, size=21)  # FWHM = 7
    #convolved_data = convolve(data, kernel)
    
    segm = detect_sources(
        convolved_data, threshold, npixels=npixels_value, connectivity=8
    )
    
    print("segm: ", segm)
    print(type(segm))
    print("detect_sources done.")
    finder = SourceFinder(npixels=npixels_value, progress_bar=False, deblend=False)
    segment_map = finder(convolved_data, threshold= threshold)

 
    canvas = np.zeros(shape=(len(data), len(data[0])))
    if segm is not None:
        #cat = SourceCatalog(data_orig, segment_map, convolved_data=convolved_data)
        cat = SourceCatalog(data_orig, segment_map, convolved_data=convolved_data)
        
        # Extract positions(x/ycentroid, flux, and eliptical aperture information),
        columns = ["label", "xcentroid", "ycentroid", "segment_flux", "kron_flux", "kron_aperture"]
        tbl = cat.to_table(columns=columns)
        # Convert to pandas dataframe 

        #df = tbl.to_pandas()
        #df_dropped = df[(df["kron_aperture"].a < 500) & (df["kron_aperture"].a < 500)]

        norm = simple_norm(data, "sqrt")
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 12.5))
        ax1.imshow(data_orig, origin="lower", norm=colors.LogNorm(vmin=1e-4, vmax=5))
        #ax1.imshow(data_orig, origin="lower", norm=colors.LogNorm(vmin=0.1, vmax=1300))
        ax1.set_title(fn_current, fontsize=14)
        ax2.imshow(np.array(segm), origin="lower", cmap=segm.cmap, interpolation = None)  #
        ax2.set_title( str(npixels_value) + '_' + str(th_coeff), fontsize=12)
        
        fig, ax = plt.subplots()
        im = ax.imshow(np.array(segm), origin="lower", cmap=segm.cmap, interpolation = None)  #
        plt.colorbar(im)

        #hdu[0].data = segm
        #hdu.writeto("/xdisk/hamden/hina0830/venv39/raw_files/2024-01-11-RA195-DEC21/" + "NGA_M51-fd-segm_color.fits", overwrite=True)
        
        final_mask = segm

        mask_only = np.zeros(segm.shape)
        
        for i in range(len(tbl)):
            print("flux: ", tbl['segment_flux'][i])
            print("kron aperture: ", tbl["kron_aperture"][i])
            a, b = tbl["kron_aperture"][i].a, tbl["kron_aperture"][i].b
            ap_size = 1.75
            

            if a > 500000 or b > 500000:
                pass
            
            else:
                aperture = EllipticalAperture(
                                tbl["kron_aperture"][i].positions,
                                a=ap_size* tbl["kron_aperture"][i].a,
                                b=ap_size* tbl["kron_aperture"][i].b,
                                theta=tbl["kron_aperture"][i].theta,
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
            
        ''' 
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
 

    #if counter == 0:
    #    final_mask = canvas
    
    elif segm is None:
        final_mask = canvas
 

    #final_mask[np.isnan(data_orig)] = 0    
    final_mask[np.isnan(convolved_data)] = 0    
    
    '''
        
    elif segm is None:
        final_mask = canvas
        print("No extended sourced detected.")
        
    final_mask[np.isnan(convolved_data)] = 0  

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
    
    #fig, ax = plt.subplots()
    #plt.imshow(segmented_bg, origin="lower")#, norm=colors.LogNorm(vmin=0.1, vmax=1300))
    #plt.title("segment_bg")
    #plt.colorbar()
    
    #fig, ax = plt.subplots()
    #plt.imshow(data, origin="lower")#, norm=colors.LogNorm(vmin=0.1, vmax=1300))
    #plt.title("skybg * rrhr")
    #plt.colorbar()
    
    # Feed background values from skybg and generate poisson distribution

    nanarry = np.where(np.isnan(segmented_bg))
    print("mask values; ", mask[nanarry[0], nanarry[1]])
    print("data values; ", data[nanarry[0], nanarry[1]])    
    
    # Drop nans by replacing with 0s
    segmented_bg[np.isnan(segmented_bg)] = 0
 
    noise_added = np.random.poisson(lam=segmented_bg)

    # Replace the pixels in cnt with the generated noise
    cnt_noise = np.where(segmented_bg != 0, noise_added, cnt)
    
    fig, ax = plt.subplots()
    plt.imshow(cnt_noise, origin="lower", norm=colors.LogNorm(vmin=1e-2, vmax=1))
    plt.title("final image, noise added")
    plt.colorbar()
    plt.close()

    # Return the final image
    return cnt_noise

#def starfinder(file_tuple):
    #cnt_file, flags_file = file_tuple
    
    #return mask_data
