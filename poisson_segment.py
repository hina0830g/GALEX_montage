from astropy.io import fits
from photutils.datasets import make_100gaussians_image
from photutils.background import Background2D, MedianBackground
from astropy.convolution import convolve
from photutils.segmentation import detect_sources, make_2dgaussian_kernel
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import simple_norm
from photutils.segmentation import deblend_sources
from photutils.segmentation import SourceFinder
from photutils.segmentation import SourceCatalog
from astropy.stats import sigma_clipped_stats

from photutils.aperture import EllipticalAperture
from astropy.coordinates import Angle

from matplotlib.patches import Ellipse
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture
import pickle
from astropy.table import QTable
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

def segmentation_ver2(fn_current, data_orig, npixels_value, th_coeff):
    """
    ### input(s)
    
    fn_current(string): name of the file for labeling purposes
    
    data_orig(2d array, float64): preprocessed (convolved) intensity file; a Guassian filter has been applied

    npixel_value(int64): number of pixels connected. Default = 10,000 for large objects

    th_coeff(float64): coefficient for determing the threshold of pixel value for detection.
    threshold is determined by this vlaue * standard diviation of the image

    ### output(s)
    bool_mask(0s and 1s): segmented mask. 0 = no mask, 1 = mask

    """ 
    
    convolved_data = data_orig # preprocessed intensity file is already convolved.

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
    canvas = np.zeros(shape=(len(convolved_data), len(convolved_data[0])))
    
    fig, ax = plt.subplots()
    plt.imshow(data_orig, origin="lower", norm=colors.LogNorm(vmin=1e-4, vmax=1e-2))
    plt.title(fn_current + '_' + str(npixels_value) + '_' + str(th_coeff)) 
    plt.colorbar()
    
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
                              width=12.5 * tbl['semimajor_sigma'][i].value,  
                              height=12.5 * tbl['semiminor_sigma'][i].value, 
                              angle=tbl['orientation'][i].value, edgecolor = 'red', facecolor= 'none')
            print("ellipse: ", ellipse)
            print("ellipse parameters: ", tbl['xcentroid'][i], tbl['ycentroid'][i], 12.5*tbl['semimajor_sigma'][i].value, 12.5*tbl['semiminor_sigma'][i].value, tbl['orientation'][i].value)
            
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
        
    fig, ax = plt.subplots()
    plt.imshow(mask, origin="lower")
    plt.title("generated mask") 
    plt.colorbar()
    plt.show()

    gc.collect()
    return mask



def segmentation(fn, data_orig, npixels_value=500, th_coeff=6):
    """
    ### input(s)
    data(2d array, float64): preprocessed int files

    npixel_value(int64): number of pixels connected. Default = 10,000 for large objects

    th_coeff(float64): coefficient for determing the threshold of pixel value for detection.
    threshold is determined by this vlaue * standard diviation of the image

    ### output(s)
    bool_mask(0s and 1s): segmented mask. 0 = no mask, 1 = mask

    """
    #hdu = fits.open(fn)
    fn_current = fn.removesuffix('-cnt_nan.fits')
    print(fn_current)
    
    convolved_data = data_orig

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
        columns = ["label", "xcentroid", "ycentroid", "segment_flux", "kron_flux", "kron_aperture"]
        tbl = cat.to_table(columns=columns)

        norm = simple_norm(data, "sqrt")
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 12.5))
        ax1.imshow(data_orig, origin="lower", norm=colors.LogNorm(vmin=1e-4, vmax=1e-2))
        ax1.set_title(fn_current, fontsize=14)
        ax2.imshow(np.array(segm), origin="lower", interpolation = None)  #
        ax2.set_title( str(npixels_value) + '_' + str(th_coeff), fontsize=12)
        
        #hdu[0].data = segm
        #hdu.writeto("/xdisk/hamden/hina0830/venv39/raw_files/2024-01-11-RA195-DEC21/" + "NGA_M51-fd-segm_color.fits", overwrite=True)
        
        final_mask = segm

        mask_only = np.zeros(segm.shape)
        
        for i in range(len(tbl)):
            a, b = tbl["kron_aperture"][i].a, tbl["kron_aperture"][i].b
            ap_size = 1.00
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
    
    
    # Drop nans by replacing with 0s
    segmented_bg[np.isnan(segmented_bg)] = 0
    # Exclude (unmask) pixels where the data has nans
    segmented_bg[np.isnan(data)] = 0
 
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

def psfinder(divided_data, flags_data, th_coeff, DAO_fwhm, mask_size):
    
    """
    ### input(s)
    divided_data(2d array, float64): flagged intensity image 
    
    flags_data(2d array, float64): wcs transformed flag image
    
    th_coeff(float64): coefficient for determing the threshold of pixel value for star finder/2D gaussian fit
    threshold is determined by this vlaue * standard diviation of the image
    
    DAO_fwhm(int64): the target value of fwhm (full width half maxima) for star finder/2D gaussian fit
    
    mask_size(int64): the radius of a mask for each point source

    ### output(s)
    bool_mask(0s and 1s): segmented mask. 0 = no mask, 1 = mask

    """
    # Determine the standard dev.
    mean, median, std = sigma_clipped_stats(divided_data)  
    print((mean, median, std)) 
    
    daofind = DAOStarFinder(fwhm=DAO_fwhm, threshold=th_coeff*std) # fwhm=3.0, 3.0*std 
    sources = daofind(divided_data - median)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=4.0)
    norm = ImageNormalize(stretch=SqrtStretch())
    
    # Define custom radii for the stars 
    star_radii = np.ones(len(sources)) * mask_size

    # Mask stars by replacing pixels with 0s
    mask_data = np.copy(divided_data)

    # create an epmty array, same size as the image
    canvas = np.zeros(shape=(divided_data.shape[0], divided_data.shape[1]))

    y_indices, x_indices = np.indices(divided_data.shape)
    x_cent, y_cent = int(round(len(divided_data) / 2)), int(round(len(divided_data) / 2))
    coord_lis = []
    
    for star, radius in zip(sources, star_radii):
        x, y = int(star['xcentroid']), int(star['ycentroid'])

        # Distances between the star and every pixel in the image
        dist = np.sqrt((x_indices - x) ** 2 + (y_indices - y) ** 2) 

        # Find the distance between the center of the image and the star
        star_dist = np.sqrt((x_cent - x) ** 2 + (y_cent - y) ** 2)

        # Ignore the stars near the edge (outside the inner radius)
        if star_dist <= 1400: 
            # fill in the pixel around the star and mask the star in the image and the mask
            # the size of masking for each star is mask_size defined at the top
            # masking the image; masked region = 0
            mask_data[dist <= radius] = 0

            # mask file; masked = 1 and nonmasked = 0
            canvas[dist <= radius] = 1

            coord = [x, y]
            coord_lis.append(coord)

    masked_data = mask_data
    mask_data = np.where(mask_data==0, 1, canvas)
    print(len(coord_lis), " stars detected")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 12.5))

    ax1.imshow(divided_data, origin="lower", norm=colors.LogNorm(vmin=1e-4, vmax=1e-2))
    ax1.set_title("Image Before", fontsize=16)
    ax2.imshow(divided_data, origin="lower", norm=colors.LogNorm(vmin=1e-4, vmax=1e-2))
    ax2.set_title("Image with apertures", fontsize=16)

    apertures.plot(ax=ax2, color='red', lw=1.5, alpha=1.)
    
    ######### Flags segmentation

    ### Parameters below are fixed to detect artifacts
    npixels_value_wcs = 25
    threshold_wcs = 50
    
    segm = detect_sources(
    flags_data, threshold=threshold_wcs, npixels=npixels_value_wcs, connectivity=8
    )
    
    finder = SourceFinder(npixels=npixels_value_wcs, progress_bar=False, deblend=False)
    segment_map = finder(flags_data, threshold=threshold_wcs)

    canvas = np.zeros(shape=(len(flags_data), len(flags_data[0])))

    if segm is not None:
        cat = SourceCatalog(flags_data, segment_map)

        # Extract positions(x/ycentroid, flux, and eliptical aperture information),
        columns = ["label", "xcentroid", "ycentroid", "semimajor_sigma", "semiminor_sigma"]

        tbl = cat.to_table(columns=columns)

        # Assuming tbl is your QTable, you can convert it to a DataFrame if needed
        df = tbl.to_pandas()

        # Find the index of the row with the largest value of "semimajor_sigma"
        indices_to_drop = df[df['semimajor_sigma'] > 50].index

        # Drop the row using the index
        df.drop(indices_to_drop, inplace=True)

        # Convert back to QTable if needed
        tbl_new = QTable.from_pandas(df)

        mask_shape = flags_data.shape
        mask = np.zeros(mask_shape, dtype=bool)

        # Loop over ellipses and draw them onto the mask
        for i in range(len(tbl_new)):
            x, y = int(tbl_new['xcentroid'][i]), int(tbl_new['ycentroid'][i])
            # Find the distance between the center of the image and the star
            star_dist = np.sqrt((x_cent - x) ** 2 + (y_cent - y) ** 2)

            if star_dist < 1400: # Save the coordinates if the coordinates are inside the circle
                coord = [int(tbl_new['xcentroid'][i]), int(tbl_new['ycentroid'][i])]
                coord_lis.append(coord)
    
    fig, ax = plt.subplots()
    plt.imshow(flags_data, origin="lower")
    plt.title("flag") 
    plt.colorbar()  
    
    positions = np.transpose((tbl_new['xcentroid'], tbl_new['ycentroid']))
    apertures = CircularAperture(positions, r=4.0)
    apertures.plot(ax=ax, color='red', lw=1.5, alpha=1.)
    print(len(tbl_new),  " flagged coordinates.")
    
    print(len(coord_lis), "coordinates in total stored.")
    
    return masked_data, mask_data, coord_lis
