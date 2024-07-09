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

# Scaling factor for EllipticalAperture. 2 * kron radius should cover > 90% of the flux (Kron et al 1980).
ap_size = 2.00

def segmentation(fn, data_orig, npixels_value=500, th_coeff=6):
    
    """
    Detects large sources ( extended sources and galaxies ) in an image and 
    determiens the regions to be be infilled.
    
    Parameters
    ----------
    fn : str
         filename of the image. Only used for labeling. 
        
    data_orig : numpy array (float 64)
        The image to perform segmentation on. It should be 
        intesity file that has been preprocessed (convolved).
    
    npixels_value : int
        The minimum number of connected pixels in order to be
        detedcted as a segment.
    
    th_coeff : float
        Used to determine the minimum value of pixels in order to be detected as a segment.
        The threshold is th_coeff * std of the image

    Returns
    -------
    bool_mask : numpy array ( bool )
        The boolean mask generated via segmentation. 
        Indicates the regions where large sources are detected.
        0 = not masked, 1 = masked

    Description
    -----------
    This functions detects large sources and creates a mask for
    detected regions. The detection process is through segmentation, 
    which segments regions that have a minimum number of connected pixels
    above our threhsold value.
    
    """  
    
    fn_current = fn.removesuffix('-cnt_nan.fits')
    print(fn_current)
    
    convolved_data = data_orig # preprocessed, thus already convolved.

    mean, _, std = sigma_clipped_stats(convolved_data) # use sigma-clipped statistics to (roughly) estimate the background noise levels
    data = convolved_data - mean # subtract the background

    threshold = th_coeff * std  # Defines the threhsold for segmentation

    # Performs segmentation 
    finder = SourceFinder(npixels=npixels_value, progress_bar=False, deblend=False)
    segment_map = finder(convolved_data, threshold=threshold)
    print("detect_sources done.")

    canvas = np.zeros(shape=(len(data), len(data[0])))
    if segment_map is not None:
        cat = SourceCatalog(data_orig, segment_map, convolved_data=convolved_data)
        
        # Extract positions(x/ycentroid, flux, and eliptical aperture information),
        columns = ["label", "xcentroid", "ycentroid", "segment_flux", "kron_flux", "kron_aperture"]
        tbl = cat.to_table(columns=columns)

        norm = simple_norm(data, "sqrt")
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 12.5))
        ax1.imshow(data_orig, origin="lower", norm=colors.LogNorm(vmin=1e-4, vmax=1e-2))
        ax1.set_title(fn_current, fontsize=14)
        ax2.imshow(np.array(segment_map), origin="lower", interpolation = None)  
        ax2.set_title( str(npixels_value) + '_' + str(th_coeff), fontsize=12)

        final_mask = canvas 
        
        for i in range(len(tbl)):
            
            # Empty canvas for individual sources
            canvas_source = np.zeros(shape=(len(data), len(data[0])))

            a, b = tbl["kron_aperture"][i].a, tbl["kron_aperture"][i].b
            
            aperture = EllipticalAperture(
                            tbl["kron_aperture"][i].positions,
                            a=ap_size* tbl["kron_aperture"][i].a,
                            b=ap_size* tbl["kron_aperture"][i].b,
                            theta=tbl["kron_aperture"][i].theta,
                        )

            mask = aperture.to_mask()
            ap_patches = aperture.plot(color="white", lw=1, ax=ax1)
            ap_patches = aperture.plot(color="white", lw=1, ax=ax2)


            x, y = (
                tbl["kron_aperture"][i].positions[1],
                tbl["kron_aperture"][i].positions[0],
            )
            width, height = np.shape(mask)[0], np.shape(mask)[1]

            top_left_x = int(x - (width / 2))
            top_left_y = int(y - (height / 2))
            bottom_right_x = top_left_x + width
            bottom_right_y = top_left_y + height
            
            try:
                canvas_source[top_left_x:bottom_right_x, top_left_y:bottom_right_y] = mask
                canvas += canvas_source

            except ValueError:
                print("--- !!! VALUE ERROR WARNING !!! ---")
                print(top_left_x, top_left_y, bottom_right_x, bottom_right_y)

                top_left_x_new, top_left_y_new, bottom_right_x_new, bottom_right_y_new = max(0, min(top_left_x, 2840)), max(0, min(top_left_y, 2840)), max(0, min(bottom_right_x, 2840)), max(0, min(bottom_right_y, 2840))
                print("new values: ", top_left_x_new, top_left_y_new, bottom_right_x_new, bottom_right_y_new)

                top_left_x_diff = top_left_x - top_left_x_new
                top_left_y_diff = top_left_y - top_left_y_new
                bottom_right_x_diff = bottom_right_x - bottom_right_x_new
                bottom_right_y_diff = bottom_right_y - bottom_right_y_new

                trimmed_mask  = mask.data
                
                if top_left_x_diff != 0: 
                    trimmed_mask = trimmed_mask[abs(top_left_x_diff):, :] # Bottom cut

                if top_left_y_diff != 0:
                    trimmed_mask = trimmed_mask[:, abs(top_left_y_diff):] # Left cut 

                if bottom_right_x_diff != 0:
                    trimmed_mask = trimmed_mask[:-abs(bottom_right_x_diff), :] # Top cut

                if bottom_right_y_diff != 0:
                    trimmed_mask = trimmed_mask[:, :-abs(bottom_right_y_diff)] # Right cut

                canvas_source[top_left_x_new:bottom_right_x_new, top_left_y_new:bottom_right_y_new] = trimmed_mask
                canvas += canvas_source

        final_mask = canvas
        
    elif segment_map is None:
        final_mask = canvas
        print("No extended sourced detected.")

    bool_mask = np.where(final_mask != 0, 1, final_mask) # make the mask 0s and 1s
    
    fig, ax = plt.subplots()
    plt.imshow(bool_mask, origin="lower")
    plt.title("final mask, aperture + segment")
    plt.clim(0, 1)
    plt.colorbar()

    return bool_mask


def poisson_noise(cnt, rrhr, mask, skybg):

    """
    Poisson infills the segmented regions and 
    removes the large sources.
    
    Parameters
    ----------
    cnt : numpy array (float 64)
        The preprocessed cnt file 
        
    rrhr : numpy array (float 64)
        The preprocessed rrhr file 
        
    mask : numpy array (float 64)
        The mask generated from segmentation
        
    skybg : numpy array (float 64)
        The preprocessed skybg file 

    Returns
    -------
    bool_mask : numpy array ( bool )
        The boolean mask generated via segmentation. 
        Indicates the regions where large sources are detected.
        0 = not masked, 1 = masked

    Description
    -----------
    This functions generates Poisson noise from the background 
    for the masked regins ( where large soures are ).
    """  
 
    # Define the sampple data to draw Poisson distribution 
    data = skybg * rrhr  # background * rrhr so it's in counts
    
    # Clean the sample data. Set unmaksed regions and nans to be 0
    segmented_bg = np.where(mask == 0, 0, data)  
    segmented_bg[np.isnan(segmented_bg)] = 0

    segmented_bg[np.isnan(data)] = 0
 
    # Draw Poisson distrubtion here
    noise_added = np.random.poisson(lam=segmented_bg)

    # Now fill in the large souces with drawn data ( Poisson noise )
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
    Detect the point sources ( stars ) in the image via DAOStarFinder
    and artifacts in flags via segmentation.
    
    Parameters
    ----------
    divided_data : numpy array (float 64)
        Poinsson infilled cnt data divided by preprocessed rrhr.
        
    flags_data : numpy array (float 64)
        The preprocessed skybg file.
        
    th_coeff : int
        This value * std of the image is the threshold 
        for source detection.
        
    DAO_fwhm : int
        The full-width half maximum of the major axis of
        the Gaussian kernel in units of pixels.
        
    mask_size : int
        The radius ( in pixel ) of individual mask for each point source
        The default = 16 pix. Determined from the instrument sigma 
        and the kernel size from Gaussian filter

    Returns
    -------
    masked_data : numpy array ( float 64 )
        Boolean mask generated via segmentation. 
        Indicates the regions where large sources are detected.
        0 = not masked, 1 = masked
        
    mask_data : numpy array ( bool )
        Boolean mask where 0 = not masked, 
        1 = masked
        
    coord_lis : list
        list of coordinates (x, y) of the detected sources.

    Description
    -----------
    1. This functions runs DAOStarFinder on the image to detect 
    and mask point sources.
    
    2. This function runs segmentation on the flags file to detect
    approximate centers of flagged regions and also flags the image.
    
    3. This functions stores the coordinates of point sources from 1 
    and flagged regions from 2 in a list ( coord_lis ). 
    
    """      
    
    # Determine the standard dev.
    mean, median, std = sigma_clipped_stats(divided_data)  
    
    daofind = DAOStarFinder(fwhm=DAO_fwhm, threshold=th_coeff*std) 
    sources = daofind(divided_data - median)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=4.0)
    norm = ImageNormalize(stretch=SqrtStretch())
    
    # Define custom radii for the stars 
    star_radii = np.ones(len(sources)) * mask_size

    mask_data = np.copy(divided_data)

    # create an epmty array, same size as the image; 
    canvas = np.zeros(shape=(divided_data.shape[0], divided_data.shape[1]))

    y_indices, x_indices = np.indices(divided_data.shape)
    x_cent, y_cent = int(round(len(divided_data) / 2)), int(round(len(divided_data) / 2))
    coord_lis = []
    
    # Loop over detected sources and calculate the distance from the center
    for star, radius in zip(sources, star_radii):
        x, y = int(star['xcentroid']), int(star['ycentroid'])
        star_dist = np.sqrt((x_cent - x) ** 2 + (y_cent - y) ** 2)

         # Store the coordinate and mask the point source 
         # if the distance is less than 1400 ( the boundary of the circle )
        if star_dist <= 1400: 
            mask_data[dist <= radius] = 0

            # mask the star; masked = 1 and nonmasked = 0
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

    finder = SourceFinder(npixels=npixels_value_wcs, progress_bar=False, deblend=False)
    segment_map = finder(flags_data, threshold=threshold_wcs)

    canvas = np.zeros(shape=(len(flags_data), len(flags_data[0])))

    if segment_map is not None:
        cat = SourceCatalog(flags_data, segment_map)

        # Specify columns to extract 
        columns = ["label", "xcentroid", "ycentroid", "semimajor_sigma", "semiminor_sigma"]

        tbl = cat.to_table(columns=columns)
        df = tbl.to_pandas()

        # Drop the sources with semimajor axis larger than 50 since we don't want to infill the "dents"
        indices_to_drop = df[df['semimajor_sigma'] > 50].index
        df.drop(indices_to_drop, inplace=True)

        tbl_new = QTable.from_pandas(df)

        mask_shape = flags_data.shape
        mask = np.zeros(mask_shape, dtype=bool)

        # Loop over detected sources and calculate the distance from the center
        for i in range(len(tbl_new)):
            x, y = int(tbl_new['xcentroid'][i]), int(tbl_new['ycentroid'][i])
            star_dist = np.sqrt((x_cent - x) ** 2 + (y_cent - y) ** 2)
            
            if star_dist < 1400: # Store the coordinate if the distance is less than 1400 ( the boundary of the circle )
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
