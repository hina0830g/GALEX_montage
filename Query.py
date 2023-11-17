"""
Created on Thu May 12 14:23:17 2022

@author: Hina Goto
"""
import numpy as np
import matplotlib.pyplot as py
from IPython.display import Image
from astroquery.mast import Catalogs
from astroquery.mast import Observations
import pandas as pd
from astropy.table import QTable
from astropy import units as u
from astropy.coordinates import SkyCoord
import os
import shutil
import sys
import glob
import gzip
from sh import gunzip

# Description: This program runs coordinate query and downloads cnt, rrhr,
# skybg, objmask files needed to run MontagePy on. 
# To run the program, specify the path to run (=path_init) below and 
# coordinate of your target (=c)

# Specify the initial path to start the program; a working directory will be made in this location
path_init = '/mnt/c/Users/Amiti/Documents/Spring2022/GALEX/raw_files/Test_py'
os.chdir( path_init )
os.getcwd()

# Specify the coordinate below; it is set to be RA 195 degree and DEC 21 degree 
c = SkyCoord(ra=195*u.degree, dec=21*u.degree, frame='icrs')
r = float(input('What search radius? '))

# Runs GALEX coordinate query 
result_table = Observations.query_criteria(coordinates=c, obs_collection='GALEX', radius=r, filters=["*FUV"])
obsid = result_table ['obsid']
print(len(obsid), "obsids available")

data_products = Observations.get_product_list(obsid)
data_products =  data_products.to_pandas()

# Only keeps the fd and AIS files;
data_products = data_products[(data_products['productFilename'].str.startswith('AIS'))]
data_products = data_products[(data_products['productFilename'].str.contains('-fd'))]

# Chooses files of the types we want (cnt, rrhr, objmask, and skybg)
DP_cnt = data_products[(data_products['productFilename'].str.endswith('-cnt.fits.gz'))]

print(len(DP_cnt), 'files per type (cnt, rrhr, objmask, skybg) available to download')
answer = input("Would you like to download the files? Type yes or no: ")
if answer.lower().startswith("y"):
    print('\n', 'Ok, downloading', len(DP_cnt), 'cnt, rrhr, skybg, and objmask files.')
    
    dir_name = str(input('Create a dirctory download the files in; give the directory a name: '))
    os.mkdir(dir_name)
    path_new = path_init + '/' + dir_name
    os.chdir(path_new)
    print('\n')
    print('Downloading files to the following directory: ', os.getcwd())
     
    answer3 = str(input('proceed? yes or no: ')) 
    if answer3.lower().startswith("y"):
        data_products = Observations.get_product_list(obsid)
        print(type(data_products))
        data_products =  data_products.to_pandas() # Converts data products to panda data frame to execute suffix search
    
        # Only keeps the fd and AIS files;
        data_products = data_products[(data_products['productFilename'].str.startswith('AIS'))]
        data_products = data_products[(data_products['productFilename'].str.contains('-fd'))]
    
        # Finds cnt, rrhr, and skybg files and stores as DP_cnt, DP_rrhr, and DP_skybg (=suffix search)
        DP_cnt = data_products[(data_products['productFilename'].str.endswith('-cnt.fits.gz'))] #| (data_products['productFilename'].str.endswith('-rrhr.fits.gz')) | (data_products['productFilename'].str.endswith('-skybg.fits.gz'))]
        DP_rrhr = data_products[(data_products['productFilename'].str.endswith('-rrhr.fits.gz'))] #| (data_products['productFilename'].str.endswith('-rrhr.fits.gz')) | (data_products['productFilename'].str.endswith('-skybg.fits.gz'))]
        DP_skybg = data_products[(data_products['productFilename'].str.endswith('-skybg.fits.gz'))] #| (data_products['productFilename'].str.endswith('-rrhr.fits.gz')) | (data_products['productFilename'].str.endswith('-skybg.fits.gz'))]
        DP_mask = data_products[(data_products['productFilename'].str.endswith('-objmask.fits.gz'))] #| (data_products['productFilename'].str.endswith('-rrhr.fits.gz')) | (data_products['productFilename'].str.endswith('-skybg.fits.gz'))]
    
        # Convert data frame back to QTable
        DP_cnt = QTable.from_pandas(DP_cnt)
        DP_rrhr = QTable.from_pandas(DP_rrhr)
        DP_skybg = QTable.from_pandas(DP_skybg)
        DP_mask = QTable.from_pandas(DP_mask)
        print(type(data_products))
    
        # Download selected files
        manifest = Observations.download_products(DP_cnt)
        manifest = Observations.download_products(DP_rrhr)
        manifest = Observations.download_products(DP_skybg)
        manifest = Observations.download_products(DP_mask)
      
        print(manifest)
        print('installation completed: ')
      
    elif answer3.lower().startswith("n"):
        print("ok, finishing the program. Start over with a new search radius")
        sys.exit()
     
     
elif answer.lower().startswith("n"):
   print("ok, finishing the program. Start over with a new search radius")
   #sys.exit()


# ----- Below organizes the downloaded files that ----- #
# ----- are in several foulders and put them all  ----- #
# ----- in one directory                          ----- #   
# ------------------------------------------------------#
# ----- In the end, the program removes the empty ----- #
# ----- directories /mastDownload/GALEX           ----- #

ext = 'mastDownload/GALEX'
os.chdir( path_new + '/' + ext)
print(os.getcwd())

current_folder = os.getcwd() 
list_dir = os.listdir(current_folder) 
  
# enumerate on list_dir to get the content of all the folders 
# and store it in a dictionary
content_list = {}
for index, val in enumerate(list_dir):
    path = os.path.join(current_folder, val)
    content_list[ list_dir[index] ] = os.listdir(path)
print(content_list)

# Distination of the files; default is the newly created directory 
merge_folder_path = path_new 
print('\n')
print('Transferring downloaded files to the following directory: ', merge_folder_path)
answer4 = str(input('Please confirm (yes or no): '))

if answer4.lower().startswith("y"):
    # loop through the folders under /mastDownload/GALEX
    for sub_dir in content_list:
    
        # loop through the contents of the
        # list of folders
        for contents in content_list[sub_dir]:
    
            # make the path of the content to move 
            path_to_content = sub_dir + "/" + contents  
    
            # make the path with the current folder
            dir_to_move = os.path.join(current_folder, path_to_content )
  
            # move the file
            shutil.move(dir_to_move, path_new)
        os.rmdir(sub_dir)

elif answer4.lower().startswith("n"):
    print("ok, finishing the program without transferring the files")

# Removes empty folders GALEX and mastDownload
os.chdir( path_new + '/mastDownload' )
os.rmdir('GALEX')
os.chdir( path_new )
print(os.getcwd())
os.rmdir('mastDownload') 

merge_folder_path = path_new
sourcepath = merge_folder_path
sourcefiles = os.listdir(sourcepath)
destinationpath = sourcepath # Moves all the files to the raw directory

answer5 = str(input('Gunzip all the download files; please confirm (yes or no): '))
if answer5.lower().startswith("y"):
    for file in sourcefiles:
        if file.endswith('.gz'):
            shutil.move(os.path.join(sourcepath,file), os.path.join(destinationpath,file))
    os.chdir(destinationpath)
    for file in glob.glob("*.gz"):
        try:
            gunzip(file)
        except:
            print('already exists')
    
elif answer4.lower().startswith("n"):
    print("ok, finishing the program without gunzipping the files")
    sys.exit()
