import os, time
from datetime import date
import argparse
import ast
import shutil
from IPython.display import Image
from astroquery.mast import Catalogs
from astroquery.mast import Observations
import pandas as pd
from astropy.table import QTable
from astropy import units as u
from astropy.coordinates import SkyCoord
import glob
from sh import gunzip

# Input here:
r = 5
project_lis = ['AIS', 'GII', 'CAI', 'NGS', 'DIS', 'CAI', 'MIS']

# Swicth to the raw file directory
path_init = os.getcwd() + '/raw_files'
os.chdir( path_init )
print('initial path (path_init): ', path_init)

today = date.today()
print("Today's date:", today)

# Create an ArgumentParser object
parser = argparse.ArgumentParser()

# Add an argument for the list of directory and filename pairs
# nargs='+' indicates that the argument can take one or more values
parser.add_argument('--pair-list', nargs=2, type=str, help='Pair of integers (int1 int2)')

# Parse the command-line arguments
args = parser.parse_args()

# Use the list of directory and filename pairs in your script
if args.pair_list:
    # Split a pair into a directory and a filename
    RA, DEC = map(int, args.pair_list)
    print("RA:", RA)
    print("DEC:", DEC)
    try:
        # Specify the coordinates using the input information
        c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='icrs')
        dir_name = str(today) + "-RA" + str(RA) + "-DEC" + str(DEC)
        
        # Create a new directory under raw_files
        os.makedirs(dir_name, exist_ok=True)
        
        # Switch to the new directory where files
        # will be downloaded
        path_new = path_init + '/' + dir_name
        os.chdir( path_new )
        
        print("directory name:", path_new)
    except FileNotFoundError:
        print("Error")
 
project_lis = ['AIS', 'GII', 'CAI', 'NGS', 'DIS', 'CAI', 'MIS']

def coord_query(coord, radius, project):
    # Runs GALEX coordinate query 
    result_table = Observations.query_criteria(coordinates=c, obs_collection='GALEX', radius=r*u.degree, project=project_lis, filters=["*FUV"])
    obsid = result_table['obsid']
    print(len(obsid), "obsids available")
    data_products = Observations.get_product_list(obsid)
    data_products =  data_products.to_pandas()
    
    # Only keeps the fd files;
    data_products = data_products[(data_products['productFilename'].str.contains('-fd'))]
    DataProduct_final = data_products[
    ((data_products['productFilename'].str.contains('_000' + r'\d{1}' + '-' )) == False) &
    ((data_products['productFilename'].str.contains('_000' + r'\d{1}' + '_' )) == False) &
    (data_products['productFilename'].str.endswith(('-cnt.fits.gz', '-rrhr.fits.gz', '-skybg.fits.gz', '-objmask.fits.gz', 'fd-ncat.fits.gz')))
    ]    

    print(len(DataProduct_final), 'files total per type available to download')
    print(DataProduct_final['productFilename'])    

    return DataProduct_final

def download_dp(data_products):

    # Convert data frame back to QTable
    data_products = QTable.from_pandas(data_products)
    # Download selected MIS files
    manifest = Observations.download_products(data_products)
    print("Installed: " , data_products['productFilename'])
    return path_new

def download_dp(data_products):

    # Convert data frame back to QTable
    data_products = QTable.from_pandas(data_products)
    # Download selected MIS files
    manifest = Observations.download_products(data_products)
    print("Installed: " , data_products['productFilename'])
    return path_new


def move_files(path_new):
    ext = '/mastDownload/GALEX'
    os.chdir( path_new + ext)
    print(os.getcwd())
    
    current_folder = os.getcwd() 
    list_dir = os.listdir(current_folder)  # All the subdirectories
    
    # enumerate on list_dir to get the content of all the folders 
    # and store it in a dictionary
    content_dict = {}
    for index, val in enumerate(list_dir):
        path = os.path.join(current_folder, val)
        # os.listdir lists all everything (files/dir) in the specified directory
        content_dict [ list_dir[index] ] = os.listdir(path) # Assigns the files as keys in the dictionary 
    print(content_dict)
    
    # Distination of the files; default is the newly created directory 
    merge_folder_path = path_new 
    
    # loop through the folders under /mastDownload/GALEX
    for sub_dir in content_dict:
    
        # loop through the contents inside
        for contents in content_dict[sub_dir]:
    
            # make the path of the content to move 
            file_path = path_new + ext + "/" + sub_dir + "/" + contents
    
            # move the file
            shutil.move(file_path, path_new)

        os.rmdir(sub_dir) # remove the empty subdirectories
    
    # Removes empty folders: GALEX and mastDownload
    os.chdir( path_new + '/mastDownload' )
    os.rmdir('GALEX')
    os.chdir( path_new )
    print(os.getcwd())
    os.rmdir('mastDownload') 
    
    merge_folder_path = path_new
    sourcepath = merge_folder_path
    sourcefiles = os.listdir(sourcepath)
    destinationpath = sourcepath # Moves all the files to the raw directory
    
    for file in sourcefiles:
        if file.endswith('.gz'):
            shutil.move(os.path.join(sourcepath,file), os.path.join(destinationpath,file))
    os.chdir(destinationpath)
    
    for file in glob.glob("*.gz"):
        try:
            gunzip(file)
        except:
            print('already exists')
            
    return file


if __name__=="__main__":
    # Run query & extract data product
    DataProduct_final = coord_query(c, r, project_lis)

    path_new = download_dp(DataProduct_final)
        
    # Organize the files
    output = move_files( path_new )

