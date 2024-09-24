#!/usr/bin/python3
 
#################################################################
#  STELLAPY.DATA.MAGNETICFIELD AND STELLAPY.PLOT.MAGNETICFIELD
#################################################################

# Load the modules
import sys, os, pathlib, time, copy, h5py
import netCDF4 as nc4 
import numpy as np 

# Load the personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)
from stellapy.utils import get_filesInFolder
from stellapy.utils.decorators.verbose import indent

#======================
# Read the NETCDF file
#======================

# Print the size of the vectors
print_size = False

# Get the folder
netcdf_folder = pathlib.Path(os.getcwd())
h5_files = get_filesInFolder(netcdf_folder, end=".out.h5")[0]

# Open the "out.h5" file
try: h5_netcdf = h5py.File(h5_files, 'r')
except: print("Could not read the h5 file: \n   ", h5_files); sys.exit(0)

# Safe the data to a dictionary
print(list(h5_netcdf.attrs.keys()))
print(list(h5_netcdf.keys()))
