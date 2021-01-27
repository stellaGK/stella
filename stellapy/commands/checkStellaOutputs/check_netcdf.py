#!/usr/bin/python3
 
#################################################################
#  STELLAPY.DATA.MAGNETICFIELD AND STELLAPY.PLOT.MAGNETICFIELD
#################################################################

# Load the modules
import sys, os, pathlib
import netCDF4 as nc4 
import numpy as np

# Load the personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.utils import get_filesInFolder
from stellapy.utils.decorators.verbose_wrapper import indent

#======================
# Read the NETCDF file
#======================

# Get the folder
netcdf_folder = pathlib.Path(os.getcwd())

# Load the file
if get_filesInFolder(netcdf_folder, end="out.nc") or get_filesInFolder(netcdf_folder, start="restart.nc")[0]:
    if get_filesInFolder(netcdf_folder, end="out.nc"):
        netcdf_paths = get_filesInFolder(netcdf_folder, end="out.nc")
    elif get_filesInFolder(netcdf_folder, start="restart.nc")[0]:
        netcdf_paths = get_filesInFolder(netcdf_folder, start="restart.nc") 
        netcdf_paths = netcdf_paths[:1]
        
    # Read the file
    for netcdf_path in netcdf_paths:
        netcdf_file = nc4.Dataset(netcdf_path)
        
        # Print the dimensions
        print()
        print(indent, '######################################################')
        print(indent, '  Read netcdf file: ', netcdf_path)
        print(indent, '######################################################', '\n')
        print(indent, '======================')
        print(indent, '  NETCDF DIMENSIONS')
        print(indent, '======================', '\n')
        print(indent, ' ', '{0:16}'.format("Dimensions"), '{0:16}'.format("  Shape"), '{0:35}'.format("     Value"), '\n')
        for name, value in netcdf_file.dimensions.items():
            variable = name
            shape = str(value).split('size =')[1].split('\n')[0]
            try:   
                value = list(netcdf_file.variables[variable])
            except: 
                value = [np.nan]
        
            if np.array(value).size > 1:   
                value = np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[:30] + ' ... ' + \
                        np.array2string(np.array(value), formatter={'float_kind':lambda x: "%.2f" % x})[-30:]  
                value = value.replace("\n", "") 
            else:  
                try:                          
                    value = str(value[0])
                except: 
                    value = 'nan'
            if value == 'nan':               
                value = " /  "
        
            # Print the dimensions
            print( indent, '    ', '{0:16}'.format(variable+':'), '{0:16}'.format(shape), '{0:50}'.format(value))
        
        # Print the netcdf variables
        print()
        print(indent, '======================')
        print(indent, '   NETCDF VARIABLES')
        print(indent, '======================', '\n')
        print(indent, '    ', '{0:16}'.format("Variable"), '{0:20}'.format("Shape"), '{0:10}'.format("      Size"), '{0:30}'.format("     Function"), '{0:60}'.format("Long name"), '\n')
        for name, value in netcdf_file.variables.items():
            variable            = name
            function            = variable + str(value).split(" "+variable)[1].split('\n')[0].replace(" ", "")     
            shape               = str(value).split("shape = ")[1].split('\n')[0].replace(" ", "") 
            try:    long_name   = str(value).split("long_name: ")[1].split('\n')[0]
            except: long_name   = ""
            shapes = shape.split("(")[-1].split(")")[0].split(",")
            shapes = [ int(s) for s in shapes if s!=""]; size=1
            for i in range(len(shapes)):
                size = size*shapes[i]
            
            # Print the netcdf variables
            print( indent, '    ', '{0:16}'.format(variable+':'), '{0:20}'.format(shape), '{0:10}'.format(size), '{0:30}'.format("   "+function), '{0:60}'.format(long_name))
else:
    print("No 'out.nc' file was found.")