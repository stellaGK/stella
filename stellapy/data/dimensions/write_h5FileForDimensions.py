""" 

#===============================================================================
#                     Write *.dimensions based on *.out.nc                     #
#=============================================================================== 

Write the (t,s,z,mu,vpa,kx,ky) dimensions to an h5 file. If multiple linear 
simulations have been launched separately, then combine their (kx,ky) dimensions
into a dummy dimensions file, e.g. input_v1_dummy.dimensions, to save memory. 
     >> write_stellapyDataFiles -s dim 

Hanne Thienpondt
20/01/2023

"""

#!/usr/bin/python3
import h5py
import os, sys
import pathlib
import numpy as np
from datetime import datetime, timedelta 
from stellapy.data.input.read_inputFile import read_modeFromInputFile

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.input.write_listOfMatchingInputFiles import read_inputFilesInDummyInputFiles
from stellapy.data.input.write_listOfMatchingInputFiles import read_inputFilesInDummyInputFile
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables
from stellapy.data.utils.print_status import print_status, get_nestingDepth
from stellapy.utils.files.get_filesInFolder import get_filesInFolder     

#===============================================================================
#                     Write *.dimensions based on *.out.nc                     #
#=============================================================================== 
 
def write_h5FileForDimensions(folder):      

    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))] 
    nesting_depth = get_nestingDepth(input_files)
    
    # For linear flux tube simulations, read the dummy input file instead 
    dummy_input_files = get_filesInFolder(folder, end="_dummy.in") 
    
    # Remove the input files listed in the dummy inputs from the list of input files
    input_files_in_dummy_files = read_inputFilesInDummyInputFiles(dummy_input_files)
    input_files = [i for i in input_files if i not in input_files_in_dummy_files]
    input_files += dummy_input_files
    if input_files==[]: return  
    
    # Initiate the dimensions
    dimensions = {}
 
    # Iterate through the input files   
    for input_file in input_files: 
            
        # Path of the new dimensions file and the netcdf file 
        dimensions_path = input_file.with_suffix(".dimensions")  
        netcdf_path = input_file.with_suffix(".out.nc")
        
        # For a dummy in put file, grab one of the netcdf files
        if "_dummy.in" in str(input_file): 
            input_files_in_dummy_file = read_inputFilesInDummyInputFile(input_file) 
            input_files_in_dummy_file_with_outnc = [i for i in input_files_in_dummy_file if os.path.isfile(i.with_suffix('.out.nc'))]
            if len(input_files_in_dummy_file_with_outnc)==0: continue
            netcdf_path = input_files_in_dummy_file_with_outnc[0].with_suffix(".out.nc")
        
        # Check whether the txt file is older than the simulation 
        if os.path.isfile(dimensions_path):  
            if datetime.fromtimestamp(netcdf_path.stat().st_mtime)<datetime.fromtimestamp(dimensions_path.stat().st_mtime)+timedelta(seconds=5):
                if "_dummy.in" not in str(input_file): 
                    continue 
                if "_dummy.in" in str(input_file): 
                    for i in read_inputFilesInDummyInputFile(input_file):
                        if datetime.fromtimestamp(i.with_suffix('.out').stat().st_mtime)>datetime.fromtimestamp(dimensions_path.stat().st_mtime)+timedelta(seconds=5):
                            break 
                    else: 
                        continue 
                    
        # Read the dimensions from the output file
        netcdf_path = read_outputFile(netcdf_path)  
        dimensions['vec_z']   = read_netcdfVariables('vec_z',   netcdf_path)
        dimensions['vec_kx']  = read_netcdfVariables('vec_kx',  netcdf_path)
        dimensions['vec_ky']  = read_netcdfVariables('vec_ky',  netcdf_path)
        dimensions['vec_mu']  = read_netcdfVariables('vec_mu',  netcdf_path)
        dimensions['vec_vpa'] = read_netcdfVariables('vec_vpa', netcdf_path)
        dimensions['vec_time'] = read_netcdfVariables('vec_time', netcdf_path)
        dimensions['vec_species'] = read_netcdfVariables('vec_species',  netcdf_path) 
        netcdf_path.close()
        
        # Get the size of the dimensions
        dimensions['dim_z']   = len(dimensions['vec_z'])
        dimensions['dim_kx']  = len(dimensions['vec_kx'])
        dimensions['dim_ky']  = len(dimensions['vec_ky'])
        dimensions['dim_mu']  = len(dimensions['vec_mu'])
        dimensions['dim_vpa'] = len(dimensions['vec_vpa'])
        dimensions['dim_time'] = len(dimensions['vec_time'])
        dimensions['dim_species'] = len(dimensions['vec_species'])
        
        # For the dummy input file, create vec_kx and vec_ky
        if "_dummy.in" in str(input_file): 
            dimensions['vec_kx'] = []; dimensions['vec_ky'] = []; dim_time = 0
            for input_file_in_dummy_file in input_files_in_dummy_file:  
                if os.path.isfile(input_file_in_dummy_file.with_suffix(".out.nc")):
                    netcdf_path = read_outputFile(input_file_in_dummy_file.with_suffix(".out.nc"))  
                    dimensions['vec_kx'] += list(read_netcdfVariables('vec_kx', netcdf_path))
                    dimensions['vec_ky'] += list(read_netcdfVariables('vec_ky', netcdf_path))
                    dim_time = np.max([dim_time, len(read_netcdfVariables('vec_time', netcdf_path))])
                    netcdf_path.close() 
                if not os.path.isfile(input_file_in_dummy_file.with_suffix(".out.nc")):
                    kx, ky = read_modeFromInputFile(input_file_in_dummy_file)
                    dimensions['vec_kx'] += [kx]
                    dimensions['vec_ky'] += [ky]
            dimensions['vec_kx'] = np.array(sorted(list(set(dimensions['vec_kx']))))
            dimensions['vec_ky'] = np.array(sorted(list(set(dimensions['vec_ky']))))
            dimensions['dim_kx'] = len(dimensions['vec_kx'])
            dimensions['dim_ky'] = len(dimensions['vec_ky'])
            dimensions['dim_time'] = dim_time 
        
        # Create the new dimensions h5 file 
        with h5py.File(dimensions_path, 'w') as h5_file:
            for key, value in dimensions.items():    
                h5_file.create_dataset(key, data=value) 
        print_status(dimensions_path, "dimensions", input_file, input_files, nesting_depth)  
           
    return 
