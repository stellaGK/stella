
import os, sys 
import h5py
import pathlib
import numpy as np
from datetime import datetime, timedelta
from stellapy.data.input.read_inputFile import read_inputFile 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder 
from stellapy.data.input.write_listOfMatchingInputFiles import read_inputFilesInDummyInputFiles
from stellapy.data.geometry.read_geometry import read_geometryFromTxtFile, read_geometryFromNetcdfFile
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.data.geometry.read_wout import read_woutFromNcFile, read_millerParameters
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_gridDivisionsAndSize 
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_muWeights, calculate_vpaWeights
from stellapy.data.geometry.read_output import calculate_dlOverB
from stellapy.utils.decorators.exit_program import exit_program
 
#===============================================================================
#                          WRITE THE GEOMETRY FILE
#===============================================================================
 
def write_h5FileForGeometry(folder): 
    ''' Get the geometry data from the netcdf file, the wout file and the geometry
    txt file and combine the data in a single '.geo' file. '''    

    # Get the input_files inside <folder> that have a *.out.nc file
    input_files = get_filesInFolder(folder, end='.in') 
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]

    # For linear flux tube simulations, read the dummy input file instead 
    dummy_input_files = get_filesInFolder(folder, end='_dummy.in') 
    
    # Remove the input files in the dummy inputs from the list of input files 
    input_files_in_dummy_files = read_inputFilesInDummyInputFiles(dummy_input_files)
    input_files = [i for i in input_files if i not in input_files_in_dummy_files]
    input_files += dummy_input_files
    
    # Collect the parent directories (Ignore OLD directories)
    parent_directories = list(set([i.parent if 'OLD' not in i.parent.name else i.parent.parent for i in input_files]))
    
    # Only look at the input files in each parent folder
    for directory in parent_directories:     
        
        # Get the files in <directory>
        input_files_directory = [i for i in input_files if (i.parent==directory) or (i.parent.parent==directory and 'OLD' in i.parent.name)]
        
        # Iterate over the input files
        for input_file in input_files_directory:
            
            # Geometry file name
            geometry_path = input_file.with_suffix('.geo')
            
            # Only write the *.geometry file if it doesn't exist or it's older than the *.out.nc file                                                         
            if os.path.isfile(geometry_path):
                if '_dummy' not in str(input_file):
                    if datetime.fromtimestamp(input_file.stat().st_mtime) < (datetime.fromtimestamp(geometry_path.stat().st_mtime)+timedelta(minutes=5)):
                        continue
                elif '_dummy' in str(input_file):
                    from stellapy.data.input.write_listOfMatchingInputFiles import read_inputFilesInDummyInputFile
                    for i in read_inputFilesInDummyInputFile(input_file):
                        if datetime.fromtimestamp(i.stat().st_mtime) > (datetime.fromtimestamp(geometry_path.stat().st_mtime)+timedelta(minutes=5)):
                            break
                    else:
                        continue
            
            # Write the *.geo file
            write_geometryFile(input_file, geometry_path)   
            
    return

#===============================================================================
#                           SAVE THE GEOMETRY OBJECT                           #
#=============================================================================== 

def write_geometryFile(input_file, geometry_path):
    
    # Grab a real input file if we have a dummy input file
    if '_dummy' in str(input_file):
        from stellapy.data.input.write_listOfMatchingInputFiles import read_inputFilesInDummyInputFile
        input_files = read_inputFilesInDummyInputFile(input_file) 
        input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
        if len(input_files)==0:
            exit_reason = "fThe geometry file is missing and there are no netcdf files.\n" 
            exit_reason += f"Therefore, the geometry file can not be written for:\n" 
            exit_reason += f"     {input_file}\n" 
            exit_program(exit_reason, write_geometryFile, sys._getframe().f_lineno)   
        input_file = input_files[0] 
         
    # Read the input parameters
    input_parameters = read_inputFile(input_file) 
    y0 = input_parameters['kt_grids_box_parameters']['y0'] 
    svalue = input_parameters['vmec_parameters']['torflux'] 
    vmec_filename = input_parameters['vmec_parameters']['vmec_filename'] 
    nfield_periods = input_parameters['vmec_parameters']['nfield_periods'] 
    
    # Load all the geometry data  
    geometry = read_geometryFromTxtFile(input_file.with_suffix('.geometry'))
    geometry.update(read_geometryFromNetcdfFile(input_file.with_suffix('.out.nc')))
    
    # Load the VMEC or Miller data 
    path = create_dummyPathObject(input_file, vmec_filename)
    if vmec_filename=='wout*.nc': geometry.update(read_millerParameters(path)); geometry['b0'] = 1; geometry['nfp'] = 1 
    if vmec_filename!='wout*.nc': geometry.update(read_woutFromNcFile(path.vmec))
    geometry['sign_B'] = np.sign(geometry['b0'])

    # Calculate the grid divisions set in stella based on the geometry data
    geometry.update(calculate_gridDivisionsAndSize(y0, nfield_periods, geometry, svalue))
    
    # For miller, set the parameters to nan
    if path.vmec_filename=='wout*.nc': 
        for key in ['b0', 'aminor', 'rmajor', 'volume']: geometry[key] = np.nan 
        
    # Calculate the integration weights 
    geometry['vec_z'] = geometry['zed']; geometry['dim_z'] = len(geometry['vec_z'])
    geometry['dl_over_B'] = calculate_dlOverB(geometry['dim_z'], geometry['vec_z'], geometry['jacob'])
    geometry['mu_weights'] = calculate_muWeights(input_file, bmag_psi0=geometry['bmag'])[0,:,:] 
    geometry['vpa_weights'] = calculate_vpaWeights(input_file) 

    # Save the geometry data to the '*.geo.h5' file   
    with h5py.File(geometry_path, 'w') as f:
        for key, value in geometry.items():  
            if isinstance(value, pathlib.Path): value = str(value)
            f.create_dataset(key, data=value)  
                       
    # Track progress
    print('    ----> Saved the geometry file as ' + geometry_path.name)  
    return 
               

 

