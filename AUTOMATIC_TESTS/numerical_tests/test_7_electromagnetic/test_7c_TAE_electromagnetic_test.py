
################################################################################
#               Check the TAE physics for Electromagentic stella               #
################################################################################
# These test are designed to test that the code correctly capture the TAE
# (Toroidal Alfven Eigenmode) physics in Electromagnetic stella.
# This input has been benchmarked against GS2. The resolution has been greatly
# reduced for the efficiency of automatic testing, but the physics is still
# captured.
################################################################################

# Python modules
import os, sys
import pathlib 
import numpy as np
import xarray as xr  

# Package to run stella 
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Global variables

input_filename_stem = 'EM'  
stella_local_run_directory = 'Not/Run/Yet' 

#-------------------------------------------------------------------------------
#           Check whether the potential data for the TAE mode matches          #
#-------------------------------------------------------------------------------
def test_TAE_physics_for_Electromagnetic_stella(tmp_path):
    
    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global stella_local_run_directory
    stella_local_run_directory = tmp_path
    
    input_filename = input_filename_stem + '_TAE.in'
    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path)
     
    # File names  
    local_netcdf_file = tmp_path / (input_filename_stem + '_TAE.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   
    compare_local_potential_with_expected_potential_em(local_netcdf_file, expected_netcdf_file, error=False, check_bpar=False)
            
    print(f'  --> The TAE (Kinetic Ballooning Mode) physics is not captured correctly.')

    #-------------------------------------------------------------------------------
    #   Check whether the growth rate and frequency data for the TAE mode matches  #
    #-------------------------------------------------------------------------------
    
    # Check whether the netCDF data matches
    keys = ['omega']
    for key in keys: compare_local_netcdf_quantity_to_expected_netcdf_quantity(local_netcdf_file, expected_netcdf_file, key=key, error=False)

    # If we made it here the test was run correctly
    print(f'  -->  The growth rate and frequency for the TAE are the same.')

    return
