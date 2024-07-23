
################################################################################
#                     Check the evolution of the potential                     #
################################################################################
# These test are designed to test the evolution of the potential, starting
# with a linear and nonlinear simulation using Miller geometry.
################################################################################

# Python modules
import os, sys
import pathlib 
import numpy as np
import xarray as xr  

# Package to run stella 
module_path = str(pathlib.Path(__file__).parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())
     
################################################################################
#                Run stella locally and check the output files                 #
################################################################################
# The argument of any test function is the temporary path where the test is 
# performed, e.g., test_local_stella_run(tmp_path) executes stella in <tmp_path>. 
################################################################################

#-------------------------------------------------------------------------------
#                Check whether miller linear evolves correctly                 #
#-------------------------------------------------------------------------------
def test_whether_miller_linear_evolves_correctly(tmp_path):

    # Input file name  
    input_filename = 'miller_geometry_linear.in'  
    
    # Run a local stella simulation
    run_local_stella_simulation(input_filename, tmp_path)
     
    # File names  
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   
    
    # Check whether the potential data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
    
        # Check whether all the keys are present
        assert set(local_netcdf.keys()) == set(expected_netcdf.keys()), f'The netcdf file contains different quantities.'
        
        # Read the time axis
        local_time = local_netcdf['t']
        expected_time = expected_netcdf['t']
        
        # Read the potential axis
        local_phi2 = local_netcdf['phi2']
        expected_phi2 = expected_netcdf['phi2'] 
                     
        # Check whether we have the same time and potential data
        assert (np.allclose(local_time, expected_time, equal_nan=True)), f'The time axis does not match in the netcdf files.'
        assert (np.allclose(local_phi2, expected_phi2, equal_nan=True)), f'The potential data does not match in the netcdf files.' 
                
    print(f'\n  -->  The potential is evolving correctly in a linear flux-tube simulation using Miller geometry ({int(local_netcdf["nproc"])} CPUs).')
    return
    
#-------------------------------------------------------------------------------
#               Check whether miller nonlinear evolves correctly               #
#-------------------------------------------------------------------------------
def test_whether_miller_nonlinear_evolves_correctly(tmp_path):

    # Input file name  
    input_filename = 'miller_geometry_nonlinear.in'   
    
    # Run a local stella simulation
    run_local_stella_simulation(input_filename, tmp_path)
     
    # File names  
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   
    
    # Check whether the potential data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
    
        # Check whether all the keys are present
        assert set(local_netcdf.keys()) == set(expected_netcdf.keys()), f'The netcdf file contains different quantities.'
        
        # Read the time axis
        local_time = local_netcdf['t']
        expected_time = expected_netcdf['t']
        
        # Read the potential axis
        local_phi2 = local_netcdf['phi2']
        expected_phi2 = expected_netcdf['phi2'] 
                     
        # Check whether we have the same time and potential data
        assert (np.allclose(local_time, expected_time, equal_nan=True)), f'The time axis does not match in the netcdf files.'
        assert (np.allclose(local_phi2, expected_phi2, equal_nan=True)), f'The potential data does not match in the netcdf files.' 
                
    print(f'\n  -->  The potential is evolving correctly in a nonlinear flux-tube simulation using Miller geometry ({int(local_netcdf["nproc"])} CPUs).')
    return
    
