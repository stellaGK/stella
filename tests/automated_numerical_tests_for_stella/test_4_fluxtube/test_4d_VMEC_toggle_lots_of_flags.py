
################################################################################
#                     Check the evolution of the potential                     #
################################################################################
# These test are designed to test the evolution of the potential, checking
# a nonlinear simulation using VMEC geometry, checking all flags at once.
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
vmec_filename = 'wout_w7x_kjm.nc'

################################################################################
#                Run stella locally and check the output files                 #
################################################################################
# The argument of any test function is the temporary path where the test is 
# performed, e.g., test_local_stella_run(tmp_path) executes stella in <tmp_path>. 
################################################################################

#-------------------------------------------------------------------------------
#               Check whether VMEC nonlinear evolves correctly                 #
#-------------------------------------------------------------------------------
def offtest_whether_a_mistake_was_introduced_in_an_extra_flag_1(tmp_path, error=False):
    # TODO-HT Turn test back on

    # Input file name  
    input_filename = 'test_all_flags1.in'  
    
    # Run a local stella simulation
    run_local_stella_simulation(input_filename, tmp_path, vmec_filename)
     
    # File names  
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   
    
    # Check whether the potential data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        
        # Read the time axis
        local_time = local_netcdf['t']
        expected_time = expected_netcdf['t']
        
        # Read the potential axis
        local_phi2 = local_netcdf['phi2']
        expected_phi2 = expected_netcdf['phi2'] 
        
        # Check whether we have the same time and potential data
        if not (np.allclose(local_time, expected_time, equal_nan=True)):
            print('\nERROR: The time axis does not match in the netCDF files.'); error = True
            print('\nCompare the time arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_time, expected_time)  
        if not (np.allclose(local_phi2, expected_phi2, equal_nan=True)):
            print('\nERROR: The potential data does not match in the netCDF files.'); error = True 
            print('\nCompare the potential arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_phi2, expected_phi2) 
        assert (not error), f'The potential data does not match in the netCDF files.'  
                
    print(f'\n  -->  The potential is evolving correctly in a nonlinear flux-tube simulation using VMEC geometry when toggling usefull flags ({int(local_netcdf["nproc"])} CPUs).')
    return
 
#-------------------------------------------------------------------------------
def offtest_whether_a_mistake_was_introduced_in_an_extra_flag_2(tmp_path, error=False):
    # TODO-HT Turn test back on

    # Input file name  
    input_filename = 'test_all_flags2.in'  
    
    # Run a local stella simulation
    run_local_stella_simulation(input_filename, tmp_path, vmec_filename)
     
    # File names  
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   
    
    # Check whether the potential data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        
        # Read the time axis
        local_time = local_netcdf['t']
        expected_time = expected_netcdf['t']
        
        # Read the potential axis
        local_phi2 = local_netcdf['phi2']
        expected_phi2 = expected_netcdf['phi2'] 
        
        # Check whether we have the same time and potential data
        if not (np.allclose(local_time, expected_time, equal_nan=True)):
            print('\nERROR: The time axis does not match in the netCDF files.'); error = True
            print('\nCompare the time arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_time, expected_time)  
        if not (np.allclose(local_phi2, expected_phi2, equal_nan=True)):
            print('\nERROR: The potential data does not match in the netCDF files.'); error = True 
            print('\nCompare the potential arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_phi2, expected_phi2) 
        assert (not error), f'The potential data does not match in the netCDF files.'  
                
    print(f'\n  -->  The potential is evolving correctly in a nonlinear flux-tube simulation using VMEC geometry when toggling usefull flags ({int(local_netcdf["nproc"])} CPUs).')
    return
