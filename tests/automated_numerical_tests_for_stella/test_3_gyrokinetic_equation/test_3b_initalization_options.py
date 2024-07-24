
################################################################################
#               Check all the terms in the gyrokinetic equation                #
################################################################################
# These test are designed to test the terms in the gyrokinetic equation one 
# by one. As a sanity check we start by turning of all terms to check that
# the potential is not being evolved in time. Moreover, the different 
# initialization options for the distribution function # will be tested, to 
# ensure that their definitions has not been changed. To ensure this we need:
#       &physics_flags 
#         nonlinear = .false.
#         include_parallel_streaming = .false.
#         include_mirror = .false.
#       /
#       &time_advance_knobs
#         xdriftknob = 0
#         ydriftknob = 0
#         wstarknob = 0 
#       /
#       &dissipation
#         hyper_dissipation = .false.
#       /   
# And we check the following initialization options:
#  
#       &physics_flags 
#         nonlinear = .true.  # Note that we need nonlinear = .true. or 'noise' will be changed to 'default'
#       /
#      &init_g_knobs
#         ginit_option = 'noise'
#         phiinit = 0.01 
#       /
#
#       &physics_flags 
#         nonlinear = .false.   
#       /
#      &init_g_knobs
#        ginit_option = 'default'
#        phiinit = 0.01 
#        width0 = 1.0
#      /
################################################################################

# Python modules
import os, sys
import pathlib 
import numpy as np
import xarray as xr  

# Package to run stella 
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

#-------------------------------------------------------------------------------
#           Check whether the the 'noise' initialization is the same           #
#-------------------------------------------------------------------------------
def test_whether_init_noise_option_is_the_same(tmp_path, error=False):

    # Input file
    input_filename = 'init_noise.in'  
    
    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path)
     
    # File names  
    local_netcdf_file = tmp_path / input_filename.replace('.in','.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   

    # Check whether the geometry data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        
        # Read the potential axis
        local_phi2 = local_netcdf['phi2']
        expected_phi2 = expected_netcdf['phi2']
                     
        # Check whether the potential is initialized the same
        compare_local_array_with_expected_array(local_phi2, expected_phi2) 
        if not (np.allclose(local_phi2, expected_phi2, equal_nan=True)):
            print('\nERROR: The potential data does not match in the netCDF files.'); error = True
            print('       It is likely that something changed in the initialization of <g> of <phi>.');
            print('\nCompare the potential arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_phi2, expected_phi2) 
        assert (not error), f'The potential data does not match in the netCDF files.'  
                
    print('  -->  The potential is initialized the same as in the previous run (noise).')
    return
    

#-------------------------------------------------------------------------------
#          Check whether the the 'default' initialization is the same          #
#-------------------------------------------------------------------------------
def test_whether_init_default_option_is_the_same(tmp_path, error=False):

    # Input file
    input_filename = 'init_default.in'  
    
    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path)
     
    # File names  
    local_netcdf_file = tmp_path / input_filename.replace('.in','.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   

    # Check whether the geometry data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        
        # Read the potential axis
        local_phi2 = local_netcdf['phi2']
        expected_phi2 = expected_netcdf['phi2']
                     
        # Check whether the potential is initialized the same
        compare_local_array_with_expected_array(local_phi2, expected_phi2) 
        if not (np.allclose(local_phi2, expected_phi2, equal_nan=True)):
            print('\nERROR: The potential data does not match in the netCDF files.'); error = True
            print('       It is likely that something changed in the initialization of <g> of <phi>.');
            print('\nCompare the potential arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_phi2, expected_phi2) 
        assert (not error), f'The potential data does not match in the netCDF files.'  
                
    print('  -->  The potential is initialized the same as in the previous run (default).')
    return
    
