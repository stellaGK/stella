
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
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Global variables
input_filename = 'miller_nonlinear.in'  
input_file_stem = input_filename.replace(".in","")
stella_local_run_directory = 'Not/Run/Yet'
local_netcdf_file = 'Not/Run/Yet'
expected_netcdf_file = 'Not/Run/Yet'

#-------------------------------------------------------------------------------
#                  Check whether the diagnostics are working                   #
#-------------------------------------------------------------------------------
def test_whether_potential_diagnostics_are_correct(tmp_path):    

    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global stella_local_run_directory
    stella_local_run_directory = tmp_path
    
    # Run a local stella simulation
    run_local_stella_simulation(input_filename, stella_local_run_directory)
     
    # File names  
    global local_netcdf_file, expected_netcdf_file
    local_netcdf_file = stella_local_run_directory / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file_stem}.out.nc'    
    
    # Check whether the potential data matches in the netcdf file
    keys = ['t', 'phi2', 'phi2_vs_kxky', 'phi_vs_t', 'apar_vs_t', 'bpar_vs_t']
    for key in keys: compare_local_netcdf_quantity_to_expected_netcdf_quantity(local_netcdf_file, expected_netcdf_file, key=key, error=False)
    print(f'\n  -->  The potential diagnostics are working.')
    return

#-------------------------------------------------------------------------------
def test_whether_final_fields_diagnostics_are_correct(error=False):   
    
    # Txt file names 
    local_file = stella_local_run_directory / f'{input_file_stem}.final_fields' 
    expected_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file_stem}.final_fields'  
     
    # Check whether the txt files match  
    compare_local_txt_with_expected_txt(local_file, expected_file, name='Final fields', error=False)

    # If we made it here the test was run correctly 
    print(f'  -->  The final fields diagnostics are working.')
    return 
    
#-------------------------------------------------------------------------------
def test_whether_fluxes_diagnostics_are_correct(error=False):   
    
    # Txt file names 
    local_file = stella_local_run_directory / f'{input_file_stem}.fluxes' 
    expected_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file_stem}.fluxes'  
     
    # Check whether the txt files match  
    compare_local_txt_with_expected_txt(local_file, expected_file, name='Fluxes txt', error=False)
    
    # Check whether the netCDF data matches 
    keys = ['pflx', 'vflx', 'qflx', 'pflx_kxky', 'qflx_kxky', 'vflx_kxky', 'pflx_vs_kxky', 'vflx_vs_kxky', 'qflx_vs_kxky', 'pflux_x', 'vflux_x', 'qflux_x']
    for key in keys: compare_local_netcdf_quantity_to_expected_netcdf_quantity(local_netcdf_file, expected_netcdf_file, key=key, error=False) 

    # If we made it here the test was run correctly 
    print(f'  -->  The final fields diagnostics are working.')
    return 
    
#-------------------------------------------------------------------------------
def test_whether_omega_diagnostics_are_correct(error=False):   
    
    # Txt file names 
    local_file = stella_local_run_directory / f'{input_file_stem}.omega' 
    expected_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file_stem}.omega'  
     
    # Check whether the txt files match  
    compare_local_txt_with_expected_txt(local_file, expected_file, name='Omega txt', error=False)
    
    # Check whether the netCDF data matches 
    keys = ['omega']
    for key in keys: compare_local_netcdf_quantity_to_expected_netcdf_quantity(local_netcdf_file, expected_netcdf_file, key=key, error=False) 

    # If we made it here the test was run correctly 
    print(f'  -->  The omega diagnostics are working.')
    return 
    
#-------------------------------------------------------------------------------
def test_whether_moments_diagnostics_are_correct(error=False):    
    
    # Check whether the netCDF data matches 
    keys = ['density', 'upar', 'temperature', 'spitzer2', 'dens_x', 'upar_x', 'temp_x']
    for key in keys: compare_local_netcdf_quantity_to_expected_netcdf_quantity(local_netcdf_file, expected_netcdf_file, key=key, error=False) 

    # If we made it here the test was run correctly 
    print(f'  -->  The moments diagnostics are working.')
    return 
    
#-------------------------------------------------------------------------------
def test_whether_distribution_function_diagnostics_are_correct(error=False):    
    
    # Check whether the netCDF data matches 
    keys = ['gvmus', 'gzvs']
    for key in keys: compare_local_netcdf_quantity_to_expected_netcdf_quantity(local_netcdf_file, expected_netcdf_file, key=key, error=False) 

    # If we made it here the test was run correctly 
    print(f'  -->  The distribution function diagnostics are working.')
    return 
