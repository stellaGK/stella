
################################################################################
#                     Check the evolution of the potential                     #
################################################################################
# These test are designed to test the evolution of the potential, starting
# with a linear and nonlinear simulation using Miller geometry.
################################################################################

# Python modules
import pytest
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
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                  Check whether the diagnostics are working                   #
#-------------------------------------------------------------------------------
def test_whether_potential_diagnostics_are_correct(tmp_path, stella_version):

    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global stella_local_run_directory
    stella_local_run_directory = tmp_path
    
    # Run a local stella simulation
    run_local_stella_simulation(input_filename, stella_local_run_directory, stella_version)
     
    # File names  
    global local_netcdf_file, expected_netcdf_file
    local_netcdf_file = stella_local_run_directory / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file_stem}.out.nc'    
    
    # Check whether the potential data matches in the netcdf file
    # We only test the electrostatic quantities, EM stella will be tested later
    keys = ['t', 'phi2', 'phi2_vs_kxky', 'phi_vs_t'] 
    for key in keys: compare_local_netcdf_quantity_to_expected_netcdf_quantity(local_netcdf_file, expected_netcdf_file, key=key, error=False)
    print(f'\n  -->  The potential diagnostics are working.')
    return

#-------------------------------------------------------------------------------
def test_whether_final_fields_diagnostics_are_correct(error=False):
    
    # Txt file names 
    local_file = stella_local_run_directory / f'{input_file_stem}.final_fields' 
    expected_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file_stem}.final_fields'  

    local_fields_data = np.loadtxt(local_file, dtype='float', skiprows=2)
    expected_fields_data = np.loadtxt(expected_file, dtype='float', skiprows=2)

    # Check whether the txt files match  
    for i in range(len(local_fields_data[0,:])):
        # Add in a tolerance here to allow for a small differences in the final fields data due to numerical noise
        # This is needd because there are operations with zero's in the final fields for zed = 0.0, which can lead to
        # which can lead to different values for zeros for different operating systems.
        if not np.allclose(local_fields_data[:,i], expected_fields_data[:,i], rtol = 0.0, atol = 1e-15, equal_nan=True):
            print(f'\nERROR: The final fields do not match in the txt files.'); error = True
            print(f'Compare the final fields in the local and expected txt files:')
            for j in range(len(local_fields_data[:,i])):
                print(f'    column {i}: {local_fields_data[j,i]:14.6e}, {expected_fields_data[j,i]:14.6e}')

    assert (not error), f'The final fields data does not match in the txt files.'
                 
    # If we made it here the test was run correctly 
    print(f'  -->  The final fields diagnostics are working.')
    return 
    
#-------------------------------------------------------------------------------
def test_whether_fluxes_diagnostics_are_correct(error=False):  

    # Txt file names 
    local_file = stella_local_run_directory / f'{input_file_stem}.fluxes' 
    expected_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file_stem}.fluxes'  
    
    # TODO-HT TODO-GA the current master branch has broken momentum flux (its using h instead of f)
    # So test fluxes against cookie's branch, where the heat flux and particle flux do match the master branch
    expected_netcdf_file_fixed_fluxes = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT_FIXEDFLUXES.{input_file_stem}.out.nc'  
    expected_file_fixed_fluxes = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT_FIXEDFLUXES.{input_file_stem}.fluxes'  
    
    # Test is run with 2 species
    dim_species = 2
     
    # Read the fluxes text files
    local_flux_data = np.loadtxt(local_file, dtype='float').reshape(-1, 1+3*dim_species) 
    expected_flux_data = np.loadtxt(expected_file_fixed_fluxes, dtype='float').reshape(-1, 1+3*dim_species)  
    
    # At time step 0 the fluxes are calculated but they are not defined yet
    local_flux_data[0,:] = 0.0
    expected_flux_data[0,:] = 0.0
    
    # Compare the columns (time; pflux_i, pflux_e, vflux_i, vflux_e, qflux_i, qflux_e)
    # Note that the sign of the fluxes has been fixed in newer stellas
    for i in range(len(local_flux_data[0,:])):
        if not np.allclose(local_flux_data[:,i], expected_flux_data[:,i], rtol = 0.0, atol = 1e-15, equal_nan=True):
             if not np.allclose(local_flux_data[:,i], -expected_flux_data[:,i], rtol = 0.0, atol = 1e-15, equal_nan=True): 
                print(f'\nERROR: The fluxes arrays do not match in the txt files.'); error = True
                print(f'Compare the fluxes arrays in the local and expected txt files:')
                for j in range(len(local_flux_data[:,i])):
                 print(f'    column {i}: {local_flux_data[j,i]:14.6e}, {expected_flux_data[j,i]:14.6e}')
    assert (not error), f'The fluxes data does not match in the txt files.'   
    
    # Check whether the netCDF data matches 
    # Note that 'pflx_vs_kxky', 'qflx_vs_kxky', 'vflx_vs_kxky' wont match since old stella calculated
    # flux(kx,ky) as the sum over z, while new stella takes the field line average
    keys = ['pflx', 'qflx', 'vflx', 'pflx_kxky', 'qflx_kxky', 'vflx_kxky', 'pflux_x', 'vflux_x', 'qflux_x'] # Old keys
    keys = ['pflux_vs_s', 'qflux_vs_s', 'vflux_vs_s', 'pflux_x', 'vflux_x', 'qflux_x'] # New keys
    keys += ['pflux_vs_kxkyzs', 'qflux_vs_kxkyzs', 'vflux_vs_kxkyzs'] # New keys
    for key in keys: compare_local_netcdf_quantity_to_expected_netcdf_quantity(local_netcdf_file, \
                            expected_netcdf_file_fixed_fluxes, key=key, error=False) 

    # If we made it here the test was run correctly 
    print(f'  -->  The final fields diagnostics are working.')
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
    
    # TODO-HT Write tests for the new distribution function diagnostics

    # If we made it here the test was run correctly 
    print(f'  -->  The distribution function diagnostics are working.')
    return 
    
#-------------------------------------------------------------------------------
def test_whether_omega_diagnostics_are_correct(tmp_path, stella_version, error=False):

    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global stella_local_run_directory
    stella_local_run_directory = tmp_path
    input_filename = 'miller_linear.in'
    input_file_stem = input_filename.replace(".in","")
    
    # Run a local stella simulation
    run_local_stella_simulation(input_filename, stella_local_run_directory, stella_version)
     
    # Check whether the txt files match (Skip the header line)
    # Compare the columns (time, ky, kx, Re[om], Im[om], Re[omavg], Im[omavg])
    # For the frequency, without nonlinear interactions, we have a lot of noise on the zonal modes
    local_file = stella_local_run_directory / f'{input_file_stem}.omega' 
    expected_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file_stem}.omega'  
    local_omega_data = np.loadtxt(local_file, dtype='float').reshape(-1, 7) 
    expected_omega_data = np.loadtxt(expected_file, dtype='float').reshape(-1, 7)  
    for i in range(len(local_omega_data[0,:])):
        if not np.allclose(local_omega_data[:,i], expected_omega_data[:,i], rtol = 0.0, atol = 1e-12, equal_nan=True):
            print(f'\nERROR: The omega arrays do not match in the txt files.'); error = True
            print(f'Compare the omega arrays in the local and expected txt files:')
            for j in range(len(local_omega_data[:,i])):
                print(f'    column {i}: {local_omega_data[j,i]:14.6e}, {expected_omega_data[j,i]:14.6e}')
    assert (not error), f'The omega data does not match in the txt files.'   
    
    # Check whether the netCDF data matches 
    keys = ['omega']
    local_netcdf_file = stella_local_run_directory / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file_stem}.out.nc'
    for key in keys: compare_local_netcdf_quantity_to_expected_netcdf_quantity(local_netcdf_file, expected_netcdf_file, key=key, error=False) 

    # If we made it here the test was run correctly 
    print(f'  -->  The omega diagnostics are working.')
    return 
