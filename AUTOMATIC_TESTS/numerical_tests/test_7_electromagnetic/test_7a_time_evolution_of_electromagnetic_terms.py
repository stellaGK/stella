
################################################################################
#               Check all the terms in the gyrokinetic equation                #
################################################################################
# These test are designed to test the terms in the gyrokinetic equation one 
# by one for the case of electromagnetic stella.
# As a sanity check we start by turning of all terms to 
# check that the potential is not being evolved in time. Then we turn the terms 
# on one by one linearly, then all together linearly, and finally we test them 
# when every term is turned on nonlinearly. The flags we toggle are:
#       &physics_flags 
#         nonlinear = .false. / .true.
#         include_parallel_streaming = .false. / .true.
#         include_mirror = .false. / .true.
#       /
#       &time_advance_knobs
#         xdriftknob = 0.0 / 1.0
#         ydriftknob = 0.0 / 1.0
#         wstarknob = 0.0 / 1.0
#       /
#       &knobs
#         stream_implicit = .false. / .true.
#         mirror_implicit = .false. / .true.
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

input_filename_stem = 'EM'
stella_local_run_directory = 'Not/Run/Yet'

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#           Check whether the potential data does not evolve in time           #
#-------------------------------------------------------------------------------
def test_each_gyrokinetic_term_for_electromagnetic(tmp_path, stella_version):
    
    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global stella_local_run_directory
    stella_local_run_directory = tmp_path
    
    input_filename = input_filename_stem + '_no_time_evolution.in'
    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
     
    # File names  
    local_netcdf_file = tmp_path / (input_filename_stem + '_no_time_evolution.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   
    compare_local_potential_with_expected_potential_em(local_netcdf_file, expected_netcdf_file, error=False)
            
    print(f'  -->  Without gyrokinetic terms the potential in EM stella does not evolve in time.')
    return 
#-------------------------------------------------------------------------------
#     Check whether the implicit parallel streaming term evolves correctly     #
#------------------------------------------------------------------------------
def test_each_gyrokinetic_term_for_electromagnetic_parallel_streaming_implicit(tmp_path, stella_version):

    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_parallel_streaming_implicit.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_parallel_streaming_implicit.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   
    compare_local_potential_with_expected_potential_em(local_netcdf_file, expected_netcdf_file, error=False)

    print('  -->  The parallel streaming term for EM stella, implemented implicitly, is evolving correctly.') 

#-------------------------------------------------------------------------------
#     Check whether the explicit parallel streaming term evolves correctly     #
#-------------------------------------------------------------------------------
def test_each_gyrokinetic_term_for_electromagnetic_parallel_streaming_explicit(tmp_path, stella_version):

    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_parallel_streaming_explicit.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_parallel_streaming_explicit.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   
    compare_local_potential_with_expected_potential_em(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  The parallel streaming term for EM stella, implemented explicitly, is evolving correctly.')

#-------------------------------------------------------------------------------
#           Check whether the implicit mirror term evolves correctly           #
#-------------------------------------------------------------------------------
def test_each_gyrokinetic_term_for_electromagnetic_mirror_implicit(tmp_path, stella_version):

    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_mirror_implicit.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_mirror_implicit.out.nc')  
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    compare_local_potential_with_expected_potential_em(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  The mirror term for EM stella, implemented implicitly, is evolving correctly.')

#-------------------------------------------------------------------------------
#           Check whether the explicit mirror term evolves correctly           #
#-------------------------------------------------------------------------------
def test_each_gyrokinetic_term_for_electromagnetic_mirror_explicit(tmp_path, stella_version):

    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_mirror_explicit.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_mirror_explicit.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    compare_local_potential_with_expected_potential_em(local_netcdf_file, expected_netcdf_file, error=False)  
    print('  -->  The mirror term for EM stella, implemented explicitly, is evolving correctly.')

#-------------------------------------------------------------------------------
#          Check whether the diagmagnetic drift term evolves correctly         #
#-------------------------------------------------------------------------------
def test_each_gyrokinetic_term_for_electromagnetic_diagmagnetic_drift(tmp_path, stella_version):

    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_diagmagnetic_drift.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_diagmagnetic_drift.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    compare_local_potential_with_expected_potential_em(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  The diagmagnetic drift term for EM stella is evolving correctly.')


#-------------------------------------------------------------------------------
#          Check whether the magnetic drifts term evolves correctly         #
#-------------------------------------------------------------------------------
def test_each_gyrokinetic_term_for_electromagnetic_magnetic_drifts(tmp_path, stella_version):

    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_magnetic_drifts.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_magnetic_drifts.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    compare_local_potential_with_expected_potential_em(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  The magnetic drift terms for EM is evolving correctly.')
    
#-------------------------------------------------------------------------------
#              Check whether the nonlinear term evolves correctly              #
#-------------------------------------------------------------------------------
def test_each_gyrokinetic_term_for_electromagnetic_nonlinear_term(tmp_path, stella_version):

    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_nonlinear_term.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)

    # Compare phi2(t) in the netCDF files  
    local_netcdf_file = tmp_path / (input_filename_stem + '_nonlinear_term.out.nc')
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'
    compare_local_potential_with_expected_potential_em(local_netcdf_file, expected_netcdf_file, error=False)
    print('  -->  The nonlinear term is being evolved correctly .')

#-------------------------------------------------------------------------------
#       Check whether the all the terms combined term evolves correctly        #
#-------------------------------------------------------------------------------
def test_each_gyrokinetic_term_for_electromagnetic_kxky_grid_box_linear(tmp_path, stella_version):

    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_kxky_grid_box_linear.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_kxky_grid_box_linear.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    compare_local_potential_with_expected_potential_em(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  All terms are running as they should when simulated together for EM when running linearly in box mode.') 

#-------------------------------------------------------------------------------
#              Check whether the nonlinear term evolves correctly              #
#-------------------------------------------------------------------------------
def test_each_gyrokinetic_term_for_electromagnetic_kxky_grid_box_nonlinear(tmp_path, stella_version):

    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_kxky_grid_box_nonlinear.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_kxky_grid_box_nonlinear.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    compare_local_potential_with_expected_potential_em(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  All terms are running as they should when simulated together for EM stella when running nonlinearly in box mode.') 

    return
