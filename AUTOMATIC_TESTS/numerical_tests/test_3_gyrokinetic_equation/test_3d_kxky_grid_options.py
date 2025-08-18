
################################################################################
#               Check all the terms in the gyrokinetic equation                #
################################################################################
# Test the options for the (kx,ky) grid which are 'box' and 'range'
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

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                                BOX NONLINEAR                                 #
#-------------------------------------------------------------------------------
def test_kxky_grid_box_nonlinear(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('kxky_grid_box_nonlinear.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files   
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)    
    print('  -->  The grid_option="box" is evolving correctly during a nonlinear simulation.')
    return
    
#-------------------------------------------------------------------------------
#                                 BOX LINEAR                                   #
#-------------------------------------------------------------------------------
def test_kxky_grid_box_linear(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('kxky_grid_box_linear.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files   
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)    
    print('  -->  The grid_option="box" is evolving correctly during a linear simulation.')
    return
    
#-------------------------------------------------------------------------------
#                                RANGE LINEAR                                  #
#-------------------------------------------------------------------------------
def test_kxky_grid_range_linear(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('kxky_grid_range_linear.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files   
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)    
    print('  -->  The grid_option="range" is evolving correctly during a linear simulation.')
    return
    
#-------------------------------------------------------------------------------
#                                RANGE LINEAR                                  #
#-------------------------------------------------------------------------------
def test_kxky_grid_singlemode_linear(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('kxky_grid_singlemode_linear.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files   
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)    
    print('  -->  The grid_option="range" choosing a single (kx,ky) mode is evolving correctly during a linear simulation.')
    return
    
