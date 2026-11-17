################################################################################
#                          Check the z-pinch geometry                          #
################################################################################
# Test all the geometry arrays when using a z-pinch equilibrium. 
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
input_filename = 'zpinch_geometry.in'

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")
    
#-------------------------------------------------------------------------------
#                         Run local stella simulation                          #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def stella_run(tmp_path_factory, stella_version):
    '''Run a local stella simulation, which is shared by all the tests in this 
    module, so that each test can also be run on its own.'''
    if stella_version!='master': pytest.skip('The z-pinch geometry did not exist in stella versions 0.5, 0.6 and 0.7.')
    tmp_path = tmp_path_factory.mktemp('zpinch_geometry')
    return run_local_stella_simulation(input_filename, tmp_path, stella_version)

#-------------------------------------------------------------------------------
#                    Check whether output files are present                    #
#-------------------------------------------------------------------------------
def test_whether_zpinch_output_files_are_present(stella_run, error=False):
    stella_local_run_directory = stella_run['tmp_path']
    
    # Gather the output files generated during the local stella run inside <tmp_path>
    local_files = os.listdir(stella_local_run_directory)
    
    # Create a list of the output files we expect when stella has been run 
    expected_files = ['zpinch_geometry.geometry']
    
    # Check whether all these output files are present
    for expected_file in expected_files:
        if not (expected_file in local_files):
            print(f'ERROR: The "{expected_file}" output file was not generated when running stella.'); error = True

    # The <pytest> module will check whether all <assert> statements are true,
    # if it runs into a statement which is false, the test will be labeled as
    # "Failed" and the string in the second argument of the <assert> statement 
    # will be printed to the command prompt as "AssertionError: <error message>"
    assert (not error), f'Some output files were not generated when running stella.'
    
    # The test will stop as soon as an <assert> error was triggered, hence we
    # will only reach this final line of code if the test ran successfully
    print(f'  -->  All the expected files are generated.')
    return

#-------------------------------------------------------------------------------
#                    Check whether zpinch output files match                   #
#-------------------------------------------------------------------------------
def test_whether_zpinch_output_files_are_correct(stella_run, stella_version):  
    '''Check that the results are identical to a previous run.'''
    stella_local_run_directory = stella_run['tmp_path']
    

    # File names 
    local_geometry_file = stella_local_run_directory / 'zpinch_geometry.geometry' 
    expected_geometry_file = get_stella_expected_run_directory() / 'EXPECTED_OUTPUT.zpinch_geometry.geometry' 

    # Check whether the txt files match  
    compare_geometry_files(local_geometry_file, expected_geometry_file, error=False)

    # If we made it here the test was run correctly 
    print(f'  -->  Geometry output file matches.')
    return 
    
#-------------------------------------------------------------------------------
#              Check whether the data in the netcdf file matches               #
#-------------------------------------------------------------------------------
def test_whether_geometry_data_in_netcdf_file_is_correct(stella_run, error=False):
    compare_geometry_in_netcdf_files(stella_run, error=False)
    print('  -->  All geometry data in the netcdf file matches the expected output.')
    return
