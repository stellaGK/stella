################################################################################
#                           Check the Miller geometry                          #
################################################################################
# Test all the geometry arrays when using a Miller equilibrium. 
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
input_filename = 'miller_geometry.in'
stella_local_run_directory = 'Not/Run/Yet'
input_file = input_filename.replace('.in','')
miller_file_name = 'geometry_miller'
run_data = {}

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                    Check whether output files are present                    #
#-------------------------------------------------------------------------------
def test_whether_miller_output_files_are_present(tmp_path, stella_version, error=False):
    
    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global stella_local_run_directory, miller_file_name, input_file, run_data
    stella_local_run_directory = tmp_path
    
    # Run stella inside of <tmp_path> based on <input_filename>
    run_data = run_local_stella_simulation(input_filename, tmp_path, stella_version)
    
    # Gather the output files generated during the local stella run inside <tmp_path>
    local_files = os.listdir(stella_local_run_directory)
    
    # Create a list of the output files we expect when stella has been run 
    expected_files = [f'{miller_file_name}.{input_file}.input', f'{miller_file_name}.{input_file}.output', f'{input_file}.geometry']; new_names = True
    
    # Check whether all these output files are present
    for expected_file in expected_files:
        if not (expected_file in local_files):
            print(f'ERROR: The "{expected_file}" output file was not generated when running stella.'); new_names = False
            
    # Old stella 
    if new_names==False:
        expected_files = [f'millerlocal.{input_file}.input', f'millerlocal.{input_file}.output', '{input_file}.geometry'] 
        for expected_file in expected_files:
            if not (expected_file in local_files):
                print(f'ERROR: The "{expected_file}" output file was not generated when running stella.'); error = True
        miller_file_name = 'millerlocal'
            
    # The <pytest> module will check whether all <assert> statements are true,
    # if it runs into a statement which is false, the test will be labeled as
    # "Failed" and the string in the second argument of the <assert> statement 
    # will be printed to the command prompt as "AssertionError: <error message>"
    assert (not error), f'Some output files were not generated when running stella.'
    
    # The test will stop as soon as an <assert> error was triggered, hence we
    # will only reach this final line of code if the test ran successfully
    print(f'  -->  All the expected files (millerlocal.input, millerlocal.output, .geometry) are generated.')

    #-------------------------------------------------------------------------------
    #                    Check whether Miller output files match                   #
    #-------------------------------------------------------------------------------
    '''Check that the results are identical to a previous run.'''
    
    # File names 
    local_geometry_file = stella_local_run_directory / f'{input_file}.geometry' 
    expected_geometry_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file}.geometry' 
    local_miller_input_file = stella_local_run_directory / f'{miller_file_name}.{input_file}.input' 
    expected_miller_input_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file}.millerlocal.input' 
    local_miller_output_file = stella_local_run_directory / f'{miller_file_name}.{input_file}.output' 
    expected_miller_output_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file}.millerlocal.output' 
    
    # Check whether the txt files match
    # Note that only 1 digit matches in the geometry file, but all digits match in the miller output file
    compare_geometry_files(local_geometry_file, expected_geometry_file, error=False, digits=1)
    shat = compare_miller_input_files(local_miller_input_file, expected_miller_input_file, error=False)
    compare_miller_output_files(local_miller_output_file, expected_miller_output_file, shat=shat, error=False)
    
    # If we made it here the test was run correctly 
    print(f'  -->  Geometry output file matches.')

    #-------------------------------------------------------------------------------
    #              Check whether the data in the netcdf file matches               #
    #-------------------------------------------------------------------------------
    compare_geometry_in_netcdf_files(run_data, error=False)
    print('  -->  All Miller geometry data in the netcdf file matches the expected output.')

    return
    
     
    

