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
    filename = input_filename if stella_version=='master' else input_filename.replace('.in', f'_v{stella_version}.in')
    run_data = run_local_stella_simulation(filename, tmp_path_factory.mktemp('miller_geometry'), stella_version)
    run_data['input_file_stem'] = filename.replace('.in','')
    # Older stella versions named the Miller files <millerlocal.*> instead of <geometry_miller.*>
    new_name_exists = (run_data['tmp_path'] / f'geometry_miller.{run_data["input_file_stem"]}.input').exists()
    run_data['miller_file_name'] = 'geometry_miller' if new_name_exists else 'millerlocal'
    return run_data

#-------------------------------------------------------------------------------
#                    Check whether output files are present                    #
#-------------------------------------------------------------------------------
def test_whether_miller_output_files_are_present(stella_run, error=False):  
    
    # Gather the output files generated during the local stella run
    local_files = os.listdir(stella_run['tmp_path'])
    input_file, miller_file_name = stella_run['input_file_stem'], stella_run['miller_file_name']
    
    # Check whether all the output files we expect are present
    expected_files = [f'{miller_file_name}.{input_file}.input', f'{miller_file_name}.{input_file}.output', f'{input_file}.geometry']
    for expected_file in expected_files:
        if not (expected_file in local_files):
            print(f'ERROR: The "{expected_file}" output file was not generated when running stella.'); error = True
    assert (not error), f'Some output files were not generated when running stella.'
    print(f'  -->  All the expected files ({miller_file_name}.input, {miller_file_name}.output, .geometry) are generated.')
    return 

#-------------------------------------------------------------------------------
#                    Check whether Miller output files match                   #
#-------------------------------------------------------------------------------
def test_whether_miller_output_files_are_correct(stella_run):
    '''Check that the results are identical to a previous run.'''
    
    # File names
    local_directory = stella_run['tmp_path']
    input_file, miller_file_name = stella_run['input_file_stem'], stella_run['miller_file_name']
    local_geometry_file = local_directory / f'{input_file}.geometry' 
    expected_geometry_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.miller_geometry.geometry' 
    local_miller_input_file = local_directory / f'{miller_file_name}.{input_file}.input' 
    expected_miller_input_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.miller_geometry.millerlocal.input' 
    local_miller_output_file = local_directory / f'{miller_file_name}.{input_file}.output' 
    expected_miller_output_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.miller_geometry.millerlocal.output'
    
    # Compare text files (first check the input file to save <shat>)
    compare_geometry_files(local_geometry_file, expected_geometry_file, error=False)
    shat = compare_miller_input_files(local_miller_input_file, expected_miller_input_file, error=False)
    compare_miller_output_files(local_miller_output_file, expected_miller_output_file, shat=shat, error=False)
    print(f'  -->  Geometry output file matches.')
    return

#-------------------------------------------------------------------------------
#              Check whether the data in the netcdf file matches               #
#-------------------------------------------------------------------------------
def test_whether_miller_geometry_data_in_netcdf_file_is_correct(stella_run, error=False): 
    compare_geometry_in_netcdf_files(stella_run, error=False)  
    print('  -->  All Miller geometry data in the netcdf file matches the expected output.')
    return
    

