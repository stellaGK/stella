
################################################################################
#               Check all the terms in the gyrokinetic equation                #
################################################################################
# These test are designed to test the terms in the gyrokinetic equation one 
# by one. Here we test the implicit parallel streaming term:
#       &physics_flags 
#         nonlinear = .true.
#         include_parallel_streaming = .true.
#         include_mirror = .false.
#       /
#       &time_advance_knobs
#         xdriftknob = 0
#         ydriftknob = 0
#         wstarknob = 0 
#       /
#       &knobs 
#         stream_implicit = .false.
#         mirror_implicit = .false.
#         drifts_implicit = .false.
#       /
#       &dissipation
#         hyper_dissipation = .false.
#       /  
################################################################################

import numpy as np
import os
import pathlib
import shutil
import subprocess
import xarray as xr
import pytest
import difflib

# Global variables  
input_filename = 'parallel_streaming_explicit.in'  

################################################################################
#                 Routines to launch a local stella simulation                 #
################################################################################

def get_test_directory():
    '''Get the directory of this test file.'''
    return pathlib.Path(__file__).parent

def get_stella_path():
    '''Returns the absolute path to the stella executable.
    Can be controlled by setting the STELLA_EXE_PATH environment variable.'''
    default_path = get_test_directory() / '../../../stella'
    stella_path = pathlib.Path(os.environ.get('STELLA_EXE_PATH', default_path))
    return stella_path.absolute()

def run_stella(stella_path, input_filename):
    '''Run stella with a given input file.'''
    subprocess.run([stella_path, input_filename], check=True)

def copy_input_file(input_filename: str, destination):
    '''Copy input_filename to destination directory.'''
    shutil.copyfile(get_test_directory() / input_filename, destination / input_filename)
    
def convert_byte_array(array):
    '''Tool to convert netcdf text arrays, to a text string that we can compare.'''
    return '\n'.join((str(line, encoding='utf-8').strip() for line in array.data))

################################################################################
#                Run stella locally and check the output files                 #
################################################################################
# The argument of any test function is the temporary path where the test is 
# performed, e.g., test_local_stella_run(tmp_path) executes stella in <tmp_path>.
# The test are performed from top to bottom, so we use the first test to run
# a local stella simulation in a temporary folder, and later we can access 
# the temporary data to run tests on it.
################################################################################

#-------------------------------------------------------------------------------
#     Check whether the explicit parallel streaming term evolves correctly     #
#-------------------------------------------------------------------------------
def test_whether_parallel_streaming_explicit_evolves_correctly(tmp_path):

    # Set the paths to correct stella netcdf file and the newly generated netcdf file
    stella_local_run_directory = tmp_path
    stella_expected_run_directory = get_test_directory()  
    input_file_stem = input_filename.replace('.in','')
    
    # Run a local stella simulation
    copy_input_file(input_filename, tmp_path)
    os.chdir(tmp_path)
    run_stella(get_stella_path(), input_filename)
     
    # File names  
    local_netcdf_file = stella_local_run_directory / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = stella_expected_run_directory / f'EXPECTED_OUTPUT.{input_file_stem}.out.nc'   

    # Check whether the geometry data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
    
        # Check again whether all the keys are present
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
                
    print('  -->  The parallel streaming term, implemented explicitly, is evolving correctly.')
    return
    
