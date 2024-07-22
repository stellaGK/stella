#!/usr/bin/env python3

import numpy as np
import os
import pathlib
import shutil
import subprocess
import xarray as xr
import pytest
import difflib

# Global variables 
input_filename = "miller_geometry.in"  
stella_local_run_directory = 'Not/Run/Yet'
stella_expected_run_directory = 'Not/Run/Yet'

################################################################################
#                 Routines to launch a local stella simulation                 #
################################################################################

def get_test_directory():
    """Get the directory of this test file."""
    return pathlib.Path(__file__).parent

def get_stella_path():
    """Returns the absolute path to the stella executable.
    Can be controlled by setting the STELLA_EXE_PATH environment variable."""
    default_path = get_test_directory() / "../../../stella"
    stella_path = pathlib.Path(os.environ.get("STELLA_EXE_PATH", default_path))
    return stella_path.absolute()

def run_stella(stella_path, input_filename):
    """Run stella with a given input file."""
    subprocess.run([stella_path, input_filename], check=True)

def copy_input_file(input_filename: str, destination):
    """Copy input_filename to destination directory."""
    shutil.copyfile(get_test_directory() / input_filename, destination / input_filename)
    
def convert_byte_array(array):
    """Tool to convert netcdf text arrays, to a text string that we can compare."""
    return "\n".join((str(line, encoding="utf-8").strip() for line in array.data))

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
#                         Run local stella simulation                          #
#-------------------------------------------------------------------------------
def test_run_local_stella_simulation(tmp_path):
    """Run a local stella simulation in a temporary folder""" 
    
    # Set the paths to correct stella netcdf file and the newly generated netcdf file
    global stella_expected_run_directory, stella_local_run_directory 
    stella_local_run_directory = tmp_path  
    stella_expected_run_directory = get_test_directory() 
    
    # Run a local stella simulation
    copy_input_file(input_filename, tmp_path)
    os.chdir(tmp_path)
    run_stella(get_stella_path(), input_filename)
    
    # Always pass this test
    assert True
    print('\n  -->  Successfully ran a local stella simulation.')
    return 
    
#-------------------------------------------------------------------------------
#                    Check whether output files are present                    #
#-------------------------------------------------------------------------------
def test_whether_miller_output_files_are_present():  
    extensions = []
    stem = input_filename.split('.in')[0]
    files = os.listdir(stella_local_run_directory) 
    expected_files = ['millerlocal.miller_geometry.input', 'millerlocal.miller_geometry.output', 'miller_geometry.geometry']
    for file in expected_files: 
        assert (file in files), 'Some output files were not generated when running stella.'
    print(f'  -->  All the expected files (millerlocal.input, millerlocal.output, .geometry) are generated.')
    return 

#-------------------------------------------------------------------------------
#                    Check whether Miller output files match                   #
#-------------------------------------------------------------------------------
def test_whether_miller_output_files_are_correct():  

    # Initialize
    geometry_files_match = True
    miller_input_files_match = True
    miller_output_files_match = True
    
    # File names 
    local_geometry_file = stella_local_run_directory / 'miller_geometry.geometry' 
    expected_geometry_file = stella_expected_run_directory / 'EXPECTED_OUTPUT.geometry' 
    local_miller_input_file = stella_local_run_directory / 'millerlocal.miller_geometry.input' 
    expected_miller_input_file = stella_expected_run_directory / 'EXPECTED_OUTPUT.millerlocal.input' 
    local_miller_output_file = stella_local_run_directory / 'millerlocal.miller_geometry.output' 
    expected_miller_output_file = stella_expected_run_directory / 'EXPECTED_OUTPUT.millerlocal.output' 
     
    # Check whether the millerlocal input files match  
    if os.path.getsize(local_miller_input_file) != os.path.getsize(expected_miller_input_file): miller_input_files_match = False; print('ERROR: Miller input files do not match.')
    elif open(local_miller_input_file,'r').read() != open(expected_miller_input_file,'r').read(): miller_input_files_match = False; print('ERROR: Miller input files do not match.')

    # Check whether the millerlocal output files match  
    if os.path.getsize(local_miller_output_file) != os.path.getsize(expected_miller_output_file): miller_output_files_match = False; print('ERROR: Miller output files do not match.')
    elif open(local_miller_output_file,'r').read() != open(expected_miller_output_file,'r').read(): miller_output_files_match = False; print('ERROR: Miller output files do not match.')

    # Check whether the geometry files match  
    if os.path.getsize(local_geometry_file) != os.path.getsize(expected_geometry_file): geometry_files_match = False; print('ERROR: Geometry files do not match.')
    elif open(local_geometry_file,'r').read() != open(expected_geometry_file,'r').read(): geometry_files_match = False; print('ERROR: Geometry files do not match.')
    
    # If the files don't match, print the differences
    if miller_input_files_match==False:
        print_differences_in_text_files(local_miller_input_file, expected_miller_input_file, name='MILLER INPUT')
        assert False, 'Miller input files do not match.'
    if miller_output_files_match==False:
        print_differences_in_text_files(local_miller_output_file, expected_miller_output_file, name='MILLER OUTPUT')
        assert False, 'Miller output files do not match.'
    if geometry_files_match==False:
        print_differences_in_text_files(local_geometry_file, expected_geometry_file, name='GEOMETRY')
        assert False, 'Geometry files do not match.'

    # If we made it here the test was run correctly
    assert True
    print(f'  -->  Geometry output file matches.')
    return 

#-------------------------------------------------------------------------------
def print_differences_in_text_files(file1, file2, name='', maxlines=10): 
    print(f'\nDIFFERENCE BETWEEN {name} FILES:\n') 
    with open(file1) as f1:
        local_file = f1.readlines()
    with open(file2) as f2:
        expected_file = f2.readlines()  
    if len(local_file)>maxlines:
        print(f'    Warning: the compared files are very long, we limit them to {maxlines} lines')
        local_file = local_file[:maxlines]; expected_file = expected_file[:maxlines]
    for line in difflib.unified_diff(local_file, expected_file, fromfile=str(file1), tofile=str(file2), lineterm=''): 
        print('    ', line) 

#-------------------------------------------------------------------------------
#              Check whether the data in the netcdf file matches               #
#-------------------------------------------------------------------------------
def test_whether_geometry_data_in_netcdf_file_is_correct():
    '''Check that the results are identical to a previous run.
    Note that if the results are expected to change, due to a valid change to stella, 
    the EXPECTED_OUTPUT.out.nc file should be updated on Github. ''' 
     
    # File names  
    local_netcdf_file = stella_local_run_directory / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = stella_expected_run_directory / 'EXPECTED_OUTPUT.out.nc'   

    # Check whether the geometry data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
    
        # Check again whether all the keys are present
        assert set(local_netcdf.keys()) == set(expected_netcdf.keys())
        
        # Relevant keys for the geometry
        geometry_keys = ["bmag", "b_dot_grad_z", "gradpar", "gbdrift", "gbdrift0", "cvdrift", "cvdrift0", "kperp2", \
            "gds2", "gds21", "gds22", "grho", "jacob", "djacdrho", "beta", "q", "shat", "d2qdr2", "drhodpsi", "d2psidr2", "jtwist"]  
        for key in geometry_keys:
        
            # Compare integers and floats
            if expected_netcdf[key].shape == ():
                if key=='nproc': continue
                assert (local_netcdf[key] == expected_netcdf[key]), f'The quantity <{key}> does not match in the netcdf files.'
                continue
                
            # Compare texts (code_info and input_file)
            if expected_netcdf[key].dtype.kind == 'S':
                expected_netcdf_str = convert_byte_array(expected_netcdf[key])
                local_netcdf_str = convert_byte_array(local_netcdf[key])  
                if key!='input_file':
                    assert (local_netcdf_str == expected_netcdf_str), f'The quantity <{key}> does not match in the netcdf files.'
                    
            # Compare data arrays 
            else: 
                assert (np.allclose(local_netcdf[key], expected_netcdf[key], equal_nan=True)), f'The quantity <{key}> does not match in the netcdf files.'
                
    print('  -->  All the geometry data in the netcdf file matches the expected output.')
    return
    

