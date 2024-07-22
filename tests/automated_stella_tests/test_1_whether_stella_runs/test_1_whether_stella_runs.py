#!/usr/bin/env python3

import numpy as np
import os
import pathlib
import shutil
import subprocess
import xarray as xr
import pytest

# Global variables 
input_filename = "small_stella_run.in" 
correct_netcdf_file = "EXPECTED_OUTPUT.out.nc"
stella_run_directory = 'Not/Run/Yet'

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
    global local_netcdf_file, correct_netcdf_file, stella_run_directory 
    stella_run_directory = tmp_path
    correct_netcdf_file = get_test_directory() / correct_netcdf_file 
    local_netcdf_file = stella_run_directory / input_filename.replace("in", "out.nc")
    
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
def test_check_whether_all_output_files_are_present():  
    extensions = []
    stem = input_filename.split('.in')[0]
    for file in os.listdir(stella_run_directory): 
        if file.startswith('.'): continue
        elif file.startswith(stem): extensions.append('.'.join(file.split('.')[1:]))
        else: extensions.append(file)
    expected_files = ['in', 'out.nc', 'geometry', 'fluxes', 'omega', 'out', 'final_fields', 'species.input', 'error']
    for file_extension in expected_files:
        assert (file_extension in extensions), 'Some output files were not generated when running stella.'
    print(f'  -->  All the expected files (.in, .geometry, .out.nc, .fluxes, .omega, ...) are generated.')
    return 
        
#-------------------------------------------------------------------------------
#         Check whether all quantities are present in the netcdf file          #
#-------------------------------------------------------------------------------
def test_correct_quantities_are_present_in_netcdf_file():
    """Run a short test case to check that the output is generated with
    the expected fields, and the expected number of timesteps. """ 
    with xr.open_dataset(local_netcdf_file) as local_netcdf:
        expected_keys = ["kx", "ky", "tube", "zed", "alpha", "vpa", "mu", "species", "t", "char10", "char200", "ri"]  
        for key in expected_keys:
            assert (key in local_netcdf), f'The quantity {key} could not be found in the netcdf file.'
        print(f'  -->  The netcdf file contains the following dimensions: {expected_keys}.')  
        expected_keys = ["charge", "mass", "dens", "temp", "tprim", "fprim", "vnew", "type_of_species"] 
        for key in expected_keys:
            assert (key in local_netcdf), f'The quantity {key} could not be found in the netcdf file.'
        print(f'  -->  The netcdf file contains the following species information: {expected_keys}.')
        expected_keys = ["bmag", "b_dot_grad_z", "gradpar", "gbdrift", "gbdrift0", "cvdrift", "cvdrift0", "kperp2", \
            "gds2", "gds21", "gds22", "grho", "jacob", "djacdrho", "beta", "q", "shat", "d2qdr2", "drhodpsi", "d2psidr2", "jtwist"] 
        for key in expected_keys:
            assert (key in local_netcdf), f'The quantity {key} could not be found in the netcdf file.'
        print(f'  -->  The netcdf file contains the following geometry information: {expected_keys}.')
        expected_keys = ["t", "omega", "phi2", "apar2", "bpar2", "pflx", "vflx", "phi_vs_t", "phi2_vs_kxky", "pflx_vs_kxky", "vflx_vs_kxky", "qflx_vs_kxky",\
            "pflux_x", "vflux_x", "qflux_x", "dens_x", "upar_x", "temp_x", "pflx_kxky", "vflx_kxky", "qflx_kxky", "gvmus", "gzvs"] 
        for key in expected_keys:
            assert (key in local_netcdf), f'The quantity {key} could not be found in the netcdf file.' 
        print(f'  -->  The netcdf file contains the following diagnostics information: {expected_keys}.')
    return 
