################################################################################
#                       Check the non-twisting flux tube                       #
################################################################################
# These tests check the non-twisting flux tube (boundary_option = 'non-twisting',
# Ball and Brunner 2021) in a nonlinear simulation using VMEC geometry, and in a
# linear electromagnetic simulation with the implicit drifts and parallel streaming.
# The expected output was made with stella-fixes-13 (d4c457095d), before OpenMP.
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
vmec_filename = 'wout_w7x_kjm.nc'

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                Compare the fields and the full state of stella               #
#-------------------------------------------------------------------------------
def check_non_twisting_run(input_filename, tmp_path, error=False):

    # File names
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc')
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'

    # Check whether the time traces of the fields match in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        if not (np.allclose(local_netcdf['t'], expected_netcdf['t'], equal_nan=True)):
            print('\nERROR: The time axis does not match in the netCDF files.'); error = True
            compare_local_array_with_expected_array(local_netcdf['t'], expected_netcdf['t'])
        for key in ['phi2', 'apar2', 'bpar2']:
            if key not in expected_netcdf.variables: continue
            if not (np.allclose(local_netcdf[key], expected_netcdf[key], equal_nan=True)):
                print(f'\nERROR: The {key} data does not match in the netCDF files.'); error = True
                compare_local_array_with_expected_array(local_netcdf[key], expected_netcdf[key], name=key)
        nproc = int(local_netcdf["nproc"])
    assert (not error), f'The fields do not match in the netCDF files.'
    compare_full_state_with_expected(local_netcdf_file, expected_netcdf_file)
    return nproc

#-------------------------------------------------------------------------------
#          Nonlinear VMEC simulation in the non-twisting flux tube             #
#-------------------------------------------------------------------------------
def test_whether_vmec_nonlinear_non_twisting_evolves_correctly(tmp_path, stella_version):

    # Input file name
    input_filename = 'nontwisting_vmec_nonlinear.in'

    # Run a local stella simulation
    run_local_stella_simulation(input_filename, tmp_path, stella_version, vmec_filename)
    nproc = check_non_twisting_run(input_filename, tmp_path)

    print(f'\n  -->  The potential is evolving correctly in a nonlinear simulation using VMEC geometry in the non-twisting flux tube ({nproc} CPUs).')
    return

#-------------------------------------------------------------------------------
#    Linear electromagnetic simulation with implicit drifts, non-twisting      #
#-------------------------------------------------------------------------------
def test_whether_electromagnetic_implicit_non_twisting_evolves_correctly(tmp_path, stella_version):

    # Input file name
    input_filename = 'nontwisting_miller_EM_implicit.in'

    # Run a local stella simulation
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    nproc = check_non_twisting_run(input_filename, tmp_path)

    print(f'\n  -->  The fields are evolving correctly in an electromagnetic simulation with implicit drifts in the non-twisting flux tube ({nproc} CPUs).')
    return
