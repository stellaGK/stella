
################################################################################
#                Check the new explicit time integrators (regression)           #
################################################################################
# These tests regression-check the new explicit time-advance schemes added in
# the 'numerical upgrades' port:
#     SSP(3,2), SSP(4,2), SSP(5,2), SSP(4,3), SSP(9,3) and RK6-4.
#
# Each scheme is run on a single, minimal linear problem (only the explicit
# parallel-streaming term active, single-ky-mode Gaussian init). The resulting
# potential <phi2>(t) is compared against a stored golden netCDF file
# (EXPECTED_OUTPUT.<input>.out.nc). This catches any change to the scheme
# coefficients or the time-advance wiring.
#
# The relevant parameters in the input files are:
#       &gyrokinetic_terms
#         include_parallel_streaming = .true.
#         include_mirror = .false.
#         include_nonlinear = .false.
#       /
#       &numerical_algorithms
#         explicit_algorithm = '<scheme>'
#       /
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
#                       SSP(3,2) explicit integrator                           #
#-------------------------------------------------------------------------------
def test_whether_ssp32_integrator_reproduces_expected_output(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('integrator_ssp32.in', tmp_path, stella_version)

    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The SSP(3,2) explicit integrator reproduces the expected potential.')
    return

#-------------------------------------------------------------------------------
#                       SSP(4,2) explicit integrator                           #
#-------------------------------------------------------------------------------
def test_whether_ssp42_integrator_reproduces_expected_output(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('integrator_ssp42.in', tmp_path, stella_version)

    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The SSP(4,2) explicit integrator reproduces the expected potential.')
    return

#-------------------------------------------------------------------------------
#                       SSP(5,2) explicit integrator                           #
#-------------------------------------------------------------------------------
def test_whether_ssp52_integrator_reproduces_expected_output(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('integrator_ssp52.in', tmp_path, stella_version)

    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The SSP(5,2) explicit integrator reproduces the expected potential.')
    return

#-------------------------------------------------------------------------------
#                       SSP(4,3) explicit integrator                           #
#-------------------------------------------------------------------------------
def test_whether_ssp43_integrator_reproduces_expected_output(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('integrator_ssp43.in', tmp_path, stella_version)

    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The SSP(4,3) explicit integrator reproduces the expected potential.')
    return

#-------------------------------------------------------------------------------
#                       SSP(9,3) explicit integrator                           #
#-------------------------------------------------------------------------------
def test_whether_ssp93_integrator_reproduces_expected_output(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('integrator_ssp93.in', tmp_path, stella_version)

    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The SSP(9,3) explicit integrator reproduces the expected potential.')
    return

#-------------------------------------------------------------------------------
#                         RK6-4 explicit integrator                            #
#-------------------------------------------------------------------------------
def test_whether_rk64_integrator_reproduces_expected_output(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('integrator_rk64.in', tmp_path, stella_version)

    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The RK6-4 explicit integrator reproduces the expected potential.')
    return

#-------------------------------------------------------------------------------
#                   'fourierblob' initialization option                        #
#-------------------------------------------------------------------------------
def test_whether_init_fourierblob_option_is_the_same(tmp_path, stella_version, error=False):

    # Input file
    input_filename = 'init_fourierblob.in'

    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path, stella_version)

    # File names
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc')
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'

    # Check whether the potential at t=0 matches the stored file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:

        # Read the potential axis
        local_phi2 = local_netcdf['phi2']
        expected_phi2 = expected_netcdf['phi2']

        compare_local_array_with_expected_array(local_phi2, expected_phi2)
        if not (np.allclose(local_phi2, expected_phi2, equal_nan=True)):
            print('\nERROR: The potential data does not match in the netCDF files.'); error = True
            print('       It is likely that something changed in the fourierblob initialization of <g> of <phi>.')
            print('\nCompare the potential arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_phi2, expected_phi2)
        assert (not error), f'The potential data does not match in the netCDF files.'

    print('  -->  The potential is initialized the same as in the previous run (fourierblob).')
    return
