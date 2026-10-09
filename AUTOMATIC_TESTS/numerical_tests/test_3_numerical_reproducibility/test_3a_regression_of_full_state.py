################################################################################
#           Check that the full state of the simulation has not changed        #
################################################################################
# Most of the other numerical tests only compare |phi|^2(t), which is a single
# number per time step and which is insensitive to e.g. the phase of phi, to a
# permutation of the (kx,ky) modes, or to errors in the distribution function
# which have not (yet) affected the potential. Here we compare the fields on the
# full (t, tube, z, kx, ky) grid, the moments, the fluxes, the frequency and the
# distribution functions on the (z, vpa, mu) grid, with a tight tolerance.
#
# Each input file switches on a different numerical scheme (see the top of each
# input file), so that optimisations of the memory usage or speed of stella can
# be checked against all of them. The simulations are run on a single processor,
# the dependence on the number of processors is tested in test_3b.
#
# If a change of the numerics is intended, recreate the expected output with:
#     python3 create_expected_output.py
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

# Settings of the reproducibility tests
module_path = str(pathlib.Path(__file__).parent / 'reproducibility_settings.py')
with open(module_path, 'r') as file: exec(file.read())

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                 Compare the full state with the expected output              #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('input_file', input_files)
def test_whether_full_state_matches_expected_output(input_file, tmp_path, stella_version):

    # These input files only exist for the current stella version
    if stella_version != 'master': pytest.skip('Only implemented for the master branch of stella.')

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation(input_file, tmp_path, stella_version, nproc=nproc_expected_output)
    local_netcdf_file = tmp_path / input_file.replace('.in', '.out.nc')
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file.replace(".in","")}.out.nc'

    # Only compare the keys which are present in the expected output
    with xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        keys = [key for key in regression_keys if key in expected_netcdf.variables]

    # Compare the full arrays
    compare_netcdf_quantities_normwise(local_netcdf_file, expected_netcdf_file, keys, rtol=rtol, atol=atol, label='expected')
    print(f'  -->  The full state of {input_file} matches the expected output.')
    return
