################################################################################
#          Quickly check the code paths of the KBM and TAE benchmarks          #
################################################################################
# The KBM and TAE benchmarks (test_8b and test_8c) are slow, so they are skipped
# by default. They test a few options which are not covered by the other tests:
#     - KBM: a single (kx,ky) mode with grid_option = 'range' and nperiod = 2 in
#       electromagnetic stella, with zero upwinding.
#     - TAE: three species (including an energetic ion species), electromagnetic
#       stella with apar only (include_bpar = .false.), betaprim, scale_to_phiinit,
#       and stopping the simulation at <tend> instead of <nstep>.
# Here we run the same input files for only a few time steps, and compare the
# full fields, moments, fluxes and frequency with the expected
# output. These short runs do not test whether the KBM and TAE physics is
# captured, run the benchmarks for that with:
#     make numerical-tests-slow
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

# Quantities which are compared on the full grid
keys = ['t', 'phi_vs_t', 'apar_vs_t', 'bpar_vs_t', 'omega', 'density', 'upar', 'temperature',
        'pflux_vs_kxkys', 'qflux_vs_kxkys']

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                      Short runs of the KBM and TAE inputs                    #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('input_file, check_bpar', [('EM_KBM_short.in', True), ('EM_TAE_short.in', False)])
def test_short_KBM_and_TAE_runs(input_file, check_bpar, tmp_path, stella_version):

    # These input files only exist for the current stella version
    if stella_version != 'master': pytest.skip('Only implemented for the master branch of stella.')

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation(input_file, tmp_path, stella_version)
    local_netcdf_file = tmp_path / input_file.replace('.in', '.out.nc')
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file.replace(".in","")}.out.nc'

    # Compare the time traces of |phi|^2, |apar|^2 and |bpar|^2
    compare_local_potential_with_expected_potential_em(local_netcdf_file, expected_netcdf_file, check_bpar=check_bpar)

    # Compare the full arrays, for the keys which are present in the expected output
    with xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        keys_to_compare = [key for key in keys if key in expected_netcdf.variables]
    compare_netcdf_quantities_normwise(local_netcdf_file, expected_netcdf_file, keys_to_compare, rtol=1e-8)
    print(f'  -->  The short {input_file.split("_")[1]} run matches the expected output.')
    return
