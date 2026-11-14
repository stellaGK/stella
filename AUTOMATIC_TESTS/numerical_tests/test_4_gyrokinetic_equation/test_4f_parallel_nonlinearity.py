################################################################################
#                       Check the parallel nonlinearity                        #
################################################################################
# The parallel nonlinearity (include_parallel_nonlinearity = .true.) is the
# acceleration along the field line by the fluctuating potential, proportional
# to rhostar, so it is only tested with rhostar > 0 (it vanishes for the default
# rhostar = 0). The input file parallel_nonlinearity.in has only the ExB and the
# parallel nonlinearity, with rhostar = 0.05 and adiabatic electrons, and starts
# from noise that is even in vpa.
#
# Without the linear terms, the parallel nonlinearity is the only term that can
# create a parallel flow from an initial distribution that is even in vpa, and
# the parallel flow is the moment it changes the most. Therefore:
#     - the full state, and the parallel flow relative to its own size, are
#       compared with the expected output (the parallel flow is much smaller than
#       the density, so the common tolerance of the moments would not see it);
#     - with rhostar = 0, the parallel flow has to stay zero up to round-off.
# Both bugs fixed in November 2026 (the advection speed taken at vpa = -vpa_max
# for every vpa, and the vpa-sign dependent d<phi>/dz at the ends of the chains)
# change the parallel flow of this test by 960% and 6%.
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

# Input file of the tests
input_filename = 'parallel_nonlinearity.in'

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                    Compare with the expected output                          #
#-------------------------------------------------------------------------------
def test_whether_parallel_nonlinearity_evolves_correctly(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation(input_filename, tmp_path, stella_version)

    # Compare phi2(t) and the full state in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)

    # Compare the parallel flow relative to its own size
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc')
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'
    compare_netcdf_quantities_normwise(local_netcdf_file, expected_netcdf_file, ['upar'], rtol=1e-6, atol=0.0)
    print('  -->  The parallel nonlinearity is evolving correctly.')
    return

#-------------------------------------------------------------------------------
#             Only the parallel nonlinearity drives a parallel flow            #
#-------------------------------------------------------------------------------
def test_whether_parallel_nonlinearity_drives_the_parallel_flow(tmp_path, stella_version):

    # The input file with rhostar = 0.05 and a copy with rhostar = 0
    upar_norm = {}
    for rhostar in ['0.05', '0.0']:
        run_path = tmp_path / f'rhostar{rhostar}'
        run_path.mkdir()
        path_input_file = copy_input_file(input_filename, run_path)
        text = path_input_file.read_text()
        assert '  rhostar = 0.05\n' in text, 'The input file should set rhostar = 0.05.'
        path_input_file.write_text(text.replace('  rhostar = 0.05\n', f'  rhostar = {rhostar}\n'))
        run_local_input_file(path_input_file, stella_version)
        with xr.open_dataset(run_path / input_filename.replace('.in', '.out.nc')) as netcdf:
            upar_norm[rhostar] = float(np.abs(netcdf['upar'].values).max())
            phi_norm = float(np.abs(netcdf['phi_vs_t'].values).max())

    # With rhostar = 0 the parallel flow stays zero up to round-off; with rhostar = 0.05 it grows
    print(f'\n    max|upar| = {upar_norm["0.05"]:.3e} (rhostar = 0.05), {upar_norm["0.0"]:.3e} (rhostar = 0), max|phi| = {phi_norm:.3e}')
    assert upar_norm['0.05'] > 1e-12 * phi_norm, 'The parallel nonlinearity does not drive a parallel flow.'
    assert upar_norm['0.0'] < 1e-6 * upar_norm['0.05'], 'A parallel flow appears without the parallel nonlinearity.'
    print('  -->  Only the parallel nonlinearity drives a parallel flow.')
    return

#-------------------------------------------------------------------------------
def run_local_input_file(path_input_file, stella_version):
    ''' Run stella with the input file <path_input_file> inside its own folder (it may differ from the
    input file of the test folder, which run_local_stella_simulation() would copy). '''
    if stella_version!='master':
        pytest.skip('The parallel nonlinearity test only runs on the current stella version.')
    os.chdir(path_input_file.parent)
    run_stella(get_stella_path(stella_version), path_input_file.name)
    return
