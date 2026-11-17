################################################################################
#       Check that writing diagnostics does not change the time evolution      #
################################################################################
# The diagnostics are a typical target when reducing the memory usage of stella
# (e.g. by reusing work arrays). They should never modify the distribution
# function or the fields. Here we run the same simulation writing the
# diagnostics every 10 time steps (nwrite = 10) and every time step
# (nwrite = 1), and check that the results agree at the common time steps.
# The relevant parameters are:
#       &diagnostics
#         nwrite = 10
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
@pytest.mark.parametrize('input_file', ['es_nonlinear.in', 'em_nonlinear.in'])
def test_whether_writing_diagnostics_changes_the_time_evolution(input_file, tmp_path, stella_version):

    # These input files only exist for the current stella version
    if stella_version != 'master': pytest.skip('Only implemented for the master branch of stella.')

    # Run stella with nwrite = 10 (as in the input file) and nwrite = 1
    netcdf_files = {}
    for nwrite in [10, 1]:
        folder = tmp_path / f'nwrite{nwrite}'; folder.mkdir()
        path_input_file = copy_input_file(input_file, folder)
        text = open(path_input_file, 'r').read()
        assert 'nwrite = 10' in text, f'Expected "nwrite = 10" in {input_file}.'
        open(path_input_file, 'w').write(text.replace('nwrite = 10', f'nwrite = {nwrite}'))
        os.chdir(folder)
        run_stella(get_stella_path(stella_version), input_file, nproc=nproc_expected_output)
        netcdf_files[nwrite] = folder / input_file.replace('.in', '.out.nc')

    # Compare the fields at the common time steps
    with xr.open_dataset(netcdf_files[10]) as sparse, xr.open_dataset(netcdf_files[1]) as dense:
        dense = dense.isel(t=slice(None, None, 10))
        assert np.allclose(sparse['t'], dense['t'], rtol=1e-12), 'The time axes do not match.'
        failed = []
        for key in ['phi_vs_t', 'apar_vs_t', 'bpar_vs_t', 'density', 'upar', 'temperature', 'g2_vs_zvpamus']:
            if key not in sparse.variables: continue
            scale = np.max(np.abs(sparse[key].values))
            difference = np.max(np.abs(sparse[key].values - dense[key].values))
            print(f'    {key:<18s} max|diff|/max = {difference/scale:.3e}')
            if difference > rtol * scale + atol: failed.append(key)
    assert not failed, f'Writing the diagnostics more often changes {failed}.'
    print(f'  -->  Writing the diagnostics does not change the time evolution of {input_file}.')
    return
