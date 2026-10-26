################################################################################
#   Check that implicit collisions use the time step after a CFL time step change  #
################################################################################
# The matrices of the implicit collision algorithms contain the time step code_dt.
# When the user-specified time step is larger than the CFL condition allows, stella
# reduces code_dt and calls reset_dt, which has to rebuild these matrices; otherwise
# the collision frequency is effectively multiplied by delt/code_dt.
# Here we run the same linear simulation twice:
#     (A) with delt = 1.0, which stella reduces to the CFL limit at the start
#     (B) with delt equal to the reduced time step of run (A), so no reset is needed
# Both runs advance with the same time step, so they should give the same results.
# The relevant parameters are:
#       &dissipation_and_collisions_options
#         include_collisions = .true.
#         collisions_implicit = .true.
#       /
#       &time_step
#         delt = 1.0
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

# The runs (A) and (B) differ in their time step by a relative 1e-9 (see below), so the
# results agree to about 1e-8; a stale collision matrix gives differences of order 1 or larger
rtol_time_step = 1e-5

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
def write_linear_input_file(input_file, folder, delt):
    """ Copy <input_file> to <folder> as a linear simulation with a stronger collision frequency
    (vnew_ref = 0.1, so that the collisions clearly affect the evolution) and time step <delt>.
    Only the ions are kept (adiabatic electrons): with kinetic electrons the explicit terms are
    numerically unstable at the CFL-limited time step of the ions, even without collisions. """
    path_input_file = copy_input_file(input_file, folder)
    text = open(path_input_file, 'r').read()
    for old, new in (('include_nonlinear = .true.', 'include_nonlinear = .false.'),
                     ('nspec = 2', 'nspec = 1'),
                     ('vnew_ref = 0.01', 'vnew_ref = 0.1'),
                     ('delt = 0.01', f'delt = {delt!r}')):
        assert old in text, f'Expected "{old}" in {input_file}.'
        text = text.replace(old, new)
    open(path_input_file, 'w').write(text)
    return path_input_file

#-------------------------------------------------------------------------------
@pytest.mark.parametrize('input_file', ['collisions_dougherty_implicit.in', 'collisions_fokker_planck_implicit.in'])
def test_whether_implicit_collisions_follow_time_step_changes(input_file, tmp_path, stella_version):

    # These input files only exist for the current stella version
    if stella_version != 'master': pytest.skip('Only implemented for the master branch of stella.')

    # (A) delt = 1.0, which stella reduces to the CFL limit before the first time step
    folder_a = tmp_path / 'delt_reduced_by_cfl'; folder_a.mkdir()
    write_linear_input_file(input_file, folder_a, 1.0)
    os.chdir(folder_a)
    run_stella(get_stella_path(stella_version), input_file, nproc=nproc_expected_output)
    netcdf_a = folder_a / input_file.replace('.in', '.out.nc')

    # The reduced time step: t is written every nwrite = 10 time steps of a constant code_dt
    with xr.open_dataset(netcdf_a) as netcdf:
        code_dt = float(netcdf['t'].values[1] - netcdf['t'].values[0]) / 10
    assert code_dt < 0.9, f'The time step was not reduced by the CFL condition (code_dt = {code_dt}).'

    # (B) delt = the reduced time step (slightly smaller, so that stella does not reduce it again)
    folder_b = tmp_path / 'delt_cfl'; folder_b.mkdir()
    write_linear_input_file(input_file, folder_b, code_dt * (1 - 1e-9))
    os.chdir(folder_b)
    run_stella(get_stella_path(stella_version), input_file, nproc=nproc_expected_output)
    netcdf_b = folder_b / input_file.replace('.in', '.out.nc')

    # Compare the time evolution
    with xr.open_dataset(netcdf_a) as reduced, xr.open_dataset(netcdf_b) as direct:
        assert np.allclose(reduced['t'], direct['t'], rtol=1e-6), 'The time axes do not match.'
        failed = []
        for key in ['phi_vs_t', 'density', 'upar', 'temperature', 'g2_vs_zvpamus']:
            if key not in reduced.variables: continue
            scale = np.max(np.abs(direct[key].values))
            difference = np.max(np.abs(reduced[key].values - direct[key].values))
            assert np.isfinite(scale) and np.isfinite(difference), f'{key} is not finite: the simulation is unstable.'
            print(f'    {key:<18s} max|diff|/max = {difference/scale:.3e}')
            if not difference <= rtol_time_step * scale: failed.append(key)
    assert not failed, f'The implicit collisions do not follow the reduced time step: {failed} differ.'
    print(f'  -->  The implicit collisions of {input_file} use the reduced time step.')
    return
