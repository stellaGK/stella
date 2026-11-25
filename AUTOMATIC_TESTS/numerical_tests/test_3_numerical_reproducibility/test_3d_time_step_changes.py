################################################################################
#     Check that every numerical scheme follows a change of the time step      #
################################################################################
# Several parts of stella store matrices which contain the time step code_dt,
# e.g. the response matrix of the implicit parallel streaming, the implicit
# mirror and drift terms, and the implicit Dougherty and Fokker-Planck collision
# operators. When the CFL condition changes code_dt, reset_dt has to rebuild all
# of them; a matrix which is not rebuilt keeps advancing its term with the old
# time step (e.g. the implicit collision frequency was effectively multiplied by
# dt_old/dt_new until this was fixed by R. Nies in October 2026).
#
# Here we run every input file of test 3 twice, as a linear simulation:
#     (A) with delt = 10, which stella reduces to the CFL limit before the first
#         time step, so that reset_dt rebuilds everything for the new time step
#     (B) with delt equal to the reduced time step of run (A), so that the
#         matrices are only built once, directly with this time step
# Both runs advance with the same time step, so they should give the same
# results. A term which does not follow the time step change gives differences
# of order 1. The time step changes during nonlinear simulations go through the
# same reset_dt routine, so they are covered as well.
#
# When a new term with a time step dependent matrix is added to stella, add an
# input file which switches it on to <input_files> in reproducibility_settings.py.
# The relevant parameters are:
#       &time_step
#         delt = 10.0
#       /
################################################################################

# Python modules
import pytest
import os, sys
import re
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
# results agree to about 1e-8; a term which keeps the old time step gives differences of order 1
rtol_time_step = 1e-5

# Quantities which are compared on the full grid
keys_time_step = ['phi_vs_t', 'apar_vs_t', 'bpar_vs_t', 'density', 'upar', 'temperature', 'g2_vs_zvpamus']

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
def write_linear_input_file(input_file, folder, delt):
    '''Copy <input_file> to <folder> as a linear simulation with time step <delt>. The
    nonlinear term is switched off, since its CFL condition would change the time step
    again during the simulation, at slightly different times in the runs (A) and (B).'''
    path_input_file = copy_input_file(input_file, folder)
    text = open(path_input_file, 'r').read()
    for old, new in (('include_nonlinear = .true.', 'include_nonlinear = .false.'),):
        assert old in text, f'Expected "{old}" in {input_file}.'
        text = text.replace(old, new)
    # Replace the time step of the input file, whatever its value
    text, count = re.subn(r'^(\s*delt\s*=\s*)\S+', lambda m: f'{m.group(1)}{delt!r}', text, flags=re.M)
    assert count == 1, f'Expected one "delt = ..." line in {input_file}.'
    open(path_input_file, 'w').write(text)
    return path_input_file

#-------------------------------------------------------------------------------
def run_linear_simulation(input_file, folder, stella_version, delt):
    '''Run <input_file> as a linear simulation with time step <delt> inside <folder>.'''
    folder.mkdir()
    write_linear_input_file(input_file, folder, delt)
    os.chdir(folder)
    run_stella(get_stella_path(stella_version), input_file, nproc=nproc_expected_output)
    return folder / input_file.replace('.in', '.out.nc')

#-------------------------------------------------------------------------------
@pytest.mark.parametrize('input_file', input_files)
def test_whether_the_time_step_change_is_followed_by_every_term(input_file, tmp_path, stella_version):

    # These input files only exist for the current stella version
    if stella_version != 'master': pytest.skip('Only implemented for the master branch of stella.')

    # (A) delt = 10, which stella reduces to the CFL limit before the first time step
    netcdf_a = run_linear_simulation(input_file, tmp_path / 'delt_reduced_by_cfl', stella_version, 10.0)

    # The reduced time step: t is written every nwrite = 10 time steps of a constant code_dt
    with xr.open_dataset(netcdf_a) as netcdf:
        code_dt = float(netcdf['t'].values[1] - netcdf['t'].values[0]) / 10
    assert code_dt < 9, f'The time step was not reduced by the CFL condition (code_dt = {code_dt}).'

    # (B) delt = the reduced time step (slightly smaller, so that stella does not reduce it again)
    netcdf_b = run_linear_simulation(input_file, tmp_path / 'delt_cfl', stella_version, code_dt * (1 - 1e-9))

    # Compare the time evolution
    with xr.open_dataset(netcdf_a) as reduced, xr.open_dataset(netcdf_b) as direct:
        assert np.allclose(reduced['t'], direct['t'], rtol=1e-6), 'The time axes do not match.'
        failed = []
        for key in keys_time_step:
            if key not in reduced.variables: continue
            scale = np.max(np.abs(direct[key].values))
            difference = np.max(np.abs(reduced[key].values - direct[key].values))
            assert np.isfinite(scale) and np.isfinite(difference), f'{key} is not finite: the simulation is unstable.'
            print(f'    {key:<18s} max|diff|/max = {difference/scale:.3e}')
            if not difference <= rtol_time_step * scale + atol: failed.append(key)
    assert not failed, f'{input_file} does not follow the time step change of the CFL condition: {failed} differ.'
    print(f'  -->  All terms of {input_file} follow the time step change (code_dt = {code_dt:.4g}).')
    return
