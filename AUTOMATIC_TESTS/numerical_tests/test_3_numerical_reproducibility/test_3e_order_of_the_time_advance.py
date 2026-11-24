################################################################################
#          Check the order of accuracy of the time advance of stella           #
################################################################################
# A bug in the time advance often keeps the results close to the correct ones,
# but lowers the order of accuracy of the scheme, e.g. a term which is evaluated
# at the wrong Runge-Kutta stage, or a boundary condition which is only imposed
# at the end of a time step. Such bugs are invisible in a comparison with the
# expected output, but they show up in the convergence with the time step.
#
# Here we run a short linear simulation up to the same time with the time steps
# dt, dt/2 and dt/4, and estimate the order p of the time advance from
#     |g(dt) - g(dt/2)| / |g(dt/2) - g(dt/4)| = 2^p
# for the potential and the distribution function at the final time. We check
#     - the explicit algorithms rk2, rk3 and rk4 with the parallel streaming and
#       mirror terms treated explicitly, which should converge at order 2, 3, 4
#     - the default algorithm (rk3 with implicit parallel streaming and mirror),
#       whose operator splitting should converge at least at order 1
# The relevant parameters are:
#       &numerical_algorithms
#         explicit_algorithm = 'rk3'
#         stream_implicit = .false.
#         mirror_implicit = .false.
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

# The measured order may be slightly below the order of the scheme, since dt is not
# infinitely small (e.g. rk4 gives 3.75 for g, while rk2 and rk3 give 2.0-2.1 and 2.9-3.0)
tolerance_order = 0.35

# (input file, explicit algorithm, largest time step, number of time steps, expected order)
# The time steps of the fully explicit runs are limited by the stability of the explicit
# parallel streaming of the kinetic electrons, well below the CFL limit of 6.8e-3.
time_advance_schemes = [
    ('stream_and_mirror_explicit.in', 'rk2', 0.0008, 25, 2),
    ('stream_and_mirror_explicit.in', 'rk3', 0.0008, 25, 3),
    ('stream_and_mirror_explicit.in', 'rk4', 0.0008, 25, 4),
    ('es_nonlinear.in', 'rk3', 0.01, 10, 1),
]

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
def run_linear_simulation(input_file, folder, stella_version, explicit_algorithm, delt, nstep):
    '''Run <input_file> inside <folder> as a linear simulation with the given explicit
    algorithm, time step and number of time steps, writing only the initial and final state.'''
    folder.mkdir()
    path_input_file = copy_input_file(input_file, folder)
    text = open(path_input_file, 'r').read()
    for old, new in (('include_nonlinear = .true.', 'include_nonlinear = .false.'),
                     ('nstep = 50', f'nstep = {nstep}'),
                     ('nwrite = 10', f'nwrite = {nstep}'),
                     ("explicit_algorithm = 'rk3'", f"explicit_algorithm = '{explicit_algorithm}'")):
        assert old in text, f'Expected "{old}" in {input_file}.'
        text = text.replace(old, new)
    # Replace the time step of the input file, whatever its value
    text, count = re.subn(r'^(\s*delt\s*=\s*)\S+', lambda m: f'{m.group(1)}{delt!r}', text, flags=re.M)
    assert count == 1, f'Expected one "delt = ..." line in {input_file}.'
    open(path_input_file, 'w').write(text)
    os.chdir(folder)
    run_stella(get_stella_path(stella_version), input_file, nproc=nproc_expected_output)
    with xr.open_dataset(folder / input_file.replace('.in', '.out.nc')) as netcdf:
        assert np.isclose(netcdf['t'].values[-1], delt * nstep, rtol=1e-6), \
            f'The simulation did not reach t = {delt * nstep} with delt = {delt} (did the CFL condition change it?).'
        return {key: netcdf[key].values[-1].copy() for key in ['phi_vs_t', 'g2_vs_zvpamus']}

#-------------------------------------------------------------------------------
@pytest.mark.parametrize('input_file, explicit_algorithm, delt, nstep, expected_order', time_advance_schemes)
def test_order_of_the_time_advance(input_file, explicit_algorithm, delt, nstep, expected_order, tmp_path, stella_version):

    # These input files only exist for the current stella version
    if stella_version != 'master': pytest.skip('Only implemented for the master branch of stella.')

    # Run up to the same time t = delt * nstep with the time steps delt, delt/2 and delt/4
    states = [run_linear_simulation(input_file, tmp_path / f'dt{level}', stella_version,
              explicit_algorithm, delt / 2**level, nstep * 2**level) for level in range(3)]

    # Estimate the order from the differences between consecutive time steps
    failed = []
    for key in ['phi_vs_t', 'g2_vs_zvpamus']:
        scale = np.max(np.abs(states[-1][key]))
        difference_coarse = np.max(np.abs(states[0][key] - states[1][key])) / scale
        difference_fine = np.max(np.abs(states[1][key] - states[2][key])) / scale
        order = np.log2(difference_coarse / difference_fine)
        print(f'    {key:<14s} |dt - dt/2| = {difference_coarse:.3e}   |dt/2 - dt/4| = {difference_fine:.3e}   order = {order:.2f}')
        if not order >= expected_order - tolerance_order: failed.append(f'{key} (order {order:.2f})')
    assert not failed, f'The time advance of {input_file} with {explicit_algorithm} converges at a lower order than {expected_order}: {failed}.'
    print(f'  -->  The time advance of {input_file} with {explicit_algorithm} converges at order {expected_order} or higher.')
    return
