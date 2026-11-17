################################################################################
#               Check whether simulation is restarted correctly                #
################################################################################
# Run a nonlinear simulation for 200 time steps, and run the same simulation for
# 100 time steps, followed by a restart for another 100 time steps. The restarted
# simulation should reproduce the continuous simulation, hence we compare the full
# state (the fields on the full grid, the moments and the distribution functions)
# after the restart. This is tested both when each processor writes its own
# restart file (save_many = .true.) and when a single restart file is written.
#
# The restart is controlled by renaming two namelists in the input file:
#       &time_step_temp                    ->  &time_step
#         delt_option = 'check_restart'          (read the time step from the restart file)
#       &initialise_distribution_temp      ->  &initialise_distribution
#         initialise_distribution_option = 'many'  (read g from the restart files)
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

# Quantities which are compared after the restart
keys = ['phi2', 'phi_vs_t', 'density', 'upar', 'temperature', 'g2_vs_zvpamus', 'h2_vs_zvpamus']

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
def write_input_file(folder, input_filename, nstep, restart):
    '''Copy <input_filename> to <folder>, with <nstep> time steps, and switch on the
    restart namelists if <restart> is True. Returns the name of the input file.'''
    text = open(get_stella_expected_run_directory() / input_filename, 'r').read()
    assert 'nstep = 200' in text, f'Expected "nstep = 200" in {input_filename}.'
    text = text.replace('nstep = 200', f'nstep = {nstep}')
    if restart:
        assert '&time_step_temp' in text and '&initialise_distribution_temp' in text
        text = text.replace('&time_step_temp', '&time_step')
        text = text.replace('&initialise_distribution_temp', '&initialise_distribution')
    folder.mkdir(parents=True, exist_ok=True)
    open(folder / input_filename, 'w').write(text)
    return input_filename

#-------------------------------------------------------------------------------
def first_time_step_printed_by_stella(stella_path, input_filename, nproc=None):
    '''Run stella and return the first time step it prints to the command prompt. A 
    restarted simulation continues from the time step at which it was stopped, while
    a simulation that (silently) starts from scratch would print time step 0.'''
    if not nproc: nproc = read_nproc()
    result = subprocess.run(['mpirun', '--oversubscribe', '-np', f'{nproc}', stella_path, input_filename], 
                            check=True, capture_output=True, text=True)
    print(result.stdout)
    for line in result.stdout.split('\n'):
        values = line.split()
        if len(values) == 4 and values[0].isdigit(): return int(values[0])
    pytest.fail('Could not find the time steps in the output of stella.')

#-------------------------------------------------------------------------------
#                Compare the restarted and the continuous simulation           #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('input_filename', ['input.in', 'input_singlerestartfile.in'])
def test_whether_restarted_simulation_matches_continuous_simulation(input_filename, tmp_path, stella_version):

    # The restart namelists only exist for the current stella version
    if stella_version != 'master': pytest.skip('Only implemented for the master branch of stella.')
    stella_path = get_stella_path(stella_version)
    netcdf_filename = input_filename.replace('.in', '.out.nc')

    # Run the continuous simulation for 200 time steps
    continuous_folder = tmp_path / 'continuous'
    write_input_file(continuous_folder, input_filename, nstep=200, restart=False)
    os.chdir(continuous_folder); run_stella(stella_path, input_filename)

    # Run 100 time steps, and restart the simulation for another 100 time steps
    restarted_folder = tmp_path / 'restarted'
    write_input_file(restarted_folder, input_filename, nstep=100, restart=False)
    os.chdir(restarted_folder); run_stella(stella_path, input_filename)
    with xr.open_dataset(restarted_folder / netcdf_filename) as netcdf: time_restart = float(netcdf['t'][-1])
    write_input_file(restarted_folder, input_filename, nstep=200, restart=True)
    os.chdir(restarted_folder); first_time_step = first_time_step_printed_by_stella(stella_path, input_filename)
    assert first_time_step > 0, 'The simulation was not restarted, it started again at time step 0.'

    # Compare the time steps after the restart (the netcdf file of the restarted
    # simulation also contains the time steps from before the restart)
    failed = []
    with xr.open_dataset(continuous_folder / netcdf_filename) as continuous, xr.open_dataset(restarted_folder / netcdf_filename) as restarted:
        assert_same_shape(restarted['t'], continuous['t'], 't')
        assert np.allclose(restarted['t'], continuous['t'], rtol=1e-12), 'The time axes do not match.'
        after_restart = (continuous['t'].values > time_restart + 1e-10)
        assert after_restart.sum() > 0, 'The netcdf file contains no time steps after the restart.'
        print(f'\n    Compare {after_restart.sum()} time steps after the restart at t = {time_restart}:')
        for key in keys:
            a = continuous[key].values[after_restart]; b = restarted[key].values[after_restart]
            relative_difference = np.max(np.abs(b - a)) / np.max(np.abs(a))
            print(f'    {key:<18s} max|diff|/max = {relative_difference:.3e}')
            if not relative_difference <= 1e-8: failed.append(key)
    assert not failed, f'The restarted simulation does not match the continuous simulation for {failed}.'
    print(f'  -->  The restarted simulation matches the continuous simulation ({input_filename}).')
    return
