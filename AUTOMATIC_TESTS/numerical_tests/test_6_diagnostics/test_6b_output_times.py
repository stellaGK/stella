################################################################################
#              Check the output times set by <dt_write> (issue #182)           #
################################################################################
# With <dt_write> > 0, stella writes the diagnostics to the text files at the
# first time step after every multiple n*<dt_write> of the time, and to the NetCDF
# file after every multiple of <dt_write>*<nc_mult>, instead of every <nwrite> time
# steps. Every multiple is written exactly once, also when the time step changes.
#   - Run the linear test of test 6a (delt = 0.01, nstep = 100) with dt_write = 0.034,
#     which is not a multiple of the time step, and nc_mult = 2.
#   - Run the nonlinear test of test 9 with delt = 0.02 and dt_write = 0.05, and restart
#     it with delt = 0.007, so that the time step changes during the simulation.
################################################################################

# Python modules
import pytest
import os, re
import pathlib
import subprocess
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
def latest_multiple(time, dt_write):
    '''The index n of the latest multiple n*dt_write which <time> has reached, as in stella.'''
    return np.floor(np.asarray(time) / dt_write + 1e-4).astype(int)

#-------------------------------------------------------------------------------
def times_after_multiples(times, dt_write, nc_mult=1):
    '''The times of the time steps at which stella writes to the text files (nc_mult = 1) or to the
    NetCDF file, i.e. the first time step after every multiple of dt_write (or dt_write*nc_mult).'''
    written, last = [], 0
    for time in times[1:]:
        n = latest_multiple(time, dt_write)
        if n // nc_mult > last // nc_mult: written.append(time)
        last = max(n, last)
    return np.array(written)

#-------------------------------------------------------------------------------
def test_whether_diagnostics_are_written_every_dt_write(tmp_path, stella_version):
    if stella_version!='master': pytest.skip('Older stella versions do not have <dt_write>.')
    input_filename = 'miller_linear.in'
    delt, nstep, dt_write, nc_mult = 0.01, 100, 0.034, 2

    # Write the input file with dt_write and nc_mult in the <diagnostics> namelist
    path_input_file = copy_input_file(input_filename, tmp_path)
    text = path_input_file.read_text()
    assert f'delt = {delt}' in text and f'nstep = {nstep}' in text and '  nwrite = 10\n' in text
    assert 'dt_write' not in text and 'nc_mult' not in text
    path_input_file.write_text(text.replace('  nwrite = 10\n', f'  nwrite = 10\n  dt_write = {dt_write}\n  nc_mult = {nc_mult}\n'))

    # Run stella, and keep its output to check the message about the number of time slices
    os.chdir(tmp_path)
    result = subprocess.run(['mpirun', '--oversubscribe', '-np', f'{read_nproc()}', get_stella_path(stella_version), input_filename],
        check=True, capture_output=True, text=True)
    assert 'The diagnostics are written every dt_write' in result.stdout, 'stella did not report the output interval.'

    # The time of every time step (the time step stays at <delt> in this linear run)
    times = np.insert(np.cumsum(np.full(nstep, delt)), 0, 0.0)

    # The NetCDF file contains t = 0 and the first time step after every multiple of dt_write*nc_mult
    with xr.open_dataset(tmp_path / input_filename.replace('.in', '.out.nc')) as data:
        netcdf_times = data['t'].values
    expected_netcdf_times = np.insert(times_after_multiples(times, dt_write, nc_mult), 0, 0.0)
    assert netcdf_times.shape == expected_netcdf_times.shape and np.allclose(netcdf_times, expected_netcdf_times, rtol=0, atol=1e-10), \
        f'The NetCDF file was written at t = {netcdf_times}, instead of {expected_netcdf_times}.'

    # The text files are written at the first time step after every multiple of dt_write (the omega file skips t = 0)
    omega_times = np.unique(np.loadtxt(tmp_path / input_filename.replace('.in', '.omega'))[:, 0])
    expected_text_times = times_after_multiples(times, dt_write)
    assert omega_times.shape == expected_text_times.shape and np.allclose(omega_times, expected_text_times, rtol=0, atol=1e-7), \
        f'The omega file was written at t = {omega_times}, instead of {expected_text_times}.'
    print(f'\n  -->  The diagnostics are written at the first time step after every multiple of dt_write = {dt_write}.')
    return

#-------------------------------------------------------------------------------
def test_whether_every_multiple_of_dt_write_is_written_once_when_the_time_step_changes(tmp_path, stella_version):
    if stella_version!='master': pytest.skip('Older stella versions do not have <dt_write>.')
    input_filename = 'input.in'
    dt_write, delt_before, delt_after = 0.05, 0.02, 0.007

    # Read the input file of the restart test (test 9), which saves the restart files every 10 time steps
    text = (get_automatic_tests_directory() / 'numerical_tests' / 'test_9_restarting' / input_filename).read_text()
    assert 'nstep = 200' in text and '  nwrite = 10\n' in text and 'dt_write' not in text
    text = text.replace('  nwrite = 10\n', f'  nwrite = 10\n  dt_write = {dt_write}\n')
    def write_input_file(nstep, delt, restart):
        new_text = text.replace('nstep = 200', f'nstep = {nstep}')
        new_text = re.sub(r"&time_step_temp\n.*?\n/\n", f"&time_step\n  delt_option = 'set_by_hand'\n  delt = {delt}\n/\n", new_text, flags=re.S)
        if restart: new_text = new_text.replace('&initialise_distribution_temp', '&initialise_distribution')
        else: new_text = re.sub(r'&initialise_distribution_temp\n.*?\n/\n', '', new_text, flags=re.S)
        (tmp_path / input_filename).write_text(new_text)

    # Run 100 time steps with delt = 0.02, and restart for another 100 time steps with delt = 0.007
    # (delt_option = 'set_by_hand', since by default a restarted simulation continues with the saved time step)
    os.chdir(tmp_path)
    write_input_file(nstep=100, delt=delt_before, restart=False); run_stella(get_stella_path(stella_version), input_filename)
    with xr.open_dataset(tmp_path / input_filename.replace('.in', '.out.nc')) as data: time_restart = float(data['t'][-1])
    write_input_file(nstep=200, delt=delt_after, restart=True); run_stella(get_stella_path(stella_version), input_filename)
    with xr.open_dataset(tmp_path / input_filename.replace('.in', '.out.nc')) as data: times = data['t'].values
    assert times[-1] > time_restart + 10 * delt_after, 'The simulation was not restarted.'

    # Every multiple n*dt_write is written exactly once (dt_write is larger than the time step), within a time step after it
    n = latest_multiple(times, dt_write)
    assert np.array_equal(n, np.arange(len(times))), f'The multiples n of dt_write that were written are {n}, instead of 0, 1, 2, ...'
    delt = np.where(times <= time_restart + 1e-10, delt_before, delt_after)
    assert np.all(times - n * dt_write >= -1e-4 * dt_write) and np.all(times - n * dt_write < delt + 1e-10), \
        f'The diagnostics were not written at the first time step after every multiple of dt_write: t = {times}'
    print(f'\n  -->  Every multiple of dt_write = {dt_write} is written once, before ({delt_before}) and after ({delt_after}) the restart.')
    return
