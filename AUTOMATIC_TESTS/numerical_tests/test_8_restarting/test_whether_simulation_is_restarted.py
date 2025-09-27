################################################################################
#               Check whether simulation is restarted correctly                #
################################################################################
# Test whether a stella simulation is restarted correctly.
################################################################################

# Python modules
import pytest
import shutil
import os, sys
import pathlib
import numpy as np
import xarray as xr

# Package to run stella 
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Global variables
input_filename = 'input.in'
local_stella_run_directory = 'Not/Run/Yet'

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                         Run local stella simulation                          #
#-------------------------------------------------------------------------------
def test_whether_we_can_run_a_local_stella_simulation(tmp_path, stella_version):
    '''Run a local stella simulation in a temporary folder <tmp_path>.'''

    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global local_stella_run_directory, input_filename
    local_stella_run_directory = tmp_path
    os.chdir(tmp_path)

    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    
    # Check the potential at the last time step (nstep = 200)
    output_file = tmp_path / input_filename.replace('.in','.out')
    data = np.loadtxt(output_file,skiprows=2,dtype='float').reshape(-1, 5)
    steps = data[:,0]; time = data[:,1]; potential = data[:,2]
    last_potential_point_test1 = potential[-1]
    print(f'\nRan a stella simulation for nstep = {int(steps[-1])} giving |phi|^2 = {potential[-1]} at t = {time[-1]}.\n')
    
    # Now create two input files, the first one runs upto nstep = 100
    input_filename_restart_simulation = input_filename.replace('.in','_restarted.in')
    shutil.copyfile(tmp_path / input_filename, tmp_path / input_filename_restart_simulation)
    input_file = tmp_path / input_filename_restart_simulation
    with open(input_file, 'r') as file: filedata = file.read()
    filedata = filedata.replace('nstep = 200', 'nstep = 100')
    with open(input_file, 'w') as file: file.write(filedata)
      
    # Run the first input file
    run_stella(get_stella_path(stella_version), input_filename_restart_simulation)
    
    # Now modify the input file to restart the simulation upto nstep = 200
    with open(input_file, 'r') as file: filedata = file.read()
    filedata = filedata.replace('nstep = 100', 'nstep = 200')
    filedata = filedata.replace('&time_step_temp', '&time_step')
    filedata = filedata.replace('&initialise_distribution_temp', '&initialise_distribution')
    with open(input_file, 'w') as file: file.write(filedata)
      
    # Restart the first simulation
    run_stella(get_stella_path(stella_version), input_filename_restart_simulation)
    
    # Check the potential at the last time step of the restarted simulation (nstep = 200)
    output_file = tmp_path / input_filename_restart_simulation.replace('.in','.out')
    data = np.loadtxt(output_file,skiprows=2,dtype='float').reshape(-1, 5)
    steps = data[:,0]; time = data[:,1]; potential = data[:,2]
    last_potential_point_test2 = potential[-1]
    print(f'\nRan a stella simulation for nstep = {int(steps[-1])} giving |phi|^2 = {potential[-1]} at t = {time[-1]}.\n')
    
    # Check whether it worked
    if (last_potential_point_test1!=last_potential_point_test2):
        print(f'ERROR: The simulation was not succesfully restarted.') 
        assert False, f'The simulation was not succesfully restarted.'
    else:
        print(f'The simulation was succesfully restarted.\n') 
    
    return
    
    
