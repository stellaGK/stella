################################################################################
#              Check whether pyrokinetics can plot stella output               #
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

# Pyrokinetics package
from pyrokinetics import Pyro

# Global variables  
input_filename = 'miller_linear_CBC.in'
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
def test_whether_we_can_load_pyrokinetics(tmp_path, stella_version):
    '''Run a local stella simulation in a temporary folder <tmp_path>.'''  
    
    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global local_stella_run_directory
    local_stella_run_directory = tmp_path
    
    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    
    # Load the pyro data
    local_input_file = local_stella_run_directory / input_filename
    pyro = Pyro(gk_file=local_input_file, gk_code="STELLA")
    print('\n  -->  Successfully loaded pyrokinetics.')
    return
