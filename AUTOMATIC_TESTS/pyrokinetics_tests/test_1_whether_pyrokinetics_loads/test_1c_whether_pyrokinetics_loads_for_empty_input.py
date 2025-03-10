################################################################################
#              Check whether pyrokinetics can plot stella output               #
################################################################################

# Python modules
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
input_filename = 'empty.in'  
local_stella_run_directory = 'Not/Run/Yet'

#-------------------------------------------------------------------------------
#                         Run local stella simulation                          #
#-------------------------------------------------------------------------------
def TODO_test_whether_we_can_load_pyrokinetics(tmp_path):
    '''Run a local stella simulation in a temporary folder <tmp_path>.'''  
    
    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global local_stella_run_directory
    local_stella_run_directory = tmp_path
    
    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path) 
    
    # Load the pyro data
    local_input_file = local_stella_run_directory / input_filename
    pyro = Pyro(gk_file=local_input_file, gk_code="STELLA")
    print('\n  -->  Successfully loaded pyrokinetics.')
    return 

