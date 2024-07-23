
################################################################################
#               Check all the terms in the gyrokinetic equation                #
################################################################################
# These test are designed to test the terms in the gyrokinetic equation one 
# by one. As a sanity check we start by turning of all terms to check that
# the potential is not being evolved in time. To ensure this we need:
#       &physics_flags 
#         nonlinear = .false.
#         include_parallel_streaming = .false.
#         include_mirror = .false.
#       /
#       &time_advance_knobs
#         xdriftknob = 0
#         ydriftknob = 0
#         wstarknob = 0 
#       /
#       &dissipation
#         hyper_dissipation = .false.
#       /   
################################################################################

# Python modules
import os, sys
import pathlib 
import numpy as np
import xarray as xr  

# Package to run stella 
module_path = str(pathlib.Path(__file__).parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Global variables  
input_filename = 'no_time_evolution.in'  
stella_local_run_directory = 'Not/Run/Yet' 

#-------------------------------------------------------------------------------
#           Check whether the potential data does not evolve in time           #
#-------------------------------------------------------------------------------
def test_whether_potential_data_in_netcdf_file_remains_constant(tmp_path):
    
    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global stella_local_run_directory
    stella_local_run_directory = tmp_path
    
    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path)
     
    # File names  
    local_netcdf_file = stella_local_run_directory / input_filename.replace('.in','.out.nc') 

    # Check whether the geometry data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf:
        
        # Read the potential axis
        local_phi2 = local_netcdf['phi2']
        
        # Check whether the potential did not evolve in time
        if not (np.all(local_phi2 == local_phi2[0])):
            print(f'ERROR: The potential is evolving in time, while it should not.')
            for i in range(len(local_phi2)):
                print(f'phi2 = {float(local_phi2.data[i]):.5e}')
            assert False, f'The potential is evolving in time, while it should not.'
                
    print(f'  -->  Without gyrokinetic terms the potential does not evolve in time ({int(local_netcdf["nproc"])} CPUs).')
    return
