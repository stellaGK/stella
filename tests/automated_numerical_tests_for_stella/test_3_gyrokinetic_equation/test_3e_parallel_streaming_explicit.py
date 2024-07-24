
################################################################################
#               Check all the terms in the gyrokinetic equation                #
################################################################################
# These test are designed to test the terms in the gyrokinetic equation one 
# by one. Here we test the implicit parallel streaming term:
#       &physics_flags 
#         nonlinear = .true.
#         include_parallel_streaming = .true.
#         include_mirror = .false.
#       /
#       &time_advance_knobs
#         xdriftknob = 0
#         ydriftknob = 0
#         wstarknob = 0 
#       /
#       &knobs 
#         stream_implicit = .false.
#         mirror_implicit = .false.
#         drifts_implicit = .false.
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
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Global variables  
input_filename = 'parallel_streaming_explicit.in'  

#-------------------------------------------------------------------------------
#     Check whether the explicit parallel streaming term evolves correctly     #
#-------------------------------------------------------------------------------
def test_whether_parallel_streaming_explicit_evolves_correctly(tmp_path, error=False):

    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   
    compare_local_potential_with_expected_potential(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  The parallel streaming term, implemented explicitly, is evolving correctly.')
    return
    
