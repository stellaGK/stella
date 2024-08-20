
################################################################################
#               Check all the terms in the gyrokinetic equation                #
################################################################################
# These test are designed to test the terms in the gyrokinetic equation one 
# by one for the case of Full Flux Surface with modified adiabatic electrons 
# (Boltzmann electrons). As a sanity check we start by turning of all terms to 
# check that the potential is not being evolved in time. Then we turn the terms 
# on one by one linearly, then all together linearly, and finally we test them 
# when every term is turned on nonlinearly. The flags we toggle are:
#       &physics_flags 
#         nonlinear = .false. / .true.
#         include_parallel_streaming = .false. / .true.
#         include_mirror = .false. / .true.
#       /
#       &time_advance_knobs
#         xdriftknob = 0.0 / 1.0
#         ydriftknob = 0.0 / 1.0
#         wstarknob = 0.0 / 1.0
#       /
#       &knobs
#         stream_implicit = .false. / .true.
#         mirror_implicit = .false. / .true.
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

input_filename_stem = 'FFS_adiabatic'  
stella_local_run_directory = 'Not/Run/Yet' 

#-------------------------------------------------------------------------------
#           Check whether the potential data does not evolve in time           #
#-------------------------------------------------------------------------------
def test_each_gyrokinetic_term_for_modified_adiabatic_full_flux_surface(tmp_path):
    
    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global stella_local_run_directory
    stella_local_run_directory = tmp_path
    
    input_filename = input_filename_stem + '_no_time_evolution.in'
    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path)
     
    # File names  
    local_netcdf_file = stella_local_run_directory / (input_filename_stem + '_no_time_evolution.out.nc')

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
                
    print(f'  -->  Without gyrokinetic terms the potential in FFS does not evolve in time ({int(local_netcdf["nproc"])} CPUs).')

#-------------------------------------------------------------------------------
#     Check whether the implicit parallel streaming term evolves correctly     #
#------------------------------------------------------------------------------
    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_parallel_streaming_implicit.in'
    run_local_stella_simulation(input_filename, tmp_path)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_parallel_streaming_implicit.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   
    compare_local_potential_with_expected_potential(local_netcdf_file, expected_netcdf_file, error=False)

    print('  -->  The parallel streaming term for adiabatic FFS, implemented implicitly, is evolving correctly.') 

#-------------------------------------------------------------------------------
#           Check whether the implicit mirror term evolves correctly           #
#-------------------------------------------------------------------------------
    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_mirror_implicit.in'
    run_local_stella_simulation(input_filename, tmp_path)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_mirror_implicit.out.nc')  
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    compare_local_potential_with_expected_potential(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  The mirror term for adiabatic FFS, implemented implicitly, is evolving correctly.')

#-------------------------------------------------------------------------------
#     Check whether the explicit parallel streaming term evolves correctly     #
#-------------------------------------------------------------------------------
    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_parallel_streaming_explicit.in'
    run_local_stella_simulation(input_filename, tmp_path)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_parallel_streaming_explicit.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'   
    compare_local_potential_with_expected_potential(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  The parallel streaming term for adiabatic FFS, implemented explicitly, is evolving correctly.')

#-------------------------------------------------------------------------------
#           Check whether the explicit mirror term evolves correctly           #
#-------------------------------------------------------------------------------
    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_mirror_explicit.in'
    run_local_stella_simulation(input_filename, tmp_path) 
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_mirror_explicit.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    compare_local_potential_with_expected_potential(local_netcdf_file, expected_netcdf_file, error=False)  
    print('  -->  The mirror term for adiabatic FFS, implemented explicitly, is evolving correctly.')

#-------------------------------------------------------------------------------
#          Check whether the diagmagnetic drift term evolves correctly         #
#-------------------------------------------------------------------------------
    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_diagmagnetic_drift.in'
    run_local_stella_simulation(input_filename, tmp_path)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_diagmagnetic_drift.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    compare_local_potential_with_expected_potential(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  The diagmagnetic drift term for adiabatic FFS is evolving correctly.')


#-------------------------------------------------------------------------------
#          Check whether the magnetic drifts term evolves correctly         #
#-------------------------------------------------------------------------------
    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_magnetic_drifts.in'
    run_local_stella_simulation(input_filename, tmp_path)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_magnetic_drifts.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    compare_local_potential_with_expected_potential(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  The magnetic drift terms for adiabatic FFS is evolving correctly.')


#-------------------------------------------------------------------------------
#       Check whether the all the terms combined term evolves correctly        #
#-------------------------------------------------------------------------------
    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_all.in'
    run_local_stella_simulation(input_filename, tmp_path)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_all.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    compare_local_potential_with_expected_potential(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  All terms are running as they should when simulated together for adiabatic FFS when running linearly.') 

#-------------------------------------------------------------------------------
#              Check whether the nonlinear term evolves correctly              #
#-------------------------------------------------------------------------------
    # Run stella inside of <tmp_path> based on <input_filename>
    input_filename = input_filename_stem + '_all_nonlinear.in'
    run_local_stella_simulation(input_filename, tmp_path)
     
    # Compare phi2(t) in the netCDF files
    local_netcdf_file = tmp_path / (input_filename_stem + '_all_nonlinear.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'    
    compare_local_potential_with_expected_potential(local_netcdf_file, expected_netcdf_file, error=False)    
    print('  -->  All terms are running as they should when simulated together for adiabatic FFS when running nonlinearly.') 

    return
