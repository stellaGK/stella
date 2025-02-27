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

# Global variables  
input_filename = 'miller_linear_CBC.in'  
local_stella_run_directory = 'Not/Run/Yet'

#-------------------------------------------------------------------------------
#                         Run local stella simulation                          #
#-------------------------------------------------------------------------------
def test_whether_we_can_run_a_local_stella_simulation(tmp_path):
    '''Run a local stella simulation in a temporary folder <tmp_path>.'''  
    
    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global local_stella_run_directory
    local_stella_run_directory = tmp_path
    
    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path)
    print('\n  -->  Successfully ran a local stella simulation.')
    return 
    
#-------------------------------------------------------------------------------
#                    Check whether output files are present                    #
#-------------------------------------------------------------------------------
def test_whether_all_output_files_are_gerenated_when_running_stella(error=False):  
    '''To ensure that stella ran correctly, check that all output files are generated.'''
    
    # Gather the output files generated during the local stella run inside <tmp_path>
    local_files = os.listdir(local_stella_run_directory)
    
    # Create a list of the output files we expect when stella has been run
    stem = input_filename.replace('.in','')
    expected_files = ['in', 'out.nc', 'geometry', 'fluxes', 'omega', 'out', 'final_fields', 'species.input', 'error']
    expected_files = [stem+'.'+extension for extension in expected_files]
    
    # Check whether all these output files are present
    for expected_file in expected_files:
        if not (expected_file in local_files):
            print(f'ERROR: The "{expected_file}" output file was not generated when running stella.'); error = True
            
    # The <pytest> module will check whether all <assert> statements are true,
    # if it runs into a statement which is false, the test will be labeled as
    # "Failed" and the string in the second argument of the <assert> statement 
    # will be printed to the command prompt as "AssertionError: <error message>"
    assert (not error), f'Some output files were not generated when running stella.'
    
    # The test will stop as soon as an <assert> error was triggered, hence we
    # will only reach this final line of code if the test ran successfully
    print(f'  -->  All the expected files (.in, .geometry, .out.nc, .fluxes, .omega, ...) are generated.')
    return 
        
#-------------------------------------------------------------------------------
#         Check whether all quantities are present in the netcdf file          #
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#         Check whether all quantities are present in the netcdf file          #
#-------------------------------------------------------------------------------
def test_whether_correct_quantities_are_present_in_netcdf_file(error=False):
    '''Check whether the correct quantities are present in the netcdf output file.''' 
    
    # Find the netcdf output file which was generated during our local stella run
    local_netcdf_file = local_stella_run_directory / input_filename.replace('.in', '.out.nc')
    
    # Open the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf:
    
        # Check whether all the dimensions are present in the netcdf file
        expected_keys = ['kx', 'ky', 'tube', 'zed', 'alpha', 'vpa', 'mu', 'species', 't', 'char10', 'char200', 'ri']  
        for key in expected_keys:
            if (key not in local_netcdf): print(f'ERROR: The quantity {key} could not be found in the netcdf file.'); error = True
        if not error: print(f'  -->  The netcdf file contains the following dimensions: {expected_keys}.')  
    
        # Check whether all the species information is present in the netcdf file
        expected_keys = ['charge', 'mass', 'dens', 'temp', 'tprim', 'fprim', 'vnew', 'type_of_species'] 
        for key in expected_keys:
            if (key not in local_netcdf): print(f'ERROR: The quantity {key} could not be found in the netcdf file.'); error = True
        if not error: print(f'  -->  The netcdf file contains the following species information: {expected_keys}.')
         
        # Check whether all the geometry data is present in the netcdf file
        expected_keys = ['bmag', 'b_dot_grad_z', 'gradpar', 'gbdrift', 'gbdrift0', 'cvdrift', 'cvdrift0', 'kperp2', \
            'gds2', 'gds21', 'gds22', 'grho', 'jacob', 'djacdrho', 'beta', 'q', 'shat', 'd2qdr2', 'drhodpsi', 'd2psidr2', 'jtwist'] 
        for key in expected_keys:
            if (key not in local_netcdf): print(f'ERROR: The quantity {key} could not be found in the netcdf file.'); error = True
        if not error: print(f'  -->  The netcdf file contains the following geometry data: {expected_keys}.')
        
        # Check whether all the diagnostics data is present in the netcdf file
        dist_extensions = ['_vs_vpamus', '_vs_zvpas', '_vs_zmus', '_vs_zvpamus']
        expected_keys = ['t', 'omega', 'phi2', 'apar2', 'bpar2', 'omega']; new_keys = True
        expected_keys += [i+j for i in ['pflux', 'qflux', 'vflux'] for j in ['_vs_s', '_vs_kxkyzs']]
        expected_keys += [i+j for i in ['g2', 'h2', 'f2', 'g2nozonal', 'h2nozonal', 'f2nozonal'] for j in dist_extensions]
        expected_keys += ['g2_vs_zkykxs', 'h2_vs_zkykxs', 'f2_vs_zkykxs']
        for key in expected_keys:
            if (key not in local_netcdf): 
                print(f'ERROR: The quantity {key} could not be found in the netcdf file.'); new_keys = False
            
        # Try old stella names
        if new_keys==False:
            print('--------------------------------------------------')
            print(' The new stella keys werent found, try old keys.')
            print('--------------------------------------------------')
            expected_keys = ['t', 'omega', 'phi2', 'apar2', 'bpar2', 'pflx', 'vflx', 'phi_vs_t', 'phi2_vs_kxky', \
               'pflx_vs_kxky', 'vflx_vs_kxky', 'qflx_vs_kxky', 'pflux_x', 'vflux_x', 'qflux_x', \
               'dens_x', 'upar_x', 'temp_x', 'pflx_kxky', 'vflx_kxky', 'qflx_kxky', 'gvmus', 'gzvs']  
            for key in expected_keys:
               if (key not in local_netcdf): 
                   print(f'ERROR: The quantity {key} could not be found in the netcdf file.'); error = True
                       
        if not error: print(f'  -->  The netcdf file contains the following diagnostics data: {expected_keys}.')
            
        # The <pytest> module will check whether all <assert> statements are true,
        # if it runs into a statement which is false, the test will be labeled as
        # "Failed" and the string in the second argument of the <assert> statement 
        # will be printed to the command prompt as "AssertionError: <error message>"
        assert (not error), f'Some quantities were not found in the netcdf file.' 
        
        # The test will stop as soon as an <assert> error was triggered, hence we
        # will only reach this final line of code if the test ran successfully
        print(f'  -->  The netcdf file contains all the expected data.')
        
    return 
     
