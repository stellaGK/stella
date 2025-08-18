################################################################################
#                            Check the VMEC geometry                           #
################################################################################
# Test all the geometry arrays when using a VMEC equilibrium. 
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

# Global variables 
vmec_filename = 'wout_w7x_kjm.nc'
input_filename = 'vmec_geometry.in'  
stella_local_run_directory = 'Not/Run/Yet' 

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")
    
#-------------------------------------------------------------------------------
#                    Check whether output files are present                    #
#-------------------------------------------------------------------------------
def test_whether_VMEC_output_files_are_present(tmp_path, stella_version, error=False):  
    
    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global stella_local_run_directory, input_filename
    stella_local_run_directory = tmp_path
    
    # Run stella inside of <tmp_path> based on <input_filename>
    if stella_version!='master': 
        input_filename = input_filename.replace('.in', f'_v{stella_version}.in') 
        print('WARNING: TODO: Not implemented yet for stella versions 0.5, 0.6, 0.7')
        return 
    run_local_stella_simulation(input_filename, tmp_path, stella_version, vmec_filename)
    
    # Gather the output files generated during the local stella run inside <tmp_path>
    local_files = os.listdir(stella_local_run_directory)
    
    # Create a list of the output files we expect when stella has been run 
    expected_files = ['vmec_geometry.vmec.geo', 'vmec_geometry.geometry']
    
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
    print(f'  -->  All the expected files (.vmec.geo, .geometry) are generated.')
    return 

#-------------------------------------------------------------------------------
#                     Check whether VMEC output files match                    #
#-------------------------------------------------------------------------------
def offtest_whether_vmec_output_files_are_correct():  
    # TODO-HT Fix this test

    # Initialize
    geometry_files_match = True
    vmec_files_match = True 
    
    # File names 
    local_geometry_file = stella_local_run_directory / 'vmec_geometry.geometry' 
    expected_geometry_file = get_stella_expected_run_directory() / 'EXPECTED_OUTPUT.vmec_geometry.geometry' 
    local_vmec_file = stella_local_run_directory / 'vmec_geometry.vmec.geo' 
    expected_vmec_file = get_stella_expected_run_directory() / 'EXPECTED_OUTPUT.vmec_geometry.vmec.geo'  
    
    # Check whether the txt files match  
    # compare_local_txt_with_expected_txt(local_geometry_file, expected_geometry_file, name='Geometry txt', error=False)
    # compare_local_txt_with_expected_txt(local_vmec_file, expected_vmec_file, name='vmec.geo txt', error=False)  

    # For new stella, we print <flux_fac> instead of <exb_nonlin_p>, and we do not print <btor> 
    compare_local_txt_with_expected_txt_newstella(local_geometry_file, expected_geometry_file, name='Geometry txt', error=False) 
    # TODO-HT Add vmec.geo

    # If we made it here the test was run correctly 
    print(f'  -->  The geometry output file matches.')
    return 
        
#-------------------------------------------------------------------------------  
def compare_local_txt_with_expected_txt_newstella(local_file, expected_file, name, error=False, maxlines=10):
     
    # If both tests were run with new stella it will be fine
    # If the files match return without giving an error
    if os.path.getsize(local_file) == os.path.getsize(expected_file) and open(local_file,'r').read() == open(expected_file,'r').read(): 
        return 
    
    # Cut-off <flux_fac> or <exb_nonlin_p> in the last column of the variables
    # and cut-off <btor> in the last column of the arrays
    with open(local_file) as f1, open(expected_file) as f2: 
        local_file_txt = f1.readlines() 
        expected_file_txt = f2.readlines()  
        
    # Modify the geometry files slightly so they match whether they come from old or new stella
    if '.geometry' in str(local_file):
        for file_txt in [local_file_txt, expected_file_txt]: 
            for i, line in enumerate(file_txt):
                if len(line)==181: file_txt[i] = file_txt[i][:168] # Cut off <btor> in the old stella files 
                if len(line)==169: file_txt[i] = file_txt[i][:168] # Cut off <btor> in the old stella files 
                if len(line)==142: file_txt[i] = file_txt[i][:129] # Cut off <flux_fac> or <exb_nonlin_p> in the last column of the variables
                if 'dxdXcoord' in line: file_txt[i] = file_txt[i].replace('dxdXcoord','   dxdpsi')

    # Check whether the txt files match now 
    if local_file_txt!=expected_file_txt:
        error = True
    
    # If they do not match show the differences
    if error == True:
        if len(local_file_txt)>maxlines:
            print(f'    Warning: the compared files are very long, we limit them to {maxlines} lines')
            local_file_txt = local_file_txt[:maxlines]; expected_file_txt = expected_file_txt[:maxlines]
        print(f'\nDIFFERENCE BETWEEN {name} FILES:\n')  
        for line in difflib.unified_diff(local_file_txt, expected_file_txt, fromfile=str(local_file), tofile=str(expected_file), lineterm=''): 
            print('    ', line)  
        
    # If the files don't match, print the differences 
    if (error==True): print(f'\nERROR: {name} files do not match.') 
    assert (not error), f'{name} files do not match.' 
    return error
    
#-------------------------------------------------------------------------------
#              Check whether the data in the netcdf file matches               #
#-------------------------------------------------------------------------------
def test_whether_VMEC_geometry_data_in_netcdf_file_is_correct(error=False):
    '''Check that the results are identical to a previous run.'''
    
    # Turn off tests for older stella versions for now
    if stella_version!='master':  
        print('WARNING: TODO: Not implemented yet for stella versions 0.5, 0.6, 0.7')
        return 
     
    # File names  
    local_netcdf_file = stella_local_run_directory / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / 'EXPECTED_OUTPUT.vmec_geometry.out.nc'   

    # Check whether the geometry data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        
        # Relevant keys for the geometry
        geometry_keys = ["bmag", "b_dot_grad_z", "gradpar", "gbdrift", "gbdrift0", "cvdrift", "cvdrift0", "kperp2", \
            "gds2", "gds21", "gds22", "grho", "jacob", "djacdrho", "q", "shat", "d2qdr2", "drhodpsi", "d2psidr2", "jtwist"]  
        for key in geometry_keys:
        
            # Compare integers and floats
            if expected_netcdf[key].shape == ():
                if key=='nproc': continue # The number of processors is allowed to differ
                if (local_netcdf[key] != expected_netcdf[key]):
                    print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                    print(f'    LOCAL:    {local_netcdf[key].data}')
                    print(f'    EXPECTED: {expected_netcdf[key].data}')
                
            # Compare texts (code_info and input_file)
            if expected_netcdf[key].dtype.kind == 'S':
                expected_netcdf_str = convert_byte_array(expected_netcdf[key])
                local_netcdf_str = convert_byte_array(local_netcdf[key])  
                if key=='input_file': continue # The input file is allowed to differ
                if (local_netcdf_str != expected_netcdf_str):
                    print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                    print(f'    LOCAL:    {local_netcdf_str}')
                    print(f'    EXPECTED: {expected_netcdf_str}') 
                    
            # Compare data arrays 
            else: 
                if not (np.allclose(local_netcdf[key], expected_netcdf[key], equal_nan=True)):
                    print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                    print('\nCompare the {key} arrays in the local and expected netCDF files:')
                    compare_local_array_with_expected_array(local_netcdf[key], expected_netcdf[key])  
                    
        # Print "AssertionError: <error message>" if an error was encountered
        assert (not error), f'Some Miller geometry arrays in the netCDF file did not match the previous run.' 
                
    print('  -->  All VMEC geometry data in the netCDF file matches the expected output.')
    return
