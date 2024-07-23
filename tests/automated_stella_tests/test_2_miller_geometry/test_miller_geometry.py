################################################################################
#                           Check the Miller geometry                          #
################################################################################
# Test all the geometry arrays when using a Miller equilibrium. 
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
input_filename = 'miller_geometry.in'  
stella_local_run_directory = 'Not/Run/Yet'

#-------------------------------------------------------------------------------
#                    Check whether output files are present                    #
#-------------------------------------------------------------------------------
def test_whether_miller_output_files_are_present(tmp_path, error=False):  
    
    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global stella_local_run_directory
    stella_local_run_directory = tmp_path
    
    # Run stella inside of <tmp_path> based on <input_filename>
    run_local_stella_simulation(input_filename, tmp_path)
    
    # Gather the output files generated during the local stella run inside <tmp_path>
    local_files = os.listdir(stella_local_run_directory)
    
    # Create a list of the output files we expect when stella has been run 
    expected_files = ['millerlocal.miller_geometry.input', 'millerlocal.miller_geometry.output', 'miller_geometry.geometry']
    
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
    print(f'  -->  All the expected files (millerlocal.input, millerlocal.output, .geometry) are generated.')
    return 

#-------------------------------------------------------------------------------
#                    Check whether Miller output files match                   #
#-------------------------------------------------------------------------------
def test_whether_miller_output_files_are_correct():  
    '''Check that the results are identical to a previous run.'''

    # Initialize
    geometry_files_match = True
    miller_input_files_match = True
    miller_output_files_match = True
    
    # File names 
    local_geometry_file = stella_local_run_directory / 'miller_geometry.geometry' 
    expected_geometry_file = get_stella_expected_run_directory() / 'EXPECTED_OUTPUT.geometry' 
    local_miller_input_file = stella_local_run_directory / 'millerlocal.miller_geometry.input' 
    expected_miller_input_file = get_stella_expected_run_directory() / 'EXPECTED_OUTPUT.millerlocal.input' 
    local_miller_output_file = stella_local_run_directory / 'millerlocal.miller_geometry.output' 
    expected_miller_output_file = get_stella_expected_run_directory() / 'EXPECTED_OUTPUT.millerlocal.output' 
     
    # Check whether the millerlocal input files match  
    if os.path.getsize(local_miller_input_file) != os.path.getsize(expected_miller_input_file): miller_input_files_match = False; print('ERROR: Miller input files do not match.')
    elif open(local_miller_input_file,'r').read() != open(expected_miller_input_file,'r').read(): miller_input_files_match = False; print('ERROR: Miller input files do not match.')

    # Check whether the millerlocal output files match  
    if os.path.getsize(local_miller_output_file) != os.path.getsize(expected_miller_output_file): miller_output_files_match = False; print('ERROR: Miller output files do not match.')
    elif open(local_miller_output_file,'r').read() != open(expected_miller_output_file,'r').read(): miller_output_files_match = False; print('ERROR: Miller output files do not match.')

    # Check whether the geometry files match  
    if os.path.getsize(local_geometry_file) != os.path.getsize(expected_geometry_file): geometry_files_match = False; print('ERROR: Geometry files do not match.')
    elif open(local_geometry_file,'r').read() != open(expected_geometry_file,'r').read(): geometry_files_match = False; print('ERROR: Geometry files do not match.')
    
    # If the files don't match, print the differences
    if miller_input_files_match==False:
        print_differences_in_text_files(local_miller_input_file, expected_miller_input_file, name='MILLER INPUT')
        assert False, 'Miller input files do not match.'
    if miller_output_files_match==False:
        print_differences_in_text_files(local_miller_output_file, expected_miller_output_file, name='MILLER OUTPUT')
        assert False, 'Miller output files do not match.'
    if geometry_files_match==False:
        print_differences_in_text_files(local_geometry_file, expected_geometry_file, name='GEOMETRY')
        assert False, 'Geometry files do not match.'

    # If we made it here the test was run correctly 
    print(f'  -->  Geometry output file matches.')
    return 

#-------------------------------------------------------------------------------
#              Check whether the data in the netcdf file matches               #
#-------------------------------------------------------------------------------
def test_whether_miller_geometry_data_in_netcdf_file_is_correct(error=False):
    '''Check that the results are identical to a previous run.'''
     
    # File names  
    local_netcdf_file = stella_local_run_directory / input_filename.replace('.in', '.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / 'EXPECTED_OUTPUT.out.nc'   

    # Check whether the geometry data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
    
        # Check whether all the keys are present
        if set(local_netcdf.keys()) != set(expected_netcdf.keys()):
            print('ERROR: The quantities in the netcdf file differ from the EXPECTED_OUTPUT.out.nc file.')
            assert False, 'The quantities in the netcdf file differ from the EXPECTED_OUTPUT.out.nc file.'
        
        # Relevant keys for the geometry
        geometry_keys = ["bmag", "b_dot_grad_z", "gradpar", "gbdrift", "gbdrift0", "cvdrift", "cvdrift0", "kperp2", \
            "gds2", "gds21", "gds22", "grho", "jacob", "djacdrho", "beta", "q", "shat", "d2qdr2", "drhodpsi", "d2psidr2", "jtwist"]  
        for key in geometry_keys:
        
            # Compare integers and floats
            if expected_netcdf[key].shape == ():
                if key=='nproc': continue # The number of processors is allowed to differ
                if (local_netcdf[key] != expected_netcdf[key]):
                    print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                    print('    LOCAL:    {local_netcdf[key]}')
                    print('    EXPECTED: {expected_netcdf[key]}')
                
            # Compare texts (code_info and input_file)
            if expected_netcdf[key].dtype.kind == 'S':
                expected_netcdf_str = convert_byte_array(expected_netcdf[key])
                local_netcdf_str = convert_byte_array(local_netcdf[key])  
                if key=='input_file': continue # The input file is allowed to differ
                if (local_netcdf_str != expected_netcdf_str):
                    print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                    print('    LOCAL:    {local_netcdf[key]}')
                    print('    EXPECTED: {expected_netcdf[key]}') 
                    
            # Compare data arrays 
            else: 
                if not (np.allclose(local_netcdf[key], expected_netcdf[key], equal_nan=True)):
                    print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                    print('    LOCAL:    {local_netcdf[key]}')
                    print('    EXPECTED: {expected_netcdf[key]}')
                    
        # Print "AssertionError: <error message>" if an error was encountered
        assert (not error), f'Some Miller geometry arrays in the netcdf file did not match the previous run.' 
                
    print('  -->  All Miller geometry data in the netcdf file matches the expected output.')
    return
    

