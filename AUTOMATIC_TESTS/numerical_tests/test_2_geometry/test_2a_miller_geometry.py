################################################################################
#                           Check the Miller geometry                          #
################################################################################
# Test all the geometry arrays when using a Miller equilibrium. 
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
input_filename = 'miller_geometry.in'
input_file = input_filename.replace('.in','')
stella_local_run_directory = 'Not/Run/Yet'
miller_file_name = 'geometry_miller'
run_data = {}

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                    Check whether output files are present                    #
#-------------------------------------------------------------------------------
def test_whether_miller_output_files_are_present(tmp_path, stella_version, error=False):  
    
    # Save the temporary folder <tmp_path> as a global variable so the
    # other tests can access the output files from the local stella run.
    global stella_local_run_directory, miller_file_name, run_data, input_filename, input_file
    stella_local_run_directory = tmp_path
    
    # Run stella inside of <tmp_path> based on <input_filename>
    if stella_version!='master': 
       input_filename = input_filename.replace('.in', f'_v{stella_version}.in')
       input_file = input_filename.replace('.in','')
       miller_file_name = 'millerlocal'
    run_data = run_local_stella_simulation(input_filename, tmp_path, stella_version)
    
    # Gather the output files generated during the local stella run inside <tmp_path>
    local_files = os.listdir(stella_local_run_directory)
    
    # Create a list of the output files we expect when stella has been run 
    expected_files = [f'{miller_file_name}.{input_file}.input', f'{miller_file_name}.{input_file}.output', f'{input_file}.geometry']; new_names = True
    
    # Check whether all these output files are present
    for expected_file in expected_files:
        if not (expected_file in local_files):
            print(f'ERROR: The "{expected_file}" output file was not generated when running stella.'); new_names = False
            
    # Old stella 
    if new_names==False:
        expected_files = [f'millerlocal.{input_file}.input', f'millerlocal.{input_file}.output', '{input_file}.geometry'] 
        for expected_file in expected_files:
            if not (expected_file in local_files):
                print(f'ERROR: The "{expected_file}" output file was not generated when running stella.'); error = True
        miller_file_name = 'millerlocal'
            
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
    
    # File names 
    local_geometry_file = stella_local_run_directory / f'{input_file}.geometry' 
    expected_geometry_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.miller_geometry.geometry' 
    local_miller_input_file = stella_local_run_directory / f'{miller_file_name}.{input_file}.input' 
    expected_miller_input_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.miller_geometry.millerlocal.input' 
    local_miller_output_file = stella_local_run_directory / f'{miller_file_name}.{input_file}.output' 
    expected_miller_output_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.miller_geometry.millerlocal.output' 
    
    # Check whether the txt files match for old stella
    if miller_file_name=='millerlocal':
        compare_local_txt_with_expected_txt_newstella(local_geometry_file, expected_geometry_file, name='Geometry txt', error=False)
        compare_local_txt_with_expected_txt(local_miller_input_file, expected_miller_input_file, name='Miller input txt', error=False)
        compare_old_stella_with_new_stella_geometry_output(local_miller_output_file, expected_miller_output_file, name='Miller output txt', error=False)
    
    # For new stella, we print <flux_fac> instead of <exb_nonlin_p>, and we do not print <btor> 
    if miller_file_name=='miller_geometry':
        compare_local_txt_with_expected_txt_newstella(local_geometry_file, expected_geometry_file, name='Geometry txt', error=False) 
        compare_local_txt_with_expected_txt_newstella(local_miller_input_file, expected_miller_input_file, name='Miller input txt', error=False) 
        compare_local_txt_with_expected_txt_newstella(local_miller_output_file, expected_miller_output_file, name='Miller output txt', error=False) 
    
    # If we made it here the test was run correctly 
    print(f'  -->  Geometry output file matches.')
    return 
    
#-------------------------------------------------------------------------------
def compare_old_stella_with_new_stella_geometry_output(local_file, expected_file, name, error=False):

    # Read the files
    with open(local_file) as f1, open(expected_file) as f2: 
        local_file_txt = f1.readlines()
        expected_file_txt = f2.readlines()
        
    # The first two lines should match exactly
    for i in [0,1]:
        if local_file_txt[0]!=expected_file_txt[0]: error = True
        if local_file_txt[1]!=expected_file_txt[1]: error = True
    if (error==True): 
        print(f'\nERROR: First two lines in the geometry file do not match:') 
        print(local_file_txt[0])
        print(local_file_txt[1])
        print('-------------------')
        print(expected_file_txt[0])
        print(expected_file_txt[1])
    assert (not error), f'{name} files do not match.' 
    
    # The data columns should match up to a certain number of digits
    # Testing this shows that 7 digits do not pass, but 6 do.
    number_of_lines = len(expected_file_txt)
    number_of_digits = 6; ibad = 0
    for i in range(2, number_of_lines):
        numbers_local = [np.round(float(number), number_of_digits) for number in local_file_txt[i].split(' ') if number!='']
        numbers_expected = [np.round(float(number), number_of_digits)  for number in expected_file_txt[i].split(' ') if number!='']  
        if numbers_local!=numbers_expected: error = True; ibad = i
    if (error==True): 
        print(f'\nERROR: Arrays in the geometry file do not match:') 
        print(local_file_txt[ibad])
        print(expected_file_txt[ibad])
    assert (not error), f'{name} files do not match.' 
    return error

#-------------------------------------------------------------------------------  
def compare_local_txt_with_expected_txt_newstella(local_file, expected_file, name, error=False):
     
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
    if error==True: 
        print(f'\nDIFFERENCE BETWEEN {name} FILES:\n')  
        for line in difflib.unified_diff(local_file_txt, expected_file_txt, fromfile=str(local_file), tofile=str(expected_file), lineterm=''): 
            print('    ', line)  
            break
        print(local_file)
        print(expected_file)
        
    # If the files don't match, print the differences 
    if (error==True): print(f'\nERROR: {name} files do not match.') 
    assert (not error), f'{name} files do not match.' 
    return error
     
#-------------------------------------------------------------------------------
#              Check whether the data in the netcdf file matches               #
#-------------------------------------------------------------------------------
def test_whether_miller_geometry_data_in_netcdf_file_is_correct(error=False): 
    compare_geometry_in_netcdf_files(run_data, error=False)  
    print('  -->  All Miller geometry data in the netcdf file matches the expected output.')
    return
    

