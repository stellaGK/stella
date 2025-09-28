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
    
    # Compare text files
    compare_miller_geometry_files(local_geometry_file, expected_geometry_file, error=False)
    compare_miller_input_files(local_miller_input_file, expected_miller_input_file, name='Miller input txt', error=False)
    
    # If we made it here the test was run correctly 
    print(f'  -->  Geometry output file matches.')
    return 
    
#-------------------------------------------------------------------------------
def process_error(variable):
    print(f'\nERROR: {variable} does not match in the *.geometry file.'); 
    return True
    
#-------------------------------------------------------------------------------
def compare_miller_geometry_files(local_geometry_file, expected_geometry_file, error=False):

    # Read variables in old *.geometry file
    with open(expected_geometry_file, "r") as f:
       lines = f.readlines()
       variables = lines[1]
    variables = variables.replace('#','').replace('\n','').split(' ')
    variables = [i for i in variables if i!='']
    rhoc = float(variables[0])
    qinp = float(variables[1])
    shat = float(variables[2])
    rhotor = float(variables[3])
    aref = float(variables[4])
    bref = float(variables[5])
    dxdpsi = float(variables[6])
    dydalpha = float(variables[7])
    exb_nonlin = float(variables[8])
    fluxfac = float(variables[9])
    
    # Read arrays in old *.geometry file
    data = np.loadtxt(expected_geometry_file,skiprows=4,dtype='float').reshape(-1, 15)
    alpha_old = data[:,0]
    zed_old = data[:,1]
    zeta_old = data[:,2]
    bmag_old = data[:,3]
    b_dot_gradz_old = data[:,4]
    gds2_old = data[:,5]
    gds21_old = data[:,6]
    gds22_old = data[:,7]
    gds23_old = data[:,8]
    gds24_old = data[:,9]
    gbdrift_old = data[:,10]
    cvdrift_old = data[:,11]
    gbdrift0_old = data[:,12]
    bmag_psi0_old = data[:,13]
    btor_old = data[:,14]
    
    # Read variables in new *.geometry file
    with open(local_geometry_file, "r") as f:
       lines = f.readlines()
       variables = lines[1]
    variables = variables.replace('#','').replace('\n','').split(' ')
    variables = [i for i in variables if i!='']
    rhoc = float(variables[0])
    qinp = float(variables[1])
    shat = float(variables[2])
    rhotor = float(variables[3])
    aref = float(variables[4])
    bref = float(variables[5])
    dxdpsi = float(variables[6])
    dydalpha = float(variables[7])
    exb_nonlin = float(variables[8])
    fluxfac = float(variables[9])
    
    # Read arrays in new *.geometry file
    data = np.loadtxt(local_geometry_file,skiprows=4,dtype='float').reshape(-1, 15)
    alpha_new = data[:,0]
    zed_new = data[:,1]
    zeta_new = data[:,2]
    bmag_new = data[:,3]
    b_dot_gradz_new = data[:,4]
    gds2_new = data[:,5]
    gds21_new = data[:,6]
    gds22_new = data[:,7]
    gds23_new = data[:,8]
    gds24_new = data[:,9]
    gbdrift_new = data[:,10]
    cvdrift_new = data[:,11]
    B_times_gradB_dot_gradx_new = data[:,12]
    bmag_psi0_new = data[:,13]
    btor_new = data[:,14]
    
    # New definitions
    digits = 2
    gbdrift0_new = B_times_gradB_dot_gradx_new * 2 * shat
    gbdrift0_new = np.round(gbdrift0_new, digits)
    gbdrift0_old = np.round(gbdrift0_old, digits)
    
    # Compare arrays
    error = False
    if not (np.allclose(alpha_old, alpha_new, equal_nan=True)): error = process_error('alpha')
    if not (np.allclose(zed_old, zed_new, equal_nan=True)): error = process_error('zed')
    if not (np.allclose(zeta_old, zeta_new, equal_nan=True)): error = process_error('zeta')
    if not (np.allclose(bmag_old, bmag_new, equal_nan=True)): error = process_error('bmag')
    if not (np.allclose(b_dot_gradz_new, b_dot_gradz_new, equal_nan=True)): error = process_error('b_dot_gradz')
    if not (np.allclose(gds2_old, gds2_new, equal_nan=True)): error = process_error('gds2')
    if not (np.allclose(gds21_old, gds21_new, equal_nan=True)): error = process_error('gds21')
    if not (np.allclose(gds23_old, gds23_new, equal_nan=True)): error = process_error('gds23')
    if not (np.allclose(gds24_old, gds24_new, equal_nan=True)): error = process_error('gds24')
    if not (np.allclose(gbdrift_old, gbdrift_new, equal_nan=True)): error = process_error('gbdrift')
    if not (np.allclose(cvdrift_old, cvdrift_new, equal_nan=True)): error = process_error('cvdrift')
    if not (np.allclose(gbdrift0_old, gbdrift0_new, equal_nan=True)): error = process_error('gbdrift0')
    if not (np.allclose(bmag_psi0_old, bmag_psi0_new, equal_nan=True)): error = process_error('bmag_psi0')
    if not (np.allclose(btor_old, btor_new, equal_nan=True)): error = process_error('btor')
    assert (not error), f'The geometry data does not match in the *.geometry file.'
    return

#-------------------------------------------------------------------------------  
def compare_miller_input_files(local_file, expected_file, name, error=False):
     
    # If both tests were run with new stella it will be fine
    # If the files match return without giving an error
    if os.path.getsize(local_file) == os.path.getsize(expected_file) and open(local_file,'r').read() == open(expected_file,'r').read(): 
        return

    # Check whether the text files match
    with open(local_file) as f1, open(expected_file) as f2: 
        local_file_txt = f1.readlines()
        expected_file_txt = f2.readlines()
    if local_file_txt!=expected_file_txt:
        error = True
        
    for i in range(len(local_file_txt)):
        print(i, local_file_txt[i] == expected_file_txt[i])
        
    # If the files don't match, print the differences 
    if (error==True):
      print(f'\nERROR: {name} files do not match.')
      print('---------------------')
      print(local_file_txt)
      print('---------------------')
      print(expected_file_txt)
    assert (not error), f'{name} files do not match.'
    return error
     
#-------------------------------------------------------------------------------
#              Check whether the data in the netcdf file matches               #
#-------------------------------------------------------------------------------
def off_test_whether_miller_geometry_data_in_netcdf_file_is_correct(error=False): 
    compare_geometry_in_netcdf_files(run_data, error=False)  
    print('  -->  All Miller geometry data in the netcdf file matches the expected output.')
    return
    

