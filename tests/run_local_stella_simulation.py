#!/usr/bin/python3   
import os
import shutil
import pathlib
import subprocess   

################################################################################
#                 Routines to launch a local stella simulation                 #
################################################################################
# Note that the argument of any test function is the temporary path where the 
# test is performed, i.e., test_*(tmp_path) executes stella in <tmp_path>. 
################################################################################
 
#-------------------------------------------------------------------------------
def get_stella_expected_run_directory():
    '''Get the directory of this test file.'''
    return pathlib.Path(__file__).parent

#-------------------------------------------------------------------------------
def get_stella_path():
    '''Returns the absolute path to the stella executable.
    Can be controlled by setting the STELLA_EXE_PATH environment variable.'''
    default_path = get_stella_expected_run_directory() / '../../../stella'
    stella_path = pathlib.Path(os.environ.get('STELLA_EXE_PATH', default_path))
    return stella_path.absolute()

#-------------------------------------------------------------------------------
def run_stella(stella_path, input_filename):
    '''Run stella with a given input file.''' 
    subprocess.run(['mpirun', '-np', '2', stella_path, input_filename], check=True)

#-------------------------------------------------------------------------------
def copy_input_file(input_filename: str, destination):
    '''Copy input_filename to destination directory.'''
    shutil.copyfile(get_stella_expected_run_directory() / input_filename, destination / input_filename)

#-------------------------------------------------------------------------------
def copy_vmec_file(vmec_filename: str, destination):
    '''Copy input_filename to destination directory.'''
    shutil.copyfile(get_stella_expected_run_directory() / vmec_filename, destination / vmec_filename)
    
#-------------------------------------------------------------------------------
def run_local_stella_simulation(input_filename, tmp_path, vmec_filename=None):
    ''' Run a local stella simulation in <tmp_path>. '''
    copy_input_file(input_filename, tmp_path)
    if vmec_filename: copy_vmec_file(vmec_filename, tmp_path)
    os.chdir(tmp_path); run_stella(get_stella_path(), input_filename)
    return 
    
################################################################################
#                         Routines to compare txt files                        #
################################################################################
    
#-------------------------------------------------------------------------------  
def compare_local_txt_with_expected_txt(local_file, expected_file, name, error=False):
     
    # Check whether the files match  
    if os.path.getsize(local_file) != os.path.getsize(expected_file): error = True
    elif open(local_file,'r').read() != open(expected_file,'r').read(): error = True 
    
    # If the files don't match, print the differences
    if error==True:
        print(f'ERROR: {name} files do not match.')
        print_differences_in_text_files(local_file, expected_file, name=name.upper())
        assert False, f'{name} files do not match.' 
     
#-------------------------------------------------------------------------------
def print_differences_in_text_files(file1, file2, name='', maxlines=10): 
    print(f'\nDIFFERENCE BETWEEN {name} FILES:\n') 
    with open(file1) as f1:
        local_file = f1.readlines()
    with open(file2) as f2:
        expected_file = f2.readlines()  
    if len(local_file)>maxlines:
        print(f'    Warning: the compared files are very long, we limit them to {maxlines} lines')
        local_file = local_file[:maxlines]; expected_file = expected_file[:maxlines]
    for line in difflib.unified_diff(local_file, expected_file, fromfile=str(file1), tofile=str(file2), lineterm=''): 
        print('    ', line) 
        
################################################################################
#                       Routines to compare netCDF files                       #
################################################################################
    
#-------------------------------------------------------------------------------
def compare_local_netcdf_quantity_to_expected_netcdf_quantity(local_netcdf_file, expected_netcdf_file, key, error=False):

    # Check whether the potential data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
    
        # Check whether the key is present
        if key not in local_netcdf.keys() and key not in expected_netcdf.keys():
            print(f'\nERROR: The key "{key}" does not exist in the netCDF files.')
            assert False, f'The key "{key}" does not exist in the netCDF files.'
        elif key not in local_netcdf.keys() and key in expected_netcdf.keys():
            print(f'\nERROR: The key "{key}" is present in the expected netCDF file, but not')
            print(f'         in the local netcdf file. Check whether the diagnostics has been')
            print(f'         printed and whether the key of the diagnostics has changed.')
            assert False, f'The key "{key}" does not exist in the local netCDF file.'
        elif key in local_netcdf.keys() and key not in expected_netcdf.keys():
            print(f'\nERROR: The key "{key}" is present in the local netCDF file, but not')
            print(f'         in the expected netcdf file. This mightbe a new diagnostic,')
            print(f'         hence the expected output file should be updated.')
            assert False, f'The key "{key}" does not exist in the expected netCDF file.'
        
        # Read the quantity
        local_quantity = local_netcdf[key]
        expected_quantity = expected_netcdf[key] 
                     
        # Check whether the quantity matches
        if not (np.allclose(local_quantity, expected_quantity, equal_nan=True)):
            print(f'\nERROR: The {key} arrays do not match in the netCDF files.'); error = True
            print('\nCompare the {key} arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_quantity, expected_quantity)   
        assert (not error), f'The {key} data does not match in the netCDF files.' 
        
    return error
        
#-------------------------------------------------------------------------------  
def compare_local_potential_with_expected_potential(local_netcdf_file, expected_netcdf_file, error=False): 
    
    # Check whether the potential data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        
        # Read the time axis
        local_time = local_netcdf['t']
        expected_time = expected_netcdf['t']
        
        # Read the potential axis
        local_phi2 = local_netcdf['phi2']
        expected_phi2 = expected_netcdf['phi2'] 
                     
        # Check whether we have the same time and potential data
        if not (np.allclose(local_time, expected_time, equal_nan=True)):
            print('\nERROR: The time axis does not match in the netCDF files.'); error = True
            print('\nCompare the time arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_time, expected_time)  
        if not (np.allclose(local_phi2, expected_phi2, equal_nan=True)):
            print('\nERROR: The potential data does not match in the netCDF files.'); error = True 
            print('\nCompare the potential arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_phi2, expected_phi2) 
        assert (not error), f'The potential data does not match in the netCDF files.' 
    
    return error
    
#-------------------------------------------------------------------------------
def convert_byte_array(array):
    '''Tool to convert netcdf text arrays, to a text string that we can compare.'''
    return '\n'.join((str(line, encoding='utf-8').strip() for line in array.data))
      
#-------------------------------------------------------------------------------  
def compare_local_array_with_expected_array(local_array, expected_array): 
    print('                                       ')
    print('                LOCAL       EXPECTED   ')
    print('                -----       --------   ')
    for i in range(np.min([len(local_array.data),10])):
        print(f' phi2[{i}]  =  {local_array.data[i]:8.5e}   {expected_array.data[i]:8.5e}')
    print('                                       ')
