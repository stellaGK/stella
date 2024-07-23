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
#                           Routines to compare data                           #
################################################################################
    
#-------------------------------------------------------------------------------
def convert_byte_array(array):
    '''Tool to convert netcdf text arrays, to a text string that we can compare.'''
    return '\n'.join((str(line, encoding='utf-8').strip() for line in array.data))
     
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
      
#-------------------------------------------------------------------------------  
def compare_local_array_with_expected_array(local_array, expected_array): 
    print('                                       ')
    print('                LOCAL       EXPECTED   ')
    print('                -----       --------   ')
    for i in range(np.min([len(local_array.data),10])):
        print(f' phi2[{i}]  =  {local_array.data[i]:8.5e}   {expected_array.data[i]:8.5e}')
    print('                                       ')
