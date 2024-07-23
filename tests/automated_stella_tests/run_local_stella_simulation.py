#!/usr/bin/python3   
import os
import shutil
import pathlib
import subprocess   

################################################################################
#                 Routines to launch a local stella simulation                 #
################################################################################
 
def get_stella_expected_run_directory():
    '''Get the directory of this test file.'''
    return pathlib.Path(__file__).parent

def get_stella_path():
    '''Returns the absolute path to the stella executable.
    Can be controlled by setting the STELLA_EXE_PATH environment variable.'''
    default_path = get_stella_expected_run_directory() / '../../../stella'
    stella_path = pathlib.Path(os.environ.get('STELLA_EXE_PATH', default_path))
    return stella_path.absolute()

def run_stella(stella_path, input_filename):
    '''Run stella with a given input file.''' 
    subprocess.run(['mpirun', '-np', '2', stella_path, input_filename], check=True)

def copy_input_file(input_filename: str, destination):
    '''Copy input_filename to destination directory.'''
    shutil.copyfile(get_stella_expected_run_directory() / input_filename, destination / input_filename)

def copy_vmec_file(vmec_filename: str, destination):
    '''Copy input_filename to destination directory.'''
    shutil.copyfile(get_stella_expected_run_directory() / vmec_filename, destination / vmec_filename)
    
def convert_byte_array(array):
    '''Tool to convert netcdf text arrays, to a text string that we can compare.'''
    return '\n'.join((str(line, encoding='utf-8').strip() for line in array.data))
    
def run_local_stella_simulation(input_filename, tmp_path, vmec_filename=None):
    ''' Run a local stella simulation in <tmp_path>. '''
    copy_input_file(input_filename, tmp_path)
    if vmec_filename: copy_vmec_file(vmec_filename, tmp_path)
    os.chdir(tmp_path); run_stella(get_stella_path(), input_filename)
    return 
