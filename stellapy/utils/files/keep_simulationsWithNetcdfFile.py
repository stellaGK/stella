 

import os, sys   
from stellapy.utils.decorators.exit_program import exit_program

#===============================================================================
#                 Keep only input files which have a netcdf file          
#=============================================================================== 

def keep_simulationsWithNetcdfFile(input_files, forceexit=False, verbose=False):
    ''' Returns a list of input_files for which a *out.nc exists. '''
    
    # Message 
    if verbose: print("Checking existence of *out file corresponding inputs.")
    
    # Initiate 
    not_valid = []
    count = 0

    # Select the cases without an *out file
    for input_file in input_files: 
        if os.path.isfile(input_file.with_suffix('.out.nc')):     
            count = count + 1   
        else:
            if verbose: print("    File "+input_file.name+" removed from input list.")
            not_valid.append(input_file)
            
    # Message
    if verbose: print("    Number of found input files with corresponding output *.nc = "+str(count))  

    # Remove the selected cases
    for input_file in not_valid:
        i = input_files.index(input_file)
        input_files.pop(i) 
        
    # Stop the program if no input files were found 
    if input_files == []: 
        if verbose: 
            print("There were no input files found, from the following selection:")
            for i in input_files: print("      ", str(i))
        if forceexit:
            exit_reason = "There were no input files found, from the following selection:\n" 
            for i in input_files:  exit_reason += "      ", str(i)+"\n" 
            exit_program(exit_reason, keep_simulationsWithNetcdfFile, sys._getframe().f_lineno)
        
    # Return the input files that have corresponding *.out.nc files
    return input_files

