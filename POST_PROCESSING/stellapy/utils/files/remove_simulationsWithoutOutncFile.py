
#===============================================================
# Remove cases without an *.out.nc file
#===============================================================
# TODO: Remove all mentions of <remove_simulationsWithoutOutncFile>

import os, sys  
from stellapy.utils.decorators import printv
from stellapy.utils.decorators.verbose import noverbose 
from stellapy.utils.decorators.exit_program import exit_program

@noverbose
def remove_simulationsWithoutOutncFile(input_files, forceexit=False):
    ''' Analysis can only be performed if *out.nc file or *out.h5 file is present.

    Returns
    -------       
    input_files : list of PosixPaths
        List of input_files of which the existence of *out.nc or *out.h5 files is checked.
    '''

    # Select the cases without an *out file
    printv("Checking existence of *out file corresponding inputs.")
    count=0; not_valid=[]
    for input_file in input_files: 
        if os.path.isfile(input_file.with_suffix('.out.nc')):     
            count = count + 1  
        elif os.path.isfile(input_file.with_suffix('.out.h5')):                 
            count = count + 1
        elif os.path.isfile(input_file.with_suffix('.dimensions')):                 
            count = count + 1
        elif os.path.isfile(input_file.with_suffix('.dt1.fluxes_vs_t')):                 
            count = count + 1
        elif os.path.isfile(input_file.with_suffix('.dt1.omega_vs_t')):            
            count = count + 1
        else:
            printv("    File "+input_file.name+" removed from input list.")
            not_valid.append(input_file)
    printv("    Number of found input files with corresponding output *.nc = "+str(count))  

    # Remove the selected cases
    for input_file in not_valid:
        i = input_files.index(input_file)
        input_files.pop(i) 
        
    # Stop the program if no input files were found 
    if input_files == []: 
        printv("There were no input files found, from the following selection:")
        for i in input_files: printv("      ", str(i))
        if forceexit:
            exit_reason = "There were no input files found, from the following selection:\n" 
            for i in input_files:  exit_reason += "      ", str(i)+"\n" 
            exit_program(exit_reason, remove_simulationsWithoutOutncFile, sys._getframe().f_lineno)
        
    # Return the input files that have corresponding *.out.nc files
    return input_files
