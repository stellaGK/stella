#!/usr/bin/python3  
import os, sys 
import pathlib
import numpy as np

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.data.input import get_differencesInStellaInputs 
from stellapy.utils.commandprompt.bash import Bash
        
#===============================================================================
#                      CHECK DIFFERENCES IN STELLA INPUTS                      #
#===============================================================================

def check_differencesInStellaInputs(folder, input_files=None, verbose=True):   
    
    # Get the differences in the stella input files
    differences = get_differencesInStellaInputs(folder, input_files)
             
    # Check the differences
    if verbose:
        print()
        for key in differences.keys():
            if key in ["delt"]:
                if len(differences[key])==2:
                    if (not np.isnan(differences[key][0])) and (not np.isnan(differences[key][1])):
                        print("{0:>15}".format(key), " ", "{0:<15}".format(str(differences[key])))  
                if len(differences[key])!=2:
                    print("{0:>15}".format(key), " ", "{0:<15}".format(str(differences[key])))    
            elif key not in ["ginit_option", "delt_option"]:
                print("{0:>15}".format(key), " ", "{0:<15}".format(str(differences[key])))
        print()
    return differences

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    bash = Bash(check_differencesInStellaInputs, __doc__)   
    check_differencesInStellaInputs(**bash.get_arguments())   