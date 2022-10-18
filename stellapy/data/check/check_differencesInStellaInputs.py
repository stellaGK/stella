#!/usr/bin/python3  
import numpy as np
import os, sys
import configparser 

# Stellapy package
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])   
from stellapy.utils.files.get_filesInFolder import get_filesInFolder
from stellapy.data.input.read_inputFile import read_inFile
from stellapy.utils.commandprompt.bash import Bash
        
#===============================================================================
#                      CHECK DIFFERENCES IN STELLA INPUTS                      #
#===============================================================================

def check_differencesInStellaInputs(folder, verbose=True):   
    
    # Get the ini files
    input_files = get_filesInFolder(folder, end=".in", subfolder=True)  

    # Store the differences in a dictionary
    differences = {}
    
    # Find the differences$
    for in_file1 in input_files:
        
        # Read the files
        inputParameters1 = read_inFile(in_file1) 
        
        # Find the differences
        for in_file2 in input_files: 
            
            # Read the files
            inputParameters2 = read_inFile(in_file2) 
             
            # Go through the section 
            for knob in inputParameters1.keys():
                dict1 = set(inputParameters1[knob].items())
                dict2 = set(inputParameters2[knob].items())
                
                # Save the differences
                diff = set(dict1) ^ set(dict2) 
                for item in diff: 
                    key = item[0]
                    value = item[1] 
                    if key=="tprim" and "1" in knob: key = "tiprim"
                    if key=="tprim" and "2" in knob: key = "teprim"
                    if key not in differences:
                        differences[key] = [] 
                    if key in differences: 
                        differences[key].append(value)
                
    # Differences
    keys = sorted(differences.keys())
    if "nstep" in keys: keys.remove("nstep")
    if "tend" in keys: keys.remove("tend")
            
    # Clean up the differences
    for key in keys:   
        differences[key] = sorted(list(set(differences[key])))
                
    # Check the differences
    if verbose:
        print()
        for key in keys:
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
        
#-----------------------------------
def check_differencesInStellaIniFiles(folder):   
    
    # Get the ini files
    ini_files = get_filesInFolder(folder, end=".ini") 
    
    # Only keep the input files for which a simulation was run
    ini_files = [ f for f in ini_files if os.path.isfile(f.with_suffix(".out"))]
    
    # Store the differences in a dictionary
    differences = {}
    
    # Find the differences$
    for ini_file1 in ini_files:
        for ini_file2 in ini_files:
            
            # Read the files
            config1 = configparser.ConfigParser()
            config2 = configparser.ConfigParser()
            config1.read(ini_file1)
            config2.read(ini_file2)
            
            # Go through the section 
            for section in config1.sections():
                dict1 = config1.items(section)
                dict2 = config2.items(section)
                
                # Save the differences
                diff = set(dict1) ^ set(dict2)
                for item in diff: 
                    key = item[0]
                    value = item[1] 
                    if key not in differences:
                        differences[key] = [] 
                    if key in differences: 
                        differences[key].append(value)
                
    # Differences
    keys = sorted(differences.keys())
            
    # Clean up the differences
    for key in keys:
        differences[key] = sorted(list(set(differences[key])))
                
    # Check the differences
    for key in keys:
        print("{0:>15}".format(key), " ", "{0:<15}".format(str(differences[key])))
                
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    bash = Bash(check_differencesInStellaInputs, __doc__)   
    check_differencesInStellaInputs(**bash.get_arguments())   