#!/usr/bin/python3
import os, sys, h5py

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.utils.files.get_filesInFolder import get_filesInFolder   
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                               CHECK GEOMETRY                                 #
#===============================================================================

def check_geometry(folder, geometry_files = []):
 
    # Get the geometry files
    geometry_files1 = get_filesInFolder(folder, start="input.unique.geometry") 
    geometry_files2 = get_filesInFolder(folder, end="rho70.geometry")  
    geometry_files3 = get_filesInFolder(folder, end="rho25.geometry")  
    if geometry_files1: geometry_files += geometry_files1
    if geometry_files2: geometry_files += geometry_files2
    if geometry_files3: geometry_files += geometry_files3

    # Check whether the geometry file contains "nfp"
    for geometry_file in geometry_files:
        with h5py.File(geometry_file, 'r') as f: 
            if "nfp" not in f.keys(): 
                print()
                print("    GEOMETRY DATA WAS MISSING INSIDE:")
                print("    "+str(geometry_file))
                sys.exit()
            else:
                print("OK: "+str(geometry_file))

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    bash = Bash(check_geometry)   
    check_geometry(**bash.get_arguments())   