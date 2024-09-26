#!/usr/bin/python3  
import os, sys 
import pathlib

# Personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder    
from stellapy.utils.commandprompt.bash import Bash
 
#===============================================================================
#               REPLACE THE NETCDF FILE BY THE REDUCED NETCDF FILE             #
#===============================================================================

def replace_netcdfFile(folder): 
  
    # Make sure we have a folder 
    if folder is None:
        return
    
    # Get the netcdf and reduced netcdf files inside this folder 
    netcdf_files_full = get_filesInFolder(folder, end="out.nc") 
    netcdf_files_redu = get_filesInFolder(folder, end="out.nc.reduced")
    
    # Only keep the files that have a reduced file
    if netcdf_files_full: 
        netcdf_files_full = [f for f in netcdf_files_full if f.with_suffix(".nc.reduced") in netcdf_files_redu]
    
    # Only continue if we found files
    if netcdf_files_full is None:
        print("   EXIT: No netcdf file found."); sys.exit(0) 
    if len(netcdf_files_full)==0:
        print("   EXIT: No netcdf file found."); sys.exit(0) 
    if netcdf_files_redu is None:
        print("   EXIT: No reduced netcdf file found."); sys.exit(0) 
    if len(netcdf_files_redu)==0:
        print("   EXIT: No reduced netcdf file found."); sys.exit(0) 
        
    # Replace the *.out.nc file by the *.out.nc.reduced file
    for netcdf_file in netcdf_files_full:
        reduced_file = netcdf_file.with_suffix(".nc.reduced")
        print("   --> remove "+str(netcdf_file.name), "and replace by "+str(reduced_file.name))
        os.system("rm "+str(netcdf_file))  
        os.system("mv "+str(reduced_file)+" "+str(netcdf_file)) 
    return    

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    bash = Bash(replace_netcdfFile, __doc__)  
    replace_netcdfFile(**bash.get_arguments())
