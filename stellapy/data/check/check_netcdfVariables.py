#!/usr/bin/python3
import os, sys
import netCDF4 as nc

# Personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0]) 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder 
from stellapy.utils.commandprompt.bash import Bash

def check_netcdfVariables(folder):
    ''' Reduce the time dimension of the netcdf file. '''
    
    # Get the netcdf files inside this folder
    netcdf_files = get_filesInFolder(folder, end="out.nc")
     
    # Reduce the time axis of each netcdf file
    if netcdf_files: 
        for netcdf_file in netcdf_files:
                 
            # New files
            old_file = str(netcdf_file)
            print("=======================================================")
            print(old_file)
            print("=======================================================") 
             
            # Constuct a new netcdf file with less data
            with nc.Dataset(old_file) as src:
                 
                # Copy all attributes (only <title>)
                print("\nATTRIBUTES")
                for name in src.ncattrs():
                    print("   ...  ", name)
                     
                # Copy all dimensions (kx,ky,tube,theta0,zed,alpha,vpa,mu,species,t,char10,char200,nlines,ri)
                print("\nDIMENSIONS")
                for name, dimension in src.dimensions.items():
                    if dimension.isunlimited():
                        print("{0:>10}".format(name), " ", "unlimited")
                    else:
                        print("{0:>10}".format(name), " ", len(dimension))
                         
                # Copy all file data for variables that are included in the toinclude list
                print("\nVARIABLES")
                for name, variable in src.variables.items():
                    print("{0:>15}".format(name), "   ", variable.datatype, variable.dimensions)
                    if name=="qflx_kxky": saved_dimension = variable.dimensions
                print()

    return    
        
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    bash = Bash(check_netcdfVariables, __doc__)    
    check_netcdfVariables(**bash.get_arguments())
