#!/usr/bin/python3
import copy
import pathlib
import os, sys 
import numpy as np
import time as timer
import netCDF4 as nc
from scipy.io import netcdf as scnetcdf  
from datetime import datetime, timedelta

# Personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile  
from stellapy.utils.commandprompt.print_progressbar import print_progressbar
from stellapy.utils.files.get_filesInFolder import get_filesInFolder 
from stellapy.data.output.read_outFile import read_outFile
from stellapy.utils.commandprompt.bash import Bash

def write_netcdfFileOverSaturatedPhase(folder, skip=1, verbose=False):
    ''' Reduce the time dimension of the netcdf file. '''
    
    # Get the netcdf files inside this folder
    netcdf_files = get_filesInFolder(folder, end="out.nc")
     
    # Reduce the time axis of each netcdf file
    if netcdf_files: 
        for netcdf_file in netcdf_files:

            # Check whether we have a linear or nonlinear simulation
            linear, nonlinear = read_linearNonlinearFromInputFile(str(netcdf_file).replace(".out.nc",".in")) 
             
            # Start timer
            start = timer.time()
            
            # Check the time trace
            input_file = pathlib.Path(str(netcdf_file).replace(".out.nc", ".in"))
            time = np.array(read_outFile(input_file)["time"])
            tstart = int(time[-1]/2); tend = int(time[-1]) 
                 
            # New files
            time_frame = "t"+str(tstart)+"-"+str(tend)
            time_frame = time_frame+"-skip"+str(skip) if skip!=1 else time_frame
            old_file = pathlib.Path(netcdf_file)
            if linear: new_file = pathlib.Path(str(netcdf_file).replace(".out.nc", ".out.nc.tlast"))
            if nonlinear: new_file = pathlib.Path(str(netcdf_file).replace(".out.nc", ".out.nc."+time_frame))
            finished_file = copy.deepcopy(new_file); new_file = pathlib.Path(str(new_file)+"_in_progress")
            if verbose: print("=======================================================")
            if verbose: print(old_file)
            if verbose: print(new_file)
            if verbose: print("=======================================================") 
            
            # Check whether the saturated file contains all information 
            if os.path.isfile(finished_file): 
                with scnetcdf.netcdf_file(finished_file,'r') as netcdf_data: 
                    variables_new = list(netcdf_data.variables.keys()) 
                if "input_file" in variables_new:
                    if verbose: print("The netcdf file averaged over the saturated phase already exists. ")
                    if (not verbose) and nonlinear: print("   The saturated netcdf already exists: " +  finished_file.parent.name+"/"+finished_file.name)  
                    continue
                    
            # Constuct a new netcdf file averaged over the saturated time trace
            print("\nCreate a new netcdf file averaged over t = ["+str(tstart)+", "+str(tend)+"] with skip = "+str(skip)+" for "+finished_file.parent.name+".")
            with nc.Dataset(old_file) as src, nc.Dataset(new_file, "w", format="NETCDF3_CLASSIC") as dst:
                 
                # Copy all attributes (only <title>)
                if nonlinear:
                    if verbose: print("\nSAVE ATTRIBUTES")
                    for name in src.ncattrs():
                        if verbose: print("   ...  ", name)
                        dst.setncattr(name, src.getncattr(name))
                     
                # Copy all dimensions (kx,ky,tube,theta0,zed,alpha,vpa,mu,species,t,char10,char200,nlines,ri)
                if verbose: print("\nSAVE DIMENSIONS")
                for name, dimension in src.dimensions.items():
                    if name=="t": 
                        if verbose: print("{0:>10}".format(name), " ", "Put to 1!")
                        dst.createDimension(name, 1)
                    else: 
                        if verbose: print("{0:>10}".format(name), " ", len(dimension))
                        dst.createDimension(name, len(dimension))
     
                # Read the time axis
                if nonlinear:
                    time = src.variables["t"][:]
                    istart = np.where(time==time[time>=tstart][0])[0][0]; iend = len(time)
                         
                # Copy all file data for variables that are included in the toinclude list
                keys = src.variables.keys()
                if verbose: print("\nSAVE VARIABLES")
                print_progressbar(0, len(keys), prefix = '   Progress:', suffix = 'Complete', length = 50)
                keys = src.variables.keys(); count = 1
                for name, variable in src.variables.items():
                    if name!="dt": 
                         
                        # Create the variable 
                        if verbose: print("{0:>15}".format(name), "   ", variable.datatype, variable.dimensions)
                        dst.createVariable(name, variable.datatype, variable.dimensions)
                         
                        # Fill in the variable, but first cut off time positions   
                        if 't' in variable.dimensions and linear: dst.variables[name][:] = src.variables[name][-1]
                        if 't' in variable.dimensions and nonlinear: dst.variables[name][:] = np.mean(src.variables[name][istart:iend:skip], axis=0)
                        if 't' not in variable.dimensions and nonlinear: dst.variables[name][:] = src.variables[name][:]
                        
                    # Show progress
                    print_progressbar(count, len(keys), prefix = '   Progress:', suffix = 'Complete', length = 50); count += 1
                        
                # Time the calculation of the new time indices 
                elapsed_time = str(timedelta(seconds=round(timer.time()-start,0)))
                         
            # Final message
            os.system("mv "+str(new_file)+" "+str(finished_file))
            if verbose: print("\nSaved the new netcdf file averaged over t = ["+str(tstart)+", "+str(tend)+".")
            if verbose: print("   ---> Reduced from "+str(len(time))+" to 1 time point.")
            if verbose: print("   ---> Time reduction of a factor "+str(round(len(time)/len([1]),2))+" which took "+elapsed_time+" hours.")  
            if verbose: print()
            if not verbose: print("   ---> The saturated netcdf is saved as " +  finished_file.parent.name+"/"+finished_file.name)  
    
    return    
        
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    bash = Bash(write_netcdfFileOverSaturatedPhase) 
    #bash.add_option('dt', 'float', 't', 'Only keep the time values [0, dt, 2*dt, 3*dt, ...]') 
    write_netcdfFileOverSaturatedPhase(**bash.get_arguments())
