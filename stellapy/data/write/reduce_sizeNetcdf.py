#!/usr/bin/python3
import pickle
import os, sys
import pathlib
import datetime
import time as timer
import netCDF4 as nc

# Personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder 
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                       REDUCE THE SIZE OF THE NETCDF FILE                     #
#===============================================================================

def reduce_sizeNetcdfBigFiles(folder, dt=1, verbose=True):
    ''' Reduce the time dimension of the netcdf file. '''
    
    # Get the netcdf files inside this folder
    netcdf_files = get_filesInFolder(folder, end="out.nc")
     
    # Reduce the time axis of each netcdf file
    if netcdf_files: 
        for netcdf_file in netcdf_files:
            
            # Start timer
            start = timer.time()
            
            # Check which files exist
            file_step1 = os.path.isfile(str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step1"))
            file_step2 = os.path.isfile(str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step2"))
            file_step3 = os.path.isfile(str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step3"))
            file_step4 = os.path.isfile(str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step4"))
            file_step5 = os.path.isfile(str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step5"))
            
            # Remove the in progress files
            progress_files = get_filesInFolder(folder, end="_in_progress")
            if progress_files: os.system("rm "+str(netcdf_file.parent)+"/*_in_progress")
                            
            #=================================================================
            #                       Save the dimensions
            #=================================================================
            
            if (not file_step1) and (not file_step2) and (not file_step3) and (not file_step4) and (not file_step5):
            
                # Start timer
                start_step1 = timer.time()
                     
                # New files
                old_file = str(netcdf_file)
                new_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step1_in_progress")  
                finished_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step1")  
                
                # Print a message
                print("=======================================================")
                print(old_file if ("SCRATCH" not in old_file) else old_file.split("SCRATCH")[-1])
                print(new_file if ("SCRATCH" not in new_file) else new_file.split("SCRATCH")[-1])
                print("=======================================================") 
    
                # Constuct a new netcdf file with less data
                with nc.Dataset(old_file, 'r') as src, nc.Dataset(new_file, "w", format="NETCDF4") as dst:
                     
                    # Save the chosen time step
                    dst.createVariable('dt', 'f8', ())
                    dst.variables['dt'][:] = dt
                     
                    # Copy all attributes (only <title>)
                    if verbose: print("\nSAVE ATTRIBUTES")
                    for name in src.ncattrs():
                        if verbose: print("   ...  ", name)
                        dst.setncattr(name, src.getncattr(name))
                         
                    # Copy all dimensions (kx,ky,tube,theta0,zed,alpha,vpa,mu,species,t,char10,char200,nlines,ri)
                    if verbose: print("\nSAVE DIMENSIONS")
                    for name, dimension in src.dimensions.items():
                        if dimension.isunlimited():
                            if verbose: print("{0:>10}".format(name), " ", "unlimited")
                            dst.createDimension(name, None)
                        else:
                            if verbose: print("{0:>10}".format(name), " ", len(dimension))
                            dst.createDimension(name, len(dimension))
                             
                    # To reduce the time axis, calculate the length of the new time axis
                    # And the indixes of the data points that need to be kept 
                    if verbose: print("\nREAD TIME AXIS")
                    time = src.variables["t"][:]
                    time_reduced = []
                    current_time = 0 
                    print("   ---> Calculate time indices for dt = "+str(dt)+".")
                    for i in range(len(time)):
                        if time[i]>=current_time:
                            time_reduced.append(time[i])
                            current_time += dt 
                    indexes = [list(time).index(t) for t in time_reduced]
                    
                    # Save the indexes and time reduction as a pickle
                    with open(netcdf_file.parent / 'indexes.pkl','wb') as f:
                        pickle.dump(indexes, f) 
                    with open(netcdf_file.parent / 'time_reduction.pkl','wb') as f:
                        pickle.dump([len(time), len(time_reduced)], f) 
                    
                    # If the file is already reduced, don't write the file
                    if len(time)==len(time_reduced):
                        print("   ---> The file already has data at every dt = "+str(dt)+".\n")
                        os.system("rm "+str(new_file))
                        continue
                    
                    # Otherwise start the reduction 
                    print("   ---> Reduce from "+str(len(time))+" to "+str(len(time_reduced))+" time points.")  
                             
                    # Copy all file data for variables that are included in the toinclude list
                    if verbose: print("\nSAVE VARIABLES")
                    for name, variable in src.variables.items():
                        if name!="dt": 
                            
                            # Stop when we reach the time dependent quantities
                            if 't' in variable.dimensions and name!="t": break
                             
                            # Create the variable 
                            if verbose: print("{0:>15}".format(name), "   ", variable.datatype, variable.dimensions)
                            dst.createVariable(name, variable.datatype, variable.dimensions)
                             
                            # Fill in the variable, but first cut off time positions  
                            if 't' not in variable.dimensions: dst.variables[name][:] = src.variables[name][:] 
                            if 't' in variable.dimensions: dst.variables[name][:] = src.variables[name][indexes] 
                    
                    # Time the calculation of the new time indices 
                    elapsed_time = str(datetime.timedelta(seconds=round(timer.time()-start_step1,0)))
                    
                os.system("mv "+new_file+" "+finished_file); file_step1 = True
                print("\nCompleted step 1 of reducing new netcdf file at every dt = "+str(dt)+".")
                print("   ---> Reduced from "+str(len(time))+" to "+str(len(time_reduced))+" time points.")
                print("   ---> Time reduction of a factor "+str(round(len(time)/len(time_reduced),2))+" which took "+elapsed_time+" hours.")  
                print()
                
            #=================================================================
            #                       Save the phi(t)
            #=================================================================
            
            if (file_step1) and (not file_step2) and (not file_step3) and (not file_step4) and (not file_step5):
            
                # Start timer
                start_step2 = timer.time()
                     
                # New files
                old_file = str(netcdf_file)
                previous_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step1") 
                new_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step2_in_progress")  
                finished_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step2")  
                os.system("cp -r "+previous_file+" "+new_file)
                
                # Print a message
                print("=======================================================")
                print(old_file if ("SCRATCH" not in old_file) else old_file.split("SCRATCH")[-1])
                print(new_file if ("SCRATCH" not in new_file) else new_file.split("SCRATCH")[-1])
                print("=======================================================\n") 
    
                # Constuct a new netcdf file with less data
                with nc.Dataset(old_file, 'r') as src, nc.Dataset(new_file, "r+") as dst:
                     
                    # Save the indexes and time reduction as a pickle
                    with open(netcdf_file.parent / 'indexes.pkl','rb') as f:
                        indexes = pickle.load(f) 
                    with open(netcdf_file.parent / 'time_reduction.pkl','rb') as f:
                        len_times = pickle.load(f)
                             
                    # Copy all file data for variables that are included in the toinclude list
                    for name, variable in src.variables.items():
                        if 't' in variable.dimensions and name!="t": 
                             
                            # Create the variable  
                            if verbose: print("{0:>15}".format(name), "   ", variable.datatype, variable.dimensions)
                            dst.createVariable(name, variable.datatype, variable.dimensions)
                             
                            # Fill in the variable, but first cut off time positions 
                            if 't' in variable.dimensions: dst.variables[name][:] = src.variables[name][indexes] 
                            
                            # Stop after phi_vs_t
                            if name=="phi_vs_t": break
                    
                    # Time the calculation of the new time indices 
                    elapsed_time = str(datetime.timedelta(seconds=round(timer.time()-start_step2,0)))
                
                # Final message
                os.system("mv "+new_file+" "+finished_file); file_step2 = True
                os.system("rm "+previous_file); file_step1 = False
                print("\nCompleted step 2 of reducing new netcdf file at every dt = "+str(dt)+".")
                print("   ---> Reduced from "+str(len_times[0])+" to "+str(len_times[1])+" time points.")
                print("   ---> Time reduction of a factor "+str(round(len_times[0]/len_times[1],2))+" which took "+elapsed_time+" hours.")  
                print()
                
            #=================================================================
            #                       Save phi2_vs_kxky
            #=================================================================
            
            if (not file_step1) and (file_step2) and (not file_step3) and (not file_step4) and (not file_step5):
            
                # Start timer
                start_step3 = timer.time()
                     
                # New files
                old_file = str(netcdf_file)
                previous_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step2") 
                new_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step3_in_progress")  
                finished_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step3")  
                os.system("cp "+previous_file+" "+new_file)
                
                # Check whether we have the quantities versus (kx,ky)
                keys_to_write = ["phi2_vs_kxky", "pflx_kxky", "vflx_kxky", "qflx_kxky"]  
                with nc.Dataset(old_file) as src:
                    keys = src.variables.keys()
                write_step_spectra = False
                for key in keys_to_write:
                    if key in keys:
                        write_step_spectra = True
                
                # Skip if we don't have these variables
                if not write_step_spectra:
                    os.system("mv "+new_file+" "+finished_file)
                    print()
                    print(" ---> Skipping step 3 since neither phi2_vs_kxky nor pflx_kxky are present.")
                    
                # Print a message
                if write_step_spectra:
                    print("=======================================================")
                    print(old_file if ("SCRATCH" not in old_file) else old_file.split("SCRATCH")[-1])
                    print(new_file if ("SCRATCH" not in new_file) else new_file.split("SCRATCH")[-1])
                    print("=======================================================\n") 
        
                    # Constuct a new netcdf file with less data
                    with nc.Dataset(old_file, "r") as src, nc.Dataset(new_file, "r+") as dst:
                     
                        # Save the indexes and time reduction as a pickle
                        with open(netcdf_file.parent / 'indexes.pkl','rb') as f:
                            indexes = pickle.load(f) 
                        with open(netcdf_file.parent / 'time_reduction.pkl','rb') as f:
                            len_times = pickle.load(f)
                                 
                        # Copy all file data for variables that are included in the toinclude list
                        for name, variable in src.variables.items():
                            if name in keys_to_write: 
                                 
                                # Create the variable 
                                if verbose: print("{0:>15}".format(name), "   ", variable.datatype, variable.dimensions)
                                dst.createVariable(name, variable.datatype, variable.dimensions)
                                 
                                # Fill in the variable, but first cut off time positions 
                                dst.variables[name][:] = src.variables[name][indexes] 
                        
                        # Time the calculation of the new time indices 
                        elapsed_time = str(datetime.timedelta(seconds=round(timer.time()-start_step3,0)))
                    
                    # Final message
                    os.system("mv "+new_file+" "+finished_file); file_step3 = True
                    os.system("rm "+previous_file); file_step2 = False
                    print("\nCompleted step 3 of reducing new netcdf file at every dt = "+str(dt)+".")
                    print("   ---> Reduced from "+str(len_times[0])+" to "+str(len_times[1])+" time points.")
                    print("   ---> Time reduction of a factor "+str(round(len_times[0]/len_times[1],2))+" which took "+elapsed_time+" hours.")  
                    print()
                    
            #=================================================================
            #                       Save the moments
            #=================================================================
            
            if (not file_step1) and (not file_step2) and (file_step3) and (not file_step4) and (not file_step5):
                
                # Start timer
                start_step4 = timer.time()
                     
                # New files
                old_file = str(netcdf_file)
                previous_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step3") 
                new_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step4_in_progress")  
                finished_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step4")  
                
                # Check whether we have "density" or "spitzer2"
                keys_to_write = ["density", "upar"]  
                with nc.Dataset(old_file) as src:
                    keys = src.variables.keys()
                write_step_moments = False
                for key in keys_to_write:
                    if key in keys:
                        write_step_moments = True
                
                # Skip if we don't have these variables
                if not write_step_moments:
                    os.system("mv "+previous_file+" "+finished_file); file_step3 = False; file_step4 = True
                    print("\nSkipping step 4 since neither the density nor upar are present.\n")
                    
                # Print a message
                if write_step_moments:
                    os.system("cp "+previous_file+" "+new_file)
                    print("=======================================================")
                    print(old_file if ("SCRATCH" not in old_file) else old_file.split("SCRATCH")[-1])
                    print(new_file if ("SCRATCH" not in new_file) else new_file.split("SCRATCH")[-1])
                    print("=======================================================\n") 
        
                    # Constuct a new netcdf file with less data
                    with nc.Dataset(old_file, "r") as src, nc.Dataset(new_file, "r+") as dst:
                     
                        # Save the indexes and time reduction as a pickle
                        with open(netcdf_file.parent / 'indexes.pkl','rb') as f:
                            indexes = pickle.load(f) 
                        with open(netcdf_file.parent / 'time_reduction.pkl','rb') as f:
                            len_times = pickle.load(f)
                                 
                        # Copy all file data for variables that are included in the toinclude list
                        for name, variable in src.variables.items():
                            if name in keys_to_write: 
                                 
                                # Create the variable 
                                if verbose: print("{0:>15}".format(name), "   ", variable.datatype, variable.dimensions)
                                dst.createVariable(name, variable.datatype, variable.dimensions)
                                 
                                # Fill in the variable, but first cut off time positions 
                                dst.variables[name][:] = src.variables[name][indexes] 
                        
                        # Time the calculation of the new time indices 
                        elapsed_time = str(datetime.timedelta(seconds=round(timer.time()-start_step4,0)))
                    
                    # Final message
                    os.system("mv "+new_file+" "+finished_file); file_step4 = True
                    os.system("rm "+previous_file); file_step3 = False
                    print("\nCompleted step 4 of reducing new netcdf file at every dt = "+str(dt)+".")
                    print("   ---> Reduced from "+str(len_times[0])+" to "+str(len_times[1])+" time points.")
                    print("   ---> Time reduction of a factor "+str(round(len_times[0]/len_times[1],2))+" which took "+elapsed_time+" hours.")  
                    print()
                    
            #=================================================================
            #                       Save the moments
            #=================================================================
            
            if (not file_step1) and (not file_step2) and (not file_step3) and (file_step4) and (not file_step5):
                
                # Start timer
                start_step5 = timer.time()
                     
                # New files
                old_file = str(netcdf_file)
                previous_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step4") 
                new_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step5_in_progress")  
                finished_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step5")  
                
                # Check whether we have "density" or "spitzer2"
                keys_to_write = ["temperature", "spitzer2"]  
                with nc.Dataset(old_file) as src:
                    keys = src.variables.keys()
                write_step_moments = False
                for key in keys_to_write:
                    if key in keys:
                        write_step_moments = True
                
                # Skip if we don't have these variables
                if not write_step_moments:
                    os.system("mv "+previous_file+" "+finished_file); file_step3 = False; file_step4 = True
                    print("\nSkipping step 5 since neither the temperature nor spitzer2 are present.\n")
                    
                # Print a message
                if write_step_moments:
                    os.system("cp "+previous_file+" "+new_file)
                    print("=======================================================")
                    print(old_file if ("SCRATCH" not in old_file) else old_file.split("SCRATCH")[-1])
                    print(new_file if ("SCRATCH" not in new_file) else new_file.split("SCRATCH")[-1])
                    print("=======================================================\n") 
        
                    # Constuct a new netcdf file with less data
                    with nc.Dataset(old_file, "r") as src, nc.Dataset(new_file, "r+") as dst:
                     
                        # Save the indexes and time reduction as a pickle
                        with open(netcdf_file.parent / 'indexes.pkl','rb') as f:
                            indexes = pickle.load(f) 
                        with open(netcdf_file.parent / 'time_reduction.pkl','rb') as f:
                            len_times = pickle.load(f)
                                 
                        # Copy all file data for variables that are included in the toinclude list
                        for name, variable in src.variables.items():
                            if name in keys_to_write: 
                                 
                                # Create the variable 
                                if verbose: print("{0:>15}".format(name), "   ", variable.datatype, variable.dimensions)
                                dst.createVariable(name, variable.datatype, variable.dimensions)
                                 
                                # Fill in the variable, but first cut off time positions 
                                dst.variables[name][:] = src.variables[name][indexes] 
                        
                        # Time the calculation of the new time indices 
                        elapsed_time = str(datetime.timedelta(seconds=round(timer.time()-start_step5,0)))
                    
                    # Final message
                    os.system("mv "+new_file+" "+finished_file); file_step5 = True
                    os.system("rm "+previous_file); file_step4 = False
                    print("\nCompleted step 5 of reducing new netcdf file at every dt = "+str(dt)+".")
                    print("   ---> Reduced from "+str(len_times[0])+" to "+str(len_times[1])+" time points.")
                    print("   ---> Time reduction of a factor "+str(round(len_times[0]/len_times[1],2))+" which took "+elapsed_time+" hours.")  
                    print()

            #=================================================================
            #                       Save the remainder
            #=================================================================
            
            if (not file_step1) and (not file_step2) and (not file_step3) and (not file_step4) and (file_step5):
                
                # Start timer
                start_step6 = timer.time()
                     
                # New files
                old_file = str(netcdf_file)
                previous_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step5") 
                new_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced.step6_in_progress")  
                finished_file = str(netcdf_file).replace(".out.nc", ".out.nc.reduced")  
                os.system("cp "+previous_file+" "+new_file)
                
                # Check whether we have any missing variables
                with nc.Dataset(old_file) as src:
                    keys_old = src.variables.keys()
                with nc.Dataset(new_file) as src:
                    keys_new = src.variables.keys()
                keys_missing = [key for key in keys_old if key not in keys_new]

                # Print a message
                print("=======================================================")
                print(old_file if ("SCRATCH" not in old_file) else old_file.split("SCRATCH")[-1])
                print(new_file if ("SCRATCH" not in new_file) else new_file.split("SCRATCH")[-1])
                print("=======================================================\n") 
    
                # Constuct a new netcdf file with less data
                with nc.Dataset(old_file, "r") as src, nc.Dataset(new_file, "r+") as dst:
                     
                    # Save the indexes and time reduction as a pickle
                    with open(netcdf_file.parent / 'indexes.pkl','rb') as f:
                        indexes = pickle.load(f) 
                    with open(netcdf_file.parent / 'time_reduction.pkl','rb') as f:
                        len_times = pickle.load(f)
                             
                    # Copy all file data for variables that are included in the toinclude list
                    for name, variable in src.variables.items():
                        if name in keys_missing: 
                             
                            # Create the variable 
                            if verbose: print("{0:>15}".format(name), "   ", variable.datatype, variable.dimensions)
                            dst.createVariable(name, variable.datatype, variable.dimensions)
                                 
                            # Fill in the variable, but first cut off time positions 
                            if 't' in variable.dimensions: dst.variables[name][:] = src.variables[name][indexes] 
                            if 't' not in variable.dimensions:  dst.variables[name][:] = src.variables[name][:] 
                    
                    # Time the calculation of the new time indices 
                    elapsed_time = str(datetime.timedelta(seconds=round(timer.time()-start_step6,0)))
                    
                    # Final message
                    os.system("mv "+new_file+" "+finished_file)
                    os.system("rm "+previous_file) 
                    os.system("rm "+str(netcdf_file.parent / 'indexes.pkl'))
                    os.system("rm "+str(netcdf_file.parent / 'time_reduction.pkl'))
                    print("\nCompleted step 6 of reducing new netcdf file at every dt = "+str(dt)+".")
                    print("   ---> Reduced from "+str(len_times[0])+" to "+str(len_times[1])+" time points.")
                    print("   ---> Time reduction of a factor "+str(round(len_times[0]/len_times[1],2))+" which took "+elapsed_time+" hours.")  
                    print()
                    
                    # Final message
                    elapsed_time = str(datetime.timedelta(seconds=round(timer.time()-start,0)))
                    print("=================================================================")
                    print("                 FINISHED THE LAST STEP ")
                    print("=================================================================")
                    print("\nSaved the new netcdf file with data at every dt = "+str(dt)+".")
                    print("   ---> Reduced from "+str(len_times[0])+" to "+str(len_times[1])+" time points.")
                    print("   ---> Time reduction of a factor "+str(round(len_times[0]/len_times[1],2))+" which took "+elapsed_time+" hours.")  
                    print()
                    print()
    return    
        
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    bash = Bash(reduce_sizeNetcdfBigFiles, __doc__) 
    bash.add_option('dt', 'float', 't', '', 'Only keep the time values [0, dt, 2*dt, 3*dt, ...]') 
    reduce_sizeNetcdfBigFiles(**bash.get_arguments())
