
import os
import numpy as np  
from stellapy.utils.files.get_filesInFolder import get_filesInFolder
from stellapy.data.input.read_inputFile import read_nspecFromInputFile  

#===============================================================================
#                      REPLACE *.FLUXES BY *.DT1.FLUXES
#===============================================================================
# Make sure to run this bash command on marconi after every simulation.
#     > reduce_sizeFluxes -t 1
#===============================================================================

def write_txtFileForFluxesVsTime(folder, dt=1, verbose=False):
    
    # Time step
    dt = int(dt) if int(dt)==dt else dt   
    
    # Suffix of the new file
    suffix = ".dt"+str(dt)+".fluxes_vs_t" 
    
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.fluxes'))]
    if input_files==[]: return 
    
    # Go through the input files
    for input_file in input_files:  
        
        # Processing status
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   " 
        
        try:
                 
            # Check whether the reduced file is older than the netcdf file
            already_written = False
            if os.path.isfile(input_file.with_suffix(suffix)):
                time_h5File = input_file.with_suffix(suffix).stat().st_mtime
                time_fluxesFile = input_file.with_suffix(".fluxes").stat().st_mtime
                if time_h5File>time_fluxesFile:
                    already_written = True
    
            # Only write it if we have to 
            if already_written==False:
    
                # Get dim_species 
                dim_species = read_nspecFromInputFile(input_file) 
    
                # Read the fluxes file 
                flux_file = input_file.with_suffix(".fluxes") 
                flux_data = np.loadtxt(flux_file, dtype='float').reshape(-1, 1+3*dim_species) 
            
                # Header for the file 
                header = "      #time            "
                for s in range(dim_species): header += "pflux (s="+str(s)+")        "    
                for s in range(dim_species): header += "vflux (s="+str(s)+")        "    
                for s in range(dim_species): header += "qflux (s="+str(s)+")        "    
                header += "\n"
                
                # Create the new file and add the header
                fluxes_file = open(input_file.with_suffix(suffix),'w') 
                fluxes_file.write(header)
                np.savetxt(fluxes_file, [flux_data[0,:]], fmt='%16.8e', delimiter="   ")
                
                # Keep track of the current time
                current_time = 0; count = 0
                
                # Save the fluxes data  
                for i,line in enumerate(flux_data):
                    if flux_data[i,0] >= current_time+dt:
                        np.savetxt(fluxes_file, [line], fmt='%16.8e', delimiter="   ")
                        current_time = current_time+dt; count += 1
                        
                # Close the file
                fluxes_file.close()
                
                # Write a message
                print(status+"   ---> The fluxes(t) file is saved as " + input_file.with_suffix(suffix).name)  
                
            else: 
                if verbose: print(status+"The fluxes(t) file already exists:", input_file.with_suffix(suffix).parent.name+"/"+input_file.with_suffix(suffix).name)
        except:
            print(status+"Something went wrong for:", input_file.with_suffix(suffix).parent.parent.name+"/"+input_file.with_suffix(suffix).parent.name+"/"+input_file.with_suffix(suffix).name)
            import sys; sys.exit()
    return
