

import os
import numpy as np
from datetime import datetime, timedelta
from stellapy.utils.files.get_filesInFolder import get_filesInFolder     
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables

################################################################################
#                        TIME EVOLUTION OF THE POTENTIAL
################################################################################ 

def write_txtFileForPotentialVsZ(folder):  
    
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc')) or os.path.isfile(i.with_suffix('.final_fields'))]
    if input_files==[]: return 

    # Iterate through the input files 
    for input_file in input_files: 
        
        # Processing status
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
        
        # Only write the final fields file for linear simulations 
        if read_linearNonlinearFromInputFile(input_file)[0]:
                
            # Path of the new file   
            netcdf_file = input_file.with_suffix(".out.nc")
            fields_file = input_file.with_suffix(".final_fields")
            potential_file = input_file.with_suffix(".phi_vs_z")  
            
            # Check whether the txt file is older than the simulation
            simulationWasntFinished = False
            if os.path.isfile(potential_file): 
                time_txtFile = datetime.fromtimestamp(potential_file.stat().st_mtime)
                if os.path.isfile(netcdf_file): time_simulation = datetime.fromtimestamp(netcdf_file.stat().st_mtime) 
                if os.path.isfile(potential_file): time_simulation = datetime.fromtimestamp(potential_file.stat().st_mtime) 
                if time_simulation>time_txtFile+timedelta(minutes=5):
                    simulationWasntFinished = True 
            
            # Create text files for dphiz(t)
            if not os.path.isfile(potential_file) or simulationWasntFinished:   
                
                # Read the final fields from the netcdf file
                if os.path.isfile(netcdf_file):
                    
                    # Read the potential from the output file
                    netcdf_data = read_outputFile(netcdf_file)  
                    vec_time = read_netcdfVariables('vec_time', netcdf_data) 
                    phi_vs_tzkxkyri = read_netcdfVariables('phi_vs_tzkxkyri', netcdf_data) 
                    if "/mnt/lustre" in os.path.abspath(__file__): phi_vs_tz = np.apply_along_axis(lambda args: [complex(*args)], 2, phi_vs_tzkxkyri[:,:,0,0,:])[:,:,0]
                    else: phi_vs_tz = phi_vs_tzkxkyri[:,:,0,0,0] + 1j*phi_vs_tzkxkyri[:,:,0,0,1]
                    netcdf_data.close()
                    
                    # Get the last time step without nans
                    for it in range(len(vec_time)-1, -1, -1):
                        if not np.any(np.isnan(phi_vs_tz[it,:])):
                            phi_vs_z = phi_vs_tz[it,:] 
                            break 
                    
                    # Get the potential squared
                    phi2_vs_z = np.abs(phi_vs_z)**2 
                    
                # For old simulations, hope that the field file doesn't contain nans
                elif os.path.isfile(fields_file): 
                    try:    data = np.loadtxt(fields_file, dtype='float').reshape(-1, 10) # Old code has 10 columns
                    except: data = np.loadtxt(fields_file, dtype='float').reshape(-1, 11) # New code has 11 columns
                    phi_vs_z = data[:,4] + 1j*data[:,5]
                    phi2_vs_z = np.abs(phi_vs_z)**2 
             
                # Write the potential data to a text file  
                phi_file = open(potential_file,'w') 
                phi_file.write(" "*7+"|phi|^2"+" "*12+"Re(phi)"+" "*12+"Im(phi)"+"\n")
                for i in range(len(phi2_vs_z)): 
                    row = [[phi2_vs_z[i], phi_vs_z[i].real, phi_vs_z[i].imag]]   
                    np.savetxt(phi_file, row, fmt='%16.8e', delimiter="   ") 
                phi_file.close() 
                    
                # Notify that we finished creating the file
                print(status+"   ---> The phi(z) file is saved as " + potential_file.parent.name+"/"+potential_file.name)    
             
            # Notify that the file already existed   
            else:
                print(status+"The phi(z) file already exists:", potential_file.parent.name+"/"+potential_file.name)
            
    return  

################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    import pathlib
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")
    write_txtFileForPotentialVsZ(folder)
    
