

import os
import numpy as np
from stellapy.utils.files.get_filesInFolder import get_filesInFolder   
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile 
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep

################################################################################
#                    DIFFERENCE IN THE SHAPE OF PHI(Z,T)
################################################################################
# This file contains the maximum difference of the shape of phi(z,t) compared 
# to the shape at the final time step. This is a good measure for convergence.
# If this quantity does not reach numerical noise, we might encounter a jumper.
#       phi(z,t)/(sum_z phi(z,t) - phi(z,tend)/(sum_z phi(z,tend)
################################################################################

def write_txtFileForDPhiZVsTime(folder, dt=0.1):  
    
    # Time step
    dt = int(dt) if int(dt)==dt else dt    
    
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 

    # Iterate through the input files 
    for input_file in input_files:
        
        # Processing status
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
        
        # Only write this data for linear simulations 
        if read_linearNonlinearFromInputFile(input_file)[0]: 
            
            try:
        
                # Path of the new file  
                txt_file = input_file.with_suffix(".dt"+str(dt)+".dphiz_vs_t") 
                netcdf_file = input_file.with_suffix(".out.nc")
                
                # Check whether the txt file is older than the simulation
                outputFileHasChanged = False
                if os.path.isfile(txt_file): 
                    if netcdf_file.stat().st_mtime > txt_file.stat().st_mtime:
                        outputFileHasChanged = True 
                
                # Notify that the file already existed   
                if os.path.isfile(txt_file) and not outputFileHasChanged:
                    print(status+"The dphiz(t) file already exists:", txt_file.parent.name+"/"+txt_file.name)
                    continue
                
                # Create the text file for dphiz(t)
                elif (not os.path.isfile(txt_file)) or (os.path.isfile(txt_file) and outputFileHasChanged):    
                    
                    # Read the potential versus time from the output file
                    vec_time, dphiz_vs_t = read_dphizVsTime(dt, netcdf_file)  
                         
                # Write the dphiz(t) data to a text file  
                dphiz_file = open(txt_file,'w') 
                dphiz_file.write("       time               |dphi(z)|\n")
                for i in range(len(dphiz_vs_t)):    
                    np.savetxt(dphiz_file, [[vec_time[i], dphiz_vs_t[i]]], fmt='%16.8e', delimiter="   ") 
                dphiz_file.close() 
                        
                # Notify that we finished creating the file
                if outputFileHasChanged: print(status+"   ---> The dphiz(t) file is updated as " +  txt_file.parent.name+"/"+txt_file.name)   
                if not outputFileHasChanged:  print(status+"   ---> dphiz phi(t) file is saved as " +  txt_file.parent.name+"/"+txt_file.name)  

            except:
                print(status+"   ---> The dphiz(t) file couldn't be written for " +  txt_file.parent.name+"/"+txt_file.name)   
    return  

#---------------------------------------------
def read_dphizVsTime(dt, netcdf_file):
    
    # Read the potential from the output file
    netcdf_data = read_outputFile(netcdf_file)   
    vec_time = read_netcdfVariables('vec_time', netcdf_data, netcdf_path=netcdf_file) 
    phi_vs_tzkxkyri = read_netcdfVariables('phi_vs_tzkxkyri', netcdf_data, netcdf_path=netcdf_file) 
    netcdf_data.close()
    
    # Get the data at every <dt> timestep and assume we have one mode per file
    indices = get_indicesAtFixedStep(vec_time, dt)  
    if "/mnt/lustre" in os.path.abspath(__file__): phi_vs_tz = np.apply_along_axis(lambda args: [complex(*args)], 2, phi_vs_tzkxkyri[indices,:,0,0,:])[:,:,0]
    else: phi_vs_tz = phi_vs_tzkxkyri[indices,:,0,0,0] + 1j*phi_vs_tzkxkyri[indices,:,0,0,1]
    vec_time = vec_time[indices]
    
    # At a certain point a NaN might occur in the (z,kx,ky) grid
    # Remove this data, as well as data with infinite values
    phi_vs_tz[phi_vs_tz==np.inf] = np.NaN 
    for it in range(len(vec_time)): vec_time[it] = np.nan if np.any(np.isnan(phi_vs_tz[it,:])) else vec_time[it] 
    phi_vs_tz = phi_vs_tz[~np.isnan(vec_time),:] 
    vec_time = vec_time[~np.isnan(vec_time)]
                    
    # Find the maximum difference of the shape of phi(z,t) compared to phi(z,tend)
    phi_vs_tz = phi_vs_tz/(np.sum(phi_vs_tz, axis=1)[:, np.newaxis])  
    dphiz_vs_tz = np.abs(phi_vs_tz - phi_vs_tz[-1,:]) 
    dphiz_vs_t = np.nanmax(dphiz_vs_tz,axis=1)
    return vec_time, dphiz_vs_t
    
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    import pathlib
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")
    write_txtFileForDPhiZVsTime(folder)
    
