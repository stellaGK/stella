
"""

#================================================================================
# Write the maximum difference of the shape phi(z,t) compared with phi(z,tlast) #
#================================================================================

This txt file contains the maximum difference of the shape of phi(z,t) compared 
to the shape at the final time step. For a linear converged simulation, the shape
of phi(z,t), or thus phi(z,t)/sum_z(phi(z,t)) is a constant. Therefore, the 
shape of phi(z,t) is a good measure for convergence. If this quantity does not 
reach numerical noise (around 1e-19), we might encounter a jumper: this is a 
mode (kx,ky) for which the growth rate seems converged, however, after a certain
time the growth rate "jumps" to a different value, since numerical convergence
had not yet been achieved.

The quantity which is written is:
      np.max(phi(z,t)/(sum_z phi(z,t)) - phi(z,tend)/(sum_z phi(z,tend)))

If multiple modes were launched together for a linear flux tube simulation, then
the file is an h5 file instead, since a txt file would be too long.

Hanne Thienpondt
20/10/2022

"""

#!/usr/bin/python3
import h5py
import os, sys
import pathlib
import numpy as np

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.input.read_inputFile import read_numberOfModesFromInputFile 
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep
from stellapy.utils.files.get_filesInFolder import get_filesInFolder   

#================================================================================
# Write the maximum difference of the shape phi(z,t) compared with phi(z,tlast) #
#================================================================================

def write_txtFileForDPhiZVsTime(folder, dt=1):  
    
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
                dphiz_path = input_file.with_suffix(".dt"+str(dt)+".dphiz_vs_t") 
                netcdf_path = input_file.with_suffix(".out.nc")
                
                # Check whether the txt file is older than the simulation
                outputFileHasChanged = False
                if os.path.isfile(dphiz_path): 
                    if netcdf_path.stat().st_mtime > dphiz_path.stat().st_mtime:
                        outputFileHasChanged = True 
                
                # Notify that the file already existed   
                if os.path.isfile(dphiz_path) and not outputFileHasChanged:
                    print(status+"The dphiz(t) file already exists:", dphiz_path.parent.name+"/"+dphiz_path.name)
                    continue
                
                # Read the number of modes simulated in the input file
                nakx, naky = read_numberOfModesFromInputFile(input_file); nakxnaky = nakx*naky    
                
                # Read the potential versus time from the output file
                vec_time, dphiz_vs_tkxky = read_dphizVsTime(dt, netcdf_path)  
                
                # h5 file for multiple modes per simulation
                if nakxnaky>1:  
            
                    # Print dphiz to an h5 file
                    with h5py.File(dphiz_path, 'w') as h5_file:  
                        h5_file.create_dataset("vec_time", data=vec_time) 
                        h5_file.create_dataset("dphiz_vs_tkxky", data=dphiz_vs_tkxky)    
                                
                # txt file for 1 mode per simulation
                elif nakxnaky==1:
                
                    # Write the dphiz(t) data to a text file  
                    dphiz_file = open(dphiz_path,'w') 
                    dphiz_file.write("       time               |dphi(z)|\n")
                    for i in range(len(vec_time)):    
                        np.savetxt(dphiz_file, [[vec_time[i], dphiz_vs_tkxky[i,0,0]]], fmt='%16.8e', delimiter="   ") 
                    dphiz_file.close() 
                        
                # Notify that we finished creating the file
                print(status+"   ---> dphiz phi(t) file is saved as " +  dphiz_path.parent.name+"/"+dphiz_path.name)  

            except:
                print(status+"   ---> The dphiz(t) file couldn't be written for " +  dphiz_path.parent.name+"/"+dphiz_path.name)   
    return  

#---------------------------------------------
def read_dphizVsTime(dt, netcdf_path):
    
    # Read the potential from the output file
    netcdf_data = read_outputFile(netcdf_path)   
    vec_time = read_netcdfVariables('vec_time', netcdf_data, netcdf_path=netcdf_path) 
    time_indices = get_indicesAtFixedStep(vec_time, dt)  
    phi_vs_tzkxkyri = read_netcdfVariables('phi_vs_tzkxkyri', netcdf_data, time_indices, netcdf_path=netcdf_path) 
    vec_time = vec_time[time_indices]
    netcdf_data.close()
    
    # Construct the potential from its real and imaginary part (issue on xula solved with code below)
    if "/mnt/lustre" in os.path.abspath(__file__): phi_vs_tkxkyz = np.apply_along_axis(lambda args: [complex(*args)], 2, phi_vs_tzkxkyri)[:,:,:,:,0]
    else: phi_vs_tkxkyz = phi_vs_tzkxkyri[:,:,:,:,0] + 1j*phi_vs_tzkxkyri[:,:,:,:,1]
    
    # At a certain point a NaN might occur in the (z,kx,ky) grid
    # Remove this data, as well as data with infinite values
    phi_vs_tkxkyz[phi_vs_tkxkyz==np.inf] = np.NaN 
    for it in range(len(vec_time)): vec_time[it] = np.nan if np.any(np.isnan(phi_vs_tkxkyz[it,:,:,:])) else vec_time[it] 
    phi_vs_tkxkyz = phi_vs_tkxkyz[~np.isnan(vec_time),:] 
    vec_time = vec_time[~np.isnan(vec_time)]
                    
    # Find the maximum difference of the shape of phi(z,t) compared to phi(z,tend)
    phi_vs_tkxkyz = phi_vs_tkxkyz/(np.sum(phi_vs_tkxkyz, axis=1)[:, np.newaxis])  
    dphiz_vs_tkxkyz = np.abs(phi_vs_tkxkyz - phi_vs_tkxkyz[-1,:,:,:]) 
    dphiz_vs_tkxky = np.nanmax(dphiz_vs_tkxkyz,axis=1)
    return vec_time, dphiz_vs_tkxky

