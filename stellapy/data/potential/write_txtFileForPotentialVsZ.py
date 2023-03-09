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

#===============================================================================
#                  PARALLEL MODE STRUCTURE OF THE POTENTIAL
#===============================================================================

def write_txtFileForPotentialVsZ(folder):  
    
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 

    # Iterate through the input files 
    for input_file in input_files: 
        
        # Processing status
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
        
        # Only write the final fields file for linear simulations  
        nonlinear = read_linearNonlinearFromInputFile(input_file)[1]
        if nonlinear: continue
                
        # Path of the new file   
        netcdf_path = input_file.with_suffix(".out.nc") 
        potential_path = input_file.with_suffix(".phi_vs_z")  
        
        # Check whether the txt file is older than the simulation
        outputFileHasChanged = False
        if os.path.isfile(potential_path):  
            if netcdf_path.stat().st_mtime > potential_path.stat().st_mtime:
                outputFileHasChanged = True 
    
        # Notify that the file already existed   
        if os.path.isfile(potential_path) and not outputFileHasChanged:
            print(status+"The phi(z) file already exists:", potential_path.parent.name+"/"+potential_path.name)
            continue 
    
        # Read the number of modes simulated in the input file
        nakx, naky = read_numberOfModesFromInputFile(input_file); nakxnaky = nakx*naky
        
        # Read the potential versus time from the output file
        phi_vs_zkxky, phi2_vs_zkxky = read_potentialVsZ(netcdf_path) 
    
        # h5 file for multiple modes per simulation
        if nakxnaky>1:  
            
            # Print omega to an h5 file
            with h5py.File(potential_path, 'w') as h5_file:   
                h5_file.create_dataset("phi_vs_zkxky", data=phi_vs_zkxky) 
                h5_file.create_dataset("phi2_vs_zkxky", data=phi2_vs_zkxky)  
            
        # txt file for 1 mode per simulation
        elif nakxnaky==1:
            
            # Data to be written
            header = " "*7+"|phi|^2"+" "*12+"Re(phi)"+" "*12+"Im(phi)"+"\n"
            potential_data = np.concatenate((phi2_vs_zkxky[:,0,0,np.newaxis], phi_vs_zkxky[:,0,0,np.newaxis].real, phi_vs_zkxky[:,0,0,np.newaxis].imag), axis=1)
             
            # Write the potential data to a text file  
            phi_file = open(potential_path,'w') 
            phi_file.write(header) 
            np.savetxt(phi_file, potential_data, fmt='%16.8e', delimiter="   ") 
            phi_file.close() 
            
        # Notify that we finished creating the file
        print(status+"   ---> The phi(z) file is saved as " + potential_path.parent.name+"/"+potential_path.name)    
        
    return  

#---------------------------------------------
def read_potentialVsZ(netcdf_path):
                
    # Read the potential from the output file
    netcdf_data = read_outputFile(netcdf_path)  
    vec_time = read_netcdfVariables('vec_time', netcdf_data) 
    time_indices = get_indicesAtFixedStep(vec_time, dx=1)  
    phi_vs_tzkxkyri = read_netcdfVariables('phi_vs_tzkxkyri', netcdf_data, time_indices) 
    vec_time = vec_time[time_indices]
    netcdf_data.close()
    
    # Construct the potential from its real and imaginary part (code below fixes bug on xula)
    if "/mnt/lustre" in os.path.abspath(__file__): phi_vs_tzkxky = np.apply_along_axis(lambda args: [complex(*args)], 2, phi_vs_tzkxkyri)[:,:,:,:,0]
    else: phi_vs_tzkxky = phi_vs_tzkxkyri[:,:,:,:,0] + 1j*phi_vs_tzkxkyri[:,:,:,:,1]
    
    # Get the last time step without nans for each mode (kx,ky)
    phi_vs_zkxky = np.ones(np.shape(phi_vs_tzkxky)[1:], dtype="complex")*np.nan
    dim_kx = np.shape(phi_vs_zkxky)[1]
    dim_ky = np.shape(phi_vs_zkxky)[2]
    nakxnaky = dim_kx*dim_ky; count = 0
    for it in range(len(vec_time)-1, -1, -1):
        for ikx in range(dim_kx):
            for iky in range(dim_ky):
                if not np.any(np.isnan(phi_vs_tzkxky[it,:,ikx,iky])):
                    phi_vs_zkxky[:,ikx,iky] = phi_vs_tzkxky[it,:,ikx,iky]; count += 1
                if count==nakxnaky:
                    break
    
    # Get the potential squared
    phi2_vs_zkxky = np.abs(phi_vs_zkxky)**2 
    return phi_vs_zkxky, phi2_vs_zkxky
