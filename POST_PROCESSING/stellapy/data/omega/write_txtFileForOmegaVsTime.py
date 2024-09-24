""" 

#===============================================================================
#                      Write *.dt1.omega based on *.omega
#=============================================================================== 

Write omega(t) and gamma(t) to a text file if the simulation contains one
(kx,ky) mode, or to an h5 file if there are multiple (kx,ky) modes.
The purpose is to reduce the number of time points in the input.omega file.
     >> write_stellapyDataFiles -s omega -t 1

Hanne Thienpondt
20/01/2023

"""
 
#!/usr/bin/python3
import pathlib
import numpy as np   
import os, sys, h5py

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.utils.print_status import print_status, print_statusSkipped, get_nestingDepth
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables 
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.input.read_inputFile import read_numberOfModesFromInputFile 
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep  
from stellapy.utils.files.get_filesInFolder import get_filesInFolder

#===============================================================================
#                      Write *.dt1.omega based on *.omega
#=============================================================================== 

def write_txtFileForOmegaVsTime(folder, dt=1): 
    
    # Time step
    dt = int(dt) if int(dt)==dt else dt    
    
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = sorted([i for i in input_files if os.path.isfile(i.with_suffix('.omega'))])
    nesting_depth = get_nestingDepth(input_files)
    if input_files==[]: return  
    
    # Go through the input files  
    for input_file in input_files: 
        
        # Only write the growth rate data for linear simulations 
        nonlinear = read_linearNonlinearFromInputFile(input_file)[1]
        if nonlinear: continue 
        
        # Path of the new omega file and the netcdf file    
        omega_path = input_file.with_suffix(".dt"+str(dt)+".omega_vs_t")   
        netcdf_path = input_file.with_suffix(".out.nc")  
        
        # Check whether the reduced file has already been created 
        if os.path.isfile(omega_path): 
            time_omegaFile = omega_path.stat().st_mtime
            time_stellaOmegaFile = input_file.with_suffix(".omega").stat().st_mtime
            if time_omegaFile>time_stellaOmegaFile:
                continue  
        
        # Read the number of modes simulated in the input file
        nakx, naky = read_numberOfModesFromInputFile(input_file); nakxnaky = nakx*naky
        
        # h5 file for multiple modes per simulation
        if nakxnaky>1:
            
            # Read the fluxes from the output file
            netcdf_data = read_outputFile(netcdf_path)  
            vec_time = read_netcdfVariables('vec_time', netcdf_data)
            indices = get_indicesAtFixedStep(vec_time, dt)
            omega_vs_tkxky = read_netcdfVariables("omega_vs_tkxky", netcdf_data, indices) 
            vec_time = vec_time[indices]
            netcdf_data.close()   
            
            # Print omega to an h5 file
            with h5py.File(omega_path, 'w') as h5_file:  
                h5_file.create_dataset("vec_time", data=vec_time) 
                h5_file.create_dataset("omega_vs_tkxky", data=omega_vs_tkxky) 
            print_status(omega_path, "omega(t)", input_file, input_files, nesting_depth)
            
        # txt file for 1 mode per simulation
        elif nakxnaky==1: 
            
            # Check how many modes are written 
            dim_kx, dim_ky = read_numberOfModesFromInputFile(input_file)
            multiple_modes_per_file = True if (dim_kx*dim_ky>1) else False

            # Read the omega file and return omega_data[kx,ky,time,{0:kx,1:ky,2:time,3:omega,4:gamma}]
            omega_data = read_omega(input_file) 
            
            # If the simulation didn't work the omega time axis is empty
            if np.shape(omega_data)[2]==0:
                print_statusSkipped(omega_path, "Omega file was empty, skipping.", input_file, input_files, nesting_depth) 
                continue
            
            # Find the indices at which to keep the data (time axis) 
            indices = [0]; current_time = omega_data[0,0,0,2]; time_dim = len(omega_data[0,0,:,2])
            for it in range(time_dim):
                if omega_data[0,0,it,2] >= current_time+dt:
                    current_time = omega_data[0,0,it,2]
                    indices.append(it) 

            # Data and header if there is one mode in the file
            if not multiple_modes_per_file: 
                header = "      #time              omega            gamma  \n" 
                omega_data = omega_data[0,0,indices,2:]
            
            # Data and header if there are multiple modes in the file
            if multiple_modes_per_file: 
                header = "      #time               omega              gamma               kx                ky  \n"
                header = "      # kx                  ky               time               omega              gamma               \n"
                omega_data = omega_data[:,:,indices,:].reshape(-1, 5)

            # Create the new file and add the header
            omega_file = open(omega_path,'w') 
            omega_file.write(header)
            np.savetxt(omega_file, omega_data, fmt='%16.8e', delimiter="   ") 
            omega_file.close()
            print_status(omega_path, "omega(t)", input_file, input_files, nesting_depth)
            
    return 

#------------------------------
def read_omega(input_file):
    ''' Read the *.omega file: [time,  ky,  kx,  Re(om),  Im(om),  Re(omavg), Im(omavg)]
    and assume we only have one mode (kx, ky) for this reduction!'''
    
    # Find the omega file corresponding to the input file
    omega_file = input_file.with_suffix(".omega")
    
    # Read the data in the omega file
    dim_kx, dim_ky = read_numberOfModesFromInputFile(input_file)
    omega_data = np.loadtxt(omega_file, dtype='float').reshape(-1, dim_kx, dim_ky, 7)[:,:,:,:]  
    
    # In the fourth axis we have [time, ky, kx, omega, gamma, omega_avg, gamma_avg], trow away the last two columns
    omega_data = omega_data[:,:,:,:5]
    
    # Move the data so that in the fourth columns we have [kx, ky, time, omega, gamma]
    omega_data = np.concatenate((omega_data[:,:,:,np.newaxis,2], omega_data[:,:,:,np.newaxis,1], omega_data[:,:,:,np.newaxis,0], omega_data[:,:,:,np.newaxis,3], omega_data[:,:,:,np.newaxis,4]), axis=3)
    omega_data = np.swapaxes(omega_data,0,2)
    omega_data = np.swapaxes(omega_data,0,1)

    # Return the data
    return omega_data

