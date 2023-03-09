
"""

#===============================================================================
#                  Write the time evolution of the potential                   #
#=============================================================================== 

Create a text file that contains the following columns:
    Linear simulations:       {time, |phi|^2, Re(phi), Im(phi)}
    Nonlinear simulations:    {time, |phi|^2, zonal|phi|^2, nozonal|phi|^2, Re(phi), Im(phi)}
     
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
from stellapy.utils.files.get_filesInFolder import get_filesInFolder    
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile
from stellapy.data.input.read_inputFile import read_numberOfModesFromInputFile 
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry 

#===============================================================================
#                        TIME EVOLUTION OF THE POTENTIAL
#===============================================================================

def write_txtFileForPotentialVsTime(folder, dt=1):   
    
    # Segmentation fault (core dumped) on xula
    if "/mnt/lustre" in os.path.abspath(__file__): return 
    
    # Time step
    dt = int(dt) if int(dt)==dt else dt    
    
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 
    
    # Iterate through the input files 
    for input_file in input_files: 
        
        # Check progress
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
    
        # Check whether we have a linear or nonlinear simulation   
        nonlinear = read_linearNonlinearFromInputFile(input_files[0])[1]
        
        # Depending on whether the simulation is linear or nonlinear write different data
        try:
            if not nonlinear: write_txtFileForPotentialVsTimeLinearSimulations(input_file, dt, status)
            if nonlinear: write_txtFileForPotentialVsTimeNonlinearSimulations(input_file, dt, status)
        except:
            print(status+"   ---> Something went wrong writing phi(t) for " +  input_file.parent.parent.name+"/"+input_file.parent.name+"/"+input_file.name)   
            sys.exit()
        
    return 

#===============================================================================
#                            NONLINEAR SIMULATIONS                             #
#===============================================================================
 
def write_txtFileForPotentialVsTimeNonlinearSimulations(input_file, dt, status):    
    
    # Path of the new file  
    potential_path = input_file.with_suffix(".dt"+str(dt)+".phi_vs_t") 
    netcdf_path = input_file.with_suffix(".out.nc") 
    
    # Check whether the txt file is older than the simulation
    outputFileHasChanged = False
    if os.path.isfile(potential_path):  
        if netcdf_path.stat().st_mtime > potential_path.stat().st_mtime:
            outputFileHasChanged = True 
    
    # Notify that the file already existed   
    if os.path.isfile(potential_path) and not outputFileHasChanged:
        print(status+"The phi(t) file already exists:", potential_path.parent.name+"/"+potential_path.name)
        return
    
    # If the output file changed, then append to the txt file
    elif os.path.isfile(potential_path) and outputFileHasChanged:
        
        # Check whether the output file (*.out.nc) contains extra time points
        netcdf_data  = read_outputFile(netcdf_path)   
        vec_time     = read_netcdfVariables('vec_time', netcdf_data, netcdf_path=netcdf_path) 
        netcdf_data.close()  
        
        # Edit the time vector in the output file
        indices  = get_indicesAtFixedStep(vec_time, dt)
        vec_time = vec_time[indices]
        vec_time = [round(n, 8) for n in vec_time]
        tlast_outputfile = vec_time[-1]  

        # Read the potential versus time from the txt file
        data = np.loadtxt(potential_path,skiprows=1,dtype='float').reshape(-1, 6)
        vec_time_txtfile  = data[:,0]
        phi2_vs_t_txtfile = data[:,1]
        phi2_vs_t_zonal_txtfile = data[:,2]
        phi2_vs_t_nozonal_txtfile = data[:,3]
        phi_vs_t_txtfile  = data[:,4] + 1j*data[:,5] 
        tlast_txtfile = vec_time_txtfile[-1]
        
        # If the output file (*.out.nc) contains extra time steps, append to the text file 
        if tlast_outputfile > tlast_txtfile: 
            index_tlast = np.argwhere(tlast_txtfile < vec_time)[0][0] 
            phi_vs_t_outputFile, phi2_vs_t_outputFile, phi2_vs_t_zonal_outputFile, \
            phi2_vs_t_nozonal_outputFile, vec_time_outputFile = read_potentialVsTimeNonlinearSimulations(dt, input_file, netcdf_path)  
            vec_time = np.append(vec_time_txtfile, vec_time_outputFile[index_tlast:], axis=0)  
            phi_vs_t = np.append(phi_vs_t_txtfile, phi_vs_t_outputFile[index_tlast:], axis=0)  
            phi2_vs_t = np.append(phi2_vs_t_txtfile, phi2_vs_t_outputFile[index_tlast:], axis=0) 
            phi2_vs_t_zonal = np.append(phi2_vs_t_zonal_txtfile, phi2_vs_t_zonal_outputFile[index_tlast:], axis=0)
            phi2_vs_t_nozonal = np.append(phi2_vs_t_nozonal_txtfile, phi2_vs_t_nozonal_outputFile[index_tlast:], axis=0) 
            
        # If the output file has the same time steps, touch the text file
        elif tlast_outputfile <= tlast_txtfile:
            print(status+"The phi(t) file is up to date:", potential_path.parent.name+"/"+potential_path.name)  
            os.system("touch "+str(potential_path))
            return
        
    # Create the text file for dphiz(t)
    elif not os.path.isfile(potential_path):    
        
        # Read the potential versus time from the output file
        phi_vs_t, phi2_vs_t, phi2_vs_t_zonal, phi2_vs_t_nozonal, vec_time = read_potentialVsTimeNonlinearSimulations(dt, input_file, netcdf_path) 
             
    # Write the potential data to a text file  
    potential_file = open(potential_path,'w') 
    potential_file.write(" "*7+"time"+" "*13+"|phi|^2"+" "*10+"zonal|phi|^2"+" "*6+"nozonal|phi|^2"+" "*8+"Re(phi)"+" "*12+"Im(phi)"+"\n")
    for i in range(len(phi_vs_t)): 
        row = [[vec_time[i], phi2_vs_t[i], phi2_vs_t_zonal[i], phi2_vs_t_nozonal[i], phi_vs_t[i].real, phi_vs_t[i].imag]]   
        np.savetxt(potential_file, row, fmt='%16.8e', delimiter="   ") 
    potential_file.close() 
        
    # Notify that we finished creating the file
    if outputFileHasChanged: print(status+"   ---> The phi(t) file is updated as " +  potential_path.parent.name+"/"+potential_path.name)   
    if not outputFileHasChanged:  print(status+"   ---> The phi(t) file is saved as " +  potential_path.parent.name+"/"+potential_path.name)  
    return  

#---------------------------------------------
def read_potentialVsTimeNonlinearSimulations(dt, input_file, netcdf_path):

    # Read the geometry data in the output file
    vmec_filename = read_vmecFileNameFromInputFile(input_file) 
    path = create_dummyPathObject(input_file, vmec_filename)
    geometry = read_outputFileForGeometry(path) 
    dl_over_B = geometry["dl_over_B"]  
    
    # Read the potential from the output file
    netcdf_data = read_outputFile(netcdf_path)  
    vec_time = read_netcdfVariables('vec_time', netcdf_data) 
    phi_vs_tzkxkyri = read_netcdfVariables('phi_vs_tzkxkyri', netcdf_data) 
    netcdf_data.close()
    
    # Get the data at every <dt> timestep 
    indices = get_indicesAtFixedStep(vec_time, dt) 
    phi_vs_tzkxky = phi_vs_tzkxkyri[indices,:,:,:,0] + 1j*phi_vs_tzkxkyri[indices,:,:,:,1] 
    vec_time = vec_time[indices]
    
    # Sum away the (kx,ky) modes and the (z) dimension  
    phi_vs_tz = np.sum(phi_vs_tzkxky[:,:,:,0],axis=2) + 2*np.sum(phi_vs_tzkxky[:,:,:,1:],axis=(2,3)) 
    phi_vs_t = np.sum(phi_vs_tz[:,:]*dl_over_B[np.newaxis,:], axis=1) 
    
    # Carefull: for phi2 we need to first take abs(phi)**2 and then sum away the dimensions (kx,ky,z)!
    phi2_vs_tzkxky = np.abs(phi_vs_tzkxky)**2
    phi2_vs_tz_zonal = np.sum(phi2_vs_tzkxky[:,:,:,0],axis=2)
    phi2_vs_tz_nozonal = 2*np.sum(phi2_vs_tzkxky[:,:,:,1:],axis=(2,3)) 
    phi2_vs_tz = phi2_vs_tz_zonal + phi2_vs_tz_nozonal 
    phi2_vs_t = np.sum(phi2_vs_tz[:,:]*dl_over_B[np.newaxis,:], axis=1)
    phi2_vs_t_zonal = np.sum(phi2_vs_tz_zonal[:,:]*dl_over_B[np.newaxis,:], axis=1)
    phi2_vs_t_nozonal = np.sum(phi2_vs_tz_nozonal[:,:]*dl_over_B[np.newaxis,:], axis=1)
            
    # Return the data
    return phi_vs_t, phi2_vs_t, phi2_vs_t_zonal, phi2_vs_t_nozonal, vec_time

#===============================================================================
#                             LINEAR SIMULATIONS                               #
#===============================================================================

def write_txtFileForPotentialVsTimeLinearSimulations(input_file, dt, status):   
    
    # Path of the new file  
    potential_path = input_file.with_suffix(".dt"+str(dt)+".phi_vs_t") 
    netcdf_path = input_file.with_suffix(".out.nc") 
    
    # Check whether the txt file is older than the simulation
    outputFileHasChanged = False
    if os.path.isfile(potential_path):  
        if netcdf_path.stat().st_mtime > potential_path.stat().st_mtime:
            outputFileHasChanged = True 
    
    # Notify that the file already existed   
    if os.path.isfile(potential_path) and not outputFileHasChanged:
        print(status+"The phi(t) file already exists:", potential_path.parent.name+"/"+potential_path.name)
        return
    
    # Read the number of modes simulated in the input file
    nakx, naky = read_numberOfModesFromInputFile(input_file); nakxnaky = nakx*naky
        
    # Read the potential versus time from the output file
    vec_time, phi_vs_tkxky, phi2_vs_tkxky = read_potentialVsTimeLinearSimulations(dt, input_file, netcdf_path) 
    
    # h5 file for multiple modes per simulation
    if nakxnaky>1:  
        
        # Print omega to an h5 file
        with h5py.File(potential_path, 'w') as h5_file:  
            h5_file.create_dataset("vec_time", data=vec_time) 
            h5_file.create_dataset("phi_vs_tkxky", data=phi_vs_tkxky) 
            h5_file.create_dataset("phi2_vs_tkxky", data=phi2_vs_tkxky)  
            
    # txt file for 1 mode per simulation
    elif nakxnaky==1:
        
        # Data to be written
        header = " "*7+"time"+" "*13+"|phi|^2"+" "*10+"Re(phi)"+" "*12+"Im(phi)"+"\n"
        potential_data = np.concatenate((vec_time[:,np.newaxis], phi2_vs_tkxky[:,0,0,np.newaxis], phi_vs_tkxky[:,0,0,np.newaxis].real, phi_vs_tkxky[:,0,0,np.newaxis].imag), axis=1)
             
        # Write the potential data to a text file   
        potential_file = open(potential_path,'w')  
        potential_file.write(header)
        np.savetxt(potential_file, potential_data, fmt='%16.8e', delimiter="   ")  
        potential_file.close() 
        
    # Notify that we finished creating the file
    print(status+"   ---> The phi(t) file is saved as " +  potential_path.parent.name+"/"+potential_path.name)  
    return  

#---------------------------------------------
def read_potentialVsTimeLinearSimulations(dt, input_file, netcdf_path):

    # Read the geometry data in the output file
    vmec_filename = read_vmecFileNameFromInputFile(input_file) 
    path = create_dummyPathObject(input_file, vmec_filename)
    geometry = read_outputFileForGeometry(path) 
    dl_over_B = geometry["dl_over_B"]  
    
    # Read the potential from the output file
    netcdf_data = read_outputFile(netcdf_path)  
    vec_time = read_netcdfVariables('vec_time', netcdf_data) 
    time_indices = get_indicesAtFixedStep(vec_time, dt) 
    phi_vs_tzkxkyri = read_netcdfVariables('phi_vs_tzkxkyri', netcdf_data, time_indices)  
    vec_time = vec_time[time_indices]
    netcdf_data.close()
    
    # Construct the potential from its real and imaginary part
    phi_vs_tzkxky = phi_vs_tzkxkyri[:,:,:,:,0] + 1j*phi_vs_tzkxkyri[:,:,:,:,1]  
    
    # Sum away the (kx,ky) modes and the (z) dimension   
    phi_vs_tkxky = np.sum(phi_vs_tzkxky[:,:,:,:]*dl_over_B[np.newaxis,:,np.newaxis,np.newaxis], axis=1) 
    
    # Carefull: for phi2 we need to first take abs(phi)**2 and then sum away the dimensions (kx,ky,z)!
    phi2_vs_tzkxky = np.abs(phi_vs_tzkxky)**2
    phi2_vs_tkxky = np.sum(phi2_vs_tzkxky[:,:,:,:]*dl_over_B[np.newaxis,:,np.newaxis,np.newaxis], axis=1)
 
    # Return the data
    return vec_time, phi_vs_tkxky, phi2_vs_tkxky
 
