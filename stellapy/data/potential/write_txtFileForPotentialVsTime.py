

import os
import numpy as np 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder    
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry 
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile
from stellapy.data.input.read_inputFile import read_numberOfModesFromInputFile 
from stellapy.data.input.read_inputFile import read_modeFromInputFile
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables

################################################################################
#                        TIME EVOLUTION OF THE POTENTIAL
################################################################################ 

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
        if not nonlinear: write_txtFileForPotentialVsTimeLinearSimulations(input_file, dt, status)
        if nonlinear: write_txtFileForPotentialVsTimeNonlinearSimulations(input_file, dt, status)
        
        
#===============================================================================
#                            NONLINEAR SIMULATIONS                             #
#===============================================================================
 
def write_txtFileForPotentialVsTimeNonlinearSimulations(input_file, dt, status):   
        
    try:
    
        # Path of the new file  
        txt_file = input_file.with_suffix(".dt"+str(dt)+".phi_vs_t") 
        netcdf_file = input_file.with_suffix(".out.nc") 
        
        # Check whether the txt file is older than the simulation
        outputFileHasChanged = False
        if os.path.isfile(txt_file):  
            if netcdf_file.stat().st_mtime > txt_file.stat().st_mtime:
                outputFileHasChanged = True 
        
        # Notify that the file already existed   
        if os.path.isfile(txt_file) and not outputFileHasChanged:
            print(status+"The phi(t) file already exists:", txt_file.parent.name+"/"+txt_file.name)
            return
        
        # If the output file changed, then append to the txt file
        elif os.path.isfile(txt_file) and outputFileHasChanged:
            
            # Check whether the output file (*.out.nc) contains extra time points
            netcdf_data  = read_outputFile(netcdf_file)   
            vec_time     = read_netcdfVariables('vec_time', netcdf_data, netcdf_path=netcdf_file) 
            netcdf_data.close()  
            
            # Edit the time vector in the output file
            indices  = get_indicesAtFixedStep(vec_time, dt)
            vec_time = vec_time[indices]
            vec_time = [round(n, 8) for n in vec_time]
            tlast_outputfile = vec_time[-1]  
    
            # Read the potential versus time from the txt file
            data = np.loadtxt(txt_file,skiprows=1,dtype='float').reshape(-1, 6)
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
                phi2_vs_t_nozonal_outputFile, vec_time_outputFile = read_potentialVsTimeNonlinearSimulations(dt, input_file, netcdf_file)  
                vec_time = np.append(vec_time_txtfile, vec_time_outputFile[index_tlast:], axis=0)  
                phi_vs_t = np.append(phi_vs_t_txtfile, phi_vs_t_outputFile[index_tlast:], axis=0)  
                phi2_vs_t = np.append(phi2_vs_t_txtfile, phi2_vs_t_outputFile[index_tlast:], axis=0) 
                phi2_vs_t_zonal = np.append(phi2_vs_t_zonal_txtfile, phi2_vs_t_zonal_outputFile[index_tlast:], axis=0)
                phi2_vs_t_nozonal = np.append(phi2_vs_t_nozonal_txtfile, phi2_vs_t_nozonal_outputFile[index_tlast:], axis=0) 
                
            # If the output file has the same time steps, touch the text file
            elif tlast_outputfile <= tlast_txtfile:
                print(status+"The phi(t) file is up to date:", txt_file.parent.name+"/"+txt_file.name)  
                os.system("touch "+str(txt_file))
                return
            
        # Create the text file for dphiz(t)
        elif not os.path.isfile(txt_file):    
            
            # Read the potential versus time from the output file
            phi_vs_t, phi2_vs_t, phi2_vs_t_zonal, phi2_vs_t_nozonal, vec_time = read_potentialVsTimeNonlinearSimulations(dt, input_file, netcdf_file) 
                 
        # Write the potential data to a text file  
        dphiz_file = open(txt_file,'w') 
        dphiz_file.write(" "*7+"time"+" "*13+"|phi|^2"+" "*10+"zonal|phi|^2"+" "*6+"nozonal|phi|^2"+" "*8+"Re(phi)"+" "*12+"Im(phi)"+"\n")
        for i in range(len(phi_vs_t)): 
            row = [[vec_time[i], phi2_vs_t[i], phi2_vs_t_zonal[i], phi2_vs_t_nozonal[i], phi_vs_t[i].real, phi_vs_t[i].imag]]   
            np.savetxt(dphiz_file, row, fmt='%16.8e', delimiter="   ") 
        dphiz_file.close() 
            
        # Notify that we finished creating the file
        if outputFileHasChanged: print(status+"   ---> The phi(t) file is updated as " +  txt_file.parent.name+"/"+txt_file.name)   
        if not outputFileHasChanged:  print(status+"   ---> The phi(t) file is saved as " +  txt_file.parent.name+"/"+txt_file.name)  
    
    except:
        print(status+"   ---> Something went wrong writing phi(t) for " +  txt_file.parent.name+"/"+txt_file.name)   
    return  

#---------------------------------------------
def read_potentialVsTimeNonlinearSimulations(dt, input_file, netcdf_file):

    # Read the geometry data in the output file
    vmec_filename = read_vmecFileNameFromInputFile(input_file)
    nonlinear = read_linearNonlinearFromInputFile(input_file)[1]
    path = create_dummyPathObject(input_file, vmec_filename, nonlinear)
    geometry = read_outputFileForGeometry(path) 
    dl_over_B = geometry["dl_over_B"]  
    
    # Read the potential from the output file
    netcdf_data = read_outputFile(netcdf_file)  
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
    
    # For a linear mode, we need to sum the non-zonal modes twice
    if read_linearNonlinearFromInputFile(input_file)[0]:  
        if read_modeFromInputFile(input_file)[1]!=0:
            phi_vs_t = phi_vs_t*2
            phi2_vs_t = phi2_vs_t*2
            
    # Return the data
    return phi_vs_t, phi2_vs_t, phi2_vs_t_zonal, phi2_vs_t_nozonal, vec_time

#===============================================================================
#                             LINEAR SIMULATIONS                               #
#===============================================================================

def write_txtFileForPotentialVsTimeLinearSimulations(input_file, dt, status):   
        
#     try:
    
        # Path of the new file  
        txt_file = input_file.with_suffix(".dt"+str(dt)+".phi_vs_t") 
        netcdf_file = input_file.with_suffix(".out.nc") 
        
        # Check the (kx,ky) dimensions
        dim_kx, dim_ky = read_numberOfModesFromInputFile(input_file)
        multiple_modes_per_file = True if (dim_kx*dim_ky>1) else False
        
        # Check whether the txt file is older than the simulation
        outputFileHasChanged = False
        if os.path.isfile(txt_file):  
            if netcdf_file.stat().st_mtime > txt_file.stat().st_mtime:
                outputFileHasChanged = True 
        
        # Notify that the file already existed   
        if os.path.isfile(txt_file) and not outputFileHasChanged:
            print(status+"The phi(t) file already exists:", txt_file.parent.name+"/"+txt_file.name)
            return
            
        # Create the text file for dphiz(t)
        else:    
            
            # Read the potential versus time from the output file
            potential_data = read_potentialVsTimeLinearSimulations(dt, input_file, netcdf_file) 
            
            # Change the data according to whether we have one or multiple (kx,ky)
            if not multiple_modes_per_file: 
                header = " "*7+"time"+" "*13+"|phi|^2"+" "*10+"Re(phi)"+" "*12+"Im(phi)"+"\n"
                potential_data = potential_data[0,0,:,2:]
            if multiple_modes_per_file: 
                potential_data = potential_data.reshape(-1, 6)
                header = " "*8+"kx"+" "*17+"ky"+" "*16+"time"+" "*13+"|phi|^2"+" "*12+"Re(phi)"+" "*12+"Im(phi)"+"\n"
                 
        # Write the potential data to a text file  
        dphiz_file = open(txt_file,'w')  
        dphiz_file.write(header)
        np.savetxt(dphiz_file, potential_data, fmt='%16.8e', delimiter="   ")  
        dphiz_file.close() 
            
        # Notify that we finished creating the file
        print(status+"   ---> The phi(t) file is saved as " +  txt_file.parent.name+"/"+txt_file.name)  
    
#     except:
#         print(status+"   ---> Something went wrong writing phi(t) for " +  txt_file.parent.name+"/"+txt_file.name)   
#     return  

#---------------------------------------------
def read_potentialVsTimeLinearSimulations(dt, input_file, netcdf_file):

    # Read the geometry data in the output file
    vmec_filename = read_vmecFileNameFromInputFile(input_file)
    nonlinear = read_linearNonlinearFromInputFile(input_file)[1]
    dim_kx, dim_ky = read_numberOfModesFromInputFile(input_file)
    path = create_dummyPathObject(input_file, vmec_filename, nonlinear)
    geometry = read_outputFileForGeometry(path) 
    dl_over_B = geometry["dl_over_B"]  
    
    # Read the potential from the output file
    netcdf_data = read_outputFile(netcdf_file)  
    vec_time = read_netcdfVariables('vec_time', netcdf_data) 
    phi_vs_tzkxkyri = read_netcdfVariables('phi_vs_tzkxkyri', netcdf_data) 
    vec_kx = read_netcdfVariables('vec_kx', netcdf_data) 
    vec_ky = read_netcdfVariables('vec_ky', netcdf_data) 
    netcdf_data.close()
    
    # Get the data at every <dt> timestep 
    indices = get_indicesAtFixedStep(vec_time, dt) 
    phi_vs_tzkxky = phi_vs_tzkxkyri[indices,:,:,:,0] + 1j*phi_vs_tzkxkyri[indices,:,:,:,1] 
    vec_time = vec_time[indices]
    
    # Sum away the (kx,ky) modes and the (z) dimension   
    phi_vs_tkxky = np.sum(phi_vs_tzkxky[:,:,:,:]*dl_over_B[np.newaxis,:,np.newaxis,np.newaxis], axis=1) 
    
    # Carefull: for phi2 we need to first take abs(phi)**2 and then sum away the dimensions (kx,ky,z)!
    phi2_vs_tzkxky = np.abs(phi_vs_tzkxky)**2
    phi2_vs_tkxky = np.sum(phi2_vs_tzkxky[:,:,:,:]*dl_over_B[np.newaxis,:,np.newaxis,np.newaxis], axis=1)
    
    # Swap the axis so that we have (kx,ky,t) instead of (t,kx,ky)
    phi_vs_kxtky = np.swapaxes(phi_vs_tkxky,0,1)
    phi_vs_kxkyt = np.swapaxes(phi_vs_kxtky,1,2)
    phi2_vs_kxtky = np.swapaxes(phi2_vs_tkxky,0,1)
    phi2_vs_kxkyt = np.swapaxes(phi2_vs_kxtky,1,2)  
    
    # Create an array to hold the (kx,ky) information versus (kx,ky,t)
    kx_ky_time_array = np.zeros((dim_kx, dim_ky, len(vec_time), 3))
    for i in range(dim_kx):
        for j in range(dim_ky):
            kx_ky_time_array[i,j,:,0] = vec_kx[i]
            kx_ky_time_array[i,j,:,1] = vec_ky[j]
            kx_ky_time_array[i,j,:,2] = vec_time
            
    # Make one big data array holding (kx,ky,t,phi2,phireal,phiimag) in the fourth column as a function of (kx,ky,t)
    potential_data = np.concatenate((kx_ky_time_array[:,:,:,np.newaxis,0], kx_ky_time_array[:,:,:,np.newaxis,1], kx_ky_time_array[:,:,:,np.newaxis,2], phi2_vs_kxkyt[:,:,:,np.newaxis], phi_vs_kxkyt[:,:,:,np.newaxis].real, phi_vs_kxkyt[:,:,:,np.newaxis].imag), axis=3)
    
    # Return the data
    return potential_data
            
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    import pathlib
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI/test")
    write_txtFileForPotentialVsTime(folder)
    
