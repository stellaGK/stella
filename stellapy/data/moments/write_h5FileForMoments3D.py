
import copy
import os, h5py
import numpy as np
import time as timer  
from scipy.io import netcdf as scnetcdf  
from datetime import datetime, timedelta 
from stellapy.utils.decorators.verbose import noverbose 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder     
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry  
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables 
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile, read_modeFromInputFile  
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.stella.load_stellaVariables import stella_variables

################################################################################
#                       WRITE THE DIMENSIONS TO AN H5 FILE
################################################################################

@noverbose
def write_h5FileForMoments3D(folder, dt=10, automatic=False):     
     
    # Time step
    dt = int(dt) if int(dt)==dt else dt     
     
    # Get the input files for which a netcdf file exists
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 
 
    # Iterate through the input files  
    for input_file in input_files:  
            
        # Start timer
        start = timer.time()
    
        # Processing status
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
        nonlinear = read_linearNonlinearFromInputFile(input_file)[1] 
        if (not nonlinear) and automatic: continue
         
        # Path of the new file
        netcdf_path = input_file.with_suffix(".out.nc")
        moments_path = input_file.with_suffix(".dt"+str(dt)+".moments3D")   
                
        # Check whether the potential file is older than the simulation
        outputFileHasChanged = False
        if os.path.isfile(moments_path):  
            if datetime.fromtimestamp(netcdf_path.stat().st_mtime) > (datetime.fromtimestamp(moments_path.stat().st_mtime)+timedelta(minutes=5)):
                outputFileHasChanged = True
            if nonlinear:
                with h5py.File(moments_path, 'r') as h5_file:   
                    if "upar2_vs_tsz" not in h5_file.keys():
                        print("Remove the 3D moments file since <upar2_vs_tsz> was missing.")
                        os.system("rm "+str(moments_path)) 
                
        # Notify that the file already existed   
        if os.path.isfile(moments_path) and not outputFileHasChanged:
            print(status+"The 3D moments file already exists:", moments_path.parent.name+"/"+moments_path.name)
            continue
             
        # If the output file changed, then append to the h5 file
        elif os.path.isfile(moments_path) and outputFileHasChanged:
            
            #  Read the variables and the time axis
            netcdf_data = read_outputFile(netcdf_path)  
            vec_time = read_netcdfVariables('vec_time', netcdf_data)  
            variables = list(netcdf_data.variables.keys())
            netcdf_data.close()    
            
            # Check whether the moments are written
            if "density" not in variables:
                print(status+"The moments were not written to the netcdf file.")
                continue 
            
            # Edit the time vector in the output file
            indices  = get_indicesAtFixedStep(vec_time, dt)
            vec_time = vec_time[indices]
            vec_time = [round(n, 8) for n in vec_time]
            tlast_outputfile = vec_time[-1]  
    
            # Read the data in the h5 file            
            with h5py.File(moments_path, 'r') as f:  
                vec_time_h5file = f['vec_time'][()]
                upar_vs_tsz_h5file = f['upar_vs_tsz'][()]
                dens_vs_tsz_h5file = f['dens_vs_tsz'][()]
                temp_vs_tsz_h5file = f['temp_vs_tsz'][()]
                tlast_h5file = vec_time_h5file[-1] 
                if nonlinear:
                    upar_vs_tskx_h5file = f['upar_vs_tskx'][()]
                    dens_vs_tskx_h5file = f['dens_vs_tskx'][()]
                    temp_vs_tskx_h5file = f['temp_vs_tskx'][()]
                    upar_vs_tsky_h5file = f['upar_vs_tsky'][()]
                    dens_vs_tsky_h5file = f['dens_vs_tsky'][()]
                    temp_vs_tsky_h5file = f['temp_vs_tsky'][()]
                    upar2_vs_tsz_h5file = f['upar2_vs_tsz'][()]
                    dens2_vs_tsz_h5file = f['dens2_vs_tsz'][()]
                    temp2_vs_tsz_h5file = f['temp2_vs_tsz'][()]
            
            # If the output file (*.out.nc) contains extra time steps, append to the h5 file 
            if tlast_outputfile > tlast_h5file: 
            
                # Get the integration weigths and the time indices
                dl_over_B = get_integrationWeightsAlongZ(input_file) 
    
                # Get the data at every <dt> timestep 
                index_tlast = np.argwhere(tlast_h5file < vec_time)[0][0] 
                time_indices = get_indicesAtFixedStep(vec_time, dt)
                vec_time = np.array(vec_time)[time_indices]  
                
                # Read the potential versus time from the output file ['upar_vs_tszkxkyri', 'dens_vs_tszkxkyri', 'temp_vs_tszkxkyri']
                upar_vs_tsz, upar2_vs_tsz, upar_vs_tsky, upar_vs_tskx = read_moments3D(netcdf_path, 'upar_vs_tszkxkyri', time_indices, dl_over_B, input_file)
                dens_vs_tsz, dens2_vs_tsz, dens_vs_tsky, dens_vs_tskx = read_moments3D(netcdf_path, 'dens_vs_tszkxkyri', time_indices, dl_over_B, input_file)
                temp_vs_tsz, temp2_vs_tsz, temp_vs_tsky, temp_vs_tskx = read_moments3D(netcdf_path, 'temp_vs_tszkxkyri', time_indices, dl_over_B, input_file)         
        
                # Combine the old and new arrays
                vec_time = np.append(vec_time_h5file, vec_time[index_tlast:], axis=0)  
                upar_vs_tsz = np.append(upar_vs_tsz_h5file, upar_vs_tsz[index_tlast:,:,:], axis=0)  
                dens_vs_tsz = np.append(dens_vs_tsz_h5file, dens_vs_tsz[index_tlast:,:,:], axis=0)  
                temp_vs_tsz = np.append(temp_vs_tsz_h5file, temp_vs_tsz[index_tlast:,:,:], axis=0)  
                upar_vs_tskx = np.append(upar_vs_tskx_h5file, upar_vs_tskx[index_tlast:,:,:], axis=0)  
                dens_vs_tskx = np.append(dens_vs_tskx_h5file, dens_vs_tskx[index_tlast:,:,:], axis=0)  
                temp_vs_tskx = np.append(temp_vs_tskx_h5file, temp_vs_tskx[index_tlast:,:,:], axis=0)  
                upar_vs_tsky = np.append(upar_vs_tsky_h5file, upar_vs_tsky[index_tlast:,:,:], axis=0)  
                dens_vs_tsky = np.append(dens_vs_tsky_h5file, dens_vs_tsky[index_tlast:,:,:], axis=0)  
                temp_vs_tsky = np.append(temp_vs_tsky_h5file, temp_vs_tsky[index_tlast:,:,:], axis=0)  
                upar2_vs_tsz = np.append(upar2_vs_tsz_h5file, upar2_vs_tsz[index_tlast:,:,:], axis=0)  
                dens2_vs_tsz = np.append(dens2_vs_tsz_h5file, dens2_vs_tsz[index_tlast:,:,:], axis=0)  
                temp2_vs_tsz = np.append(temp2_vs_tsz_h5file, temp2_vs_tsz[index_tlast:,:,:], axis=0)
                
            # If the output file has the same time steps, touch the text file
            elif tlast_outputfile <= tlast_h5file:
                print(status+"The 3D moments file is up to date:", moments_path.parent.name+"/"+moments_path.name)  
                os.system("touch "+str(moments_path))
                continue
            
        # Otherwise read the 3D moments data
        elif not os.path.isfile(moments_path): 
            
            # Read the time axis and variables in the netcdf file
            netcdf_data = read_outputFile(netcdf_path) 
            vec_time = read_netcdfVariables('vec_time', netcdf_data)   
            variables = list(netcdf_data.variables.keys())
            netcdf_data.close()   
            
            # Check whether the moments are written
            if "density" not in variables:
                print(status+"The moments were not written to the netcdf file.")
                continue    
            
            # Get the integration weigths and the time indices
            dl_over_B = get_integrationWeightsAlongZ(input_file) 

            # Get the data at every <dt> timestep 
            time_indices = get_indicesAtFixedStep(vec_time, dt)
            vec_time = vec_time[time_indices]  
            
            # Read the potential versus time from the output file ['upar_vs_tszkxkyri', 'dens_vs_tszkxkyri', 'temp_vs_tszkxkyri']
            upar_vs_tsz, upar2_vs_tsz, upar_vs_tsky, upar_vs_tskx = read_moments3D(netcdf_path, 'upar_vs_tszkxkyri', time_indices, dl_over_B, input_file)
            dens_vs_tsz, dens2_vs_tsz, dens_vs_tsky, dens_vs_tskx = read_moments3D(netcdf_path, 'dens_vs_tszkxkyri', time_indices, dl_over_B, input_file)
            temp_vs_tsz, temp2_vs_tsz, temp_vs_tsky, temp_vs_tskx = read_moments3D(netcdf_path, 'temp_vs_tszkxkyri', time_indices, dl_over_B, input_file)         
        
        # Create the new h5 file 
        with h5py.File(moments_path, 'w') as h5_file:   
            h5_file.create_dataset("vec_time", data=vec_time) 
            h5_file.create_dataset("upar_vs_tsz", data=upar_vs_tsz) 
            h5_file.create_dataset("dens_vs_tsz", data=dens_vs_tsz)  
            h5_file.create_dataset("temp_vs_tsz", data=temp_vs_tsz)  
            if nonlinear:
                h5_file.create_dataset("upar2_vs_tsz", data=upar2_vs_tsz) 
                h5_file.create_dataset("dens2_vs_tsz", data=dens2_vs_tsz)  
                h5_file.create_dataset("temp2_vs_tsz", data=temp2_vs_tsz)  
                h5_file.create_dataset("upar_vs_tskx", data=upar_vs_tskx) 
                h5_file.create_dataset("dens_vs_tskx", data=dens_vs_tskx)  
                h5_file.create_dataset("temp_vs_tskx", data=temp_vs_tskx)  
                h5_file.create_dataset("upar_vs_tsky", data=upar_vs_tsky) 
                h5_file.create_dataset("dens_vs_tsky", data=dens_vs_tsky)  
                h5_file.create_dataset("temp_vs_tsky", data=temp_vs_tsky)  
                
        # Notify that we finished creating the file
        end = timer.time(); elapsed_time = str(timedelta(seconds=round(end-start,0))) 
        if outputFileHasChanged: print(status+"   ---> The 3D moments file is updated as " +  moments_path.parent.name+"/"+moments_path.name+" ("+elapsed_time+")")   
        if not outputFileHasChanged:  print(status+"   ---> The 3D moments file is saved as " +  moments_path.parent.name+"/"+moments_path.name+" ("+elapsed_time+")")   
 
    return 

#---------------------------------------------
def get_integrationWeightsAlongZ(input_file):
    vmec_filename = read_vmecFileNameFromInputFile(input_file)
    nonlinear = read_linearNonlinearFromInputFile(input_file)[1]
    path = create_dummyPathObject(input_file, vmec_filename, nonlinear)
    geometry = read_outputFileForGeometry(path)  
    dl_over_B = geometry["dl_over_B"]   
    return dl_over_B

#---------------------------------------------
def read_moments3D(netcdf_path, variable, time_indices, dl_over_B, input_file):    
     
    # Get the dimensions and the stella key of the variable 
    key = stella_variables[variable][0]
     
    # Open the netcdf file
    netcdf_file = scnetcdf.netcdf_file(netcdf_path,'r') 
    variable_vs_tszkxkyri = copy.deepcopy(netcdf_file.variables[key][time_indices,:,0,:,:,:,:])
    netcdf_file.close()  
 
    # Sum away the (kx,ky) dimensions
    variable_vs_tszri = np.sum(variable_vs_tszkxkyri[:,:,:,:,0,:],axis=3) + 2*np.sum(variable_vs_tszkxkyri[:,:,:,:,1:,:],axis=(3,4))  
 
    # Sum away the (z) dimension
    variable_vs_tskxkyri = np.sum(variable_vs_tszkxkyri[:,:,:,:,:,:]*dl_over_B[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis,np.newaxis], axis=2) 

    # Sum away the (kx) dimension
    variable_vs_tskyri = np.sum(variable_vs_tskxkyri[:,:,:,:,:],axis=2)  
     
    # Sum away the (ky) dimension
    variable_vs_tskxri = variable_vs_tskxkyri[:,:,:,0,:] + 2*np.sum(variable_vs_tskxkyri[:,:,:,1:,:],axis=3)   

    # Combine the real and imaginary part to a complex number
    variable_vs_tsz = variable_vs_tszri[:,:,:,0] + 1j*variable_vs_tszri[:,:,:,1]
    variable_vs_tsky = variable_vs_tskyri[:,:,:,0] + 1j*variable_vs_tskyri[:,:,:,1]
    variable_vs_tskx = variable_vs_tskxri[:,:,:,0] + 1j*variable_vs_tskxri[:,:,:,1]
    
    # Clean up the memory 
    del variable_vs_tskxkyri
    del variable_vs_tskyri
    del variable_vs_tskxri
    del variable_vs_tszri
    
    # Calculate variable^2(z)
    variable_vs_tszkxky = variable_vs_tszkxkyri[:,:,:,:,:,0] + 1j*variable_vs_tszkxkyri[:,:,:,:,:,1]
    variable2_vs_tszkxky = np.abs(variable_vs_tszkxky)**2
    variable2_vs_tsz = np.sum(variable2_vs_tszkxky[:,:,:,0,:],axis=3) + 2*np.sum(variable2_vs_tszkxky[:,:,:,1:,:],axis=(3,4))

    # Clean up the memory 
    del variable_vs_tszkxkyri
    del variable_vs_tszkxky
    del variable2_vs_tszkxky
    
    # For a linear mode, we need to sum the non-zonal modes twice
    if read_linearNonlinearFromInputFile(input_file)[0]: 
        if read_modeFromInputFile(input_file)[1]!=0: 
            variable_vs_tsz = variable_vs_tsz*2 
            variable2_vs_tsz = variable2_vs_tsz*2 
    
    # Return the data 
    return variable_vs_tsz, variable2_vs_tsz, variable_vs_tsky, variable_vs_tskx

    