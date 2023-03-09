
import os, h5py
import numpy as np
import time as timer  
from datetime import datetime, timedelta 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder     
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry  
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables 
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile, read_modeFromInputFile 
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile

################################################################################
#                       WRITE THE DIMENSIONS TO AN H5 FILE
################################################################################

def write_h5FileForMoments2D(folder, dt=1):     
     
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
        linear, _ = read_linearNonlinearFromInputFile(input_file)  
        if linear: continue
         
        # Path of the new file
        netcdf_path = input_file.with_suffix(".out.nc")
        moments_path = input_file.with_suffix(".dt"+str(dt)+".moments2D")   
                
        # Check whether the potential file is older than the simulation
        outputFileHasChanged = False
        if os.path.isfile(moments_path):  
            if datetime.fromtimestamp(netcdf_path.stat().st_mtime) > (datetime.fromtimestamp(moments_path.stat().st_mtime)+timedelta(minutes=5)):
                outputFileHasChanged = True
                
        # Notify that the file already existed   
        if os.path.isfile(moments_path) and not outputFileHasChanged:
            print(status+"The moments 2D file already exists:", moments_path.parent.name+"/"+moments_path.name)
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
            
            # Read the moments versus time from the txt file      
            with h5py.File(moments_path, 'r') as f:  
                vec_time_h5file = f['vec_time'][()]
                upar_vs_ts_h5file = f['upar_vs_ts'][()]
                dens_vs_ts_h5file = f['dens_vs_ts'][()]
                temp_vs_ts_h5file = f['temp_vs_ts'][()]
                upar2_vs_ts_h5file = f['upar2_vs_ts'][()]
                dens2_vs_ts_h5file = f['dens2_vs_ts'][()]
                temp2_vs_ts_h5file = f['temp2_vs_ts'][()]
                tlast_h5file = vec_time_h5file[-1] 
            
            # If the output file (*.out.nc) contains extra time steps, append to the h5 file 
            if tlast_outputfile > tlast_h5file: 
            
                # Get the integration weigths and the time indices
                dl_over_B = get_integrationWeightsAlongZ(input_file) 
    
                # Get the data at every <dt> timestep 
                index_tlast = np.argwhere(tlast_h5file < vec_time)[0][0] 
                time_indices = get_indicesAtFixedStep(vec_time, dt)
                vec_time = np.array(vec_time)[time_indices]  
                
                # Read the potential versus time from the output file ['upar_vs_tszkxkyri', 'dens_vs_tszkxkyri', 'temp_vs_tszkxkyri']
                upar_vs_ts, upar2_vs_ts = read_moments2D(netcdf_path, 'upar_vs_tszkxkyri', time_indices, dl_over_B, input_file)
                dens_vs_ts, dens2_vs_ts = read_moments2D(netcdf_path, 'dens_vs_tszkxkyri', time_indices, dl_over_B, input_file)
                temp_vs_ts, temp2_vs_ts = read_moments2D(netcdf_path, 'temp_vs_tszkxkyri', time_indices, dl_over_B, input_file)         
        
                # Combine the old and new arrays
                vec_time = np.append(vec_time_h5file, vec_time[index_tlast:], axis=0)  
                upar_vs_ts = np.append(upar_vs_ts_h5file, upar_vs_ts[index_tlast:,:], axis=0)  
                dens_vs_ts = np.append(dens_vs_ts_h5file, dens_vs_ts[index_tlast:,:], axis=0)  
                temp_vs_ts = np.append(temp_vs_ts_h5file, temp_vs_ts[index_tlast:,:], axis=0)   
                upar2_vs_ts = np.append(upar2_vs_ts_h5file, upar2_vs_ts[index_tlast:,:], axis=0)  
                dens2_vs_ts = np.append(dens2_vs_ts_h5file, dens2_vs_ts[index_tlast:,:], axis=0)  
                temp2_vs_ts = np.append(temp2_vs_ts_h5file, temp2_vs_ts[index_tlast:,:], axis=0)
                
            # If the output file has the same time steps, touch the text file
            elif tlast_outputfile <= tlast_h5file:
                print(status+"The moments 2D file is up to date:", moments_path.parent.name+"/"+moments_path.name)  
                os.system("touch "+str(moments_path))
                continue
            
        # Otherwise read the 2D moments data from the netcdf file and create a new txt file
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
            time_indices = get_indicesAtFixedStep(vec_time, dt)
            vec_time = vec_time[time_indices]  
            
            # Read the potential versus time from the output file ['upar_vs_tszkxkyri', 'dens_vs_tszkxkyri', 'temp_vs_tszkxkyri']
            upar_vs_ts, upar2_vs_ts = read_moments2D(netcdf_path, 'upar_vs_tszkxkyri', time_indices, dl_over_B, input_file)
            dens_vs_ts, dens2_vs_ts = read_moments2D(netcdf_path, 'dens_vs_tszkxkyri', time_indices, dl_over_B, input_file)
            temp_vs_ts, temp2_vs_ts = read_moments2D(netcdf_path, 'temp_vs_tszkxkyri', time_indices, dl_over_B, input_file)    
            
        # Create the new h5 file 
        with h5py.File(moments_path, 'w') as h5_file:   
            h5_file.create_dataset("vec_time", data=vec_time) 
            h5_file.create_dataset("upar_vs_ts", data=upar_vs_ts) 
            h5_file.create_dataset("dens_vs_ts", data=dens_vs_ts)  
            h5_file.create_dataset("temp_vs_ts", data=temp_vs_ts)  
            h5_file.create_dataset("upar2_vs_ts", data=upar2_vs_ts) 
            h5_file.create_dataset("dens2_vs_ts", data=dens2_vs_ts)  
            h5_file.create_dataset("temp2_vs_ts", data=temp2_vs_ts)  
                
        # Notify that we finished creating the file
        end = timer.time(); elapsed_time = str(timedelta(seconds=round(end-start,0))) 
        if outputFileHasChanged: print(status+"   ---> The moments 2D file is updated as " +  moments_path.parent.name+"/"+moments_path.name+" ("+elapsed_time+")")   
        if not outputFileHasChanged:  print(status+"   ---> The moments 2D file is saved as " +  moments_path.parent.name+"/"+moments_path.name+" ("+elapsed_time+")")   
 
    return 

#---------------------------------------------
def get_integrationWeightsAlongZ(input_file):
    vmec_filename = read_vmecFileNameFromInputFile(input_file) 
    path = create_dummyPathObject(input_file, vmec_filename)
    geometry = read_outputFileForGeometry(path)  
    dl_over_B = geometry["dl_over_B"]   
    return dl_over_B

#---------------------------------------------
def read_moments2D(netcdf_path, variable, time_indices, dl_over_B, input_file):     
     
    # Open the netcdf file
    netcdf_data = read_outputFile(netcdf_path)  
    variable_vs_tszkxkyri = read_netcdfVariables(variable, netcdf_data, time_indices)   
    netcdf_data.close()    
 
    # Sum away the (kx,ky, z) dimensions
    variable_vs_tszri = np.sum(variable_vs_tszkxkyri[:,:,:,:,0,:],axis=3) + 2*np.sum(variable_vs_tszkxkyri[:,:,:,:,1:,:],axis=(3,4))  
    variable_vs_tsri = np.sum(variable_vs_tszri*dl_over_B[np.newaxis,np.newaxis,:,np.newaxis], axis=2); del variable_vs_tszri

    # Combine the real and imaginary part to a complex number
    variable_vs_ts = variable_vs_tsri[:,:,0] + 1j*variable_vs_tsri[:,:,1]; del variable_vs_tsri
    
    # Calculate variable^2(z)
    variable_vs_tszkxky = variable_vs_tszkxkyri[:,:,:,:,:,0] + 1j*variable_vs_tszkxkyri[:,:,:,:,:,1]; del variable_vs_tszkxkyri
    variable2_vs_tszkxky = np.abs(variable_vs_tszkxky)**2; del variable_vs_tszkxky
    variable2_vs_tsz = np.sum(variable2_vs_tszkxky[:,:,:,0,:],axis=3) + 2*np.sum(variable2_vs_tszkxky[:,:,:,1:,:],axis=(3,4)); del variable2_vs_tszkxky
    variable2_vs_ts = np.sum(variable2_vs_tsz*dl_over_B[np.newaxis,np.newaxis,:], axis=2); del variable2_vs_tsz
    
    # For a linear mode, we need to sum the non-zonal modes twice
    if read_linearNonlinearFromInputFile(input_file)[0]: 
        if read_modeFromInputFile(input_file)[1]!=0: 
            variable_vs_ts = variable_vs_ts*2 
            variable2_vs_ts = variable2_vs_ts*2 
    
    # Return the data 
    return variable_vs_ts, np.real(variable2_vs_ts)


  

    
