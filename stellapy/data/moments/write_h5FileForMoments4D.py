  
import os, h5py
import numpy as np
from datetime import datetime, timedelta
from stellapy.utils.decorators.verbose import noverbose 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder     
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry  
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables 
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile  
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile

################################################################################
#                       WRITE THE MOMENTS TO AN H5 FILE
################################################################################

@noverbose
def write_h5FileForMoments4D(folder, dt=10):    
     
    # Time step
    dt = int(dt) if int(dt)==dt else dt     
     
    # Get the input files for which a netcdf file exists
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 
 
    # Iterate through the input files  
    for input_file in input_files:  
            
        # Only write the moments file for nonlinear simulations 
        if read_linearNonlinearFromInputFile(input_file)[1]:
        
            # Processing status
            status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
            
            # Path of the new file
            netcdf_file = input_file.with_suffix(".out.nc")
            moments_file = input_file.with_suffix(".dt"+str(dt)+".moments4D")   
                    
            # Check whether the potential file is older than the simulation
            outputFileHasChanged = False
            if os.path.isfile(moments_file):  
                if datetime.fromtimestamp(netcdf_file.stat().st_mtime) > (datetime.fromtimestamp(moments_file.stat().st_mtime)+timedelta(minutes=5)):
                    outputFileHasChanged = True 
                with h5py.File(moments_file, 'r') as h5_file:   
                    if "dens_vs_tskxky_zeta0" not in h5_file.keys():
                        print("Remove the 4D moments file since <dens_vs_tskxky_zeta0> was missing.")
                        os.system("rm "+str(moments_file))
                    elif "dens2_vs_tskxky" not in h5_file.keys():
                        print("Remove the 4D moments file since <dens2_vs_tskxky> was missing.")
                        os.system("rm "+str(moments_file))
                    
            # Notify that the file already existed   
            if os.path.isfile(moments_file) and not outputFileHasChanged:
                print(status+"The 4D moments file already exists:", moments_file.parent.name+"/"+moments_file.name)
                continue
                 
            # If the output file changed, then append to the h5 file
            elif os.path.isfile(moments_file) and outputFileHasChanged:
                
                # Check whether the moments are written
                netcdf_data  = read_outputFile(netcdf_file)   
                if "density" not in netcdf_file.variables.keys():
                    print(status+"The moments were not written to the netcdf file.")
                
                # Check whether the output file (*.out.nc) contains extra time points 
                vec_time     = read_netcdfVariables('vec_time', netcdf_data) 
                netcdf_data.close()  
                
                # Edit the time vector in the output file
                indices  = get_indicesAtFixedStep(vec_time, dt)
                vec_time = vec_time[indices]
                vec_time = [round(n, 8) for n in vec_time]
                tlast_outputfile = vec_time[-1]  
        
                # Read the data in the h5 file            
                with h5py.File(moments_file, 'r') as f:  
                    vec_time_h5file = f['vec_time'][()]
                    upar_vs_tskxky_h5file = f['upar_vs_tskxky'][()]
                    dens_vs_tskxky_h5file = f['dens_vs_tskxky'][()]
                    temp_vs_tskxky_h5file = f['temp_vs_tskxky'][()]
                    upar2_vs_tskxky_h5file = f['upar2_vs_tskxky'][()]
                    dens2_vs_tskxky_h5file = f['dens2_vs_tskxky'][()]
                    temp2_vs_tskxky_h5file = f['temp2_vs_tskxky'][()]
                    upar_vs_tskxky_zeta0_h5file = f['upar_vs_tskxky_zeta0'][()]
                    dens_vs_tskxky_zeta0_h5file = f['dens_vs_tskxky_zeta0'][()]
                    temp_vs_tskxky_zeta0_h5file = f['temp_vs_tskxky_zeta0'][()]
                    tlast_h5file = vec_time_h5file[-1] 
                
                # If the output file (*.out.nc) contains extra time steps, append to the h5 file 
                if tlast_outputfile > tlast_h5file: 
                    index_tlast = np.argwhere(tlast_h5file < vec_time)[0][0] 
                    vec_time, upar_vs_tskxky, dens_vs_tskxky, temp_vs_tskxky, \
                    upar_vs_tskxky_zeta0, dens_vs_tskxky_zeta0, temp_vs_tskxky_zeta0, \
                    upar2_vs_tskxky, dens2_vs_tskxky, temp2_vs_tskxky = read_moments4D(dt, netcdf_file, input_file) 
                    vec_time = np.append(vec_time_h5file, vec_time[index_tlast:], axis=0)  
                    upar_vs_tskxky = np.append(upar_vs_tskxky_h5file, upar_vs_tskxky[index_tlast:,:,:,:], axis=0)  
                    dens_vs_tskxky = np.append(dens_vs_tskxky_h5file, dens_vs_tskxky[index_tlast:,:,:,:], axis=0)  
                    temp_vs_tskxky = np.append(temp_vs_tskxky_h5file, temp_vs_tskxky[index_tlast:,:,:,:], axis=0)  
                    upar2_vs_tskxky = np.append(upar2_vs_tskxky_h5file, upar2_vs_tskxky[index_tlast:,:,:,:], axis=0)  
                    dens2_vs_tskxky = np.append(dens2_vs_tskxky_h5file, dens2_vs_tskxky[index_tlast:,:,:,:], axis=0)  
                    temp2_vs_tskxky = np.append(temp2_vs_tskxky_h5file, temp2_vs_tskxky[index_tlast:,:,:,:], axis=0)  
                    upar_vs_tskxky_zeta0 = np.append(upar_vs_tskxky_zeta0_h5file, upar_vs_tskxky_zeta0[index_tlast:,:,:,:], axis=0)  
                    dens_vs_tskxky_zeta0 = np.append(dens_vs_tskxky_zeta0_h5file, dens_vs_tskxky_zeta0[index_tlast:,:,:,:], axis=0)  
                    temp_vs_tskxky_zeta0 = np.append(temp_vs_tskxky_zeta0_h5file, temp_vs_tskxky_zeta0[index_tlast:,:,:,:], axis=0)  
                    
                # If the output file has the same time steps, touch the text file
                elif tlast_outputfile <= tlast_h5file:
                    print(status+"The 4D moments file is up to date:", moments_file.parent.name+"/"+moments_file.name)  
                    os.system("touch "+str(moments_file))
                    continue
                
            # Otherwise read the 4D moments data
            elif not os.path.isfile(moments_file):    
                
                # Read the potential versus time from the output file
                vec_time, upar_vs_tskxky, dens_vs_tskxky, temp_vs_tskxky, \
                upar_vs_tskxky_zeta0, dens_vs_tskxky_zeta0, temp_vs_tskxky_zeta0, \
                upar2_vs_tskxky, dens2_vs_tskxky, temp2_vs_tskxky = read_moments4D(dt, netcdf_file, input_file) 
                
            # Create the new h5 file 
            with h5py.File(moments_file, 'w') as h5_file:   
                h5_file.create_dataset("vec_time", data=vec_time) 
                h5_file.create_dataset("upar_vs_tskxky", data=upar_vs_tskxky) 
                h5_file.create_dataset("dens_vs_tskxky", data=dens_vs_tskxky)  
                h5_file.create_dataset("temp_vs_tskxky", data=temp_vs_tskxky)  
                h5_file.create_dataset("upar2_vs_tskxky", data=upar2_vs_tskxky) 
                h5_file.create_dataset("dens2_vs_tskxky", data=dens2_vs_tskxky)  
                h5_file.create_dataset("temp2_vs_tskxky", data=temp2_vs_tskxky)  
                h5_file.create_dataset("upar_vs_tskxky_zeta0", data=upar_vs_tskxky_zeta0) 
                h5_file.create_dataset("dens_vs_tskxky_zeta0", data=dens_vs_tskxky_zeta0)  
                h5_file.create_dataset("temp_vs_tskxky_zeta0", data=temp_vs_tskxky_zeta0)  
                                     
            # Notify that we finished creating the file
            if outputFileHasChanged: print(status+"   ---> The 4D moments file is updated as " +  moments_file.parent.name+"/"+moments_file.name)   
            if not outputFileHasChanged:  print(status+"   ---> The 4D moments file is saved as " +  moments_file.parent.name+"/"+moments_file.name) 
            
    return 

#---------------------------------------------
def read_moments4D(dt, netcdf_file, input_file):              

    # Read the geometry data in the output file
    vmec_filename = read_vmecFileNameFromInputFile(input_file)
    nonlinear = read_linearNonlinearFromInputFile(input_file)[1]
    path = create_dummyPathObject(input_file, vmec_filename, nonlinear)
    geometry = read_outputFileForGeometry(path)  
    dl_over_B = geometry["dl_over_B"]   
    
    # Get the moments from the output file
    netcdf_file = read_outputFile(netcdf_file)  
    vec_time = read_netcdfVariables('vec_time', netcdf_file)  
    vec_z = read_netcdfVariables('vec_z', netcdf_file) 
    indices = get_indicesAtFixedStep(vec_time, dt)
    upar_vs_tszkxkyri = read_netcdfVariables('upar_vs_tszkxkyri', netcdf_file, indices)
    dens_vs_tszkxkyri = read_netcdfVariables('dens_vs_tszkxkyri', netcdf_file, indices)
    temp_vs_tszkxkyri = read_netcdfVariables('temp_vs_tszkxkyri', netcdf_file, indices) 
    netcdf_file.close()
    
    # Get the data at every <dt> timestep 
    upar_vs_tszkxky = upar_vs_tszkxkyri[:,:,:,:,:,0] + 1j*upar_vs_tszkxkyri[:,:,:,:,:,1]
    dens_vs_tszkxky = dens_vs_tszkxkyri[:,:,:,:,:,0] + 1j*dens_vs_tszkxkyri[:,:,:,:,:,1]
    temp_vs_tszkxky = temp_vs_tszkxkyri[:,:,:,:,:,0] + 1j*temp_vs_tszkxkyri[:,:,:,:,:,1]
    vec_time = vec_time[indices]
    
    # Sum away the (z) dimension
    upar_vs_tskxky = np.sum(upar_vs_tszkxky[:,:,:,:,:]*dl_over_B[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis], axis=2)
    dens_vs_tskxky = np.sum(dens_vs_tszkxky[:,:,:,:,:]*dl_over_B[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis], axis=2)
    temp_vs_tskxky = np.sum(temp_vs_tszkxky[:,:,:,:,:]*dl_over_B[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis], axis=2) 
    
    # Get the data where zeta=0
    zeta0 = np.where(vec_z==0.0)[0][0]
    upar_vs_tskxky_zeta0 = upar_vs_tszkxky[:,:,zeta0,:,:]
    dens_vs_tskxky_zeta0 = dens_vs_tszkxky[:,:,zeta0,:,:]
    temp_vs_tskxky_zeta0 = temp_vs_tszkxky[:,:,zeta0,:,:]
    
    # Careful: for dens^2 we need to first take abs(dens)**2 and then sum away the dimensions (z)!
    upar2_vs_tszkxky = np.abs(upar_vs_tszkxky)**2
    dens2_vs_tszkxky = np.abs(dens_vs_tszkxky)**2
    temp2_vs_tszkxky = np.abs(temp_vs_tszkxky)**2
    upar2_vs_tskxky = np.sum(upar2_vs_tszkxky[:,:,:,:,:]*dl_over_B[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis], axis=1)
    dens2_vs_tskxky = np.sum(dens2_vs_tszkxky[:,:,:,:,:]*dl_over_B[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis], axis=1)
    temp2_vs_tskxky = np.sum(temp2_vs_tszkxky[:,:,:,:,:]*dl_over_B[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis], axis=1)
    return vec_time, upar_vs_tskxky, dens_vs_tskxky, temp_vs_tskxky, \
        upar_vs_tskxky_zeta0, dens_vs_tskxky_zeta0, temp_vs_tskxky_zeta0, \
        upar2_vs_tskxky, dens2_vs_tskxky, temp2_vs_tskxky
    
    