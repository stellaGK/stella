  
import numpy as np
import os, h5py, pathlib 
from datetime import datetime, timedelta 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder   
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile, read_fullFluxSurfaceFromInputFile
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables 
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry  
from stellapy.data.input.read_inputFile import read_modeFromInputFile
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile

################################################################################
#                       WRITE THE POTENTIAL TO AN H5 FILE
################################################################################
 
def write_h5FileForPotential3D(folder, dt=10, automatic=False):   
    
    # Time step
    dt = int(dt) if int(dt)==dt else dt  
     
    # Get the input files for which a netcdf file exists
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 
 
    # Iterate through the input files  
    for input_file in input_files:  
        
        try:
        
            # Processing status
            status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
            nonlinear = read_linearNonlinearFromInputFile(input_file)[1] 
            if (not nonlinear) and automatic: continue
                
            # Path of the new file 
            potential_file = input_file.with_suffix(".dt"+str(dt)+".potential3D") 
            netcdf_file = input_file.with_suffix(".out.nc")  
            
            # Check whether the potential file is older than the simulation
            outputFileHasChanged = False
            if os.path.isfile(potential_file):  
                if datetime.fromtimestamp(netcdf_file.stat().st_mtime) > (datetime.fromtimestamp(potential_file.stat().st_mtime)+timedelta(minutes=5)):
                    outputFileHasChanged = True 
                    
            # Notify that the file already existed   
            if os.path.isfile(potential_file) and not outputFileHasChanged:
                print(status+"The 3D potential file already exists:", potential_file.parent.name+"/"+potential_file.name)
                continue
                 
            # If the output file changed, then append to the h5 file
            elif os.path.isfile(potential_file) and outputFileHasChanged:
                
                # Check whether the output file (*.out.nc) contains extra time points
                netcdf_data  = read_outputFile(netcdf_file)   
                vec_time     = read_netcdfVariables('vec_time', netcdf_data) 
                netcdf_data.close()  
                
                # Edit the time vector in the output file
                indices  = get_indicesAtFixedStep(vec_time, dt)
                vec_time = vec_time[indices]
                vec_time = [round(n, 8) for n in vec_time]
                tlast_outputfile = vec_time[-1]  
        
                # Read the data in the h5 file            
                with h5py.File(potential_file, 'r') as f: 
                    try:
                        vec_time_h5file = f['vec_time'][()]
                        phi_vs_tz_h5file = f['phi_vs_tz'][()]
                        phi2_vs_tz_h5file = f['phi2_vs_tz'][()]
                        tlast_h5file = vec_time_h5file[-1]
                        if nonlinear:
                            phi_vs_tz_zonal_h5file = f['phi_vs_tz_zonal'][()]
                            phi_vs_tz_nozonal_h5file = f['phi_vs_tz_nozonal'][()]
                            phi2_vs_tz_zonal_h5file = f['phi2_vs_tz_zonal'][()]
                            phi2_vs_tz_nozonal_h5file = f['phi2_vs_tz_nozonal'][()]
                            phi2_vs_tkx_h5file = f['phi2_vs_tkx'][()]
                            phi2_vs_tky_h5file = f['phi2_vs_tky'][()]
                            phi_vs_tkx_h5file = f['phi_vs_tkx'][()]
                            phi_vs_tky_h5file = f['phi_vs_tky'][()]
                    except:
                        print(status+"Data was missing in the already existing potential 3D file: ")
                        print("                 "+str(potential_file))
                        print("                    ---> Removed this file, please rerun <write_dataFiles>.")
                        os.system("rm "+str(potential_file))
                        continue 
                
                # If the output file (*.out.nc) contains extra time steps, append to the h5 file 
                if tlast_outputfile > tlast_h5file: 
                    index_tlast = np.argwhere(tlast_h5file < vec_time)[0][0] 
                    vec_time, phi_vs_tz, phi_vs_tz_zonal, phi_vs_tz_nozonal, phi2_vs_tz, phi2_vs_tz_zonal,\
                    phi2_vs_tz_nozonal, phi2_vs_tkx, phi2_vs_tky, phi_vs_tkx, phi_vs_tky = read_potential3D(dt, netcdf_file, input_file) 
                    vec_time = np.append(vec_time_h5file, vec_time[index_tlast:], axis=0)  
                    phi_vs_tz = np.append(phi_vs_tz_h5file, phi_vs_tz[index_tlast:,:], axis=0)  
                    phi_vs_tz_zonal = np.append(phi_vs_tz_zonal_h5file, phi_vs_tz_zonal[index_tlast:,:], axis=0)  
                    phi_vs_tz_nozonal = np.append(phi_vs_tz_nozonal_h5file, phi_vs_tz_nozonal[index_tlast:,:], axis=0) 
                    phi2_vs_tz = np.append(phi2_vs_tz_h5file, phi2_vs_tz[index_tlast:,:], axis=0)  
                    phi2_vs_tz_zonal = np.append(phi2_vs_tz_zonal_h5file, phi2_vs_tz_zonal[index_tlast:,:], axis=0)  
                    phi2_vs_tz_nozonal = np.append(phi2_vs_tz_nozonal_h5file, phi2_vs_tz_nozonal[index_tlast:,:], axis=0) 
                    phi2_vs_tkx = np.append(phi2_vs_tkx_h5file, phi2_vs_tkx[index_tlast:,:], axis=0)  
                    phi2_vs_tky = np.append(phi2_vs_tky_h5file, phi2_vs_tky[index_tlast:,:], axis=0) 
                    phi_vs_tkx = np.append(phi_vs_tkx_h5file, phi_vs_tkx[index_tlast:,:], axis=0)  
                    phi_vs_tky = np.append(phi_vs_tky_h5file, phi_vs_tky[index_tlast:,:], axis=0)   
                    
                # Otherwise read the 3D potential data
                elif tlast_outputfile <= tlast_h5file:
                    print(status+"The 3D potential file is up to date:", potential_file.parent.name+"/"+potential_file.name)  
                    os.system("touch "+str(potential_file))
                    continue
                
            # Create the h5 file for the 3D potential data
            elif not os.path.isfile(potential_file):    
                
                # Read the potential versus time from the output file
                vec_time, phi_vs_tz, phi_vs_tz_zonal, phi_vs_tz_nozonal, phi2_vs_tz, phi2_vs_tz_zonal,\
                phi2_vs_tz_nozonal, phi2_vs_tkx, phi2_vs_tky, phi_vs_tkx, phi_vs_tky = read_potential3D(dt, netcdf_file, input_file) 
                     
            # Create the new h5 file 
            with h5py.File(potential_file, 'w') as h5_file: 
                h5_file.create_dataset("vec_time", data=vec_time) 
                h5_file.create_dataset("phi_vs_tz", data=phi_vs_tz) 
                h5_file.create_dataset("phi2_vs_tz", data=phi2_vs_tz)  
                if nonlinear:
                    h5_file.create_dataset("phi_vs_tz_zonal", data=phi_vs_tz_zonal)  
                    h5_file.create_dataset("phi_vs_tz_nozonal", data=phi_vs_tz_nozonal)  
                    h5_file.create_dataset("phi2_vs_tz_zonal", data=phi2_vs_tz_zonal)  
                    h5_file.create_dataset("phi2_vs_tz_nozonal", data=phi2_vs_tz_nozonal)
                    h5_file.create_dataset("phi2_vs_tkx", data=phi2_vs_tkx)  
                    h5_file.create_dataset("phi2_vs_tky", data=phi2_vs_tky)  
                    h5_file.create_dataset("phi_vs_tkx", data=phi_vs_tkx)  
                    h5_file.create_dataset("phi_vs_tky", data=phi_vs_tky)  
                                     
            # Notify that we finished creating the file
            if outputFileHasChanged: print(status+"   ---> The 3D potential file is updated as " +  potential_file.parent.name+"/"+potential_file.name)   
            if not outputFileHasChanged:  print(status+"   ---> The 3D potential file is saved as " +  potential_file.parent.name+"/"+potential_file.name)  
        
        except: 
            print(status+"Something went wrong for:", input_file.parent.parent.name+"/"+input_file.parent.name+"/"+input_file.name)
            import sys; sys.exit()
    return 



#---------------------------------------------
def read_potential3D(dt, netcdf_file, input_file):
    
    # Read the geometry data in the output file
    vmec_filename = read_vmecFileNameFromInputFile(input_file) 
    path = create_dummyPathObject(input_file, vmec_filename)
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
    
    # Get phi(t,z) by summing away the (kx,ky) modes  
    phi_vs_tz_zonal = np.sum(phi_vs_tzkxky[:,:,:,0],axis=2)
    phi_vs_tz_nozonal = 2*np.sum(phi_vs_tzkxky[:,:,:,1:],axis=(2,3))
    phi_vs_tz = phi_vs_tz_zonal + phi_vs_tz_nozonal
    
    # Carefull: for phi2 we need to first take abs(phi)**2 and then sum away the dimensions (kx,ky)!
    phi2_vs_tzkxky = np.abs(phi_vs_tzkxky)**2
    phi2_vs_tz_zonal = np.sum(phi2_vs_tzkxky[:,:,:,0],axis=2)
    phi2_vs_tz_nozonal = 2*np.sum(phi2_vs_tzkxky[:,:,:,1:],axis=(2,3))  
    phi2_vs_tz = phi2_vs_tz_zonal + phi2_vs_tz_nozonal
    
    # Get phi2(t,kx) and phi2(t,ky)
    phi2_vs_tkxky = np.sum(phi2_vs_tzkxky[:,:,:,:]*dl_over_B[np.newaxis,:,np.newaxis,np.newaxis], axis=1)
    phi2_vs_tkx = phi2_vs_tkxky[:,:,0] + 2*np.sum(phi2_vs_tkxky[:,:,1:],axis=2)  
    phi2_vs_tky = np.sum(phi2_vs_tkxky[:,:,:],axis=1)  
    
    # Also get phi(t,kx) and phi(t,ky)
    phi_vs_tkxky = np.sum(phi_vs_tzkxky[:,:,:,:]*dl_over_B[np.newaxis,:,np.newaxis,np.newaxis], axis=1)
    phi_vs_tkx = phi_vs_tkxky[:,:,0] + 2*np.sum(phi_vs_tkxky[:,:,1:],axis=2)  
    phi_vs_tky = np.sum(phi_vs_tkxky[:,:,:],axis=1) 
    
    # For a linear mode, we need to sum the non-zonal modes twice
    if read_linearNonlinearFromInputFile(input_file)[0]: 
        if read_modeFromInputFile(input_file)[1]!=0: 
            phi_vs_tz = phi_vs_tz*2
            phi2_vs_tz = phi2_vs_tz*2 
            
    return vec_time, phi_vs_tz, phi_vs_tz_zonal, phi_vs_tz_nozonal, phi2_vs_tz, \
        phi2_vs_tz_zonal, phi2_vs_tz_nozonal, phi2_vs_tkx, phi2_vs_tky, phi_vs_tkx, phi_vs_tky
                     
    
