  
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
#                       WRITE THE POTENTIAL TO AN H5 FILE
################################################################################

@noverbose
def write_h5FileForPotential4D(folder, dt=10):   
     
    # Time step
    dt = int(dt) if int(dt)==dt else dt     
     
    # Get the input files for which a netcdf file exists
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 
 
    # Iterate through the input files  
    for input_file in input_files:  
        
        # Only write the phi(t,kx,ky) file for nonlinear simulations 
        if read_linearNonlinearFromInputFile(input_file)[1]: 
        
            # Processing status
            status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
            
            # Path of the new file
            netcdf_file = input_file.with_suffix(".out.nc")
            potential_file = input_file.with_suffix(".dt"+str(dt)+".potential4D")   
                    
            # Check whether the potential file is older than the simulation
            outputFileHasChanged = False
            if os.path.isfile(potential_file):  
                if datetime.fromtimestamp(netcdf_file.stat().st_mtime) > (datetime.fromtimestamp(potential_file.stat().st_mtime)+timedelta(minutes=5)):
                    outputFileHasChanged = True 
                with h5py.File(potential_file, 'r') as h5_file:   
                    if "phi_vs_tkxky_zeta0" not in h5_file.keys():
                        print("Remove the 4D potential file since <phi_vs_tkxky_zeta0> was missing.")
                        os.system("rm "+str(potential_file))
                    
            # Notify that the file already existed   
            if os.path.isfile(potential_file) and not outputFileHasChanged:
                print(status+"The 4D potential file already exists:", potential_file.parent.name+"/"+potential_file.name)
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
                    vec_time_h5file = f['vec_time'][()]
                    phi_vs_tkxky_h5file = f['phi_vs_tkxky'][()]
                    phi2_vs_tkxky_h5file = f['phi2_vs_tkxky'][()] 
                    phi_vs_tkxky_zeta0_h5file = f['phi_vs_tkxky_zeta0'][()] 
                    tlast_h5file = vec_time_h5file[-1] 
                
                # If the output file (*.out.nc) contains extra time steps, append to the h5 file 
                if tlast_outputfile > tlast_h5file: 
                    index_tlast = np.argwhere(tlast_h5file < vec_time)[0][0] 
                    vec_time, phi_vs_tkxky, phi_vs_tkxky_zeta0, phi2_vs_tkxky = read_potential4D(dt, netcdf_file, input_file) 
                    vec_time = np.append(vec_time_h5file, vec_time[index_tlast:], axis=0)  
                    phi_vs_tkxky = np.append(phi_vs_tkxky_h5file, phi_vs_tkxky[index_tlast:,:,:], axis=0)  
                    phi2_vs_tkxky = np.append(phi2_vs_tkxky_h5file, phi2_vs_tkxky[index_tlast:,:,:], axis=0)   
                    phi_vs_tkxky_zeta0 = np.append(phi_vs_tkxky_zeta0_h5file, phi_vs_tkxky_zeta0[index_tlast:,:,:], axis=0)   
                    
                # If the output file has the same time steps, touch the text file
                elif tlast_outputfile <= tlast_h5file:
                    print(status+"The 4D potential file is up to date:", potential_file.parent.name+"/"+potential_file.name)  
                    os.system("touch "+str(potential_file))
                    continue
                
            # Otherwise read the 4D potential data
            elif not os.path.isfile(potential_file):    
                
                # Read the potential versus time from the output file
                vec_time, phi_vs_tkxky, phi_vs_tkxky_zeta0, phi2_vs_tkxky = read_potential4D(dt, netcdf_file, input_file) 
                
            # Create the new h5 file 
            with h5py.File(potential_file, 'w') as h5_file:   
                h5_file.create_dataset("vec_time", data=vec_time) 
                h5_file.create_dataset("phi_vs_tkxky", data=phi_vs_tkxky) 
                h5_file.create_dataset("phi2_vs_tkxky", data=phi2_vs_tkxky)  
                h5_file.create_dataset("phi_vs_tkxky_zeta0", data=phi_vs_tkxky_zeta0)  
                                     
            # Notify that we finished creating the file
            if outputFileHasChanged: print(status+"   ---> The 4D potential file is updated as " +  potential_file.parent.name+"/"+potential_file.name)   
            if not outputFileHasChanged:  print(status+"   ---> The 4D potential file is saved as " +  potential_file.parent.name+"/"+potential_file.name)  
                        
    return 

#---------------------------------------------
def read_potential4D(dt, netcdf_file, input_file):                 

    # Read the geometry data in the output file
    vmec_filename = read_vmecFileNameFromInputFile(input_file)
    nonlinear = read_linearNonlinearFromInputFile(input_file)[1]
    path = create_dummyPathObject(input_file, vmec_filename, nonlinear)
    geometry = read_outputFileForGeometry(path)  
    dl_over_B = geometry["dl_over_B"]  
    
    # Read the potential from the output file
    netcdf_data = read_outputFile(netcdf_file)  
    vec_time = read_netcdfVariables('vec_time', netcdf_data) 
    vec_z = read_netcdfVariables('vec_z', netcdf_data) 
    indices = get_indicesAtFixedStep(vec_time, dt)
    phi_vs_tzkxkyri = read_netcdfVariables('phi_vs_tzkxkyri', netcdf_data, indices) 
    netcdf_data.close()
    
    # Get the data at every <dt> timestep 
    phi_vs_tzkxky = phi_vs_tzkxkyri[:,:,:,:,0] + 1j*phi_vs_tzkxkyri[:,:,:,:,1]
    vec_time = vec_time[indices]
    
    # Sum away the (z) dimension
    phi_vs_tkxky = np.sum(phi_vs_tzkxky[:,:,:,:]*dl_over_B[np.newaxis,:,np.newaxis,np.newaxis], axis=1)
    
    # Get the data where zeta=0
    zeta0 = np.where(vec_z==0.0)[0][0]
    phi_vs_tkxky_zeta0 = phi_vs_tzkxky[:,zeta0,:,:]
    
    # Careful: for phi2 we need to first take abs(phi)**2 and then sum away the dimensions (z)!
    phi2_vs_tzkxky = np.abs(phi_vs_tzkxky)**2
    phi2_vs_tkxky = np.sum(phi2_vs_tzkxky[:,:,:,:]*dl_over_B[np.newaxis,:,np.newaxis,np.newaxis], axis=1)
    return vec_time, phi_vs_tkxky, phi_vs_tkxky_zeta0, phi2_vs_tkxky
    
    