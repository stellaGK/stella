  
import numpy as np
import os, h5py, pathlib 
from stellapy.utils.decorators.verbose import noverbose 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder   
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep  
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables  
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile 

################################################################################
#                       WRITE THE FLUXES TO AN H5 FILE
################################################################################

@noverbose
def write_fluxesVsZ(folder, dt=1):  
    
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
        
        # Only write the fluxes(t,kx,ky) file for nonlinear simulations 
        if read_linearNonlinearFromInputFile(input_file)[1]:    
            
            # Path of the new file 
            fluxes_file = input_file.with_suffix(".dt"+str(dt)+".fluxes2D") 
            netcdf_file = input_file.with_suffix(".out.nc")  
            
            # If the file doesn't exist, check whether the fluxes were written to the netcdf file
            # Also the fluxes need to contain the z-dimension, which has been added in december 2021
            if not os.path.isfile(fluxes_file): 
                netcdf_data  = read_outputFile(netcdf_file)   
                if "qflx_kxky" not in netcdf_data.variables.keys():
                    qflux_vs_tskxky = read_netcdfVariables('qflux_vs_tskxky', netcdf_data) 
                    if len(np.shape(qflux_vs_tskxky))!=6: 
                        continue 
                
            # Check whether the fluxes file is older than the simulation
            outputFileHasChanged = False
            if os.path.isfile(fluxes_file):  
                if netcdf_file.stat().st_mtime > fluxes_file.stat().st_mtime:
                    outputFileHasChanged = True 
                    
            # Notify that the file already existed   
            if os.path.isfile(fluxes_file) and not outputFileHasChanged:
                print(status+"The 2D fluxes file already exists:", fluxes_file.parent.name+"/"+fluxes_file.name)
                continue
                 
            # If the output file changed, then append to the h5 file
            elif os.path.isfile(fluxes_file) and outputFileHasChanged:
                
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
                with h5py.File(fluxes_file, 'r') as f: 
                    vec_time_h5file = f['vec_time'][()]
                    qflux_vs_tsz_h5file = f['qflux_vs_tsz'][()]
                    pflux_vs_tsz_h5file = f['pflux_vs_tsz'][()]
                    vflux_vs_tsz_h5file = f['vflux_vs_tsz'][()] 
                    tlast_h5file = vec_time_h5file[-1]
                
                # If the output file (*.out.nc) contains extra time steps, append to the h5 file 
                if tlast_outputfile > tlast_h5file: 
                    index_tlast = np.argwhere(tlast_h5file < vec_time)[0][0] 
                    vec_time, qflux_vs_tsz, pflux_vs_tsz, vflux_vs_tsz = read_fluxes2D(dt, netcdf_file) 
                    vec_time = np.append(vec_time_h5file, vec_time[index_tlast:], axis=0)  
                    qflux_vs_tsz = np.append(qflux_vs_tsz_h5file, qflux_vs_tsz[index_tlast:,:,:], axis=0)  
                    pflux_vs_tsz = np.append(pflux_vs_tsz_h5file, pflux_vs_tsz[index_tlast:,:,:], axis=0)  
                    vflux_vs_tsz = np.append(vflux_vs_tsz_h5file, vflux_vs_tsz[index_tlast:,:,:], axis=0)  
                    
                # If the output file has the same time steps, touch the text file
                elif tlast_outputfile <= tlast_h5file:
                    print(status+"The 2D fluxes file is up to date:", fluxes_file.parent.name+"/"+fluxes_file.name)  
                    os.system("touch "+str(fluxes_file))
                    continue
                
            # Create the text file for dphiz(t)
            elif not os.path.isfile(fluxes_file):    
                
                # Read the fluxes versus time from the output file
                vec_time, qflux_vs_tsz, pflux_vs_tsz, vflux_vs_tsz = read_fluxes2D(dt, netcdf_file) 
                     
            # Create the new h5 file 
            with h5py.File(fluxes_file, 'w') as h5_file: 
                h5_file.create_dataset("vec_time", data=vec_time) 
                h5_file.create_dataset("qflux_vs_tsz", data=qflux_vs_tsz) 
                h5_file.create_dataset("pflux_vs_tsz", data=pflux_vs_tsz) 
                h5_file.create_dataset("vflux_vs_tsz", data=vflux_vs_tsz)  
                                     
            # Notify that we finished creating the file
            if outputFileHasChanged: print(status+"   ---> The 2D fluxes file is updated as " +  fluxes_file.parent.name+"/"+fluxes_file.name)   
            if not outputFileHasChanged:  print(status+"   ---> The 2D fluxes file is saved as " +  fluxes_file.parent.name+"/"+fluxes_file.name)  
        
    return 

#---------------------------------------------
def read_fluxes2D(dt, netcdf_file):
                 
    # Read the fluxes from the output file
    netcdf_data = read_outputFile(netcdf_file)  
    vec_time = read_netcdfVariables('vec_time', netcdf_data) 
    qflux_vs_tszkxky = read_netcdfVariables('qflux_vs_tszkxky', netcdf_data) 
    pflux_vs_tszkxky = read_netcdfVariables('pflux_vs_tszkxky', netcdf_data) 
    vflux_vs_tszkxky = read_netcdfVariables('vflux_vs_tszkxky', netcdf_data)  
    netcdf_data.close()  
        
    # Get the data at every <dt> timestep 
    indices = get_indicesAtFixedStep(vec_time, dt)
    qflux_vs_tszkxky = qflux_vs_tszkxky[indices,:,:,:,:] 
    pflux_vs_tszkxky = pflux_vs_tszkxky[indices,:,:,:,:]
    vflux_vs_tszkxky = vflux_vs_tszkxky[indices,:,:,:,:]
    vec_time = vec_time[indices]
    
    # Get flux(t,z) by summing over the (kx,ky) components
    qflux_vs_tsz = np.sum(qflux_vs_tszkxky[:,:,:,:,0] + 2*np.sum(qflux_vs_tszkxky[:,:,:,:,1:],axis=4),axis=3)   
    pflux_vs_tsz = np.sum(pflux_vs_tszkxky[:,:,:,:,0] + 2*np.sum(pflux_vs_tszkxky[:,:,:,:,1:],axis=4),axis=3)   
    vflux_vs_tsz = np.sum(vflux_vs_tszkxky[:,:,:,:,0] + 2*np.sum(vflux_vs_tszkxky[:,:,:,:,1:],axis=4),axis=3)   
    return vec_time, qflux_vs_tsz, pflux_vs_tsz, vflux_vs_tsz
                     
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")   
    write_fluxesVsZ(folder) 
    