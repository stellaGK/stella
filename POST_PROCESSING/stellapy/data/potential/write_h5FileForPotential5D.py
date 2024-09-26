   
import numpy as np
import os, h5py, copy, pathlib
from datetime import datetime, timedelta
from stellapy.utils.decorators.verbose import noverbose 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder   
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables  
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep

################################################################################
#                       WRITE THE POTENTIAL TO AN H5 FILE
################################################################################

@noverbose
def write_h5FileForPotential5D(folder, dt=100):  
      
    # Time step
    dt = int(dt) if int(dt)==dt else dt      
     
    # Get the input files for which a netcdf file exists
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 
 
    # Iterate through the input files  
    for input_file in get_filesInFolder(folder, end=".in"):  
        
        # Only write the phi(t,z,kx,ky) file for nonlinear simulations 
        if read_linearNonlinearFromInputFile(input_file)[1]: 
        
            # Processing status
            status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
            
            # Path of the new file
            netcdf_file = input_file.with_suffix(".out.nc")
            potential_file = input_file.with_suffix(".dt"+str(dt)+".potential5D")     
                    
            # Check whether the potential file is older than the simulation
            outputFileHasChanged = False
            if os.path.isfile(potential_file):  
                if datetime.fromtimestamp(netcdf_file.stat().st_mtime) > (datetime.fromtimestamp(potential_file.stat().st_mtime)+timedelta(minutes=5)):
                    outputFileHasChanged = True 
                    
            # Notify that the file already existed   
            if os.path.isfile(potential_file) and not outputFileHasChanged:
                print(status+"The 5D potential file already exists:", potential_file.parent.name+"/"+potential_file.name)
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
                    phi_vs_tzkxky_h5file = f['phi_vs_tzkxky'][()]
                    phi_vs_tzkxky_tavg_h5file = f['phi_vs_tzkxky_tavg'][()] 
                    tlast_h5file = vec_time_h5file[-1] 
                
                # If the output file (*.out.nc) contains extra time steps, append to the h5 file 
                if tlast_outputfile > tlast_h5file: 
                    index_tlast = np.argwhere(tlast_h5file < vec_time)[0][0] 
                    vec_time, phi_vs_tzkxky, phi_vs_tzkxky_tavg = read_potential5D(dt, netcdf_file) 
                    vec_time = np.append(vec_time_h5file, vec_time[index_tlast:], axis=0)  
                    phi_vs_tzkxky = np.append(phi_vs_tzkxky_h5file, phi_vs_tzkxky[index_tlast:,:,:,:], axis=0)  
                    phi_vs_tzkxky_tavg = np.append(phi_vs_tzkxky_tavg_h5file, phi_vs_tzkxky_tavg[index_tlast:,:,:,:], axis=0)   
                    
                # If the output file has the same time steps, touch the text file
                elif tlast_outputfile <= tlast_h5file:
                    print(status+"The 5D potential file is up to date:", potential_file.parent.name+"/"+potential_file.name)  
                    os.system("touch "+str(potential_file))
                    continue
                
            # Otherwise read the 5D potential data
            elif not os.path.isfile(potential_file):    
                
                # Read the potential versus time from the output file
                vec_time, phi_vs_tzkxky, phi_vs_tzkxky_tavg = read_potential5D(dt, netcdf_file) 
                
            # Create the new h5 file 
            with h5py.File(potential_file, 'w') as h5_file:   
                h5_file.create_dataset("vec_time", data=vec_time) 
                h5_file.create_dataset("phi_vs_tzkxky", data=phi_vs_tzkxky) 
                h5_file.create_dataset("phi_vs_tzkxky_tavg", data=phi_vs_tzkxky_tavg)  
                                     
            # Notify that we finished creating the file
            if outputFileHasChanged: print(status+"   ---> The 5D potential file is updated as " +  potential_file.parent.name+"/"+potential_file.name)   
            if not outputFileHasChanged:  print(status+"   ---> The 5D potential file is saved as " +  potential_file.parent.name+"/"+potential_file.name)  

    return 

#---------------------------------------------
def read_potential5D(dt, netcdf_file):      
    
    # Read the potential from the output file
    netcdf_data = read_outputFile(netcdf_file)  
    vec_time_full = read_netcdfVariables('vec_time', netcdf_data) 
    phi_vs_tzkxkyri = read_netcdfVariables('phi_vs_tzkxkyri', netcdf_data)  
    phi_vs_tzkxky_full = phi_vs_tzkxkyri[:,:,:,:,0] + 1j*phi_vs_tzkxkyri[:,:,:,:,1]
    netcdf_data.close()
    
    # Get the data at every <dt> timestep 
    indices = get_indicesAtFixedStep(vec_time_full, dt) 
    phi_vs_tzkxky = copy.deepcopy(phi_vs_tzkxky_full)[indices,:,:,:] 
    vec_time = copy.deepcopy(vec_time_full)[indices] 
    
    # Average the data over [(0,dt), (dt,2dt), (2dt, 3dt), ...]
    phi_vs_tzkxky_tavg = []; vec_time[-1] += 1000
    for i in range(len(vec_time)-1):
        time_filter = (vec_time_full >= vec_time[i]) & (vec_time_full < vec_time[i+1])
        phi_vs_tzkxky_tavg.append(np.nanmean(phi_vs_tzkxky_full[time_filter,:,:,:], axis=0))
        
    return vec_time, phi_vs_tzkxky, phi_vs_tzkxky_tavg
    
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")   
    write_h5FileForPotential5D(folder) 
    