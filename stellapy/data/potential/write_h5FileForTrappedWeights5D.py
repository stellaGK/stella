  
import numpy as np
import os, h5py, pathlib 
from stellapy.utils.decorators.verbose import noverbose 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder    
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables 
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep 

################################################################################
#                WRITE THE TRAPPED PARTICLE WEIGHTS TO AN H5 FILE
################################################################################

@noverbose
def write_h5FileForTrappedWeights5D(folder, dt=1):  
    
    # Time step
    dt = int(dt) if int(dt)==dt else dt   
     
    # Get the input files for which a netcdf file exists
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 
 
    # Iterate through the input files  
    for input_file in input_files:   
        
        # Processing status
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
            
        # Path of the new file 
        potential_file = input_file.with_suffix(".dt"+str(dt)+".trappedWeights5D") 
        netcdf_file = input_file.with_suffix(".out.nc")  
        
        # If the file doesn't exist, check whether the weights were written
        if not os.path.isfile(potential_file): 
            netcdf_data  = read_outputFile(netcdf_file)  
            if "trappedw_vs_kxkyz" not in netcdf_data.variables.keys():
                continue
        
        # Check whether the potential file is older than the simulation
        outputFileHasChanged = False
        if os.path.isfile(potential_file):  
            if netcdf_file.stat().st_mtime > potential_file.stat().st_mtime:
                outputFileHasChanged = True 
                
        # Notify that the file already existed   
        if os.path.isfile(potential_file) and not outputFileHasChanged:
            print(status+"The trapped weights file already exists:", potential_file.parent.name+"/"+potential_file.name)
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
                trappedw_vs_tkxkyz_h5file = f['trappedw_vs_tszkxky'][()] 
                tlast_h5file = vec_time_h5file[-1]
            
            # If the output file (*.out.nc) contains extra time steps, append to the h5 file 
            if tlast_outputfile > tlast_h5file: 
                index_tlast = np.argwhere(tlast_h5file < vec_time)[0][0] 
                vec_time, trappedw_vs_tszkxky = read_trappedWeights5D(dt, netcdf_file) 
                vec_time = np.append(vec_time_h5file, vec_time[index_tlast:], axis=0)  
                trappedw_vs_tszkxky = np.append(trappedw_vs_tkxkyz_h5file, trappedw_vs_tszkxky[index_tlast:,:,:,:], axis=0)   
                
            # If the output file has the same time steps, touch the text file
            elif tlast_outputfile <= tlast_h5file:
                print(status+"The trapped weights file is up to date:", potential_file.parent.name+"/"+potential_file.name)  
                os.system("touch "+str(potential_file))
                continue
            
        # Create the text file for dphiz(t)
        elif not os.path.isfile(potential_file):    
            
            # Read the trapped particle weights versus time from the output file
            vec_time, trappedw_vs_tszkxky = read_trappedWeights5D(dt, netcdf_file) 
                 
        # Create the new h5 file 
        with h5py.File(potential_file, 'w') as h5_file: 
            h5_file.create_dataset("vec_time", data=vec_time) 
            h5_file.create_dataset("trappedw_vs_tszkxky", data=trappedw_vs_tszkxky)  
                                 
        # Notify that we finished creating the file
        if outputFileHasChanged: print(status+"   ---> The trapped weights file is updated as " +  potential_file.parent.name+"/"+potential_file.name)   
        if not outputFileHasChanged:  print(status+"   ---> The trapped weights file is saved as " +  potential_file.parent.name+"/"+potential_file.name)  
        
    return 

#---------------------------------------------
def read_trappedWeights5D(dt, netcdf_file):
                 
    # Read the potential from the output file
    netcdf_data = read_outputFile(netcdf_file)  
    vec_time = read_netcdfVariables('vec_time', netcdf_data) 
    trappedw_vs_tszkxky = read_netcdfVariables('trappedw_vs_tszkxky', netcdf_data) 
    netcdf_data.close()
    
    # Get the data at every <dt> timestep 
    indices = get_indicesAtFixedStep(vec_time, dt)
    trappedw_vs_tszkxky = trappedw_vs_tszkxky[indices,:,:,:,:] 
    vec_time = vec_time[indices]

    return vec_time, trappedw_vs_tszkxky
                     
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")   
    write_h5FileForTrappedWeights5D(folder) 
    