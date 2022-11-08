   
import os, h5py, pathlib
from datetime import datetime, timedelta
from stellapy.utils.decorators.verbose import noverbose 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder     
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables

################################################################################
#                       WRITE THE DIMENSIONS TO AN H5 FILE
################################################################################

@noverbose
def write_h5FileForMoments5D(folder, dt=10):   
     
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
            output_file = input_file.with_suffix(".out")
            netcdf_file = input_file.with_suffix(".out.nc")
            moments_file = input_file.with_suffix(".dt"+str(dt)+".moments5D")   
        
            # Check whether the txt file is older than the simulation
            if os.path.isfile(moments_file): 
                time_txtFile = datetime.fromtimestamp(moments_file.stat().st_mtime)
                time_simulation = datetime.fromtimestamp(output_file.stat().st_mtime) 
                if time_simulation>time_txtFile+timedelta(minutes=5):
                    simulationWasntFinished = True 
                
            # Create an h5 file that contains the moments  
            if not os.path.isfile(moments_file) or simulationWasntFinished:   
                
                # Check whether the moments are written
                netcdf_file  = read_outputFile(netcdf_file)   
                if "density" not in netcdf_file.variables.keys():
                    print(status+"The moments were not written to the netcdf file.")
                
                # Get the moments from the output file
                vec_time = read_netcdfVariables('vec_time', netcdf_file)  
                upar_vs_tszkxkyri = read_netcdfVariables('upar_vs_tszkxkyri', netcdf_file)
                dens_vs_tszkxkyri = read_netcdfVariables('dens_vs_tszkxkyri', netcdf_file)
                temp_vs_tszkxkyri = read_netcdfVariables('temp_vs_tszkxkyri', netcdf_file) 
                netcdf_file.close()
                
                # Get the data at every <dt> timestep 
                indices = get_indicesAtFixedStep(vec_time, dt)
                upar_vs_tszkxky = upar_vs_tszkxkyri[indices,:,:,:,:,0] + 1j*upar_vs_tszkxkyri[indices,:,:,:,:,1]
                dens_vs_tszkxky = dens_vs_tszkxkyri[indices,:,:,:,:,0] + 1j*dens_vs_tszkxkyri[indices,:,:,:,:,1]
                temp_vs_tszkxky = temp_vs_tszkxkyri[indices,:,:,:,:,0] + 1j*temp_vs_tszkxkyri[indices,:,:,:,:,1]
                vec_time = vec_time[indices]
                
                # Create the new h5 file 
                with h5py.File(moments_file, 'w') as h5_file:
                    h5_file.create_dataset("vec_time", data=vec_time) 
                    h5_file.create_dataset("upar_vs_tszkxky", data=upar_vs_tszkxky) 
                    h5_file.create_dataset("dens_vs_tszkxky", data=dens_vs_tszkxky) 
                    h5_file.create_dataset("temp_vs_tszkxky", data=temp_vs_tszkxky) 
        
                # Notify that we finished creating the file
                print("       ---> The 5D moments file is saved as " + moments_file.name) 
                       
            # Notify that the file already existed   
            else:
                print("    The 5D moments file already exists:", moments_file.parent.name+"/"+moments_file.name)
    return 
    
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")   
    write_h5FileForMoments5D(folder) 
    