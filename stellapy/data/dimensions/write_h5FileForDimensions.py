  
import os, h5py, pathlib
from datetime import datetime, timedelta
from stellapy.utils.decorators.verbose import noverbose 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder     
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables

################################################################################
#                       WRITE THE DIMENSIONS TO AN H5 FILE
################################################################################

@noverbose
def write_h5FileForDimensions(folder, dimensions={}):       
    
    # Get the input files
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix('.out.nc'))]
    if input_files==[]: return 
 
    # Iterate through the input files   
    for input_file in input_files: 
        
        # Processing status
        status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
            
        # Path of the new file 
        output_file = input_file.with_suffix(".out") if os.path.isfile(input_file.with_suffix(".out")) else input_file.with_suffix(".out.nc")
        netcdf_file = input_file.with_suffix(".out.nc")
        dimensions_file = input_file.with_suffix(".dimensions")   
        
        # Check whether the txt file is older than the simulation
        simulationWasntFinished = False
        if os.path.isfile(dimensions_file): 
            time_txtFile = datetime.fromtimestamp(dimensions_file.stat().st_mtime)
            time_simulation = datetime.fromtimestamp(output_file.stat().st_mtime) 
            if time_simulation>time_txtFile+timedelta(minutes=5):
                simulationWasntFinished = True 
        
        # Create an h5 file that contains the dimensions (t,s,z,mu,vpa,kx,ky)
        if not os.path.isfile(dimensions_file) or simulationWasntFinished:   
            
            # Get the dimensions from the output file
            netcdf_file = read_outputFile(netcdf_file)  
            dimensions['vec_z']   = read_netcdfVariables('vec_z',   netcdf_file)
            dimensions['vec_kx']  = read_netcdfVariables('vec_kx',  netcdf_file)
            dimensions['vec_ky']  = read_netcdfVariables('vec_ky',  netcdf_file)
            dimensions['vec_mu']  = read_netcdfVariables('vec_mu',  netcdf_file)
            dimensions['vec_vpa'] = read_netcdfVariables('vec_vpa', netcdf_file)
            dimensions['vec_time'] = read_netcdfVariables('vec_time', netcdf_file)
            dimensions['vec_species'] = read_netcdfVariables('vec_species',  netcdf_file) 
            netcdf_file.close()
            
            # Get the size of the dimensions
            dimensions['dim_z']   = len(dimensions['vec_z'])
            dimensions['dim_kx']  = len(dimensions['vec_kx'])
            dimensions['dim_ky']  = len(dimensions['vec_ky'])
            dimensions['dim_mu']  = len(dimensions['vec_mu'])
            dimensions['dim_vpa'] = len(dimensions['vec_vpa'])
            dimensions['dim_time'] = len(dimensions['vec_time'])
            dimensions['dim_species'] = len(dimensions['vec_species'])
            
            # Create the new h5 file 
            with h5py.File(dimensions_file, 'w') as h5_file:
                for key, value in dimensions.items():    
                    h5_file.create_dataset(key, data=value) 

            # Notify that we finished creating the file
            print(status+"   ---> The dimensions file is saved as " + dimensions_file.name)    
    return 
    
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")   
    write_h5FileForDimensions(folder) 
    