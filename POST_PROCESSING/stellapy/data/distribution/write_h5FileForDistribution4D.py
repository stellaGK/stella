  
import os, h5py, pathlib
from datetime import datetime, timedelta
from stellapy.utils.decorators.verbose import noverbose 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder   
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables 
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep

################################################################################
#                     WRITE THE DISTRIBUTION TO AN H5 FILE                     #
################################################################################

@noverbose
def write_h5FileForDistribution4D(folder, dt=10, simulationWasntFinished=False):      
    
    # Check whether we have any simulations
    if not get_filesInFolder(folder, end=".out.nc"):
        return
 
    # Iterate through the input files  
    for input_file in get_filesInFolder(folder, end=".in"):  
            
        # Path of the new file
        output_file = input_file.with_suffix(".out")
        netcdf_file = input_file.with_suffix(".out.nc")
        distribution_file = input_file.with_suffix(".dt"+str(dt)+".distribution4D")   
    
        # If the file doesn't exist, check whether gvmus and gzvs were written
        if not os.path.isfile(distribution_file): 
            netcdf_data  = read_outputFile(netcdf_file)  
            if "gvmus" not in netcdf_data.variables.keys() or "gzvs" not in netcdf_data.variables.keys():
                continue
    
        # Check whether the txt file is older than the simulation
        if os.path.isfile(distribution_file): 
            time_txtFile = datetime.fromtimestamp(distribution_file.stat().st_mtime)
            time_simulation = datetime.fromtimestamp(output_file.stat().st_mtime) 
            if time_simulation>time_txtFile+timedelta(minutes=5):
                simulationWasntFinished = True 
            
        # Create an h5 file that contains the distribution  
        if not os.path.isfile(distribution_file) or simulationWasntFinished:    
            
            # Get the distribution data from the output file  
            netcdf_data  = read_outputFile(netcdf_file)  
            g_vs_tsmuvpa = read_netcdfVariables('g_vs_tsmuvpa', netcdf_data) 
            g_vs_tsvpaz  = read_netcdfVariables('g_vs_tsvpaz', netcdf_data) 
            vec_time     = read_netcdfVariables('vec_time', netcdf_data) 
            netcdf_data.close()
            
            # Get the data at every <dt> timestep
            indices = get_indicesAtFixedStep(vec_time, dt)
            g_vs_tsmuvpa = g_vs_tsmuvpa[indices] 
            g_vs_tsvpaz = g_vs_tsvpaz[indices] 
            vec_time = vec_time[indices]
                
            # Create the new h5 file 
            with h5py.File(distribution_file, 'w') as h5_file: 
                h5_file.create_dataset("vec_time", data=vec_time) 
                h5_file.create_dataset("g_vs_tsvpaz", data=g_vs_tsvpaz) 
                h5_file.create_dataset("g_vs_tsmuvpa", data=g_vs_tsmuvpa)  
    
            # Notify that we finished creating the file
            print("       ---> The distribution file is saved as " + distribution_file.name) 
                   
        # Notify that the file already existed   
        else:
            print("    The distribution file already exists:", distribution_file.parent.name+"/"+distribution_file.name)
    return 

    