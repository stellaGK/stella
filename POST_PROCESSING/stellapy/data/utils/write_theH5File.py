
import os, h5py
def write_theH5File(file_path, netcdf_file, data_keys):

    # Check whether the reduced file has already been created
    already_written = False
    if os.path.isfile(file_path):
        
        # Check whether the reduced file is older than the netcdf file
        time_h5File = file_path.stat().st_mtime
        time_netcdfFile = netcdf_file.stat().st_mtime
        if time_h5File>time_netcdfFile:
            already_written = True
         
        # Check whether the existing file contains the correct data listed above
        h5_netcdf = h5py.File(file_path, 'r')  
        for key in data_keys: 
            if key not in h5_netcdf and key!="dimensions":    
                already_written = False 
        h5_netcdf.close()
        
    # If the file is missing or the data inside needs to change, then write the h5 file 
    return (not os.path.isfile(file_path) or (already_written==False))
