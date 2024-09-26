

from stellapy.data.utils.save_memory import save_memory
def create_datasets(h5_netcdf, netcdf_data, h5_keys, data_keys, memory):
    
    # Check whether the h5 keys were read from the netcdf file 
    if any(item in h5_keys for item in netcdf_data):
        
        # Iterate over the h5 keys to save them to the <h5_netcdf> file
        for key in h5_keys: 
            
            # Only save the quantity if it was read from the netcdf file
            if (key in netcdf_data) and (key in data_keys):  
                
                # Save the quantity to the h5 file   
                h5_netcdf.create_dataset(key, data=netcdf_data[key])  
                
                # If we are doing a memory check, save the name and size of the quantity
                if memory['check']==True: memory = save_memory(memory, key, netcdf_data[key]) 
    return 