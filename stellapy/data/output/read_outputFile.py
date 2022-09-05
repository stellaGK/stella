
#!/usr/bin/python3 
import copy
import os, sys
import numpy as np
from scipy.io import netcdf as scnetcdf     

# Stellapy package
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])  
from stellapy.data.stella.load_stellaVariables import stella_variables  
 
#===============================================================================
#                               READ THE OUT.NC FILE
#===============================================================================

def read_outputFile(netcdf_path):
    ''' Normally we read with the 'scnetcdf' package, however when the netcdf file
    is bugged/corrupted this will fail and we need to use the 'nc4' package instead. '''
    try: 
        netcdf_file = scnetcdf.netcdf_file(netcdf_path,'r') 
    except: 
        import netCDF4 as nc4 
        print("WARNING: THE NETCDF FILE IS CORRUPTED")
        print("    ", netcdf_path)
        netcdf_file = nc4.Dataset(netcdf_path)  
    return netcdf_file


#===============================================================================
#                 READ STELLA VARIABLES FROM THE NETCDF FILE
#===============================================================================

def read_netcdfVariables(key, netcdf_file, time_indices=None, netcdf_path=None): 
    ''' Read variables from the netcdf file based on their dimensions. '''

    # Derived variables that only exist in h5 files
    derived_variables = ['dimensions', 'g_vs_smu', 'g_vs_svpa', 'g_vs_ts', \
        'small', 'medium', 'large', 'phi_vs_tri', 'phi_vs_tkxkyri', 'phi_vs_tzri']
    
    # Check whether the key exists
    if key not in stella_variables:
        if key not in derived_variables:
            if "h5file" not in key and "reduced" not in key and 'final' not in key: 
                print("THE KEY ", key, " IS NOT A STELLAPY VARIABLE")
        return None

    # Get the dimensions and the stella key of the variable
    dim = stella_variables[key][1] 
    key = stella_variables[key][0]
    
    # Check whether the variable is in the netcdf file
    if key not in netcdf_file.variables.keys():
        if netcdf_path:
            print("WARNING:", key, "was not found in the netcdf file.")
            print(netcdf_path)
        return None
    
    # Make sure we have time indices
    if not time_indices: 
        time_indices = range(len(netcdf_file.variables["t"][:]))
        
    # Read values
    if len(dim)==1 and dim[0]=="-":
        variable = copy.deepcopy(netcdf_file.variables[key][:])
    
    # Read 1D vectors
    elif len(dim)==1 and dim[0]!="-":
        if dim[0]!="t": variable = copy.deepcopy(netcdf_file.variables[key][:])
        if dim[0]=="t": variable = copy.deepcopy(netcdf_file.variables[key][time_indices])
    
    # Read 2D
    elif len(dim)==2:
        if dim[0]!="t": variable = copy.deepcopy(netcdf_file.variables[key][:,:])
        if dim[0]=="t": variable = copy.deepcopy(netcdf_file.variables[key][time_indices,:])
    
    # Read 3D
    elif len(dim)==3:
        if dim[0]!="t": variable = copy.deepcopy(netcdf_file.variables[key][:,:,:])
        if dim[0]=="t": variable = copy.deepcopy(netcdf_file.variables[key][time_indices,:,:])
        
    # Read 4D
    elif len(dim)==4:
        if dim[0]!="t": variable = copy.deepcopy(netcdf_file.variables[key][:,:,:,:])
        if dim[0]=="t": variable = copy.deepcopy(netcdf_file.variables[key][time_indices,:,:,:])
    
    # Read 5D
    elif len(dim)==5:
        if dim[0]!="t": variable = copy.deepcopy(netcdf_file.variables[key][:,:,:,:,:])
        if dim[0]=="t": variable = copy.deepcopy(netcdf_file.variables[key][time_indices,:,:,:,:])
        
    # Read 6D
    elif len(dim)==6:
        if dim[0]!="t": variable = copy.deepcopy(netcdf_file.variables[key][:,:,:,:,:,:])
        if dim[0]=="t": variable = copy.deepcopy(netcdf_file.variables[key][time_indices,:,:,:,:,:])
            
    # Read 7D
    elif len(dim)==7:
        if dim[0]!="t": variable = copy.deepcopy(netcdf_file.variables[key][:,:,:,:,:,:,:])
        if dim[0]=="t": variable = copy.deepcopy(netcdf_file.variables[key][time_indices,:,:,:,:,:,:])
    
    # Remove the <tube> dimension since for now stella only simulates one tube
    if "tube" in dim and len(dim)>1:
        index = dim.index("tube")
        while index!=0: variable = np.swapaxes(variable, index-1, index); index -= 1
        if len(dim)==2: variable = variable[0,:]
        if len(dim)==3: variable = variable[0,:,:]
        if len(dim)==4: variable = variable[0,:,:,:]
        if len(dim)==5: variable = variable[0,:,:,:,:]
        if len(dim)==6: variable = variable[0,:,:,:,:,:] 
        if len(dim)==7: variable = variable[0,:,:,:,:,:,:]  
    
    # Return the variable 
    return variable
