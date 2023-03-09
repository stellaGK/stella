
#!/usr/bin/python3 
import copy
import os, sys
import pathlib
import numpy as np
import netCDF4 as nc4  
from scipy.io import netcdf as scnetcdf     

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.stella.load_stellaVariables import stella_variables 
from stellapy.utils.decorators.exit_program import exit_program
 
#===============================================================================
#                               READ THE OUT.NC FILE
#===============================================================================

def read_outputFile(netcdf_path):
    ''' Since a stella update around summer 2022, there are two unlimited variables  
    in the netcdf file and <scnetcdf> will no longer work.'''
    try:
        try: netcdf_file = nc4.Dataset(netcdf_path)  
        except: netcdf_file = scnetcdf.netcdf_file(netcdf_path,'r') 
    except:
        if not os.path.isfile(netcdf_path):
            exit_reason = "The netcdf file could not be found:\n"
            exit_reason += "     "+str(netcdf_path)+"\n"
            exit_program(exit_reason, read_outputFile, sys._getframe().f_lineno)
        if os.path.isfile(netcdf_path):
            exit_reason = "The netcdf file could not be read with nc4 or scnetcdf,\n"
            exit_reason += "there must be something wrong with it, please delete it and rerun:\n"
            exit_reason += "     "+str(netcdf_path)+"\n"
            exit_program(exit_reason, read_outputFile, sys._getframe().f_lineno)
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
        try: time_indices = range(len(netcdf_file.variables["t"][:]))
        except:
            exit_reason = "No time dimension was found in the following netcdf file,\n"
            exit_reason += "something went wrong with the simulation, relaunch it:\n"
            exit_reason += "     "+str(netcdf_file.filepath())+"\n"
            exit_program(exit_reason, read_netcdfVariables, sys._getframe().f_lineno)
    
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

 
