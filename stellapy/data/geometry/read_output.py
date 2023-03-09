
import h5py
import os, sys
import numpy as np
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.data.output.read_outputFile import read_netcdfVariables
from stellapy.data.output.read_outputFile import read_outputFile as read_outputFileFromNetcdf
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_muWeights, calculate_vpaWeights
from stellapy.data.dimensions.read_dimensionsAndVectors import read_dimensions
 
#===============================================================================
#                    READ THE *.OUT.NC OR *.OUT.H5 FILE
#===============================================================================
 
def read_outputFile(path): 
    
    # Read from the *rho70.geometry file
    if os.path.isfile(path.geometry): 
        return read_outFromGeometryFile(path.geometry)
         
    # Read from the *.out.nc file
    elif os.path.isfile(path.output_stella): 
        return read_outFromNcFile(path.output_stella, path.input_file)
    
    # Read from the *.out.h5 file
    elif os.path.isfile(path.output): 
        return read_outFromH5File(path.output, path.input_file)
    
    # Critical error if we didn't find any data
    exit_reason = "The output data can not be found to read the geometry data.\n"
    exit_reason += "For linear simulations with 1 mode per simulation, make sure \n"
    exit_reason += "the geometry data for the dummy input is written by doing: \n" 
    exit_reason += ">>> write_dataFiles"
    exit_program(exit_reason, read_outputFile, sys._getframe().f_lineno)    

    # Return nans if we didn't find any data 
    dim = read_dimensions(path); z = dim["dim_z"]; mu = dim["dim_mu"]; vpa = dim["dim_vpa"] 
    return {"jacob" : np.ones((z))*np.nan, "bmag" : np.ones((z))*np.nan,\
            "dl_over_B" : np.ones((z))*np.nan, "mu_weights" : np.ones((z,mu))*np.nan,\
            "vpa_weights" : np.ones((vpa))*np.nan, "dim_z" : z, "vec_z" : np.ones((z))*np.nan}
    return 

#-------------------------------------
def read_outFromGeometryFile(path):
    
    # Recall what file we are reading from
    out_data = {"source" : path.name}
  
    # Read the output data from the *.rho70.geometry file
    with h5py.File(path, 'r') as f:  
        for key in ["jacob", "bmag", "dim_z", "vec_z", "dl_over_B", "mu_weights", "vpa_weights"]:
            if key in f.keys():
                out_data[key] = f[key][()]  
                
    # Return the geometry data
    return out_data

#-------------------------------------
def read_outFromH5File(path, input_file):
     
    # Recall what file we are reading from
    out_data = {"source" : path.name}
     
    # Read the Jacobian and the magnetic field strength from the output file
    with h5py.File(path, 'r') as h5_file: 
        if 'jacob' in h5_file.keys(): out_data['jacob'] = h5_file['jacob'][:]
        if 'bmag'  in h5_file.keys(): out_data['bmag']  = h5_file['bmag'][:]
        if 'vec_z' in h5_file.keys(): out_data['vec_z'] = h5_file['vec_z'][:]
        if 'vec_z' in h5_file.keys(): out_data['dim_z'] = len(out_data['vec_z'])
        
    # Calculate the integration weights <dl_over_B>
    out_data['dl_over_B'] = calculate_dlOverB(out_data['dim_z'], out_data['vec_z'], out_data['jacob'])
    out_data['mu_weights'] = calculate_muWeights(input_file, bmag_psi0=out_data['bmag'])[0,:,:] 
    out_data['vpa_weights'] = calculate_vpaWeights(input_file) 
    
    # Return the geometry data that is read from the output file
    return out_data
        
#-------------------------------------      
def read_outFromNcFile(path, input_file):
    
    # Open the "out.nc" file   
    netcdf_file = read_outputFileFromNetcdf(path)  
    
    # Recall what file we are reading from
    out_data = {"source" : path.name}
    
    # Read the stella variables
    out_data['jacob'] = read_netcdfVariables('jacob', netcdf_file)  # jacob(zed,alpha)
    out_data['bmag']  = read_netcdfVariables('bmag', netcdf_file)   # bmag(zed,alpha)
    out_data['vec_z'] = read_netcdfVariables('vec_z', netcdf_file)  # zed(zed)
    out_data['dim_z'] = len(out_data['vec_z'])
    
    # Close the netcdf file
    netcdf_file.close()
    
    # Calculate the integration weights <dl_over_B>, <vpa_weights> and <mu_weights>
    out_data['dl_over_B'] = calculate_dlOverB(out_data['dim_z'], out_data['vec_z'], out_data['jacob'])
    out_data['mu_weights'] = calculate_muWeights(input_file, bmag_psi0=out_data['bmag'])[0,:,:] 
    out_data['vpa_weights'] = calculate_vpaWeights(input_file) 
    
    # Return the geometry data that is read from the output file
    return out_data


#-------------------------------------
def calculate_dlOverB(dim_z, vec_z, jacobian):
    ''' Get the integration weights for z. ''' 

    # Calculate the step in <z>
    delzed = np.zeros((int(dim_z))) 
    for i in range(dim_z-1): delzed[i] = vec_z[i+1]-vec_z[i] 
    delzed[-1] = delzed[0] 
    
    # The integration weights are dz*J = dz/B
    dl_over_B = delzed*jacobian[:,0] 
    
    # Avoid double counting the end points for ky = 0 modes
    dl_over_B[-1] = 0
    
    # Normalize the dl/B factor by int dl/B = sum_z dl/B
    dl_over_B = dl_over_B / np.nansum(dl_over_B) 
    
    # Return the integration weights
    return dl_over_B

 
#===============================================================================
#               ATTACH THE OUTPUT DATA TO THE SIMULATION OBJECT
#===============================================================================

def get_outputData(self):
    """ Attach the data inside the output file to the simulation object. """

    # Read the output file
    out_data = read_outputFile(self.path)

    # Attach the data   
    self.bmag = out_data['bmag']
    self.jacob = out_data['jacob']
    self.vec_z = out_data['vec_z']
    self.dl_over_B = out_data['dl_over_B']
    self.mu_weights = out_data['mu_weights']
    self.vpa_weights = out_data['vpa_weights']
    return


