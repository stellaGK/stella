
import os, sys
import numpy as np   
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.data.potential.write_txtFileForPotentialVsZ import write_txtFileForPotentialVsZ
from stellapy.data.dimensions.get_dimensionsAndVectors import read_dimensions

#===============================================================================
#                 READ THE 1D POTENTIAL DATA FROM THE TXT FILE
#===============================================================================    

def read_potentialVsZ(path, potential_data={}):
    
    # Make sure the *.finalphiz file exists
    if not os.path.isfile(path.phi_vs_z) and os.path.isfile(path.output_stella): 
        write_txtFileForPotentialVsZ(path.folder)

    # Read the *.finalphiz file  
    if os.path.isfile(path.phi_vs_z): 
        data = np.loadtxt(path.phi_vs_z,skiprows=1,dtype='float').reshape(-1, 3) 
        potential_data["phi2_vs_z"] = data[:,0]  
        potential_data["phi_vs_z"] = data[:,1] + 1j*data[:,2]
        return potential_data 

    # Read the *.final_fields file  
    if os.path.isfile(path.finalphi_stella): 
        try:    data = np.loadtxt(path.finalphi_stella, dtype='float').reshape(-1, 10) # Old code has 10 columns
        except: data = np.loadtxt(path.finalphi_stella, dtype='float').reshape(-1, 11) # New code has 11 columns
        potential_data["phi_vs_z"] = data[:,4] + 1j*data[:,5]
        potential_data["phi2_vs_z"] = np.abs(potential_data["phi_vs_z"])**2
        return potential_data
    
    # Return nans if we didn't find any data
    print("The potential(z) can not be found for:")
    print("      "+str(path.input_file)) 
    dim = read_dimensions(path); z = dim["dim_z"] 
    return {"phi_vs_z" : np.ones((z))*np.nan, "phi2_vs_z" : np.ones((z))*np.nan}

    # Critical error if we didn't find any data
    exit_program("The 1D potential data can not be found.", read_potentialVsZ, sys._getframe().f_lineno) 
    return 

#===============================================================================
#             ATTACH THE POTENTIAL DATA TO THE SIMULATION OBJECT            #
#===============================================================================

def get_potentialVsZ(self): 
    
    # Read the data in the fluxes file
    potential = read_potentialVsZ(self.path)  
    
    # Sometimes phi is still good but phi2 is too big
    if np.any(~np.isfinite(potential["phi2_vs_z"])):
        potential["phi2_vs_z"] = np.abs(potential["phi_vs_z"]/1e300)**2

    # Save the potential data  
    self.phi_vs_z = Data(["phi", "z"], potential["phi_vs_z"], self.vec.z)
    self.phi2_vs_z = Data(["phi2", "z"], potential["phi2_vs_z"], self.vec.z)
    return 