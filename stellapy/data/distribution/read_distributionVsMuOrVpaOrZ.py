
import os, sys
import numpy as np   
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.data.distribution.write_txtFileForDistributionVsMuOrVpaOrZ import write_txtFileForDistributionVsMuOrVpaOrZ
from stellapy.data.dimensions.get_dimensionsAndVectors import read_dimensions

#===============================================================================
#                 READ THE 1D DISTRIBUTION DATA FROM THE TXT FILE
#===============================================================================    

def read_distributionVsMuOrVpaOrZ(path, dim_species, distribution_data={}):
    
    # Make sure the *.finalgmu and *.finalgvpa files exist
    if not os.path.isfile(path.g_vs_z) and os.path.isfile(path.output_stella): 
        write_txtFileForDistributionVsMuOrVpaOrZ(path.folder)

    # Read the *.finalgmu and *.finalgvpa file
    if os.path.isfile(path.g_vs_z): 
        distribution_data["g_vs_zs"] = np.loadtxt(path.g_vs_z,skiprows=1,dtype='float').reshape(-1, dim_species)
        distribution_data["g_vs_mus"] = np.loadtxt(path.g_vs_mu,skiprows=1,dtype='float').reshape(-1, dim_species)
        distribution_data["g_vs_vpas"] = np.loadtxt(path.g_vs_vpa,skiprows=1,dtype='float').reshape(-1, dim_species)   
        return distribution_data

    # Return nans if we didn't find any data
    print("The 1D distribution data can not be found for:")
    print("      "+str(path.input_file))
    dim = read_dimensions(path); z = dim["dim_z"]; mu = dim["dim_mu"]; vpa = dim["dim_vpa"]; s = dim["dim_species"];  
    return {"g_vs_zs" : np.ones((z,s))*np.nan, "g_vs_mus" : np.ones((mu,s))*np.nan, "g_vs_vpas" : np.ones((vpa,s))*np.nan}

    # Critical error if we didn't find any data
    exit_program("The 1D distribution data can not be found.", read_distributionVsMuOrVpaOrZ, sys._getframe().f_lineno) 
    return 

#===============================================================================
#             ATTACH THE DISTRIBUTION DATA TO THE SIMULATION OBJECT            #
#===============================================================================

def get_distributionVsMuOrVpaOrZ(self): 

    # Read the data in the fluxes file
    distribution = read_distributionVsMuOrVpaOrZ(self.path, self.dim.species)  
    
    # Save the distribution data 
    self.g_vs_sz = Data(["g","s","z"], np.swapaxes(distribution["g_vs_zs"],0,1), range(self.dim.species), self.vec.z)
    self.g_vs_smu = Data(["g","s","mu"], np.swapaxes(distribution["g_vs_mus"],0,1), range(self.dim.species), self.vec.mu)
    self.g_vs_svpa = Data(["g","s","vpa"], np.swapaxes(distribution["g_vs_vpas"],0,1), range(self.dim.species), self.vec.vpa)
    return 