
import numpy as np
import os, sys, h5py  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program  
from stellapy.data.paths.load_pathObject import get_distribution3DPath
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry
from stellapy.data.distribution.write_h5FileForDistribution3D import write_h5FileForDistribution3D

#===============================================================================
#                        ATTACH THE DISTRIBUTION DATA
#===============================================================================

def read_distribution3D(path): 
    ''' Read the distribution from *.distribution3D or *.out.nc ''' 
               
    # Make sure that the *.dt10.distribution file exists
    if not os.path.isfile(path.distribution3D) and os.path.isfile(path.output_stella): 
        write_h5FileForDistribution3D(path.folder)
        get_distribution3DPath(path)  
                              
    # Read from the *.dt10.distribution file 
    if os.path.isfile(path.distribution3D):
        return read_fromdistributionFile(path.distribution3D)    
    
    # This file doesn't exist for old simulations
    else: 
        distribution_data = read_fromH5File(path)
        return distribution_data 
    
    # Critical error if we didn't find any data
    exit_reason = "The distribution data could not be found."
    exit_program(exit_reason, read_distribution3D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromdistributionFile(path, distribution={}):   
    with h5py.File(path, 'r') as f: 
        for key in ["vec_time", "g_vs_tsz", "g_vs_tsmu", "g_vs_tsvpa"]: 
            if key in f.keys():   
                distribution[key] = f[key][()]  
    return distribution 
    
#===============================================================================
#                            READ FROM OLD H5 FILE                             #
#===============================================================================

def read_fromH5File(path, distribution_data={}):
            
    # Read the geometry data in the output file
    vmec_filename = read_vmecFileNameFromInputFile(path.input_file)
    nonlinear = read_linearNonlinearFromInputFile(path.input_file)[1]
    path = create_dummyPathObject(path.input_file, vmec_filename, nonlinear)
    geometry = read_outputFileForGeometry(path) 
    vpa_weights = geometry["vpa_weights"] 
    mu_weights = geometry["mu_weights"] 
    dl_over_B = geometry["dl_over_B"] 
    
    # Read the data   
    with h5py.File(path.output, 'r') as f:    
        if "vec_time_reduced" in f.keys():  
            distribution_data["vec_time"] = f["vec_time_reduced"][()] 
            g_vs_tsmuvpa = f["g_vs_tsmuvpa_reduced"][()] 
            g_vs_tsvpaz = f["g_vs_tsvpaz_reduced"][()] 
        elif "vec_time_reduced" not in f.keys():  
            distribution_data["vec_time"] = f["vec_time"][()]  
            step_size = int(len(distribution_data["vec_time"])/100)
            if step_size==0: step_size = 1
            distribution_data["vec_time"] = distribution_data["vec_time"][::step_size]
            g_vs_tsmuvpa = f["gvmus"][()] 
            g_vs_tsvpaz = f["gzvs"][()][:,:,:,:,0]
    
    # Get the dimensions  
    dim_time, dim_species, dim_mu, dim_vpa = np.shape(g_vs_tsmuvpa)  
    dim_time, dim_species, dim_vpa, dim_z = np.shape(g_vs_tsvpaz) 
     
    # Calculate the 3D distribution functions
    distribution_data["g_vs_tsz"] = np.zeros((dim_time, dim_species, dim_z)) 
    distribution_data["g_vs_tsmu"] = np.zeros((dim_time, dim_species, dim_mu)) 
    distribution_data["g_vs_tsvpa"] = np.zeros((dim_time, dim_species, dim_vpa))  
    
    # Sum over (vpa,mu,z)
    for vpa in range(dim_vpa): 
        for mu in range(dim_mu):
            for z in range(dim_z): 
                product = g_vs_tsmuvpa[:,:,mu,vpa]*mu_weights[z,mu]*vpa_weights[vpa]*dl_over_B[z] 
                distribution_data["g_vs_tsz"][:,:,z] += product  
                distribution_data["g_vs_tsmu"][:,:,mu] += product  
                distribution_data["g_vs_tsvpa"][:,:,vpa] += product  
                
    return distribution_data
       
#===============================================================================
#                  ATTACH THE distribution TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_distribution3D(self): 
    
    # Read the distribution
    distribution = read_distribution3D(self.path)  
    
    # Save the distribution
    self.vec_time_3D = distribution["vec_time"]
    self.g_vs_tsz = Data(["g","t","z"], distribution["g_vs_tsz"], distribution["vec_time"], self.vec.z) 
    self.g_vs_tsmu = Data(["g","t","mu"], distribution["g_vs_tsmu"], self.g_vs_tsz.t, self.vec.mu) 
    self.g_vs_tsvpa = Data(["g","t","vpa"], distribution["g_vs_tsvpa"], self.g_vs_tsz.t, self.vec.vpa)  
    return 
