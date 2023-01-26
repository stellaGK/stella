
import h5py
import os, sys
import numpy as np   
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.data.paths.load_pathObject import create_dummyPathObject
from stellapy.data.input.read_inputFile import read_vmecFileNameFromInputFile
from stellapy.data.input.read_inputFile import read_linearNonlinearFromInputFile
from stellapy.data.geometry.read_output import read_outputFile as read_outputFileForGeometry
from stellapy.data.distribution.write_txtFileForDistributionVsTime import write_txtFileForDistributionVsTime
from stellapy.data.utils.get_indicesAtFixedStep import get_indicesAtFixedStep

#===============================================================================
#                 READ THE 1D DISTRIBUTION DATA FROM THE TXT FILE
#===============================================================================    

def read_distributionVsTime(path, dim_species, distribution_data={}):
    
    # Make sure the *.g_vs_t file exists
    if not os.path.isfile(path.g_vs_t) and os.path.isfile(path.output_stella): 
        write_txtFileForDistributionVsTime(path.folder)

    # Read the *.g_vs_t file
    if os.path.isfile(path.g_vs_t): 
        data = np.loadtxt(path.g_vs_t,skiprows=1, dtype='float').reshape(-1, dim_species+1)
        distribution_data["g_vs_ts"] = data[:,1:]
        distribution_data["vec_time"] = data[:,0]
        return distribution_data
    
    # This file doesn't exist for old simulations
    else: 
        distribution_data = read_fromH5File(path)
        return distribution_data 
    
    # Critical error if we didn't find any data
    exit_reason = "The distribution versus time can not be found.\n"
    exit_reason += "      "+str(path.folder)+"\n"
    exit_reason += "      "+str(path.g_vs_t)
    exit_program(exit_reason, read_distributionVsTime, sys._getframe().f_lineno) 
    return 

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
        elif "vec_time_reduced" not in f.keys(): 
            distribution_data["vec_time"] = f["vec_time"][()]  
            step_size = int(len(distribution_data["vec_time"])/100)
            if step_size==0: step_size = 1
            distribution_data["vec_time"] = distribution_data["vec_time"][::step_size]
            g_vs_tsmuvpa = f["gvmus"][()] 
    
    # Get the (t, species, mu, vpa) dimensions  
    dim_time, dim_species, dim_mu, dim_vpa = np.shape(g_vs_tsmuvpa) 
    
    # Calculate the final g(vpa) and g(mu)
    distribution_data["g_vs_ts"] = np.zeros((dim_time, dim_species))  
    
    # Sum over (vpa,mu,z) in g_vs_smuvpa(s,mu,vpa) 
    for vpa in range(dim_vpa): 
        for mu in range(dim_mu):
            for z in range(len(dl_over_B)): 
                product = g_vs_tsmuvpa[:,:,mu,vpa]*mu_weights[z,mu]*vpa_weights[vpa]*dl_over_B[z] 
                distribution_data["g_vs_ts"][:,:] += product 
                
    return distribution_data

#===============================================================================
#             ATTACH THE DISTRIBUTION DATA TO THE SIMULATION OBJECT            #
#===============================================================================

def get_distributionVsTime(self): 
    
    # Read the data in the fluxes file
    distribution = read_distributionVsTime(self.path, self.dim.species)  
    
    # Save the distribution data 
    self.g_vs_ts = Data(["g","t","s"], distribution["g_vs_ts"], distribution["vec_time"], range(self.dim.species))
    return 