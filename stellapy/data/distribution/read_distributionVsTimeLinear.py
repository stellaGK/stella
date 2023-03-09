""" 

#===============================================================================
#                        Read the 2D distribution data                         #
#===============================================================================
 
Read the distribution function "g" from the *.g2_vs_t or *.out.nc files.

Returns
-------
    {vec_time, g2_vs_tskxky}

Hanne Thienpondt
27/01/2023

"""
 
#!/usr/bin/python3 
import numpy as np   
import os, sys, pathlib 

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.data.distribution.write_txtFileForDistributionVsTime import write_txtFileForDistributionVsTime  
from stellapy.data.distribution.read_distributionVsMuOrVpaOrZ import print_warningForRangedSimulations
from stellapy.data.input.read_inputFile import read_modeFromInputFile
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.data.utils import Data

#===============================================================================
#                        Read the 2D distribution data                         #
#===============================================================================  

def read_distributionVsTimeLinear(path, dim, vec, nakxnaky):
    
    # Make sure the *.g2_vs_t file exists
    if not os.path.isfile(path.g2_vs_t) and os.path.isfile(path.output_stella): 
        write_txtFileForDistributionVsTime(path.folder)
        
    # Read <range> simulations
    if os.path.isfile(path.g2_vs_t): 
        if nakxnaky>1:  return read_distributionFromH5File(path, dim) 
        if nakxnaky==1: return read_distributionFromTxtFiles(dim, vec, path.dummy_paths) 

    # Read one-mode-per-file simulations
    if path.dummy_input_file:
        return read_distributionFromTxtFiles(dim, vec, path.dummy_paths) 
    
    # Critical error if we didn't find any data
    exit_reason = "The distribution versus time can not be found.\n"
    exit_reason += "      "+str(path.folder)+"\n"
    exit_reason += "      "+str(path.g2_vs_t)
    exit_program(exit_reason, read_distributionVsTimeLinear, sys._getframe().f_lineno) 
    return 

#-------------------------------------
def read_distributionFromH5File(path, dim): 
    """ For <range> simulations we do not have g(kx,ky) but only a single g() for all modes. """
    print_warningForRangedSimulations(path) 
    dummy_tskxky = np.ones((10, dim.species, dim.kx, dim.ky))*np.nan 
    dummy_tkxky = np.ones((10, dim.kx, dim.ky))*np.nan 
    return {"g2_vs_tskxky" : dummy_tskxky, "time_vs_tkxky" : dummy_tkxky} 

#-----------------------------
def read_distributionFromTxtFiles(dim, vec, paths): 
    
    # Create matrices (kx,ky) 
    g2_vs_tskxky = np.ones((0, dim.species, dim.kx, dim.ky))*np.nan 
    time_vs_tkxky = np.ones((0, dim.kx, dim.ky))*np.nan  
    
    # Iterate over the simulations (1 mode per simulation)
    for path in paths:
        
        # Get the mode (kx,ky)
        kx, ky = read_modeFromInputFile(path.input_file)
        ikx = list(vec.kx).index(kx)
        iky = list(vec.ky).index(ky)
        
        # Get the potential data 
        distribution_data = read_distributionFromTxtFile(path, dim) 
        
        # Make sure we have enough time points 
        dim_time = len(distribution_data["vec_time"][:]); dim_time_old = np.shape(time_vs_tkxky)[0]
        if dim_time_old < dim_time:
            time_vs_tkxky = np.append(time_vs_tkxky, np.ones((dim_time-dim_time_old, dim.kx, dim.ky))*np.nan, axis=0)
            g2_vs_tskxky  = np.append(g2_vs_tskxky, np.ones((dim_time-dim_time_old, dim.specie, dim.kx, dim.ky))*np.nan, axis=0) 
            
        # Put the distribution data in the matrices    
        dim_time = len(distribution_data["vec_time"][:]) 
        time_vs_tkxky[:dim_time,ikx,iky] = distribution_data["vec_time"]  
        g2_vs_tskxky[:dim_time,dim.species,ikx,iky] = distribution_data["g2_vs_ts"]  
 
    return {"g2_vs_tskxky" : g2_vs_tskxky, "time_vs_tkxky" : time_vs_tkxky}

#-----------------------------
def read_distributionFromTxtFile(path, dim): 
        
    # Read the *.g2_vs_t file 
    if os.path.isfile(path.g2_vs_t): 
        data = np.loadtxt(path.g2_vs_t,skiprows=1, dtype='float').reshape(-1, dim.species+1)
        g2_vs_ts = data[:,1:]; vec_time = data[:,0]
        return {"g2_vs_ts" : g2_vs_ts, "vec_time" : vec_time}

#===============================================================================
#           ATTACH THE DISTRIBUTION SQUARED TO THE SIMULATION OBJECT           #
#===============================================================================

def get_distributionVsTimeLinear(self): 
    
    # Read the data in the fluxes file
    distribution = read_distributionVsTimeLinear(self.path, self.dim, self.vec, self.nakxnaky)  
    
    # Save the distribution data 
    self.g2_vs_tskxky = Data(["g2","t","s","kx","ky"], distribution["g2_vs_tskxky"], distribution["time_vs_tkxky"], range(self.dim.species), self.vec.kx, self.vec.ky)
    return 