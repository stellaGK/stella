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

def read_distributionVsTime(path, dim):
    
    # Make sure the *.g2_vs_t file exists
    if not os.path.isfile(path.g2_vs_t) and os.path.isfile(path.output_stella): 
        write_txtFileForDistributionVsTime(path.folder)
        
    # Read the text file
    if os.path.isfile(path.g2_vs_t): 
        return read_distributionFromTxtFile(path, dim) 
    elif os.path.isfile(path.g_vs_t): 
        path.g2_vs_t = path.g_vs_t
        return read_distributionFromTxtFile(path, dim) 
    
    # Critical error if we didn't find any data
    exit_reason = "The distribution versus time can not be found.\n"
    exit_reason += "      "+str(path.folder)+"\n"
    exit_reason += "      "+str(path.g2_vs_t)
    exit_program(exit_reason, read_distributionVsTime, sys._getframe().f_lineno) 
    return 

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

def get_distributionVsTime(self): 
    
    # Read the data in the fluxes file
    distribution = read_distributionVsTime(self.path, self.dim)  
    
    # Save the distribution data 
    self.g2_vs_ts = Data(["g2","t","s"], distribution["g2_vs_ts"], distribution["vec_time"], range(self.dim.species))
    return 