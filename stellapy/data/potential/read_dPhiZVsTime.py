
import os, sys
import numpy as np   
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.data.potential.write_txtFileForDPhiZVsTime import write_txtFileForDPhiZVsTime

#===============================================================================
#                    READ THE DPHI(Z) DATA FROM THE TXT FILE                   #
#===============================================================================    

def read_dPhiZVsTime(path, potential_data={}):
    
    # Make sure the *.dphiz_vs_t file exists
    if not os.path.isfile(path.dphiz_vs_t) and os.path.isfile(path.output_stella): 
        write_txtFileForDPhiZVsTime(path.folder)

    # Read the *.dphiz_vs_t file
    if os.path.isfile(path.dphiz_vs_t): 
        data = np.loadtxt(path.dphiz_vs_t,skiprows=1, dtype='float').reshape(-1, 2)
        potential_data["vec_time"] = data[:,0]
        potential_data["dphiz_vs_t"] = data[:,1]
        return potential_data
    
    # Otherwise return nan
    potential_data["vec_time"] = [np.nan]
    potential_data["dphiz_vs_t"] = [np.nan]
    return potential_data

    # Critical error if we didn't find any data
    exit_program("The dphi(z) data versus time can not be found.", read_dPhiZVsTime, sys._getframe().f_lineno) 
    return 

#===============================================================================
#                ATTACH THE DPHI(Z) DATA TO THE SIMULATION OBJECT              #
#===============================================================================

def get_dPhiZVsTime(self): 
    
    # Read the data in the fluxes file
    potential = read_dPhiZVsTime(self.path)  
    
    # Save the distribution data  
    self.dphiz_vs_t = Data(["dphiz", "t"], potential["dphiz_vs_t"], potential["vec_time"])
    return 