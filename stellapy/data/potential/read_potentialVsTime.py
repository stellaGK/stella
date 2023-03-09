
#!/usr/bin/python3  
import os, sys
import pathlib
import numpy as np  

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.data.potential.write_txtFileForPotentialVsTime import write_txtFileForPotentialVsTime  
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.data.utils import Data 

#===============================================================================
#                              READ POTENTIAL(TIME)                            #
#===============================================================================    

def read_potentialVsTime(path):  
      
    # Make sure the *.phi_vs_t file exists 
    if not os.path.isfile(path.phi_vs_t) and os.path.isfile(path.output_stella):   
        write_txtFileForPotentialVsTime(path.folder) 
        
    # Read the txt file for a nonlinear simulation
    if path.nonlinear: 
        return read_potentialFromTxtFile(path) 
    
    # For a linear simulation, you should read phi2_vs_tkxky (potential 4D) instead
    if path.linear:
        exit_reason = "For linear simulations read phi2_vs_tkxky instead of phi2_vs_t." 
        exit_program(exit_reason, read_potentialVsTime, sys._getframe().f_lineno)   
    
    # Critical error if we didn't find any data
    exit_reason = "The potential data could not be found for:\n"
    exit_reason += "     "+str(path.input_file)
    exit_program(exit_reason, read_potentialVsTime, sys._getframe().f_lineno)   
    return

#-----------------------------
def read_potentialFromTxtFile(path): 
    data = np.loadtxt(path.phi_vs_t,skiprows=1,dtype='float').reshape(-1, 6)
    potential_data = {} 
    potential_data["vec_time"]  = data[:,0]
    potential_data["phi2_vs_t"] = data[:,1]
    potential_data["phi2_vs_t_zonal"] = data[:,2]
    potential_data["phi2_vs_t_nozonal"] = data[:,3]
    potential_data["phi_vs_t"]  = data[:,4] + 1j*data[:,5]
    return potential_data

#===============================================================================
#                ATTACH THE DPHI(Z) DATA TO THE SIMULATION OBJECT              #
#===============================================================================

def get_potentialVsTime(self): 
    
    # Read the data in the potential file 
    potential = read_potentialVsTime(self.path)  
    
    # Save the distribution data 
    self.phi_vs_t = Data(["phi", "t"], potential["phi_vs_t"], potential["vec_time"]) 
    self.phi2_vs_t = Data(["phi2", "t"], potential["phi2_vs_t"], self.phi_vs_t.t)  
    self.phi2_vs_t_zonal = Data(["phi2_zonal", "t"], potential["phi2_vs_t_zonal"], self.phi_vs_t.t)  
    self.phi2_vs_t_nozonal = Data(["phi2_nozonal", "t"], potential["phi2_vs_t_nozonal"], self.phi_vs_t.t) 
    return 