
import numpy as np
import os, sys, h5py  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program  
from stellapy.data.potential.write_h5FileForPotential4D import write_h5FileForPotential4D
from stellapy.data.paths.load_pathObject import get_potential4DPath

#===============================================================================
#                        ATTACH THE POTENTIAL DATA
#===============================================================================

def read_potential4D(path): 
    ''' Read the potential from *.potential4D or *.out.nc ''' 
               
    # Make sure that the *.dt10.potential4D file exists
    if not os.path.isfile(path.potential4D) and os.path.isfile(path.output_stella): 
        write_h5FileForPotential4D(path.folder)
        get_potential4DPath(path)  
                              
    # Read from the *.dt10.potential4D file 
    if os.path.isfile(path.potential4D):
        return read_fromPotentialFile(path.potential4D)    

    # Critical error if we didn't find any data
    exit_reason = "The potential data could not be found for:\n"
    exit_reason += "     "+str(path.input_file)
    exit_program(exit_reason, read_potential4D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromPotentialFile(path):
    potential={}   
    with h5py.File(path, 'r') as f: 
        for key in ["vec_time", "phi_vs_tkxky", "phi2_vs_tkxky", "phi_vs_tkxky_zeta0"]:  
            if key in f.keys():  
                potential[key] = f[key][()]  
    return potential 
            
#===============================================================================
#                  ATTACH THE POTENTIAL TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_potential4D(self): 
    
    # Read the potential
    potential = read_potential4D(self.path)  
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
    sorted_indexes = np.argsort(self.vec.kx_stella)  
    potential["phi_vs_tkxky"] = potential["phi_vs_tkxky"][:,sorted_indexes,:]
    potential["phi2_vs_tkxky"] = potential["phi2_vs_tkxky"][:,sorted_indexes,:] 
    if "phi_vs_tkxky_zeta0" in potential:
        potential["phi_vs_tkxky_zeta0"] = potential["phi_vs_tkxky_zeta0"][:,sorted_indexes,:] 
    
    # Save the potential
    self.phi_vs_tkxky = Data(["phi","t","kx","ky"], potential["phi_vs_tkxky"], potential["vec_time"], self.vec.kx, self.vec.ky) 
    self.phi2_vs_tkxky = Data(["phi2","t","kx","ky"], potential["phi2_vs_tkxky"], self.phi_vs_tkxky.t, self.vec.kx, self.vec.ky)    
    if "phi_vs_tkxky_zeta0" in potential:
        self.phi_vs_tkxky_zeta0 = Data(["phi","t","kx","ky"], potential["phi_vs_tkxky_zeta0"], self.phi_vs_tkxky.t, self.vec.kx, self.vec.ky)    
    else:
        self.phi_vs_tkxky_zeta0 = None
    return 
