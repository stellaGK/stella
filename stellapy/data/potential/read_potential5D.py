
import numpy as np
import os, sys, h5py  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program  
from stellapy.data.paths.load_pathObject import get_potential5DPath
from stellapy.data.potential.write_h5FileForPotential5D import write_h5FileForPotential5D

#===============================================================================
#                        ATTACH THE POTENTIAL DATA
#===============================================================================

def read_potential5D(path): 
    ''' Read the potential from *.potential5D or *.out.nc ''' 
               
    # Make sure that the *.dt10.potential5D file exists
    if not os.path.isfile(path.potential5D) and os.path.isfile(path.output_stella): 
        write_h5FileForPotential5D(path.folder)
        get_potential5DPath(path)  
                              
    # Read from the *.dt10.potential5D file 
    if os.path.isfile(path.potential5D):
        return read_fromPotentialFile(path.potential5D)    

    # Critical error if we didn't find any data
    exit_reason = "The potential data could not be found."
    exit_program(exit_reason, read_potential5D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromPotentialFile(path):
    potential={}   
    with h5py.File(path, 'r') as f: 
        for key in ["vec_time", "phi_vs_tzkxky", "phi_vs_tzkxky_tavg"]:  
            if key in f.keys():  
                potential[key] = f[key][()]  
    return potential 
            
#===============================================================================
#                  ATTACH THE POTENTIAL TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_potential5D(self): 
    
    # Read the potential
    potential = read_potential5D(self.path)  
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
    sorted_indexes = np.argsort(self.vec.kx_stella)  
    potential["phi_vs_tzkxky"] = potential["phi_vs_tzkxky"][:,:,sorted_indexes,:]
    potential["phi_vs_tzkxky_tavg"] = potential["phi_vs_tzkxky_tavg"][:,:,sorted_indexes,:] 
    
    # Save the potential
    self.phi_vs_tzkxky = Data(["phi","t","z","kx","ky"], potential["phi_vs_tzkxky"], potential["vec_time"], self.vec.z, self.vec.kx, self.vec.ky) 
    self.phi_vs_tzkxky_tavg = Data(["phi_tavg","t","z","kx","ky"], potential["phi_vs_tzkxky_tavg"], self.phi_vs_tzkxky.t, self.vec.z, self.vec.kx, self.vec.ky)     
    return 
