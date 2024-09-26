
import numpy as np
import os, sys, h5py  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program  
from stellapy.data.paths.load_pathObject import get_potential5DPath
from stellapy.data.potential.write_h5FileForPotential5D import write_h5FileForPotential5D

#===============================================================================
#                        ATTACH THE POTENTIAL DATA
#===============================================================================

def read_trappedWeights5D(path): 
    ''' Read the trapped particle weights from *.trappedWeights5D or *.out.nc ''' 
               
    # Make sure that the *.dt10.trappedWeights5D file exists
    if not os.path.isfile(path.trappedWeights5D) and os.path.isfile(path.output_stella): 
        write_h5FileForPotential5D(path.folder)
        get_potential5DPath(path)  
                              
    # Read from the *.dt10.trappedWeights5D file 
    if os.path.isfile(path.trappedWeights5D):
        return read_fromPotentialFile(path.trappedWeights5D)    

    # Critical error if we didn't find any data
    exit_reason = "The potential data could not be found."
    exit_program(exit_reason, read_trappedWeights5D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromPotentialFile(path, potential={}):   
    with h5py.File(path, 'r') as f: 
        for key in ["vec_time", "trappedw_vs_tszkxky"]:  
            if key in f.keys():  
                potential[key] = f[key][()] 
    return potential 
            
#===============================================================================
#                  ATTACH THE POTENTIAL TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_trappedWeights5D(self): 
    
    # Read the potential
    potential = read_trappedWeights5D(self.path)  
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
    sorted_indexes = np.argsort(self.vec.kx_stella)  
    potential["trappedw_vs_tszkxky"] = potential["trappedw_vs_tszkxky"][:,:,:,sorted_indexes,:] 
    
    # Save the potential
    self.trappedw_vs_tszkxky = Data(["trappedw","t","s","z","kx","ky"], potential["trappedw_vs_tszkxky"], potential["vec_time"], self.vec.z, self.vec.kx, self.vec.ky) 
    return 
