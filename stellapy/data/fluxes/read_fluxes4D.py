
import numpy as np
import os, sys, h5py  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program  
from stellapy.data.paths.load_pathObject import get_fluxes4DPath
from stellapy.data.fluxes.write_h5FileForFluxes4D import write_h5FileForFluxes4D

#===============================================================================
#                        ATTACH THE FLUXES DATA
#===============================================================================

def read_fluxes4D(path): 
    ''' Read the fluxes from *.fluxes4D or *.out.nc ''' 
               
    # Make sure that the *.dt10.write_h5FileForfluxes4D file exists
    if not os.path.isfile(path.fluxes4D) and os.path.isfile(path.output_stella): 
        write_h5FileForFluxes4D(path.folder)
        get_fluxes4DPath(path)  
                              
    # Read from the *.dt10.fluxes4D file 
    if os.path.isfile(path.fluxes4D):
        return read_fromfluxesFile(path.fluxes4D)    
    
    # Critical error if we didn't find any data 
    exit_reason = "The fluxes data could not be found for:\n"
    exit_reason += "     "+str(path.input_file)
    exit_program(exit_reason, read_fluxes4D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromfluxesFile(path, fluxes={}):
    with h5py.File(path, 'r') as f:  
        for key in ["vec_time", "qflux_vs_ts", "pflux_vs_ts", "vflux_vs_ts"]:
            for extra in ["", "kxky"]: 
                if key+extra in f.keys(): fluxes[key+extra] = f[key+extra][()]   
                else: fluxes[key+extra] = None 
    return fluxes 
      
#===============================================================================
#                  ATTACH THE fluxes TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_fluxes4D(self): 
    
    # Read the fluxes
    fluxes = read_fluxes4D(self.path)   
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
    sorted_indexes = np.argsort(self.vec.kx_stella)  
    fluxes["qflux_vs_tskxky"] = fluxes["qflux_vs_tskxky"][:,:,sorted_indexes,:]
    fluxes["pflux_vs_tskxky"] = fluxes["pflux_vs_tskxky"][:,:,sorted_indexes,:]
    fluxes["vflux_vs_tskxky"] = fluxes["vflux_vs_tskxky"][:,:,sorted_indexes,:]

    # Save the fluxes
    self.qflux_vs_tskxky = Data(["qflux","t","s","kx","ky"], fluxes["qflux_vs_tskxky"], fluxes["vec_time"], range(self.dim.species), self.vec.kx, self.vec.ky) 
    self.pflux_vs_tskxky = Data(["pflux","t","s","kx","ky"], fluxes["pflux_vs_tskxky"], self.qflux_vs_tskx.t, range(self.dim.species), self.vec.kx, self.vec.ky) 
    self.vflux_vs_tskxky = Data(["vflux","t","s","kx","ky"], fluxes["vflux_vs_tskxky"], self.qflux_vs_tskx.t, range(self.dim.species), self.vec.kx, self.vec.ky)
    return 
