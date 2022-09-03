
import numpy as np
import os, sys, h5py  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program  
from stellapy.data.paths.load_pathObject import get_fluxes3DPath
from stellapy.data.fluxes.write_h5FileForFluxes3D import write_h5FileForFluxes3D

#===============================================================================
#                        ATTACH THE FLUXES DATA
#===============================================================================

def read_fluxes3D(path): 
    ''' Read the fluxes from *.fluxes3D or *.out.nc ''' 
               
    # Make sure that the *.dt10.write_h5FileForfluxes3D file exists
    if not os.path.isfile(path.fluxes3D) and os.path.isfile(path.output_stella): 
        write_h5FileForFluxes3D(path.folder)
        get_fluxes3DPath(path)  

    # Read from the *.dt10.fluxes3D file 
    if os.path.isfile(path.fluxes3D):
        return read_fromfluxesFile(path.fluxes3D)    
    
    # This file doesn't exist for old simulations
    elif os.path.isfile(path.output): 
        return read_fromH5File(path)
    
    # Critical error if we didn't find any data 
    exit_reason = "The fluxes data could not be found."
    exit_program(exit_reason, read_fluxes3D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromfluxesFile(path, fluxes={}):
    with h5py.File(path, 'r') as f:  
        for key in ["vec_time", "qflux_vs_ts", "pflux_vs_ts", "vflux_vs_ts"]:
            for extra in ["", "kx", "ky","z"]:  
                if key+extra in f.keys(): 
                    fluxes[key+extra] = f[key+extra][()]   
    return fluxes 

#===============================================================================
#                            READ FROM OLD H5 FILE                             #
#===============================================================================

def read_fromH5File(path, fluxes={}):   
    with h5py.File(path.output, 'r') as f: 
        fluxes["vec_time"] = f["vec_time"][()]
        qflux_vs_tskxky = f["qflux_vs_tskxky"][()] if ("qflux_vs_tskxky" in f.keys()) else f["qflx_kxky"][()]   
        fluxes["qflux_vs_tskx"] = qflux_vs_tskxky[:,:,:,0] + 2*np.sum(qflux_vs_tskxky[:,:,:,1:],axis=3) 
        fluxes["qflux_vs_tsky"] = np.sum(qflux_vs_tskxky[:,:,:,:],axis=2)
        fluxes["pflux_vs_tskx"] = np.ones((np.shape(fluxes["qflux_vs_tskx"])))*np.nan
        fluxes["vflux_vs_tskx"] = np.ones((np.shape(fluxes["qflux_vs_tskx"])))*np.nan
        fluxes["pflux_vs_tsky"] = np.ones((np.shape(fluxes["qflux_vs_tsky"])))*np.nan
        fluxes["vflux_vs_tsky"] = np.ones((np.shape(fluxes["qflux_vs_tsky"])))*np.nan
    return fluxes
      
#===============================================================================
#                  ATTACH THE fluxes TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_fluxes3D(self): 
    
    # Read the fluxes
    fluxes = read_fluxes3D(self.path)   
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
    if "qflux_vs_tskx" in fluxes:
        sorted_indexes = np.argsort(self.vec.kx_stella)  
        fluxes["qflux_vs_tskx"] = fluxes["qflux_vs_tskx"][:,:,sorted_indexes]
        fluxes["pflux_vs_tskx"] = fluxes["pflux_vs_tskx"][:,:,sorted_indexes]
        fluxes["vflux_vs_tskx"] = fluxes["vflux_vs_tskx"][:,:,sorted_indexes]

    # Save the fluxes
    if "qflux_vs_tskx" in fluxes:
        self.qflux_vs_tskx = Data(["qflux","t","s","kx"], fluxes["qflux_vs_tskx"], fluxes["vec_time"], range(self.dim.species), self.vec.kx) 
        self.pflux_vs_tskx = Data(["pflux","t","s","kx"], fluxes["pflux_vs_tskx"], self.qflux_vs_tskx.t, range(self.dim.species), self.vec.kx) 
        self.vflux_vs_tskx = Data(["vflux","t","s","kx"], fluxes["vflux_vs_tskx"], self.qflux_vs_tskx.t, range(self.dim.species), self.vec.kx) 
        self.qflux_vs_tsky = Data(["qflux","t","s","ky"], fluxes["qflux_vs_tsky"], self.qflux_vs_tskx.t, range(self.dim.species), self.vec.ky) 
        self.pflux_vs_tsky = Data(["pflux","t","s","ky"], fluxes["pflux_vs_tsky"], self.qflux_vs_tskx.t, range(self.dim.species), self.vec.ky) 
        self.vflux_vs_tsky = Data(["vflux","t","s","ky"], fluxes["vflux_vs_tsky"], self.qflux_vs_tskx.t, range(self.dim.species), self.vec.ky) 
    else:
        self.qflux_vs_tskx = np.nan
        self.pflux_vs_tskx = np.nan
        self.vflux_vs_tskx = np.nan
        self.qflux_vs_tsky = np.nan
        self.pflux_vs_tsky = np.nan
        self.vflux_vs_tsky = np.nan
    if "qflux_vs_tsz" in fluxes:
        self.qflux_vs_tsz = Data(["qflux","t","s","z"], fluxes["qflux_vs_tsz"], fluxes["vec_time"], range(self.dim.species), self.vec.z) 
        self.pflux_vs_tsz = Data(["pflux","t","s","z"], fluxes["pflux_vs_tsz"], self.qflux_vs_tsz.t, range(self.dim.species), self.vec.z) 
        self.vflux_vs_tsz = Data(["vflux","t","s","z"], fluxes["vflux_vs_tsz"], self.qflux_vs_tsz.t, range(self.dim.species), self.vec.z) 
    else:
        dummy = np.ones((len(self.qflux_vs_tskx.t), self.dim.species, self.dim.z))*np.nan
        self.qflux_vs_tsz = Data(["qflux","t","s","z"], dummy, self.qflux_vs_tskx.t, range(self.dim.species), self.vec.z) 
        self.pflux_vs_tsz = Data(["pflux","t","s","z"], dummy, self.qflux_vs_tskx.t, range(self.dim.species), self.vec.z) 
        self.vflux_vs_tsz = Data(["vflux","t","s","z"], dummy, self.qflux_vs_tskx.t, range(self.dim.species), self.vec.z)
    return 
