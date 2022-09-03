
import numpy as np
import os, sys, h5py  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program 
from stellapy.data.paths.load_pathObject import get_moments3DPath
from stellapy.data.moments.write_h5FileForMoments3D import write_h5FileForMoments3D

#===============================================================================
#                        ATTACH THE DIMENSIONS DATA
#===============================================================================

def read_moments3D(path):
    ''' Read the moments from *.moments or *.out.nc ''' 
               
    # Make sure that the *.moments3D file exists
    if not os.path.isfile(path.moments3D) and os.path.isfile(path.output_stella): 
        write_h5FileForMoments3D(path.folder)
        get_moments3DPath(path)   
                           
    # Read from the *.moments3D file 
    if os.path.isfile(path.moments3D):
        return read_fromMomentsFile(path.moments3D)  

    # Critical error if we didn't find any data
    exit_reason = "The 3D moments data could not be found."
    exit_program(exit_reason, read_moments3D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromMomentsFile(path):
    moments={}
    with h5py.File(path, 'r') as f: 
        for key in ["vec_time", "upar_vs_ts", "dens_vs_ts", "temp_vs_ts", "upar2_vs_ts", "dens2_vs_ts", "temp2_vs_ts"]: 
            for extra in ["", "z", "kx","ky"]: 
                if key+extra in f.keys():   
                    moments[key+extra] = f[key+extra][()]   
    return moments 
            
#===============================================================================
#                  ATTACH THE MOMENTS TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_moments3D(self): 
    
    # Read the moments
    moments = read_moments3D(self.path) 
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
    if "upar_vs_tskx" in moments:
        sorted_indexes = np.argsort(self.vec.kx_stella)  
        moments["upar_vs_tskx"] = moments["upar_vs_tskx"][:,:,sorted_indexes]
        moments["dens_vs_tskx"] = moments["dens_vs_tskx"][:,:,sorted_indexes]
        moments["temp_vs_tskx"] = moments["temp_vs_tskx"][:,:,sorted_indexes] 
    
    # Save the moments
    self.upar_vs_tsz = Data(["upar","t","s","z"], moments["upar_vs_tsz"], moments["vec_time"], self.vec.species, self.vec.z) 
    self.dens_vs_tsz = Data(["dens","t","s","z"], moments["dens_vs_tsz"], moments["vec_time"], self.vec.species, self.vec.z) 
    self.temp_vs_tsz = Data(["temp","t","s","z"], moments["temp_vs_tsz"], moments["vec_time"], self.vec.species, self.vec.z) 
    if "upar_vs_tskx" in moments:
        self.upar_vs_tskx = Data(["upar","t","s","kx"], moments["upar_vs_tskx"], moments["vec_time"], self.vec.species, self.vec.kx) 
        self.dens_vs_tskx = Data(["dens","t","s","kx"], moments["dens_vs_tskx"], moments["vec_time"], self.vec.species, self.vec.kx) 
        self.temp_vs_tskx = Data(["temp","t","s","kx"], moments["temp_vs_tskx"], moments["vec_time"], self.vec.species, self.vec.kx)
        self.upar_vs_tsky = Data(["upar","t","s","ky"], moments["upar_vs_tsky"], moments["vec_time"], self.vec.species, self.vec.ky) 
        self.dens_vs_tsky = Data(["dens","t","s","ky"], moments["dens_vs_tsky"], moments["vec_time"], self.vec.species, self.vec.ky) 
        self.temp_vs_tsky = Data(["temp","t","s","ky"], moments["temp_vs_tsky"], moments["vec_time"], self.vec.species, self.vec.ky)
    else:
        self.upar_vs_tskx = np.nan
        self.dens_vs_tskx = np.nan
        self.temp_vs_tskx = np.nan
        self.upar_vs_tsky = np.nan
        self.dens_vs_tsky = np.nan
        self.temp_vs_tsky = np.nan
    if "upar2_vs_tsz" in moments:
        self.upar2_vs_tsz = Data(["upar2","t","s","z"], moments["upar2_vs_tsz"], moments["vec_time"], self.vec.species, self.vec.z) 
        self.dens2_vs_tsz = Data(["dens2","t","s","z"], moments["dens2_vs_tsz"], moments["vec_time"], self.vec.species, self.vec.z) 
        self.temp2_vs_tsz = Data(["temp2","t","s","z"], moments["temp2_vs_tsz"], moments["vec_time"], self.vec.species, self.vec.z) 
    else:
        self.upar2_vs_tsz = np.nan
        self.dens2_vs_tsz = np.nan
        self.dens2_vs_tsz = np.nan
    return 

