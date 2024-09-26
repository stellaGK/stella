
import numpy as np
import os, sys, h5py  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program 
from stellapy.data.paths.load_pathObject import get_moments4DPath
from stellapy.data.moments.write_h5FileForMoments4D import write_h5FileForMoments4D

#===============================================================================
#                        ATTACH THE DIMENSIONS DATA
#===============================================================================

def read_moments4D(path):
    ''' Read the moments from *.moments or *.out.nc ''' 
               
    # Make sure that the *.moments4D file exists
    if not os.path.isfile(path.potential4D) and os.path.isfile(path.output_stella): 
        write_h5FileForMoments4D(path.folder)
        get_moments4DPath(path)   
                           
    # Read from the *.moments4D file 
    if os.path.isfile(path.moments4D):
        return read_fromMomentsFile(path.moments4D)  

    # Critical error if we didn't find any data
    exit_reason = "The 4D moments data could not be found."
    exit_program(exit_reason, read_moments4D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromMomentsFile(path): 
    moments={}
    with h5py.File(path, 'r') as f: 
        for key in ["vec_time", "upar_vs_tskxky", "dens_vs_tskxky", "temp_vs_tskxky", "upar2_vs_tskxky", "dens2_vs_tskxky", "temp2_vs_tskxky"]: 
            for extra in ["", "_zeta0"]:
                if key+extra in f.keys():  
                    moments[key+extra] = f[key+extra][()]  
    return moments 
            
#===============================================================================
#                  ATTACH THE MOMENTS TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_moments4D(self): 
    
    # Read the moments
    moments = read_moments4D(self.path) 
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
    sorted_indexes = np.argsort(self.vec.kx_stella)  
    moments["upar_vs_tskxky"] = moments["upar_vs_tskxky"][:,:,sorted_indexes,:]
    moments["dens_vs_tskxky"] = moments["dens_vs_tskxky"][:,:,sorted_indexes,:]
    moments["temp_vs_tskxky"] = moments["temp_vs_tskxky"][:,:,sorted_indexes,:]
    if "upar2_vs_tskxky" in moments:
        moments["upar2_vs_tskxky"] = moments["upar2_vs_tskxky"][:,:,sorted_indexes,:]
        moments["dens2_vs_tskxky"] = moments["dens2_vs_tskxky"][:,:,sorted_indexes,:]
        moments["temp2_vs_tskxky"] = moments["temp2_vs_tskxky"][:,:,sorted_indexes,:]
    if "upar_vs_tskxky_zeta0" in moments:
        moments["upar_vs_tskxky_zeta0"] = moments["upar_vs_tskxky_zeta0"][:,:,sorted_indexes,:]
        moments["dens_vs_tskxky_zeta0"] = moments["dens_vs_tskxky_zeta0"][:,:,sorted_indexes,:]
        moments["temp_vs_tskxky_zeta0"] = moments["temp_vs_tskxky_zeta0"][:,:,sorted_indexes,:]
         
    # Save the moments
    self.upar_vs_tskxky = Data(["upar","t","s","kx","ky"], moments["upar_vs_tskxky"], moments["vec_time"], self.vec.species, self.vec.kx, self.vec.ky) 
    self.dens_vs_tskxky = Data(["dens","t","s","kx","ky"], moments["dens_vs_tskxky"], moments["vec_time"], self.vec.species, self.vec.kx, self.vec.ky) 
    self.temp_vs_tskxky = Data(["temp","t","s","kx","ky"], moments["temp_vs_tskxky"], moments["vec_time"], self.vec.species, self.vec.kx, self.vec.ky) 
    if "upar_vs_tskxky_zeta0" in moments:
        self.upar_vs_tskxky_zeta0 = Data(["upar","t","s","kx","ky"], moments["upar_vs_tskxky_zeta0"], moments["vec_time"], self.vec.species, self.vec.kx, self.vec.ky) 
        self.dens_vs_tskxky_zeta0 = Data(["dens","t","s","kx","ky"], moments["dens_vs_tskxky_zeta0"], moments["vec_time"], self.vec.species, self.vec.kx, self.vec.ky) 
        self.temp_vs_tskxky_zeta0 = Data(["temp","t","s","kx","ky"], moments["temp_vs_tskxky_zeta0"], moments["vec_time"], self.vec.species, self.vec.kx, self.vec.ky) 
    else:
        self.upar_vs_tskxky_zeta0 = None
        self.dens_vs_tskxky_zeta0 = None
        self.temp_vs_tskxky_zeta0 = None    
    if "upar2_vs_tskxky" in moments:
        self.upar2_vs_tskxky = Data(["upar2","t","s","kx","ky"], moments["upar2_vs_tskxky"], moments["vec_time"], self.vec.species, self.vec.kx, self.vec.ky) 
        self.dens2_vs_tskxky = Data(["dens2","t","s","kx","ky"], moments["dens2_vs_tskxky"], moments["vec_time"], self.vec.species, self.vec.kx, self.vec.ky) 
        self.temp2_vs_tskxky = Data(["temp2","t","s","kx","ky"], moments["temp2_vs_tskxky"], moments["vec_time"], self.vec.species, self.vec.kx, self.vec.ky) 
    else:
        self.upar2_vs_tskxky = None
        self.dens2_vs_tskxky = None
        self.temp2_vs_tskxky = None
    return 

