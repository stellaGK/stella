
import numpy as np
import os, sys, h5py  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program 
from stellapy.data.output.read_outputFile import read_outputFile, read_netcdfVariables 

#===============================================================================
#                        ATTACH THE DIMENSIONS DATA
#===============================================================================

def read_moments5D(path):
    ''' Read the moments from *.moments or *.out.nc ''' 
                                     
    # Read from the *.moments file 
    if os.path.isfile(path.moments5D):
        return read_fromMomentsFile(path.moments5D)  
    
    # Read from the *.out.nc file 
    if os.path.isfile(path.output_stella):
        return read_fromNcFile(path.output_stella)   

    # Critical error if we didn't find any data
    exit_reason = "The 5D moments data could not be found."
    exit_program(exit_reason, read_moments5D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromMomentsFile(path):  
    moments={}
    with h5py.File(path, 'r') as f: 
        for key in ["vec_time", "upar_vs_tszkxky", "dens_vs_tszkxky", "temp_vs_tszkxky"]: 
            if key in f.keys():  
                moments[key] = f[key][()]  
    return moments 

#-------------------------------------
def read_fromNcFile(path):  
    moments={}
    netcdf_file = read_outputFile(path)   
    moments['vec_time'] = read_netcdfVariables('vec_time', netcdf_file)
    moments['upar_vs_tszkxky'] = read_netcdfVariables('upar_vs_tszkxky', netcdf_file)
    moments['dens_vs_tszkxky'] = read_netcdfVariables('dens_vs_tszkxky', netcdf_file)
    moments['temp_vs_tszkxky'] = read_netcdfVariables('temp_vs_tszkxky', netcdf_file)  
    netcdf_file.close()
    return moments
            
#===============================================================================
#                  ATTACH THE MOMENTS TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_moments5D(self): 
    
    # Read the moments
    moments = read_moments5D(self.path) 
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
    sorted_indexes = np.argsort(self.vec.kx_stella)  
    moments["upar_vs_tszkxky"] = moments["upar_vs_tszkxky"][:,:,:,sorted_indexes,:]
    moments["dens_vs_tszkxky"] = moments["dens_vs_tszkxky"][:,:,:,sorted_indexes,:]
    moments["temp_vs_tszkxky"] = moments["temp_vs_tszkxky"][:,:,:,sorted_indexes,:] 
    
    # Save the moments
    self.upar_vs_tszkxky = Data(["upar","t","s","z","kx","ky"], moments["upar_vs_tszkxky"], moments["vec_time"], self.vec.species, self.vec.z, self.vec.kx, self.vec.ky) 
    self.dens_vs_tszkxky = Data(["dens","t","s","z","kx","ky"], moments["dens_vs_tszkxky"], moments["vec_time"], self.vec.species, self.vec.z, self.vec.kx, self.vec.ky) 
    self.temp_vs_tszkxky = Data(["temp","t","s","z","kx","ky"], moments["temp_vs_tszkxky"], moments["vec_time"], self.vec.species, self.vec.z, self.vec.kx, self.vec.ky) 
    return 

