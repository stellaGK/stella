#!/usr/bin/python3 
import os
import sys 
import h5py
import pathlib
import numpy as np

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.utils.decorators.exit_program import exit_program    
from stellapy.data.utils import Data

#===============================================================================
#                        READ THE 4D POTENTIAL DATA
#===============================================================================

def read_potential4DLinear(path, dim, vec, nakxnaky): 
    ''' Read the potential from *.potential4D or *.out.nc '''  
        
    # Read the potential files from multiple simulations 
    if path.multiple_input_files: 
        return read_potentialFromMultipleFiles(dim, vec, path.paths) 
    
    # Read the h5 or txt file      
    if os.path.isfile(path.omega):  
        if nakxnaky>1: return read_potentialFromH5File(path.phi_vs_z)  
        if nakxnaky==1: return read_potentialFromTxtFile(path.phi_vs_z)  

    # Critical error if we didn't find any data
    exit_reason = "The potential data versus (z,kx,ky) could not be found for:\n"
    exit_reason += "     "+str(path.input_file)
    exit_program(exit_reason, read_potential4DLinear, sys._getframe().f_lineno)   
    return 

#-------------------------------------
def read_potentialFromH5File(path_phi_vs_z): 
    with h5py.File(path_phi_vs_z, 'r') as f:   
        phi_vs_zkxky = f["phi_vs_zkxky"][()]  
        phi2_vs_zkxky = f["phi2_vs_zkxky"][()]  
    return {"phi_vs_zkxky" : phi_vs_zkxky, "phi2_vs_zkxky" : phi2_vs_zkxky} 
   
#-----------------------------
def read_potentialFromTxtFile(path_phi_vs_z):  
    data = np.loadtxt(path_phi_vs_z,skiprows=1,dtype='float').reshape(-1, 3) 
    phi_vs_z = data[:,1] + 1j*data[:,2]; phi2_vs_z = data[:,0]  
    return {"phi_vs_zkxky" : phi_vs_z[:,np.newaxis,np.newaxis] , "phi2_vs_zkxky" : phi2_vs_z[:,np.newaxis,np.newaxis] }

#-----------------------------
def read_potentialFromMultipleFiles(dim, vec, paths): 
        
    # Create matrices (kx,ky)
    phi_vs_zkxky = np.ones((dim.z, dim.kx, dim.ky), dtype=complex)*np.nan
    phi2_vs_zkxky = np.ones((dim.z, dim.kx, dim.ky))*np.nan 
    
    # Iterate over the simulations (1 mode per simulation)
    for path in paths:
        
        # Don't read the dummy simulation, we need to read each input file inside of it
        if "_dummy.in" not in str(path.input_file):
        
            # Read dphiz(t) from a txt or an h5 file
            if path.nakxnaky==1: data = read_potentialFromTxtFile(path.phi_vs_z)
            if path.nakxnaky>1: data = read_potentialFromH5File(path.phi_vs_z)   

            # Iterate over the modes (kx,ky) 
            for iikx, kx in enumerate(path.vec_kx):
                for iiky, ky in enumerate(path.vec_ky):  
                    
                    # Get the mode (kx,ky)   
                    try:
                        ikx = list(vec.kx).index(kx) 
                        iky = list(vec.ky).index(ky)
                    except: 
                        exit_reason = f"The mode (kx,ky) = ({kx},{ky}) could not be found inside:\n"
                        exit_reason += "     "+str(path.phi_vs_z)+"\n"
                        exit_reason += "     "+str(path.input_file)+"\n"
                        exit_reason += f"     vec_kx = {path.vec_kx}\n"
                        exit_reason += f"     vec_ky = {path.vec_ky}\n"
                        exit_program(exit_reason, read_potentialFromMultipleFiles, sys._getframe().f_lineno) 
                            
                    # Put the omega data in the matrices
                    phi_vs_zkxky[:,ikx,iky] = data["phi_vs_zkxky"][:,iikx,iiky]  
                    phi2_vs_zkxky[:,ikx,iky] = data["phi2_vs_zkxky"][:,iikx,iiky]  
 
    return {"phi_vs_zkxky" : phi_vs_zkxky, "phi2_vs_zkxky" : phi2_vs_zkxky}
            
#===============================================================================
#                  ATTACH THE POTENTIAL TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_potential4DLinear(self): 
    
    # Read the potential
    potential = read_potential4DLinear(self.path, self.dim, self.vec, self.nakxnaky)  
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
    sorted_indexes = np.argsort(self.vec.kx_stella)  
    potential["phi_vs_zkxky"] = potential["phi_vs_zkxky"][:,sorted_indexes,:]
    potential["phi2_vs_zkxky"] = potential["phi2_vs_zkxky"][:,sorted_indexes,:]  
    
    # Save the potential
    self.phi_vs_zkxky = Data(["phi","z","kx","ky"], potential["phi_vs_zkxky"], self.vec.z, self.vec.kx, self.vec.ky) 
    self.phi2_vs_zkxky = Data(["phi2","z","kx","ky"], potential["phi2_vs_zkxky"], self.vec.z, self.vec.kx, self.vec.ky)    
    return 
