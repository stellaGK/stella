
#!/usr/bin/python3  
import sys, os   
import pathlib
import numpy as np

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)    
from stellapy.data.fluxes.write_txtFileForFluxesVsZ import write_txtFileForFluxesVsZ
from stellapy.utils.decorators.exit_program import exit_program   
from stellapy.data.utils import Data

#===============================================================================
#                        ATTACH THE FLUXES DATA
#===============================================================================

def read_fluxesVsZ(path, dim, vec): 
    ''' Read the fluxes from *.fluxes_vs_z file. ''' 
               
    # Make sure that the *.fluxes_vs_z file exists
    if not os.path.isfile(path.fluxes_vs_z) and os.path.isfile(path.output_stella): 
        write_txtFileForFluxesVsZ(path.folder)
        
    # Read the fluxes files from multiple simulations 
    if path.multiple_input_files: 
        return read_fluxesVsZFromMultipleFiles(dim, vec, path.paths) 
                              
    # Read the txt file for one simulation
    if os.path.isfile(path.fluxes_vs_z):
        return read_fluxesVsZFromTxtFile(path.fluxes_vs_z, dim)    
    
    # Critical error if we didn't find any data 
    exit_reason = "The quasi-linear fluxes data could not be found for:\n"
    exit_reason += "     "+str(path.input_file)
    exit_program(exit_reason, read_fluxesVsZ, sys._getframe().f_lineno)   
    return

#-----------------------------
def read_fluxesVsZFromTxtFile(path_fluxes_vs_z, dim):  
    data = np.loadtxt(path_fluxes_vs_z,skiprows=1,dtype='float').reshape(-1, 1+3*dim.species)   
    data_dict = { 
        "phi2_vs_zkxky" : np.ones((dim.z, 1, 1))*data[:,0][:,np.newaxis,np.newaxis],\
        "pflux_vs_szkxky" : np.ones((dim.species, dim.z, 1, 1))*np.nan,\
        "qflux_vs_szkxky" : np.ones((dim.species, dim.z, 1, 1))*np.nan,\
        "vflux_vs_szkxky" : np.ones((dim.species, dim.z, 1, 1))*np.nan}
    for i in range(dim.species):
        data_dict["pflux_vs_szkxky"][i,:,0,0] = data[:,0*dim.species+1+i]
        data_dict["vflux_vs_szkxky"][i,:,0,0] = data[:,1*dim.species+1+i]
        data_dict["qflux_vs_szkxky"][i,:,0,0] = data[:,2*dim.species+1+i]  
    return data_dict

#-----------------------------
def read_fluxesVsZFromMultipleFiles(dim, vec, paths): 
        
    # Create matrices (s,kx,ky)
    phi2_vs_zkxky = np.ones((dim.z, dim.kx, dim.ky))*np.nan
    pflux_vs_szkxky = np.ones((dim.species, dim.z, dim.kx, dim.ky))*np.nan
    qflux_vs_szkxky = np.ones((dim.species, dim.z, dim.kx, dim.ky))*np.nan
    vflux_vs_szkxky = np.ones((dim.species, dim.z, dim.kx, dim.ky))*np.nan 
    
    # Iterate over the simulations (1 mode per simulation)
    for path in paths:
        
        # Don't read the dummy simulation, we need to read each input file inside of it
        if "_dummy.in" not in str(path.input_file):
        
            # Read QLfluxes(t) from a txt file
            data = read_fluxesVsZFromTxtFile(path.fluxes_vs_z, dim)

            # Iterate over the modes (kx,ky) 
            for iikx, kx in enumerate(path.vec_kx):
                for iiky, ky in enumerate(path.vec_ky): 
                    
                    # Get the mode (kx,ky)   
                    ikx = list(vec.kx).index(kx) 
                    iky = list(vec.ky).index(ky)
                        
                    # Put the quasi-linear fluxes data in the matrices
                    phi2_vs_zkxky[:,ikx,iky] = data["phi2_vs_zkxky"][:,iikx,iiky]  
                    pflux_vs_szkxky[:,:,ikx,iky] = data["pflux_vs_szkxky"][:,:,iikx,iiky]  
                    qflux_vs_szkxky[:,:,ikx,iky] = data["qflux_vs_szkxky"][:,:,iikx,iiky]  
                    vflux_vs_szkxky[:,:,ikx,iky] = data["vflux_vs_szkxky"][:,:,iikx,iiky]  

    return {"phi2_vs_zkxky" : phi2_vs_zkxky, "pflux_vs_szkxky" : pflux_vs_szkxky, "qflux_vs_szkxky" : qflux_vs_szkxky, "vflux_vs_szkxky" : vflux_vs_szkxky}
      
#===============================================================================
#                   ATTACH THE FLUXES TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_fluxesVsZ(self): 
    
    # Read the fluxes
    fluxes = read_fluxesVsZ(self.path, self.dim, self.vec)   

    # Save the fluxes
    self.phi2_vs_zkxky = Data(["phi2","z","kx","ky"], fluxes["phi2_vs_zkxky"], self.vec.z, self.vec.kx, self.vec.ky) 
    self.qflux_vs_szkxky = Data(["qflux","s","z","kx","ky"], fluxes["qflux_vs_szkxky"], range(self.dim.species), self.vec.z, self.vec.kx, self.vec.ky) 
    self.pflux_vs_szkxky = Data(["pflux","s","z","kx","ky"], fluxes["pflux_vs_szkxky"], range(self.dim.species), self.vec.z, self.vec.kx, self.vec.ky) 
    self.vflux_vs_szkxky = Data(["vflux","s","z","kx","ky"], fluxes["vflux_vs_szkxky"], range(self.dim.species), self.vec.z, self.vec.kx, self.vec.ky)
    return 
     



