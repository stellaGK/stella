
#!/usr/bin/python3  
import sys, os   
import pathlib
import numpy as np

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.data.fluxes.write_txtFileForQuasiLinearFluxes import write_txtFileForQuasiLinearFluxes 
from stellapy.utils.decorators.exit_program import exit_program   
from stellapy.data.utils import Data

#===============================================================================
#                        ATTACH THE FLUXES DATA
#===============================================================================

def read_quasilinearfluxes(path, dim, vec): 
    ''' Read the fluxes from *.QLfluxes file. ''' 
               
    # Make sure that the *.QLfluxes file exists
    if not os.path.isfile(path.QLfluxes) and os.path.isfile(path.output_stella): 
        write_txtFileForQuasiLinearFluxes(path.folder)
        
    # Read the fluxes files from multiple simulations 
    if path.multiple_input_files: 
        return read_QLfluxesFromMultipleFiles(dim, vec, path.paths) 
                              
    # Read the txt file for one simulation
    if os.path.isfile(path.QLfluxes):
        return read_QLfluxesFromTxtFile(path.QLfluxes, dim)    
    
    # Critical error if we didn't find any data 
    exit_reason = "The quasi-linear fluxes data could not be found for:\n"
    exit_reason += "     "+str(path.input_file)
    exit_program(exit_reason, read_quasilinearfluxes, sys._getframe().f_lineno)   
    return

#-----------------------------
def read_QLfluxesFromTxtFile(path_QLfluxes, dim):  
    data = np.loadtxt(path_QLfluxes,skiprows=1,dtype='float').reshape(1, 9)[0,:]  
    data_dict = {
        "kx" : data[0], "ky" : data[1], 
        "phi2_vs_kxky" : np.ones((1, 1, 1))*data[2],\
        "pflux_vs_skxky" : np.ones((dim.species, 1, 1))*np.nan,\
        "qflux_vs_skxky" : np.ones((dim.species, 1, 1))*np.nan,\
        "vflux_vs_skxky" : np.ones((dim.species, 1, 1))*np.nan}
    for i in range(dim.species):
        data_dict["pflux_vs_skxky"][i,:,:] = data[0*dim.species+3+i]
        data_dict["vflux_vs_skxky"][i,:,:] = data[1*dim.species+3+i]
        data_dict["qflux_vs_skxky"][i,:,:] = data[2*dim.species+3+i] 
    return data_dict

#-----------------------------
def read_QLfluxesFromMultipleFiles(dim, vec, paths): 
        
    # Create matrices (s,kx,ky)
    phi2_vs_kxky = np.ones((dim.kx, dim.ky))*np.nan
    pflux_vs_skxky = np.ones((dim.species, dim.kx, dim.ky))*np.nan
    qflux_vs_skxky = np.ones((dim.species, dim.kx, dim.ky))*np.nan
    vflux_vs_skxky = np.ones((dim.species, dim.kx, dim.ky))*np.nan 
    
    # Iterate over the simulations (1 mode per simulation)
    for path in paths:
        
        # Don't read the dummy simulation, we need to read each input file inside of it
        if "_dummy.in" not in str(path.input_file):
        
            # Read QLfluxes(t) from a txt file
            data = read_QLfluxesFromTxtFile(path.QLfluxes, dim)

            # Iterate over the modes (kx,ky) 
            for iikx, kx in enumerate(path.vec_kx):
                for iiky, ky in enumerate(path.vec_ky): 
                    
                    # Get the mode (kx,ky)   
                    ikx = list(vec.kx).index(kx) 
                    iky = list(vec.ky).index(ky)
                        
                    # Put the quasi-linear fluxes data in the matrices
                    phi2_vs_kxky[ikx,iky] = data["phi2_vs_kxky"][iikx,iiky]  
                    pflux_vs_skxky[:,ikx,iky] = data["pflux_vs_skxky"][:,iikx,iiky]  
                    qflux_vs_skxky[:,ikx,iky] = data["qflux_vs_skxky"][:,iikx,iiky]  
                    vflux_vs_skxky[:,ikx,iky] = data["vflux_vs_skxky"][:,iikx,iiky]  
 
    return {"phi2_vs_kxky" : phi2_vs_kxky, "pflux_vs_skxky" : pflux_vs_skxky, "qflux_vs_skxky" : qflux_vs_skxky, "vflux_vs_skxky" : vflux_vs_skxky}
      
#===============================================================================
#                   ATTACH THE FLUXES TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_quasilinearfluxes(self): 
    
    # Read the fluxes
    fluxes = read_quasilinearfluxes(self.path, self.dim, self.vec)   

    # Save the fluxes
    self.phi2_vs_kxky = Data(["phi2","kx","ky"], fluxes["phi2_vs_kxky"], self.vec.kx, self.vec.ky) 
    self.qflux_vs_skxky = Data(["qflux","s","kx","ky"], fluxes["qflux_vs_skxky"], range(self.dim.species), self.vec.kx, self.vec.ky) 
    self.pflux_vs_skxky = Data(["pflux","s","kx","ky"], fluxes["pflux_vs_skxky"], range(self.dim.species), self.vec.kx, self.vec.ky) 
    self.vflux_vs_skxky = Data(["vflux","s","kx","ky"], fluxes["vflux_vs_skxky"], range(self.dim.species), self.vec.kx, self.vec.ky)
    return 
     



