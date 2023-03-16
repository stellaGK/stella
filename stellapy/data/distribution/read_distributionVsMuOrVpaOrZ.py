""" 

#===============================================================================
#                        Read the 2D distribution data                         #
#===============================================================================
 
Read the distribution function "g" from the {*.g_vs_mu, *.g_vs_vpa, *g_vs_z} 
or the *.out.nc files. This data is only written for linear simulations.

Returns
-------
    {g2_vs_zskxky, g2_vs_muskxky, g2_vs_vpaskxky}

Hanne Thienpondt
27/01/2023

"""

#!/usr/bin/python3 
import numpy as np   
import os, sys, pathlib

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.distribution.write_txtFileForDistributionVsMuOrVpaOrZ import write_txtFileForDistributionVsMuOrVpaOrZ
from stellapy.data.input.read_inputFile import read_modeFromInputFile
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.data.utils import Data

#===============================================================================
#                        Read the 2D distribution data                         #
#===============================================================================

def read_distributionVsMuOrVpaOrZ(path, dim, vec, nakxnaky):
    
    # Make sure the {*.g_vs_mu, *.g_vs_vpa, *g_vs_z}  files exist
    if not os.path.isfile(path.g2_vs_z) and os.path.isfile(path.output_stella): 
        write_txtFileForDistributionVsMuOrVpaOrZ(path.folder)
        
    # Read the distribution files from multiple simulations 
    if path.multiple_input_files: 
        return read_distributionFromMultipleFiles(dim, vec, path.paths) 
        
    # Read the distribution from one file for <range> simulations or one-mode-per-file simulations 
    if os.path.isfile(path.g2_vs_z):  
        if nakxnaky>1: return read_distributionFromH5File(path, dim)  
        if nakxnaky==1: return read_distributionFromTxtFile(path.g2_vs_z, path.g2_vs_mu, path.g2_vs_vpa, dim.species) 

    # Critical error if we didn't find any data
    exit_reason = "The distribution data versus (vpa), (mu) and (z) does not exist:\n"
    exit_reason += "    "+str(path.g2_vs_z)+"\n"
    exit_reason += "    "+str(path.g2_vs_mu)+"\n"
    exit_reason += "    "+str(path.g2_vs_vpa)+"\n\n"
    exit_reason += "Please write the data files through: \n"
    exit_reason += "     >> write_dataFiles -s g2D \n"
    exit_program(exit_reason, read_distributionVsMuOrVpaOrZ, sys._getframe().f_lineno) 
    return 

#-------------------------------------
def read_distributionFromH5File(path, dim): 
    """ For <range> simulations we do not have g(kx,ky) but only a single g() for all modes. """
    print_warningForRangedSimulations(path)
    dummy_zskxky = np.ones((dim.z, dim.species, dim.kx, dim.ky))*np.nan
    dummy_muskxky = np.ones((dim.mu, dim.species, dim.kx, dim.ky))*np.nan
    dummy_vpaskxky = np.ones((dim.vpa, dim.species, dim.kx, dim.ky))*np.nan
    return {"g2_vs_zskxky" : dummy_zskxky, "g2_vs_muskxky" : dummy_muskxky, "g2_vs_vpaskxky" : dummy_vpaskxky} 

#-----------------------------
def read_distributionFromTxtFile(path_g2_vs_z, path_g2_vs_mu, path_g2_vs_vpa, dim_species): 
    
    # Read the {*.g_vs_mu, *.g_vs_vpa, *g_vs_z} files 
    if os.path.isfile(path_g2_vs_z): 
        g2_vs_zs = np.loadtxt(path_g2_vs_z,skiprows=1,dtype='float').reshape(-1, dim_species)
        g2_vs_mus = np.loadtxt(path_g2_vs_mu,skiprows=1,dtype='float').reshape(-1, dim_species)
        g2_vs_vpas = np.loadtxt(path_g2_vs_vpa,skiprows=1,dtype='float').reshape(-1, dim_species)   
        return {"g2_vs_zs" : g2_vs_zs, "g2_vs_mus" : g2_vs_mus, "g2_vs_vpas" : g2_vs_vpas} 
    
    # Critical error if we didn't find any data
    exit_reason = "The distribution data versus (vpa), (mu) and (z) does not exist:\n"
    exit_reason += "    "+str(path_g2_vs_z)+"\n"
    exit_reason += "    "+str(path_g2_vs_mu)+"\n"
    exit_reason += "    "+str(path_g2_vs_vpa)+"\n\n"
    exit_reason += "Please write the data files through: \n"
    exit_reason += "     >> write_dataFiles -s g2D \n"
    exit_program(exit_reason, read_distributionVsMuOrVpaOrZ, sys._getframe().f_lineno) 
    return 

#-----------------------------
def read_distributionFromMultipleFiles(dim, vec, paths): 
    
    # Create matrices (kx,ky) 
    g2_vs_zskxky = np.ones((dim.z, dim.species, dim.kx, dim.ky))*np.nan 
    g2_vs_muskxky = np.ones((dim.mu, dim.species, dim.kx, dim.ky))*np.nan 
    g2_vs_vpaskxky = np.ones((dim.vpa, dim.species, dim.kx, dim.ky))*np.nan 
    
    # Iterate over the simulations (1 mode per simulation)
    for path in paths:
        
        # Don't read the dummy simulation, we need to read each input file inside of it
        if "_dummy.in" not in str(path.input_file):
            
            # We do not have g2(kx,ky) information for range simulations
            if path.nakxnaky>1: 
                print_warningForRangedSimulations(path)
                continue
        
            # Get the mode (kx,ky)
            kx, ky = read_modeFromInputFile(path.input_file)
            ikx = list(vec.kx).index(kx)
            iky = list(vec.ky).index(ky)
            
            # Get the potential data 
            distribution_data = read_distributionFromTxtFile(path.g2_vs_z, path.g2_vs_mu, path.g2_vs_vpa, dim.species) 
                
            # Put the potential data in the matrices     
            g2_vs_zskxky[:,:,ikx,iky] = distribution_data["g2_vs_zs"]  
            g2_vs_muskxky[:,:,ikx,iky] = distribution_data["g2_vs_mus"]  
            g2_vs_vpaskxky[:,:,ikx,iky] = distribution_data["g2_vs_vpas"]  
 
    return {"g2_vs_zskxky" : g2_vs_zskxky, "g2_vs_muskxky" : g2_vs_muskxky, "g2_vs_vpaskxky" : g2_vs_vpaskxky}

#-----------------------------
def print_warningForRangedSimulations(path):
    print("\n", "".center(80,"="))
    print("", "WARNING".center(80," ")); 
    print("", "".center(80,"="))
    print("  If multiple modes are launched together in a linear simulation, ")
    print("  then we do do not have access to the g^2(kx,ky) data. Affected modes:")
    print("\n  "+str(path.input_file.parent))
    for kx in path.vec_kx:
        for ky in path.vec_ky: 
            print("         (kx,ky) = ("+str(round(kx,10))+", "+str(round(ky,10))+")")
    print()
    return

#===============================================================================
#         Attach the distribution squared to the <distribution> object         #
#===============================================================================

def get_distributionVsMuOrVpaOrZ(self): 

    # Read the data in the fluxes file
    distribution = read_distributionVsMuOrVpaOrZ(self.path, self.dim, self.vec, self.nakxnaky)  
    
    # Save the distribution data 
    self.g2_vs_szkxky = Data(["g2","s","z","kx","ky"], np.swapaxes(distribution["g2_vs_zskxky"],0,1), self.vec.species, self.vec.z, self.vec.kx, self.vec.ky)
    self.g2_vs_smukxky = Data(["g2","s","mu","kx","ky"], np.swapaxes(distribution["g2_vs_muskxky"],0,1), self.vec.species, self.vec.mu, self.vec.kx, self.vec.ky)
    self.g2_vs_svpakxky = Data(["g2","s","vpa","kx","ky"], np.swapaxes(distribution["g2_vs_vpaskxky"],0,1), self.vec.species, self.vec.vpa, self.vec.kx, self.vec.ky)
    return 
