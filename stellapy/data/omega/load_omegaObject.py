""" 

#===============================================================================
#                          Create the <omega> object                           #
#===============================================================================   

Attach the data related to the *.omega file to an <omega> object.

Thanks to the "calculate_attributeWhenReadFirstTime" wrapper, the <omega> 
attributes will not be loaded/read/calculated until they have been requested.

Attributes
----------
    omega_vs_tkxky = simulation.omega.omega_vs_tkxky.omega
    gamma_vs_tkxky = simulation.omega.gamma_vs_tkxky.gamma
    vec_time = simulation.omega.gamma_vs_tkxky.t
    vec_kx = simulation.omega.gamma_vs_tkxky.kx 
    vec_ky = simulation.omega.gamma_vs_tkxky.ky 
    
    
Hanne Thienpondt
20/01/2023

"""

#!/usr/bin/python3
import os, sys
import pathlib

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime 
from stellapy.data.omega.read_omegaFile import get_omegaAndGamma

#===============================================================================
#                          Create the <omega> object                           #
#===============================================================================  

class Omega:
    
    # Copy the data from <simulation> that is needed to construct <omega>
    def __init__(self, simulation): 
        
        # Remember the file paths and the progress
        self.path = simulation.path 
        self.Progress = simulation.Progress  
        
        # Remember the sign of the magnetic field
        self.sign_B = simulation.geometry.sign_B  
        
        # For multiple modes per simulation the omega file is an h5 file 
        self.nakxnaky = simulation.nakxnaky
        self.dim = simulation.dim 
        self.vec = simulation.vec
        return
        
    # Read the omega data    
    @calculate_attributeWhenReadFirstTime 
    def omega_vs_tkxky(self):   get_omegaAndGamma(self);        return self.omega_vs_tkxky
    @calculate_attributeWhenReadFirstTime 
    def gamma_vs_tkxky(self):   get_omegaAndGamma(self);        return self.gamma_vs_tkxky 
    
#----------------------    
def load_omegaObject(self):
    self.omega = Omega(self)
    return
     
    
    
    
