""" 

#===============================================================================
#                          Create the <vectors> object                         #
#===============================================================================   

Attach the dimensions of the simulation to a <dimension> object.

Thanks to the "calculate_attributeWhenReadFirstTime" wrapper, the <vector> 
attributes will not be loaded/read/calculated until they have been requested. 

Attributes
----------
    vec_x = simulation.vec.x
    vec_y = simulation.vec.y
    vec_z = simulation.vec.z
    vec_kx = simulation.vec.kx
    vec_ky = simulation.vec.ky
    vec_mu = simulation.vec.mu
    vec_vpa = simulation.vec.vpa
    vec_time = simulation.vec.time
    vec_species = simulation.vec.species
    vec_kx_stella = simulation.vec.kx_stella
    vec_pol = simulation.vec.pol
    vec_tor = simulation.vec.tor
    vec_zeta = simulation.vec.zeta

Hanne Thienpondt
20/01/2023

"""  

#!/usr/bin/python3
import os, sys
import pathlib

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime 
from stellapy.data.dimensions.read_dimensionsAndVectors import get_dimensionsAndVectors as get_vectors 
from stellapy.data.dimensions.read_dimensionsAndVectors import get_extraZVectors
from stellapy.data.dimensions.calculate_xyGrid import get_xyGrid

#===============================================================================
#                          Create the <vectors> object                         #
#=============================================================================== 
 
class Vectors:
     
    def __init__(self, simulation):  
        
        # Copy the data from <simulation> that is needed to construct <dimensions>
        self.input = simulation.input 
        self.path = simulation.path
        self.vec = self  
        
        # Be able to show the reading progress
        self.simulation = simulation   
        self.Progress = simulation.Progress 
        return
    
    # The dimension and vector objects are filled at the same time    
    def linkToDimension(self, simulation): 
        self.dim = simulation.dim
         
    # Read the vectors when asked for 
    @calculate_attributeWhenReadFirstTime 
    def z(self):            get_vectors(self);     return self.z
    @calculate_attributeWhenReadFirstTime 
    def kx(self):           get_vectors(self);     return self.kx
    @calculate_attributeWhenReadFirstTime 
    def ky(self):           get_vectors(self);     return self.ky
    @calculate_attributeWhenReadFirstTime 
    def mu(self):           get_vectors(self);     return self.mu
    @calculate_attributeWhenReadFirstTime 
    def vpa(self):          get_vectors(self);     return self.vpa 
    @calculate_attributeWhenReadFirstTime 
    def species(self):      get_vectors(self);     return self.species 
    @calculate_attributeWhenReadFirstTime 
    def kx_stella(self):    get_vectors(self);     return self.kx_stella  
    
    # We can express the z-axis in different units
    @calculate_attributeWhenReadFirstTime 
    def pol(self):          get_extraZVectors(self);        return self.pol
    @calculate_attributeWhenReadFirstTime 
    def tor(self):          get_extraZVectors(self);        return self.tor
    @calculate_attributeWhenReadFirstTime 
    def zeta(self):         get_extraZVectors(self);        return self.zeta 
    
    # Read the real space dimensions
    @calculate_attributeWhenReadFirstTime 
    def x(self):            get_xyGrid(self);         return self.x
    @calculate_attributeWhenReadFirstTime 
    def y(self):            get_xyGrid(self);         return self.y
            
#-------------------------- 
def load_vectorsObject(self):   
    self.vec = Vectors(self)
    self.vec.linkToDimension(self)
    return  
     
            
            
            
