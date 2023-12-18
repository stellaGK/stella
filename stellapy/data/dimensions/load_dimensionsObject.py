""" 

#===============================================================================
#                        Create the <dimensions> object                        #
#===============================================================================   

Attach the dimensions of the simulation to a <dimension> object. 

Thanks to the "calculate_attributeWhenReadFirstTime" wrapper, the <dimension> 
attributes will not be loaded/read/calculated until they have been requested.

Attributes
----------
    dim_x = simulation.dim.x
    dim_y = simulation.dim.y
    dim_z = simulation.dim.z
    dim_kx = simulation.dim.kx
    dim_ky = simulation.dim.ky
    dim_mu = simulation.dim.mu
    dim_vpa = simulation.dim.vpa
    dim_time = simulation.dim.time
    dim_species = simulation.dim.species 
    
Note that if multiple linear simulations have been combined, dim_time represents
the maximum of the time dimensions of the individual simulations.

Hanne Thienpondt
20/01/2023

"""

#!/usr/bin/python3
import os, sys
import pathlib

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime
from stellapy.data.dimensions.read_dimensionsAndVectors import get_dimensionsAndVectors as get_dimensions
from stellapy.data.dimensions.calculate_xyGrid import get_xyGrid

#===============================================================================
#                        Create the <dimensions> object                        #
#===============================================================================   
       
class Dimensions:
    
    def __init__(self, simulation): 
        
        # Copy the data from <simulation> that is needed to construct <dimensions>
        self.input = simulation.input 
        self.path = simulation.path
        self.dim = self  
        
        # Be able to show the reading progress
        self.simulation = simulation   
        self.Progress = simulation.Progress 
        return
        
    # The dimension and vector objects are filled at the same time    
    def linkToVector(self, simulation): 
        self.vec = simulation.vec
        
    # Read the dimensions when asked for
    @calculate_attributeWhenReadFirstTime 
    def z(self):            get_dimensions(self);     return self.z
    @calculate_attributeWhenReadFirstTime 
    def kx(self):           get_dimensions(self);     return self.kx
    @calculate_attributeWhenReadFirstTime 
    def ky(self):           get_dimensions(self);     return self.ky
    @calculate_attributeWhenReadFirstTime 
    def mu(self):           get_dimensions(self);     return self.mu
    @calculate_attributeWhenReadFirstTime 
    def vpa(self):          get_dimensions(self);     return self.vpa
    @calculate_attributeWhenReadFirstTime 
    def time(self):         get_dimensions(self);     return self.time
    @calculate_attributeWhenReadFirstTime 
    def species(self):      get_dimensions(self);     return self.species  
    
    # Read the real space dimensions
    @calculate_attributeWhenReadFirstTime 
    def x(self):            get_xyGrid(self);         return self.x
    @calculate_attributeWhenReadFirstTime 
    def y(self):            get_xyGrid(self);         return self.y

#-------------------------- 
def load_dimensionsObject(self):  
    self.dim = Dimensions(self)
    self.dim.linkToVector(self) 
    return 


