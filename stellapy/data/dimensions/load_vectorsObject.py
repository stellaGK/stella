  
from stellapy.data.dimensions.get_dimensionsAndVectors import get_extraZVectors
from stellapy.data.dimensions.get_dimensionsAndVectors import get_dimensionsAndVectors as get_vectors 
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime 

#===============================================================================
#                      CREATE THE VECTORS OBJECT
#===============================================================================
 
class Vectors:
     
    def __init__(self, simulation):  
        
        # Copy the data from <simulation> that is needed to construct <dimensions>
        self.input = simulation.input 
        self.path = simulation.path
        self.vec = self  
        
        # Be able to show the reading progress
        self.simulation = simulation if simulation.object=="Simulation" else simulation.simulation
        self.mode = None if simulation.object=="Simulation" else simulation 
        self.Progress = simulation.Progress
        self.object = simulation.object
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
            
#-------------------------- 
def load_vectorsObject(self):  
    
    # Add the dimensions for a simulation
    if self.object=="Simulation": 
        self.vec = Vectors(self)
        self.vec.linkToDimension(self)
        return 
    
    # For a linear simulation, add the input per mode. Since most inputs are
    # identical, only read each unique input once
    if self.object=="Mode":
        load_onlyReferenceVectors(self)
        return 
    
#-------------------------- 
def load_onlyReferenceVectors(self):
 
    # Initialize 
    loaded_list=[]
    loaded_vectors={}
    
    # Only read sets of vectors/dimensions per unique input
    for imode, mode in enumerate(self.simulation.modes): 
        
        # If the input is already read, create a reference to the existing input
        if mode.path.input in loaded_list: 
            mode.vec = self.simulation.modes[loaded_vectors[mode.path.input]].vec
        if mode.path.input not in loaded_list:
            loaded_list.append(mode.path.input) 
            loaded_vectors[mode.path.input] = imode
            mode.vec = Vectors(mode) 
            mode.vec.linkToDimension(mode)
            
            
            
            
