
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime
from stellapy.data.dimensions.get_dimensionsAndVectors import get_dimensionsAndVectors as get_dimensions

#===============================================================================
#                      CREATE THE DIMENSIONS OBJECT
#===============================================================================
       
class Dimensions:
    
    def __init__(self, simulation): 
        
        # Copy the data from <simulation> that is needed to construct <dimensions>
        self.input = simulation.input 
        self.path = simulation.path
        self.dim = self  
        
        # Be able to show the reading progress
        self.simulation = simulation if simulation.object=="Simulation" else simulation.simulation
        self.mode = None if simulation.object=="Simulation" else simulation
        self.Progress = simulation.Progress
        self.object = simulation.object
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
    

#-------------------------- 
def load_dimensionsObject(self): 
    
    # Add the dimensions for a simulation
    if self.object=="Simulation": 
        self.dim = Dimensions(self)
        self.dim.linkToVector(self)
        return 
    
    # For a linear simulation, add the input per mode. Since most inputs are
    # identical, only read each unique input once (self==mode)
    if self.object=="Mode":
        load_onlyReferenceDimensions(self) 
        return 

#-------------------------- 
def load_onlyReferenceDimensions(self): 
    
    # Initialize 
    loaded_list=[]
    loaded_dimensions={}

    # Only read sets of vectors/dimensions per unique input
    for imode, mode in enumerate(self.simulation.modes): 
        
        # If the input is already read, create a reference to the existing input
        if mode.path.input in loaded_list:   
            mode.dim = self.simulation.modes[loaded_dimensions[mode.path.input]].dim
        if mode.path.input not in loaded_list:
            loaded_list.append(mode.path.input) 
            loaded_dimensions[mode.path.input] = imode
            mode.dim = Dimensions(mode)  
            mode.dim.linkToVector(mode)   
            
    # Clean up the namespace
    del loaded_dimensions
    del loaded_list 
    return
    
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    
    import pathlib
    from stellapy.simulations.Simulation import create_simulations 
    folder = pathlib.Path("/home/hanne/CIEMAT/PREVIOUSRUNS/LINEARMAPS/W7Xstandard_rho0.7_aLTe0/LinearMap/fprim4tprim4") 
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")    
    simulations = create_simulations(folders=folder, ignore_resolution=True, number_variedVariables=5)
    print("We have "+str(len(simulations))+" simulations") 
    print("Test the Dimensions class.\n")
    for simulation in simulations:  
        for mode in simulation.modes:  
            name = "("+str(mode.kx)+", "+str(mode.ky)+")"  
            print("{:<15}".format(name), mode.dim.z, mode.dim.species, mode.vec.z[:3]) 


