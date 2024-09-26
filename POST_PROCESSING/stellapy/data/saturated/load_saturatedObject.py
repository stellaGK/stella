 
from stellapy.data.saturated.read_saturatedNetcdf import get_saturatedData 
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime 

#===============================================================================
#                        CREATE THE SATURATED OBJECT
#===============================================================================
  
class Saturated:
    
    # Copy the data from <simulation> that is needed to construct <distribution>
    def __init__(self, simulation):
        
        # Remember the file paths and the progress
        self.path = simulation.path  
        self.vec = simulation.vec
        self.dim = simulation.dim
        return  
    
    # Get the simulation data over t = [tend/2, tend]
    @calculate_attributeWhenReadFirstTime 
    def trange(self):                 get_saturatedData(self);      return self.trange 
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_kxky(self):           get_saturatedData(self);      return self.phi2_vs_kxky 
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_zkxky(self):           get_saturatedData(self);      return self.phi_vs_zkxky 
    @calculate_attributeWhenReadFirstTime 
    def pflux_vs_szkxky(self):        get_saturatedData(self);      return self.pflux_vs_szkxky 
    @calculate_attributeWhenReadFirstTime 
    def vflux_vs_szkxky(self):        get_saturatedData(self);      return self.vflux_vs_szkxky 
    @calculate_attributeWhenReadFirstTime 
    def qflux_vs_szkxky(self):        get_saturatedData(self);      return self.qflux_vs_szkxky 
    @calculate_attributeWhenReadFirstTime 
    def dens_vs_szkxky(self):         get_saturatedData(self);      return self.dens_vs_szkxky 
    @calculate_attributeWhenReadFirstTime 
    def upar_vs_szkxky(self):         get_saturatedData(self);      return self.upar_vs_szkxky 
    @calculate_attributeWhenReadFirstTime 
    def temp_vs_szkxky(self):         get_saturatedData(self);      return self.temp_vs_szkxky 
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_svpaz(self):            get_saturatedData(self);      return self.g2_vs_svpaz 
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_smuvpa(self):           get_saturatedData(self);      return self.g2_vs_smuvpa 
    
#-------------------------
def load_saturatedObject(self): 
    self.saturated = Saturated(self) 
            