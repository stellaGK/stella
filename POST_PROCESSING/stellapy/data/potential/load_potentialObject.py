 
#!/usr/bin/python3   
import sys, os
import pathlib

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)    
from stellapy.data.potential.read_potential3D import get_potential3D
from stellapy.data.potential.read_potential4D import get_potential4D
from stellapy.data.potential.read_potential5D import get_potential5D
from stellapy.data.potential.read_dPhiZVsTime import get_dPhiZVsTime  
from stellapy.data.potential.read_potentialVsZ import get_potentialVsZ
from stellapy.data.potential.read_potentialVsTime import get_potentialVsTime
from stellapy.data.potential.read_trappedWeights5D import get_trappedWeights5D
from stellapy.data.potential.read_potential4DLinear import get_potential4DLinear
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime 

#===============================================================================
#                        CREATE THE DISTRIBUTION OBJECT
#===============================================================================
  
class Potential:
    
    # Copy the data from <simulation> that is needed to construct <distribution>
    def __init__(self, simulation):
        
        # Remember the file paths and the progress
        self.path = simulation.path  
        self.vec = simulation.vec
        self.dim = simulation.dim
        self.linear = simulation.linear
        self.nonlinear = simulation.nonlinear
        self.nakxnaky = simulation.nakxnaky
        
        # Remember the mode
        self.ikx = simulation.ikx if hasattr(simulation, "ikx") else 0
        self.iky = simulation.iky if hasattr(simulation, "iky") else 0
        return  
    
    # Get the maximum difference of the shape of phi(z,t) compared to phi(z,tend)
    @calculate_attributeWhenReadFirstTime 
    def dphiz_vs_tkxky(self):       get_dPhiZVsTime(self);          return self.dphiz_vs_tkxky 
    
    # Get the final fields phi(z)
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_z(self):             get_potentialVsZ(self);         return self.phi_vs_z
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_z(self):            get_potentialVsZ(self);         return self.phi2_vs_z
    
    # Get the potential versus time 
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_t(self):             get_potentialVsTime(self);      return self.phi_vs_t
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_t(self):            get_potentialVsTime(self);      return self.phi2_vs_t
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_t_zonal(self):      get_potentialVsTime(self);      return self.phi2_vs_t_zonal
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_t_nozonal(self):    get_potentialVsTime(self);      return self.phi2_vs_t_nozonal
 
    # Get the 3D potential data 
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_tz(self):            get_potential3D(self);          return self.phi_vs_tz
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_tz(self):           get_potential3D(self);          return self.phi2_vs_tz
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_tz_zonal(self):      get_potential3D(self);          return self.phi_vs_tz_zonal
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_tz_nozonal(self):    get_potential3D(self);          return self.phi_vs_tz_nozonal
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_tz_zonal(self):     get_potential3D(self);          return self.phi2_vs_tz_zonal
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_tz_nozonal(self):   get_potential3D(self);          return self.phi2_vs_tz_nozonal
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_tkx(self):          get_potential3D(self);          return self.phi2_vs_tkx
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_tky(self):          get_potential3D(self);          return self.phi2_vs_tky
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_tkx(self):           get_potential3D(self);          return self.phi_vs_tkx
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_tky(self):           get_potential3D(self);          return self.phi_vs_tky
    
    # Get the 4D potential data 
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_tkxky(self):         get_potential4D(self);          return self.phi_vs_tkxky
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_tkxky(self):        get_potential4D(self);          return self.phi2_vs_tkxky
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_tkxky_zeta0(self):   get_potential4D(self);          return self.phi_vs_tkxky_zeta0
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_zkxky(self):         get_potential4DLinear(self);    return self.phi_vs_zkxky
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_zkxky(self):        get_potential4DLinear(self);    return self.phi2_vs_zkxky
    
    # Get the 5D potential data 
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_tzkxky(self):        get_potential5D(self);          return self.phi_vs_tzkxky
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_tzkxky_tavg(self):   get_potential5D(self);          return self.phi_vs_tzkxky_tavg
    
    # Get the trapped particle weights
    @calculate_attributeWhenReadFirstTime 
    def trappedw_vs_tszkxky(self):  get_trappedWeights5D(self);     return self.trappedw_vs_tszkxky
    
#-------------------------
def load_potentialObject(self): 
    self.potential = Potential(self) 

            