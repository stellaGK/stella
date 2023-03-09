
import numpy as np
from stellapy.data.lineardata.read_linearDataPerMode import get_linearDataPerMode
from stellapy.data.lineardata.get_filterForSelectedModes import get_filterForSelectedModes
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime  

#===============================================================================
#                     CREATE THE MOSTUNSTBABLEMODE OBJECT
#===============================================================================  

class MostUnstableMode:
    
    # Copy the data from <simulation> that is needed to construct <mostunstablemode>
    def __init__(self, simulation, modes_id, kx_range, ky_range): 
        
        # Remember the data from <simulation>
        self.input_file = simulation.input_file 
        self.simulation = simulation 
        
        # Selected modes
        self.modes_id = modes_id
        self.kx_range = kx_range
        self.ky_range = ky_range
        return     
    
    # Get (stable; unstable; converged)  
    @calculate_attributeWhenReadFirstTime 
    def gamma(self):            get_mostUnstableMode(self);         return self.gamma
    @calculate_attributeWhenReadFirstTime 
    def omega(self):            get_mostUnstableMode(self);         return self.omega
    @calculate_attributeWhenReadFirstTime 
    def ky(self):               get_mostUnstableMode(self);         return self.ky
    
#----------------------    
def load_mostUnstableModeObject(self, modes_id="unstable", kx_range=[-9999,9999], ky_range=[-9999,9999]):
    
    # Check whether it's already calculated
    if hasattr(self, "mostunstablemode"):
        if self.mostunstablemode.modes_id==modes_id and self.mostunstablemode.kx_range==kx_range and self.mostunstablemode.ky_range==ky_range:
            return
        
    # Calculate the most unstable mode of the (kx,ky) modes in <simulation>
    self.mostunstablemode = MostUnstableMode(self, modes_id, kx_range, ky_range)
    return

#----------------------   
def get_mostUnstableMode(self):  
        
    # Get the modes
    vec_kx = np.array(self.simulation.vec.kx)
    vec_ky = np.array(self.simulation.vec.ky)
    
    # Get the modes which are selected (usually the unstable modes only)
    selected_modes = get_filterForSelectedModes(self.simulation, self.modes_id, self.kx_range, self.ky_range)
    
    # Attach nans if no modes were selected
    if np.sum(selected_modes)==0:
        self.kx = np.nan
        self.ky = np.nan
        self.gamma = np.nan
        self.omega = np.nan
        self.gamma_error = [np.nan, np.nan]
        self.omega_error = [np.nan, np.nan]
        return
        
    
    # Get the growth rate and remove the unwanted modes
    gamma = self.simulation.lineardata.gamma_avg_vs_kxky 
    gamma[~selected_modes] = np.nan  
    
    # Identify the most unstable mode
    index_most_unstable_mode = np.unravel_index(np.nanargmax(gamma, axis=None), gamma.shape) 
    
    # Save (gamma, omega, ky) of the most unstable mode, as well as the minimum
    # and maximum value of gamma and omega in the last 15% of the time trace 
    self.kx = vec_kx[index_most_unstable_mode[0]]
    self.ky = vec_ky[index_most_unstable_mode[1]]
    self.gamma = self.simulation.lineardata.gamma_avg_vs_kxky[index_most_unstable_mode]
    self.omega = self.simulation.lineardata.omega_avg_vs_kxky[index_most_unstable_mode]
    self.gamma_error = [self.simulation.lineardata.gamma_min_vs_kxky[index_most_unstable_mode], self.simulation.lineardata.gamma_max_vs_kxky[index_most_unstable_mode]]
    self.omega_error = [self.simulation.lineardata.omega_min_vs_kxky[index_most_unstable_mode], self.simulation.lineardata.omega_max_vs_kxky[index_most_unstable_mode]]
    return
    
 
    
    
    
    
    
