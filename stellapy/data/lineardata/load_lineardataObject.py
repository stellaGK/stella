
import numpy as np
from stellapy.data.utils import Data
from stellapy.data.lineardata.read_removedModes import read_removedModes 
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime  
from stellapy.data.lineardata.read_linearDataForFFS import get_linearDataForFFS
from stellapy.data.lineardata.read_linearDataPerMode import get_linearDataPerMode

#===============================================================================
#                        CREATE THE LINEARDATA OBJECT
#===============================================================================  
class LinearData:
    
    # Copy the data from <simulation> that is needed to construct <lineardata>
    def __init__(self, simulation):
        
        # Calculate (omega; gamma) over the following last percentage of (time)
        self.percentage = 0.85
        
        # Calculate (gamma_ffs) after the first converged simulation 
        self.error_max = 0.01
        
        # Remember the file paths and the progress 
        self.full_flux_surface = simulation.full_flux_surface 
        self.path = simulation.path
        
        # We need the following data to calculate the linear data for e.g. gamma(ky) 
        self.input_file = simulation.input_file 
        self.simulation = simulation
        self.input = simulation.input
        self.path = simulation.path  
        self.dim = simulation.dim 
        self.vec = simulation.vec  
        self.omega = simulation.omega
        self.Progress = simulation.Progress
        self.input_file = simulation.input_file  
        self.sign_B = simulation.geometry.sign_B  
        self.kx_range = None
        self.ky_range = None
        return    
    
    # To the <simulation> attach the vector of removed modes
    @calculate_attributeWhenReadFirstTime 
    def removedKxModes(self):          read_removedModes(self);             return self.removedKxModes
    @calculate_attributeWhenReadFirstTime 
    def removedKyModes(self):          read_removedModes(self);             return self.removedKyModes
    
    # Get (stable; unstable; converged)  
    @calculate_attributeWhenReadFirstTime 
    def stable_vs_kxky(self):          get_linearDataPerMode(self);         return self.stable_vs_kxky
    @calculate_attributeWhenReadFirstTime 
    def unstable_vs_kxky(self):        get_linearDataPerMode(self);         return self.unstable_vs_kxky
    @calculate_attributeWhenReadFirstTime 
    def converged_vs_kxky(self):       get_linearDataPerMode(self);         return self.converged_vs_kxky
    
    # Calculate (stable; omega_last; gamma_last) when it's asked for 
    @calculate_attributeWhenReadFirstTime 
    def omega_last_vs_kxky(self):      get_linearDataPerMode(self);         return self.omega_last_vs_kxky
    @calculate_attributeWhenReadFirstTime 
    def omega_avg_vs_kxky(self):       get_linearDataPerMode(self);         return self.omega_avg_vs_kxky
    @calculate_attributeWhenReadFirstTime 
    def omega_min_vs_kxky(self):       get_linearDataPerMode(self);         return self.omega_min_vs_kxky
    @calculate_attributeWhenReadFirstTime 
    def omega_max_vs_kxky(self):       get_linearDataPerMode(self);         return self.omega_max_vs_kxky
    @calculate_attributeWhenReadFirstTime 
    def gamma_last_vs_kxky(self):      get_linearDataPerMode(self);         return self.gamma_last_vs_kxky
    @calculate_attributeWhenReadFirstTime 
    def gamma_avg_vs_kxky(self):       get_linearDataPerMode(self);         return self.gamma_avg_vs_kxky
    @calculate_attributeWhenReadFirstTime 
    def gamma_min_vs_kxky(self):       get_linearDataPerMode(self);         return self.gamma_min_vs_kxky
    @calculate_attributeWhenReadFirstTime 
    def gamma_max_vs_kxky(self):       get_linearDataPerMode(self);         return self.gamma_max_vs_kxky
    
    # Calculate (gamma, omega) for a full flux surface simulation
    @calculate_attributeWhenReadFirstTime 
    def ffs_gamma(self):               get_linearDataForFFS(self);          return self.ffs_gamma
    @calculate_attributeWhenReadFirstTime 
    def ffs_omega(self):               get_linearDataForFFS(self);          return self.ffs_omega
    @calculate_attributeWhenReadFirstTime 
    def ffs_kymin(self):               get_linearDataForFFS(self);          return self.ffs_kymin
    @calculate_attributeWhenReadFirstTime 
    def ffs_gamma_min(self):           get_linearDataForFFS(self);          return self.ffs_gamma_min
    @calculate_attributeWhenReadFirstTime 
    def ffs_gamma_max(self):           get_linearDataForFFS(self);          return self.ffs_gamma_max
    @calculate_attributeWhenReadFirstTime 
    def ffs_omega_min(self):           get_linearDataForFFS(self);          return self.ffs_omega_min
    @calculate_attributeWhenReadFirstTime 
    def ffs_omega_max(self):           get_linearDataForFFS(self);          return self.ffs_omega_max
    
    #----------------------   
    def get_linearDataPerSimulation(self, modes, plot):
    
        # Get the quantity on the x-axis
        x = [getattr(mode, plot.xdim) for mode in modes]  
         
        # Get the quantity on the y-axis 
        y = [(getattr(mode.lineardata, plot.ydim) if mode.lineardata.unstable else 0) for mode in modes]
 
        # Divide by k**2 if we want too
        if "/ky**2" in plot.yname: y = [y[i]/(x[i]**2) for i in range(len(y))]
        
        # Save the y-data under the correct name: e.g. "gamma_avg_vs_ky"
        dim0, dim1 = plot.quantity.split("_vs_")[0], plot.quantity.split("_vs_")[-1]
        setattr(self, plot.quantity, Data([dim0, dim1], y, x))   
        
        # Get the error
        if plot.ydim=="kx" or plot.ydim=="ky":
            minimum = np.ones((len(modes)))*np.NaN
            maximum = np.ones((len(modes)))*np.NaN
        if plot.ydim!="kx" and plot.ydim!="ky":
            minimum = [getattr(mode.lineardata, plot.ydim.split("_")[0]+"_min") if mode.lineardata.unstable else 0 for mode in modes]
            maximum = [getattr(mode.lineardata, plot.ydim.split("_")[0]+"_max") if mode.lineardata.unstable else 0 for mode in modes]
            if "/ky**2" in plot.yname: minimum = [minimum[i]/x[i]**2 for i in range(len(y))]
            if "/ky**2" in plot.yname: maximum = [maximum[i]/x[i]**2 for i in range(len(y))]
        setattr(self, plot.quantity+"_error", Data([dim0, dim1], [minimum, maximum] , x))      
        return

#----------------------    
def load_linearDataObject(self):
    self.lineardata = LinearData(self)
    return
    
    
    
    
