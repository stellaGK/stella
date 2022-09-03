
import numpy as np
from stellapy.data.utils import Data
from stellapy.data.lineardata.read_removedModes import get_removedModes
from stellapy.data.lineardata.read_linearDataPerMode import get_linearDataPerMode
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime  

#===============================================================================
#                        CREATE THE LINEARDATA OBJECT
#===============================================================================  

class LinearData:
    
    # Copy the data from <simulation> that is needed to construct <distribution>
    def __init__(self, mode):
        
        # Calculate (omega; gamma) over the following last percentage of (time)
        self.percentage = 0.85
        
        # Remember the file paths and the progress
        self.object = mode.object 
        self.path = mode.path
        
        # To get the removed modes, attach them to the simulation
        if self.object=="Simulation":
            self.modes = mode.modes 
            self.path = mode.path 
            self.kx_range = None
            self.ky_range = None
            
        # To decide whether a mode needs to be removed, access the simulation
        if self.object=="Mode":
            self.kx = mode.kx
            self.ky = mode.ky 
            self.dim = mode.dim 
            self.vec = mode.vec  
            self.omega = mode.omega
            self.Progress = mode.Progress
            self.input_file = mode.input_file 
            self.simulation = mode.simulation
            self.sign_B = mode.geometry.sign_B 
        return    
    
    # To the <simulation> attach the vector of removed modes
    @calculate_attributeWhenReadFirstTime 
    def removedKxModes(self):  get_removedModes(self);              return self.removedKxModes
    @calculate_attributeWhenReadFirstTime 
    def removedKyModes(self):  get_removedModes(self);              return self.removedKyModes
    @calculate_attributeWhenReadFirstTime 
    def removeMode(self):      get_removedModes(self);              return self.removeMode 
    
    # Calculate (stable; omega_last; gamma_last) when it's asked for 
    @calculate_attributeWhenReadFirstTime 
    def stable(self):          get_linearDataPerMode(self);         return self.stable
    @calculate_attributeWhenReadFirstTime 
    def unstable(self):        get_linearDataPerMode(self);         return self.unstable
    @calculate_attributeWhenReadFirstTime 
    def converged(self):       get_linearDataPerMode(self);         return self.converged
    @calculate_attributeWhenReadFirstTime 
    def omega_last(self):      get_linearDataPerMode(self);         return self.omega_last
    @calculate_attributeWhenReadFirstTime 
    def omega_avg(self):       get_linearDataPerMode(self);         return self.omega_avg
    @calculate_attributeWhenReadFirstTime 
    def omega_min(self):       get_linearDataPerMode(self);         return self.omega_min
    @calculate_attributeWhenReadFirstTime 
    def omega_max(self):       get_linearDataPerMode(self);         return self.omega_max
    @calculate_attributeWhenReadFirstTime 
    def gamma_last(self):      get_linearDataPerMode(self);         return self.gamma_last
    @calculate_attributeWhenReadFirstTime 
    def gamma_avg(self):       get_linearDataPerMode(self);         return self.gamma_avg
    @calculate_attributeWhenReadFirstTime 
    def gamma_min(self):       get_linearDataPerMode(self);         return self.gamma_min
    @calculate_attributeWhenReadFirstTime 
    def gamma_max(self):       get_linearDataPerMode(self);         return self.gamma_max
    
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
    
    
    
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    
    from stellapy.simulations.Simulation import create_simulations
    import timeit, pathlib; start = timeit.timeit()
    
    folder = pathlib.Path("/home/hanne/CIEMAT/PREVIOUSRUNS/LINEARMAPS/W7Xstandard_rho0.7_aLTe0/LinearMap/fprim4tprim4")   
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")    
    simulations = create_simulations(folders=folder, input_files=None, ignore_resolution=True, number_variedVariables=5) 
    print("\nWe have "+str(len(simulations))+" simulations.") 
    print("Test the LinearData class.\n")
    for simulation in simulations: 
        for mode in simulation.modes:   
            name = "("+str(mode.kx)+", "+str(mode.ky)+")"   
            gamma = mode.lineardata.gamma_avg
            unstable = mode.lineardata.unstable
            removeMode = mode.lineardata.removeMode
            print("{:<15}".format(name), gamma, unstable, removeMode) 
    
    
    
    
    
    
    
