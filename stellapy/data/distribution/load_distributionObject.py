 
from stellapy.data.distribution.read_distribution3D import get_distribution3D
from stellapy.data.distribution.read_distribution4D import get_distribution4D
from stellapy.data.distribution.read_distributionVsTime import get_distributionVsTime
from stellapy.data.distribution.read_distributionVsMuOrVpaOrZ import get_distributionVsMuOrVpaOrZ
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime 

#===============================================================================
#                        CREATE THE DISTRIBUTION OBJECT
#===============================================================================
  
class Distribution:
    
    # Copy the data from <simulation> that is needed to construct <distribution>
    def __init__(self, simulation):
        
        # Remember the file paths and the progress
        self.path = simulation.path 
        self.Progress = simulation.Progress
        
        # Remember the dimensions and vectors
        self.dim = simulation.dim 
        self.vec = simulation.vec    
        return
    
    # Get the distribution versus time
    @calculate_attributeWhenReadFirstTime 
    def g_vs_ts(self):          get_distributionVsTime(self);        return self.g_vs_ts 
    
    # Get the 3D distribution data
    @calculate_attributeWhenReadFirstTime 
    def g_vs_tsz(self):         get_distribution3D(self);           return self.g_vs_tsz
    @calculate_attributeWhenReadFirstTime 
    def g_vs_tsmu(self):        get_distribution3D(self);           return self.g_vs_tsmu 
    @calculate_attributeWhenReadFirstTime 
    def g_vs_tsvpa(self):       get_distribution3D(self);           return self.g_vs_tsvpa  
    
    # Get the 4D distribution data
    @calculate_attributeWhenReadFirstTime 
    def g_vs_tsvpaz(self):      get_distribution4D(self);           return self.g_vs_tsvpaz
    @calculate_attributeWhenReadFirstTime 
    def g_vs_tsmuvpa(self):     get_distribution4D(self);           return self.g_vs_tsmuvpa  

    # Get the distribution along (vpa) and (mu) at the final time step for each species
    @calculate_attributeWhenReadFirstTime 
    def g_vs_sz(self):          get_distributionVsMuOrVpaOrZ(self);     return self.g_vs_sz
    @calculate_attributeWhenReadFirstTime 
    def g_vs_smu(self):         get_distributionVsMuOrVpaOrZ(self);     return self.g_vs_smu
    @calculate_attributeWhenReadFirstTime 
    def g_vs_svpa(self):        get_distributionVsMuOrVpaOrZ(self);     return self.g_vs_svpa
    
#-------------------------
def load_distributionObject(self): 
    self.distribution = Distribution(self) 


################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    
    from stellapy.simulations.Simulation import create_simulations
    import timeit, pathlib; start = timeit.timeit()
    import numpy as np
     
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")    
    simulations = create_simulations(folders=folder, input_files=None, ignore_resolution=True, number_variedVariables=5) 
    print("\nWe have "+str(len(simulations))+" simulations.") 
    print("Test the Distribution class.\n")
    for simulation in simulations: 
        for mode in simulation.modes:   
            name = "("+str(mode.kx)+", "+str(mode.ky)+")"   
            gmu = mode.distribution.g_vs_smu[0,:]
            g_vs_ts = mode.distribution.g_vs_ts
            vec_time = mode.distribution.vec_time_1D
            g_vs_tsz = mode.distribution.g_vs_tsz
            vec_time_3D = mode.distribution.vec_time_3D
            g_vs_tsvpaz = mode.distribution.g_vs_tsvpaz
            vec_time_3D = mode.distribution.vec_time_4D
            g_vs_sz = mode.distribution.g_vs_sz
            g_vs_smu = mode.distribution.g_vs_smu
            g_vs_svpa = mode.distribution.g_vs_svpa
            print("{:<15}".format(name), np.shape(g_vs_svpa),  g_vs_svpa[:,10]) 
            