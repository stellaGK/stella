 
from stellapy.data.potential.read_potential3D import get_potential3D
from stellapy.data.potential.read_potential4D import get_potential4D
from stellapy.data.potential.read_potential5D import get_potential5D
from stellapy.data.potential.read_dPhiZVsTime import get_dPhiZVsTime  
from stellapy.data.potential.read_potentialVsZ import get_potentialVsZ
from stellapy.data.potential.read_potentialVsTime import get_potentialVsTime
from stellapy.data.potential.read_trappedWeights5D import get_trappedWeights5D
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
        
        # Remember the mode
        self.ikx = simulation.ikx if hasattr(simulation, "ikx") else 0
        self.iky = simulation.iky if hasattr(simulation, "iky") else 0
        return  
    
    # Get the maximum difference of the shape of phi(z,t) compared to phi(z,tend)
    @calculate_attributeWhenReadFirstTime 
    def dphiz_vs_t(self):           get_dPhiZVsTime(self);          return self.dphiz_vs_t 
    
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


################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    
    from stellapy.simulations.Simulation import create_simulations
    import timeit, pathlib; start = timeit.timeit() 
      
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")    
    simulations = create_simulations(folders=folder, input_files=None, ignore_resolution=True, number_variedVariables=5) 
    print("\nWe have "+str(len(simulations))+" simulations.") 
    print("Test the Potential class.\n")
    for simulation in simulations:  
        for i, mode in enumerate(simulation.modes):   
            name = "("+str(mode.kx)+", "+str(mode.ky)+")"     
            phi2_vs_z = mode.potential.phi2_vs_z
            print("{:<15}".format(name),phi2_vs_z[-3:])  
            
    
if __name__ == "__main__" and False:  
    
    from stellapy.simulations.Simulation import create_simulations
    import timeit, pathlib; start = timeit.timeit()
     
    import numpy as np
    import matplotlib.pyplot as plt 
    import matplotlib.gridspec as gridspec  
    fig = plt.figure(figsize=(18,9))  
    gs1 = gridspec.GridSpec(1,1) 
    ax = plt.subplot(gs1[0,0])   
    ax.set_xlabel("$t\, v_{\mathrm{th},r}/a$") 
    ax.set_ylabel("log($\\Delta |\\hat{\\varphi}_{\\mathbf{k}}(z)|$)") 
      
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI")    
    simulations = create_simulations(folders=folder, input_files=None, ignore_resolution=True, number_variedVariables=5) 
    print("\nWe have "+str(len(simulations))+" simulations.") 
    print("Test the Potential class.\n")
    for simulation in simulations: 
        color = plt.cm.jet(np.linspace(0,1,len(simulation.modes))) #@undefinedvariable
        for i, mode in enumerate(simulation.modes):   
            name = "("+str(mode.kx)+", "+str(mode.ky)+")"    
            dphiz_vs_t = mode.potential.dphiz_vs_t 
            vec_time = mode.potential.dphiz_vs_t_time 
            print("{:<15}".format(name),dphiz_vs_t[-3:]) 
            dphiz_vs_t[dphiz_vs_t==0] = 1.e-20   
            ax.plot(vec_time, np.log10(dphiz_vs_t),label="$k_y\\rho_i=$"+str(mode.ky),color=color[i])
    
    plt.legend(labelspacing=0, handlelength=1)
    plt.show()
            