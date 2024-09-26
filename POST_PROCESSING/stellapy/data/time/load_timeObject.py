 
from stellapy.data.time.get_peakTime import get_peakTime
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime   
from stellapy.data.time.read_timeFrames import read_timeFrames
from stellapy.data.time.read_timeFrames import update_timeFrame as update_timeFrame2
from stellapy.data.time.get_saturatedFlux import get_saturatedFluxes
from stellapy.data.time.get_saturatedPotential import get_saturatedPotential
from stellapy.data.time.get_saturatedDistribution import get_saturatedDistribution

#===============================================================================
#                        CREATE THE TIMEFRAME OBJECT
#===============================================================================  
# All attached quantities depend on the chosen time frame, therefore, when we 
# change the time frame, re-initiate the entire class!
#=============================================================================== 

class Time:
    
    # Copy the data from <simulation> that is needed to construct <distribution>
    def __init__(self, simulation):  
        self.dim = simulation.dim
        self.path = simulation.path 
        self.fluxes = simulation.fluxes
        self.potential = simulation.potential
        self.distribution = simulation.distribution
        self.simulation = simulation
        return    
    
    # Update the time frame
    def update_timeFrame(self, trange):
        update_timeFrame2(self, trange)
        return
    
    # Data concerning the "timeframe" file
    @calculate_attributeWhenReadFirstTime 
    def file(self):             read_timeFrames(self);      return self.file 
    @calculate_attributeWhenReadFirstTime 
    def section(self):          read_timeFrames(self);      return self.section 
    @calculate_attributeWhenReadFirstTime 
    def timeFrame(self):        read_timeFrames(self);      return self.timeFrame 
    @calculate_attributeWhenReadFirstTime 
    def tstart(self):           read_timeFrames(self);      return self.tstart 
    @calculate_attributeWhenReadFirstTime 
    def tend(self):             read_timeFrames(self);      return self.tend      
    
    # Data calculated for the selected time frame
    @calculate_attributeWhenReadFirstTime 
    def peakTime(self):         get_peakTime(self);         return self.peakTime 
    
    # Saturated fluxes
    @calculate_attributeWhenReadFirstTime 
    def saturatedFluxes(self):  get_saturatedFluxes(self);  return self.saturatedFluxes 
    @calculate_attributeWhenReadFirstTime 
    def satFluxMinimum(self):   get_saturatedFluxes(self);  return self.satFluxMinimum 
    @calculate_attributeWhenReadFirstTime 
    def satFluxMaximum(self):   get_saturatedFluxes(self);  return self.satFluxMaximum 
    @calculate_attributeWhenReadFirstTime 
    def satFluxStdErrors(self): get_saturatedFluxes(self);  return self.satFluxStdErrors 
    
    # Saturated potential
    @calculate_attributeWhenReadFirstTime  
    def saturatedPotential(self):   get_saturatedPotential(self);   return self.saturatedPotential 
    @calculate_attributeWhenReadFirstTime 
    def satPotMinimum(self):        get_saturatedPotential(self);   return self.satPotMinimum 
    @calculate_attributeWhenReadFirstTime 
    def satPotMaximum(self):        get_saturatedPotential(self);   return self.satPotMaximum 
    @calculate_attributeWhenReadFirstTime 
    def satPotStdErrors(self):      get_saturatedPotential(self);   return self.satPotStdErrors 
    
    # Saturated distribution
    @calculate_attributeWhenReadFirstTime  
    def saturatedDistribution(self): get_saturatedDistribution(self);   return self.saturatedDistribution 
    @calculate_attributeWhenReadFirstTime 
    def satDistMinimum(self):        get_saturatedDistribution(self);   return self.satDistMinimum 
    @calculate_attributeWhenReadFirstTime 
    def satDistMaximum(self):        get_saturatedDistribution(self);   return self.satDistMaximum 
    @calculate_attributeWhenReadFirstTime 
    def satDistStdErrors(self):      get_saturatedDistribution(self);   return self.satDistStdErrors 
    
#----------------------    
def load_timeObject(self):
    self.time = Time(self)
    return
    
    
    
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    
    from stellapy.simulations.Simulation import create_simulations
    import timeit, pathlib; start = timeit.timeit()
     
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/nonlinearSimulationWithNCFile")    
    simulations = create_simulations(folders=folder, input_files=None, ignore_resolution=True, number_variedVariables=5) 
    print("\nWe have "+str(len(simulations))+" simulations.") 
    print("Test the TimeFrame class.\n")
    for simulation in simulations:   
        print(simulation.input_file.name, simulation.time.peakTime, simulation.time.saturatedFluxes["qflux"])
    
    
    
    
    
    
