 
from stellapy.data.referenceunits.read_referenceUnits import get_referenceUnits
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime
 
#===============================================================================
#                     CREATE THE REFERENCE UNITS OBJECT
#===============================================================================

class ReferenceUnits:
    
    # Copy the data from <simulation> that is needed to construct <fluxes>
    def __init__(self, simulation):
        
        # Remember the file paths and the progress
        self.prof_n = 1
        self.prof_T = 1
        self.nspec = simulation.dim.species
        self.ref_a = simulation.geometry.aref
        self.ref_B = simulation.geometry.bref 
        self.inputParameters = simulation.inputParameters
        return

    # Read the fluxes data
    @calculate_attributeWhenReadFirstTime 
    def ref(self):      get_referenceUnits(self);       return self.ref   
    
#-------------------------
def load_referenceObject(self):
    self.referenceunits = ReferenceUnits(self)
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
    print("Test the ReferenceUnits class.\n")
    for simulation in simulations: 
        for mode in simulation.modes:   
            name = "("+str(mode.kx)+", "+str(mode.ky)+")"   
            ref = mode.referenceunits.ref['flux_heat']
            print("{:<15}".format(name), ref) 