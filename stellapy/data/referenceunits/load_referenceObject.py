 
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

