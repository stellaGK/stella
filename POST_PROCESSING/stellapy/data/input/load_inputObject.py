 
from stellapy.data.input.get_gridParameters import get_gridParameters
from stellapy.data.input.get_modeParameters import get_modeParameters
from stellapy.data.input.get_inputParameters import get_inputParameters
from stellapy.data.input.get_basicParameters import get_basicParameters 
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime

#===============================================================================
#                        CREATE THE INPUT OBJECT
#===============================================================================

class Input:
    
    # Copy the data from <simulation> that is needed to construct <input>
    def __init__(self, simulation):
        
        # Remember the path object of <simulation> or <mode>   
        self.path = simulation.path
        return
        
    # Load all the data for the input object
    def load_data(self):
        get_inputParameters(self)
        get_basicParameters(self)
        get_gridParameters(self)
        get_modeParameters(self)
        return 
    
    # Read all input parameters
    @calculate_attributeWhenReadFirstTime
    def inputParameters(self):  get_inputParameters(self);    return self.inputParameters 
    @calculate_attributeWhenReadFirstTime
    def input_parameters(self): get_inputParameters(self);    return self.input_parameters 
    
    # Read basic input parameters
    @calculate_attributeWhenReadFirstTime
    def nspec(self):            get_basicParameters(self);    return self.nspec
    @calculate_attributeWhenReadFirstTime
    def linear(self):           get_basicParameters(self);    return self.linear
    @calculate_attributeWhenReadFirstTime
    def nonlinear(self):        get_basicParameters(self);    return self.nonlinear
    @calculate_attributeWhenReadFirstTime
    def full_flux_surface(self):get_basicParameters(self);    return self.full_flux_surface
    
    # Read basic input parameters for the equilibrium
    @calculate_attributeWhenReadFirstTime
    def y0(self):               get_basicParameters(self);    return self.y0
    @calculate_attributeWhenReadFirstTime
    def rho(self):              get_basicParameters(self);    return self.rho
    @calculate_attributeWhenReadFirstTime
    def svalue(self):           get_basicParameters(self);    return self.svalue
    @calculate_attributeWhenReadFirstTime
    def vmec(self):             get_basicParameters(self);    return self.vmec 
    @calculate_attributeWhenReadFirstTime
    def miller(self):           get_basicParameters(self);    return self.miller
    @calculate_attributeWhenReadFirstTime
    def vmec_filename(self):    get_basicParameters(self);    return self.vmec_filename
    @calculate_attributeWhenReadFirstTime
    def nperiod(self):          get_basicParameters(self);    return self.nperiod
    @calculate_attributeWhenReadFirstTime
    def nfield_periods(self):   get_basicParameters(self);    return self.nfield_periods
    
    # Read the grid parameters
    @calculate_attributeWhenReadFirstTime
    def nzed(self):             get_gridParameters(self);     return self.nzed
    @calculate_attributeWhenReadFirstTime
    def nzgrid(self):           get_gridParameters(self);     return self.nzgrid
    @calculate_attributeWhenReadFirstTime
    def nvgrid(self):           get_gridParameters(self);     return self.nvgrid
    @calculate_attributeWhenReadFirstTime
    def nmu(self):              get_gridParameters(self);     return self.nmu  
    
    # Read the nonlinear mode parameters
    @calculate_attributeWhenReadFirstTime
    def nx(self):               get_modeParameters(self);     return self.nx
    @calculate_attributeWhenReadFirstTime
    def ny(self):               get_modeParameters(self);     return self.ny
    @calculate_attributeWhenReadFirstTime
    def nakx(self):             get_modeParameters(self);     return self.nakx
    @calculate_attributeWhenReadFirstTime
    def naky(self):             get_modeParameters(self);     return self.naky  
    
#===============================================================================
#                                  LOAD INPUT                                  #
#===============================================================================

def load_inputObject(self):     
    self.input = Input(self)
    return

