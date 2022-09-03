 
import numpy as np 

#===============================================================================
#                    GET THE BASIC INPUT PARAMETERS
#===============================================================================

def get_basicParameters(self): 
    
    # Basic input parameters
    self.y0 = self.inputParameters["kt_grids_box_parameters"]["y0"]
    self.nspec = int(self.inputParameters["species_knobs"]["nspec"])  
    
    # Linear or nonlinear simulations
    self.linear = True if (self.inputParameters["physics_flags"]["nonlinear"]!=True) else False
    self.nonlinear = True if (self.inputParameters["physics_flags"]["nonlinear"]==True) else False
    
    # Read the number of periods or field periods
    self.nperiod = self.inputParameters['zgrid_parameters']['nperiod']
    self.nfield_periods = self.inputParameters['vmec_parameters']['nfield_periods'] 
    
    # Radial location in the vmec or miller file
    self.vmec_filename = self.inputParameters["vmec_parameters"]["vmec_filename"]
    if self.vmec_filename=='wout*.nc':
        self.vmec = False
        self.miller = True
        self.rho = self.inputParameters["millergeo_parameters"]["rhoc"]
        self.svalue = self.rho*self.rho 
    if self.vmec_filename!='wout*.nc':
        self.vmec = True
        self.miller = False
        self.rho = np.sqrt(self.inputParameters["vmec_parameters"]["torflux"])  
        self.svalue = self.inputParameters["vmec_parameters"]["torflux"]     
    return