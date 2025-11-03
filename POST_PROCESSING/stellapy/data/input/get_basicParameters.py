 
import numpy as np 

#===============================================================================
#                    GET THE BASIC INPUT PARAMETERS
#===============================================================================

def get_basicParameters(self): 
    
    # Basic input parameters
    self.y0 = self.inputParameters['kxky_grid_box']['y0']
    self.nspec = int(self.inputParameters['species_options']['nspec'])
    
    # Linear or nonlinear simulations
    self.linear = True if (self.inputParameters['gyrokinetic_terms']['include_nonlinear']!=True) else False
    self.nonlinear = True if (self.inputParameters['gyrokinetic_terms']['include_nonlinear']==True) else False
    self.full_flux_surface = True if (self.inputParameters['gyrokinetic_terms']['include_full_flux_annulus']==True) else False
    
    # Read the number of periods or field periods
    self.nperiod = self.inputParameters['z_grid']['nperiod']
    self.nfield_periods = self.inputParameters['geometry_vmec']['nfield_periods']
    
    # Radial location in the vmec or miller file
    self.vmec_filename = self.inputParameters['geometry_vmec']['vmec_filename']
    if self.vmec_filename=='wout*.nc':
        self.vmec = False
        self.miller = True
        self.rho = self.inputParameters['geometry_miller']['rhoc']
        self.svalue = self.rho*self.rho 
    if self.vmec_filename!='wout*.nc':
        self.vmec = True
        self.miller = False
        self.rho = np.sqrt(self.inputParameters['geometry_vmec']['torflux'])
        self.svalue = self.inputParameters['geometry_vmec']['torflux']
    return
