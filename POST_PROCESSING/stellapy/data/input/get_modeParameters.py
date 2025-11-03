
import sys
 
#===============================================================================
#                    GET THE BASIC MODE PARAMETERS
#===============================================================================

def get_modeParameters(self):
    
    # Save the (kx, ky) grid for a nonlinear simulation
    if self.inputParameters['gyrokinetic_terms']['include_nonlinear']==True:
        self.nx = int(self.inputParameters['kxky_grid_box']['nx'])
        self.ny = int(self.inputParameters['kxky_grid_box']['ny'])
        self.naky = int(self.inputParameters['kxky_grid_box']['naky'])
        self.nakx = int(self.inputParameters['kxky_grid_box']['nakx'])
        
    # Save the (kx, ky) mode for a linear simulation
    if self.inputParameters['gyrokinetic_terms']['include_nonlinear']==False:
        self.nx = 1
        self.ny = 1
        self.nakx = int(self.inputParameters['kxky_grid_box']['naky'])
        self.naky = int(self.inputParameters['kxky_grid_box']['nakx'])
    return
