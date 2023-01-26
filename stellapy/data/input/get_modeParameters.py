 
import sys 
 
#===============================================================================
#                    GET THE BASIC MODE PARAMETERS
#===============================================================================

def get_modeParameters(self): 
    
    # Save the (kx, ky) grid for a nonlinear simulation
    if self.inputParameters["physics_flags"]["nonlinear"]==True:
        self.nx = int(self.inputParameters["kt_grids_box_parameters"]["nx"])         
        self.ny = int(self.inputParameters["kt_grids_box_parameters"]["ny"])     
        self.naky = int(self.inputParameters["kt_grids_box_parameters"]["naky"]) 
        self.nakx = int(self.inputParameters["kt_grids_box_parameters"]["nakx"])  
        
    # Save the (kx, ky) mode for a linear simulation
    if self.inputParameters["physics_flags"]["nonlinear"]==False:
        self.nx = 1
        self.ny = 1
        self.nakx = int(self.inputParameters["kt_grids_box_parameters"]["naky"]) 
        self.naky = int(self.inputParameters["kt_grids_box_parameters"]["nakx"])            
    return
