 
#===============================================================================
#                    GET THE BASIC GRID PARAMETERS
#===============================================================================

def get_gridParameters(self):
    
    # Save the (z, mu, vpa) grid
    self.nzed = int(self.inputParameters["zgrid_parameters"]["nzed"])
    self.nzgrid = int(self.inputParameters["zgrid_parameters"]["nzgrid"])
    self.nvgrid = int(self.inputParameters["vpamu_grids_parameters"]["nvgrid"])
    self.nmu = int(self.inputParameters["vpamu_grids_parameters"]["nmu"])
    return
