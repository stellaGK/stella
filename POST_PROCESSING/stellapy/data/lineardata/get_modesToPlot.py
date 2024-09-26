

#===============================================================================
#                               GET MODES TO PLOT                              #
#=============================================================================== 
  
def get_modesToPlot(self, kx_range, ky_range):
    ''' Reduces the vec_kx and vec_ky to fall with kx_range and ky_range. '''
 
    # Check with the previous parameters to avoid recalculating 
    if self.kx_range != kx_range or self.ky_range != ky_range: 
        
        # Save the current range of modes
        self.kx_range = kx_range
        self.ky_range = ky_range 
        
        # Gather the modes that lie within (kx_range, ky_range)
        modes = [mode for mode in self.modes if (round(mode.ky,10) >= ky_range[0] and round(mode.ky,10) <= ky_range[1])]
        modes = [mode for mode in modes if (round(mode.kx,10) >= kx_range[0] and round(mode.kx,10) <= kx_range[1])]

        # Gather the modes that are: stable, unstable, converged and unconverged
        self.allModes = modes 
        self.stableModes = [mode for mode in modes if mode.lineardata.stable]
        self.unstableModes = [mode for mode in modes if mode.lineardata.unstable]
        self.convergedModes = [mode for mode in modes if mode.lineardata.converged]
        self.unconvergedModes = [mode for mode in modes if not mode.lineardata.converged]
        
    return self
