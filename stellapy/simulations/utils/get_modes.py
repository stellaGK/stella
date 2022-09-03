from stellapy.data.lineardata.get_modesToPlot import get_modesToPlot

def get_modes(simulation, kx_range, ky_range, modes_id): 
    
    # Update the vectors of modes 
    simulation.lineardata = get_modesToPlot(simulation.lineardata, kx_range, ky_range)
    
    # Only return the modes we want to plot 
    if modes_id=="all":         return simulation.lineardata.allModes
    if modes_id=="stable":      return simulation.lineardata.stableModes
    if modes_id=="unstable":    return simulation.lineardata.unstableModes
    if modes_id=="converged":   return simulation.lineardata.convergedModes
    if modes_id=="unconverged": return simulation.lineardata.unconvergedModes
    return
