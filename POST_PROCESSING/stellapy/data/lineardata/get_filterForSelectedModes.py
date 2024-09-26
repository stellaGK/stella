
import numpy as np

#===============================================================================
#                               GET MODES TO PLOT                              #
#=============================================================================== 
  
def get_filterForSelectedModes(simulation, modes_id="all", kx_range=[-9999,9999], ky_range=[-9999,9999]): 

    # Get either the stable or the unstable modes  
    if modes_id=="unstable": selected_modes = simulation.lineardata.unstable_vs_kxky
    if modes_id=="stable": selected_modes = simulation.lineardata.stable_vs_kxky
    if modes_id=="all": selected_modes = np.ones((simulation.dim.kx, simulation.dim.ky))
    
    # Give a warning if none of the modes got selected
    if np.all(selected_modes==0): 
        print("\n", "".center(80,"="))
        print("", "WARNING".center(80," ")); 
        print("", "".center(80,"="))
        print("  No modes (kx,ky) were selected to be plotted for: ")
        print("     ", str(simulation.input_file), "\n") 
        print("  Try changing the selected modes to 'all' by adding the flag '-m all'. ")
        print("", "".center(80,"="), "\n")
     
    # Get the modes (kx,ky) within kx_range and ky_range  
    selected_modes[:, (np.array(simulation.vec.ky)<ky_range[0]) | (np.array(simulation.vec.ky)>ky_range[1])] = 0
    selected_modes[(np.array(simulation.vec.kx)<kx_range[0]) | (np.array(simulation.vec.kx)>kx_range[1]), :] = 0 
    return selected_modes.astype(bool)
