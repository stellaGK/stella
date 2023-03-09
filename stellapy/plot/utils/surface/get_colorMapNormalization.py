"""

#===============================================================================
#                       PERSONALIZED COLOR MAP NORMALIZATION                   #
#===============================================================================

In order to plot a surface we generally need {vmin, vmax, norm} which sets the 
color scale of the color map. This script sets {vmin, vmax, norm}. If <z_range> 
is not None, it overwrites the range of the actual data z(x,y).

Hanne Thienpondt
01/09/2022

"""

#!/usr/bin/python3
import numpy as np
import matplotlib.colors as mcolors 

#===============================================================================
#                       PERSONALIZED COLOR MAP NORMALIZATION                   #
#===============================================================================

def get_colorMapNormalization(z, z_range=None, start_at_zero=True, center_at_zero=False): 
         
    # Get the range of the data
    zmin = np.nanmin(np.nanmin(z, axis=0)) if z_range==None else z_range[0]
    zmax = np.nanmax(np.nanmax(z, axis=0)) if z_range==None else z_range[1]
    
    # Move the minimum to zero
    if start_at_zero and zmin>0: zmin = 0 
    
    # Create the normalization for the color scale
    norm = mcolors.TwoSlopeNorm(vmin=zmin, vcenter=0, vmax=zmax) if center_at_zero else None
    return zmin, zmax, norm
