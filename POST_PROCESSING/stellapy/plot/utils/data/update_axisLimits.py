#!/usr/bin/python3
import numpy as np

#===============================================================================
#                                UPDATE AXIS LIMITS                            #
#===============================================================================

def update_axisLimits(x, y, xlims, ylims, percentage=None, logx=False, logy=False, overshoot_x=1, overshoot_y=1): 
        
    # Grab the arrays
    x = np.array(x)
    y = np.array(y)   
    
    # Grab the limit of the y-axis based on the last part 
    if percentage!=None:
        y = y[int(np.size(y)*percentage):]  
    
    # Calculate the limits
    if not np.isnan(y).all(): 
        
        # Don't take into account zeros for logaritmic scales
        if logx==True: x = abs(x[x != 0])
        if logy==True: y = abs(y[y != 0])
        
        # Get the x- and y-limits
        xmin = overshoot_x*np.nanmin(x) if np.nanmin(x)<0 else np.nanmin(x)
        ymin = overshoot_y*np.nanmin(y) if np.nanmin(y)<0 else np.nanmin(y)
        xmax = overshoot_x*np.nanmax(x) if np.nanmax(x)>0 else np.nanmax(x)
        ymax = overshoot_y*np.nanmax(y) if np.nanmax(y)>0 else np.nanmax(y)
        
        # Remember the limits
        xlims[0] = np.nanmin([xmin, xlims[0]])
        xlims[1] = np.nanmax([xmax, xlims[1]]) 
        ylims[0] = np.nanmin([ymin, ylims[0]])
        ylims[1] = np.nanmax([ymax, ylims[1]]) 
        
    return xlims, ylims



