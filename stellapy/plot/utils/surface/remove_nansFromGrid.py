"""

#===============================================================================
#              INTERPOLATE/EXTRAPOLATE TO REMOVE NANS FROM Z(X,Y)              #
#===============================================================================

In order to interpolate 2D data such as z(x,y), it is important that there are 
no nan's in the grid z(x,y). Use <interpolate.griddata> to interpolate or 
extrapolate to remove all nan's from the grid.

Hanne Thienpondt
01/09/2022

"""

#!/usr/bin/python3
import numpy as np
from scipy import interpolate

#===============================================================================
#              INTERPOLATE/EXTRAPOLATE TO REMOVE NANS FROM Z(X,Y)              #
#===============================================================================
 
def remove_nansFromGrid(x, y, z):  
    if len(np.argwhere(np.isnan(z)))>0:
        z = np.ma.masked_invalid(z) 
        x, y = np.meshgrid(x, y)
        xx = x[~z.mask]
        yy = y[~z.mask]
        zz = z[~z.mask]
        z = interpolate.griddata((xx, yy), zz.ravel(), (x, y), method='cubic')  
    return z
