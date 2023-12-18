#!/usr/bin/python3
import copy
import numpy as np
from scipy import interpolate
from scipy.interpolate import interp2d

# Stellapy package     
from stellapy.plot.utils.surface.interpolate_booleanMap import interpolate_booleanMap
from stellapy.plot.utils.surface.remove_nansFromGrid import remove_nansFromGrid

#===============================================================================
#                         INTERPOLATE A SURFACE Z(X,Y)                         #
#===============================================================================

def interpolate_surface(x, y, z, step=4, mode_maps=None, removeBoundaries=False): 

    # Recall where the data is zero or nan
    filter_zero_and_nans = (np.abs(z) < 1.E-25) | np.isnan(z)
    
    # Remove the nans from the grid
    z = remove_nansFromGrid(x, y, z)
    
    # Create new (x,y) grid with higher resolution
    xnew = np.linspace(x[0], x[-1], int(len(x)+1)*step+1)
    ynew = np.linspace(y[0], y[-1], int(len(y)+1)*step+1)  
    xnew = list(xnew) + list(x); xnew = np.array(sorted(list(set(xnew))))  
    ynew = list(ynew) + list(y); ynew = np.array(sorted(list(set(ynew))))  
    xnew = np.array([round(x,12) for x in xnew])
    ynew = np.array([round(y,12) for y in ynew])
    
    # Interpolate z(x,y) to znew(xnew,ynew) for the entire map  
    if not mode_maps:  
        function = interp2d(x, y, z.T, kind='linear')
        filter_zero_and_nans = ~interpolate_booleanMap(y, x, ynew, xnew, ~filter_zero_and_nans.T)
        return xnew, ynew, (function(xnew,ynew)).T, filter_zero_and_nans.T
        
    # Interpolate z(x,y) to znew(xnew,ynew) for each piece of the map seperatly
    if mode_maps:
        
        # Initiate the znew(xnew,ynew) to np.nan
        znew = np.ones((len(xnew),len(ynew)))*np.nan 
        
        # Iterate over the pieces of the map
        for imap in range(len(mode_maps)):  
            
            # Initiate the znew_map(xnew,ynew) to np.nan, which only has the a piece of the grid
            znew_map = np.ones((len(xnew),len(ynew)))*np.nan 
            
            # Get the map of booleans that defined which points (fprim, tiprim) belong to the mode
            # And fill <z_map> which the grid points of <z> that need to be retained according to <mode_map>
            mode_map = mode_maps[imap]
            z_map = copy.deepcopy(z)  
            z_map[~mode_map] = np.nan    
            
            # First interpolate black squares and then extrapolate (x=ky and y=a/Ln)
            z_map = np.ma.masked_invalid(z_map)  
            xmesh, ymesh = np.meshgrid(y, x)
            xx = xmesh[~z_map.mask]
            yy = ymesh[~z_map.mask]
            zz = z_map[~z_map.mask]
            try: z_map = interpolate.griddata((xx, yy), zz.ravel(), (xmesh, ymesh), method='cubic') 
            except: pass   
            
            for iy in range(len(y)):
                interpolate_filter = ~np.isnan(z_map[:,iy])
                if len(np.array(x)[interpolate_filter])!=0 and len(np.array(x)[interpolate_filter])!=1:
                    f = interpolate.interp1d(np.array(x)[interpolate_filter], z_map[:,iy][interpolate_filter], fill_value="extrapolate")
                    z_map[~interpolate_filter, iy] = f(np.array(x)[~interpolate_filter])
                if len(np.array(x)[interpolate_filter])==1:
                    z_map[~interpolate_filter, iy] = z_map[interpolate_filter, iy] 
            for ix in range(len(x)):
                interpolate_filter = ~np.isnan(z_map[ix,:])
                if len(np.array(y)[interpolate_filter])!=0 and len(np.array(y)[interpolate_filter])!=1:
                    f = interpolate.interp1d(np.array(y)[interpolate_filter], z_map[ix,:][interpolate_filter], fill_value="extrapolate")
                    z_map[ix, ~interpolate_filter] = f(np.array(y)[~interpolate_filter])
                if len(np.array(y)[interpolate_filter])==1:
                    z_map[ix, ~interpolate_filter] = z_map[ix, interpolate_filter]
 
            # Fill in the high resolution grid <znew_map> basde on the low resolution <z_map>
            for i in range(len(x)):
                for j in range(len(y)):
                    a = list(xnew).index(x[i])
                    b = list(ynew).index(y[j])
                    znew_map[a,b] = z_map[i,j]
                    
            # Mask the nan values to be invalid
            znew_map = np.ma.masked_invalid(znew_map)  
            
            # Transform (x,y) to a mesh grid and apply the mask 
            xmesh, ymesh = np.meshgrid(ynew, xnew)
            xx = xmesh[~znew_map.mask]
            yy = ymesh[~znew_map.mask]
            zz = znew_map[~znew_map.mask]
            
            # Interpolate the <znew_map> which should interpolate between the nan values 
            znew_map = interpolate.griddata((xx, yy), zz.ravel(), (xmesh, ymesh), method='cubic')  
                        
            # Only save the piece of <znew_map> where the interpolation is considered correct 
            mode_map = interpolate_booleanMap(x,y,xnew,ynew,mode_map,removeBoundaries) 
            znew[mode_map & ~np.isnan(znew_map)] = znew_map[mode_map & ~np.isnan(znew_map)]   
             
    # Get where we had zeros (stable modes)  
    filter_gamma_zero = ~interpolate_booleanMap(x,y,xnew,ynew,~filter_gamma_zero.astype(np.bool), method1=True)
    return xnew, ynew, znew, filter_gamma_zero
