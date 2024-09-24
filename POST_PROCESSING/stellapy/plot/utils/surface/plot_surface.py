
import copy
import numpy as np
import matplotlib.pyplot as plt 
from scipy.interpolate import interp2d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm, LinearSegmentedColormap

def plot_surface(
        ax, x, y, z, 
        zlabel=None, 
        cmap='jet with white start', \
        ordersOfMagnitude=None,
        interpolate=True,
        crange=None, 
        log=False,
        step=20):
    ''' Plot a surface plot z(x,y) with a colormap. ''' 
    
    # Interpolate the data  
    if interpolate:  
        
        # At ky,kx=0,0 the values are nan's so replace them with data
        index_nans = np.argwhere(np.isnan(z))
        for index in index_nans: 
            if ( index[0]<=len(x) or index[0]==len(z[:,1])-1 ) and (not index[0]==0): 
                z[index[0], index[1]] = z[index[0]-1, index[1]]  
            if ( index[0]>=len(x) or index[0]==0 ) and (not index[0]==len(z[:,1])-1):           
                z[index[0], index[1]] = z[index[0]+1, index[1]]
                
        # Interpolate the data now that nans are removed
        function = interp2d(y, x, z, kind='linear')
        xnew = np.linspace(y[0], y[-1], int(len(y))*step)
        ynew = np.linspace(x[0], x[-1], int(len(x))*step)
        z_interp = function(xnew,ynew)
        y_interp, x_interp = np.meshgrid(xnew, ynew)
        y_interp = y_interp[0,:]; x_interp = x_interp[:,0]
        x = x_interp; y = y_interp; z = z_interp
        
        # Remove again the broken tiles 
        z[np.repeat(np.repeat(index_nans, interpolate, axis=1), interpolate, axis=0)] = np.nan
        z[np.repeat(np.repeat(index_nans, interpolate, axis=1), interpolate, axis=0)] = np.nan

    # For z(x,y) the columns refer to x and the rows to y
    zdata = np.transpose(z)

    # Value kx=0 is plotted from kx=0 to the next kx, correct this by shifting xdata half a tile left
    x = x-(x[-1]-x[-2])/2; 
    x = list(x); x.append(x[-1]+(x[-1]-x[-2]))
    y = list(y); y.append(y[-1]+(y[-1]-y[-2]))

    # Set the axis limits
    ax.axis([x[0], x[np.size(x)-1], y[0], y[np.size(y)-1]])

    # Add white to the start of the jet colormap
    if cmap=='jet with white start':
        cdict = {'red': ((0., 1, 1),
                     (0.05, 1, 1),
                     (0.11, 0, 0),
                     (0.66, 1, 1),
                     (0.89, 1, 1),
                     (1, 0.5, 0.5)),
             'green': ((0., 1, 1),
                       (0.05, 1, 1),
                       (0.11, 0, 0),
                       (0.375, 1, 1),
                       (0.64, 1, 1),
                       (0.91, 0, 0),
                       (1, 0, 0)),
             'blue': ((0., 1, 1),
                      (0.05, 1, 1),
                      (0.11, 1, 1),
                      (0.34, 1, 1),
                      (0.65, 0, 0),
                      (1, 0, 0))}
        cmap = LinearSegmentedColormap('my_colormap',cdict,256)
        
    # Or load the standard jet map
    else:
        cmap = copy.copy(plt.get_cmap(cmap))

    # Show nan values in black 
    cmap.set_bad(color='black')
    
    # Put the z-axis in logaritmic scales
    if log:
        zdata[zdata < 1.E-25] = np.NaN
        if not crange:
            crange = [ np.nanmin(abs(zdata)), np.nanmax(abs(zdata)) ]  
            if ordersOfMagnitude:
                crange[0] = 10**(np.log10(crange[1])-ordersOfMagnitude)

    # Plot the surface plotz(x,y) 
    if log == True:   
        if not crange: norm = LogNorm(vmin=np.nanmin(zdata), vmax=np.nanmax(zdata))   
        if crange:     norm = LogNorm(vmin=crange[0], vmax=crange[1])  
        img  = ax.pcolormesh(x, y, np.real(zdata), cmap=cmap, norm=norm)   
    if log == False:  
        if not crange: img  = ax.pcolormesh(x, y, np.real(zdata), cmap=cmap)   
        if crange:     img  = ax.pcolormesh(x, y, np.real(zdata), cmap=cmap, vmin=crange[0], vmax=crange[1])

    # Add colorbar 
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    cbar = plt.colorbar(img, cax=cax) 
    cbar.set_label(zlabel)
    cbar.update_ticks()
        
    # Change the axis style
    if not log:
        cbar.formatter.set_powerlimits((0,0))
        cbar.formatter.set_scientific(True)

    # Set axis
    ax.set_xlim(xmin=np.min(x), xmax=np.max(x))
    ax.set_ylim(ymin=np.min(y), ymax=np.max(y)) 
    return ax, cbar, crange

 
