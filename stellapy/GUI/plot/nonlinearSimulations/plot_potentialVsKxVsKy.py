
#==========================================================================
# Plot |phi(kx,ky)^2| averaged over the time range with steady state fluxes
#===========================================================================

# Load the modules
import copy
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors
from scipy.interpolate import interp2d

# Personal modules 
from stellapy.GUI.utils import DisplayProgressGUI 
from stellapy.GUI.plot.utils import create_figure 
from stellapy.utils.decorators import verbose
from stellapy.simulations.utils.get_simulations import get_simulations
from stellapy.simulations.utils.get_experiments import get_experiments   
from stellapy.calculations.calculate_inverseFourierTransform import calculate_inverseFourierTransform
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_nx, calculate_ny

@verbose
def plot_potentialVsKxVsKy(
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
            # Specify data range
                species=[0],     
                x_quantity="kx",\
                y_quantity="kx",\
                z_quantity="qflux",\
                x_range=None,\
                y_range=None,\
                t_specific=None,\
                z_specific=None,\
                c_range=None,\
            # Labels
                x_label=None,\
                y_label=None,\
                z_label=None,\
                title=None,\
            # For the GUI the figure object already exists 
                show_figure = True,\
                ax=None,\
                Progress=None,\
            # Toggles
                ordersOfMagnitude=None,\
                units="normalized",\
                interpolate=False,\
                log=False,\
                step=4):
    ''' Plot |z(x,y)^2| averaged over (t,z) or at specific (t,z). '''
 
    # Group the input data for easier passing to functions
    data = [research, experiment_id, simulation_id, species]
    
    # Update the progress bar of the GUI 
    DisplayProgressGUI(Progress, "nonlinear map", data) 
     
    # Create the figure and keep track of the axis limits and the legend labels
    ax, z_label = create_figure(ax, units, title, x_quantity, x_label, y_quantity, y_label, z_quantity, z_label, species=species)  

    # Plot the data 
    experiment = get_experiments(research, experiment_id)[0]
    simulation = get_simulations(experiment, simulation_id)[0]
                     
    # Get the data 
    vec_t   = simulation.vec_t  
    vec_z   = simulation.vec_z 
    vec_kx  = simulation.vec_kx
    vec_ky  = simulation.vec_ky
    
    # Show the data averaged over z: read |phi**2|(t,kx,ky) or fluxes(t,kx,ky,species)
    if z_specific is None:  
        if z_quantity == "gzvs":           vec_q = simulation.gvmus[:,species[0],:,:]
        if z_quantity == "gvmus":          vec_q = simulation.gvmus[:,species[0],:,:]
        if z_quantity == "phi2":           vec_q = simulation.phi2_kxky[:,:,:]
        if z_quantity == "phi2_nozonal":   vec_q = simulation.phi2_kxky[:,:,:]
        if z_quantity == "qflux":          vec_q = simulation.qflx_kxky[:,:,:,species[0]]
        if z_quantity == "vflux":          vec_q = simulation.vflx_kxky[:,:,:,species[0]]
        if z_quantity == "pflux":          vec_q = simulation.pflx_kxky[:,:,:,species[0]]
        if z_quantity in ["phi_real", "phi_imag"]:
            print("WARNING: HAVE NOT IMPLEMENED THE SPECTRUM OF REAL AND IMAG PHI YET.")
            return 
        
    # Show the data at a specific z: read |phi**2|(t,z,kx,ky) but we don't have this for fluxes
    if z_specific is not None:  
        if z_quantity in ["qflux", "vflux", "pflux"]: 
            print("WARNING: CANT LOOK AT FLUXES SPECTRUM OF SPECIFIC Z.")
            print("         LOOK AT THE FIELD LINE AVERAGE INSTEAD PLEASE.")
            return
        
        # Get |phi**2|(t,kx,ky) from |phi**2|(t,z,kx,ky)
        if z_quantity in ["phi2", "phi2_nozonal"]:
        
            # Calculate where |phi**2|(z) is maximum along z
            if z_specific=="max":
                phi2_vs_z = simulation.phi2_vs_z
                
                # Get |phi**2|(z) by averging over a time frame
                if not t_specific: 
                    t_last  = vec_t[(~np.isnan(vec_t)).sum() - 1]
                    t_start = simulation.t_range[species[0]][0]
                    t_stop  = np.nanmin([t_last, simulation.t_range[species[0]][1]])
                    phi2_vs_z = phi2_vs_z[(vec_t > t_start) & (vec_t < t_stop),:]
                    phi2_vs_z = np.nanmean(phi2_vs_z, axis=0)
                    
                # Get |phi**2|(z) by selecting a specific t  
                if t_specific:
                    index_time = np.argwhere(vec_t > t_specific)[0][0]
                    phi2_vs_z = phi2_vs_z[index_time,:]
                
                # Calculate where |phi**2|(z) is maximum along z
                z_index = np.argmax(phi2_vs_z)
                
            # If we have a value for z, get the index
            if not z_specific=="max": 
                z_index = np.argwhere(vec_z > z_specific)[0][0] 

            # Get |phi**2|(t,z,kx,ky) to calculate |phi**2|(kx,ky) at
            # the point along z where |phi**2|(z) is maximum 
            vec_q = np.abs(simulation.phi)**2   
            vec_q = vec_q[:,z_index,:,:]
            
    # Get the quantity at a specific t or as an average over a time frame
    if not t_specific: 
        t_last  = vec_t[(~np.isnan(vec_t)).sum() - 1]
        t_start = simulation.t_range[species[0]][0]
        t_stop  = np.nanmin([t_last, simulation.t_range[species[0]][1]])
        vec_q = vec_q[(vec_t > t_start) & (vec_t < t_stop),:,:]
        vec_q = np.nanmean(vec_q, axis=0)
    if t_specific:
        index_time = np.argwhere(vec_t > t_specific)[0][0]
        vec_q = vec_q[index_time,:,:]
        
    # If we want the data in real space, do an Inverse Fourier Transform
    if x_quantity=="x":  
        
        # Change the (kx,ky) vectors to (x,y) 
        nkx = len(simulation.vec.kx)
        nky = len(simulation.vec.ky)
        simulation.nx = calculate_nx(nkx)  
        simulation.ny = calculate_ny(nky)   
        Lx = 2*np.pi/(vec_kx[1]-vec_kx[0])
        Ly = 2*np.pi/(vec_ky[1]-vec_ky[0]) 
        vec_kx = np.linspace(-Lx/2,Lx/2,simulation.nx) 
        vec_ky = np.linspace(-Ly/2,Ly/2,simulation.ny)  
        
        # For direct quantities (fluxes) we can immediately take the Fourier transform 
        if not z_quantity in ["phi2", "phi2_nozonal"]:
            vec_q = calculate_inverseFourierTransform(simulation, vec_q, axis_kx=0, axis_ky=1) 
            
        # For phi2 we need to fourier transform phi and then take the square!
        if z_quantity in ["phi2", "phi2_nozonal"]:
            
            # Get phi(t,z,kx,ky)
            vec_phi_vs_tzkxky = simulation.phi
            
            # Get the quantity at a specific t or as an average over a time frame
            if not t_specific: 
                vec_phi_vs_zkxky = vec_phi_vs_tzkxky[(vec_t > t_start) & (vec_t < t_stop),:,:,:]
                vec_phi_vs_zkxky = np.nanmean(vec_phi_vs_zkxky, axis=0)
            if t_specific: 
                vec_phi_vs_zkxky = vec_phi_vs_tzkxky[index_time,:,:]
            
            # Remove the zonal modes 
            if z_quantity=="phi2_nozonal": 
                vec_phi_vs_zkxky[:,:,0] = 0
            
            # Take the inverse Fourier
            vec_phi_vs_zxy = calculate_inverseFourierTransform(simulation, vec_phi_vs_zkxky, axis_kx=1, axis_ky=2) 
            
            # Now take the square and add (nx*ny) to justify the Parseval theorem
            vec_phi2_vs_zxy = np.abs(vec_phi_vs_zxy)**2*simulation.nx*simulation.ny
            
            # Average over z or take the spectrum at a specific z
            if z_specific is None: vec_q = calculate_fieldLineAverage(simulation, vec_phi2_vs_zxy, z_axis=0)
            if z_specific is not None: vec_q = vec_phi2_vs_zxy[z_index,:,:]
            
            # Make sure the end result is real (we have +Oj components)
            vec_q = np.real(vec_q)
            
    # Get the spectrum along vpa or vmu
    if z_quantity in ["gzvs","gvmus"]:  
        vec_ky = simulation.vec_vpa 
        vec_kx = simulation.vec_mu   
            
    # Interpolate the data 
    if interpolate:  
        
        # At ky,kx=0,0 the values are nan's so replace them with data
        index_nans = np.argwhere(np.isnan(vec_q))
        for index in index_nans: 
            if ( index[0]<=len(vec_kx) or index[0]==len(vec_q[:,1])-1 ) and (not index[0]==0): 
                vec_q[index[0], index[1]] = vec_q[index[0]-1, index[1]]  
            if ( index[0]>=len(vec_kx) or index[0]==0 ) and (not index[0]==len(vec_q[:,1])-1):           
                vec_q[index[0], index[1]] = vec_q[index[0]+1, index[1]]
                
        # Interpolate the data now that nans are removed
        function = interp2d(vec_ky, vec_kx, vec_q, kind='linear')
        xnew = np.linspace(vec_ky[0], vec_ky[-1], int(len(vec_ky))*step)
        ynew = np.linspace(vec_kx[0], vec_kx[-1], int(len(vec_kx))*step)
        vec_q_interp = function(xnew,ynew)
        vec_ky_interp, vec_kx_interp = np.meshgrid(xnew, ynew)
        vec_ky_interp = vec_ky_interp[0,:]; vec_kx_interp = vec_kx_interp[:,0]
        vec_kx = vec_kx_interp; vec_ky = vec_ky_interp; vec_q = vec_q_interp
        
        # Remove again the broken tiles 
        vec_q_interp[np.repeat(np.repeat(index_nans, interpolate, axis=1), interpolate, axis=0)] = np.nan
        vec_q_interp[np.repeat(np.repeat(index_nans, interpolate, axis=1), interpolate, axis=0)] = np.nan
    
    # Remove the zonal modes 
    if z_quantity=="phi2_nozonal" and x_quantity=="kx":  
        vec_q[:,0]=0
        
    # Plot quantity(kx,ky)
    ax, cbar = create_surfacePlot(xdata=vec_kx, ydata=vec_ky, zdata=vec_q, ax=ax, cmap='jet', \
                crange=c_range, zlabel=z_label, log=log, x_quantity=x_quantity, z_quantity=z_quantity, ordersOfMagnitude=ordersOfMagnitude)

    # Show the figure
    if show_figure: plt.show()
    if True: return cbar
     

#======================
# Make a surface plot
#======================

# Plot a surface plot z(x,y) with a colormap
def create_surfacePlot(xdata=None, ydata=None, zdata=None, cmap='jet', ax=None,\
        crange=None, zlabel=None, log=False, x_quantity="kx", z_quantity="phi", \
        ordersOfMagnitude=None, center_color_map_around_zero=False):

    # Load the modules
    from matplotlib.colors import LogNorm, LinearSegmentedColormap
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # For z(x,y) the columns refer to x and the rows to y
    zdata = np.transpose(zdata)

    # Value kx=0 is plotted from kx=0 to the next kx, correct this by shifting xdata half a tile left
    xdata = xdata-(xdata[-1]-xdata[-2])/2; 
    xdata = list(xdata); xdata.append(xdata[-1]+(xdata[-1]-xdata[-2]))
    ydata = list(ydata); ydata.append(ydata[-1]+(ydata[-1]-ydata[-2]))

    # Set the labels    
    ax.axis([xdata[0], xdata[np.size(xdata)-1], ydata[0], ydata[np.size(ydata)-1]])

    # Add white to the start of the jet colormap
    if x_quantity=="kx" and z_quantity=="phi2" and log==False:
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
        
    # Map centered on zero
    if center_color_map_around_zero:
        cdict_seismic = {'red': 
                [[0.  , 0.  , 0.  ],
                [0.25,  0.  , 0.  ],
                [0.495, 1.  , 0.  ],     # White
                [0.5,   0.  , 0.  ],     # Black
                [0.505, 0.  , 1.  ],     # White
                [0.75,  1.  , 1.  ],
                [1.  ,  0.5 , 0.5 ]], 
            'green': 
                [[0.  , 0.  , 0.  ],
                [0.25, 0.  , 0.  ],
                [0.495, 1.  , 0.  ],     # White
                [0.5,   0.  , 0.  ],     # Black
                [0.505, 0.  , 1.  ],     # White
                [0.75, 0.  , 0.  ],
                [1.  , 0.  , 0.  ]], 
            'blue': 
                [[0. , 0.3 , 0.3 ],     # Darkbue
                [0.25, 1.  , 1.  ],     # Blue
                [0.495, 1.  , 0.  ],     # White
                [0.5,   0.  , 0.  ],     # Black
                [0.505, 0.  , 1.  ],     # White
                [0.75, 0.  , 0.  ],     # Red
                [1.  , 0.  , 0.  ]]}    # Darkred
        cdict_seismic = LinearSegmentedColormap('my_colormap2',cdict_seismic,256)
        cdict_seismic = copy.copy(plt.get_cmap("seismic"))
        cdict_seismic.set_bad(color='black')
        cmap = cdict_seismic

    # Show nan values in black 
    cmap.set_bad(color='black')
    
    # Put the z-axis in logaritmic scales
    if log:
        zdata[zdata < 1.E-25] = np.NaN
        crange = [ np.nanmin(abs(zdata)), np.nanmax(abs(zdata)) ]  
        if ordersOfMagnitude:
            crange[0] = 10**(np.log10(crange[1])-ordersOfMagnitude)

    # Plot the surface plotz(x,y) 
    if not center_color_map_around_zero:
        if log == True:   
            if not crange: norm = LogNorm(vmin=np.nanmin(zdata), vmax=np.nanmax(zdata))   
            if crange:     norm = LogNorm(vmin=crange[0], vmax=crange[1])  
            img  = ax.pcolormesh(xdata, ydata, np.real(zdata), cmap=cmap, norm=norm)   
        if log == False:  
            if not crange: img  = ax.pcolormesh(xdata, ydata, np.real(zdata), cmap=cmap)   
            if crange:     img  = ax.pcolormesh(xdata, ydata, np.real(zdata), cmap=cmap, vmin=crange[0], vmax=crange[1])
    if center_color_map_around_zero:
        vmax = np.nanmax(np.nanmax(zdata, axis=0))
        vmin = np.nanmin(np.nanmin(zdata, axis=0)) 
        if abs(vmin)>abs(vmax): vmax = abs(vmin)
        if abs(vmin)<abs(vmax): vmin = -abs(vmax)
        #vmin = 0 if (vmin>0) else vmin
        #vmin = -0.1 if vmin==0 else vmin 
        #vmax = -vmin if vmax<=0 else vmax
        norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax) 
        img = ax.pcolormesh(xdata, ydata, zdata, cmap=cmap, norm=norm)   
        
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
        cbar.ax.yaxis.set_offset_position('left')
        cbar.update_ticks()  

    # Set axis
    ax.set_xlim(xmin=np.min(xdata), xmax=np.max(xdata))
    ax.set_ylim(ymin=np.min(ydata), ymax=np.max(ydata)) 
    return ax, cbar

 
