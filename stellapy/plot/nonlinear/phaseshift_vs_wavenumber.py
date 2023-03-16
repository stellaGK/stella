#!/usr/bin/python3 
import copy, math
import sys, os 
import pathlib
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec 
from scipy.interpolate import interp2d, interp1d
from matplotlib.colors import LinearSegmentedColormap  

# Personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)     
from stellapy.simulations.Research import create_research  
from stellapy.plot.utils.labels import standardLabels
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                        PLOT FFT_AMPLITUDE(KY, OMEGA)                         #
#=============================================================================== 

def plot_phaseshift_vs_wavenumber(
        folder,
        # Quantities to be plotted
        y1_quantity="n0",
        y2_quantity="phi", 
        # Time range over which to plot the phase shift distribution
        tstart=None,    
        # Other 
        step_angle=5,
        normalize_to_one=True,  
        highlight_dominant_angle=True,
        show=True): 
           
    # Create <simulations> based on the given <folder>
    research = create_research(folders=folder) 
    simulations = [s for experiment in research.experiments for s in experiment.simulations] 
    
    # Iterate over the simulations
    for simulation in simulations: 
     
        # Create a figure
        fig = plt.figure(figsize=(18,6))
        fig.grid_specifications = gridspec.GridSpec(1, 2)
        fig.grid_specifications.update(top=0.9, left=0.1, right=0.9, bottom=0.15, wspace=0.25, hspace=0.25)
        ax_vs_kx = plt.subplot(fig.grid_specifications[0])  
        ax_vs_ky = plt.subplot(fig.grid_specifications[1])  
        
        # Plot the phase shifts 
        subplot_phaseshift_vs_wavenumber(ax_vs_kx, ax_vs_ky, simulation, y1_quantity, y2_quantity,  
                tstart=tstart, step_angle=step_angle, normalize_to_one=normalize_to_one, highlight_dominant_angle=highlight_dominant_angle)
        
    # Show the figure   
    if show: plt.show() 
    return 
 
#---------------------------------------------
def subplot_phaseshift_vs_wavenumber(ax_vs_kx, ax_vs_ky, simulation, y1_quantity="n0", y2_quantity="phi", 
        interpolate=True, tstart=None, step_angle=5, normalize_to_one=True, highlight_dominant_angle=True):
    
    # Sum every 5 degrees together to reduce flakyness of plot
    if step_angle==1: angles = np.linspace(-180, 180, 361)
    if step_angle==5: angles = np.linspace(-180, 180, 72)
    if step_angle==10: angles = np.linspace(-180, 180, 36)
    if step_angle==20: angles = np.linspace(-180, 180, 18)
    #angles = angles + np.abs(angles[1]-angles[0])/2
    
    # Get the kx and ky vectors
    vec_kx = np.array(simulation.vec.kx)
    vec_ky = np.array(simulation.vec.ky)
    
    # Get the time index where the saturated phase starts
    if tstart==None: tstart = simulation.time.tstart
    tstarts = simulation.moments.tstarts
    itstart = np.where(tstarts>=tstart)[0][0]
    
    # Read the phase shifts over the saturated phase
    variable_name1kx = "phase_shifts_"+y1_quantity+"_and_"+y2_quantity+"_vs_tanglekx"
    variable_name2kx = "phase_shifts_"+y2_quantity+"_and_"+y1_quantity+"_vs_tanglekx"
    variable_name1ky = "phase_shifts_"+y1_quantity+"_and_"+y2_quantity+"_vs_tangleky"
    variable_name2ky = "phase_shifts_"+y2_quantity+"_and_"+y1_quantity+"_vs_tangleky"
    if hasattr(simulation.moments, variable_name1kx): variable_namekx = variable_name1kx; factor = 1
    if hasattr(simulation.moments, variable_name2kx): variable_namekx = variable_name2kx; factor = -1
    if hasattr(simulation.moments, variable_name1ky): variable_nameky = variable_name1ky; factor = 1
    if hasattr(simulation.moments, variable_name2ky): variable_nameky = variable_name2ky; factor = -1
    phaseshift_vs_tanglekx = getattr(simulation.moments, variable_namekx).phaseshift[itstart:,:,:] 
    phaseshift_vs_tangleky = getattr(simulation.moments, variable_nameky).phaseshift[itstart:,:,:] 
    weightskx_vs_t = getattr(simulation.moments, variable_namekx.split("shifts_")[-1].split("_vs_")[0]+"_weights")[itstart:] 
    weightsky_vs_t = getattr(simulation.moments, variable_nameky.split("shifts_")[-1].split("_vs_")[0]+"_weights")[itstart:] 
    
    # If no phaseshifts were found for a certain angle, it is set to nan
    phaseshift_vs_tanglekx[np.isnan(phaseshift_vs_tanglekx)] = 0
    phaseshift_vs_tangleky[np.isnan(phaseshift_vs_tangleky)] = 0

    # Remove the weight of each time-piece, sum the pieces, and divide by the total weight
    phaseshift_vs_anglekx = np.sum(phaseshift_vs_tanglekx*weightskx_vs_t[:,np.newaxis,np.newaxis]/np.sum(weightskx_vs_t), axis=0)
    phaseshift_vs_angleky = np.sum(phaseshift_vs_tangleky*weightsky_vs_t[:,np.newaxis,np.newaxis]/np.sum(weightsky_vs_t), axis=0)
    
    # Sum every x degrees together
    if step_angle!=1:
        phaseshift_vs_anglekx = np.array([np.sum(phaseshift_vs_anglekx[i*step_angle:(i+1)*step_angle,:], axis=0) for i in range(len(angles))])
        phaseshift_vs_angleky = np.array([np.sum(phaseshift_vs_angleky[i*step_angle:(i+1)*step_angle,:], axis=0) for i in range(len(angles))]) 

    # Along kx, sum -kx with +kx 
    ikx_zero = np.where(np.array(simulation.vec.kx)==0.0)[0][0]
    phaseshift_vs_anglekx = phaseshift_vs_anglekx[:,ikx_zero:] + phaseshift_vs_anglekx[:,:ikx_zero+1][:,::-1]
    phaseshift_vs_anglekx[:,0] = phaseshift_vs_anglekx[:,0]/2
    vec_kx = vec_kx[ikx_zero:]
    
    # Normalise so that the maximum weight/count is 1 
    if not normalize_to_one:
        phaseshift_vs_anglekx = phaseshift_vs_anglekx/np.max(np.abs(phaseshift_vs_anglekx),axis=(0,1)) 
        phaseshift_vs_angleky = phaseshift_vs_angleky/np.max(np.abs(phaseshift_vs_angleky),axis=(0,1)) 
    if normalize_to_one: 
        phaseshift_vs_anglekx = phaseshift_vs_anglekx/np.max(np.abs(phaseshift_vs_anglekx),axis=0)[np.newaxis,:]
        phaseshift_vs_angleky = phaseshift_vs_angleky/np.max(np.abs(phaseshift_vs_angleky),axis=0)[np.newaxis,:]
    
    # Highlight the dominant angle
    if highlight_dominant_angle:
        if ax_vs_kx: dominantangle_vs_kx = np.array([factor*angles[i] for i in np.nanargmax(phaseshift_vs_anglekx, axis=0)])
        if ax_vs_ky: dominantangle_vs_ky = np.array([factor*angles[i] for i in np.nanargmax(phaseshift_vs_angleky, axis=0)])

    # Interpolate and plot weight(angles, kx) or weight(angles, ky)
    if interpolate:
        if ax_vs_kx: angles_new, vec_kx_new, phaseshift_vs_anglekx = interpolate_data(angles, vec_kx, phaseshift_vs_anglekx)
        if ax_vs_ky: angles_new, vec_ky_new, phaseshift_vs_angleky = interpolate_data(angles, vec_ky, phaseshift_vs_angleky)
    if not interpolate:
        if ax_vs_kx: angles_new, vec_kx_new, phaseshift_vs_anglekx = angles, vec_kx, phaseshift_vs_anglekx
        if ax_vs_ky: angles_new, vec_ky_new, phaseshift_vs_angleky = angles, vec_ky, phaseshift_vs_angleky
        
    # Plot the phaseshift
    if ax_vs_kx: create_surface_plot(ax_vs_kx, xdata=vec_kx_new, ydata=factor*angles_new, zdata=phaseshift_vs_anglekx.T)
    if ax_vs_ky: create_surface_plot(ax_vs_ky, xdata=vec_ky_new, ydata=factor*angles_new, zdata=phaseshift_vs_angleky.T)
    
    # Highlight the dominant angle
    if highlight_dominant_angle: 
        print("\n Dominant angle versus ky")
        print(dominantangle_vs_ky)
        if ax_vs_kx: ax_vs_kx.plot(vec_kx[dominantangle_vs_kx>0], dominantangle_vs_kx[dominantangle_vs_kx>0], color="red", lw=0, ms=5, marker="o")
        if ax_vs_kx: ax_vs_kx.plot(vec_kx[dominantangle_vs_kx<0], dominantangle_vs_kx[dominantangle_vs_kx<0], color="blue", lw=0, ms=5, marker="o")
        if ax_vs_ky: ax_vs_ky.plot(vec_ky[dominantangle_vs_ky>0], dominantangle_vs_ky[dominantangle_vs_ky>0], color="red", lw=0, ms=1, marker="o")
        if ax_vs_ky: ax_vs_ky.plot(vec_ky[dominantangle_vs_ky<0], dominantangle_vs_ky[dominantangle_vs_ky<0], color="blue", lw=0, ms=1, marker="o")

    # Finish the plot
    if ax_vs_kx: ax_vs_kx.set_ylim([-180,180])    
    if ax_vs_ky: ax_vs_ky.set_ylim([-180,180])             
    if ax_vs_kx: ax_vs_kx.set_xlabel(standardLabels["normalized"]["kx"])
    if ax_vs_ky: ax_vs_ky.set_xlabel(standardLabels["normalized"]["ky"])
    if ax_vs_kx: ax_vs_kx.yaxis.set_ticks([-180, -90, 0, 90, 180])
    if ax_vs_ky: ax_vs_ky.yaxis.set_ticks([-180, -90, 0, 90, 180])
    if ax_vs_kx: ax_vs_kx.set_xlim([0, np.max(simulation.vec.kx)])
    if ax_vs_ky: ax_vs_ky.set_xlim([simulation.vec.ky[1], np.max(simulation.vec.ky)])
    
    # Label y-axis
    labels = {"phi" : "\\varphi", "n0" : "\\delta n_0", "n1" : "\\delta n_1", "T0" : "\\delta T_0", "T1" : "\\delta T_1"}
    if ax_vs_kx: ax_vs_kx.set_ylabel("angle$("+labels[y1_quantity]+"/"+labels[y2_quantity]+") (^\\circ)$", labelpad=-5)
    if ax_vs_ky: ax_vs_ky.set_ylabel("angle$("+labels[y1_quantity]+"/"+labels[y2_quantity]+") (^\\circ)$", labelpad=-5)
    return  

#---------------------------------------------
def average_angles(weights, angles):
    """Average (mean) of angles

    Return the average of an input sequence of angles. The result is between
    ``0`` and ``2 * math.pi``.
    If the average is not defined (e.g. ``average_angles([0, math.pi]))``,
    a ``ValueError`` is raised.
    """

    weights = weights/np.sum(weights)
    angles = angles/180*np.pi
    x = sum([ weights[i]*math.cos(angles[i]) for i in range(len(angles)) ])
    y = sum([ weights[i]*math.sin(angles[i]) for i in range(len(angles)) ])

    if x == 0 and y == 0:
        raise ValueError(
            "The angle average of the inputs is undefined: %r" % angles)

    return math.atan2(y, x)/np.pi*180

#---------------------------
def get_phase_shift(y1, y2):
    vec_phaseshifts = np.angle(y1) - np.angle(y2)
    vec_phaseshifts[vec_phaseshifts>np.pi] = vec_phaseshifts[vec_phaseshifts>np.pi]-2*np.pi
    vec_phaseshifts[vec_phaseshifts<-np.pi] = vec_phaseshifts[vec_phaseshifts<-np.pi]+2*np.pi
    return np.array(vec_phaseshifts)/np.pi*180

#---------------------------------------------
def create_surface_plot(ax=None, xdata=None, ydata=None, zdata=None, zlabel=None, remove_color_bar=True, cmap=None):

    # Load the modules
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # For z(x,y) the columns refer to x and the rows to y
    zdata = np.transpose(zdata)

    # Set the labels    
    ax.axis([xdata[0], xdata[np.size(xdata)-1], ydata[0], ydata[np.size(ydata)-1]])

    # Load the standard jet map
    if cmap==None: cmap = copy.copy(plt.get_cmap('afmhot_r'))
    else: cmap = copy.copy(plt.get_cmap(cmap))
    
    # Make the beginning pure white
    if map=="Reds" or map=="Blues":
        cdict = cmap._segmentdata
        cdict['red'][0] = (0, 1, 1)
        cdict['green'][0] = (0, 1, 1)
        cdict['blue'][0] = (0, 1, 1)
        cmap = LinearSegmentedColormap('my_cmap', cdict)

    # Show nan values in black 
    cmap.set_bad(color='black')

    # Plot the surface plotz(x,y)  
    img = ax.pcolormesh(xdata, ydata, np.real(zdata), cmap=cmap, vmin=0, vmax=1.2)   
        
    # Add colorbar 
    if remove_color_bar: cbar = None
    if not remove_color_bar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.1)
        cbar = plt.colorbar(img, cax=cax) 
        cbar.set_label(zlabel)
        cbar.update_ticks() 
        cbar.formatter.set_powerlimits((0,0))
        cbar.formatter.set_scientific(True)

    # Set axis
    ax.set_xlim(xmin=np.min(xdata), xmax=np.max(xdata))
    ax.set_ylim(ymin=np.min(ydata), ymax=np.max(ydata)) 
    return ax, cbar

#---------------------------------------------
def interpolate_data(x,y,z,step=20,interpolate=True): 
    if interpolate: 
        function = interp2d(x, y, z.T, kind='linear')
        xnew = np.linspace(x[0], x[-1], int(len(x))*step)
        ynew = np.linspace(y[0], y[-1], int(len(y))*step)
        z_interp = function(xnew,ynew)
        x_interp, y_interp = np.meshgrid(xnew, ynew)
        x_interp = x_interp[0,:]; y_interp = y_interp[:,0] 
    return x_interp, y_interp, z_interp.T

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface
    bash = Bash(plot_phaseshift_vs_wavenumber, __doc__)  
    
    # Launch the script
    plot_phaseshift_vs_wavenumber(**bash.get_arguments())  


