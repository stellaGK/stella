"""

#===============================================================================
#                          Plot phi2(kx) and phi2(ky)                          #
#===============================================================================
 
Based on <folder> create a <research> to plot phi2(kx) and phi2(ky).
    
    
Average potential squared in real space
---------------------------------------
In stella the Fourier components hat{phi}_{kxky} are related to the convential
Fourier components through hat{phi}_{kxky} = hat{PHI}_{kxky}/NxNy. In other words,
the inverse Fourier transformation is defined as
    
    phi_{xy} = 1/NxNy sum_{kxky} hat{PHI}_{kxky} exp(i*kx*x+i*ky*y)
             = sum_{kxky} hat{phi}_{kxky} exp(i*kx*x+i*ky*y)
             
As a result, Parseval's theorem is given by

    sum_{xy} |phi_{xy}|^2 = 1/NxNy sum_{kxky} |hat{PHI}_{kxky}|^2
                          = NxNy sum_{kxky} |hat{phi}_{kxky}|^2  
                          
Therefore, the average potential squared in real space, is simply given by the sum
over the (kx,ky) Fourier components as calculated by stella. 

    |phi|^2 = 1/NxNy sum_{xy} |phi(x,y)|^2 
            = 1/(NxNy)^2 sum_{kxky} |hat{PHI}_{kxky}|^2
            = sum_{kxky} |hat{phi}_{kxky}|^2
            
Note that phi refers to the perturbed potential delta phi_1 in the stella bible.

x_quantity
----------
Choose between {kx, ky}.
    >> plot_potential_spectra --ky
    >> plot_potential_vs_ky 
    
folderIsExperiment
------------------
To ignore the automatic sorting of experiment and simulation, you can force stellapy
to create an <experiment> for each <subfolder>, the folder name will be the label
    >> plot_potential_spectra -f (or --folderIsExperiment) 

Hanne Thienpondt 
15/12/2022

"""

#!/usr/bin/python3  
import sys, os
import pathlib
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec

# Personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import get_styleForLinesAndMarkers
from stellapy.plot.utils.labels.get_timeFrameString import get_timeFrameString 
from stellapy.plot.utils.style.create_figure import update_figure_style
from stellapy.plot.utils.labels.standardLabels import standardLabels   
from stellapy.simulations.Research import create_research   
from stellapy.plot.utils.style import Axis, Legend, Plot 
from stellapy.utils.commandprompt.bash import Bash 

#===============================================================================
#                          Plot phi2(kx) and phi2(ky)                          #
#===============================================================================

def plot_potential_vs_wavenumber(
        folder, 
        # Quantities to be plotted
        x_quantity="both",
        y_quantity="phi2",  
        # Plotting options
        log=False, 
        folderIsExperiment=False, 
        normalize_to_one=False, ): 
    
    # Create a <research> based on the given <folder>
    try: research = create_research(folders=folder, folderIsExperiment=folderIsExperiment)
    except: research = create_research(folders=folder, folderIsExperiment=True)
    
    # Get the x-quantity
    nk = 2 if x_quantity=="both" else 1
    x1_quantity = "kx" if x_quantity=="both" else x_quantity
    
    # Create a figure
    if y_quantity=="phi2":  title = "Potential squared spectra"
    fig = plt.figure(figsize=(18, 9)); axes = []
    grid_specifications = gridspec.GridSpec(nk, 1)
    grid_specifications.update(top=0.93, left=0.08, right=0.97, bottom=0.1, wspace=0.2, hspace=0.3)
    for i in range(nk): axes.append(plt.subplot(grid_specifications[i]))
    update_figure_style(fig, axes)  
    fig.suptitle(title)
    
    # Add the data to the plot
    subplot_potential_vs_wavenumber(axes[0], research, x1_quantity, y_quantity, log=log, normalize_to_one=normalize_to_one)
    if nk==2: subplot_potential_vs_wavenumber(axes[1], research, "ky", y_quantity, log=log, normalize_to_one=normalize_to_one)
    
    # Appearance  
    for ax in axes: (ax.ticklabel_format(style='plain', useOffset=False) if (axes[0].get_ylim()[1]<1000) else ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))) if not log else ""
    for ax in axes: ax.xaxis.labelpad = 10
    for ax in axes: ax.yaxis.labelpad = 15
    
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    plt.show()
    return

#===============================================================================
#                           Plot phi2(kx) or phi2(ky)                          #
#===============================================================================

def subplot_potential_vs_wavenumber(
        ax, research, 
        # Quantities to be plotted
        x_quantity="ky", 
        y_quantity="phi2",  
        # Plotting options
        log=True, 
        normalize_to_one=False):
    
    # Automate the axis limits and legend
    plot = Plot(); legend = Legend(ax, plot) 
    axis = Axis(ax, plot, xbot_pos=0, ytop_neg=0, overshoot_y=1.1) 
    plot.process_plottingVariables(research)
    tstarts = [np.nan,np.nan]; tends = [np.nan,np.nan] 
    
    # Iterate over the experiments and simulations
    for experiment in research.experiments:
        for simulation in experiment.simulations:
            
            # Get the lines and marker styles for the experiments and simulations
            style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation)
            
            # Get flux and kx or ky
            x, y, tstarts, tends = get_quantity_data(simulation, x_quantity, y_quantity, tstarts, tends)
            if normalize_to_one: y = y/np.max(np.abs(y))
            if np.any(y<0): log = False
            
            # Plot qflux(kx) or qflux(ky) or pflux(kx) or vflux(kx) or ...
            ax.plot(x, y, marker="o", **style)  
        
            # Keep track of the axis limits
            axis.update_axisLimits(x, y)
    
    # Axis labels
    extra = "$\\sum_{k_x}$" if x_quantity=="ky" else "$\\sum_{k_y}$"
    ylabel = standardLabels["normalized"][y_quantity+"_fourier"]
    ylabel = extra+"$\\langle$" + ylabel + "$\\rangle_{z, t="+get_timeFrameString(tstarts, tends)+"}$"
    ax.set_xlabel(standardLabels["normalized"][x_quantity]) 
    ax.set_ylabel(ylabel) 

    # Automatically set the axis limits and legends 
    if normalize_to_one: ax.set_ylim([0,1])
    if log: ax.set_yscale("log")
    legend.add_legend()
    axis.rescale_axis()  
    return 

#===============================================================================
#                           Read phi2(kx) or phi2(ky)                          #
#===============================================================================

def get_quantity_data(simulation, x_quantity, y_quantity, tstarts, tends):
    
    # Data for the x-axis
    if x_quantity=="kx": x = simulation.vec.kx 
    if x_quantity=="ky": x = simulation.vec.ky 
    
    # Get the data for the y-axis  
    try: y, t0, t1 = get_quantity_data_from3DFiles(simulation, x_quantity, y_quantity)
    except: y, t0, t1 = get_quantity_data_fromSaturatedFiles(simulation, x_quantity, y_quantity)

    # Keep track of the time frames
    tstarts[0] = np.nanmin([tstarts[0], t0]) 
    tstarts[1] = np.nanmax([tstarts[1], t0]) 
    tends[0] = np.nanmin([tends[0], t1]) 
    tends[1] = np.nanmax([tends[1], t1]) 
    return x, y, tstarts, tends 

#-------------------------------
def get_quantity_data_from3DFiles(simulation, x_quantity, y_quantity,): 
    if x_quantity=="kx": 
        if y_quantity=="phi2":  y = simulation.potential.phi2_vs_tkx.phi2[:,:]; t = simulation.potential.phi2_vs_tkx.t
    if x_quantity=="ky": 
        if y_quantity=="phi2":  y = simulation.potential.phi2_vs_tky.phi2[:,:]; t = simulation.potential.phi2_vs_tky.t
    y = np.nanmean(y[t>simulation.time.tstart,:],axis=0) 
    return y, simulation.time.tstart, t[-1]
        
#-------------------------------
def get_quantity_data_fromSaturatedFiles(simulation, x_quantity, y_quantity):
    trange = simulation.saturated.trange
    if y_quantity=="phi2":  y = simulation.saturated.phi2_vs_kxky.phi2[:,:]
    if y_quantity!="phi2":  y = np.sum(y[:,:,:]*simulation.geometry.dl_over_B[:,np.newaxis,np.newaxis], axis=0) 
    return y, trange[0], trange[1]

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash interface
    # Toggles are defined through (name, value, shortoption, longoption, explanation)
    # Options are defined through (name, datatype, shortoption, longoption, explanation)
    bash = Bash(plot_potential_vs_wavenumber, __doc__)     
    
    # Data toggles for the x-quantities
    bash.add_toggleheader("x-quantity")
    bash.add_toggle('x_quantity', 'kx', '', '', 'Plot the kx-spectra.') 
    bash.add_toggle('x_quantity', 'ky', '', '', 'Plot the ky-spectra.') 
    bash.add_toggle('x_quantity', 'both', '', '', 'Plot the kx- and ky-spectra.') 
    bash.add_togglespace()
    
    # Research options      
    bash.add_toggleheader("other")
    bash.add_toggle('folderIsExperiment', True, 'f', 'folderIsExperiment', 'Each folder is an experiment.')    
    
    # Plot options  
    bash.add_toggle('normalize_to_one', True, 'n', 'normalize_to_one', 'Rescale the data to have the maximum at one.')  
    bash.add_toggle('log', True, 'l', 'log', 'Use logaritmic scales.')  
    
    # Get the arguments and extract <x_quantity>
    args = bash.get_arguments() 
    
    # Launch the script
    plot_potential_vs_wavenumber(**args)    


