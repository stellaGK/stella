"""

#===============================================================================
#                               Plot phi2(time)                                #
#===============================================================================

For <research>, created for <folder>, the time evolution of the potential squared
is plotted. If there is only one <simulation>, the zonal contributions are shown
as well. If the simulations are linear, then plot/linear/potential_vs_time.py is
called instead.

Hanne Thienpondt 
07/09/2022

"""


#!/usr/bin/python3  
import sys, os
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt

# Stellapy package
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])  
from stellapy.plot.linear.potential_vs_time import plot_potential_vs_time as linear_plot_potential_vs_time
from stellapy.plot.utils.labels.standardLabels import standardLabels 
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.simulations.Research import create_research
from stellapy.utils.commandprompt.bash import Bash
from stellapy.GUI.plot.utils import Axis, Plot

#===============================================================================
#                               Plot phi2(time)                                #
#===============================================================================

def plot_potential_vs_time(folder, trange=None, folderIsExperiment=False, zonal=None, log=False): 
  
    # Initiate the simulation objects
    simulations = [] 
    
    # Create <simulations> based on the given <folder>
    research = create_research(folders=folder, resolutionScan=False, ignoreResolution=False, folderIsExperiment=folderIsExperiment) 
    
    # If we have linear simulations, call the linear <potential_vs_time>
    if research.experiments[0].simulations[0].linear: linear_plot_potential_vs_time(folder, trange, folderIsExperiment, research); return 
        
    # Collect the simulations
    for experiment in research.experiments:
        simulations += [s for s in experiment.simulations ] 

    # Update the time frame
    if trange!=None:   
        for simulation in experiment.simulations:
            if trange=="default": simulation.time.update_timeFrame(["50%%", "100%%"])
            if trange!="default": simulation.time.update_timeFrame(trange) 
     
    # Create a figure
    xlabel = standardLabels["normalized"]["t"]
    ylabel = standardLabels["normalized"]["phi2"]
    title = "Time evolution of the potential squared"
    ax = create_figure(xlabel, ylabel, title)    
    
    # Plot phi2(t)
    if len(simulations)==1: subplot_potential_vs_time_zonal_contributions(ax, simulations[0], log)
    if len(simulations)!=1: subplot_potential_vs_time(ax, simulations, zonal, log)   
    
    # Appearance
    if log: ax.set_yscale("log")         
            
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    plt.show()
    return

#----------------------------------------
def subplot_potential_vs_time_zonal_contributions(ax, simulation, log=False):

    # Get the data
    phi2 = simulation.potential.phi2_vs_t.phi2 
    phi2_zonal = simulation.potential.phi2_vs_t_zonal.phi2_zonal
    phi2_nozonal = simulation.potential.phi2_vs_t_nozonal.phi2_nozonal
    time = simulation.potential.phi2_vs_t.t
    
    # Plot the data
    ax.plot(time, phi2[:], color="black", label="$\\langle \\varphi^2 \\rangle_z$")  
    ax.plot(time, phi2_zonal[:], color="navy", label="$\\langle \\varphi^2_\\text{Z} \\rangle_z$")  
    ax.plot(time, phi2_nozonal[:], color="crimson", label="$\\langle \\varphi^2_\\text{NZ} \\rangle_z$")  
    
    # Appearance
    ax.legend(labelspacing=0.5, handlelength=1, shadow=True)
    ax.set_xlim(left=0, right=np.max(time))
    if log==True:  ax.set_ylim(bottom=np.min(phi2), top=np.max(phi2))
    if log==False: ax.set_ylim(bottom=0, top=np.max(phi2)*1.1)
    return

#----------------------------------------
def subplot_potential_vs_time(ax, simulations, zonal=None, log=False):
    
    # Automate the axis limits  
    axis = Axis(ax, Plot(),  xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.1, log=log, percentage=0.5)
    
    # Colors based on the number of simulations
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,len(simulations))) 
    
    # Plot phi2(t)
    for i, simulation in enumerate(simulations):  
        
        # Get the data
        if zonal==None: phi2 = simulation.potential.phi2_vs_t.phi2 
        if zonal==True: phi2 = simulation.potential.phi2_vs_t_zonal.phi2_zonal 
        if zonal==False: phi2 = simulation.potential.phi2_vs_t_nozonal.phi2_nozonal 
        time = simulation.potential.phi2_vs_t.t
        ax.plot(time, phi2[:], color=colors[i], label=simulation.line_label)  
        axis.update_axisLimits(time, phi2)  
    
    # Appearance
    if zonal==None: ax.set_ylabel(standardLabels["normalized"]["phi2"])
    if zonal==True: ax.set_ylabel(standardLabels["normalized"]["phi2_zonal"])
    if zonal==False: ax.set_ylabel(standardLabels["normalized"]["phi2_nozonal"])
    ax.legend(labelspacing=0.0, handlelength=1, shadow=True)
    
    # Axis limits
    axis.rescale_axis()
    return
            
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
# This script is redirects to "stellapy/linear/potential_vs_time" for linear simulations.
 
if __name__ == "__main__":
    bash = Bash(plot_potential_vs_time, __doc__)        
    
    # Data toggles for the x- and y-quantities
    bash.add_toggle('y_quantity', 'phi2', 'Plot the time evolution of the potential squared.') 
    
    # Logarithmic axis
    bash.add_option('log', 'True', 'l', 'Use a logarithmic y-axis.') 
    
    # Change the time range for calculation of the saturated quantities
    bash.add_option('trange', 'range', 't', 'Select the time frame for the simulations.') 
    
    # Research options
    bash.add_option('folderIsExperiment', 'True', 'f', 'Create an experiment for each folder') 
    
    # Plot only the zonal modes or the non-zonal modes
    bash.add_option('zonal', 'True', 'z', 'Plot only the zonal modes.') 
    bash.add_option('zonal', 'False', 'n', 'Plot only the non-zonal modes.') 
    
    # Get the arguments and execute the script
    plot_potential_vs_time(**bash.get_arguments())   


