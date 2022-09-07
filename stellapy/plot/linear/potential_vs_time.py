#!/usr/bin/python3  
import sys, os
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt

# Stellapy package
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])   
from stellapy.plot.utils.labels.standardLabels import standardLabels
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.simulations.Research import create_research 

#===============================================================================
#                               Plot phi2(time)                                #
#===============================================================================

def plot_potential_vs_time(folder, folderIsExperiment=False, research=None): 
    
    # Initiate the mode objects 
    modes=[] 
  
    # Create <simulations> based on the given <folder>
    if not research: research = create_research(folders=folder, folderIsExperiment=folderIsExperiment) 
    
    # Collect the modes
    for experiment in research.experiments: 
        for simulation in experiment.simulations:
            modes += simulation.modes
     
    # Create a figure
    xlabel = standardLabels["normalized"]["t"]
    ylabel = standardLabels["normalized"]["phi2"]
    title = "Time evolution of the potential squared"
    ax = create_figure(xlabel, ylabel, title)    
    
    # Plot phi2(t)
    subplot_potential_vs_time(ax, modes) 
    
    # Appearance
    ax.set_yscale("log")         
            
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder 
    plt.show()
    return

#----------------------------------------
def subplot_potential_vs_time(ax, modes):
    
    # For linear simulations, iterate over the modes 
    for mode in modes: mode.line_label = "$k_y\\rho_i = "+str(round(mode.ky,2))+"$" 
    
    # Colors based on the number of simulations
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,len(modes))) 
    
    # Plot phi2(t)
    for i, mode in enumerate(modes):  
        
        # Get the data
        phi2 = mode.potential.phi2_vs_t.phi2  
        time = mode.potential.phi2_vs_t.t
        ax.plot(time, phi2[:], color=colors[i], label=mode.line_label)  
    
    # Appearance
    ax.set_ylabel(standardLabels["normalized"]["phi2"]) 
    ax.legend(labelspacing=0.0, handlelength=1, shadow=True)
    ax.set_xlim(left=0)
    return
            
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
# This script is redirected from "stellapy/nonlinear/potential_vs_time" for linear simulations.  


