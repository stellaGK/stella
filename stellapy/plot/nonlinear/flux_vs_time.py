"""

#===============================================================================
#                    Plot qflux(t) or pflux(t) or vflux(t)                     #
#===============================================================================

If <folder> contains only one simulation, then the time evolution for the fluxes 
is shown for both the ions and electrons. If <folder> contains multiple simulations, 
the time evolution of the fluxes is shown for <specie>.

y_quantity
----------
Choose between {qflux, pflux, vflux}.
    >> plot_flux_vs_time --qflux (or --pflux or --vflux)

specie
------
Choose an integer X refering to species_parameters_X.
    >> plot_flux_vs_time -s 0 (or -s 1, -s 2, ... or --specie 0, --specie 1, ...)

parameter
---------
For <research> a new <experiment> is created for each unique value of <parameter>. 
    >> plot_flux_vs_time -p rho (or --parameter rho)

t_range
-------
Change the time range where the simulation is considered saturated. 
    >> plot_saturatedflux_vs_parameter -t "[500, 1000]" (or -t "500, 1000" or --trange "[500, 1000]")
    
folderIsExperiment
------------------
To ignore the automatic sorting of experiment and simulation, you can force stellapy
to create an <experiment> for each <subfolder>, the folder name will be the label
    >> plot_flux_vs_time -f (or --folderIsExperiment)

Hanne Thienpondt 
01/09/2022

"""

#!/usr/bin/python3  
import sys, os
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt  

# Stellapy package
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])   
from stellapy.plot.utils.labels.standardParameters import standardParameters
from stellapy.plot.utils.species.recognize_species import recognize_species
from stellapy.plot.utils.labels.standardLabels import standardLabels
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.simulations.Research import create_research
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                    Plot qflux(t) or pflux(t) or vflux(t)                     #
#===============================================================================

def plot_flux_vs_time(folder, y_quantity="qflux", specie=0, trange=None, 
        folderIsExperiment=False, parameter="-", logy=False):  
    
    # Create <experiments> and <simulations> based on the given <folder>
    knob = standardParameters[parameter]["knob"]; key = standardParameters[parameter]["key"]
    try: research = create_research(folders=folder, knob1=knob, key1=key, resolutionScan=False, ignoreResolution=True, folderIsExperiment=folderIsExperiment)
    except: research = create_research(folders=folder, knob1=knob, key1=key, resolutionScan=False, ignoreResolution=True, folderIsExperiment=True)

    # Update the time frame
    if trange!=None:    
        for experiment in research.experiments:
            for simulation in experiment.simulations: 
                if trange=="default": simulation.time.update_timeFrame(["50%%", "100%%"])
                if trange!="default": simulation.time.update_timeFrame(trange) 
     
    # Create a figure 
    if y_quantity=="qflux": title = "Time evolution of the heat flux"
    if y_quantity=="pflux": title = "Time evolution of the particle flux"
    if y_quantity=="vflux": title = "Time evolution of the momentum flux"
    ax = create_figure(title=title)   
    
    # Plot qflux(time)
    if research.numberOfSimulations==1: subplot_flux_vs_time_bothspecies(ax, research, y_quantity)
    if research.numberOfSimulations!=1: subplot_flux_vs_time(ax, research, y_quantity, specie)
    
    # Appearance    
    ax.yaxis.labelpad = 15; ax.xaxis.labelpad = 10
    if logy: ax.set_yscale("log"); ax.set_ylim([10e-12, 100])
        
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    plt.show()
    return

#----------------------------------------
def subplot_flux_vs_time(ax, research, y_quantity="qflux", specie=0):
    
    # Colors based on the number of simulations
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,research.numberOfSimulations)) 
    max_time = 0; min_flux = 0; max_flux = 0; count = 0
 
    # Plot qflux(t) for one of the species
    for experiment in research.experiments:
        for simulation in experiment.simulations: 
        
            # Get the data
            if y_quantity=="qflux": flux = simulation.fluxes.qflux_vs_ts.qflux[:,specie]  
            if y_quantity=="pflux": flux = simulation.fluxes.pflux_vs_ts.pflux[:,specie]  
            if y_quantity=="vflux": flux = simulation.fluxes.vflux_vs_ts.vflux[:,specie]    
            time = simulation.fluxes.qflux_vs_ts.t
            
            # Get the saturated value
            tstart = simulation.time.tstart
            tend = simulation.time.tend
            satflux = np.mean(flux[(time>=tstart) & (time<=tend)]) 
            
            # Get the label
            s = recognize_species(research, specie)
            label = "Q" if y_quantity=="qflux" else ("\\Gamma" if y_quantity=="pflux" else "\\Pi")
            label = "$"+label+"_{\\text{sat},"+s+"}/"+label+"_{\\text{gB},i} = "+str("{:.2f}".format(satflux))+"$"
            label = simulation.line_label+"; "+label if research.numberOfSimulations!=research.numberOfExperiments else experiment.line_label+"; "+label  
            
            # Plot the time evolution and the saturated value
            ax.plot(time, flux, color=colors[count])  
            ax.plot([tstart, tend], [satflux,satflux], color=colors[count], label=label) 
            
            # Keep track of the axis limits
            max_time = np.max([max_time, np.max(time)])
            min_flux = np.min([min_flux, np.min(flux[time>tstart])*1.2])
            max_flux = np.max([max_flux, np.max(flux[time>tstart])*1.2])
            count += 1

    # Correct the y-label
    s = recognize_species(research, specie) 
    ax.set_xlabel(standardLabels["normalized"]["t"])
    ax.set_ylabel(standardLabels["normalized"][y_quantity].replace("_{s}", "_{"+s+"}"))
    ax.legend(labelspacing=0.0, prop={'size':20}, handlelength=0.8, handletextpad=0.5)
    
    # Axis limits
    ax.set_xlim(left=0, right=max_time)  
    ax.set_ylim([min_flux, max_flux])
    return 

#----------------------------------------
def subplot_flux_vs_time_bothspecies(ax, research, y_quantity="qflux"):
    """ If we have a single simulation, plot both the ion and electron flux. """
    
    # Get the first simulation
    simulation = research.experiments[0].simulations[0]
        
    # Get the data
    if y_quantity=="qflux": flux = simulation.fluxes.qflux_vs_ts.qflux 
    if y_quantity=="pflux": flux = simulation.fluxes.pflux_vs_ts.pflux 
    if y_quantity=="vflux": flux = simulation.fluxes.vflux_vs_ts.vflux 
    time = simulation.fluxes.qflux_vs_ts.t
    
    # Get the saturated value
    tstart = simulation.time.tstart
    tend = simulation.time.tend
    satflux = np.mean(flux[(time>=tstart) & (time<=tend)], axis=0)
    
    # Colors per species
    colors = ["navy", "crimson", "black", "green", "orange"]
    species = range(simulation.dim.species)
    
    # Iterate over the species
    for ispecie, specie in enumerate(species):  
        
        # Get the label
        s = recognize_species(research, specie) 
        label = "Q" if y_quantity=="qflux" else ("\\Gamma" if y_quantity=="pflux" else "\\Pi")
        label = "$"+label+"_{\\text{sat},"+s+"}/"+label+"_{\\text{gB},i} = "+str("{:.2f}".format(satflux[specie]))+"$"
        
        # Plot the time evolution and saturated value
        ax.plot(time, flux[:,specie], color=colors[ispecie])  
        ax.plot([tstart, tend], [satflux[specie],satflux[specie]], color=colors[ispecie], label=label) 

    # Axis labels
    ax.set_xlabel(standardLabels["normalized"]["t"])
    ax.set_ylabel(standardLabels["normalized"][y_quantity])

    # Axis limits
    ax.set_xlim(left=0, right=np.max(time))  
    ax.set_ylim(bottom=np.min([0, 1.2*np.min(flux[time>tstart,:])]), top=np.max([0, 1.2*np.max(flux[time>tstart,:])]))
    return

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface
    bash = Bash(plot_flux_vs_time, __doc__)  
    
    # Data options
    bash.add_option('specie', 'int', 's', 'Select the species as "-s 0" for ions and "-s 1" for electrons.') 
    bash.add_toggle('y_quantity', 'qflux', 'Plot the time evolution of the heat flux.') 
    bash.add_toggle('y_quantity', 'pflux', 'Plot the time evolution of the particle flux.') 
    bash.add_toggle('y_quantity', 'vflux', 'Plot the time evolution of the momentum flux.') 
    
    # Research options      
    bash.add_option('folderIsExperiment', 'True', 'f', 'Each folder is an experiment.')  
    bash.add_option('parameter', 'str', 'p', 'Create an <experiment> for each unique value of <parameter>.') 
    bash.add_option('trange', 'range', 't', 'Select the time frame for the simulations.') 
    
    # Plotting options
    bash.add_option('logy', 'True', 'y', 'Logarithmic scale for the y-axis.') 
    
    # Launch the script
    plot_flux_vs_time(**bash.get_arguments())   


