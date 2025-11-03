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
15/12/2022

"""

#!/usr/bin/python3  
import sys, os
import pathlib
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt  

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.plot.utils.labels.standardParameters import standardParameters
from stellapy.plot.utils.species.recognize_species import recognize_species
from stellapy.plot.utils.labels.standardLabels import standardLabels
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.Research import create_research
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                    Plot qflux(t) or pflux(t) or vflux(t)                     #
#===============================================================================

def plot_flux_vs_time(
        folder, 
        # Quantities to be plotted
        y_quantity="qflux", 
        specie=0, 
        # Time frame to calculate saturated quantities 
        trange=None, 
        # Research
        folderIsExperiment=False, 
        folderIsSimulation=False,
        experiment_parameter="-", 
        # Plotting
        xrange=None,
        yrange=None,
        logy=False):  

    # Create <experiments> and <simulations> based on the given <folder>
    knob = standardParameters[experiment_parameter]["knob"]; key = standardParameters[experiment_parameter]["key"]
    research = create_research(folders=folder, knob1=knob, key1=key, resolutionScan=False, ignore_resolution=True, folderIsExperiment=folderIsExperiment, folderIsSimulation=folderIsSimulation)
    research.update_timeFrame(trange); experiments = research.experiments
    
    # Create a figure for each experiment
    for experiment in experiments: 
     
        # Create a figure 
        if y_quantity=="qflux": title = "Time evolution of the heat flux"
        if y_quantity=="pflux": title = "Time evolution of the particle flux"
        if y_quantity=="vflux": title = "Time evolution of the momentum flux"
        if research.numberOfExperiments>1: title += ": "+experiment.line_label
        ax = create_figure(title=title)   
        
        # Plot qflux(time)
        research.experiments = [experiment]; research.numberOfSimulations = len(experiment.simulations)
        if research.numberOfSimulations==1: subplot_flux_vs_time_bothspecies(ax, research, y_quantity)
        if research.numberOfSimulations!=1: subplot_flux_vs_time(ax, research, y_quantity=y_quantity, specie=specie)
        
        # Appearance    
        ax.yaxis.labelpad = 15; ax.xaxis.labelpad = 10
        if logy: ax.set_yscale("log"); ax.set_ylim([10e-12, 100])  
        if xrange: ax.set_xlim(xrange)    
        if yrange: ax.set_ylim(yrange)  
            
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0:  plt.show()
    return

#----------------------------------------
def subplot_flux_vs_time(
        ax, 
        research, 
        # Quantities to be plotted
        y_quantity="qflux", 
        specie=0):
    
    # Colors based on the number of simulations
    colors = plt.colormaps.get_cmap('jet')(np.linspace(0,1,research.numberOfSimulations)) 
    max_time = 0; min_flux = 0; max_flux = 0; count = 0
    
    # Check whether <specie> is a valid choice 
    if specie >= research.experiments[0].simulations[0].dim.species:
        dim_species = research.experiments[0].simulations[0].dim.species
        exit_reason = "Error: specie = "+str(specie)+" was selected but only "+str(dim_species)+" species are present in the simulation.\n"
        exit_reason += "Please choose specie from {"+", ".join([str(s) for s in range(dim_species)])+"} instead.\n"
        exit_program(exit_reason, subplot_flux_vs_time, sys._getframe().f_lineno)
        
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
            label = "$"+label+"_{\\text{sat,}"+s+"}/"+label+"_{\\text{gB},i} = "+str("{:.2f}".format(satflux))+"$"
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
    
    # Add the legend
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
        label = "$"+label+"_{\\text{sat,}"+s+"}/"+label+"_{\\text{gB},i} = "+str("{:.2f}".format(satflux[specie]))+"$"
        
        # Plot the time evolution and saturated value
        ax.plot(time, flux[:,specie], color=colors[ispecie])  
        ax.plot([tstart, tend], [satflux[specie],satflux[specie]], color=colors[ispecie], label=label) 

    # Axis labels
    ax.set_xlabel(standardLabels["normalized"]["t"])
    ax.set_ylabel(standardLabels["normalized"][y_quantity])

    # Axis limits
    ax.set_xlim(left=0, right=np.max(time))   
    ax.set_ylim(bottom=np.min([0, 1.2*np.nanmin(flux[time>tstart,:])]), top=np.max([0, 1.2*np.nanmax(flux[time>tstart,:])]))
    
    # Add the legend
    ax.legend(labelspacing=0.0, prop={'size':20}, handlelength=0.8, handletextpad=0.5)
    return

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":  
    
    # Create a bash-like interface
    # Toggles are defined through (name, value, shortoption, longoption, explanation)
    # Options are defined through (name, datatype, shortoption, longoption, explanation)
    bash = Bash(plot_flux_vs_time, __doc__)  
    
    # Quantities to be plotted
    bash.add_toggleheader("y-quantity")
    bash.add_toggle('y_quantity', 'qflux', '', '', 'Plot the time evolution of the heat flux.') 
    bash.add_toggle('y_quantity', 'pflux', '', '', 'Plot the time evolution of the particle flux.') 
    bash.add_toggle('y_quantity', 'vflux', '', '', 'Plot the time evolution of the momentum flux.') 
    bash.add_togglespace()
    
    # Species
    bash.add_toggleheader("specie")
    bash.add_toggle('specie', 0, '', 'ions', 'Plot the fluxes for the ions.') 
    bash.add_toggle('specie', 1, '', 'electrons', 'Plot the fluxes for the electrons (assuming s=1 for electrons).') 
    bash.add_toggle('specie', 2, '', 'impurity', 'Plot the fluxes for the impurities (assuming s=2 for impurities).') 
    bash.add_togglespace()
    
    # Manipulate the data
    bash.add_optionheader("simulation") 
    bash.add_toggleheader("simulation")
    bash.add_option('y_quantity', 'str', 'y', '', 'Choose the y-quantity from {qflux, pflux, vflux}.') 
    bash.add_option('specie', 'int', 's', '', 'Select the species as "-s 0" for ions and "-s 1" for electrons.') 
    bash.add_toggle('folderIsExperiment', True, 'f', 'folderIsExp', 'Each folder is an experiment.') 
    bash.add_toggle('folderIsSimulation', True, '', 'all', 'Each folder is a simulation.')  
    bash.add_option('trange', 'range', 't', '', "Select the time frame for the saturated state of the simulations, e.g. '-t [500, 1000]' or -t 500.") 
    bash.add_optionspace()
    bash.add_togglespace()  
    
    # Change the appearance of the figure
    bash.add_optionheader("figure") 
    bash.add_toggleheader("figure") 
    bash.add_toggle('logy', True, 'l', 'log', 'Use a logarithmic y-axis.')     
    bash.add_option('xrange', 'range', '', 'xrange', "Select the xrange for the figure, e.g. '[0,1000]'.") 
    bash.add_option('yrange', 'range', '', 'yrange', "Select the yrange for the figure, e.g. '[0,10]'.") 
    bash.add_optionspace()
    bash.add_togglespace()
    
    # Other options 
    bash.add_toggleheader("other")
    bash.add_optionheader("other") 
    
    # Research       
    bash.add_option('experiment_parameter', 'str', 'e', 'para', """Create an <experiment> for each unique value of <experiment_parameter1> in 
    {rho, tiprim, teprim, fprim, tite, delt, nmu, nvgrid, rk, \n\
    dvpa, dmu, nz, nzed, nzgrid, nx, ny, y0, kx max, ky max\n\
    dkx, dky, Lx, Ly, nfield, pol.turns, nperiod, D_hyper}""") 
    
    # Launch the script
    plot_flux_vs_time(**bash.get_arguments())   


