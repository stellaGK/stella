"""

#===============================================================================
#         Plot qflux(parameter) or pflux(parameter) or vflux(parameter)        #
#===============================================================================

Plot the saturated flux versus the chosen <parameter>. 

y_quantity
----------
Choose between {qflux, pflux, vflux}.
    >> plot_saturatedflux_vs_parameter --qflux (or --pflux or --vflux)

specie
------
Choose an integer X refering to species_parameters_X.
    >> plot_saturatedflux_vs_parameter -s 0 (or -s 1, -s 2, ... or --specie 0, --specie 1, ...)

parameter
---------
Quantity plotted on the x-axis, should be a stella input variable.
    >> plot_saturatedflux_vs_parameter -p fprim (or --parameter fprim)

experiment_parameter1
---------
For <research> a new <experiment> is created for each unique value of <parameter>. 
    >> plot_saturatedflux_vs_parameter -e rho (or --parameter rho)

t_range
-------
Change the time range where the simulation is considered saturated. 
    >> plot_saturatedflux_vs_parameter -t [500, 1000] (or -t 500, 1000 or --trange [500, 1000])
    
folderIsExperiment
------------------
To ignore the automatic sorting of experiment and simulation, you can force stellapy
to create an <experiment> for each <subfolder>, the folder name will be the label
    >> plot_saturatedflux_vs_parameter -f (or --folderIsExperiment)
    
resolutionScan
------------------
For a resolution scan, turn this on.
    >> plot_saturatedflux_vs_parameter -r (or --resolutionScan)

Hanne Thienpondt 
27/10/2022

"""

#!/usr/bin/python3   
import sys, os
import pathlib
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt 

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import get_styleForLinesAndMarkers 
from stellapy.plot.utils.data.recognize_varied_parameter import recognize_varied_parameter
from stellapy.plot.utils.data.recognize_resolution_scan import recognize_resolution_scan
from stellapy.plot.utils.labels.standardParameters import standardParameters
from stellapy.plot.utils.species.recognize_species import recognize_species
from stellapy.plot.utils.labels.standardLabels import standardLabels
from stellapy.plot.utils.data.get_parameters import get_parameters 
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.simulations.Research import create_research 
from stellapy.utils.commandprompt.bash import Bash
from stellapy.plot.utils.style import Axis, Legend
from stellapy.plot.utils.style.Plot import Plot

#===============================================================================
#         Plot qflux(parameter) or pflux(parameter) or vflux(parameter)        #
#===============================================================================

def plot_flux_vs_parameter(
        folder, 
        # Quantities to plot
        specie=0, 
        x_quantity=None,
        y_quantity="qflux", 
        # Time frame to calculate saturated quantities 
        trange=None,  
        # Research
        experiment_parameter1=None, 
        experiment_parameter2=None, 
        folderIsExperiment=False, 
        resolutionScan=None): 
     
    # Find the stella knob and key of the given parameter
    knob1, key1, knob2, key2 = get_knobsAndKeys(experiment_parameter1, experiment_parameter2) 
    
    # Automatically detect <x_quantity> and <resolutionScan> if they're set to None
    resolutionScan = recognize_resolution_scan(folder, resolutionScan)
    x_quantity = recognize_varied_parameter(folder, x_quantity)
    
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder, knob1=knob1, key1=key1, knob2=knob2, key2=key2, resolutionScan=resolutionScan, ignore_resolution=(not resolutionScan), folderIsExperiment=folderIsExperiment) 
    research.update_timeFrame(trange) 
     
    # Create a figure
    parameter_label = standardLabels["normalized"][standardParameters[x_quantity]["key"]]
    if y_quantity=="qflux": title = "Saturated heat flux versus "+parameter_label
    if y_quantity=="pflux": title = "Saturated particle flux versus "+parameter_label
    if y_quantity=="vflux": title = "Saturated momentum flux versus "+parameter_label
    ax = create_figure(title=title)  

    # Add the data to the plot
    subplot_flux_vs_parameter(ax, research, x_quantity=x_quantity, y_quantity=y_quantity, specie=specie)
    
    # Change the appearance
    ax.ticklabel_format(style='plain', useOffset=False)
    ax.yaxis.labelpad = 15
    ax.xaxis.labelpad = 10
    
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0: plt.show()
    return

#-----------------------------
def subplot_flux_vs_parameter(
        # Axis
        ax,
        # Research 
        research, 
        # Quantities to be plotted
        x_quantity="tiprim",
        y_quantity="qflux", 
        specie=0, 
        # Plotting
        fontsize=20):
    
    # Get the knob and key of the chosen parameter (x_quantity)
    knob = standardParameters[x_quantity]["knob"]
    key = standardParameters[x_quantity]["key"]

    # Automate the axis limits and legend
    plot = Plot(); plot.update_legend(fontsize=fontsize)
    axis = Axis(ax, plot, xbot_pos=0, ytop_neg=0, ybot_pos=0)
    legend = Legend(ax, plot)
    
    # Axis labels
    s = recognize_species(research, specie)
    ax.set_xlabel(standardLabels["normalized"][key])
    ax.set_ylabel(standardLabels["normalized"][y_quantity].replace("_{s}", "_{"+s+"}"))
    
    # Plot qflux(parameter) or pflux(parameter) or vflux(parameter)
    for experiment in research.experiments:
            
        # Get the parameters that are scanned
        parameters, simulations = get_parameters(experiment, key, knob)
    
        # Initiate the data
        satflux = [np.NaN]*len(experiment.simulations)
        stdflux = [np.NaN]*len(experiment.simulations) 

        # Get the saturated data
        for i, simulation in enumerate(simulations):  
            satflux[i] = simulation.time.saturatedFluxes[y_quantity][specie]
            stdflux[i] = simulation.time.satFluxStdErrors[y_quantity][specie]
            
        # Plot the saturated data
        style = get_styleForLinesAndMarkers(plot, legend, research, experiment)  
        ax.errorbar(parameters, satflux, yerr=stdflux, capsize=2, marker="o", **style)  
        
        # Keep track of the axis limits 
        axis.update_axisLimits(parameters, satflux)
            
    # Automatically set the axis limits and legends 
    legend.add_legend()
    axis.rescale_axis()
    return

#--------------------------------------------
def get_knobsAndKeys(experiment_parameter1, experiment_parameter2): 
    if experiment_parameter1==None: experiment_parameter1 = "-"
    if experiment_parameter2==None: experiment_parameter2 = "-"
    knob1 = standardParameters[experiment_parameter1]["knob"]
    key1 = standardParameters[experiment_parameter1]["key"]
    knob2 = standardParameters[experiment_parameter2]["knob"]
    key2 = standardParameters[experiment_parameter2]["key"] 
    return knob1, key1, knob2, key2

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    bash = Bash(plot_flux_vs_parameter, __doc__)      
    
    # Quantities to be plotted
    bash.add_option('x_quantity', 'str', 'x', '', """Choose the x-quantity from
    {rho, tiprim, teprim, fprim, tite, delt, nmu, nvgrid, rk, \n\
    dvpa, dmu, nz, nzed, nzgrid, nx, ny, y0, kx max, ky max\n\
    dkx, dky, Lx, Ly, nfield, pol.turns, nperiod, D_hyper}.""") 
    bash.add_option('y_quantity', 'str', 'y', '', 'Choose the y-quantity from {qflux, pflux, vflux}.') 
    bash.add_toggle('y_quantity', 'qflux', '', '', 'Plot the time evolution of the heat flux.') 
    bash.add_toggle('y_quantity', 'pflux', '', '', 'Plot the time evolution of the particle flux.') 
    bash.add_toggle('y_quantity', 'vflux', '', '', 'Plot the time evolution of the momentum flux.') 
    bash.add_option('specie', 'int', 's', '', 'Select the species as "-s 0" for ions and "-s 1" for electrons.') 
    
    # Time range  
    bash.add_option('trange', 'range', 't', '', 'Select the time frame for the simulations through e.g. "-t [500, 1000]".') 
    
#     # Research      
#     bash.add_option('resolutionScan', 'True', 'r', 'Create a <research> for a resolution scan.')   
#     bash.add_option('folderIsExperiment', 'True', 'f', 'Each subfolder in <folder> is an <experiment>.')   
#     bash.add_option('experiment_parameter1', 'str', 'e', """Create an <experiment> for each unique value of <experiment_parameter1> in 
#     {rho, tiprim, teprim, fprim, tite, delt, nmu, nvgrid, rk, \n\
#     dvpa, dmu, nz, nzed, nzgrid, nx, ny, y0, kx max, ky max\n\
#     dkx, dky, Lx, Ly, nfield, pol.turns, nperiod, D_hyper}""") 
#     bash.add_option('experiment_parameter2', 'str', 'E', """Create an <experiment> for each unique value of <experiment_parameter2> in
#     {rho, tiprim, teprim, fprim, tite, delt, nmu, nvgrid, rk, \n\
#     dvpa, dmu, nz, nzed, nzgrid, nx, ny, y0, kx max, ky max\n\
#     dkx, dky, Lx, Ly, nfield, pol.turns, nperiod, D_hyper}""") 
    
    # Launch the script 
    plot_flux_vs_parameter(**bash.get_arguments())   

    
    
