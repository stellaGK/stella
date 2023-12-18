"""
 
 ...

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
#                    Plot dens2(t) or temp2(t) or upar2(t)                     #
#===============================================================================

def plot_moment_vs_time(
        folder, 
        # Quantities to be plotted
        y_quantity="dens2", 
        specie=0, 
        # Time frame to calculate saturated quantities 
        trange=None, 
        # Research
        folderIsExperiment=False, 
        experiment_parameter="-", 
        # Plotting
        logy=False):  
    
    # Create <experiments> and <simulations> based on the given <folder>
    knob = standardParameters[experiment_parameter]["knob"]; key = standardParameters[experiment_parameter]["key"]
    research = create_research(folders=folder, knob1=knob, key1=key, resolutionScan=False, ignore_resolution=True, folderIsExperiment=folderIsExperiment)
    research.update_timeFrame(trange)    
    experiments = research.experiments
    
    # Still implement the real and imaginary parts
    if ("real" in y_quantity) or ("imag" in y_quantity) or ("2" not in y_quantity):
        exit_program("Not implemented yet.", plot_moment_vs_time, sys._getframe().f_lineno)
    
    # Create a figure for each experiment
    for experiment in experiments: 
     
        # Create a figure 
        if y_quantity=="dens2": title = "Time evolution of the density fluctuations squared"
        if y_quantity=="temp2": title = "Time evolution of the temperature fluctuations squared"
        if y_quantity=="upar2": title = "Time evolution of the parallel velocity fluctuations squared"
        if research.numberOfExperiments>1: title += ": "+experiment.line_label
        ax = create_figure(title=title)   
        
        # Plot qmoments(time)
        research.experiments = [experiment]; research.numberOfSimulations = len(experiment.simulations)
        if research.numberOfSimulations==1: subplot_moment_vs_time_bothspecies(ax, research, y_quantity)
        if research.numberOfSimulations!=1: subplot_moment_vs_time(ax, research, y_quantity=y_quantity, specie=specie)
        
        # Appearance    
        ax.yaxis.labelpad = 15; ax.xaxis.labelpad = 10
        if logy: ax.set_yscale("log"); ax.set_ylim([10e-12, 100])
            
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0:  plt.show()
    return

#----------------------------------------
def subplot_moment_vs_time(
        ax, 
        research, 
        # Quantities to be plotted
        y_quantity="dens2", 
        specie=0):
    
    # Colors based on the number of simulations
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,research.numberOfSimulations)) 
    max_time = 0; min_moments = 0; max_moments = 0; count = 0
 
    # Plot qmoments(t) for one of the species
    for experiment in research.experiments:
        for simulation in experiment.simulations: 
        
            # Get the data
            if y_quantity=="dens2": moments = simulation.moments.dens2_vs_ts.dens2[:,specie]  
            if y_quantity=="temp2": moments = simulation.moments.temp2_vs_ts.temp2[:,specie]  
            if y_quantity=="upar2": moments = simulation.moments.upar2_vs_ts.upar2[:,specie]    
            time = simulation.momentses.qmoments_vs_ts.t
            
            # Get the saturated value
            tstart = simulation.time.tstart
            tend = simulation.time.tend
            satmoments = np.mean(moments[(time>=tstart) & (time<=tend)]) 
            
            # Get the label
            s = recognize_species(research, specie)
            label = "|\\delta n|^2" if y_quantity=="dens2" else ("|\\delta T|^2" if y_quantity=="temp2" else "|\\delta v_\\parallel|^2")
            label = "$"+label+"_{\\text{sat},"+s+"}/"+label+"_{\\text{gB},i} = "+str("{:.2f}".format(satmoments))+"$"
            label = simulation.line_label+"; "+label if research.numberOfSimulations!=research.numberOfExperiments else experiment.line_label+"; "+label  
            
            # Plot the time evolution and the saturated value
            ax.plot(time, moments, color=colors[count])  
            ax.plot([tstart, tend], [satmoments,satmoments], color=colors[count], label=label) 
            
            # Keep track of the axis limits
            max_time = np.max([max_time, np.max(time)])
            min_moments = np.min([min_moments, np.min(moments[time>tstart])*1.2])
            max_moments = np.max([max_moments, np.max(moments[time>tstart])*1.2])
            count += 1

    # Correct the y-label
    s = recognize_species(research, specie) 
    ax.set_xlabel(standardLabels["normalized"]["t"])
    ax.set_ylabel(standardLabels["normalized"][y_quantity].replace("_{s}", "_{"+s+"}"))
    
    # Add the legend
    ax.legend(labelspacing=0.0, prop={'size':20}, handlelength=0.8, handletextpad=0.5)
    
    # Axis limits
    ax.set_xlim(left=0, right=max_time)  
    ax.set_ylim([min_moments, max_moments])
    return 

#----------------------------------------
def subplot_moment_vs_time_bothspecies(ax, research, y_quantity="qmoments"):
    """ If we have a single simulation, plot both the ion and electron moments. """
    
    # Get the first simulation
    simulation = research.experiments[0].simulations[0]
        
    # Get the data
    if y_quantity=="dens2": moments = simulation.moments.dens2_vs_ts.dens2
    if y_quantity=="temp2": moments = simulation.moments.temp2_vs_ts.temp2 
    if y_quantity=="upar2": moments = simulation.moments.upar2_vs_ts.upar2  
    time = simulation.moments.dens2_vs_ts.t
    
    # Get the saturated value
    tstart = simulation.time.tstart
    tend = simulation.time.tend
    satmoments = np.mean(moments[(time>=tstart) & (time<=tend),:], axis=0)
    
    # Colors per species
    colors = ["navy", "crimson", "black", "green", "orange"]
    species = range(simulation.dim.species)
    
    # Iterate over the species
    for ispecie, specie in enumerate(species):  
        
        # Get the label
        s = recognize_species(research, specie) 
        label = "|\\delta n|^2" if y_quantity=="dens2" else ("|\\delta T|^2" if y_quantity=="temp2" else "|\\delta v_\\parallel|^2")
        label = "$"+label+"_{\\text{sat},"+s+"}/"+label+"_{\\text{gB},i} = "+str("{:.2f}".format(satmoments[specie]))+"$"
        
        # Plot the time evolution and saturated value
        ax.plot(time, moments[:,specie], color=colors[ispecie])  
        ax.plot([tstart, tend], [satmoments[specie],satmoments[specie]], color=colors[ispecie], label=label) 

    # Axis labels
    ax.set_xlabel(standardLabels["normalized"]["t"])
    ax.set_ylabel(standardLabels["normalized"][y_quantity])

    # Axis limits
    ax.set_xlim(left=0, right=np.max(time))  
    ax.set_ylim(bottom=np.min([0, 1.2*np.min(moments[time>tstart,:])]), top=np.max([0, 1.2*np.max(moments[time>tstart,:])]))
    
    # Add the legend
    ax.legend(labelspacing=0.0, prop={'size':20}, handlelength=0.8, handletextpad=0.5)
    return

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface
    bash = Bash(plot_moment_vs_time, __doc__)  
    
    # Quantities to be plotted
    bash.add_option('y_quantity', 'str', 'y', 'Choose the y-quantity from {qmoments, pmoments, vmoments}.') 
    bash.add_toggle('y_quantity', 'dens2', 'Plot the time evolution of the heat moments.') 
    bash.add_toggle('y_quantity', 'temp2', 'Plot the time evolution of the particle moments.') 
    bash.add_toggle('y_quantity', 'upar2', 'Plot the time evolution of the momentum moments.') 
    bash.add_option('specie', 'int', 's', 'Select the species as "-s 0" for ions and "-s 1" for electrons.') 
    
    # Time range  
    bash.add_option('trange', 'range', 't', 'Select the time frame for the simulations through e.g. "-t [500, 1000]".') 
    
    # Research      
    bash.add_option('folderIsExperiment', 'True', 'f', 'Each folder is an experiment.')  
    bash.add_option('experiment_parameter', 'str', 'e', """Create an <experiment> for each unique value of <experiment_parameter1> in 
    {rho, tiprim, teprim, fprim, tite, delt, nmu, nvgrid, rk, \n\
    dvpa, dmu, nz, nzed, nzgrid, nx, ny, y0, kx max, ky max\n\
    dkx, dky, Lx, Ly, nfield, pol.turns, nperiod, D_hyper}""") 
    
    # Plotting
    bash.add_option('logy', 'True', 'y', 'Logarithmic scale for the y-axis.') 
    
    # Launch the script
    plot_moment_vs_time(**bash.get_arguments())   


