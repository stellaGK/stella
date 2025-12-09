"""



"""


#!/usr/bin/python3 
import pathlib 
import sys, os
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.data.input.read_inputFile import read_nonlinearFullFluxSurfaceFromInputFile
from stellapy.plot.utils.species.recognize_species import recognize_species
from stellapy.utils.files.get_firstInputFile import get_firstInputFile
from stellapy.plot.utils.labels.standardLabels import standardLabels 
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.Research import create_research
from stellapy.utils.commandprompt.bash import Bash
from stellapy.plot.utils.style import Axis, Plot

#===============================================================================
#                               Plot g2(time)                                #
#===============================================================================

def plot_distribution_vs_time(
        folder, 
        # Quanities to be plotted
        y_quantity="g2",
        # Time frame to calculate saturated quantities 
        trange=None, 
        # Research
        folderIsExperiment=False, 
        folderIsSimulation=False,
        # Plot
        xrange=None,
        yrange=None,
        label_units="_norm",
        zonal=None, 
        log=False):  
    
    # Create <simulations> based on the given <folder>
    research = create_research(folders=folder, resolutionScan=False, folderIsExperiment=folderIsExperiment, folderIsSimulation=folderIsSimulation, ignore_resolution=False) 
    simulations = [s for experiment in research.experiments for s in experiment.simulations] 
    research.update_timeFrame(trange); experiments = research.experiments
    
    # Create a figure for each experiment
    for experiment in experiments: 
     
        # Create a figure
        title = "Time evolution of the average distribution squared in real space"
        if research.numberOfExperiments>1: title += ": "+experiment.line_label
        ax = create_figure(title=title)    
         
        # Plot g2(t)
        research.experiments = [experiment]; research.numberOfSimulations = len(experiment.simulations)
        if len(simulations)==1: subplot_distribution_vs_time_allspecies(ax, research, y_quantity=y_quantity, log=log)
        if len(simulations)!=1: subplot_distribution_vs_time(ax, research, y_quantity=y_quantity, zonal=zonal, log=log, label_units=label_units)   
        
        # Appearance
        ax.yaxis.labelpad = 15
        ax.xaxis.labelpad = 10
        if log: ax.set_yscale("log"); ax.set_ylim([10e-9, 1000])         
        if xrange: ax.set_xlim(xrange)    
        if yrange: ax.set_ylim(yrange)  
            
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0: plt.show()
    return

#----------------------------------------
def subplot_distribution_vs_time_allspecies(ax, research, y_quantity="g2", log=False, fontsize=20):

    # Get the first simulation
    simulation = research.experiments[0].simulations[0]
    
    # Get the data
    if y_quantity!="g2": exit_program("Only y_quantity==g2 has been implemented.", subplot_distribution_vs_time_allspecies, sys._getframe().f_lineno)
    g2_vs_ts = simulation.distribution.g2_vs_ts.g2  
    time = simulation.distribution.g2_vs_ts.t
    
    # Get the saturated data
    tstart = simulation.time.tstart; tend = simulation.time.tend
    saturated_distribution = np.mean(g2_vs_ts[(time>=tstart) & (time<=tend)], axis=0)  
    
    # Colors per species
    colors = ["navy", "crimson", "black", "green", "orange"]
    species = range(simulation.dim.species)
    
    # Iterate over the species
    for ispecie, specie in enumerate(species):  
        
        # Get the label
        s = recognize_species(research, specie)  
        label = "$g^2_{\\text{sat,}"+s+"}/g^2_{\\text{gB},i} = "+str("{:.2f}".format(saturated_distribution[specie]))+"$"
        
        # Plot the time evolution and saturated value
        ax.plot(time, g2_vs_ts[:,specie], color=colors[ispecie])  
        ax.plot([tstart, tend], [saturated_distribution[specie],saturated_distribution[specie]], color=colors[ispecie], label=label) 
        
    # Axis limits
    ax.set_xlim(left=0, right=np.max(time))
    g2_vs_ts = g2_vs_ts[int(len(g2_vs_ts)/2):,:]
    if log==True:  ax.set_ylim(bottom=np.min(g2_vs_ts), top=np.max(g2_vs_ts))
    if log==False: ax.set_ylim(bottom=0, top=np.max(g2_vs_ts)*1.1)
    
    # legend and labels
    ax.legend(labelspacing=0.5, handlelength=1, shadow=True, prop={'size':fontsize})
    ax.set_xlabel(standardLabels["normalized"]["t"])
    ax.set_ylabel(standardLabels["normalized"]["g2"])
    return

#----------------------------------------
def subplot_distribution_vs_time(ax, research, y_quantity="g2", zonal=None, log=False, label_units=""):
    
    # Automate the axis limits 
    exit_program("Haven't written this script yet.", subplot_distribution_vs_time, sys._getframe().f_lineno) 
    axis = Axis(ax, Plot(),  xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.1, logy=log, percentage=0.5)
    
    # Colors based on the number of simulations
    colors = plt.colormaps.get_cmap('jet')(np.linspace(0,1,research.numberOfSimulations)); count=0
    
    # Plot g2(t)
    for experiment in research.experiments:
        for simulation in experiment.simulations:  
        
            # Get the data
            if y_quantity!="g2": exit_program("Only y_quantity==g2 has been implemented.", subplot_distribution_vs_time, sys._getframe().f_lineno)
            if zonal==None: g2 = simulation.distribution.g2_vs_t.g2 
            if zonal==True: g2 = simulation.distribution.g2_vs_t_zonal.g2_zonal 
            if zonal==False: g2 = simulation.distribution.g2_vs_t_nozonal.g2_nozonal 
            time = simulation.distribution.g2_vs_t.t 
            
            # Get the label 
            tstart = simulation.time.tstart; tend = simulation.time.tend
            saturated_distribution = np.mean(g2[(time>=tstart) & (time<=tend)]) 
            label = "$\\langle \\tilde{g}_{\\text{sat}}^2 \\rangle_z = "+str("{:.2f}".format(saturated_distribution))+"$"
            label = simulation.line_label+"; "+label if research.numberOfSimulations!=research.numberOfExperiments else experiment.line_label+"; "+label 
            
            # Plot the time evolution of the distribution squared  
            ax.plot(time, g2[:], color=colors[count])  
            ax.plot([tstart, tend], [saturated_distribution, saturated_distribution], color=colors[count], label=label) 
            axis.update_axisLimits(time, g2) ; count +=1  
    
    # Appearance
    if zonal==None: ax.set_ylabel(standardLabels["normalized"]["g2_fieldlineaverage"+label_units])
    if zonal==True: ax.set_ylabel(standardLabels["normalized"]["g2_zonal_fieldlineaverage"+label_units])
    if zonal==False: ax.set_ylabel(standardLabels["normalized"]["g2_nozonal_fieldlineaverage"+label_units])
    ax.legend(labelspacing=0.2, handlelength=1, shadow=True)
    
    # Axis limits
    axis.rescale_axis()
    return
            
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#=============================================================================== 
 
if __name__ == "__main__":
    
    # Check which script to launch      
    input_file = get_firstInputFile(pathlib.Path(os.getcwd()))  
    nonlinear, full_flux_surface = read_nonlinearFullFluxSurfaceFromInputFile(input_file)
    if not nonlinear: os.system("python3 $STELLAPY/plot/linear/distribution_vs_time.py "+" ".join(sys.argv[1:])); sys.exit()
    
    # Create a bash-like interface
    bash = Bash(plot_distribution_vs_time, __doc__)  
    
    # Choose the y-quantity
    bash.add_toggleheader("y_quantity")
    bash.add_toggle('y_quantity', 'g2', '', '', 'Plot the time evolution of the distribution squared.')  
    bash.add_toggle('zonal', True, 'z', 'zonal', 'Plot only the zonal modes.') 
    bash.add_toggle('zonal', False, 'n', 'nozonal', 'Plot only the non-zonal modes.') 
    bash.add_togglespace()  
    
    # Manipulate the data
    bash.add_optionheader("simulation") 
    bash.add_toggleheader("simulation")
    bash.add_toggle('folderIsExperiment', True, 'f', 'folderIsExp', 'Each folder is an experiment.') 
    bash.add_toggle('folderIsSimulation', True, '', 'all', 'Each folder is a simulation.')  
    bash.add_option('trange', 'range', 't', '', "Select the time frame for the saturated state of the simulations, e.g. '-t [500, 1000]'.") 
    bash.add_option('xrange', 'range', '', 'xrange', "Select the xrange for the figure, e.g. '[0,1000]'.") 
    bash.add_option('yrange', 'range', '', 'yrange', "Select the yrange for the figure, e.g. '[0,10]'.") 
    bash.add_optionspace()
    bash.add_togglespace()  
    
    # Change the appearance of the figure
    bash.add_optionheader("Figure") 
    bash.add_toggleheader("Figure") 
    bash.add_toggle('log', True, 'l', 'log', 'Use a logarithmic y-axis.')  
    bash.add_optionspace()
    bash.add_togglespace()
    
    # Other options 
    bash.add_toggleheader("other")
    bash.add_optionheader("other")    
    
    # Get the arguments and execute the script
    plot_distribution_vs_time(**bash.get_arguments())   


