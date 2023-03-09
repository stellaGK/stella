"""
 
#===============================================================================
#                            PLOT GAMMA(PARAMETER)                             #
#===============================================================================
 
A <research> is created for <folder> and for each <simulation> the most unstable
mode is identified and its growth rate is plotted as a function of the 
selected <x_quantity> to plot gamma(<x_quantity>).
 
For each <experiment> there will be one line of gamma(<x_quantity>). In order to 
sort the simulations in experiments, set <experiment_parameter1> and 
<experiment_parameter2>, for each unique value of these parameters, a new 
<experiment> will be created. 
 
Arguments
---------
    x_quantity : {fprim, tiprim, teprim, rho, ...}
    y_quantity : {omega, gamma, gamma/ky^2} 
    experiment_parameter1 : {fprim, tiprim, teprim, rho, ...}
    experiment_parameter2 : {fprim, tiprim, teprim, rho, ...}
 
Hanne Thienpondt
20/10/2022
 
"""
 
#!/usr/bin/python3  
import sys, os 
import pathlib
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt 
 
# Personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.lineardata.load_mostUnstableModeObject import load_mostUnstableModeObject
from stellapy.plot.utils.data.recognize_varied_parameter import recognize_varied_parameter
from stellapy.plot.utils.labels import standardParameters, standardLabels  
from stellapy.plot.utils.data.get_parameters import get_parameters
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy  
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.simulations.Research import create_research 
from stellapy.utils.commandprompt.bash import Bash
from stellapy.plot.utils.style import Axis, Plot
 
#===============================================================================
#                            PLOT GAMMA(PARAMETER)                             #
#===============================================================================
 
def plot_gamma_vs_parameter(
        folder, 
        # Quantities to be plotted
        x_quantity=None, 
        y_quantity="gamma", 
        # Research options
        experiment_parameter1="rho", 
        experiment_parameter2="-", 
        ignore_resolution=True,
        folderIsExperiment=False, 
        # Selected modes
        modes_id="unstable", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999]): 
    
    # Try to recognize the scanned parameter
    x_quantity = recognize_varied_parameter(folder, x_quantity) 
      
    # Find the stella knob and key of the given parameter
    knob1 = standardParameters[experiment_parameter1]["knob"]
    key1 = standardParameters[experiment_parameter1]["key"]
    knob2 = standardParameters[experiment_parameter2]["knob"]
    key2 = standardParameters[experiment_parameter2]["key"]
 
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder, knob1=knob1, key1=key1, knob2=knob2, key2=key2, ignore_resolution=ignore_resolution, folderIsExperiment=folderIsExperiment) 

    # Create a figure
    if y_quantity=="gamma": title = "Influence of "+standardLabels["normalized"][x_quantity]+" on the growth rate"
    if y_quantity=="omega": title = "Influence of "+standardLabels["normalized"][x_quantity]+" on the frequency"
    ax = create_figure(title=title)  
  
    # Add the data to the plot
    subplot_gamma_vs_parameter(ax, research, x_quantity, y_quantity, modes_id, kx_range, ky_range)
     
    # Finish figure
    ax.yaxis.labelpad = 15
    ax.xaxis.labelpad = 10
     
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder 
    if len(ax.get_lines())!=0: plt.show()
    return
 
#-----------------------------
def subplot_gamma_vs_parameter(
        ax, research, 
        # Quantities to be plotted
        x_quantity="tiprim", 
        y_quantity="gamma", 
        # Selected modes
        modes_id="unstable", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999],
        # Plotting options
        fontsize=20): 
      
    # Automate the axis limits and legend 
    axis = Axis(ax, Plot(), ytop_neg=0, ybot_pos=0, overshoot_y=1.1)  
    
    # Identify the stella (knob, key) of the parameter corresponding to <x_quantity>
    knob = standardParameters[x_quantity]["knob"]
    key = standardParameters[x_quantity]["key"] 
    
    # Plot gamma_mostunstablemode(parameter) for each experiment
    for experiment in research.experiments:
        
        # Gather the chosen parameters for the x-axis
        parameters, simulations = get_parameters(experiment, key, knob) 
        x = np.array(parameters)
        
        # Identify the most unstable mode for each simulation 
        for simulation in simulations: 
            load_mostUnstableModeObject(simulation, modes_id, kx_range, ky_range)
        
        # Gather gamma_mostunstablemode(parameter) for each experiment 
        if y_quantity=="ky":     y = np.array([s.mostunstablemode.ky for s in simulations])
        if y_quantity=="gamma":  y = np.array([s.mostunstablemode.gamma for s in simulations])
        if y_quantity=="omega":  y = np.array([s.mostunstablemode.omega for s in simulations])
        
        # The error is represented by the minimum and maximum value of gamma
        # and omega that they obtain in the last 15% of the time trace
        if y_quantity=="ky":     error = np.ones((2,len(y)))*np.nan
        if y_quantity=="gamma":  error = np.array([s.mostunstablemode.gamma_error for s in simulations])
        if y_quantity=="omega":  error = np.array([s.mostunstablemode.omega_error for s in simulations])
        
        # Remove nans 
        error = error[~np.isnan(y),:]
        x = x[~np.isnan(y)]
        y = y[~np.isnan(y)]
        if len(x)==0: continue
        
        # Plot y(x)  
        ax.errorbar(x, y, yerr=np.array(error).T, color=experiment.line_color, label=experiment.line_label, marker=experiment.marker_style) 
 
        # Keep track of the axis limits
        axis.update_axisLimits(x, y)
         
    # Only continue if we plotted data
    if len(ax.get_lines())!=0: 
        
        # Labels
        ax.set_xlabel(standardLabels["normalized"][x_quantity]) 
        ax.set_ylabel(standardLabels["normalized"][y_quantity]) 
             
        # Legend
        ax.legend(labelspacing=0.1, handlelength=1, prop={'size':fontsize})
        
        # Limits 
        axis.rescale_axis()
        
    return
 
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
  
if __name__ == "__main__":
     
    # Create a bash-like interface
    bash = Bash(plot_gamma_vs_parameter, __doc__)   
      
    # Adjust the range of wavenumbers
    bash.add_optionheader("wavenumbers")
    bash.add_option('kxmin', 'float', '', '', 'Minimum kx.')   
    bash.add_option('kxmax', 'float', '', '', 'Maximum kx.')  
    bash.add_option('kymin', 'float', '', '', 'Minimum ky.')   
    bash.add_option('kymax', 'float', 'k', '', 'Maximum ky.') 
    bash.add_optionspace()
    
    # Choose the y-quantity
    bash.add_toggleheader("x-quantity")
    bash.add_toggle('x_quantity', 'fprim', '', '', 'Plot the frequency or growthrate vs the density gradient.')  
    bash.add_toggle('x_quantity', 'tprim', '', 'tprim', 'Plot the frequency or growthrate vs the ion temperature gradient.')  
    bash.add_toggle('x_quantity', 'tiprim', '', '', 'Plot the frequency or growthrate vs the ion temperature gradient.')  
    bash.add_toggle('x_quantity', 'teprim', '', '', 'Plot the frequency or growthrate vs the electron temperature gradient.')  
    bash.add_toggle('x_quantity', 'rho', '', '', 'Plot the frequency or growthrate vs the effective minor radius.')   
    bash.add_togglespace()  
    
    # Choose the y-quantity
    bash.add_toggleheader("y_quantity")
    bash.add_toggle('y_quantity', 'omega', '', '', 'Plot the frequency vs parameter.')  
    bash.add_toggle('y_quantity', 'gamma', '', '', 'Plot the growth rate vs parameter.')   
    bash.add_togglespace()  
    
    # Other options 
    bash.add_toggleheader("other")
    bash.add_optionheader("other")   
      
    # Select the x-quantity and the experiment parameter
    bash.add_option('x_quantity', 'str', 'x', '', '{\
rho, tiprim, teprim, fprim, tite, delt, nmu, nvgrid, rk, \n\
dvpa, dmu, nz, nzed, nzgrid, nx, ny, y0, kx max, ky max\n\
dkx, dky, Lx, Ly, nfield, pol.turns, nperiod, D_hyper}') 
    bash.add_option('experiment_parameter1', 'str', 'e', '', 'See list above.')    
    bash.add_option('experiment_parameter2', 'str', '', '', 'See list above.')    
     
    # Choose whether we plot the stable, unstable or all modes
    bash.add_option('modes_id', 'str', 'm', '', 'Choose {unstable, stable, all}.')
     
    # Research options 
    bash.add_toggle('ignore_resolution', False, 'i', 'include_resolution', 'Include the resolution (delta t, nzed, ...) between simulations.')   
    bash.add_toggle('folderIsExperiment', False, 'f', 'folderIsExperiment', 'Create an experiment for each folder')
     
    # Get the bash arguments and execute the script
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_gamma_vs_parameter(**args)   
 
      
     
