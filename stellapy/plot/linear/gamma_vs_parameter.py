"""

#===============================================================================
#                            PLOT GAMMA(PARAMETER)                             #
#===============================================================================

A <research> is created for <folder> and for each <simulation> the most unstable
mode is identified and it's growth rate is plotted as a function of the 
selected <parameter> to plot gamma(<x_quantity>).

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
01/09/2022

"""

#!/usr/bin/python3  
import sys, os 
import matplotlib as mpl 
import matplotlib.pyplot as plt 

# Personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0]) 
from stellapy.plot.utils.data.get_gammaOfMostUnstableMode import get_gammaOfMostUnstableMode   
from stellapy.plot.utils.labels import standardParameters, standardLabels  
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy  
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.simulations.Research import create_research 
from stellapy.GUI.plot.utils import Axis, Legend, Plot
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                            PLOT GAMMA(PARAMETER)                             #
#===============================================================================

def plot_gamma_vs_parameter(folder, x_quantity="tiprim", y_quantity="gamma", 
        experiment_parameter1="rho", experiment_parameter2="-",
        kx_range=[-999,999], ky_range=[-999,999]): 
     
    # Find the stella knob and key of the given parameter
    knob1 = standardParameters[experiment_parameter1]["knob"]
    key1 = standardParameters[experiment_parameter1]["key"]
    knob2 = standardParameters[experiment_parameter2]["knob"]
    key2 = standardParameters[experiment_parameter2]["key"]
    
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder, knob1=knob1, key1=key1, knob2=knob2, key2=key2) 
     
    # Create a figure
    if y_quantity=="gamma": title = "Influence of "+standardLabels["normalized"][x_quantity]+" on the growth rate"
    if y_quantity=="omega": title = "Influence of "+standardLabels["normalized"][x_quantity]+" on the frequency"
    ax = create_figure(title=title)  
 
    # Add the data to the plot
    subplot_gamma_vs_parameter(ax, research, x_quantity, y_quantity, kx_range, ky_range)
    
    # Finish figure
    ax.yaxis.labelpad = 15
    ax.xaxis.labelpad = 10
    
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    plt.show()
    return

#-----------------------------
def subplot_gamma_vs_parameter(ax, research, x_quantity="ky", y_quantity="gamma", 
        kx_range=[-999,999], ky_range=[-999,999]): 
    
    # Axis labels
    ax.set_xlabel(standardLabels["normalized"][x_quantity]) 
    ax.set_ylabel(standardLabels["normalized"][y_quantity]) 
     
    # Automate the axis limits and legend
    plot = Plot(); legend = Legend(ax, plot) 
    axis = Axis(ax, plot, xtop_neg=0, xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.1) 
    plot.process_plottingVariables(research)
    
    # Read (ky, gamma, omega) of the most unstable mode for each experiments
    gamma, omega, ky, parameters = get_gammaOfMostUnstableMode(research, x_quantity, kx_range, ky_range)

    # Plot gamma_max(parameter) for each experiment
    for experiment in research.experiments:
        
        # Get the (x,y) data
        x = parameters[experiment.id]
        if y_quantity=="ky": y = ky[experiment.id]
        if y_quantity=="gamma": y = gamma[experiment.id]
        if y_quantity=="omega": y = omega[experiment.id]
        
        # Plot y(x)
        ax.plot(x, y, color=experiment.line_color, label=experiment.line_label, marker=experiment.marker_style) 
        
        # Keep track of the axis limits
        axis.update_axisLimits(x, y)
        
    # Automatically set the axis limits and legends 
    legend.add_legend()
    axis.rescale_axis()
    return

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface
    bash = Bash(plot_gamma_vs_parameter, __doc__)   
    
    # Toggle the quantities to be plotted through --kx --omega  
    bash.add_toggle('y_quantity', 'ky', 'Plot ky(parameter).') 
    bash.add_toggle('y_quantity', 'omega', 'Plot omega(parameter).')    
     
    # Select the x-quantity and the experiment parameter
    bash.add_option('x_quantity', 'str', 'p', '{\
rho, tiprim, teprim, fprim, tite, delt, nmu, nvgrid, rk, \n\
dvpa, dmu, nz, nzed, nzgrid, nx, ny, y0, kx max, ky max\n\
dkx, dky, Lx, Ly, nfield, pol.turns, nperiod, D_hyper}') 
    bash.add_option('experiment_parameter', 'str', 'e', '{\
rho, tiprim, teprim, fprim, tite, delt, nmu, nvgrid, rk, \n\
dvpa, dmu, nz, nzed, nzgrid, nx, ny, y0, kx max, ky max\n\
dkx, dky, Lx, Ly, nfield, pol.turns, nperiod, D_hyper}')     
    
    # Adjust the range of wavenumbers
    bash.add_option('kxmin', 'float', '-', 'Minimum kx.')   
    bash.add_option('kxmax', 'float', '-', 'Maximum kx.')  
    bash.add_option('kymin', 'float', '-', 'Minimum ky.')   
    bash.add_option('kymax', 'float', 'k', 'Maximum ky.')   
    
    # Get the bash arguments and execute the script
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_gamma_vs_parameter(**args)   

     
    
