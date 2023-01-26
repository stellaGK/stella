"""

#===============================================================================
#                       Plot gamma(time) or omega(time)                        #
#=============================================================================== 

For the first simulation in <research>, created for <folder>, the time evolution
of the growth rate will be shown. If <y_quantity> is "omega" then the time
evolution of the frequency will be shown for each (kx,ky)  mode. 

Arguments
---------
    y_quantity : {omega, gamma}
    modes_id : {unstable, stable, all}
    
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
from stellapy.plot.utils.labels.standardLabels import standardLabels   
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.utils.get_modes import get_modes
from stellapy.simulations.Research import create_research  
from stellapy.utils.commandprompt.bash import Bash
from stellapy.GUI.plot.utils import Axis, Plot

#===============================================================================
#                       Plot gamma(time) or omega(time)                        #
#=============================================================================== 

def plot_gamma_vs_time(folder, y_quantity="gamma", modes_id="unstable", kx_range=[-999, 999], ky_range=[0,999]):
    
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder)
    if research.numberOfSimulations>1:
        exit_reason = "Stellapy can only plot "+y_quantity+"(t) of one simulation at a time.\n" 
        exit_reason += "Please go in to a subfolder which contains a single simulation." 
        exit_program(exit_reason, plot_gamma_vs_time, sys._getframe().f_lineno)   
        
    # Create a figure
    if y_quantity=="gamma": title = "Time evolution of the growth rate"
    if y_quantity=="omega": title = "Time evolution of the frequency"
    ax = create_figure(title=title)  
    
    # Add the data to the plot
    subplot_gamma_vs_time(ax, research, y_quantity, modes_id, kx_range, ky_range)
    
    # Change the appearance
    ax.yaxis.labelpad = 15
    
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    plt.show()
    return
    
#-----------------------------
def subplot_gamma_vs_time(ax, research, y_quantity="gamma", modes_id="unstable", 
        kx_range=[-999, 999], ky_range=[0,999], fontsize=20, ):
    
    # Only plot this for the first simulation
    simulation = research.experiments[0].simulations[0] 
     
    # Axis labels
    ax.set_xlabel(standardLabels["normalized"]["t"]) 
    ax.set_ylabel(standardLabels["normalized"][y_quantity]) 
    
    # Automate the axis limits  
    axis = Axis(ax, Plot(),  xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.1, percentage=0.9)
                
    # Add the experiment label  
    if research.numberOfSimulations>1: ax.plot(-1, -1, label=simulation.line_label, color="navy") 
            
    # Color map for modes 
    simulation.plotted_modes = get_modes(simulation, kx_range, ky_range, modes_id)
    colors_modes = plt.cm.get_cmap('jet')(np.linspace(0,1,len(simulation.plotted_modes))); plot_i=0
    
    # Iterate over the modes
    for mode in simulation.plotted_modes: 
        
        # Label for the mode 
        string_ky = "{:.2f}".format(mode.ky) if round(mode.ky,1)==float(mode.ky) else str(float(mode.ky))
        label = "$k_y\\rho_i = "+ string_ky + "$"

        # Plot gamma(t) or omega(t) and keep track of the axis limits 
        if y_quantity=="gamma": y = mode.omega.gamma_vs_t.gamma; x = mode.omega.gamma_vs_t.t
        if y_quantity=="omega": y = mode.omega.omega_vs_t.omega; x = mode.omega.omega_vs_t.t
        ax.plot(x, y, color=colors_modes[plot_i], label=label)
        axis.update_axisLimits(x, y[~np.isnan(y)]); plot_i+=1 
    
    # Automatically set the axis limits and legends   
    ax.legend(loc="upper right", labelspacing=0.1, handlelength=1, prop={'size':fontsize}, ncol=np.max([int(len(simulation.plotted_modes)/10),1]))
    axis.rescale_axis()
    return  

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface
    bash = Bash(plot_gamma_vs_time, __doc__)
    
    # Toggle the quantity to be plotted through --omega     
    bash.add_toggle('y_quantity', 'omega', 'Plot the time evolution for omega.')  
    
    # Choose whether we plot the stable, unstable or all modes
    bash.add_option('modes_id', 'str', 'm', 'Choose {unstable, stable, all}.')
      
    # Adjust the range of wavenumbers
    bash.add_option('kxmin', 'float', '-', 'Minimum kx.')   
    bash.add_option('kxmax', 'float', '-', 'Maximum kx.')  
    bash.add_option('kymin', 'float', '-', 'Minimum ky.')   
    bash.add_option('kymax', 'float', 'k', 'Maximum ky.') 
    
    # Get the bash arguments and execute the script 
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_gamma_vs_time(**args) 
    
      
