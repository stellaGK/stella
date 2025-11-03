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
import pathlib
import sys, os
import numpy as np  
import matplotlib as mpl 
import matplotlib.pyplot as plt    

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.input.read_inputFile import read_nonlinearFullFluxSurfaceFromInputFile
from stellapy.utils.files.get_firstInputFile import get_firstInputFile 
from stellapy.plot.utils.labels.standardLabels import standardLabels   
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.utils.get_modes import get_modes
from stellapy.simulations.Research import create_research  
from stellapy.utils.commandprompt.bash import Bash
from stellapy.plot.utils.style import Axis, Plot

#===============================================================================
#                       Plot gamma(time) or omega(time)                        #
#=============================================================================== 

def plot_gamma_vs_time(folder, y_quantity="gamma", modes_id="unstable", 
        kx_range=[-999, 999], ky_range=[-999,999], number_of_plotted_modes=10):
    
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
    subplot_gamma_vs_time(ax, research, y_quantity, modes_id, kx_range, ky_range, number_of_plotted_modes)
    
    # Change the appearance
    ax.yaxis.labelpad = 15
    
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    plt.show()
    return
    
#-----------------------------
def subplot_gamma_vs_time(ax, research, y_quantity="gamma", modes_id="unstable", 
        kx_range=[-999, 999], ky_range=[-999,999], number_of_plotted_modes=10, fontsize=20):
    
    # Only plot this for the first simulation
    simulation = research.experiments[0].simulations[0] 
     
    # Axis labels
    ax.set_xlabel(standardLabels["normalized"]["t"]) 
    ax.set_ylabel(standardLabels["normalized"][y_quantity]) 
    
    # Automate the axis limits  
    axis = Axis(ax, Plot(),  xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.1, percentage=0.9)
                
    # Add the experiment label  
    if research.numberOfSimulations>1: ax.plot(-1, -1, label=simulation.line_label, color="navy") 
    
    # Flux tube simulations 
    if not simulation.full_flux_surface:
    
        # Don't show all modes
        simulation.plotted_modes = get_modes(simulation, kx_range, ky_range, modes_id) 
        simulation.plotted_modes = simulation.plotted_modes[::int(np.ceil(len(simulation.plotted_modes)/number_of_plotted_modes))] 
                 
        # Color map for modes 
        colors_modes = plt.colormaps.get_cmap('jet')(np.linspace(0,1,len(simulation.plotted_modes))); plot_i=0
        
        # Iterate over the modes
        for i, mode in enumerate(simulation.plotted_modes): 
            
            # Label for the mode 
            string_ky = "{:.2f}".format(mode.ky) 
            label = "$k_y\\rho_i = "+ string_ky + "$"
    
            # Plot gamma(t) or omega(t) and keep track of the axis limits 
            if y_quantity=="gamma": y = mode.omega.gamma_vs_t.gamma; x = mode.omega.gamma_vs_t.t
            if y_quantity=="omega": y = mode.omega.omega_vs_t.omega; x = mode.omega.omega_vs_t.t
            ax.plot(x, y, color=colors_modes[plot_i], label=label); plot_i+=1 
            if i>len(simulation.plotted_modes)-2: axis.update_axisLimits(x, y[~np.isnan(y)])
            
    # Full flux surface simulations 
    if simulation.full_flux_surface: 
        
        # Color map for modes 
        dim_kx = simulation.dim.kx; dim_ky = simulation.dim.ky
        colors = plt.colormaps.get_cmap('jet')(np.linspace(0,1,dim_kx*dim_ky)); plot_i=0
        
        # Iterate over the modes
        for ikx in range(dim_kx):
            for iky in range(dim_ky):
        
                # Label for the mode 
                string_kx = "{:.2f}".format(simulation.vec.ky[ikx]) 
                string_ky = "{:.2f}".format(simulation.vec.ky[iky]) 
                label = "$k_y\\rho_i = "+ string_ky + "$"
                if dim_kx>1: label = "$k_x\\rho_i = "+ string_kx + "$; "+label

                # Plot gamma(t) or omega(t) and keep track of the axis limits 
                if y_quantity=="gamma": y = simulation.omega.gamma_vs_tkxky.gamma[:,ikx,iky]; x = simulation.omega.gamma_vs_tkxky.t
                if y_quantity=="omega": y = simulation.omega.omega_vs_tkxky.omega[:,ikx,iky]; x = simulation.omega.omega_vs_tkxky.t
                ax.plot(x, y, color=colors[plot_i], label=label); plot_i+=1 
                if iky>dim_ky-2: axis.update_axisLimits(x, y[~np.isnan(y)])
    
    # Automatically set the axis limits and legends   
    ax.legend(labelspacing=0.1, handlelength=1, prop={'size':fontsize}, ncol=np.max([int(dim_kx*dim_ky/10),1]))
    axis.rescale_axis()
    return  

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Check which script to launch      
    input_file = get_firstInputFile(pathlib.Path(os.getcwd())) 
    nonlinear, full_flux_surface = read_nonlinearFullFluxSurfaceFromInputFile(input_file) 
    if not full_flux_surface: os.system("python3 $STELLAPY/plot/linear/gamma_vs_time.py "+" ".join(sys.argv[1:])); sys.exit() 
    if nonlinear: print("Plot gamma(t) makes no sense for a nonlinear simulation."); sys.exit()
    
    # Create a bash-like interface
    bash = Bash(plot_gamma_vs_time, __doc__)
    
    # Toggle the quantity to be plotted through --omega     
    bash.add_toggleheader("y_quantity")
    bash.add_toggle('y_quantity', 'omega', '', '', 'Plot the time evolution for omega.')
    bash.add_toggle('y_quantity', 'gamma', '', '', 'Plot the time evolution for gamma.')
    bash.add_togglespace()  
    
    # Choose the modes
    bash.add_toggleheader("modes")
    bash.add_toggle('modes_id', 'all', '', '', 'Plot all the modes.')  
    bash.add_toggle('modes_id', 'stable', '', '', 'Plot the stable modes or unconverged modes.')  
    bash.add_toggle('modes_id', 'unstable', '', '', 'Plot the unstable converged modes. (DEFAULT)')   
    bash.add_togglespace()  
      
    # Adjust the range of wavenumbers
    bash.add_optionheader("wavenumbers")
    bash.add_option('kxmin', 'float', '', '', 'Minimum kx.')   
    bash.add_option('kxmax', 'float', '', '', 'Maximum kx.')  
    bash.add_option('kymin', 'float', '', '', 'Minimum ky.')   
    bash.add_option('kymax', 'float', 'k', '', 'Maximum ky.') 
    bash.add_optionspace()
    
    # Other options 
    bash.add_toggleheader("other")
    bash.add_optionheader("other")
    
    # Choose whether we plot the stable, unstable or all modes
    bash.add_option('modes_id', 'str', 'm', '', 'Choose {unstable, stable, all}.') 
    
    # Research options 
    bash.add_toggle('ignore_resolution', False, 'i', 'include_resolution', 'Include the resolution (delta t, nzed, ...) between simulations.')   
    bash.add_toggle('folderIsExperiment', False, 'f', 'folderIsExperiment', 'Create an experiment for each folder')
    
    # Plotting options
    bash.add_option('number_of_decimal_places', 'int', 'd', 'decimal_places', 'Choose the number of decimal places in the ky label.')  
    
    # Get the bash arguments and execute the script 
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_gamma_vs_time(**args) 
    
      
