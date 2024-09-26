"""

#===============================================================================
#             Plot gamma(ky) or omega(ky) or gamma(kx) or omega(kx)            #
#===============================================================================

A <research> is ceated for <folder> and for each <simulation> the spectra gamma(ky); 
omega(ky), gamma(kx) or omega(kx) are plotted based on <x_quantity> and <y_quantity>.

By default, only the unstable modes are plotted, be careful that unconverged modes
can be wrongly identified as stable modes, make sure the time traces are sufficient. 

Arguments
---------
    x_quantity : {kx, ky}
    y_quantity : {omega, gamma, gamma/ky^2}
    modes_id : {unstable, stable, all} 
    ignore_resolution : {True, False}
    kx_range : [-9999, 9999]
    ky_range : [-9999, 9999]

Hanne Thienpondt
19/12/2022

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
from stellapy.data.lineardata.get_filterForSelectedModes import get_filterForSelectedModes
from stellapy.plot.utils.labels.standardLabels import standardLabels  
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.plot.utils.style.create_figure import create_figure 
from stellapy.simulations.Research import create_research   
from stellapy.plot.utils.style import Axis, Legend, Plot
from stellapy.utils.commandprompt.bash import Bash 

#===============================================================================
#             Plot gamma(ky) or omega(ky) or gamma(kx) or omega(kx)            #
#===============================================================================

def plot_gamma_vs_wavenumber(
        folder, 
        # Quantities to be plotted
        x_quantity="ky", 
        y_quantity="gamma", 
        # Selected modes
        modes_id="unstable", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999], 
        # Research options
        folderIsExperiment=False,
        ignore_resolution=True): 
    
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder, ignore_resolution=ignore_resolution, folderIsExperiment=folderIsExperiment)
    
    # Create a figure
    if y_quantity=="gamma/ky^2": title = "Growth rate spectra divided by $k_y^2$"
    if y_quantity=="gamma": title = "Growth rate spectra"
    if y_quantity=="omega": title = "Frequency spectra"
    ax = create_figure(title=title)  
    
    # Add the data to the plot
    subplot_gamma_vs_wavenumber(ax, research, x_quantity, y_quantity, modes_id, kx_range, ky_range)
    
    # Finish figure
    ax.yaxis.labelpad = 15
    ax.xaxis.labelpad = 10
    
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0: plt.show()
    return
    
#-----------------------------
def subplot_gamma_vs_wavenumber(ax, research, x_quantity="ky", y_quantity="gamma", modes_id="unstable", 
        kx_range=[-999,999], ky_range=[-999,999], highlight_most_unstable_mode=True, plotted_simulation=None): 
    
    
    # Automate the axis limits and legend
    plot = Plot(); legend = Legend(ax, plot) 
    axis = Axis(ax, plot, xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.1) 
    plot.process_plottingVariables(research) 

    # Iterate over the simulations and get the (x,y) data
    for experiment in research.experiments:
        for simulation in experiment.simulations: 

            # Sometimes we only want to plot a specific simulation
            if plotted_simulation!=None:
                if simulation!=plotted_simulation: continue

            # Get the modes
            vec_kx = np.array(simulation.vec.kx); dim_kx = simulation.dim.kx
            vec_ky = np.array(simulation.vec.ky); dim_ky = simulation.dim.ky
            
            # Get the lines and marker styles for the experiments and simulations
            style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation)
            if "marker" not in style: style["marker"] = "o"   
            
            # Get omega or gamma
            if y_quantity=="gamma": y = simulation.lineardata.gamma_avg_vs_kxky
            if y_quantity=="omega": y = simulation.lineardata.omega_avg_vs_kxky  
            if y_quantity=="gamma/ky^2": y = simulation.lineardata.gamma_avg_vs_kxky/((vec_ky[np.newaxis,:]**2))
            
            # Get the error on the time trace of gamma(t) or omega(t)
            if y_quantity=="gamma": error = np.array([simulation.lineardata.gamma_min_vs_kxky, simulation.lineardata.gamma_max_vs_kxky])
            if y_quantity=="omega": error = np.array([simulation.lineardata.omega_min_vs_kxky, simulation.lineardata.omega_max_vs_kxky])
            if y_quantity=="gamma/ky^2": error = np.array([np.ones((dim_kx, dim_ky))*np.nan, np.ones((dim_kx, dim_ky))*np.nan])
            
            # Apply the filter of which modes should be shown 
            selected_modes = get_filterForSelectedModes(simulation, modes_id, kx_range, ky_range)
            y[~selected_modes] = np.nan
            error[0, ~selected_modes] = np.nan
            error[1, ~selected_modes] = np.nan 
            
            # Get the wavenumbers
            # TO DO: iterate over the second dimension
            if x_quantity=="kx": x = vec_kx; y = y[:,0]; error = error[:,:,0]
            if x_quantity=="ky": x = vec_ky; y = y[0,:]; error = error[:,0,:] 
            
            # Plot gamma(ky) or omega(ky) or gamma(kx) or omega(kx)    
            ax.errorbar(x, y, yerr=error, **style) 
            axis.update_axisLimits(x, y)  
            
            # Highlight the most unstable mode 
            if highlight_most_unstable_mode:  
                ax.plot(x[np.nanargmax(y)], y[np.nanargmax(y)], mec="black", mfc="gold", marker="*", ms=12, zorder=10)
             
    # Only continue if we plotted data
    if len(ax.get_lines())!=0: 
                
        # Axis labels
        ax.set_xlabel(standardLabels["normalized"][x_quantity]) 
        ax.set_ylabel(standardLabels["normalized"][y_quantity]) 
        
        # Automatically set the axis limits and legends  
        ax.legend(labelspacing=0.0, prop={'size':20}, handlelength=0.8, handletextpad=0.5)
        axis.rescale_axis()
        
    return 

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":

    # Create a bash-like interface
    bash = Bash(plot_gamma_vs_wavenumber, __doc__)   
    
    # Choose the x-quantity
    bash.add_toggleheader("x_quantity")
    bash.add_toggle('x_quantity', 'ky', '', '', 'Plot the spectra as a function of ky.') 
    bash.add_toggle('x_quantity', 'kx', '', '', 'Plot the spectra as a function of kx.') 
    bash.add_togglespace()  
    
    # Choose the y-quantity
    bash.add_toggleheader("y_quantity")
    bash.add_toggle('y_quantity', 'omega', '', '', 'Plot the spectra for omega.')  
    bash.add_toggle('y_quantity', 'gamma', '', '', 'Plot the spectra for gamma.')  
    bash.add_toggle('y_quantity', 'gamma/ky^2', '', '', 'Plot the spectra for gamma/ky**2.')   
    bash.add_togglespace()  
    
    # Choose the modes
    bash.add_toggleheader("modes")
    bash.add_toggle('modes_id', 'all', '', '', 'Plot all the modes.')  
    bash.add_toggle('modes_id', 'stable', 's', '', 'Plot the stable modes or unconverged modes.')  
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
    
    # Research options 
    bash.add_toggle('ignore_resolution', False, 'i', 'include_resolution', 'Include the resolution (delta t, nzed, ...) between simulations.')   
    bash.add_toggle('folderIsExperiment', False, 'f', 'folderIsExperiment', 'Create an experiment for each folder')
    
    # Choose whether we plot the stable, unstable or all modes
    bash.add_option('modes_id', 'str', 'm', '', 'Choose {unstable, stable, all}.')
    
    # Get the bash arguments and execute the script
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_gamma_vs_wavenumber(**args)    
     


