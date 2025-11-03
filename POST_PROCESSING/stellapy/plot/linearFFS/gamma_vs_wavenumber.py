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

Hanne Thienpondt
01/09/2022

"""

#!/usr/bin/python3  
import sys, os   
import pathlib
import numpy as np
import configparser 
import matplotlib as mpl 
import matplotlib.pyplot as plt   

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import get_styleForLinesAndMarkers
from stellapy.plot.utils.labels.standardLabels import standardLabels  
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.Research import create_research   
from stellapy.plot.utils.style import Axis, Legend, Plot
from stellapy.utils.commandprompt.bash import Bash 

#===============================================================================
#             Plot gamma(ky) or omega(ky) or gamma(kx) or omega(kx)            #
#===============================================================================

def plot_gamma_vs_wavenumber(folder, x_quantity="ky", y_quantity="gamma", 
        kx_range=[-999,999], ky_range=[-999,999], modes_id="all"): 
    
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder)
    
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
    plt.show()
    return
    
#-----------------------------
def subplot_gamma_vs_wavenumber(ax, research, x_quantity="ky", y_quantity="gamma",
        modes_id="all", kx_range=[-999,999], ky_range=[-999,999], error_max=0.01, highlight_dominant_gamma=True): 
    
    # Axis labels
    ax.set_xlabel(standardLabels["normalized"][x_quantity]) 
    ax.set_ylabel(standardLabels["normalized"][y_quantity]) 
    
    # Automate the axis limits and legend
    plot = Plot(); legend = Legend(ax, plot) 
    axis = Axis(ax, plot, xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.1)  
        
    # Iterate over the experiments and simulations
    for experiment in research.experiments:
        for simulation in experiment.simulations:
            
            # Change the line label "ny = X" to "ky max = Y"
            if ('ny' in simulation.marker_label):
                simulation.marker_label = "$k_{y,max} = "+str(simulation.input.inputParameters["kxky_grid_box"]["ky max"])+"$" 
                          
            # Get the lines and marker styles for the experiments and simulations
            style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation)
            if "marker" not in style: style["marker"] = "o"  
            
            # Flux tube simulations
            if not simulation.full_flux_surface: 
                               
                # Get either the stable or the unstable modes
                if modes_id=="unstable": modes = [mode for mode in simulation.modes if mode.lineardata.unstable]
                if modes_id=="stable": modes = [mode for mode in simulation.modes if mode.lineardata.stable]
                if modes_id=="all": modes = simulation.modes  
                if len(modes)==0: exit_program("No modes were found to be plotted.", subplot_gamma_vs_wavenumber, sys._getframe().f_lineno) 
                 
                # Get the modes (kx,ky) within kx_range and ky_range 
                modes = [mode for mode in modes if (round(mode.ky,2) >= ky_range[0] and round(mode.ky,2) <= ky_range[1])]
                modes = [mode for mode in modes if (round(mode.kx,2) >= kx_range[0] and round(mode.kx,2) <= kx_range[1])]
                
                # Get the wavenumbers
                if x_quantity=="kx": x = [mode.kx for mode in modes]  
                if x_quantity=="ky": x = [mode.ky for mode in modes]  
                
                # Get omega or gamma
                if y_quantity=="gamma": y = [mode.lineardata.gamma_avg for mode in modes]
                if y_quantity=="omega": y = [mode.lineardata.omega_avg for mode in modes]  
                if y_quantity=="gamma/ky^2": y = [mode.lineardata.gamma_avg/(mode.ky**2) for mode in modes]
                
                # Get the error on the time trace of gamma(t) or omega(t)
                if y_quantity=="gamma": error = [[mode.lineardata.gamma_min, mode.lineardata.gamma_max] for mode in modes]
                if y_quantity=="omega": error = [[mode.lineardata.omega_min, mode.lineardata.omega_max] for mode in modes]
                if y_quantity=="gamma/ky^2": error = [[np.nan, np.nan]  for mode in modes]
            
            # Full flux surface simulations
            if simulation.full_flux_surface: 
                
                # Access the full flux surface gamma to make sure the data is written
                simulation.lineardata.ffs_gamma
                
                # Read the linear data file
                file = configparser.ConfigParser() 
                path = simulation.input_file.with_suffix(".lineardata.ini") 
                file.read(path)  
                
                # Get the wavenumbers
                if x_quantity=="kx": x = simulation.vec.kx
                if x_quantity=="ky": x = simulation.vec.ky
                
                # Store (gamma, omega) in a matrix
                y = np.ones((simulation.dim.kx, simulation.dim.ky))*np.nan
                y_min = np.ones((simulation.dim.kx, simulation.dim.ky))*np.nan
                y_max = np.ones((simulation.dim.kx, simulation.dim.ky))*np.nan
                for ikx in range(simulation.dim.kx):
                    for iky in range(simulation.dim.ky):
                        mode_identifier = '('+str(simulation.vec.kx[ikx])+', '+str(simulation.vec.ky[iky])+')' 
                        y[ikx,iky] = float(file[mode_identifier][y_quantity+"_avg"])
                        y_min[ikx,iky] = float(file[mode_identifier][y_quantity+"_min"])
                        y_max[ikx,iky] = float(file[mode_identifier][y_quantity+"_max"])
                
                # For now force ikx=0
                y = y[0,:]; error = np.array([y_min[0,:], y_max[0,:]]).T 
                                       
            # Plot gamma(ky) or omega(ky) or gamma(kx) or omega(kx) 
            ax.errorbar(x, y, yerr=error.T, **style) 
            
            # Highlight the dominant mode by averaging over all modes right off 
            # the first mode that converged
            if highlight_dominant_gamma and not simulation.full_flux_surface:
                ky = [mode.ky for mode in modes]  
                gamma = [mode.lineardata.gamma_avg for mode in modes]
                error = [[mode.lineardata.gamma_min, mode.lineardata.gamma_max] for mode in modes]
                relative_error = np.array(error)/np.array(gamma)[:,np.newaxis] 
                indices = [i for i in range(len(gamma)) if (relative_error[i][0]<error_max and relative_error[i][1]<error_max and error[i][0]<error_max and error[i][1]<error_max)]
                if len(indices)==0: indices = range(len(gamma)) 
                indices = [i for i in range(len(gamma)) if i>=indices[0]]
                gamma = np.array(gamma)[indices]; error = np.array(error)[indices]; ky = np.array(ky)[indices]
                average_gamma = np.average(gamma)
                ax.plot([-999,999], [average_gamma, average_gamma], lw=3, color="crimson", zorder=10)
                ax.axvline(x=ky[0], color="crimson", ls=":")
                simulation.ffs_gamma = average_gamma
                    
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
    bash = Bash(plot_gamma_vs_wavenumber, __doc__)
    
    # Toggle the quantities to be plotted through --ky --kx --omega --gamma/ky^2
    bash.add_toggle('x_quantity', 'ky', 'Plot the spectra as a function of ky.')
    bash.add_toggle('x_quantity', 'kx', 'Plot the spectra as a function of kx.')
    bash.add_toggle('y_quantity', 'omega', 'Plot the spectra for omega.')
    bash.add_toggle('y_quantity', 'gamma/ky^2', 'Plot the spectra for gamma/ky**2.')
    
    # Choose whether we plot the stable, unstable or all modes
    bash.add_option('modes_id', 'str', 'm', 'Choose {unstable, stable, all}.')
    
    # Add a quick option -s to switch to the stable modes
    bash.add_option('stable', 'True', 's', 'Plot the stable modes.')
    
    # Adjust the range of wavenumbers
    bash.add_option('kxmin', 'float', '-', 'Minimum kx.')
    bash.add_option('kxmax', 'float', '-', 'Maximum kx.')
    bash.add_option('kymin', 'float', '-', 'Minimum ky.')
    bash.add_option('kymax', 'float', 'k', 'Maximum ky.')
    
    # Get the bash arguments and execute the script
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_gamma_vs_wavenumber(**args)
     


