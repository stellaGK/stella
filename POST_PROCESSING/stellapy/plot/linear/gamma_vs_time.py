"""

#===============================================================================
#                       Plot gamma(time) or omega(time)                        #
#=============================================================================== 

For each <simulation> in <research>, created for <folder>, the time evolution
of the growth rate will be shown. If <y_quantity> is "omega" then the time
evolution of the frequency will be shown for each (kx,ky) mode. 

Arguments
---------
    y_quantity : {omega, gamma}
    modes_id : {unstable, stable, all}
    ignore_resolution : {True, False}
    kx_range : [-9999, 9999]
    ky_range : [-9999, 9999]
    
Hanne Thienpondt
19/12/2022

"""

#!/usr/bin/python3
import pathlib
import sys, os
import numpy as np  
import matplotlib as mpl 
import matplotlib.pyplot as plt    

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)    
from stellapy.data.lineardata.get_filterForSelectedModes import get_filterForSelectedModes
from stellapy.data.input.read_inputFile import read_nonlinearFullFluxSurfaceFromInputFile
from stellapy.utils.files.get_firstInputFile import get_firstInputFile 
from stellapy.plot.utils.labels.standardLabels import standardLabels   
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.plot.utils.style.create_figure import create_figure  
from stellapy.simulations.Research import create_research  
from stellapy.utils.commandprompt.bash import Bash
from stellapy.plot.utils.style import Axis, Plot

#===============================================================================
#                       Plot gamma(time) or omega(time)                        #
#=============================================================================== 

def plot_gamma_vs_time(
        folder, 
        # Quantities to be plotted
        y_quantity="gamma", 
        # Selected modes
        modes_id="unstable", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999], 
        # Research options
        folderIsExperiment=False,
        ignore_resolution=True,
        # Plotting options
        number_of_decimal_places=2):
    
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder, ignore_resolution=ignore_resolution, folderIsExperiment=folderIsExperiment) 
    simulations = [s for experiment in research.experiments for s in experiment.simulations] 

    # Create a figure for each simulation
    for simulation in simulations:  
            
        # Create a figure
        if y_quantity=="gamma": title = "Time evolution of the growth rate"
        if y_quantity=="omega": title = "Time evolution of the frequency"
        if research.numberOfSimulations>1: title += ": "+simulation.line_label
        if "k_y" in title: title = title.split("$k_y")[0]+";".join(title.split("$k_y")[-1].split(";")[1:])
        ax = create_figure(title=title)  
         
        # Add the data to the plot
        subplot_gamma_vs_time(ax, simulation, y_quantity, modes_id=modes_id, 
            kx_range=kx_range, ky_range=ky_range, number_of_decimal_places=number_of_decimal_places)
         
        # Change the appearance
        ax.yaxis.labelpad = 15
                
    # Show the figure  
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0: plt.show()
    return

#-----------------------------
def subplot_gamma_vs_time(ax, simulation, y_quantity="gamma", modes_id="unstable", 
        kx_range=[-9999,9999], ky_range=[-9999,9999], fontsize=20, number_of_decimal_places=2): 
    
    # Automate the axis limits  
    axis = Axis(ax, Plot(),  xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.1, percentage=0.9)  
                    
    # Color map for modes    
    selected_modes = get_filterForSelectedModes(simulation, modes_id, kx_range, ky_range) 
    number_of_modes = np.sum(selected_modes.astype(int), axis=(0,1)); plot_i=0
    colors = plt.colormaps.get_cmap('jet')(np.linspace(0,1,number_of_modes))

    # Iterate over the modes  
    for ikx in range(simulation.dim.kx):
        for iky in range(simulation.dim.ky):
            
            # Check whether it needs to be plotted
            if selected_modes[ikx, iky]==False: continue
        
            # Label for the mode  
            ky = simulation.vec.ky[iky]
            string_ky = "{:.2f}".format(ky) if round(ky,2)==float(ky) else ("{:."+str(number_of_decimal_places)+"f}").format(ky)
            label = "$k_y\\rho_i = "+ string_ky + "$"
     
            # Read gamma(t) or omega(t)   
            if y_quantity=="gamma": y = simulation.omega.gamma_vs_tkxky.gamma[:,ikx,iky]
            if y_quantity=="omega": y = simulation.omega.omega_vs_tkxky.omega[:,ikx,iky]
            
            # The time axis is the same for all modes for a <range> simulation 
            # and it is unique for each (kx,ky) if we launched one mode per simulation
            try: x = simulation.omega.gamma_vs_tkxky.t[:,ikx,iky]
            except: x = simulation.omega.gamma_vs_tkxky.t[:]
            
            # Remove nans from gamma(t) or omega(t)   
            y = y[~np.isnan(x)]
            x = x[~np.isnan(x)]

            # Plot gamma(t) or omega(t) and keep track of the axis limits 
            ax.plot(x, y, color=colors[plot_i], label=label)
            axis.update_axisLimits(x, y[~np.isnan(y)]); plot_i+=1 
            
            
    # Only continue if we plotted data
    if len(ax.get_lines())!=0: 
            
        # Print a warning for broken old stella
        if np.max(y)>10e10: 
            print("\n", "".center(80,"="))
            print("", "WARNING".center(80," ")); 
            print("", "".center(80,"="))
            print("  In stella versions from before May 2022, the writing of omega to")
            print("  the netcdf file was broken, it was written as a masked array and only")
            print("  the fill value of the mask can be read. Please rerun these simulations")
            print("  with an up-to-date stella version or read the omega text file instead.\n")
        
        # Labels  
        ax.set_xlabel(standardLabels["normalized"]["t"]) 
        ax.set_ylabel(standardLabels["normalized"][y_quantity]) 
        
        # Legend
        ax.legend(loc="upper right", labelspacing=0.1, handlelength=1, prop={'size':fontsize}, ncol=np.max([int(number_of_modes/10),1]))
        
        # Limits
        axis.rescale_axis()
        
    return  

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Check which script to launch      
    input_file = get_firstInputFile(pathlib.Path(os.getcwd())) 
    nonlinear, full_flux_surface = read_nonlinearFullFluxSurfaceFromInputFile(input_file)   
    if full_flux_surface: os.system("python3 $STELLAPY/plot/linearFFS/gamma_vs_time.py "+" ".join(sys.argv[1:])); sys.exit()
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
    
          
