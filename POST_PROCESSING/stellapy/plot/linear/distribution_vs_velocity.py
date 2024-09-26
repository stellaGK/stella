"""
 
#===============================================================================
#                Plot g2_s(vpa) and g2_s(mu) for each species s                #
#===============================================================================
 
Based on <folder> create a <research> and create a figure with 2*nspec subplots:
    - ax1*nspec: g2_s(kx)     (first row)
    - ax2*nspec: g2_s(ky)     (second row)
    
Arguments
--------- 
    folder : pathlib.Path directory where the command has been executed  
    modes_id : {unstable, stable, all} 
    kx_range : [-9999, 9999]
    ky_range : [-9999, 9999] 
    ignore_resolution : {True, False}
    folderIsExperiment : {False, True}
    
kx_range; ky_range
------------------
Only plot the modes for which kx*rhoi (or ky*rhoi) falls within the range [a,b].
    >> plot_distribution_vs_velocity --kymin 0.5 --kxmax 2.0
    
modes_id
--------
Only plot the "unstable" or "stable" modes, or plot all the modes with "all". 
Modes are identified as "unstable" if (1) <gamma> does not change sign in the last
20% of the time trace and (2) if the sign is positive, since for a "stable" mode 
<gamma> will consist of noise which changes it sign constantly, or <gamma> is 
negative when the potential is damped. However, if the time trace is not long
enough, modes will be mistakenly identified as stable, while they are simply 
unconverged, so make sure the time traces are sufficient. 
    >> plot_distribution_vs_velocity -m unstable
    
ignore_resolution
-----------------
If we launch one-mode-per-simulation, then we can create a smooth spectrum by
using a smaller dt for low ky, and a higher nmu for large ky, in order to 
create one spectrum we need to ignore the differences in resolution. In contrast, 
to take into account the differences in resolution, turn this off through:
    >> plot_distribution_vs_velocity -i
     
folderIsExperiment
------------------
To ignore the automatic sorting of experiment and simulation, you can force stellapy
to create an <experiment> for each <subfolder>, the folder name will be the label
    >> plot_distribution_vs_velocity -f (or --folderIsExperiment) 
 
Hanne Thienpondt 
19/12/2022
 
"""
 
#!/usr/bin/python3  
import sys, os
import pathlib
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
 
# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.data.lineardata.get_filterForSelectedModes import get_filterForSelectedModes 
from stellapy.plot.utils.species.recognize_species import recognize_species  
from stellapy.plot.utils.style.create_figure import update_figure_style
from stellapy.plot.utils.labels.standardLabels import standardLabels  
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.simulations.Research import create_research  
from stellapy.utils.commandprompt.bash import Bash
from stellapy.plot.utils.style import Axis, Plot
 
#===============================================================================
#                Plot g2_s(mu) and g2_s(vpa) for each species s                  #
#===============================================================================
 
def plot_distribution_vs_velocity(
        folder, 
        # Selected modes
        modes_id="unstable", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999], 
        # Research options
        ignore_resolution=True, 
        folderIsExperiment=False,
        # Plotting options 
        number_of_decimal_places=2): 
     
    # Create <experiments> and <simulations> based on the given <folder> 
    research = create_research(folders=folder, folderIsExperiment=folderIsExperiment, ignore_resolution=ignore_resolution)
     
    # Create a figure for each simulation
    for experiment in research.experiments:
        for simulation in experiment.simulations:   
 
            # Create a figure   
            fig = plt.figure(figsize=(18, 9)); axes = []
            grid_specifications = gridspec.GridSpec(2, simulation.input.nspec)
            if simulation.input.nspec==2: grid_specifications.update(top=0.93, left=0.08, right=0.85, bottom=0.1, wspace=0.2, hspace=0.3)
            if simulation.input.nspec==3: grid_specifications.update(top=0.92, left=0.08, right=0.85, bottom=0.1, wspace=0.3, hspace=0.3)
            for i in range(simulation.input.nspec*2): axes.append(plt.subplot(grid_specifications[i]))
             
            # Add a title and make the figure pretty
            title = "Parallel and perpendicular velocity contributions of the distribution squared"
            if research.numberOfSimulations>1: title += ": "+simulation.line_label
            if "k_y" in title: title = title.split("$k_y")[0]+";".join(title.split("$k_y")[-1].split(";")[1:])
            update_figure_style(fig, axes); fig.suptitle(title)
             
            # Add the data to the plot
            for i in range(simulation.input.nspec):
                subplot_distribution_vs_velocity(axes[i], simulation, x_quantity="mu", specie=i, modes_id=modes_id, kx_range=kx_range, ky_range=ky_range, add_legend=False, number_of_decimal_places=number_of_decimal_places)
                subplot_distribution_vs_velocity(axes[i+simulation.input.nspec], simulation, x_quantity="vpa", specie=i, modes_id=modes_id, kx_range=kx_range, ky_range=ky_range, add_legend=False, number_of_decimal_places=number_of_decimal_places)
  
            # Appearance
            for ax in axes: ax.ticklabel_format(style='plain', useOffset=False) if (axes[0].get_ylim()[1]<1000) else ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
            for ax in axes: ax.xaxis.labelpad = 10
            for ax in axes: ax.yaxis.labelpad = 15
            
            # Put the legend to the left
            axes[simulation.input.nspec-1].legend(loc="upper left", bbox_to_anchor=(1.03, 1.05), labelspacing=0.1, handlelength=1, prop={'size':20})
     
    # Show the figure  
    mpl.rcParams["savefig.directory"] = folder   
    if len(axes[0].get_lines())!=0: plt.show() 
    return 
 
#-----------------------------
def subplot_distribution_vs_velocity(
        ax, simulation,
        # Quantity to be plotted 
        x_quantity="vpa", 
        specie=0, 
        # Selected modes
        modes_id="unstable", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999],
        # Plotting
        number_of_decimal_places=2,
        add_legend=True,
        fontsize=18):
     
    # Automate the axis limits and legend 
    axis = Axis(ax, Plot(), xbot_pos=0, ytop_neg=0, overshoot_y=1.1)   
             
    # Color map for modes   
    selected_modes = get_filterForSelectedModes(simulation, modes_id, kx_range, ky_range) 
    number_of_modes = np.sum(selected_modes.astype(int), axis=(0,1)) 
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,number_of_modes)); plot_i=0 
             
    # Get the distribution versus mu or vpa
    x, y_vs_xkxky = get_distribution_data(simulation, x_quantity, specie) 
     
    # Iterate over the modes 
    for ikx in range(simulation.dim.kx):
        for iky in range(simulation.dim.ky): 
     
            # Check whether it needs to be plotted
            if selected_modes[ikx, iky]==False: continue
             
            # Label for the mode  
            ky = simulation.vec.ky[iky]
            string_ky = "{:.2f}".format(ky) if round(ky,2)==float(ky) else ("{:."+str(number_of_decimal_places)+"f}").format(ky)
            label = "$k_y\\rho_i = "+ string_ky + "$"
              
            # Plot g2(mu) or g2(vpa) normalized to one  
            y_vs_xkxky[:,ikx,iky] = y_vs_xkxky[:,ikx,iky]/np.max(np.abs(y_vs_xkxky[:,ikx,iky]))
            ax.plot(x, y_vs_xkxky[:,ikx,iky], marker="o", label=label, color=colors[plot_i]); plot_i += 1 
         
            # Keep track of the axis limits
            axis.update_axisLimits(x, y_vs_xkxky[:,ikx,iky])
             
    # Labels  
    specie_label = recognize_species(simulation, specie) 
    ylabel = standardLabels["normalized"]["g"].replace("{s}",specie_label)
    if x_quantity=="mu": ylabel = "$\\langle$" + ylabel + "$\\rangle_{z, v_\\parallel, t_\\text{last}}$"
    if x_quantity=="vpa": ylabel = "$\\langle$" + ylabel + "$\\rangle_{z, \\mu, t_\\text{last}}$"
    ax.set_xlabel(standardLabels["normalized"][x_quantity]) 
    ax.set_ylabel(ylabel)  
 
    # Limits and legends 
    if add_legend: ax.legend(loc="upper right", labelspacing=0.1, handlelength=1, prop={'size':fontsize}, ncol=np.max([int(number_of_modes/10),1]))
    axis.rescale_axis() 
    ax.set_ylim([0,1]) 
    return 
 
#-------------------------------
def get_distribution_data(simulation, x_quantity, specie):
     
    # Data for the x-axis
    if x_quantity=="mu":  x = simulation.vec.mu 
    if x_quantity=="vpa": x = simulation.vec.vpa 
     
    # Get the data for the y-axis 
    if x_quantity=="mu":  y = simulation.distribution.g2_vs_smukxky.g2[specie,:,:,:]
    if x_quantity=="vpa": y = simulation.distribution.g2_vs_svpakxky.g2[specie,:,:,:]
    return x, y
 
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
  
if __name__ == "__main__":
     
    # Create a bash interface
    bash = Bash(plot_distribution_vs_velocity, __doc__)      
      
    # Adjust the range of wavenumbers
    bash.add_optionheader("wavenumbers")
    bash.add_option('kxmin', 'float', '', '', 'Minimum kx.')   
    bash.add_option('kxmax', 'float', '', '', 'Maximum kx.')  
    bash.add_option('kymin', 'float', '', '', 'Minimum ky.')   
    bash.add_option('kymax', 'float', 'k', '', 'Maximum ky.')  
    
    # Choose the modes
    bash.add_toggleheader("modes")
    bash.add_toggle('modes_id', 'all', '', '', 'Plot all the modes.')  
    bash.add_toggle('modes_id', 'stable', 's', '', 'Plot the stable modes or unconverged modes.')  
    bash.add_toggle('modes_id', 'unstable', '', '', 'Plot the unstable converged modes. (DEFAULT)')   
    bash.add_togglespace()  
    
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
      
    # Get the arguments and execute the script
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_distribution_vs_velocity(**args)   
      
      
 
 


