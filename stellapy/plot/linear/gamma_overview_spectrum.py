"""

#===============================================================================
#               Plot Gamma(ky); Omega(ky); Gamma(t) and dphiz(t)               #
#===============================================================================

Based on <folder> create a <research> and create a figure with 4 subplots:
    - ax1: gamma(ky) based on <subplot_gamma_vs_ky>
    - ax2: omega(ky) based on <subplot_gamma_vs_ky --omega>
    - ax3: gamma(t) based on <subplot_gamma_vs_time>
    - ax4: dphiz(t) based on <subplot_gamma_vs_time>

Hanne Thienpondt
01/09/2022

"""

#!/usr/bin/python3  
import sys, os   
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Stellapy package
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])      
from stellapy.plot.utils.style.create_figure import update_figure_style
from stellapy.plot.linear.gamma_vs_time import subplot_gamma_vs_time
from stellapy.plot.linear.dphiz_vs_time import subplot_dphiz_vs_time
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.plot.linear.gamma_vs_ky import subplot_gamma_vs_ky
from stellapy.simulations.Research import create_research 
from stellapy.plot.utils.labels import standardLabels
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#               Plot Gamma(ky); Omega(ky); Gamma(t) and dphiz(t)               #
#===============================================================================

def gamma_overview_spectrum(folder, kx_range=[-999,999], ky_range=[-999,999], 
        modes_id="unstable", ignoreResolution=True, folderIsExperiment=False):  
  
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder, ignoreResolution=ignoreResolution, folderIsExperiment=folderIsExperiment)
    
    # Create a figure
    fig = plt.figure(figsize=(18, 9))
    grid_specifications = gridspec.GridSpec(2, 2)
    grid_specifications.update(top=0.95, left=0.06, right=0.95, bottom=0.1, hspace=0.25)
    ax1 = plt.subplot(grid_specifications[0])     
    ax2 = plt.subplot(grid_specifications[1])     
    ax3 = plt.subplot(grid_specifications[2])     
    ax4 = plt.subplot(grid_specifications[3])   
    update_figure_style(fig, [ax1, ax2, ax3, ax4])
    
    # Plot gamma(ky) and omega(ky)
    subplot_gamma_vs_ky(ax1, research, "ky", "gamma", modes_id, kx_range, ky_range) 
    subplot_gamma_vs_ky(ax2, research, "ky", "omega", modes_id, kx_range, ky_range) 
    
    # Plot gamma(t) for the first simulation for each ky 
    subplot_gamma_vs_time(ax3, research, "gamma", modes_id, fontsize=12, kx_range=kx_range, ky_range=ky_range)
    
    # Plot dphiz(t)
    subplot_dphiz_vs_time(ax4, research, kx_range, ky_range, fontsize=14, verbose=False) 
    
    # Set the labels
    ax3.set_xlabel(standardLabels["normalized"]["t"])
    ax3.set_ylabel(standardLabels["normalized"]["gamma"])
    
    # Show the figure  
    mpl.rcParams["savefig.directory"] = folder   
    plt.show() 
    return 

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface
    bash = Bash(gamma_overview_spectrum, __doc__)  
    
    # Choose whether we plot the stable, unstable or all modes
    bash.add_option('modes_id', 'str', 'm', 'Choose {unstable, stable, all}.')  
    
    # Add a quick option -s to switch to the stable modes
    bash.add_option('stable', 'True', 's', 'Plot the stable modes.')  
    
    # By default the resolution is ignored, include resolution through -r
    bash.add_option('ignoreResolution', 'False', 'r', 'Do not ignore the resolution.')  
    bash.add_option('folderIsExperiment', 'True', 'f', 'Create an experiment for each folder.')  
    
    # Adjust the range of wavenumbers
    bash.add_option('kxmin', 'float', '-', 'Minimum kx.')   
    bash.add_option('kxmax', 'float', '-', 'Maximum kx.')  
    bash.add_option('kymin', 'float', '-', 'Minimum ky.')   
    bash.add_option('kymax', 'float', 'k', 'Maximum ky.')   
       
    # Get the bash arguments and execute the script
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    gamma_overview_spectrum(**args)   

     
