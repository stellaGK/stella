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
import pathlib
import sys, os   
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)      
from stellapy.data.input.read_inputFile import read_nonlinearFullFluxSurfaceFromInputFile
from stellapy.plot.linear.gamma_vs_wavenumber import subplot_gamma_vs_wavenumber
from stellapy.plot.utils.style.create_figure import update_figure_style
from stellapy.utils.files.get_firstInputFile import get_firstInputFile  
from stellapy.plot.linear.gamma_vs_time import subplot_gamma_vs_time
from stellapy.plot.linear.dphiz_vs_time import subplot_dphiz_vs_time
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.simulations.Research import create_research 
from stellapy.plot.utils.labels import standardLabels
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#               Plot Gamma(ky); Omega(ky); Gamma(t) and dphiz(t)               #
#===============================================================================

def plot_gamma_overview_spectra(
        folder,  
        # Research options 
        ignore_resolution=True,
        folderIsExperiment=False, 
        # Selected modes
        modes_id="unstable", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999]): 
  
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder, ignore_resolution=ignore_resolution, folderIsExperiment=folderIsExperiment)
    simulations = [s for experiment in research.experiments for s in experiment.simulations] 
    
    # Create a figure for each simulation
    for simulation in simulations: 
    
        # Create a figure
        fig = plt.figure(figsize=(18, 9))
        grid_specifications = gridspec.GridSpec(2, 2)
        grid_specifications.update(top=0.92, left=0.06, right=0.95, bottom=0.1, hspace=0.25)
        ax1 = plt.subplot(grid_specifications[0])     
        ax2 = plt.subplot(grid_specifications[1])     
        ax3 = plt.subplot(grid_specifications[2])     
        ax4 = plt.subplot(grid_specifications[3])   
        update_figure_style(fig, [ax1, ax2, ax3, ax4]) 
        if research.numberOfSimulations>1: fig.suptitle(simulation.line_label)
        
        # Plot gamma(ky) and omega(ky) 
        subplot_gamma_vs_wavenumber(ax1, research, "ky", "gamma", modes_id=modes_id, kx_range=kx_range, ky_range=ky_range, plotted_simulation=simulation) 
        subplot_gamma_vs_wavenumber(ax2, research, "ky", "omega", modes_id=modes_id, kx_range=kx_range, ky_range=ky_range, plotted_simulation=simulation)  
        
        # Plot gamma(t) for the first simulation for each ky  
        subplot_gamma_vs_time(ax3, simulation, "gamma", modes_id=modes_id, fontsize=12, kx_range=kx_range, ky_range=ky_range)
        
        # Plot dphiz(t)
        subplot_dphiz_vs_time(ax4, simulation, modes_id=modes_id, kx_range=kx_range, ky_range=ky_range, fontsize=14, verbose=False) 
        
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
    
    # Check which script to launch      
    input_file = get_firstInputFile(pathlib.Path(os.getcwd()))  
    if input_file!=None: 
        nonlinear, full_flux_surface = read_nonlinearFullFluxSurfaceFromInputFile(input_file) 
        if full_flux_surface: os.system("python3 $STELLAPY/plot/linearFFS/gamma_overview_spectrum.py "+" ".join(sys.argv[1:])); sys.exit() 
        if nonlinear: print("Plot gamma(ky) and gamma(t) makes no sense for a nonlinear simulation."); sys.exit()
    
    # Create a bash-like interface
    bash = Bash(plot_gamma_overview_spectra, __doc__)   
      
    # Adjust the range of wavenumbers
    bash.add_optionheader("wavenumbers")
    bash.add_option('kxmin', 'float', '', '', 'Minimum kx.')   
    bash.add_option('kxmax', 'float', '', '', 'Maximum kx.')  
    bash.add_option('kymin', 'float', '', '', 'Minimum ky.')   
    bash.add_option('kymax', 'float', 'k', '', 'Maximum ky.') 
    bash.add_optionspace()
    
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
       
    # Get the bash arguments and execute the script
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_gamma_overview_spectra(**args)   

     
