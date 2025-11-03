"""

#===============================================================================
#               Plot flux_s(ky) and flux_s(kx) for each species s              #
#===============================================================================
 
Based on <folder> create a <research> and create a figure with 2*nspec subplots:
    - ax1*nspec: flux_s(kx)     (first row)
    - ax2*nspec: flux_s(ky)     (second row)
    
x_quantity
----------
Choose between {kx, ky}.
    >> plot_flux_spectra --ky
    >> plot_flux_vs_ky

y_quantity
----------
Choose between {qflux, pflux, vflux}.
    >> plot_flux_spectra --qflux (or --pflux or --vflux) 
    
folderIsExperiment
------------------
To ignore the automatic sorting of experiment and simulation, you can force stellapy
to create an <experiment> for each <subfolder>, the folder name will be the label
    >> plot_qflux_spectra -f (or --folderIsExperiment) 

Hanne Thienpondt 
15/12/2022

"""

#!/usr/bin/python3
import sys, os
import pathlib
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import get_styleForLinesAndMarkers
from stellapy.data.input.read_inputFile import read_nonlinearFullFluxSurfaceFromInputFile
from stellapy.plot.utils.labels.get_timeFrameString import get_timeFrameString
from stellapy.plot.utils.species.recognize_species import recognize_species
from stellapy.plot.utils.style.create_figure import update_figure_style
from stellapy.utils.files.get_firstInputFile import get_firstInputFile
from stellapy.plot.utils.labels.standardLabels import standardLabels
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.Research import create_research
from stellapy.plot.utils.style import Axis, Legend, Plot
from stellapy.utils.commandprompt.bash import Bash 

#===============================================================================
#               Plot flux_s(ky) and flux_s(kx) for each species s              #
#===============================================================================

def plot_flux_vs_wavenumber(
        folder,
        # Quantities to be plotted
        x_quantity="both",
        y_quantity="qflux",
        specie=None, 
        # Plotting options
        log=False,
        folderIsExperiment=False,
        normalize_to_one=False):
    
    # Create a <research> based on the given <folder>
    try: research = create_research(folders=folder, folderIsExperiment=folderIsExperiment)
    except: research = create_research(folders=folder, folderIsExperiment=True)
    simulation = research.experiments[0].simulations[0]
    
    # Get the species and the x-quantity
    nspec = simulation.input.nspec if specie==None else 1
    species = range(simulation.input.nspec) if specie==None else [specie]
    nk = 2 if x_quantity=="both" else 1
    x1_quantity = "kx" if x_quantity=="both" else x_quantity
    
    # Create a figure
    if y_quantity=="qflux": title = "Heat flux spectra"
    if y_quantity=="pflux": title = "Particle flux spectra"
    if y_quantity=="vflux": title = "Momentum flux spectra"
    fig = plt.figure(figsize=(18, 9)); axes = []
    grid_specifications = gridspec.GridSpec(nk, nspec)
    if nspec==1: grid_specifications.update(top=0.93, left=0.08, right=0.97, bottom=0.1, wspace=0.2, hspace=0.3)
    if nspec==2: grid_specifications.update(top=0.93, left=0.08, right=0.97, bottom=0.1, wspace=0.2, hspace=0.3)
    if nspec==3: grid_specifications.update(top=0.92, left=0.08, right=0.97, bottom=0.1, wspace=0.3, hspace=0.3)
    for i in range(nspec*nk): axes.append(plt.subplot(grid_specifications[i]))
    update_figure_style(fig, axes)
    fig.suptitle(title)
    
    # Add the data to the plot
    for i, specie in enumerate(species):
        subplot_flux_vs_wavenumber(axes[i], research, x1_quantity, y_quantity, specie=specie, log=log, normalize_to_one=normalize_to_one)
        if nk==2: subplot_flux_vs_wavenumber(axes[i+nspec], research, "ky", y_quantity, specie=specie, log=log, normalize_to_one=normalize_to_one)
    
    # Appearance  
    for ax in axes: (ax.ticklabel_format(style='plain', useOffset=False) if (axes[0].get_ylim()[1]<1000) else ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))) if not log else ""
    for ax in axes: ax.xaxis.labelpad = 10
    for ax in axes: ax.yaxis.labelpad = 15
    
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    plt.show()
    return

#===============================================================================
#                  Plot flux_s(ky) or flux_s(kx) for a specie s                #
#===============================================================================

def subplot_flux_vs_wavenumber(
        ax, research, 
        # Quantities to be plotted
        x_quantity="ky", 
        y_quantity="qflux", 
        specie=0, 
        # Plotting options
        log=False, 
        normalize_to_one=False):
    
    # Automate the axis limits and legend
    plot = Plot(); legend = Legend(ax, plot) 
    axis = Axis(ax, plot, xbot_pos=0, ytop_neg=0, overshoot_y=1.1, logy=log) 
    plot.process_plottingVariables(research)
    tstart = [np.nan,np.nan]; tend = [np.nan,np.nan] 
    
    # Check whether <specie> is a valid choice 
    if specie >= research.experiments[0].simulations[0].dim.species:
        dim_species = research.experiments[0].simulations[0].dim.species
        exit_reason = "Error: specie = "+str(specie)+" was selected but only "+str(dim_species)+" species are present in the simulation.\n"
        exit_reason += "Please choose specie from {"+", ".join([str(s) for s in range(dim_species)])+"} instead.\n"
        exit_program(exit_reason, subplot_flux_vs_wavenumber, sys._getframe().f_lineno)
    
    # Iterate over the experiments and simulations
    for experiment in research.experiments:
        for simulation in experiment.simulations:
            
            # Get the lines and marker styles for the experiments and simulations
            style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation)
            
            # Get flux and kx or ky
            x, y, tstart, tend = get_quantity_data(simulation, x_quantity, y_quantity, specie, tstart, tend)
            if normalize_to_one: y = y/np.max(np.abs(y))
            if np.any(y<0): log = False
            
            # Plot qflux(kx) or qflux(ky) or pflux(kx) or vflux(kx) or ...
            ax.plot(x, y, marker="o", **style)
        
            # Keep track of the axis limits
            axis.update_axisLimits(x, y)
    
    # Axis labels
    specie_label = recognize_species(research, specie) 
    extra = "$\\sum_{k_x}$" if x_quantity=="ky" else "$\\sum_{k_y}$"
    ylabel = standardLabels["normalized"][y_quantity].replace("{s}",specie_label)
    ylabel = extra+"$\\langle$" + ylabel + "$\\rangle_{z, t="+get_timeFrameString(tstart, tend)+"}$"
    ax.set_xlabel(standardLabels["normalized"][x_quantity]) 
    ax.set_ylabel(ylabel) 

    # Automatically set the axis limits and legends 
    if normalize_to_one: ax.set_ylim([0,1])
    if log: ax.set_yscale("log")
    legend.add_legend()
    axis.rescale_axis()  
    return 

#===============================================================================
#                  Read flux_s(ky) or flux_s(kx) for a specie s                #
#===============================================================================

def get_quantity_data(simulation, x_quantity, y_quantity, specie, tstart, tend):
    
    # Data for the x-axis
    if x_quantity=="kx": x = simulation.vec.kx
    if x_quantity=="ky": x = simulation.vec.ky
    
    # Get the data for the y-axis
    try:
        try: y, t0, t1 = get_quantity_data_from3DFiles(simulation, x_quantity, y_quantity, specie)
        except: y, t0, t1 = get_quantity_data_fromSaturatedFiles(simulation, x_quantity, y_quantity, specie)
    except: 
        print('WARNING: The flux(kx,ky) data could not be found.')
        return np.array([np.nan]), np.array([np.nan]), 500, 1000
    
    # Keep track of the time frames
    tstart[0] = np.nanmin([tstart[0], t0])
    tstart[1] = np.nanmax([tstart[1], t0])
    tend[0] = np.nanmin([tend[0], t1])
    tend[1] = np.nanmax([tend[1], t1])
    return x, y, tstart, tend 

#-------------------------------
def get_quantity_data_from3DFiles(simulation, x_quantity, y_quantity, specie): 
    if x_quantity=="kx": 
        if y_quantity=="qflux": y = simulation.fluxes.qflux_vs_tskx.qflux[:,specie,:]; t = simulation.fluxes.qflux_vs_tskx.t
        if y_quantity=="pflux": y = simulation.fluxes.pflux_vs_tskx.pflux[:,specie,:]; t = simulation.fluxes.pflux_vs_tskx.t
        if y_quantity=="vflux": y = simulation.fluxes.vflux_vs_tskx.vflux[:,specie,:]; t = simulation.fluxes.vflux_vs_tskx.t
    if x_quantity=="ky": 
        if y_quantity=="qflux": y = simulation.fluxes.qflux_vs_tsky.qflux[:,specie,:]; t = simulation.fluxes.qflux_vs_tsky.t
        if y_quantity=="pflux": y = simulation.fluxes.pflux_vs_tsky.pflux[:,specie,:]; t = simulation.fluxes.pflux_vs_tsky.t
        if y_quantity=="vflux": y = simulation.fluxes.vflux_vs_tsky.vflux[:,specie,:]; t = simulation.fluxes.vflux_vs_tsky.t
    y = np.nanmean(y[t>simulation.time.tstart,:],axis=0)
    return y, simulation.time.tstart, t[-1]
        
#-------------------------------
def get_quantity_data_fromSaturatedFiles(simulation, x_quantity, y_quantity, specie):
    trange = simulation.saturated.trange
    if y_quantity=="qflux": y = simulation.saturated.qflux_vs_szkxky.qflux[specie,:,:,:]
    if y_quantity=="pflux": y = simulation.saturated.pflux_vs_szkxky.pflux[specie,:,:,:]
    if y_quantity=="vflux": y = simulation.saturated.vflux_vs_szkxky.vflux[specie,:,:,:]
    if x_quantity=="kx": y = y[:,0] + 2*np.sum(y[:,1:],axis=1) 
    if x_quantity=="ky": y = np.sum(y[:,:],axis=0)
    return y, trange[0], trange[1]

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Check which script to launch      
    input_file = get_firstInputFile(pathlib.Path(os.getcwd())) 
    nonlinear, full_flux_surface = read_nonlinearFullFluxSurfaceFromInputFile(input_file)
    if not nonlinear: os.system("python3 $STELLAPY/plot/linear/qflux_vs_wavenumber.py "+" ".join(sys.argv[1:])); sys.exit()
    
    # Create a bash interface 
    bash = Bash(plot_flux_vs_wavenumber, __doc__)     
    
    # Data toggles for the x-quantities
    bash.add_toggleheader("x-quantity")
    bash.add_toggle('x_quantity', 'kx', '', '', 'Plot the kx-spectra.') 
    bash.add_toggle('x_quantity', 'ky', '', '', 'Plot the ky-spectra.') 
    bash.add_toggle('x_quantity', 'both', '', '', 'Plot the kx- and ky-spectra.') 
    bash.add_togglespace()
    
    # Data toggles for the y-quantities
    bash.add_toggleheader("y-quantity")
    bash.add_toggle('y_quantity', 'qflux', '', '', 'Plot the kx- and ky-spectra of the heat flux.') 
    bash.add_toggle('y_quantity', 'pflux', '', '', 'Plot the kx- and ky-spectra of the particle flux.') 
    bash.add_toggle('y_quantity', 'vflux', '', '', 'Plot the kx- and ky-spectra of the momentum flux.') 
    bash.add_togglespace()
    
    # Species
    bash.add_toggleheader("specie")
    bash.add_option('specie', 'int', 's', '', 'Select the species as "-s 0" for ions and "-s 1" for electrons.') 
    bash.add_toggle('specie', 0, '', 'ions', 'Plot the fluxes for the ions.') 
    bash.add_toggle('specie', 1, '', 'electrons', 'Plot the fluxes for the electrons (assuming s=1 for electrons).') 
    bash.add_toggle('specie', 2, '', 'impurity', 'Plot the fluxes for the impurities (assuming s=2 for impurities).')
    bash.add_togglespace() 
    
    # Research options      
    bash.add_toggleheader("other")
    bash.add_toggle('folderIsExperiment', True, 'f', 'folderIsExperiment', 'Each folder is an experiment.')    
    
    # Plot options  
    bash.add_toggle('normalize_to_one', True, 'n', 'normalize_to_one', 'Rescale the data to have the maximum at one.')  
    bash.add_toggle('log', True, 'l', 'log', 'Use logaritmic scales.')  
    
    # Get the arguments and extract <x_quantity>
    args = bash.get_arguments() 
    
    # Launch the script
    plot_flux_vs_wavenumber(**args)    


