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
Choose between {qflux, pflux, vflux, phi2}.
    >> plot_flux_spectra --qflux (or --pflux or --vflux or --phi2) 
    
folderIsExperiment
------------------
To ignore the automatic sorting of experiment and simulation, you can force stellapy
to create an <experiment> for each <subfolder>, the folder name will be the label
    >> plot_flux_vs_ky -f (or --folderIsExperiment) 

Hanne Thienpondt 
07/09/2022

"""

#!/usr/bin/python3  
import sys, os
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec

# Personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])   
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import get_styleForLinesAndMarkers
from stellapy.plot.utils.labels.get_timeFrameString import get_timeFrameString
from stellapy.plot.utils.species.recognize_species import recognize_species
from stellapy.plot.utils.style.create_figure import update_figure_style
from stellapy.plot.utils.labels.standardLabels import standardLabels  
from stellapy.simulations.Research import create_research   
from stellapy.GUI.plot.utils import Axis, Legend, Plot 
from stellapy.utils.commandprompt.bash import Bash 

#===============================================================================
#               Plot flux_s(ky) and flux_s(kx) for each species s              #
#===============================================================================

def plot_flux_vs_wavenumber_both(folder, y_quantity="qflux", log=False, folderIsExperiment=False, normalize_to_one=False): 
    
    # Create a <research> based on the given <folder>
    try: research = create_research(folders=folder, folderIsExperiment=folderIsExperiment)
    except: research = create_research(folders=folder, folderIsExperiment=True)
    simulation = research.experiments[0].simulations[0]
    
    # Create a figure
    if y_quantity=="qflux": title = "Heat flux spectra"
    if y_quantity=="pflux": title = "Particle flux spectra"
    if y_quantity=="vflux": title = "Momentum flux spectra"
    if y_quantity=="phi2":  title = "Potential squared spectra"
    fig = plt.figure(figsize=(18, 9)); axes = []
    grid_specifications = gridspec.GridSpec(2, simulation.input.nspec)
    if simulation.input.nspec==2: grid_specifications.update(top=0.93, left=0.08, right=0.97, bottom=0.1, wspace=0.2, hspace=0.3)
    if simulation.input.nspec==3: grid_specifications.update(top=0.92, left=0.08, right=0.97, bottom=0.1, wspace=0.3, hspace=0.3)
    for i in range(simulation.input.nspec*2): axes.append(plt.subplot(grid_specifications[i]))
    update_figure_style(fig, axes)  
    fig.suptitle(title)
    
    # Add the data to the plot
    for i in range(simulation.input.nspec):
        subplot_flux_vs_wavenumber(axes[i], research, "kx", y_quantity, specie=i, log=log, normalize_to_one=normalize_to_one)
        subplot_flux_vs_wavenumber(axes[i+simulation.input.nspec], research, "ky", y_quantity, specie=i, log=log, normalize_to_one=normalize_to_one)
    
    # Appearance  
    for ax in axes: ax.ticklabel_format(style='plain', useOffset=False) if (axes[0].get_ylim()[1]<1000) else ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    for ax in axes: ax.xaxis.labelpad = 10
    for ax in axes: ax.yaxis.labelpad = 15
    
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    plt.show()
    return

#===============================================================================
#                Plot flux_s(ky) or flux_s(kx) for each species s              #
#===============================================================================

def plot_flux_vs_wavenumber(folder, x_quantity="ky", y_quantity="qflux", log=False, folderIsExperiment=False, normalize_to_one=False): 
    
    # Create a <research> based on the given <folder>
    try: research = create_research(folders=folder, folderIsExperiment=folderIsExperiment)
    except: research = create_research(folders=folder, folderIsExperiment=True)
    simulation = research.experiments[0].simulations[0]
    
    # Create a figure
    if y_quantity=="qflux": title = "Heat flux spectra"
    if y_quantity=="pflux": title = "Particle flux spectra"
    if y_quantity=="vflux": title = "Momentum flux spectra"
    if y_quantity=="phi2":  title = "Potential squared spectra"
    fig = plt.figure(figsize=(18, 9)); axes = []
    grid_specifications = gridspec.GridSpec(1, simulation.input.nspec)
    if simulation.input.nspec==2: grid_specifications.update(top=0.93, left=0.08, right=0.97, bottom=0.1, wspace=0.2, hspace=0.3)
    if simulation.input.nspec==3: grid_specifications.update(top=0.92, left=0.08, right=0.97, bottom=0.1, wspace=0.3, hspace=0.3)
    for i in range(simulation.input.nspec): axes.append(plt.subplot(grid_specifications[i]))
    update_figure_style(fig, axes)  
    fig.suptitle(title)
    
    # Add the data to the plot
    for i in range(simulation.input.nspec):
        subplot_flux_vs_wavenumber(axes[i], research, x_quantity, y_quantity, specie=i, log=log, normalize_to_one=normalize_to_one)
    
    # Appearance
    for ax in axes: ax.ticklabel_format(style='plain', useOffset=False) 
    for ax in axes: ax.xaxis.labelpad = 10
    for ax in axes: ax.yaxis.labelpad = 15
    
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    plt.show()
    return 

#===============================================================================
#                  Plot flux_s(ky) or flux_s(kx) for a specie s                #
#===============================================================================

def subplot_flux_vs_wavenumber(ax, research, x_quantity="ky", y_quantity="qflux", specie=0, log=True, normalize_to_one=False):
    
    # Automate the axis limits and legend
    plot = Plot(); legend = Legend(ax, plot) 
    axis = Axis(ax, plot, xbot_pos=0, ytop_neg=0, overshoot_y=1.1) 
    plot.process_plottingVariables(research)
    tstart = [np.nan,np.nan]; tend = [np.nan,np.nan] 
        
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
    try: y, t0, t1 = get_quantity_data_from3DFiles(simulation, x_quantity, y_quantity, specie)
    except: y, t0, t1 = get_quantity_data_fromSaturatedFiles(simulation, x_quantity, y_quantity, specie)
    
    # Keep track of the time frames
    tstart[0] = np.nanmin([tstart[0], t0]) 
    tstart[1] = np.nanmax([tstart[1], t0]) 
    tend[0] = np.nanmin([tend[0], t1]) 
    tend[1] = np.nanmax([tend[1], t1]) 
    return x, y, tstart, tend 

#-------------------------------
def get_quantity_data_from3DFiles(simulation, x_quantity, y_quantity, specie): 
    if x_quantity=="kx": 
        if y_quantity=="phi2":  y = simulation.potential.phi2_vs_tkx.phi2[:,:]; t = simulation.potential.phi2_vs_tkx.t
        if y_quantity=="qflux": y = simulation.fluxes.qflux_vs_tskx.qflux[:,specie,:]; t = simulation.fluxes.qflux_vs_tskx.t
        if y_quantity=="pflux": y = simulation.fluxes.pflux_vs_tskx.pflux[:,specie,:]; t = simulation.fluxes.pflux_vs_tskx.t
        if y_quantity=="vflux": y = simulation.fluxes.vflux_vs_tskx.vflux[:,specie,:]; t = simulation.fluxes.vflux_vs_tskx.t
    if x_quantity=="ky": 
        if y_quantity=="phi2":  y = simulation.potential.phi2_vs_tky.phi2[:,:]; t = simulation.potential.phi2_vs_tky.t
        if y_quantity=="qflux": y = simulation.fluxes.qflux_vs_tsky.qflux[:,specie,:]; t = simulation.fluxes.qflux_vs_tsky.t
        if y_quantity=="pflux": y = simulation.fluxes.pflux_vs_tsky.pflux[:,specie,:]; t = simulation.fluxes.pflux_vs_tsky.t
        if y_quantity=="vflux": y = simulation.fluxes.vflux_vs_tsky.vflux[:,specie,:]; t = simulation.fluxes.vflux_vs_tsky.t
    y = np.mean(y[t>simulation.time.tstart,:],axis=0) 
    return y, simulation.time.tstart, t[-1]
        
#-------------------------------
def get_quantity_data_fromSaturatedFiles(simulation, x_quantity, y_quantity, specie):
    trange = simulation.saturated.trange
    if y_quantity=="phi2":  y = simulation.saturated.phi2_vs_kxky.phi2[:,:]
    if y_quantity=="qflux": y = simulation.saturated.qflux_vs_szkxky.qflux[specie,:,:,:]
    if y_quantity=="pflux": y = simulation.saturated.pflux_vs_szkxky.pflux[specie,:,:,:]
    if y_quantity=="vflux": y = simulation.saturated.vflux_vs_szkxky.vflux[specie,:,:,:]
    if y_quantity!="phi2":  y = np.sum(y[:,:,:]*simulation.geometry.dl_over_B[:,np.newaxis,np.newaxis], axis=0)
    if x_quantity=="kx": y = y[:,0] + 2*np.sum(y[:,1:],axis=1) 
    if x_quantity=="ky": y = np.sum(y[:,:],axis=0)
    return y, trange[0], trange[1]

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash interface
    bash = Bash(plot_flux_vs_wavenumber_both, __doc__)    
    
    # Data options
    bash.add_option('specie', 'int', 's', 'Select the species as "-s 0" for ions and "-s 1" for electrons.') 
    
    # Data toggles for the x- and y-quantities
    bash.add_toggle('x_quantity', 'kx', 'Plot the kx-spectra.') 
    bash.add_toggle('x_quantity', 'ky', 'Plot the ky-spectra.') 
    bash.add_toggle('y_quantity', 'qflux', 'Plot the kx- and ky-spectra of the heat flux.') 
    bash.add_toggle('y_quantity', 'pflux', 'Plot the kx- and ky-spectra  of the particle flux.') 
    bash.add_toggle('y_quantity', 'vflux', 'Plot the kx- and ky-spectra  of the momentum flux.') 
    bash.add_toggle('y_quantity', 'phi2', 'Plot the kx- and ky-spectra  of the potential squared.') 
    
    # Research options      
    bash.add_option('folderIsExperiment', 'True', 'f', 'Each folder is an experiment.')    
    
    # Plot options 
    bash.add_option('log', 'True', 'l', 'Use logaritmic scales.') 
    bash.add_option('normalize_to_one', 'True', 'n', 'Rescale the data to have the maximum at one.')  
    
    # Get the arguments and extract <x_quantity>
    args = bash.get_arguments()
    x_quantity = args['x_quantity'] if 'x_quantity' in args else None 
    
    # Launch the script
    if not x_quantity: plot_flux_vs_wavenumber_both(**args)   
    if x_quantity=="kx": plot_flux_vs_wavenumber(**args)   
    if x_quantity=="ky": plot_flux_vs_wavenumber(**args)   


