
# Load the modules 
import copy
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Personal modules
from stellapy.utils.decorators.verbose import noverbose 

# Plotting modules 
from stellapy.GUI.plot.utils.Plot import Plot
from stellapy.GUI.plot.utils import Axis, Legend
from stellapy.GUI.plot.utils import create_figure
from stellapy.GUI.utils import DisplayProgressGUI 
from stellapy.plot.utils.style import get_styleForLinesAndMarkers 
from stellapy.plot.utils.species.display_multipleSpecies import display_specieLineLabels, get_speciesLineStyle 

#===============================================================================
#              PLOT QUANTITY(KX,KY,MU,VPA) FOR NONLINEAR SIMULATIONS
#===============================================================================

@noverbose 
def plot_oneDimensionalSpectrum(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
                species=[0],\
            # Specify data range  
                x_quantity="kx",\
                y_quantity="qflux",\
                x_range=None,\
                y_range=None,\
                units="normalized",\
            # Eliminate time dimension
                percentage=50,\
                time="averaged",\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # For the GUI the figure object already exists 
                ax=None,\
                Progress=None,\
                show_figure=True,\
            # Toggles 
                fontsize=18,\
                removelegend=True,\
                cascade=False,\
                logx=False,\
                logy=True,\
                rescaleToOne=False):   
    
    # Save the plotting details to a <plot> object
    plot = Plot()  
    plot.update_xyData(x_quantity, y_quantity, x_range, y_range, units) #, logx=logx, logy=logy)  
    plot.update_toggles(rescaleToOne=rescaleToOne, cascade=cascade)
    plot.update_simulations(experiment_id, simulation_id, species=species) 
    plot.update_figure(Progress, ax, show_figure, logx=logx, logy=logy)
    plot.update_labels(title, x_label, y_label) 
    plot.update_legend(removelegend=removelegend, fontsize=fontsize) 
    plot.process_plottingVariables(research)  
    
    # Update the progress bar of the GUI 
    displayProgressGUI = DisplayProgressGUI(research, "Plot "+plot.yname+" versus "+plot.xname)  
     
    # Create the figure and a <legend> and <axis> object
    ax = create_figure(ax, plot, research) 
    legend = Legend(ax, plot)
    axis = Axis(ax, plot) 
    
    # Plot the electron and ion labels, and perhaps impurity labels, if we have multiple species
    display_specieLineLabels(ax, species) 
        
    # Iterate over the y_quantities, experiments, simulations and the species 
    for experiment in research.plotted_experiments:  
        for simulation in experiment.plotted_simulations:  
            for specie in species:
            
                # Update the progress bar of the GUI  
                displayProgressGUI.move_progressBar()

                # Get the data for the plots 
                if y_quantity=="phi2" and x_quantity=="kx":
                    vec_y = simulation.potential.phi2_vs_tkx.phi2
                    vec_x = simulation.potential.phi2_vs_tkx.kx 
                    vec_time = simulation.potential.phi2_vs_tkx.t
                if y_quantity=="phi2" and x_quantity=="ky":
                    vec_y = simulation.potential.phi2_vs_tky.phi2
                    vec_x = simulation.potential.phi2_vs_tky.ky 
                    vec_time = simulation.potential.phi2_vs_tky.t
                if y_quantity=="qflux" and x_quantity=="kx":
                    vec_y = simulation.fluxes.qflux_vs_tskx.qflux  
                    try: vec_y = vec_y.take(indices=specie, axis=1)
                    except: 
                        print("There is only one specie for: "+str(simulation.input_file))
                        vec_y = vec_y.take(indices=specie, axis=1)
                    vec_x = simulation.fluxes.qflux_vs_tskx.kx 
                    vec_time = simulation.fluxes.qflux_vs_tskx.t
                if y_quantity=="qflux" and x_quantity=="ky":
                    vec_y = simulation.fluxes.qflux_vs_tsky.qflux
                    vec_y = vec_y.take(indices=specie, axis=1)
                    vec_x = simulation.fluxes.qflux_vs_tsky.ky 
                    vec_time = simulation.fluxes.qflux_vs_tsky.t
                if y_quantity=="pflux" and x_quantity=="kx":
                    vec_y = simulation.fluxes.pflux_vs_tskx.pflux  
                    try: vec_y = vec_y.take(indices=specie, axis=1)
                    except: 
                        print("There is only one specie for: "+str(simulation.input_file))
                        vec_y = vec_y.take(indices=specie, axis=1)
                    vec_x = simulation.fluxes.pflux_vs_tskx.kx 
                    vec_time = simulation.fluxes.pflux_vs_tskx.t
                if y_quantity=="pflux" and x_quantity=="ky":
                    vec_y = simulation.fluxes.pflux_vs_tsky.pflux
                    vec_y = vec_y.take(indices=specie, axis=1)
                    vec_x = simulation.fluxes.pflux_vs_tsky.ky 
                    vec_time = simulation.fluxes.pflux_vs_tsky.t
                if y_quantity=="g" and x_quantity=="mu": 
                    vec_y = simulation.distribution.g_vs_tsmu.g
                    vec_y = vec_y.take(indices=specie, axis=1)
                    vec_x = simulation.distribution.g_vs_tsmu.mu 
                    vec_time = simulation.distribution.g_vs_tsmu.t 
                if y_quantity=="g" and x_quantity=="vpa": 
                    vec_y = simulation.distribution.g_vs_tsvpa.g
                    vec_y = vec_y.take(indices=specie, axis=1)
                    vec_x = simulation.distribution.g_vs_tsvpa.vpa
                    vec_time = simulation.distribution.g_vs_tsvpa.t
                    
                # Make sure we don't have ridiculous values
                #if logy: vec_y[vec_y<1.E-100] = np.nan
             
                # Get the lines and marker styles for the experiments and simulations
                style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation)
                style["linestyle"] = get_speciesLineStyle(species, specie);
                style["marker"] = 'o'; style["ms"] = 3; style["linewidth"] = 2
                 
                # Average over time, or take a specific time point 
                if time=="averaged" or time=="both":
                    t_index = int(len(vec_time)*percentage/100)
                    vec_y  = np.average(vec_y[t_index:,:], axis=0)  
                if time=="specific" or time=="both":
                    style["linestyle"] = (0, (1, 1)) if time=="both" else "-" 
                    vec_y = vec_y[-1,:]   
                    
                # Rescale 
                if plot.rescaleToOne: vec_y = vec_y/np.max(np.abs(vec_y))
                
                # Plot log 
                #if plot.takelogx: vec_x[vec_x<=0]=np.nan; vec_x = np.log10(vec_x)
                #if plot.takelogy: vec_y[vec_y<=0]=np.nan; vec_y = np.log10(vec_y) 
                if plot.takelogx: ax.set_xscale('log')
                if plot.takelogy: ax.set_yscale('log')
                
                # Plot 
                if not logy:
                    ax.plot(vec_x, vec_y, **style)  
                elif logy and len(research.plotted_experiments)==1 and len(experiment.plotted_simulations)==1: 
                    vec_y_pos = copy.deepcopy(vec_y); vec_y_pos[vec_y_pos<0] = np.nan
                    vec_y_neg = copy.deepcopy(vec_y); vec_y_neg[vec_y_neg>0] = np.nan
                    del style["marker"]; del style["color"]
                    ax.plot(vec_x, vec_y_pos, marker="o", color="navy", **style)  
                    ax.plot(vec_x, -vec_y_neg, marker="x", color="crimson", **style)  
                else: 
                    vec_y_pos = copy.deepcopy(vec_y); vec_y_pos[vec_y_pos<0] = np.nan
                    vec_y_neg = copy.deepcopy(vec_y); vec_y_neg[vec_y_neg>0] = np.nan
                    del style["marker"]
                    ax.plot(vec_x, vec_y_pos, marker="o", **style)  
                    ax.plot(vec_x, -vec_y_neg, marker="x", **style)  
               
                # Plot the cascade 7/3 power law
                if cascade and y_quantity not in ["gvmus", "gzvs"]:
                    # We need data_vs_ky and vec_ky
                    if x_quantity=="kx": vec_x = np.array(vec_x[1:])
                    if x_quantity=="kx": vec_y = np.real(vec_y[1:]) 
                    # Fit a linear curve to the (log,log) graph   
                    slope, intercept = linregress(np.log10(vec_x[6:]), np.log10(vec_y[6:]))[0:2]
                    # Therefore, the fitted y-data in the (log,log) graph is:
                    vec_y = 10**(np.log10(vec_x**(slope))+intercept)
                    # Ideally the exponent is 7/3, so keep the intercept and use the ideal exponent 
                    vec_y = 10**(np.log10(vec_x**(-7./3.))+intercept) 
                    # Plot the 7/3 power law
                    ax.plot(vec_x, vec_y, '--', color='black')

    # Add the legend and rescale the axis
    legend.add_legend()
    axis.rescale_axis()
    
    # Show the figure
    if show_figure: plt.show()

    # Update the progress bar of the GUI  
    displayProgressGUI.finilize_progressBar()
    if True: return ax  
        