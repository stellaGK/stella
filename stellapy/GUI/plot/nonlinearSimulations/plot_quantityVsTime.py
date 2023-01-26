 
# Load the modules 
import numpy as np
import matplotlib.pyplot as plt

# Personal modules
from stellapy.utils.decorators.verbose import noverbose 
from stellapy.GUI.plot.utils.get_dataForPlots import get_yQuantities
from stellapy.GUI.plot.utils.get_dataForPlots import get_dataForPlots 

# Plotting modules 
from stellapy.GUI.utils import DisplayProgressGUI  
from stellapy.GUI.plot.utils.Plot import Plot
from stellapy.GUI.plot.utils import Axis, Legend
from stellapy.GUI.plot.utils import create_figure
from stellapy.plot.utils.style import get_styleForLinesAndMarkers 
from stellapy.plot.utils.species.display_multipleSpecies import display_specieLineLabels, get_speciesLineStyle 

#===============================================================================
#                    PLOT QUANTITY(T) FOR NONLINEAR SIMULATIONS
#===============================================================================

@noverbose 
def plot_quantityVsTime(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
                species=[0],\
            # Specify data to plot
                x_quantity="t",\
                y_quantity="qflux",\
                x_range=None,\
                y_range=None,\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # Figure
                ax=None,\
                Progress=None,\
                show_figure=False,\
            # Legend 
                shadow=True,\
                legend_loc=None,\
                fontsize=18,\
                removelegend=False,\
            # Toggles 
                print_unconverged_sims=None,\
                linewidth=None,\
                tick_style="sci",\
                show_timeFrame=False,\
                normalize=False,\
                units="normalized"): 
    
    # If we look at the heat flux, print the time frames on the prompt
    if y_quantity=="qflux": print_header("start")

    # We can have combined graphs (phiReal, phiImag) (Phi2, Phi2_nozonal, Phi2_zonal) 
    # and (Phi2_nozonal, qflux) for which we will iterate over certain y-quantities
    y_quantities = get_yQuantities(y_quantity)  
    
    # Save the plotting details to a <plot> object
    plot = Plot()  
    plot.update_xyData(x_quantity, y_quantity, x_range, y_range, units, y_quantities=y_quantities)  
    plot.update_toggles(show_timeFrame=show_timeFrame, normalize=normalize)
    plot.update_simulations(experiment_id, simulation_id, species=species) 
    plot.update_figure(Progress, ax, show_figure)
    plot.update_labels(title, x_label, y_label) 
    plot.update_legend(fontsize=fontsize, removelegend=removelegend, loc=legend_loc, shadow=shadow) 
    plot.update_axis(tick_style)
    plot.process_plottingVariables(research)  
    
    # Update the progress bar of the GUI 
    displayProgressGUI = DisplayProgressGUI(research, "Plot "+plot.yname+" versus "+plot.xname)  
      
    # Create the figure and a <legend> and <axis> object 
    ax = create_figure(ax, plot, research) 
    legend = Legend(ax, plot); axis = Axis(ax, plot, overshoot_y=1.2) 
    axis.set_limits(xbot_pos=0, ytop_neg=0,  ybot_pos=0, percentage=0.9)  
         
    # Plot the electron and ion labels, and perhaps impurity labels, if we have multiple species
    display_specieLineLabels(ax, species)   
    
    # Iterate over the y_quantities, experiments, simulations and the species
    for y_quantity in y_quantities:     
        for experiment in research.plotted_experiments:  
            for simulation in experiment.plotted_simulations: 
                for specie in species:

                    # Update the progress bar of the GUI  
                    displayProgressGUI.move_progressBar() 
                     
                    # Get the (time,quanitity) data for the plots 
                    vec_time, vec_quantity = get_dataForPlots(plot, simulation, plot.quantity, specie)  
                    if specie=="both": vec_quantity = np.sum(vec_quantity, axis=2)
                     
                    # Get the lines and marker styles for the experiments and simulations
                    style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation)
                    style["linestyle"] = get_speciesLineStyle(species, specie)
                    if linewidth: style["linewidth"] = linewidth 
                                         
                    # Plot the data  
                    ax.plot(vec_time, vec_quantity, **style)   
                             
                    # Add a line to indicate the saturated flux mean     
                    satflux = plot_saturatedDataOverTimeFrame(ax, simulation, vec_time, vec_quantity, y_quantity, units, specie, style, print_unconverged_sims)
                     
                    # Show the time frame for which the saturated flux is calculated  
                    if show_timeFrame: 
                        tstart, tend = simulation.time.tstart, simulation.time.tend 
                        ax.axvspan(tstart, tend, alpha=0.3, color=style["color"])  
                        
                    # Keep track of the axis limits 
                    axis.update_axisLimits(vec_time, vec_quantity)  
                    
                    # For the heat flux, print some information
                    if y_quantity=="qflux": print_information(style["label"], vec_time, simulation.time.tstart, simulation.time.tend, satflux)      
 
    # Add the legend and rescale the axis
    legend.add_legend()
    axis.rescale_axis()
 
    # Show the figure 
    if show_figure: plt.show()
     
    # Update the progress bar of the GUI  
    displayProgressGUI.finilize_progressBar()
    if y_quantity=="qflux": print_header("end")
    return axis
    
#===============================================================================
#                                METHODS
#===============================================================================

def plot_saturatedDataOverTimeFrame(ax, simulation, vec_time, vec_quantity, y_quantity, units, specie, style, print_unconverged_sims):
    
    # Get the time frame
    tstart, tend = simulation.time.tstart, simulation.time.tend 
    tstart  = tstart*simulation.referenceunits.ref['t'] if (units=="SI") else tstart 
    tend    = tend*simulation.referenceunits.ref['t'] if (units=="SI") else tend 
    
    # Plot the saturated flux
    if "flux" in y_quantity:
        satflux = simulation.time.saturatedFluxes[y_quantity][specie]
        satflux = satflux*simulation.referenceunits.ref['qflux'] if (units=="SI" and y_quantity=="qflux") else satflux 
        satflux = satflux*simulation.referenceunits.ref['pflux'] if (units=="SI" and y_quantity=="pflux") else satflux 
        satflux = satflux*simulation.referenceunits.ref['vflux'] if (units=="SI" and y_quantity=="vflux") else satflux 
        ax.plot([tstart, tend], satflux*np.ones(2),  **style)   
        if print_unconverged_sims and (tend<495 or (tend<995 and tend>510)): 
            if print_unconverged_sims==True: print_unconverged_sims = 1
            ax.text(0.6, 0.95-0.05*print_unconverged_sims, style["label"]+": tend = "+str(int(tend)), horizontalalignment='center', \
                 verticalalignment='center', transform=ax.transAxes, fontsize=12, color="red")  
            print_unconverged_sims += 1
        return satflux
            
    # Plot the saturated quantity
    else:
        saturated_quantity = np.mean(vec_quantity[(vec_time>=tstart) & (vec_time<=tend)])
        ax.plot([tstart, tend], saturated_quantity*np.ones(2),  **style)   
        return saturated_quantity
                            
#-------------------------
def print_header(id_):
    if id_=="start":
        print("")
        print("     ","".center(110,"=")) 
        print("     "," HEAT FLUX ".center(110,"="))
        print("     ","".center(110,"="))
        print()
        print("    " 
              "{0:^40}".format("SIMULATION") + \
              "{0:^15}".format("LAST TIME") + \
              "{0:^15}".format("T START") + \
              "{0:^15}".format("T END") + \
              "{0:^15}".format("QSAT"))
        print("    " 
              "{0:^40}".format("----------") + \
              "{0:^15}".format("---------") + \
              "{0:^15}".format("-------") + \
              "{0:^15}".format("-----") + \
              "{0:^15}".format("---"))
    if id_=="end":
        print()
        print("     ","".center(110,"=")) 
        print()
        
def print_information(id_, time, tstart, tend, satflux):
    print("    " 
          "{0:^40}".format(id_) + \
          "{0:^15}".format(time[-1]) + \
          "{0:^15}".format(tstart) + \
          "{0:^15}".format(tend) + \
          "{0:^15}".format(satflux))
    
################################################################################
#                  USE THIS PLOTTING FUNCTION AS A MAIN SCRIPT                 #
################################################################################ 
     
if __name__ == "__main__":   
     
    import pathlib, time
    from stellapy.simulations.Research import create_research  
    start = time.time()
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/fprim2tprim3/")    
    research = create_research(folders=folder) 
    research.print_research()  
    plot_quantityVsTime(research, show_figure=False)
    print("    ---> It took "+str(time.time() - start)+" seconds to plot a quantity versus modes.")
    plt.show()    
     
     