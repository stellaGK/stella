 
# Load modules 
import numpy as np
import matplotlib.pyplot as plt

# Personal modules
from stellapy.utils.decorators.verbose import noverbose   
from stellapy.GUI.plot.utils.get_dataForPlots import get_dataForPlots  

# Plotting modules
from stellapy.GUI.plot.utils.Plot import Plot
from stellapy.GUI.plot.utils import Axis, Legend
from stellapy.GUI.plot.utils import create_figure
from stellapy.GUI.utils import DisplayProgressGUI  
from stellapy.plot.utils.style import get_styleForLinesAndMarkers 

#===============================================================================
#                    PLOT QUANTITY(T) FOR LINEAR SIMULATIONS                   #
#===============================================================================

@noverbose 
def plot_quantityVsTime(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
            # Data on the (x,y) axis  
                x_quantity="t",\
                y_quantity="omega",\
                x_range=None,\
                y_range=None,\
                units="normalized",\
            # Modes
                modes_id="unstable",\
                kx_range=[0,0],\
                ky_range=[0,100],\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # Figure
                ax=None,\
                Progress=None,\
                show_figure=False,\
            # Legend
                ncol=1,\
                fontsize=20,\
                handlelength=1,\
                removelegend=False): 

    # Save the plotting details to a <plot> object
    plot = Plot()  
    plot.update_figure(Progress, ax, show_figure)
    plot.update_simulations(experiment_id, simulation_id, species=[0]) 
    plot.update_xyData(x_quantity, y_quantity, x_range, y_range, units, logy=(y_quantity=="phi2")) 
    plot.update_legend(fontsize, handlelength, ncol, removelegend=removelegend)
    plot.update_modes(modes_id, kx_range, ky_range)
    plot.update_labels(title, x_label, y_label) 
    plot.process_plottingVariables(research)  
    
    # Update the progress bar of the GUI 
    displayProgressGUI = DisplayProgressGUI(research, "Plot "+plot.yname+" versus "+plot.xname, count="modes") 
    
    # Create the figure and a <legend> and <axis> object
    ax = create_figure(ax, plot, research) 
    legend = Legend(ax, plot)
    axis = Axis(ax, plot) 
    
    # Set some fixed axis limits depending on the plot
    if y_quantity=="phi2": axis.set_limits(xbot=0)
    if y_quantity!="phi2": axis.set_limits(percentage=0.9, ybot_pos=0, xbot=0)
    
    # Define a color for each mode that is plotted
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,research.numberOfPlottedModes)); plot_i=0
    
    # Iterate over the experiments and simulations
    for experiment in research.plotted_experiments: 
        for simulation in experiment.plotted_simulations: 
            for mode in simulation.plotted_modes:  
                
                # Update the progress bar of the GUI  
                displayProgressGUI.move_progressBar()
                    
                # Get the (time, quantity) data
                vec_time, vec_quantity = get_dataForPlots(plot, mode, plot.quantity)     
                
                # Get the lines and marker styles for the experiments and simulations
                style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation, mode)
                style['linestyle'] = '-'; style['color'] = colors[plot_i]; style['linewidth'] = 2

                # Plot quantity(t) 
                ax.plot(vec_time, vec_quantity, **style); plot_i += 1
                
                # Keep track of the axis limits 
                axis.update_axisLimits(vec_time, vec_quantity)        
                
            # If we have too many modes remove the legend 
            if len(simulation.plotted_modes)>40: plot.removelegend = True; print('removed legend!')

    # Add the legend and rescale the axis
    legend.add_legend()
    axis.rescale_axis()  
    
    # Show the figure  
    if show_figure and not Progress: plt.show()   
    return 
        
    
################################################################################
#                  USE THIS PLOTTING FUNCTION AS A MAIN SCRIPT                 #
################################################################################ 
    
if __name__ == "__main__":   
    
    import pathlib, time
    from stellapy.simulations.Research import create_research  
    start = time.time()
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI_PREVIOUS/fprim2tprim6/")    
    research = create_research(folders=folder)#, knob1="species_parameters_1", key1="tprim", knob2="species_parameters_1", key2="fprim")
    research.print_research()  
    plot_quantityVsTime(research, y_quantity="gamma", show_figure=False)
    print("    ---> It took "+str(time.time() - start)+" seconds to plot a quantity versus modes.")
    plt.show()    
    
    
