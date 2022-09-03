 
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
from stellapy.plot.utils.species.display_multipleSpecies import get_speciesLineStyle
from stellapy.plot.utils.species.display_multipleSpecies import display_specieLineLabels

@noverbose 
def plot_velocityDistribution(\
            # Simulations
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
                species=[0],\
            # Data on the (x,y) axis    
                x_quantity="vpa",\
                y_quantity="g",\
                units="normalized",\
                x_range=None,\
                y_range=None,\
            # Modes
                kx_range=[0,0],\
                ky_range=[0,100],\
                modes_id="unstable",\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # Figure
                ax=None,\
                Progress=None,\
                show_figure=False,\
            # Toggles
                removelegend=False,\
                interpolate=None, 
                logy=False):  
    
    
    # Save the plotting details to a <plot> object
    plot = Plot()  
    plot.update_figure(Progress, ax, show_figure, logy=logy)
    plot.update_simulations(experiment_id, simulation_id, species=[0]) 
    plot.update_xyData(x_quantity, y_quantity, x_range, y_range, units)  
    plot.update_toggles(interpolate=interpolate)
    plot.update_modes(modes_id, kx_range, ky_range)
    plot.update_labels(title, x_label, y_label) 
    plot.update_legend(fontsize=16, removelegend=removelegend) 
    plot.process_plottingVariables(research)  
    
    # Update the progress bar of the GUI 
    displayProgressGUI = DisplayProgressGUI(research, "Plot "+plot.yname+" versus "+plot.xname, count="modes") 
     
    # Create the figure and a <legend> and <axis> object
    ax = create_figure(ax, plot, research) 
    legend = Legend(ax, plot); axis = Axis(ax, plot)  
    axis.set_limits(ybot=0, ytop=1)
    
    # Labels for multiple species
    display_specieLineLabels(ax, species) 
      
    # Define a color for each mode that is plotted
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,research.numberOfPlottedModes)); plot_i=0
    
    # Iterate over the experiments and simulations
    for experiment in research.plotted_experiments: 
        for simulation in experiment.plotted_simulations: 
            for mode in simulation.plotted_modes:  
                for specie in species:
                 
                    # Update the progress bar of the GUI  
                    displayProgressGUI.move_progressBar()
     
                    # Get the lines and marker styles for the experiments and simulations
                    style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation, mode) 
                    style['linestyle'] = get_speciesLineStyle(species, specie) if len(species)>1 else style['linestyle'] 
                    style['color'] = colors[plot_i]; style['linewidth'] = 2
                    if 'marker' not in style: style["marker"] = 'o'
                    
                    # Get the (time, quantity) data
                    vec_velocity, vec_quantity = get_dataForPlots(plot, mode, plot.quantity, specie)  
                    
                    # Only look at the shape 
                    vec_quantity = vec_quantity/np.max(vec_quantity)
                    
                    # Plot quantity(t)  
                    ax.plot(vec_velocity, vec_quantity, **style); plot_i += 1
                    
                    # Keep track of the axis limits 
                    axis.update_axisLimits(vec_velocity, vec_quantity)    
 
    # Add the legend and rescale the axis
    legend.add_legend()
    axis.rescale_axis()  
     
    # Show the figure
    if show_figure: plt.show()
 
    # Update the progress bar of the GUI  
    displayProgressGUI.finilize_progressBar()
    return ax 
 

################################################################################
#                  USE THIS PLOTTING FUNCTION AS A MAIN SCRIPT                 #
################################################################################ 
    
if __name__ == "__main__":   
    
    import pathlib, time
    from stellapy.simulations.Research import create_research  
    start = time.time()
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI/")  
    research = create_research(folders=folder); research.print_research()  
    plot_velocityDistribution(research, show_figure=False)
    print("    ---> It took "+str(time.time() - start)+" seconds to plot a quantity versus velocity.")
    plt.show()
    