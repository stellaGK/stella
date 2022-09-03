        
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
#                    PLOT QUANTITY(Z) FOR LINEAR SIMULATIONS                   #
#===============================================================================
 
@noverbose 
def plot_quantityVsZ(\
            # Specify which simulations to plot
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
            # Data on the (x,y) axis  
                units="normalized",\
                x_quantity="z",\
                y_quantity="phi_real",\
                x_range=None,\
                y_range=None,\
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
                show_figure=True,\
            # Legend
                removelegend=False,\
                loc="upper left",\
                fontsize=16):
    ''' Plot the parallel mode structure of Re(phi), Imag(phi) or phi**2. 
    
    Note that the phi(Real) and phi(Imag) are normalized to max(abs(phi))!

    Parameters
    ----------
    units : {normalized; SI units}
    y_quantity : {phi2; phiReal; phiImag}
    x_quantity : {zeta; z; pol} 
    '''

    # Save the plotting details to a <plot> object
    plot = Plot()  
    plot.update_figure(Progress, ax, show_figure)
    plot.update_simulations(experiment_id, simulation_id, species=[0]) 
    plot.update_xyData(x_quantity, y_quantity, x_range, y_range, units) 
    plot.update_modes(modes_id, kx_range, ky_range)
    plot.update_labels(title, x_label, y_label) 
    plot.update_legend(fontsize=fontsize, loc=loc, removelegend=removelegend)
    plot.process_plottingVariables(research)  
    
    # Update the progress bar of the GUI 
    displayProgressGUI = DisplayProgressGUI(research, "Plot "+plot.yname+" versus "+plot.xname, count="modes") 

    # Create the figure and a <legend> and <axis> object
    ax = create_figure(ax, plot, research) 
    legend = Legend(ax, plot)
    axis = Axis(ax, plot) 

    # Set some fixed axis limits depending on the plot
    if y_quantity=="phi2": axis.set_limits(ybot=0, ytop=1)
    if y_quantity!="phi2": axis.set_limits(ybot=-1, ytop=1) 
    
    # Define a color for each mode that is plotted
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,research.numberOfPlottedModes)); plot_i=0

    # Iterate over the experiments and simulations
    for experiment in research.plotted_experiments: 
        for simulation in experiment.plotted_simulations: 
            for mode in simulation.plotted_modes:    
   
                # Update the progress bar of the GUI  
                displayProgressGUI.move_progressBar()
                   
                # Get the (z, quantity) data
                vec_z, vec_quantity = get_dataForPlots(plot, mode, plot.quantity)  
                
                # Look at the shape
                vec_quantity = vec_quantity/np.max(np.abs(vec_quantity))
                
                # Get the lines and marker styles for the experiments and simulations
                style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation, mode)
                style['linestyle'] = '-'; style['color'] = colors[plot_i]; style['linewidth'] = 2
                
                # Plot quantity(z)  
                ax.plot(vec_z, vec_quantity, **style); plot_i += 1
                   
                # Keep track of the axis limits  
                axis.update_axisLimits(vec_z, vec_quantity)    

    # Plot the magnetic field
    bmag = simulation.geometry.bmag
    bmag = bmag-np.min(bmag)+ax.get_ylim()[0]/ax.get_ylim()[1]
    bmag = bmag/np.max(bmag)*ax.get_ylim()[1]
    if plot_i!=0: ax.plot(vec_z, bmag, color="red", linestyle=":")
    
    # Add the legend and rescale the axis
    legend.add_legend()
    axis.rescale_axis()
    
    # Show the figure  
    if show_figure and not Progress: plt.show()   
    return
    
