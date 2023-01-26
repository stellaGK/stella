 
# Load the modules 
import sys
import numpy as np
import matplotlib.pyplot as plt

# Personal modules
from stellapy.utils.decorators.verbose import noverbose  

# Plotting modules 
from stellapy.GUI.plot.utils.Plot import Plot
from stellapy.GUI.plot.utils import Axis, Legend
from stellapy.GUI.plot.utils import create_figure
from stellapy.GUI.utils import DisplayProgressGUI  
from stellapy.plot.utils.style import get_styleForLinesAndMarkers 
from stellapy.utils.decorators.exit_program import exit_program

#===============================================================================
#                    PLOT QUANTITY(Z) FOR NONLINEAR SIMULATIONS
#===============================================================================

@noverbose 
def plot_quantityVsZ(\
            # Simulations
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
                specie=0,\
            # Data 
                x_quantity="pol",\
                y_quantity="phi2",\
                units="normalized",\
                x_range=None,\
                y_range=None,\
            # Eliminate time dimension
                percentage=50,\
                time="both",\
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
                rescaleToOne=True): 
     
    # Only implemented for specific data  
    if (y_quantity not in ["phi2", "phi2_zonal", "phi2_nozonal", "phi_real", "phi_imag", "g"]) or (x_quantity not in ["z", "zeta", "pol", "tor"]):
        exit_program("Not implemented", plot_quantityVsZ, sys._getframe().f_lineno)  

    # Save the plotting details to a <plot> object
    plot = Plot()  
    plot.update_xyData(x_quantity, y_quantity, x_range, y_range, units)  
    plot.update_simulations(experiment_id, simulation_id, species=[specie]) 
    plot.update_figure(Progress, ax, show_figure)
    plot.update_labels(title, x_label, y_label) 
    plot.update_toggles(rescaleToOne=rescaleToOne)
    plot.update_legend(fontsize=14,removelegend=removelegend) #loc='upper left'
    plot.process_plottingVariables(research)  
    
    # Update the progress bar of the GUI 
    displayProgressGUI = DisplayProgressGUI(research, "Plot "+plot.yname+" versus "+plot.xname)   
    
    # Create the figure and a <legend> and <axis> object
    ax = create_figure(ax, plot, research) 
    legend = Legend(ax, plot); axis = Axis(ax, plot, overshoot_y=1.1)      
    axis.set_limits(ybot_pos=0) 
          
    # Iterate over the experiments, simulations and the species 
    for experiment in research.plotted_experiments:  
        for simulation in experiment.plotted_simulations:   
            
            # Update the progress bar of the GUI  
            displayProgressGUI.move_progressBar()  
            
            # Get the (time,z,quanitity) data for the plots  
            if y_quantity=="phi_real": 
                vec_quantity = simulation.potential.phi_vs_tz.phi.real
                vec_time = simulation.potential.phi_vs_tz.t
                vec_z = simulation.potential.phi_vs_tz.z 
            if y_quantity=="phi_imag":
                vec_quantity = simulation.potential.phi_vs_tz.phi.imag
                vec_time = simulation.potential.phi_vs_tz.t
                vec_z = simulation.potential.phi_vs_tz.z 
            if y_quantity=="phi2":
                vec_quantity = simulation.potential.phi2_vs_tz.phi2
                vec_time = simulation.potential.phi2_vs_tz.t
                vec_z = simulation.potential.phi2_vs_tz.z 
            if y_quantity=="phi2_zonal":
                vec_quantity = simulation.potential.phi2_vs_tz_zonal.phi2_zonal
                vec_time = simulation.potential.phi2_vs_tz_zonal.t
                vec_z = simulation.potential.phi2_vs_tz_zonal.z 
            if y_quantity=="phi2_nozonal":
                vec_quantity = simulation.potential.phi2_vs_tz_nozonal.phi2_nozonal
                vec_time = simulation.potential.phi2_vs_tz_nozonal.t
                vec_z = simulation.potential.phi2_vs_tz_nozonal.z 
            if y_quantity=="g":
                vec_quantity = simulation.distribution.g_vs_tsz.g
                vec_quantity = vec_quantity.take(indices=specie, axis=1)
                vec_time = simulation.distribution.g_vs_tsz.t
                vec_z = simulation.distribution.g_vs_tsz.z  
            
            # Get the z-axis 
            if x_quantity == "z":     vec_z = simulation.vec.z 
            if x_quantity == "zeta":  vec_z = simulation.vec.zeta
            if x_quantity == "pol":   vec_z = simulation.vec.pol
            if x_quantity == "tor":   vec_z = simulation.vec.tor  
            
            # Get the lines and marker styles for the experiments and simulations
            style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation)
            if "linestyle" in style: del style["linestyle"]
            style["linewidth"] = 2
            
            # Average over time, or take a specific time point 
            if time=="averaged" or time=="both":
                t_index = int(len(vec_time)*percentage/100)
                vec_quantity_averaged = np.average(vec_quantity[t_index:,:], axis=0)
                if plot.rescaleToOne: vec_quantity_averaged = vec_quantity_averaged/np.max(np.abs(vec_quantity_averaged))  
                ax.plot(vec_z, vec_quantity_averaged, ls="-", **style) 
                axis.update_axisLimits(vec_z, vec_quantity_averaged)       
            if time=="specific" or time=="both":
                ls = (0, (1, 1)) if time=="both" else "-" 
                vec_quantity_specific1 = vec_quantity[-1,:]  
                if plot.rescaleToOne: vec_quantity_specific1 = vec_quantity_specific1/np.max(np.abs(vec_quantity_specific1))  
                ax.plot(vec_z, vec_quantity_specific1, ls=ls, **style) 
                axis.update_axisLimits(vec_z, vec_quantity_specific1)     

    # Plot the magnetic field
    bmag = simulation.geometry.bmag-np.min(simulation.geometry.bmag)
    bmag = bmag/np.max(bmag)*axis.ylims[1]
    ax.plot(vec_z, bmag, color="gray", linestyle="-", alpha=0.5)
                     
    # Add the legend and rescale the axis
    if not plot.removelegend: legend.add_legend()
    axis.rescale_axis()
    
    # Show the figure 
    if show_figure: plt.show()
     
    # Update the progress bar of the GUI  
    displayProgressGUI.finilize_progressBar()
    if True: return ax


