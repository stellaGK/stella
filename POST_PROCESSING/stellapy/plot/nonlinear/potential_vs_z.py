"""

#===============================================================================
#                                 Plot phi2(z)                                 #
#===============================================================================

Based on <folder> create a <research> and plot phi2(z).

Arguments
--------- 
    folder : pathlib.Path directory where the command has been executed
    x_quantity : {z, pol, tor}
    y_quantity : {phi2, phi, phi_real, phi_imag}
    geometry : {bmag, gradpar, gds2, gds21, gds22, gds23, gds24, cvdrift, gbdrift0, bmag_psi0} overlayed on the background 
    folderIsExperiment : {False, True} 
    interpolate : {int, False} where <int> (e.g. 20) is the interpolation step

Hanne Thienpondt 
16/12/2022

"""

#!/usr/bin/python3   
import pathlib
import sys, os
import numpy as np
import matplotlib as mpl  
import matplotlib.pyplot as plt 
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import get_styleForLinesAndMarkers 
from stellapy.data.input.read_inputFile import read_nonlinearFullFluxSurfaceFromInputFile
from stellapy.plot.utils.labels.get_timeFrameString import get_timeFrameString
from stellapy.utils.files.get_firstInputFile import get_firstInputFile
from stellapy.plot.utils.labels.standardLabels import standardLabels  
from stellapy.plot.utils.style.create_figure import create_figure 
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.Research import create_research 
from stellapy.plot.utils.style import Axis, Legend, Plot
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                                 Plot phi2(z)                                 #
#===============================================================================

def plot_potential_vs_z(
        folder, 
        # Quantities to be plotted
        x_quantity="z", 
        y_quantity="phi", 
        specie=0, 
        # Overlay a geometric quantity in gray
        geometry="bmag", 
        # Plotting options
        log=False,
        color=None, 
        interpolate=20, 
        normalize_to_one=False,
        folderIsExperiment=False):        
    
    # Create a <research> based on the given <folder>
    try: research = create_research(folders=folder, folderIsExperiment=folderIsExperiment)
    except: research = create_research(folders=folder, folderIsExperiment=True)
            
    # Create a figure
    if y_quantity=="phi": title = "Parallel mode structure of $|\\varphi|$"
    elif y_quantity=="phi_real": title = "Parallel mode structure of real($|\\varphi|$)"
    elif y_quantity=="phi_imag": title = "Parallel mode structure of imag($|\\varphi|$)"
    elif y_quantity=="phi2": title = "Parallel mode structure of the $|\\varphi|^2$"
    elif y_quantity=="phi2_zonal": title = "Parallel mode structure of the zonal modes of $|\\varphi|^2$"
    elif y_quantity=="phi2_nozonal": title = "Parallel mode structure of $|\\varphi|^2$ without zonal modes"
    elif y_quantity=="g2": title = "Parallel mode structure of the distribution $|g|^2$"
    else: exit_program("<y_quantity> = "+y_quantity+" is not a valid choice.", plot_potential_vs_z, sys._getframe().f_lineno) 
    ax = create_figure(title=title)  
    
    # Plot
    subplot_phi_vs_z(ax, research, x_quantity, y_quantity, specie=specie, log=log, color=color, geometry=geometry, normalize_to_one=normalize_to_one, interpolate=interpolate)
    
    # Appearance
    if normalize_to_one: ax.set_ylim([0,1]) 
    ax.yaxis.labelpad = 15
    
    # Show the figure   
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0: plt.show() 
    return

#-----------------------------
def subplot_phi_vs_z(
        ax, research, 
        # Quantities to be plotted
        x_quantity="z", 
        y_quantity="phi", 
        specie=0, 
        # Overlay a geometric quantity in gray
        geometry="bmag", 
        # Plotting options
        log=False,
        color=None, 
        interpolate=20, 
        normalize_to_one=False):
    
    # Automate the axis limits and legend
    plot = Plot(); plot.update_legend(loc="upper right");  legend = Legend(ax, plot) 
    axis = Axis(ax, plot, xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.02, logy=log) 
    plot.process_plottingVariables(research); maximum = np.nan; minimum = np.nan;
    tstarts = [np.nan,np.nan]; tends = [np.nan,np.nan] 
    
    # Iterate over the experiments and simulations
    for experiment in research.experiments:
        for simulation in experiment.simulations:
            
            # Style of the data
            style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation)
            style["color"] = color if color!=None else style["color"]
            
            # Read potential(z)
            x = get_field_line_coordinate(x_quantity, simulation)
            y, tstarts, tends = get_quantity_versus_z(y_quantity, simulation, specie, tstarts, tends)
            
            # Interpolate the y-data and normalize the data to have max(y)=1 
            if interpolate: x, y = interpolate_data(x, y, interpolate) 
            if normalize_to_one: y = y/np.max(np.abs(y))

            # Plot potential(z)
            ax.plot(x, y, **style)
            
            # Keep track of the axis limits
            axis.update_axisLimits(x, y)
            maximum = np.nanmax([maximum, np.nanmax(y)])
            minimum = np.nanmin([minimum, np.nanmin(y)])
    
    # Axis labels
    ax.set_xlabel(standardLabels["normalized"][x_quantity]) 
    ax.set_ylabel(get_label(y_quantity, tstarts, tends, normalize_to_one))
        
    # Automatically set the axis limits and legends 
    legend.add_legend()
    axis.rescale_axis()
            
    # Plot geometric quantity
    plot_geometric_quantity(ax, x_quantity, geometry, simulation, maximum, minimum, interpolate, plot)
    return


#-----------------------------
def interpolate_data(x, y, interpolate):
    xnew = np.linspace(x[0], x[-1], num=len(x)*interpolate, endpoint=True) 
    f = interp1d(x, y, kind='cubic')
    y = f(xnew)
    x = xnew
    return x, y

#-----------------------------
def add_extra_legend(ax, color, label, plot):
    fake_line = [Line2D([0], [0], color=color, linewidth=3, linestyle="-")]
    extra_legend = ax.legend(fake_line, [label], labelspacing=0.0, shadow=True, \
        loc="upper left", prop={'size':plot.fontsize},handlelength=plot.handlelength)  
    ax.add_artist(extra_legend) 
    return 

#-----------------------------
def plot_geometric_quantity(ax, x_quantity, geometry, simulation, maximum, minimum, interpolate, plot):
    if geometry=="bmag":
        x = get_field_line_coordinate(x_quantity, simulation)
        bmag = simulation.geometry.bmag[:,0]
        bmag = (bmag-np.min(bmag))/np.max(bmag-np.min(bmag))*(maximum-minimum)+minimum
        if interpolate: x, bmag = interpolate_data(x, bmag, interpolate) 
        ax.plot(x, bmag, color="gray")
        add_extra_legend(ax, "gray", "$B$", plot) 
    elif geometry in ["alpha", "zed", "gradpar", "gds2", "gds21", "gds22", "gds23","gds24", "cvdrift", "gbdrift0", "bmag_psi0"]: 
        x = get_field_line_coordinate(x_quantity, simulation)
        geometry_data = getattr(simulation.geometry, geometry)
        geometry_data = geometry_data/np.max(np.abs(geometry_data))*maximum
        if interpolate: x, geometry_data = interpolate_data(x, geometry_data, interpolate) 
        ax.plot(x, geometry_data, color="gray") 
        add_extra_legend(ax, "gray", geometry+'/max('+geometry+')', plot) 
    elif geometry!=None and geometry!=False:
        print("WARNING: <geometry> = "+geometry+" is not recognized.") 
    return 
         
#-----------------------------
def get_quantity_versus_z(y_quantity, simulation, specie, tstarts, tends):   
        
    # Get the quantity versus (t,z)
    if y_quantity=="phi": 
        vec_quantity = np.abs(simulation.potential.phi_vs_tz.phi)
        vec_time = simulation.potential.phi_vs_tz.t
    if y_quantity=="phi_real": 
        vec_quantity = simulation.potential.phi_vs_tz.phi.real
        vec_time = simulation.potential.phi_vs_tz.t
    if y_quantity=="phi_imag":
        vec_quantity = simulation.potential.phi_vs_tz.phi.imag
        vec_time = simulation.potential.phi_vs_tz.t
    if y_quantity=="phi2":
        vec_quantity = simulation.potential.phi2_vs_tz.phi2
        vec_time = simulation.potential.phi2_vs_tz.t
    if y_quantity=="phi2_zonal":
        vec_quantity = simulation.potential.phi2_vs_tz_zonal.phi2_zonal
        vec_time = simulation.potential.phi2_vs_tz_zonal.t
    if y_quantity=="phi2_nozonal":
        vec_quantity = simulation.potential.phi2_vs_tz_nozonal.phi2_nozonal
        vec_time = simulation.potential.phi2_vs_tz_nozonal.t
    if y_quantity=="g2":
        vec_quantity = simulation.distribution.g2_vs_tsz.g2[:,specie,:] 
        vec_time = simulation.distribution.g2_vs_tsz.t
        vec_quantity[:,-1] = vec_quantity[:,0] # Last value is zero, but the end points along z should match
        
    # Average over the saturated time frame
    tstart, tend = simulation.time.tstart, simulation.time.tend
    vec_quantity = np.mean(vec_quantity[(vec_time>=tstart)&(vec_time<=tend),:], axis=0)  
        
    # Keep track of the time frames
    tstarts[0] = np.nanmin([tstarts[0], tstart]) 
    tstarts[1] = np.nanmax([tstarts[1], tstart]) 
    tends[0] = np.nanmin([tends[0], tend]) 
    tends[1] = np.nanmax([tends[1], tend]) 
    return vec_quantity, tstarts, tends

#-----------------------------
def get_field_line_coordinate(x_quantity, simulation):
    " Get the quantity on the z-axis"
    if x_quantity == "z":     vec_z = simulation.vec.z 
    if x_quantity == "zeta":  vec_z = simulation.vec.zeta
    if x_quantity == "pol":   vec_z = simulation.vec.pol
    if x_quantity == "tor":   vec_z = simulation.vec.tor
    return vec_z

#-----------------------------
def get_label(y_quantity, tstarts, tends, normalize_to_one=False):
    if normalize_to_one==True:
        if y_quantity=="phi":           label = "$|\\varphi|/$max$(|\\varphi|)$"
        if y_quantity=="phi_real":      label = "Re$(\\hat{\\varphi}/$max$(|\\hat{\\varphi}|))$"
        if y_quantity=="phi_imag":      label = "Im$(\\hat{\\varphi}/$max$(|\\hat{\\varphi}|))$"
        if y_quantity=="phi2":          label = "${\\varphi}^2/$max$({\\varphi}^2)$"
        if y_quantity=="phi2_zonal":    label = "$\\sum_{k_x,k_y=0} \\hat{\\varphi}^2_\\mathbf{k}/$max$(\\hat{\\varphi}^2_\\mathbf{k})$"
        if y_quantity=="phi2_nozonal":  label = "$\\sum_{k_x,k_y\\neq0} \\hat{\\varphi}^2_\\mathbf{k}/$max$(\\hat{\\varphi}^2_\\mathbf{k})$"
        if y_quantity=="g2":            label = "$|g|^2/$max$(|g|^2)$" 
    if normalize_to_one==False:
        if y_quantity=="phi":           label = "$|\\varphi|$"
        if y_quantity=="phi_real":      label = "Re$(\\hat{\\varphi})$"
        if y_quantity=="phi_imag":      label = "Im$(\\hat{\\varphi})$"
        if y_quantity=="phi2":          label = "${\\varphi}^2$"
        if y_quantity=="phi2_zonal":    label = "$\\sum_{k_x,k_y=0} \\hat{\\varphi}^2_\\mathbf{k}$"
        if y_quantity=="phi2_nozonal":  label = "$\\sum_{k_x,k_y\\neq0} \\hat{\\varphi}^2_\\mathbf{k}$"
        if y_quantity=="g2":            label = "$|g|^2$" 
    label = "$\\langle$" + label + "$\\rangle_{t="+get_timeFrameString(tstarts, tends)+"}$" 
    return label 

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface
    # Toggles are defined through (name, value, shortoption, longoption, explanation)
    # Options are defined through (name, datatype, shortoption, longoption, explanation)
    bash = Bash(plot_potential_vs_z, __doc__)   
    
    # Check which script to launch      
    input_file = get_firstInputFile(pathlib.Path(os.getcwd())) 
    nonlinear, full_flux_surface = read_nonlinearFullFluxSurfaceFromInputFile(input_file)
    if not nonlinear: os.system("python3 $STELLAPY/plot/linear/potential_vs_z.py "+" ".join(sys.argv[1:])); sys.exit()

    # Data 
    bash.add_option('geometry', 'str', 'g', '', 'Plot geometric quantity from {bmag, gradpar, gds2, gds21, gds22, gds23, gds24, cvdrift, gbdrift0, bmag_psi0}.')  
    
    # Choose the y-quantity
    bash.add_toggleheader("y_quantity")
    bash.add_option('y_quantity', 'str', 'y', '', 'Choose the quantity from {phi, phi_real, phi_imag, phi2, phi2_zonal, phi2_nozonal, g2}.')  
    bash.add_toggle('y_quantity', 'phi', '', '', 'Plot the parallel mode structure of the potential.')  
    bash.add_toggle('y_quantity', 'phi_real', '', '', 'Plot the parallel mode structure of the real part of the potential.')  
    bash.add_toggle('y_quantity', 'phi_imag', '', '', 'Plot the parallel mode structure of the imaginary part of the potential.')  
    bash.add_toggle('y_quantity', 'phi2', '', '', 'Plot the parallel mode structure of the potential squared.') 
    bash.add_toggle('y_quantity', 'phi2_zonal', '', '', 'Plot the parallel mode structure of the zonal components of the potential squared.') 
    bash.add_toggle('y_quantity', 'phi2_nozonal', '', '', 'Plot the parallel mode structure of the non-zonal components of the potential squared.') 
    bash.add_toggle('y_quantity', 'g2', '', '', 'Plot the parallel mode structure of the distribution squared.') 
    bash.add_togglespace()
    
    # Plotting options 
    bash.add_toggleheader("other")
    bash.add_toggle('normalize_to_one', True, 'n', 'normalize_to_one', 'Normalize the data to a maximum of one.')  
    bash.add_toggle('log', True, 'l', 'log', 'Use logaritmic scales.')  
    
    # Research options      
    bash.add_option('folderIsExperiment', True, 'f', 'folderIsExperiment', 'Each folder is an experiment.')  
    
    # Get the arguments and execute the script
    plot_potential_vs_z(**bash.get_arguments())    





    




