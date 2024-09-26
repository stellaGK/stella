"""

#===============================================================================
#                                Plot delta n(z)                               #
#===============================================================================

Based on <folder> create a <research> and plot delta n(z).

The density fluctuations |delta n| are obtained as sqrt(|delta n|^2).
Where |delta n|^2 = sum_{kx,ky} hat{delta n}_{kx,ky}^2.
Note that dens and temp follow bmag, and upar follows cvdrift.


Average density fluctuations in real space
------------------------------------------
In stella the Fourier components hat{delta n}_{kxky} are related to the convential
Fourier components through hat{delta n}_{kxky} = hat{DELTA N}_{kxky}/NxNy. In other words,
the inverse Fourier transformation is defined as
    
    delta n_{xy} = 1/NxNy sum_{kxky} hat{DELTA N}_{kxky} exp(i*kx*x+i*ky*y)
                 = sum_{kxky} hat{delta n}_{kxky} exp(i*kx*x+i*ky*y)
             
As a result, Parseval's theorem is given by

    sum_{xy} |delta n_{xy}|^2 = 1/NxNy sum_{kxky} |hat{DELTA N}_{kxky}|^2
                              = NxNy sum_{kxky} |hat{delta n}_{kxky}|^2  
                          
Therefore, the average density fluctuations squared in real space, is simply 
given by the sum over the (kx,ky) Fourier components as calculated by stella. 

    |delta n|^2 = 1/NxNy sum_{xy} |delta n(x,y)|^2 
            = 1/(NxNy)^2 sum_{kxky} |hat{DELTA N}_{kxky}|^2
            = sum_{kxky} |hat{delta n}_{kxky}|^2

Arguments
--------- 
    folder : pathlib.Path directory where the command has been executed
    x_quantity : {z, pol, tor}
    y_quantity : {dens, temp, upar, dens2, temp2, upar2}
    geometry : {bmag, gradpar, gds2, gds21, gds22, gds23, gds24, cvdrift, gbdrift0, bmag_psi0} overlayed on the background 
    folderIsExperiment : {False, True} 
    interpolate : {int, False} where <int> (e.g. 20) is the interpolation step

Hanne Thienpondt 
19/12/2022

"""

#!/usr/bin/python3   
import sys, os
import pathlib
import numpy as np
import matplotlib as mpl  
import matplotlib.pyplot as plt  
from scipy.interpolate import interp1d

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import get_styleForLinesAndMarkers 
from stellapy.plot.utils.species.recognize_species import recognize_species
from stellapy.plot.geometry.geometry_vs_z import subplot_geometry_vs_z
from stellapy.plot.utils.labels.standardLabels import standardLabels  
from stellapy.plot.utils.style.create_figure import create_figure 
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.Research import create_research 
from stellapy.plot.utils.style import Axis, Legend, Plot
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                                Plot delta n(z)                               #
#===============================================================================

def plot_moment_vs_z(folder, 
        # Quantities to be plotted         
        x_quantity="z", 
        y_quantity="dens2",
        specie=None,
        # Time frame to calculate saturated quantities 
        trange=None,
        # Show a geometric quantity in the background 
        geometry="bmag", 
        # Research
        folderIsExperiment=False, 
        # Plotting
        normalize_to_one=False,
        interpolate=20,  
        title=""):        
    
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder, resolutionScan=False, folderIsExperiment=folderIsExperiment)
    research.update_timeFrame(trange)    
    experiments = research.experiments
    
    # Create a figure for each experiment
    for experiment in experiments: 
            
        # Create a figure
        if y_quantity=="dens": title = "Parallel mode structure of $|\\delta n|$"
        if y_quantity=="temp": title = "Parallel mode structure of $|\\delta T|$"
        if y_quantity=="upar": title = "Parallel mode structure of $|\\delta v_\\parallel|$"
        if y_quantity=="dens2": title = "Parallel mode structure of $|\\delta n|^2$"
        if y_quantity=="temp2": title = "Parallel mode structure of $|\\delta T|^2$"
        if y_quantity=="upar2": title = "Parallel mode structure of $|\\delta v_\\parallel|^2$"
        ax = create_figure(title=title)  
        
        # Plot all the species
        if specie==None:
            species = experiment.simulations[0].vec.species
            colors = ["red", "blue", "green"]; geometry = [geometry, None, None]
            labels = ["Ions", "Electrons", "Impurity"]; ylabels = ""
            for i in species: 
                subplot_moment_vs_z(ax, research, x_quantity, y_quantity, specie=i, color=colors[i], label=labels[i], geometry=geometry[i], interpolate=interpolate, normalize_to_one=normalize_to_one) 
                ylabels += ", "+ax.get_ylabel()
            ax.set_ylim([0, np.max([np.nanmax(l.get_ydata()) for l in ax.get_lines()])])
            ax.set_ylabel(ylabels[2:])
            
        # Plot one specie
        if specie!=None:
            subplot_moment_vs_z(ax, research, x_quantity, y_quantity, specie=specie, interpolate=interpolate, normalize_to_one=normalize_to_one) 
        
    # Appearance 
    ax.yaxis.labelpad = 15
        
    # Show the figure   
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0: plt.show() 
    return

#-----------------------------
def subplot_moment_vs_z(
        ax, research, 
        # Quantities to be plotted    
        x_quantity="zeta", 
        y_quantity="phi", 
        specie=0, 
        # Show a geometric quantity in the background 
        geometry="bmag", 
        # Plotting options
        normalize_to_one=False,
        interpolate=20, 
        color=None, 
        label=None):
    
    # Automate the axis limits and legend
    plot = Plot(); legend = Legend(ax, plot) 
    axis = Axis(ax, plot, xbot_pos=0, ytop_neg=0, ybot_pos=0) 
    plot.process_plottingVariables(research)
    maximum = np.nan; minimum = np.nan;
    
    # Axis labels
    s = recognize_species(research.experiments[0].simulations[0], specie)
    ax.set_xlabel(standardLabels["normalized"][x_quantity]) 
    ax.set_ylabel(standardLabels["normalized"][y_quantity].replace("_{s}","_"+s).replace("_{\\parallel,s}","_{\\parallel,"+s+"}"))
    
    # Iterate over the experiments and simulations
    for experiment in research.experiments:
        for simulation in experiment.simulations:
            
            # Style of the data
            style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation)
            style["color"] = color if color!=None else style["color"]
            style["label"] = label if label!=None else style["label"]
            
            # Read moment(z)
            x = get_field_line_coordinate(x_quantity, simulation)
            y = get_quantity_versus_z(y_quantity, simulation, specie) 
            
            # Interpolate moment(z) 
            if interpolate: x, y = interpolate_data(x, y, interpolate) 
            if normalize_to_one: rescale_factor = np.max(np.abs(y)); y = y/rescale_factor; 
            if normalize_to_one: style["label"] = style["label"]+"$\,/\,"+"{:.2f}".format(rescale_factor)+"$"
            
            # Plot moment(z)
            ax.plot(x, y, **style) 
            
            # Keep track of the axis limits
            axis.update_axisLimits(x, y) 
            maximum = np.nanmax([maximum, np.nanmax(y)])
            minimum = np.nanmin([minimum, np.nanmin(y)])
            
    # Automatically set the axis limits and legends 
    ax.legend(labelspacing=0.0, prop={'size':20}, handlelength=0.8, handletextpad=0.5)
    axis.rescale_axis()
    
    # Plot geometric quantity 
    subplot_geometry_vs_z(ax, research, x_quantity, geometry, minimum=minimum, maximum=maximum, math_mode=True, \
        extra_legend=True, handlelength=plot.handlelength, fontsize=plot.fontsize, interpolate=interpolate, color="gray") 
    return

#-----------------------------
def interpolate_data(x, y, interpolate):
    xnew = np.linspace(x[0], x[-1], num=len(x)*interpolate, endpoint=True) 
    f = interp1d(x, y, kind='cubic')
    y = f(xnew)
    x = xnew
    return x, y

#-----------------------------
def get_quantity_versus_z(y_quantity, simulation, specie):  
    
    # Check whether moments(z) have been written
    if isinstance(simulation.moments.dens2_vs_tsz, float):
        exit_reason = "The moments versus z data has not been written. \n"
        exit_reason += "Please (re-)write this data through: \n     >> write_dataFiles -s mom3D \n"
        exit_program(exit_reason, get_quantity_versus_z, sys._getframe().f_lineno) 
        
    # Get the quantity versus (t,z)
    if "dens" in y_quantity: vec_quantity = simulation.moments.dens2_vs_tsz.dens2[:,specie,:]; vec_time = simulation.moments.dens2_vs_tsz.t 
    if "temp" in y_quantity: vec_quantity = simulation.moments.temp2_vs_tsz.temp2[:,specie,:]; vec_time = simulation.moments.temp2_vs_tsz.t 
    if "upar" in y_quantity: vec_quantity = simulation.moments.upar2_vs_tsz.upar2[:,specie,:]; vec_time = simulation.moments.upar2_vs_tsz.t 
    if "2" not in y_quantity: vec_quantity = np.sqrt(vec_quantity) 
        
    # Average over the saturated time frame
    tstart, tend = simulation.time.tstart, simulation.time.tend
    vec_quantity = np.mean(vec_quantity[(vec_time>=tstart)&(vec_time<=tend),:], axis=0) 
    return vec_quantity

#-----------------------------
def get_field_line_coordinate(x_quantity, simulation):
    " Get the quantity on the z-axis"
    if x_quantity == "z":     vec_z = simulation.vec.z 
    if x_quantity == "zeta":  vec_z = simulation.vec.zeta
    if x_quantity == "pol":   vec_z = simulation.vec.pol
    if x_quantity == "tor":   vec_z = simulation.vec.tor
    return vec_z

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface 
    bash = Bash(plot_moment_vs_z, __doc__)        

    # Species
    bash.add_toggleheader("specie")
    bash.add_option('specie', 'int', 's', '', 'Select the species as "-s 0" for ions and "-s 1" for electrons.') 
    bash.add_toggle('specie', None, '', 'all', 'Plot the moments for all the species. (DEFAULT)') 
    bash.add_toggle('specie', 0, '', 'ions', 'Plot the moments for the ions.') 
    bash.add_toggle('specie', 1, '', 'electrons', 'Plot the moments for the electrons (assuming s=1 for electrons).') 
    bash.add_toggle('specie', 2, '', 'impurity', 'Plot the moments for the impurities (assuming s=2 for impurities).') 
    bash.add_togglespace()
    
    # Choose the y-quantity
    bash.add_toggleheader("y_quantity")
    bash.add_option('y_quantity', 'str', 'y', '', 'Choose the quantity from {dens, temp, upar, dens2, temp2, upar2}.')  
    bash.add_toggle('y_quantity', 'dens', '', '', 'Plot the parallel mode structure of the density fluctuations.')   
    bash.add_toggle('y_quantity', 'dens2', '', '', 'Plot the parallel mode structure of the density fluctuations squared. (DEFAULT)')   
    bash.add_toggle('y_quantity', 'temp', '', '', 'Plot the parallel mode structure of the temperature fluctuations.')   
    bash.add_toggle('y_quantity', 'temp2', '', '', 'Plot the parallel mode structure of the temperature fluctuations squared.')   
    bash.add_toggle('y_quantity', 'upar', '', '', 'Plot the parallel mode structure of the parallel velocity velocity fluctuations.')   
    bash.add_toggle('y_quantity', 'upar2', '', '', 'Plot the parallel mode structure of the parallel velocity fluctuations squared.')   
    bash.add_togglespace()
    
    # Choose the y-quantity
    bash.add_toggleheader("geometry")
    bash.add_option('geometry', 'str', 'g', '', 'Choose the quantity from {bmag, gradpar, gds2, gds21, gds22, gds23, gds24, cvdrift, gbdrift0, bmag_psi0, alpha, zed}.')
    bash.add_toggle('geometry', 'bmag', '', '', 'Plot the magnetic field strength. (DEFAULT)')
    bash.add_toggle('geometry', 'gradpar', '', '', 'Plot gradpar.')      
    bash.add_toggle('geometry', 'gds2', '', '', 'Plot gds2.')      
    bash.add_toggle('geometry', 'gds21', '', '', 'Plot gds21.')      
    bash.add_toggle('geometry', 'gds22', '', '', 'Plot gds22.')      
    bash.add_toggle('geometry', 'gds23', '', '', 'Plot gds23.')      
    bash.add_toggle('geometry', 'gds24', '', '', 'Plot gds24.')      
    bash.add_toggle('geometry', 'cvdrift', '', '', 'Plot cvdrift.')      
    bash.add_toggle('geometry', 'gbdrift0', '', '', 'Plot gbdrift0.')      
    bash.add_toggle('geometry', 'bmag_psi0', '', '', 'Plot bmag_psi0.')      
    bash.add_togglespace() 
    
    # Plotting options 
    bash.add_toggleheader("other")
    bash.add_toggle('normalize_to_one', True, 'n', 'normalize_to_one', 'Normalize the data to a maximum of one.')  
    
    # Research options      
    bash.add_toggle('folderIsExperiment', True, 'f', 'folderIsExperiment', 'Each folder is an experiment.')  
    
    # Get the arguments and execute the script
    plot_moment_vs_z(**bash.get_arguments())    





    




