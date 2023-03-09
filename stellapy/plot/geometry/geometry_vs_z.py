"""

#===============================================================================
#                               Plot geometry(z)                               #
#===============================================================================

Based on <folder> create a <research> and for the first simulation plot geometry(z).

Arguments
--------- 
    folder : pathlib.Path directory where the command has been executed
    x_quantity : {z, pol, tor}
    y_quantity : {bmag, gradpar, gds2, gds21, gds22, gds23, gds24, cvdrift, gbdrift0, bmag_psi0} 
    folderIsExperiment : {False, True} 
    interpolate : {int, False} where <int> (e.g. 20) is the interpolation step

Hanne Thienpondt 
19/12/2022

"""

#!/usr/bin/python3   
import copy
import sys, os
import pathlib
import numpy as np
import matplotlib as mpl  
import matplotlib.pyplot as plt 
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)    
from stellapy.plot.utils.data.update_axisLimits import update_axisLimits 
from stellapy.plot.utils.style.create_figure import create_figure 
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.Research import create_research  
from stellapy.plot.utils.labels import standardLabels
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                               Plot geometry(z)                               #
#===============================================================================

def plot_geometry_vs_z(folder, 
        # Quantities to be plotted         
        x_quantity="z", 
        y_quantity="bmag",    
        # Research
        folderIsExperiment=False, 
        # Rescale the data to fit between [maximum, minimum]
        maximum=None, 
        minimum=0,
        # Plotting 
        math_mode=True,
        extra_legend=False,
        handlelength=1,
        fontsize=18,
        interpolate=20, 
        color="gray", 
        title=""):     
    
    # Make sure <geometry> is a valid choice
    if y_quantity not in ["bmag", "alpha", "zed", "gradpar", "gds2", "gds21", "gds22", "gds23","gds24", "cvdrift", "gbdrift0", "bmag_psi0"]:    
        exit_reason = "Error: <y_quantity> = "+y_quantity+" is not a valid geometric quantity. \n"
        exit_reason += "Choose <y_quantity> from {bmag, gradpar, gds2, gds21, gds22, gds23, gds24, cvdrift, gbdrift0, bmag_psi0, alpha, zed}"
        exit_program(exit_reason, plot_geometry_vs_z, sys._getframe().f_lineno)
        
    # Create a <research> based on the given <folder>
    research = create_research(folders=folder, resolutionScan=False, folderIsExperiment=folderIsExperiment) 
            
    # Create a figure
    if y_quantity=="bmag": title = "Magnetic field strength $B$" 
    ax = create_figure(title=title)  
    
    # Plot geometry(z)
    subplot_geometry_vs_z(ax, research, x_quantity, y_quantity, minimum=minimum, maximum=maximum, math_mode=math_mode, \
        extra_legend=extra_legend, handlelength=handlelength, fontsize=fontsize, interpolate=interpolate, color=color)
                
    # Appearance 
    ax.yaxis.labelpad = 15
    
    # Labels 
    ax.set_xlabel(standardLabels["normalized"][x_quantity])
    ax.set_ylabel(get_labels(y_quantity, math_mode, maximum)) 
    
    # Add the legend
    if not extra_legend: ax.legend(labelspacing=0.0, prop={'size':20}, handlelength=0.8, handletextpad=0.5)
        
    # Show the figure   
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0: plt.show() 
    return


def subplot_geometry_vs_z(
        ax, research, 
        # Quantities to be plotted         
        x_quantity="z", 
        y_quantity="bmag",     
        # Rescale the data to fit between [maximum, minimum]
        maximum=None, 
        minimum=0,
        # Plotting 
        math_mode=True,
        extra_legend=False,
        handlelength=1,
        fontsize=18,
        interpolate=20, 
        loc="upper left",
        color="gray"):   
    
    # Make sure <geometry> is a valid choice
    if y_quantity not in ["bmag", "alpha", "zed", "gradpar", "gds2", "gds21", "gds22", "gds23","gds24", "cvdrift", "gbdrift0", "bmag_psi0"]:    
        if y_quantity==None: return
        exit_reason = "Error: <y_quantity> = "+str(y_quantity)+" is not a valid geometric quantity. \n"
        exit_reason += "Choose <y_quantity> from {bmag, gradpar, gds2, gds21, gds22, gds23, gds24, cvdrift, gbdrift0, bmag_psi0, alpha, zed}"
        exit_program(exit_reason, plot_geometry_vs_z, sys._getframe().f_lineno)
    
    # Keep track of the data
    plotted_geometries = {}; 
    xlims = [np.nan, np.nan]; 
    ylims = [0, np.nan]
    
    # Get the simulations
    if hasattr(research, "experiments"): simulations = [s for experiment in research.experiments for s in experiment.simulations]
    if not hasattr(research, "experiments"): simulations = [research]
        
    # Plot the geometry for each simulation
    for simulation in simulations:
            
        # Get the x-coordinate and y-coordinate
        x = get_field_line_coordinate(x_quantity, simulation) 
        y = getattr(simulation.geometry, y_quantity)
        
        # Make sure there is no alpha dimension 
        if len(np.shape(y))>1: y = y[:,0]
                
        # Only plot each unique geometry once
        for key in plotted_geometries: 
            if np.array_equal(plotted_geometries[key], y): continue
        else: plotted_geometries[simulation.id] = np.array(y)
        
        # Interpolate the data 
        if interpolate: x, y = interpolate_data(x, y, interpolate)  
        
        # Rescale to [minimum, maximum]
        if maximum!=None:  
            if np.all(y>=0):
                y = (y-np.min(y))/np.max(y-np.min(y))
                y = y*(maximum-minimum)+minimum
            else:  
                yold = copy.deepcopy(y)
                y = (y-np.min(y))/np.max(y-np.min(y))   
                y_negative = copy.deepcopy(y)
                y_negative[yold>0] = np.nan   
                ax.plot(x, y_negative, color=color, ls=":") 
                y[yold<0] = np.nan 
                
        # Label
        label = get_labels(y_quantity, math_mode, maximum) 
                
        # Plot geometry(z)
        ax.plot(x, y, color=color, label=label if extra_legend==False else "")
        
        # Put the label as an extra legend
        if extra_legend: add_extra_legend(ax, color, label, fontsize=fontsize, handlelength=handlelength, loc=loc)
        xlims, ylims = update_axisLimits(x, y, xlims, ylims) 
         
    # Axis limits 
    if not extra_legend: 
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)   
    return

#-----------------------------
def get_field_line_coordinate(x_quantity, simulation):
    " Get the quantity on the z-axis"
    if x_quantity == "z":     vec_z = simulation.vec.z 
    if x_quantity == "zeta":  vec_z = simulation.vec.zeta
    if x_quantity == "pol":   vec_z = simulation.vec.pol
    if x_quantity == "tor":   vec_z = simulation.vec.tor
    return vec_z 

#-----------------------------
def add_extra_legend(ax, color, label, fontsize=18, handlelength=1, loc="upper left"):
    fake_line = [Line2D([0], [0], color=color, linewidth=3, linestyle="-")]
    extra_legend = ax.legend(fake_line, [label], labelspacing=0.0, \
        loc=loc, prop={'size':fontsize},handlelength=handlelength)  
    ax.add_artist(extra_legend) 
    return 

#-----------------------------
def interpolate_data(x, y, interpolate):  
    xnew = np.linspace(x[0], x[-1], num=len(x)*interpolate, endpoint=True) 
    f = interp1d(x, y, kind='cubic')
    y = f(xnew)
    x = xnew
    return x, y


#-----------------------------------------------
def get_labels(geometry, math_mode, maximum): 
    if math_mode:
        labels = {
            "bmag" : "$B$",\
            "gradpar" : "$\\nabla_{\\parallel}z$",\
            "gds2" : "$|\\nabla y|^2$",\
            "gds21" : "$\\nabla x \\cdot \\nabla y$",\
            "gds22" : "$|\\nabla x|^2$",\
            "gbdrift" : "$\\mathbf{B} \\times \\nabla B \\cdot \\nabla y$",\
            "gbdrift0" : "$\\mathbf{B} \\times \\nabla B \\cdot \\nabla x$",\
            "cvdrift" : "$\\mathbf{B} \\times \\mathbf{\\kappa} \\cdot \\nabla y$",\
            "cvdrift0" : "$\\mathbf{B} \\times \\mathbf{\\kappa} \\cdot \\nabla x$",\
            "gds23" : "$((\\nabla y \\cdot \\nabla \\zeta)(\\nabla x \\cdot \\nabla y) - (\\nabla x \\cdot \\nabla \\zeta)|\\nabla y|^2)$",\
            "gds24" : "$((\\nabla y \\cdot \\nabla \\zeta)|\\nabla x|^2 - (\\nabla x \\cdot \\nabla \\zeta)(\\nabla x \\cdot \\nabla y))$",\
            }
        if maximum==1:
            labels["bmag"] = "$(B-B_{min})/B_{max}$"
    if not math_mode:
        labels = {
            "bmag" : "bmag",\
            "gradpar" : "gradpar",\
            "gds2" : "gds2",\
            "gds21" : "gds21",\
            "gds22" : "gds22",\
            "gbdrift" : "gbdrift",\
            "gbdrift0" : "gbdrift0",\
            "cvdrift" : "cvdrift",\
            "cvdrift0" : "cvdrift0",\
            "gds23" : "gds23",\
            "gds24" : "gds24",\
            }
    return labels[geometry]

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface 
    bash = Bash(plot_geometry_vs_z, __doc__)        
    
    # Choose the y-quantity
    bash.add_toggleheader("y_quantity")
    bash.add_option('y_quantity', 'str', 'y', '', 'Choose the quantity from {bmag, gradpar, gds2, gds21, gds22, gds23, gds24, cvdrift, gbdrift0, bmag_psi0, alpha, zed}.')
    bash.add_toggle('y_quantity', 'bmag', '', '', 'Plot the magnetic field strength.')
    bash.add_toggle('y_quantity', 'gradpar', '', '', 'Plot gradpar.')      
    bash.add_toggle('y_quantity', 'gds2', '', '', 'Plot gds2.')      
    bash.add_toggle('y_quantity', 'gds21', '', '', 'Plot gds21.')      
    bash.add_toggle('y_quantity', 'gds22', '', '', 'Plot gds22.')      
    bash.add_toggle('y_quantity', 'gds23', '', '', 'Plot gds23.')      
    bash.add_toggle('y_quantity', 'gds24', '', '', 'Plot gds24.')      
    bash.add_toggle('y_quantity', 'cvdrift', '', '', 'Plot cvdrift.')      
    bash.add_toggle('y_quantity', 'gbdrift0', '', '', 'Plot gbdrift0.')      
    bash.add_toggle('y_quantity', 'bmag_psi0', '', '', 'Plot bmag_psi0.')      
    bash.add_togglespace() 
    
    # Plotting options 
    bash.add_toggleheader("other")
    bash.add_option('minimum', 'float', '', '', 'Rescale the data to have the minimum at <minimum>, default is 0.')  
    bash.add_option('maximum', 'float', 'm', '', 'Rescale the data have the maximum at <maximum>, set this to 1.')  
    bash.add_option('interpolate', 'int', 'i', '', 'Choose the interpolation stop, 0 is no interpolation, default is 20.')  
    bash.add_toggle('maximum', 1, 'n', 'normalize_to_one', 'Rescale the data to range from 0 to 1.')  
    
    # Research options      
    bash.add_option('folderIsExperiment', True, 'f', 'folderIsExperiment', 'Each folder is an experiment.')  
    
    # Get the arguments and execute the script
    plot_geometry_vs_z(**bash.get_arguments())    





    




