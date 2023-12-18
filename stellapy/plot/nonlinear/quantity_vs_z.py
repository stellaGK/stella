 
#!/usr/bin/python3   
import sys, os 
import pathlib
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
from stellapy.plot.utils.labels import standardLabels

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import get_styleForLinesAndMarkers  
from stellapy.plot.utils.style.create_figure import create_figure 
from stellapy.simulations.Research import create_research
from stellapy.plot.utils.style import Axis, Legend, Plot 
from stellapy.utils.commandprompt.bash import Bash 

#===============================================================================
#                                 Plot phi(z)                                  #
#===============================================================================

def plot_quantity_vs_z(folder, x_quantity="z", y_quantity="phi2_nozonal", trange=None, 
        log=False, interpolation_step=None): 
    
    # Create <simulations> based on the given <folder>
    research = create_research(folders=folder)  
    
    # Update the time frame
    if trange!=None:    
        for experiment in research.experiments:
            for simulation in experiment.simulations:  
                if trange=="default": simulation.time.update_timeFrame(["50%%", "100%%"])
                if trange!="default": simulation.time.update_timeFrame(trange) 
            
    # Create a figure  
    ax = create_figure()    
     
    # Plot phi2(z) 
    subplot_quantity_vs_z(ax, research, x_quantity, y_quantity, log=log, interpolation_step=interpolation_step)  
                    
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    plt.show()
    return

#---------------------------------------- 
def subplot_quantity_vs_z(ax, research, x_quantity="z", y_quantity="phi2", specie=0, log=False, interpolation_step=None):
    
    # Automate the axis limits and legend
    plot = Plot(); legend = Legend(ax, plot) 
    axis = Axis(ax, plot, xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.1, log=log) 
    
    # Iterate over the simulations
    for experiment in research.experiments:
        for simulation in experiment.simulations:
            
            # Get the lines and marker styles for the experiments and simulations
            style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation)
            
            # Get the quantity along the field lines
            if x_quantity=="z":     x = simulation.vec.z
            if x_quantity=="pol":   x = simulation.vec.pol
            if x_quantity=="tor":   x = simulation.vec.tor
            
            # Get the quantity to be plotted
            if y_quantity=="phi2":          y = simulation.potential.phi2_vs_tz.phi2;                   t = simulation.potential.phi2_vs_tz.t
            if y_quantity=="phi2_nozonal":  y = simulation.potential.phi2_vs_tz_nozonal.phi2_nozonal;   t = simulation.potential.phi2_vs_tz.t
            if y_quantity=="phi2_zonal":    y = simulation.potential.phi2_vs_tz_zonal.phi2_zonal;       t = simulation.potential.phi2_vs_tz.t
            if y_quantity=="dens2":         y = simulation.moments.dens2_vs_tsz.dens2[:,specie,:];      t = simulation.moments.dens2_vs_tsz.t
            if y_quantity=="temp2":         y = simulation.moments.temp2_vs_tsz.dens2[:,specie,:];      t = simulation.moments.temp2_vs_tsz.t
            if y_quantity=="upar2":         y = simulation.moments.upar2_vs_tsz.dens2[:,specie,:];      t = simulation.moments.upar2_vs_tsz.t

            # Get rid of the time dimensions
            y = np.mean(y[(t>=simulation.time.tstart) & (t<=simulation.time.tend)], axis=0)
            
            # Interpolate the data
            if interpolation_step!=None:
                f = interp1d(x, y, kind='cubic')
                x = np.linspace(x[0], x[-1], num=len(x)*interpolation_step, endpoint=True)
                y = f(x)
                
            # Plot phi2(z)
            ax.plot(x, y, **style)
            
            # Keep track of the axis limits
            axis.update_axisLimits(x, y)
            
    # Axis labels
    ax.set_xlabel(standardLabels["normalized"]["z"])
    ax.set_ylabel(standardLabels["normalized"][y_quantity])

    # Automatically set the axis limits and legends 
    legend.add_legend()
    axis.rescale_axis()
    return


#---------------------------------------- 
def plotline_quantity_vs_z(ax, simulation, x_quantity="z", y_quantity="phi2", specie=0,
        color="navy", label="", marker="", ls="-", interpolation_step=None, normalize_to_one=False):  
    
    # Get the quantity along the field lines
    if x_quantity=="z":     x = simulation.vec.z
    if x_quantity=="pol":   x = simulation.vec.pol
    if x_quantity=="tor":   x = simulation.vec.tor
    
    # Get the quantity to be plotted
    if y_quantity=="phi2":          y = simulation.potential.phi2_vs_tz.phi2;                   t = simulation.potential.phi2_vs_tz.t
    if y_quantity=="phi2_nozonal":  y = simulation.potential.phi2_vs_tz_nozonal.phi2_nozonal;   t = simulation.potential.phi2_vs_tz.t
    if y_quantity=="phi2_zonal":    y = simulation.potential.phi2_vs_tz_zonal.phi2_zonal;       t = simulation.potential.phi2_vs_tz.t
    if y_quantity=="dens2":         y = simulation.moments.dens2_vs_tsz.dens2[:,specie,:];      t = simulation.moments.dens2_vs_tsz.t
    if y_quantity=="temp2":         y = simulation.moments.temp2_vs_tsz.dens2[:,specie,:];      t = simulation.moments.temp2_vs_tsz.t
    if y_quantity=="upar2":         y = simulation.moments.upar2_vs_tsz.dens2[:,specie,:];      t = simulation.moments.upar2_vs_tsz.t

    # Get rid of the time dimensions
    y = np.mean(y[(t>=simulation.time.tstart) & (t<=simulation.time.tend)], axis=0)
    
    # Interpolate the data
    if interpolation_step!=None:
        f = interp1d(x, y, kind='cubic')
        x = np.linspace(x[0], x[-1], num=len(x)*interpolation_step, endpoint=True)
        y = f(x)
        
    # Normalize the data
    if normalize_to_one: label += ": $y_\\text{max} = "+"{:.2e}".format(np.max(y))+"$"
    if normalize_to_one: y = y/np.max(y)
        
    # Plot phi2(z)
    ax.plot(x, y, color=color, marker=marker, label=label, ls=ls)   
    
    # Axis labels
    ax.set_xlabel(standardLabels["normalized"]["z"])
    ax.set_ylabel(standardLabels["normalized"][y_quantity])
    return x, y


#===============================================================================
#                             RUN AS BASH COMMAND                              #
#=============================================================================== 
 
if __name__ == "__main__" and False: 
    
    # Launch the bash interface
    bash = Bash(plot_quantity_vs_z, __doc__)    
    
    # Get the arguments and execute the script
    plot_quantity_vs_z(**bash.get_arguments())   

################################################################################
#                                  DEBUG MODE                                  #
################################################################################
    
if __name__ == "__main__":
    import pathlib
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/zeff1.000001_m12.011_zz6.0_fprimz6.0_3kinspecies")
    plot_quantity_vs_z(folder, y_quantity="dens2", log=False, interpolation_step=20)
    sys.exit()
     
     

