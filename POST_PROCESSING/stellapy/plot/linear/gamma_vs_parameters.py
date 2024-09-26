"""

#===============================================================================
#           PLOT GAMMA(PARAMETER,PARAMETER) OF THE MOST UNSTABLE MODE          #
#===============================================================================

Plot the stability or growth rate map, gamma(parameter1,parameter2) in <folder>.
For each set of (parameter1,parameter2) a linear scan of modes (kx,ky) should 
have been performed, of which the most unstable mode is selected. The growth rate
(gamma) and frequency (omega) of the most unstable mode is plotted as a function
of the scanned <parameter1> and <parameter2>. 

Arguments
--------- 
    z_quantity : {omega, gamma, ky} 

Hanne Thienpondt
01/09/2022

"""
# Load modules 
import sys, os
import pathlib
import numpy as np
import matplotlib.pyplot as plt  
from matplotlib.ticker import MultipleLocator

# Stellapy package   
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.plot.utils.surface.get_colorMapNormalization import get_colorMapNormalization
from stellapy.plot.utils.data.get_gammaVsParameters import get_gammaVsParameters
from stellapy.plot.utils.data.get_uniqueParameters import get_uniqueParameters
from stellapy.utils.commandprompt.print_progressbar import print_progressbar
from stellapy.plot.utils.style.create_figure import create_figure
from stellapy.plot.utils.surface.get_colorMap import get_colorMap
from stellapy.plot.utils.surface import interpolate_surface
from stellapy.simulations.Research import create_research   
from stellapy.plot.utils.labels import standardLabels 
from stellapy.utils.commandprompt.bash import Bash
 
#===============================================================================
#           PLOT GAMMA(PARAMETER,PARAMETER) OF THE MOST UNSTABLE MODE          #
#===============================================================================

def plot_gamma_vs_parameters(\
            # Simulations
                folder,
            # Parameters on the x-axis and y-axis
                knob1="species_parameters_1",
                key1="fprim",
                knob2="species_parameters_1",
                key2="tprim",
            # Data
                z_quantity='gamma',\
            # Figure 
                ax=None, \
                c_range=None,\
                show_figure=True,\
            # Interpolate 
                step=4,
                boolean_maps=None,
                interpolate=False):
    
    # Create a <research> based on the given <folder>
    print_progressbar(0, 100, prefix = '   Progress:', suffix = 'Creating research object.', length=50)
    knobs_and_keys = {"knob1" : knob1, "key1" : key1, "knob2" : knob2, "key2" : key2}
    research = create_research(folder, variables=4, **knobs_and_keys, progress_start=1, progress_stop=10)
    
    # Create a figure
    if not ax: ax = create_figure(figsize=(7.5, 7), left=0.14, right=0.97, top=0.92, bottom=0.16)   
        
    # Get the parameters
    parameters1 = get_uniqueParameters(research, key1, knob1) 
    parameters2 = get_uniqueParameters(research, key2, knob2)  
         
    # Read the {gamma, omega, ky} data     
    print_progressbar(10, 100, prefix = '   Progress:', suffix = 'Reading gamma(parameter1,parameter2).', length=50)
    gamma, omega, ky = get_gammaVsParameters(research, parameters1, parameters2, knob1, key1, knob2, key2, start=10,stop=80)
          
    # Get the z-data 
    if z_quantity=="gamma": zdata = gamma
    if z_quantity=="omega": zdata = omega
    if z_quantity=="ky":    zdata = ky   
     
    # Interpolate the data z(x,y) and replace zeros and original nans by nan
    parameters1, parameters2, zdata, filter_zero_and_nans = interpolate_data(parameters1, parameters2, zdata, boolean_maps, interpolate, step)
    zdata[filter_zero_and_nans] = np.NaN
    
    # Get the color map "black-jet" for gamma and "red-white-blue" for omega
    vmin, vmax, norm = get_colorMapNormalization(zdata, c_range, center_at_zero=(z_quantity=="omega"))
    cmap = "black-jet" if z_quantity!="omega" else "red-white-blue"
    cmap = get_colorMap(cmap)
     
    # Plot gamma(parameter1,parameter2) or omega(parameter1,parameter2)   
    print_progressbar(90, 100, prefix = '   Progress:', suffix = 'Plotting gamma(parameter1,parameter2).', length=50)
    img = ax.pcolormesh(parameters1, parameters2, zdata.T, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm, shading='auto')  
    cbar = plt.colorbar(img, ax=ax, pad=0.02); ax.fig.cbar = cbar
    if c_range: img.set_clim(c_range[0],c_range[1])  
     
    # Set a label for the colorbar
    if z_quantity=="gamma": cbar.set_label('$\\gamma a/v_{\\mathrm{th},i}$')
    if z_quantity=="omega": cbar.set_label('$\\omega a/v_{\\mathrm{th},i}$')
    if z_quantity=="ky": cbar.set_label('$k_{y}\\rho_i$')
          
    # Set axis and labels
    ax.xaxis.set_major_locator(MultipleLocator(1)) 
    ax.yaxis.set_major_locator(MultipleLocator(1)) 
    ax.set_xlim([np.min(parameters1), np.max(parameters1)])
    ax.set_ylim([np.min(parameters2), np.max(parameters2)])
    ax.set_xlabel(standardLabels["normalized"][key1])
    ax.set_ylabel(standardLabels["normalized"][key2])
        
    # Show the figure 
    print_progressbar(100, 100, prefix = '   Progress:', suffix = 'Finished.', length=50)
    if show_figure: plt.show()
    return ax.fig
 
#-------------------------------------
def interpolate_data(x, y, z, boolean_maps=None, interpolate=True, step=4):
    
    # Interpolate the data
    if interpolate:    
        x, y, z, filter_zero_and_nans = interpolate_surface(x, y, z, step, boolean_maps)  
        x = [p for p in x]
        y = [p for p in y]
          
    # Add an extra parameter because the color is shown from x to x+1, add offset the labels
    if not interpolate:
        filter_zero_and_nans = (np.abs(z) < 1.E-25) | np.isnan(z)
        x = x + [x[-1] + (x[1]-x[0])]
        y = y + [y[-1] + (y[1]-y[0])]  
        x = [p-0.5 for p in x]
        y = [p-0.5 for p in y]

    return x, y, z, filter_zero_and_nans


#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface
    bash = Bash(plot_gamma_vs_parameters, __doc__)   
    
    # Toggle the quantities to be plotted through --omega  
    bash.add_toggle('y_quantity', 'omega', 'Plot the spectra for omega.') 
    
    # Turn on the interpolation
    bash.add_option('interpolate', 'True', 'i', 'Interpolate gamma(parameter1,parameter2).')  
    bash.add_option('step', 'int', 's', 'Interpolation step: xnew = step*xold.') 
    
    # Get the bash arguments and execute the script 
    plot_gamma_vs_parameters(**bash.get_arguments())   



