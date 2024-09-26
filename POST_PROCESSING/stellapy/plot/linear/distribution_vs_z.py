"""
  
#===============================================================================
#                         Plot distribution squared(z)                         #
#=============================================================================== 
  
For each <simulation> in <research>, created for <folder>, the parallel mode 
structure of the distribution squared will be shown for each (kx,ky) mode.  
  
Arguments
--------- 
    folder : pathlib.Path directory where the command has been executed
    x_quantity : {z, pol, tor}
    y_quantity : {g2}
    geometry : {bmag} overlayed on the background
    modes_id : {unstable, stable, all}
    ignore_resolution : {True, False}
    folderIsExperiment : {False, True}
    kx_range : [-9999, 9999]
    ky_range : [-9999, 9999] 
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
from stellapy.data.lineardata.get_filterForSelectedModes import get_filterForSelectedModes
from stellapy.plot.utils.species.recognize_species import recognize_species
from stellapy.plot.geometry.geometry_vs_z import subplot_geometry_vs_z
from stellapy.plot.utils.labels.standardLabels import standardLabels  
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.plot.utils.style.create_figure import create_figure  
from stellapy.simulations.Research import create_research 
from stellapy.utils.commandprompt.bash import Bash
from stellapy.plot.utils.style import Axis, Plot
  
#===============================================================================
#                         Plot distribution squared(z)                         #
#=============================================================================== 
  
def plot_distribution_vs_z(
        folder, 
        # Quantity to be plotted
        x_quantity="z", 
        y_quantity="g2", 
        specie=0, 
        # Overlay a geometric quantity in gray
        geometry="bmag", \
        # Selected modes
        modes_id="unstable", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999], 
        # Research
        ignore_resolution=True, 
        folderIsExperiment=False, 
        # Plotting
        interpolate=20, 
        number_of_decimal_places=2):        
       
    # Create <research> based on the given <folder>
    research = create_research(folders=folder, folderIsExperiment=folderIsExperiment, ignore_resolution=ignore_resolution)   
    simulations = [s for experiment in research.experiments for s in experiment.simulations]   
     
    # Create a figure for each simulation
    for simulation in simulations: 
       
        # Create a figure 
        if y_quantity=="g2": title = "Parallel mode structure of the distribution squared"
        if y_quantity=="g":  title = "Parallel mode structure of the distribution"
        if len(simulations)>1: title += ": "+simulation.line_label
        if "k_y" in title: title = title.split("$k_y")[0]+";".join(title.split("$k_y")[-1].split(";")[1:])
        ax = create_figure(title=title)  
          
        # Plot distribution2(z)
        subplot_potential_vs_z(ax, simulation, x_quantity=x_quantity, y_quantity=y_quantity, 
            specie=specie, geometry=geometry, number_of_decimal_places=number_of_decimal_places, 
            modes_id=modes_id, kx_range=kx_range, ky_range=ky_range, interpolate=interpolate, fontsize=18)
          
        # Appearance 
        ax.yaxis.labelpad = 15
      
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0: plt.show() 
    return
  
#-----------------------------
def subplot_potential_vs_z(
        # Axis
        ax, 
        # Simulation
        simulation, 
        # Quantity to be plotted
        x_quantity="zeta", 
        y_quantity="g2", 
        specie=0, 
        # Also plot a geometric quantity
        geometry="bmag", 
        # Selected modes
        modes_id="unstable", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999], 
        # Plotting
        interpolate=20, 
        fontsize=20, 
        number_of_decimal_places=2):  
      
    # Automate the axis limits  
    axis = Axis(ax, Plot(), xbot_pos=0, ybot_pos=0, overshoot_y=1)
      
    # Color map for modes    
    selected_modes = get_filterForSelectedModes(simulation, modes_id, kx_range, ky_range) 
    number_of_modes = np.sum(selected_modes.astype(int), axis=(0,1)) 
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,number_of_modes)); plot_i=0 
      
    # Get the coordinate along the field line
    x = get_field_line_coordinate(x_quantity, simulation)
    if interpolate: xnew = np.linspace(x[0], x[-1], num=len(x)*20, endpoint=True); xold = x
              
    # Iterate over the modes 
    for ikx in range(simulation.dim.kx):
        for iky in range(simulation.dim.ky): 
             
            # Check whether it needs to be plotted
            if selected_modes[ikx, iky]==False: continue
              
            # Label for the mode  
            ky = simulation.vec.ky[iky]
            string_ky = "{:.2f}".format(ky) if round(ky,2)==float(ky) else ("{:."+str(number_of_decimal_places)+"f}").format(ky)
            label = "$k_y\\rho_i = "+ string_ky + "$"
          
            # Read distribution squared(z)
            y = simulation.distribution.g2_vs_szkxky.g2[specie,:,ikx,iky]/np.max(simulation.distribution.g2_vs_szkxky.g2[specie,:,ikx,iky])
              
            # In the distribution array y[-1]=0 but it should be y[-1] = y[0]
            y[-1] = y[0]
              
            # Interpolate the y-data  
            if interpolate: f = interp1d(xold, y, kind='cubic'); y = f(xnew); x = xnew
              
            # Plot distribution2(z) and keep track of the axis limits
            ax.plot(x, y, color=colors[plot_i], label=label)
            axis.update_axisLimits(x, y); plot_i+=1 
              
    # Only continue if we plotted data
    if len(ax.get_lines())!=0: 
                  
        # Labels
        s = recognize_species(simulation, specie) 
        ax.set_xlabel(standardLabels["normalized"][x_quantity]) 
        ax.set_ylabel(standardLabels["normalized"][y_quantity].replace("_{s}","_"+s)) 
         
        # Limits
        axis.rescale_axis()
          
        # Make sure our axis limits are [0,1] or [-1,1] 
        if y_quantity in ["phi_real", "phi_imag"]: ax.set_ylim([-1,1]) 
        if y_quantity not in ["phi_real", "phi_imag"]: ax.set_ylim([0,1]) 
        
        # Plot geometric quantity 
        subplot_geometry_vs_z(ax, simulation, x_quantity, geometry, minimum=0, maximum=1, math_mode=True, \
            extra_legend=True, interpolate=interpolate, color="gray", loc="upper right") 
          
        # Legend
        ax.legend(labelspacing=0.1, handlelength=1, prop={'size':fontsize}, ncol=int(number_of_modes/10), loc="upper left")
    
    return
 
#-----------------------------
def get_field_line_coordinate(x_quantity, simulation):
    " Get the quantity on the z-axis"
    if x_quantity == "z":     vec_z = simulation.vec.z  
    if x_quantity == "pol":   vec_z = simulation.vec.pol
    if x_quantity == "tor":   vec_z = simulation.vec.tor
    return vec_z
  
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
   
if __name__ == "__main__":
    
    # Create a bash-like interface
    bash = Bash(plot_distribution_vs_z, __doc__)   
      
    # Adjust the range of wavenumbers
    bash.add_optionheader("wavenumbers")
    bash.add_option('kxmin', 'float', '', '', 'Minimum kx.')   
    bash.add_option('kxmax', 'float', '', '', 'Maximum kx.')  
    bash.add_option('kymin', 'float', '', '', 'Minimum ky.')   
    bash.add_option('kymax', 'float', 'k', '', 'Maximum ky.') 
    bash.add_optionspace()
    
    # Species
    bash.add_toggleheader("specie")
    bash.add_toggle('specie', 0, '', 'ions', 'Plot the distribution squared for the ions. (DEFAULT)') 
    bash.add_toggle('specie', 1, '', 'electrons', 'Plot the distribution squared for the electrons (assuming s=1 for electrons).') 
    bash.add_toggle('specie', 2, '', 'impurity', 'Plot the distribution squared for the impurities (assuming s=2 for impurities).') 
    bash.add_togglespace()
    
    # Choose the x-quantity
    bash.add_toggleheader("x_quantity") 
    bash.add_toggle('x_quantity', 'z', '', 'zed', 'Plot the parallel mode structure as a function of z.') 
    bash.add_toggle('x_quantity', 'pol', '', '', 'Plot the parallel mode structure as a function of pol.') 
    bash.add_toggle('x_quantity', 'tor', '', '', 'Plot the parallel mode structure as a function of tor.')  
    bash.add_togglespace()  
    
    # Choose the y-quantity
    bash.add_toggleheader("y_quantity")
    bash.add_toggle('y_quantity', 'g2', '', '', 'Plot the parallel mode structure of the distribution squared.') 
    bash.add_togglespace()  
    
    # Choose the y-quantity
    bash.add_toggleheader("geometry")
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
    
    # Choose the modes
    bash.add_toggleheader("modes")
    bash.add_toggle('modes_id', 'all', '', '', 'Plot all the modes.')  
    bash.add_toggle('modes_id', 'stable', '', '', 'Plot the stable modes or unconverged modes.')  
    bash.add_toggle('modes_id', 'unstable', '', '', 'Plot the unstable converged modes. (DEFAULT)')   
    bash.add_togglespace()  
    
    # Other options 
    bash.add_toggleheader("other")
    bash.add_optionheader("other")  
    
    # Quantities to be plotted
    bash.add_option('specie', 'int', 's', '', 'Select the species as "-s 0" for ions and "-s 1" for electrons.')  
    bash.add_option('x_quantity', 'str', 'z', '', 'Choose the x-quantity from {z, pol, tor}.')  
    bash.add_option('y_quantity', 'str', 'y', '', 'Choose the y-quantity from {g2}.')  
    bash.add_option('geometry', 'str', 'g', '', 'Choose the geometry from {bmag, gradpar, gds2, gds21, gds22, gds23, gds24, cvdrift, gbdrift0, bmag_psi0, alpha, zed}.')
    
    # Research options 
    bash.add_toggle('ignore_resolution', False, 'i', 'include_resolution', 'Include the resolution (delta t, nzed, ...) between simulations.')   
    bash.add_toggle('folderIsExperiment', False, 'f', 'folderIsExperiment', 'Create an experiment for each folder')
    
    # Choose whether we plot the stable, unstable or all modes
    bash.add_option('modes_id', 'str', 'm', '', 'Choose {unstable, stable, all}.')
    
    # Plotting options
    bash.add_option('interpolate', 'int', '-', '', 'Interpolate the data along z, choose the interpolation step {0, 20, ...}.') 
    bash.add_option('number_of_decimal_places', 'int', 'd', 'decimal_places', 'Choose the number of decimal places in the ky label.') 
    
    # Get the bash arguments and execute the script
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_distribution_vs_z(**args)    
  
  
  
      
  
  


