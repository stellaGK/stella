"""
 
#===============================================================================
#                             Plot potential(z)                                #
#=============================================================================== 
 
For each <simulation> in <research>, created for <folder>, the parallel mode 
structure of the potential squared will be shown for each (kx,ky) mode.  
 
Arguments
--------- 
    folder : pathlib.Path directory where the command has been executed
    x_quantity : {z, pol, tor}
    y_quantity : {phi2, phi, phi_real, phi_imag}
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
from stellapy.data.input.read_inputFile import read_nonlinearFullFluxSurfaceFromInputFile
from stellapy.plot.geometry.geometry_vs_z import subplot_geometry_vs_z
from stellapy.utils.files.get_firstInputFile import get_firstInputFile 
from stellapy.plot.utils.labels.standardLabels import standardLabels  
from stellapy.plot.utils.data.get_rangesKxKy import get_rangesKxKy
from stellapy.plot.utils.style.create_figure import create_figure  
from stellapy.simulations.Research import create_research 
from stellapy.utils.commandprompt.bash import Bash
from stellapy.plot.utils.style import Axis, Plot

#===============================================================================
#                             Plot potential(z)                                #
#=============================================================================== 
 
def plot_potential_vs_z(
        folder, 
        # Quantities to be plotted
        x_quantity="z", 
        y_quantity="phi2", 
        # Overlay a geometric quantity in gray
        geometry="bmag", 
        # Selected modes
        modes_id="unstable", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999], 
        # Research options
        folderIsExperiment=False,
        ignore_resolution=True,
        # Plotting options
        interpolate=20, 
        number_of_decimal_places=2):        
      
    # Create <research> based on the given <folder>
    research = create_research(folders=folder, folderIsExperiment=folderIsExperiment, ignore_resolution=ignore_resolution)
    simulations = [s for experiment in research.experiments for s in experiment.simulations]     
             
    # Create a figure for each simulation
    for simulation in simulations: 
      
        # Create a figure
        if y_quantity=="phi": title = "Parallel mode structure of $|$potential$|$"
        if y_quantity=="phi_real": title = "Parallel mode structure of real(potential)"
        if y_quantity=="phi_imag": title = "Parallel mode structure of imag(potential)"
        if y_quantity=="phi2": title = "Parallel mode structure of the potential squared" 
        if len(simulations)>1: title += ": "+simulation.line_label
        if "k_y" in title: title = title.split("$k_y")[0]+";".join(title.split("$k_y")[-1].split(";")[1:])
        ax = create_figure(title=title)  
         
        # Plot phi2(z)
        subplot_potential_vs_z(ax, simulation, x_quantity=x_quantity, y_quantity=y_quantity, modes_id=modes_id,
            geometry=geometry, kx_range=kx_range, ky_range=ky_range, interpolate=interpolate, fontsize=18, number_of_decimal_places=number_of_decimal_places)
         
        # Appearance 
        ax.yaxis.labelpad = 15
     
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0: plt.show() 
    return
 
#-----------------------------
def subplot_potential_vs_z( 
        ax, simulation, 
        # Quantities to be plotted
        x_quantity="z", 
        y_quantity="phi2", 
        # Overlay a geometric quantity in gray
        geometry="bmag", 
        # Selected modes
        modes_id="unstable", 
        kx_range=[-9999,9999], 
        ky_range=[-9999,9999],  
        # Plotting options
        number_of_decimal_places=2,
        interpolate=20, 
        fontsize=20): 
     
    # Automate the axis limits  
    axis = Axis(ax, Plot(), xbot_pos=0, ybot_pos=0, overshoot_y=1) 
     
    # Color map for modes    
    selected_modes = get_filterForSelectedModes(simulation, modes_id, kx_range, ky_range) 
    number_of_modes = np.sum(selected_modes.astype(int), axis=(0,1)) 
    colors = plt.cm.get_cmap('jet')(np.linspace(0,1,number_of_modes)); plot_i=0
         
    # Get the coordinate along the field line
    x = get_field_line_coordinate(x_quantity, simulation)
    if interpolate: xnew = np.linspace(x[0], x[-1], num=len(x)*interpolate, endpoint=True); xold = x
             
    # Iterate over the modes 
    for ikx in range(simulation.dim.kx):
        for iky in range(simulation.dim.ky): 
            
            # Check whether it needs to be plotted
            if selected_modes[ikx, iky]==False: continue
             
            # Label for the mode  
            ky = simulation.vec.ky[iky]
            string_ky = "{:.2f}".format(ky) if round(ky,2)==float(ky) else ("{:."+str(number_of_decimal_places)+"f}").format(ky)
            label = "$k_y\\rho_i = "+ string_ky + "$"
         
            # Read phi2(z)
            y = get_quantity_versus_z(y_quantity, simulation, ikx, iky)
             
            # Interpolate the y-data
            if interpolate: f = interp1d(xold, y, kind='cubic'); y = f(xnew); x = xnew
             
            # Plot phi2(z) and keep track of the axis limits
            ax.plot(x, y, color=colors[plot_i], label=label)
            axis.update_axisLimits(x, y); plot_i+=1 
            
    # Only continue if we plotted data
    if len(ax.get_lines())!=0: 
                 
        # Labels
        ax.set_xlabel(standardLabels["normalized"][x_quantity]) 
        ax.set_ylabel(get_label(y_quantity)) 
        
        # Limits
        axis.rescale_axis()
         
        # Make sure our axis limits are [0,1] or [-1,1] 
        if y_quantity in ["phi_real", "phi_imag"]: ax.set_ylim([-1,1]) 
        if y_quantity not in ["phi_real", "phi_imag"]: ax.set_ylim([0,1]) 
        
        # Plot geometric quantity 
        subplot_geometry_vs_z(ax, simulation, x_quantity, geometry, minimum=0, maximum=1, math_mode=True, \
            extra_legend=True, interpolate=interpolate, color="gray", loc="upper right") 
         
        # Legend
        ax.legend(labelspacing=0.1, handlelength=1, prop={'size':fontsize}, ncol=np.max([int(number_of_modes/10),1]), loc="upper left")
    
    return  
 
#-----------------------------
def get_quantity_versus_z(y_quantity, simulation, ikx, iky):  
    if "2" not in y_quantity:
        phi = simulation.potential.phi_vs_zkxky.phi[:,ikx,iky]
        if y_quantity=="phi": phi = np.abs(phi)/np.max(np.abs(phi))
        if y_quantity=="phi_real": phi = (phi/np.max(np.abs(phi))).real 
        if y_quantity=="phi_imag": phi = (phi/np.max(np.abs(phi))).imag 
        return phi
    elif "2" in y_quantity:
        phi2 = simulation.potential.phi2_vs_zkxky.phi2[:,ikx,iky] 
        return phi2/np.max(phi2)
 
#-----------------------------
def get_field_line_coordinate(x_quantity, simulation):
    " Get the quantity on the z-axis"
    if x_quantity == "z":     vec_z = simulation.vec.z  
    if x_quantity == "pol":   vec_z = simulation.vec.pol
    if x_quantity == "tor":   vec_z = simulation.vec.tor
    return vec_z
 
#-----------------------------
def get_label(y_quantity):
    if y_quantity=="phi":           label = "$|\\hat{\\varphi}|/$max$(|\\hat{\\varphi}|)$"
    if y_quantity=="phi_real":      label = "Re$(\\hat{\\varphi}/$max$(|\\hat{\\varphi}|))$"
    if y_quantity=="phi_imag":      label = "Im$(\\hat{\\varphi}/$max$(|\\hat{\\varphi}|))$"
    if y_quantity=="phi2":          label = "${\\varphi}^2/$max$({\\varphi}^2)$" 
    return label
 
#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
  
if __name__ == "__main__":
    
    # Create a bash-like interface
    bash = Bash(plot_potential_vs_z, __doc__)   
     
    # Check which script to launch      
    input_file = get_firstInputFile(pathlib.Path(os.getcwd())) 
    nonlinear, full_flux_surface = read_nonlinearFullFluxSurfaceFromInputFile(input_file)
    if nonlinear: os.system("python3 $STELLAPY/plot/nonlinear/potential_vs_z.py "+" ".join(sys.argv[1:])); sys.exit() 
      
    # Adjust the range of wavenumbers
    bash.add_optionheader("wavenumbers")
    bash.add_option('kxmin', 'float', '', '', 'Minimum kx.')   
    bash.add_option('kxmax', 'float', '', '', 'Maximum kx.')  
    bash.add_option('kymin', 'float', '', '', 'Minimum ky.')   
    bash.add_option('kymax', 'float', 'k', '', 'Maximum ky.') 
    bash.add_optionspace()
    
    # Choose the x-quantity
    bash.add_toggleheader("x_quantity") 
    bash.add_toggle('x_quantity', 'z', '', 'zed', 'Plot the parallel mode structure as a function of z.') 
    bash.add_toggle('x_quantity', 'pol', '', '', 'Plot the parallel mode structure as a function of pol.') 
    bash.add_toggle('x_quantity', 'tor', '', '', 'Plot the parallel mode structure as a function of tor.')  
    bash.add_togglespace()  
    
    # Choose the y-quantity
    bash.add_toggleheader("y_quantity")
    bash.add_toggle('y_quantity', 'phi2', '', '', 'Plot the parallel mode structure of the potential squared.') 
    bash.add_toggle('y_quantity', 'phi', '', '', 'Plot the parallel mode structure of the modulus of the potential.') 
    bash.add_toggle('y_quantity', 'phi_real', '', '', 'Plot the parallel mode structure of the real part of the potential.') 
    bash.add_toggle('y_quantity', 'phi_imag', '', '', 'Plot the parallel mode structure of the imaginary part of the potential.')  
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
    bash.add_toggle('modes_id', 'stable', 's', '', 'Plot the stable modes or unconverged modes.')  
    bash.add_toggle('modes_id', 'unstable', '', '', 'Plot the unstable converged modes. (DEFAULT)')   
    bash.add_togglespace()  
    
    # Other options 
    bash.add_toggleheader("other")
    bash.add_optionheader("other")  
    
    # Quantities to be plotted
    bash.add_option('x_quantity', 'str', 'z', '', 'Choose the x-quantity from {z, pol, tor}.')  
    bash.add_option('y_quantity', 'str', 'y', '', 'Choose the y-quantity from {phi, phi_real, phi_imag, phi2}.')  
    bash.add_option('geometry', 'str', 'g', '', 'Choose the geometry from {bmag, gradpar, gds2, gds21, gds22, gds23, gds24, cvdrift, gbdrift0, bmag_psi0, alpha, zed}.')
    
    # Research options 
    bash.add_toggle('ignore_resolution', False, 'i', 'include_resolution', 'Include the resolution (delta t, nzed, ...) between simulations.')   
    bash.add_toggle('folderIsExperiment', False, 'f', 'folderIsExperiment', 'Create an experiment for each folder')
    
    # Choose whether we plot the stable, unstable or all modes
    bash.add_option('modes_id', 'str', 'm', '', 'Choose {unstable, stable, all}.')
    
    # Plotting options
    bash.add_option('number_of_decimal_places', 'int', 'd', 'decimal_places', 'Choose the number of decimal places in the ky label.') 
    
    # Get the bash arguments and execute the script
    args = bash.get_arguments()
    args = get_rangesKxKy(args)
    plot_potential_vs_z(**args)    
 
 
 
 
 
     
 
 


