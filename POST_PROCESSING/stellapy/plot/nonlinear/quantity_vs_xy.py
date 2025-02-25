"""

Plot potential(x,y).

Hanne Thienpondt 
25/02/2025

"""

#!/usr/bin/python3  
import copy
import sys, os
import pathlib
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp2d
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.calculations.calculate_inverseFourierTransform import calculate_inverseFourierTransform 
from stellapy.plot.utils.labels.add_timeAndZFramesToLabel import add_timeAndZFramesToLabel  
from stellapy.plot.utils.style.create_figure import update_figure_style
from stellapy.plot.utils.surface.get_colorMap import get_colorMap
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.Research import create_research
from stellapy.utils.files.ensure_dir import ensure_dir 
from stellapy.plot.utils.labels import standardLabels
from stellapy.utils.commandprompt.bash import Bash 

#===============================================================================
#                               Plot phi(x,y)                                #
#===============================================================================

def plot_quantity_vs_xy(folder, z_quantity="phi", specie=None, tstart=None, tend=None, z=None, 
        remove_zonal_modes=True, log=False, crange=None, ordersOfMagnitude=2, interpolation_step=20, show=True): 
    
    # Create <simulations> based on the given <folder>
    research = create_research(folders=folder) 
  
    # Plot each simulation separately 
    for experiment in research.experiments:
        for simulation in experiment.simulations: 
            
            # Moments can be plotted for each species
            nspecies = simulation.input.nspec if (specie==None and "phi" not in z_quantity) else 1
     
            # Create a figure   
            fig = plt.figure(figsize=(6*nspecies, 5)); axes = []
            grid_specifications = gridspec.GridSpec(1, nspecies)
            if nspecies==1: grid_specifications.update(top=0.9, left=0.15, right=0.85, bottom=0.15)
            if nspecies==2: grid_specifications.update(top=0.93, left=0.08, right=0.97, bottom=0.1, wspace=0.2, hspace=0.3)
            if nspecies==3: grid_specifications.update(top=0.9, left=0.07, right=0.95, bottom=0.15, wspace=0.4, hspace=0.3)
            for i in range(nspecies): axes.append(plt.subplot(grid_specifications[i]))
            update_figure_style(fig, axes)
             
            # Plot phi2(t)
            for i in range(nspecies):
                if specie!=None: i = specie
                subplot_quantity_vs_xy(axes[i], simulation, z_quantity, specie=i, t_range=[tstart, tend], z=z, interpolation_step=interpolation_step,
                    remove_zonal_modes=remove_zonal_modes, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
                    
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    if show: plt.show()
    return fig

#---------------------------------------- 
def subplot_quantity_vs_xy(ax, simulation, z_quantity, specie, t_range=None, z=None, interpolation_step=None,
        remove_zonal_modes=False, log=False, crange=None, ordersOfMagnitude=2):

    # Get the data 
    z_vs_xy, tstart, tend = get_realSpaceData(simulation, z_quantity, specie, t_range, z, remove_zonal_modes) 
    x, y = simulation.vec.x, simulation.vec.y
    
    # Interpolate the data
    if interpolation_step:    
        function = interp2d(x, y, z_vs_xy.T, kind='cubic')
        x = np.linspace(x[0], x[-1], int(len(x))*interpolation_step)
        y = np.linspace(y[0], y[-1], int(len(y))*interpolation_step)
        z_vs_xy = function(x,y).T              

    # Value x=0 is plotted from x=0 to the next x, correct this by shifting x half a tile left
    x = x-(x[-1]-x[-2])/2; y = y-(y[-1]-y[-2])/2
    x = list(x); x.append(x[-1]+(x[-1]-x[-2]))
    y = list(y); y.append(y[-1]+(y[-1]-y[-2]))
    
    # For quantities squared use a different map
    squared = True if "2" in z_quantity else False
    
    # Get the color map
    if squared: cmap = copy.copy(plt.get_cmap("jet")) 
    if not squared: cmap = get_colorMap(color_map="red-white-blue")
    cmap.set_bad(color='black')  
    
    # Get the range of the color bar
    if log:
        z_vs_xy[z_vs_xy < 1.E-25] = np.NaN
        crange = [ np.nanmin(abs(z_vs_xy)), np.nanmax(abs(z_vs_xy)) ]  
        if ordersOfMagnitude: crange[0] = 10**(np.log10(crange[1])-ordersOfMagnitude)
        norm = LogNorm(vmin=crange[0], vmax=crange[1])  
    if not log:
        norm = None
        if crange!=None: 
            vmin = crange[0]; vmax = crange[1]
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        if not squared: 
            vmin = np.nanmin(z_vs_xy, axis=(0,1)); vmax = np.nanmax(z_vs_xy, axis=(0,1)) 
            vmin = np.min([vmin, -vmax]); vmax = np.max([vmax, -vmin]) 
            norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        
    # Plot the surface z(x,y)  
    img = ax.pcolormesh(x, y, np.real(z_vs_xy).T, cmap=cmap, norm=norm) 
        
    # Add colorbar 
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    cbar = plt.colorbar(img, cax=cax) 
    cbar.update_ticks(); cbar.set_label("")
        
    # Add the colorbar label as title
    nozonal = "_nozonal" if remove_zonal_modes==True else ""
    z_label = add_timeAndZFramesToLabel(standardLabels["normalized"][z_quantity+nozonal].replace("\\hat",""), tstart, tend, z) 
    ax.set_title(z_label.replace("_s","_"+str(specie)))  
    
    # Change the axis style
    if not log:
        cbar.formatter.set_powerlimits((0,0))
        cbar.formatter.set_scientific(True)
        cbar.ax.yaxis.set_offset_position('left')
        cbar.update_ticks()  

    # Set axis
    ax.set_xlabel(standardLabels["normalized"]["x"])
    ax.set_ylabel(standardLabels["normalized"]["y"])
    ax.set_xlim(xmin=np.min(x), xmax=np.max(x))
    ax.set_ylim(ymin=np.min(y), ymax=np.max(y)) 
    return ax, cbar

#-------------------------------
def get_realSpaceData(simulation, z_quantity, specie, t_range, z, remove_zonal_modes=False): 
    
    # Try to read the data from the big 4D files
    try: 
        z_vs_kxky, tstart, tend = get_fourierSpaceDatafrom4DFiles(simulation, z_quantity, specie, z, t_range)
        
    # If this fails, read the saturated files instead.
    except:  
        if (t_range[0]!=None and t_range[0]==t_range[1]) or (t_range[0]!=None and t_range[1]==None): 
            dimension = "4D" if z==None else "5D"
            exit_reason = "A specific time instance is chosen but the data files were not available.\n"
            if "phi" in z_quantity: exit_reason += "Please write the required files through: write_dataFiles -s pot"+dimension
            elif z_quantity in ["dens", "temp", "upar"]: exit_reason += "Please write the required files through: write_dataFiles -s mom"+dimension
            exit_program(exit_reason, get_realSpaceData, sys._getframe().f_lineno)
        z_vs_kxky, tstart, tend = get_fourierSpaceDataFromSaturatedFiles(simulation, z_quantity, specie, z)
        
    # Check the time frame
    if t_range!=None:
        if t_range[0]!=None:
            if tstart!=t_range[0] or (t_range[1]!=-1 and tend!=t_range[1]):
                dimension = "4D" if z==None else "5D"
                print("WARNING: the "+dimension+" files were not written, the data is plotted from")
                print("the saturated file instead. If you want control over the time frame,")
                if "phi" in z_quantity: print("write the data files through >> write_dataFiles -s pot"+dimension)
                elif z_quantity in ["dens", "temp", "upar"]: print("write the data files through >> write_dataFiles -s mom"+dimension)
                else: exit_program("Add other options", get_realSpaceData, sys._getframe().f_lineno)
        
    # Remove the zonal modes
    if remove_zonal_modes==True: z_vs_kxky[:,0] = 0
    
    # Calculate the quantity in real space
    z_vs_xy = calculate_inverseFourierTransform(simulation, z_vs_kxky, axis_kx=0, axis_ky=1) 
    
    # If we want the quantity squared, square it after the IFFT
    if "2" in z_quantity: z_vs_xy = np.abs(z_vs_xy)**2
    return z_vs_xy, tstart, tend

#-------------------------------
def get_fourierSpaceDatafrom4DFiles(simulation, z_quantity, specie, z, t_range):
    
    # Averaged over the field line
    if z==None:
        
        # Get the potential
        if "phi" in z_quantity: 
            z_vs_tkxky = simulation.potential.phi_vs_tkxky.phi[:,:,:]
            t = simulation.potential.phi_vs_tkxky.t
            
        # Get the moments
        if z_quantity in ["dens", "temp", "upar"]:
            t = simulation.moments.dens_vs_tskxky.t
            if z_quantity=="dens": z_vs_tkxky = simulation.moments.dens_vs_tskxky.dens[:,specie,:,:]
            if z_quantity=="temp": z_vs_tkxky = simulation.moments.temp_vs_tskxky.temp[:,specie,:,:]
            if z_quantity=="upar": z_vs_tkxky = simulation.moments.upar_vs_tskxky.upar[:,specie,:,:]
            
    # Specific point along the field line
    if z!=None:
    
        # Get the potential
        if "phi" in z_quantity: 
            z_vs_tzkxky = simulation.potential.phi_vs_tzkxky.phi[:,:,:,:]
            t = simulation.potential.phi_vs_tzkxky.t
            
        # Get the moments
        if z_quantity in ["dens", "temp", "upar"]:
            t = simulation.moments.dens_vs_tskxky.t
            if z_quantity=="dens": z_vs_tzkxky = simulation.moments.dens_vs_tszkxky.dens[:,specie,:,:,:]
            if z_quantity=="temp": z_vs_tzkxky = simulation.moments.temp_vs_tszkxky.temp[:,specie,:,:,:]
            if z_quantity=="upar": z_vs_tzkxky = simulation.moments.upar_vs_tszkxky.upar[:,specie,:,:,:]
            
        # Grab the specific z-point
        iz = np.where(simulation.vec.z>=z)[0][0]
        z_vs_tkxky = z_vs_tzkxky[:,iz,:,:]
        
    # Average out the time dimension
    if t_range==None:
        t_range = simulation.time.t_range
    if t_range[1]!=None:
        z_vs_kxky = np.mean(z_vs_tkxky[(t>=t_range[0])&(t<=t_range[1]),:],axis=0) 
        
    # Grab a specific time point
    if t_range[1]==None:
        try: it = np.where(t>=t_range[0])[0][0]
        except: exit_program("The chosen time instance t = "+str(t_range[0])+" does not exist.", get_fourierSpaceDatafrom4DFiles, sys._getframe().f_lineno)
        t_range[0] = t[it]
        z_vs_kxky = z_vs_tkxky[it,:,:] 
    return z_vs_kxky, t_range[0], t_range[1]
        
#-------------------------------
def get_fourierSpaceDataFromSaturatedFiles(simulation, z_quantity, specie, z): 
    
    # Get the potential
    if "phi" in z_quantity: 
        z_vs_zkxky = simulation.saturated.phi_vs_zkxky.phi[:,:,:] 
        
    # Get the moments
    if z_quantity in ["dens", "temp", "upar"]: 
        if z_quantity=="dens": z_vs_zkxky = simulation.saturated.dens_vs_szkxky.dens[specie,:,:,:]
        if z_quantity=="temp": z_vs_zkxky = simulation.saturated.temp_vs_szkxky.temp[specie,:,:,:]
        if z_quantity=="upar": z_vs_zkxky = simulation.saturated.upar_vs_szkxky.upar[specie,:,:,:]
        
    # Average away the z-dimension or select a specific z
    if z==None:
        dl_over_B = simulation.geometry.dl_over_B 
        z_vs_kxky = np.sum(z_vs_zkxky[:,:,:]*dl_over_B[:,np.newaxis,np.newaxis], axis=0) 
    if z!=None:
        iz = np.where(simulation.vec.z>=z)[0][0]
        z_vs_kxky = z_vs_zkxky[iz,:,:]
        if z!=simulation.vec.z[iz]: print("WARNING: Choosen value of z is "+str(z)+", selected value of z is "+str(simulation.vec.z[iz]))
        
    # Read the time range over which the saturated file has been averaged
    t_range = simulation.saturated.trange
    return z_vs_kxky, t_range[0], t_range[1]    

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#=============================================================================== 
 
if __name__ == "__main__" and False: 
    
    # Launch the bash interface
    bash = Bash(plot_quantity_vs_xy, __doc__)    
    
    # Get the arguments and execute the script
    plot_quantity_vs_xy(**bash.get_arguments())   

################################################################################
#                               STELLAPY EXAMPLE                               #
################################################################################
    
if __name__ == "__main__":    
    
    # Example output file
    stella_folder = os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent.parent)
    figure_folder = stella_folder / pathlib.Path("stellapy_figures"); ensure_dir(figure_folder)
    data_folder = stella_folder / pathlib.Path("DOCUMENTATION/example_output")

    # Plot potential(x,y)
    fig = plot_quantity_vs_xy(data_folder, z_quantity='phi', show=False)
    fig.savefig(figure_folder / "potential_vs_xy.png", format='png', dpi=500)
    print("Saved potential(x,y)")  
     
     

