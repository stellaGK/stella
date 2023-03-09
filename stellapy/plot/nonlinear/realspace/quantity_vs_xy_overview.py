 
#!/usr/bin/python3   
import sys, os 
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec  

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.plot.nonlinear.potential_vs_time import subplot_potential_vs_time_zonal_contributions
from stellapy.plot.nonlinear.quantity_vs_xy import subplot_quantity_vs_xy
from stellapy.plot.nonlinear.quantity_vs_z import plotline_quantity_vs_z 
from stellapy.plot.utils.style.create_figure import update_figure_style 
from stellapy.simulations.Research import create_research 

#===============================================================================
#                               Plot phi(x,y)                                #
#===============================================================================

def plot_quantity_vs_xy_overview(folder, tspecific=150, 
        interpolation_step=20, ordersOfMagnitude=2, log=False, crange=None): 
    
    # Create <simulations> based on the given <folder>
    research = create_research(folders=folder) 
    simulation = research.experiments[0].simulations[0] 
     
    # Create a figure   
    fig = plt.figure(figsize=(18, 9)); axes = []
    grid_specifications = gridspec.GridSpec(3, 4)
    grid_specifications.update(top=0.95, left=0.05, right=0.95, bottom=0.08, wspace=0.4, hspace=0.4)
    axes.append(plt.subplot(grid_specifications[0:2]))
    axes.append(plt.subplot(grid_specifications[2:4]))
    axes.append(plt.subplot(grid_specifications[4]))
    axes.append(plt.subplot(grid_specifications[5]))
    axes.append(plt.subplot(grid_specifications[6]))
    axes.append(plt.subplot(grid_specifications[7]))
    axes.append(plt.subplot(grid_specifications[8]))
    axes.append(plt.subplot(grid_specifications[9]))
    axes.append(plt.subplot(grid_specifications[10]))
    axes.append(plt.subplot(grid_specifications[11]))
    update_figure_style(fig, axes, font_size=16)
    
    # Plot phi2(t) and its zonal contributions
    subplot_potential_vs_time_zonal_contributions(axes[0], simulation, log=False, fontsize=12)
    
    # Plot phi2(z) and dens2(z) for each species 
    x1, y1 = plotline_quantity_vs_z(axes[1], simulation, x_quantity="z", y_quantity="phi2", label="$|\\phi^2|$", color="black", normalize_to_one=True, interpolation_step=20)
    x2, y2 = plotline_quantity_vs_z(axes[1], simulation, x_quantity="z", y_quantity="dens2", specie=0, label="$(\\delta n_i)^2$", color="crimson", normalize_to_one=True, interpolation_step=20)
    x3, y3 = plotline_quantity_vs_z(axes[1], simulation, x_quantity="z", y_quantity="dens2", specie=1, label="$(\\delta n_e)^2$", color="navy", normalize_to_one=True, interpolation_step=20)
    try: x4, y4 = plotline_quantity_vs_z(axes[1], simulation, x_quantity="z", y_quantity="dens2", specie=2, label="$(\\delta n_\\text{imp})^2$", color="ForestGreen", normalize_to_one=True, interpolation_step=20)
    except: x4 = x3; y4 = y3
    axes[1].legend(labelspacing=0.0, handlelength=1, prop={'size':12})
    axes[1].set_xlim([np.min(np.array([x1,x2,x3,x4])), np.max(np.array([x1,x2,x3,x4]))])
    axes[1].set_ylim([np.min(np.array([y1,y2,y3,y4])), np.max(np.array([y1,y2,y3,y4]))*1.1])
    axes[1].set_ylim([0,1]); #axes[1].set_yscale('log')
    
    # Plot phi(x,y) 
    subplot_quantity_vs_xy(axes[2], simulation, "phi", specie=0, t_range=[tspecific,None], z=None, interpolation_step=interpolation_step,
        remove_zonal_modes=False, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
    subplot_quantity_vs_xy(axes[3], simulation, "phi", specie=0, t_range=[tspecific,None], z=0, interpolation_step=interpolation_step,
        remove_zonal_modes=False, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
    subplot_quantity_vs_xy(axes[4], simulation, "phi", specie=0, t_range=[tspecific,None], z=None, interpolation_step=interpolation_step,
        remove_zonal_modes=True, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
    subplot_quantity_vs_xy(axes[5], simulation, "phi", specie=0, t_range=[tspecific,None], z=0, interpolation_step=interpolation_step,
        remove_zonal_modes=True, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
    
    # Plot dens(x,y) for each species
    subplot_quantity_vs_xy(axes[6], simulation, "dens", specie=0, t_range=[tspecific,None], z=None, interpolation_step=interpolation_step,
        remove_zonal_modes=False, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
    subplot_quantity_vs_xy(axes[7], simulation, "dens", specie=0, t_range=[tspecific,None], z=0, interpolation_step=interpolation_step,
        remove_zonal_modes=False, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
    subplot_quantity_vs_xy(axes[8], simulation, "dens", specie=1, t_range=[tspecific,None], z=0, interpolation_step=interpolation_step,
        remove_zonal_modes=False, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
    try: subplot_quantity_vs_xy(axes[9], simulation, "dens", specie=2, t_range=[tspecific,None], z=0, interpolation_step=interpolation_step,
        remove_zonal_modes=False, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude) 
    except: pass 
                 
    # Show the figure 
    mpl.rcParams["savefig.directory"] = folder
    plt.show()
    return
 
################################################################################
#                                  DEBUG MODE                                  #
################################################################################
    
if __name__ == "__main__":
    import pathlib
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/zeff1.000001_m12.011_zz6.0_fprimz6.0_3kinspecies")
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/zeff1.000001_m12.011_zz6.0_fprimz6.0_kinions")
    plot_quantity_vs_xy_overview(folder)
    sys.exit()
     
     

