 
#!/usr/bin/python3   
import sys, os 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec 
from stellapy.plot.nonlinear.quantity_vs_xy import subplot_quantity_vs_xy

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.plot.utils.style.create_figure import update_figure_style 
from stellapy.simulations.Research import create_research 

#===============================================================================
#                               Plot phi(x,y)                                #
#===============================================================================

def plot_quantity_vs_xy(folder, z_quantity="phi", z=None, interpolation_step=20,
        remove_zonal_modes=False, log=False, crange=None, ordersOfMagnitude=2): 
    
    # Create <simulations> based on the given <folder>
    research = create_research(folders=folder) 
  
    # Plot each simulation separately 
    for experiment in research.experiments:
        for simulation in experiment.simulations:  
     
            # Create a figure   
            fig = plt.figure(figsize=(18, 9)); axes = []
            grid_specifications = gridspec.GridSpec(2, 4)
            grid_specifications.update(top=0.95, left=0.05, right=0.95, bottom=0.08, wspace=0.4, hspace=0.4)
            for i in range(8): axes.append(plt.subplot(grid_specifications[i]))
            update_figure_style(fig, axes)
             
            # Plot phi(x,y) 
            subplot_quantity_vs_xy(axes[0], simulation, z_quantity, specie=0, t_range=[None,None], z=z, interpolation_step=interpolation_step,
                remove_zonal_modes=remove_zonal_modes, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
            subplot_quantity_vs_xy(axes[1], simulation, z_quantity, specie=0, t_range=[300,None], z=z, interpolation_step=interpolation_step,
                remove_zonal_modes=remove_zonal_modes, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
            subplot_quantity_vs_xy(axes[2], simulation, z_quantity, specie=0, t_range=[400,None], z=z, interpolation_step=interpolation_step,
                remove_zonal_modes=remove_zonal_modes, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
            subplot_quantity_vs_xy(axes[3], simulation, z_quantity, specie=0, t_range=[500,None], z=z, interpolation_step=interpolation_step,
                remove_zonal_modes=remove_zonal_modes, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)    
            subplot_quantity_vs_xy(axes[4], simulation, z_quantity, specie=0, t_range=[None,None], z=0, interpolation_step=interpolation_step,
                remove_zonal_modes=remove_zonal_modes, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
            subplot_quantity_vs_xy(axes[5], simulation, z_quantity, specie=0, t_range=[300,None], z=0, interpolation_step=interpolation_step,
                remove_zonal_modes=remove_zonal_modes, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
            subplot_quantity_vs_xy(axes[6], simulation, z_quantity, specie=0, t_range=[300,None], z=2, interpolation_step=interpolation_step,
                remove_zonal_modes=remove_zonal_modes, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude)  
            subplot_quantity_vs_xy(axes[7], simulation, z_quantity, specie=0, t_range=[300,None], z=3, interpolation_step=interpolation_step,
                remove_zonal_modes=remove_zonal_modes, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude) 
                    
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
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/zeff1.5_m12.011_zz6.0_fprimz6.0")
    plot_quantity_vs_xy(folder)
    sys.exit()
     
     

