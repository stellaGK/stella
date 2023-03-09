 
#!/usr/bin/python3   
import sys, os 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec  

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.plot.nonlinear.quantity_vs_xz import subplot_quantity_vs_xz
from stellapy.plot.nonlinear.quantity_vs_yz import subplot_quantity_vs_yz
from stellapy.plot.utils.style.create_figure import update_figure_style 
from stellapy.simulations.Research import create_research 

#===============================================================================
#                               Plot phi(x,y)                                #
#===============================================================================

def quantity_vs_xz_and_xy(folder, z_quantity="phi", interpolation_step=20,
        remove_zonal_modes=True, log=False, crange=None, ordersOfMagnitude=2): 
    
    # Create <simulations> based on the given <folder>
    research = create_research(folders=folder) 
  
    # Plot each simulation separately 
    for experiment in research.experiments:
        for simulation in experiment.simulations:  
     
            # Create a figure   
            fig = plt.figure(figsize=(18, 9)); axes = []
            grid_specifications = gridspec.GridSpec(2, 2)
            grid_specifications.update(top=0.95, left=0.05, right=0.95, bottom=0.08, wspace=0.25, hspace=0.4)
            for i in range(4): axes.append(plt.subplot(grid_specifications[i]))
            update_figure_style(fig, axes)
             
            # Plot phi(x,y) 
            subplot_quantity_vs_yz(axes[0], simulation, z_quantity, specie=0, x=None, t_range=[300, None], interpolation_step=interpolation_step,
                remove_zonal_modes=remove_zonal_modes, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude) 
            subplot_quantity_vs_yz(axes[1], simulation, z_quantity, specie=0, x="max", t_range=[300, None], interpolation_step=interpolation_step,
                remove_zonal_modes=remove_zonal_modes, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude) 
            subplot_quantity_vs_xz(axes[2], simulation, z_quantity, specie=0, y=None, t_range=[300, None], interpolation_step=interpolation_step,
                remove_zonal_modes=remove_zonal_modes, log=log, crange=crange, ordersOfMagnitude=ordersOfMagnitude) 
            subplot_quantity_vs_xz(axes[3], simulation, z_quantity, specie=0, y="max", t_range=[300, None], interpolation_step=interpolation_step,
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
    quantity_vs_xz_and_xy(folder)
    sys.exit()
     
     

