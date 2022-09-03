
# Modules
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator  

# Personal modules 
from stellapy.GUI.plot.utils.Plot import Plot 

################################################################################
#                         CREATE FIGURE FOR THE GUI                            #
################################################################################

def create_figure(ax, plot, research, quantity=None, 
        top=0.95, left=0.1, right=0.9, bottom=0.12):
    ''' Function to create and manipulate the appearance of a plot.  '''

    #===========================================================================
    #                             CREATE THE FIGURE                            #
    #===========================================================================
    
    # Create a plotting class if was not given
    if not plot:
        plot = Plot(quantity)
        plot.process_plottingVariables(research)  
    
    # Create a new figure if no existing axis was given
    if not ax: 
        fig = plt.figure(figsize=plot.fig_size)
        grid_specifications = gridspec.GridSpec(1, 1, figure=fig)
        grid_specifications.update(top=top, left=left, right=right, bottom=bottom)
        ax = plt.subplot(grid_specifications[0]) 
        fig.set_tight_layout(False)

    #===========================================================================
    #                          ADD LABELS AND A TITLE                          #
    #=========================================================================== 
    
    # Add labels to the figure 
    ax.set_xlabel(plot.x_label)
    ax.set_ylabel(plot.y_label)
 
    # Add a title 
    ax.set_title(plot.title, fontsize=plot.fontsize_title)

    #===========================================================================
    #                           CHANGE THE APPEARANCE                          #
    #===========================================================================
    
    # Change the axis range (useless since this needs to be adjusted after plotting)
    ax.set_xlim(plot.x_range)
    ax.set_ylim(plot.y_range)
    
    # Change the appearance of the axis
    ax.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
    ax.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)  
    
    # Make sure we have minor ticks
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax.xaxis.set_minor_locator(AutoMinorLocator(10)) 
    
    # Set logaritmic axis or assure we have a linear axis  
    try: ax.ticklabel_format(useOffset=False)
    except: pass
    try: ax.ticklabel_format(style='plain')
    except: pass 
    
    # Return the axis object
    return ax



    
    
    
    
    
    
