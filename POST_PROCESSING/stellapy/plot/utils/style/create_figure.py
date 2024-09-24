
# Modules
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator  
from stellapy.plot.utils.style.increase_fontsize import increase_fontsize

#===============================================================================
#                                CREATE FIGURE                                 #
#===============================================================================

def create_figure(xlabel=None, ylabel=None, title=None, font_size=20, figsize=(18,9),
    scientific_axis=False,  top=0.95, left=0.08, right=0.95, bottom=0.1):
    ''' Function to create and manipulate the appearance of a plot.  '''
    
    #===========================================================================
    #                             PLOTTING PARAMETERS                          #
    #===========================================================================
    params = {'text.usetex' : True, 'font.size' : font_size, 'font.family' : 'lmodern'}
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    plt.rcParams['axes.unicode_minus'] = True
    plt.rcParams.update(params) 

    #===========================================================================
    #                             CREATE THE FIGURE                            #
    #===========================================================================
    
    # Create a new figure if no existing axis was given  
    fig = plt.figure(figsize=figsize)
    fig.grid_specifications = gridspec.GridSpec(1, 1, figure=fig)
    fig.grid_specifications.update(top=top, left=left, right=right, bottom=bottom)
    ax = plt.subplot(fig.grid_specifications[0]) 
    fig.set_tight_layout(False)
    ax.fig = fig; fig.ax = ax

    #===========================================================================
    #                          ADD LABELS AND A TITLE                          #
    #=========================================================================== 
    
    # Add labels to the figure 
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel, labelpad=10)
 
    # Add a title 
    if title: ax.set_title(title, fontsize=font_size)

    #===========================================================================
    #                           CHANGE THE APPEARANCE                          #
    #===========================================================================
    
    # Change the appearance of the axis
    ax.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
    ax.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)  
    
    # Make sure we have minor ticks
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax.xaxis.set_minor_locator(AutoMinorLocator(10)) 
    
    # Prevent scientific notation or use it 
    if scientific_axis: ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    if not scientific_axis: ax.ticklabel_format(style='plain', useOffset=False) 
    
    # Change the font size
    increase_fontsize(ax, font_size)
    
    # Return the axis object
    return ax

#------------------------------------
def update_figure_style(fig=None, axes=[], font_size=20, scientific_axis=False, lw=None):
    
    #===========================================================================
    #                             PLOTTING PARAMETERS                          #
    #===========================================================================
    params = {'text.usetex' : True, 'font.size' : font_size, 'font.family' : 'lmodern'}
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    plt.rcParams['axes.unicode_minus'] = True
    plt.rcParams.update(params) 
    if lw: mpl.rcParams['axes.linewidth'] = lw
    if fig: fig.set_tight_layout(False)             
    
    #===========================================================================
    #                           CHANGE THE APPEARANCE                          #
    #===========================================================================
    
    # Iterate over the axes
    for ax in axes:
    
        # Change the appearance of the axis
        ax.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
        ax.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)  
        
        # Make sure we have minor ticks
        ax.yaxis.set_minor_locator(AutoMinorLocator(10))
        ax.xaxis.set_minor_locator(AutoMinorLocator(10)) 
        
        # Prevent scientific notation or use it 
        if scientific_axis: ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        if not scientific_axis: ax.ticklabel_format(style='plain', useOffset=False) 
        
        # Change the font size
        increase_fontsize(ax, font_size) 
        
    return
    
    
    
#------------------------------------
def update_axis_style(ax):
    ax.grid(color='lightgray', linestyle='-', linewidth=0.2)
    ax.xaxis.set_tick_params(which='major', width=0.5, length=2, pad=2)
    ax.yaxis.set_tick_params(which='major', width=0.5, length=2, pad=2)
    ax.xaxis.set_tick_params(which='minor', width=0.2, length=1, pad=2)
    ax.yaxis.set_tick_params(which='minor', width=0.2, length=1, pad=2)
    return
    
    
    
