
import configparser
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator  

#===============================================================================
#                                 CONFIG.INI                                   #
#=============================================================================== 

def read_figsize():
    test_directory = get_automatic_tests_directory() 
    config = configparser.ConfigParser() 
    config.read(test_directory / 'config.ini')
    width = int(config['PLOTS']['width'])
    height = int(config['PLOTS']['height'])
    return (width, height)
    
def read_fontsize():
    test_directory = get_automatic_tests_directory() 
    config = configparser.ConfigParser() 
    config.read(test_directory / 'config.ini')
    font_size = int(config['PLOTS']['font_size'])
    return font_size
    
def read_latex():
    test_directory = get_automatic_tests_directory() 
    config = configparser.ConfigParser() 
    config.read(test_directory / 'config.ini')
    latex = str(config['PLOTS']['latex']) 
    if 'T' in latex.upper(): latex = True
    elif 'F' in latex.upper(): latex = False 
    return latex

def set_latex_font():
	if read_latex():
		params = {'text.usetex' : True, 'font.size' : font_size, 'font.family' : 'lmodern'}
		plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
		plt.rcParams['axes.unicode_minus'] = True
		plt.rcParams.update(params) 

#===============================================================================
#                                CREATE FIGURE                                 #
#===============================================================================

def create_figure(xlabel=None, ylabel=None, title=None, font_size=-1, figsize=(-1,-1),
    scientific_axis=False,  top=0.95, left=0.08, right=0.95, bottom=0.1):
    ''' Function to create and manipulate the appearance of a plot.  '''
    
    # Read config file
    if font_size<0:
    	font_size = read_fontsize()
    if figsize[0]<0:
    	figsize = read_figsize()

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
    if xlabel: ax.set_xlabel(xlabel, labelpad=10)
    if ylabel: ax.set_ylabel(ylabel, labelpad=20)
 
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
    
    
#------------------------------------
def increase_fontsize(ax, fontsize, cbar=None, label=None, labelpad=20, cbar_color='black'):
    
    # Increase the font size of the ticks
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize)
    ax.yaxis.get_offset_text().set_fontsize(fontsize)
    
    # Increase the font size of the title and labels
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(fontsize) 
        
    # Increase the font size of the color bar
    if cbar!=None:
        cbar.outline.set_edgecolor(cbar_color)
        cbar.ax.yaxis.set_tick_params(color=cbar_color)
        cbar.ax.tick_params(labelsize=fontsize, color=cbar_color)  
        if label: cbar.set_label(label, labelpad=labelpad, y=0.52, color=cbar_color, labelcolor=cbar_color)
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color=cbar_color)
        text = cbar.ax.yaxis.label
        font = matplotlib.font_manager.FontProperties(size=fontsize)
        text.set_font_properties(font)
        label = text.get_text()  
        cbar.ax.yaxis.get_offset_text().set_fontsize(fontsize)
    return
