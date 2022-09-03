import os
import matplotlib.pyplot as plt 

def load_styleFigures():

    # Set the style
    plt.rcParams['lines.linewidth'] = 2 
    plt.rc('font', family='serif')
    plt.rc('font', size=20)
    plt.rcParams.update({'figure.autolayout': True})  
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

    # Set the grid
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.color'] = "lightgray"

    # Try to implement this
    #ax.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
    #ax.yaxis.grid(color='grey', linestyle='-', linewidth=0.3) 
    #ax.yaxis.set_minor_locator(AutoMinorLocator(10))
    #ax.xaxis.set_minor_locator(AutoMinorLocator(10)) 
    
    # Marconi does not support the fancy latex font
    divider = '\\' if (os.name == 'nt') else '/'
    if os.getcwd().split(divider)[1] == 'marconi_work':
        plt.rc('text', usetex=False)
    if os.getcwd().split(divider)[1] == 'home':
        plt.rc('text', usetex=True)

# Force it to load before opening a pyplot interactive loop.
load_styleFigures()
