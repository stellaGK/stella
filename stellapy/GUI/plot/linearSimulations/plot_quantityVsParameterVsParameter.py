
# Load modules 
import copy
import numpy as np
import matplotlib.pyplot as plt  
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
from scipy.interpolate import interp2d
from scipy import interpolate
from matplotlib.ticker import MultipleLocator

# Personal modules
from stellapy.utils.decorators.verbose import noverbose    

# Plotting modules
from stellapy.GUI.plot.utils.Plot import Plot
from stellapy.GUI.plot.utils import Axis 
from stellapy.GUI.plot.utils import create_figure
from stellapy.GUI.utils import DisplayProgressGUI 
 
#===============================================================================
#           PLOT GAMMA(PARAMETER,PARAMETER) OF THE MOST UNSTABLE MODE          #
#===============================================================================

@noverbose
def plot_quantityVsParameterVsParameter(\
            # Simulations
                research=None,\
                experiment_id="All experiments",\
                simulation_id="All simulations",\
            # Parameters
                knob1="vmec_parameters",\
                key1="rho",\
                knob2="vmec_parameters",\
                key2="rho",\
            # Data
                z_quantity='gamma',\
                x_range=None,\
                y_range=None,\
            # Modes
                modes_id="all",\
                kx_range=-0.0, \
                ky_range=[0,100], \
                dynamickyrange=False,\
                c_range=None,\
            # Labels
                x_label=None,\
                y_label=None,\
                title=None,\
            # Figure 
                ax=None, \
                Progress=None,\
                show_figure = False,\
            # Toggles
                step=4,\
                log=False,\
                interpolate=False): 
 
    # Save the plotting details to a <plot> object
    plot = Plot()   
    plot.update_labels(title, x_label, y_label ) 
    plot.update_toggles(interpolate=interpolate, step=step)
    plot.update_figure(Progress, ax, show_figure)
    plot.update_simulations(experiment_id, simulation_id, species=[0]) 
    plot.update_xyzData(key1, key2, None, x_range, y_range) 
    plot.update_modes(modes_id, kx_range, ky_range)
    plot.update_parameter(knob1=knob1, key1=key1, knob2=knob2, key2=key2) 
    plot.process_plottingVariables(research)
    
    # Update the progress bar of the GUI   
    displayProgressGUI = DisplayProgressGUI(research, "Get linear map")   
    
    # Create the figure and a <legend> and <axis> object
    ax = create_figure(ax, plot, research) 
    axis = Axis(ax, plot, overshoot_y=1); axis.set_limits(xbot=0, ybot=0)  
        
    # Get the parameters
    parameters1 = get_parameters(research, key1, knob1)
    parameters2 = get_parameters(research, key2, knob2) 
        
    # Get the data
    gamma, omega, ky = get_linearData(research, plot, parameters1, parameters2, knob1, key1, knob2, key2)
 
    # Update the progress bar of the GUI 
    if Progress: Progress.move(50,"Plotting linear map.")
          
    # Get the z-data 
    if z_quantity=="gamma": y_data = gamma
    if z_quantity=="omega": y_data = omega
    if z_quantity=="ky":    y_data = ky  
    
    # Recall where the gamma data is zero
    filter_gamma_zero = np.abs(y_data) < 1.E-25 
             
    # Interpolate the data
    if interpolate:   
        y_data = remove_nansFromGrid(y_data, parameters1, parameters2)
        parameters1, parameters2, y_data, filter_gamma_zero = interpolate_data(step, y_data, parameters1, parameters2, filter_gamma_zero)  
        parameters1_temp = [p for p in parameters1]
        parameters2_temp = [p for p in parameters2]
        
    # Add an extra parameter because the color is shown from x to x+1, add offset the labels
    if not interpolate:
        parameters1 = parameters1 + [parameters1[-1] + (parameters1[1]-parameters1[0])]
        parameters2 = parameters2 + [parameters2[-1] + (parameters2[1]-parameters2[0])]  
        parameters1_temp = [p -0.5 for p in parameters1]
        parameters2_temp = [p -0.5 for p in parameters2]
     
    # Replace zeros by nan
    y_data[filter_gamma_zero] = np.NaN
    
    # Add black to the start of the jet colormap, to indiciate unstable modes:
    cdict_jet, cdict_seismic = get_colorMaps()
     
    # Define the color map and range 
    #cmap = copy.copy(plt.get_cmap("jet")) if z_quantity!="omega" else copy.copy(plt.get_cmap("seismic"))
    cmap = cdict_jet if z_quantity!="omega" else cdict_seismic
    vmax = np.nanmax(np.nanmax(y_data, axis=0)) if c_range==None else c_range[1]
    vmin = np.nanmin(np.nanmin(y_data, axis=0)) if c_range==None else c_range[0]
    vmin = 0 if (vmin>0 and not log) else vmin
    vmin = -0.1 if vmin==0 and  (not log and z_quantity=="omega") else vmin
    norm = LogNorm(vmin=vmin, vmax=vmax) if log else None 
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax) if (not log and z_quantity=="omega") else None
     
    # Plot y_avrg(kx,ky)    
    img = ax.pcolormesh(parameters1_temp, parameters2_temp, y_data.T, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)  
    cbar = plt.colorbar(img, ax=ax, pad=0.02)
    if c_range: img.set_clim(c_range[0],c_range[1])  
     
    # Set a label for the colorbar
    if z_quantity=="gamma": cbar.set_label('$\\gamma a/v_{\\mathrm{th},i}$')
    if z_quantity=="omega": cbar.set_label('$\\omega a/v_{\\mathrm{th},i}$')
    if z_quantity=="ky": cbar.set_label('$k_{y}\\rho_i$')
          
    # Change appearance plot 
    ax.autoscale()
    
    # Keep track of the axis limits
    ax.xaxis.set_major_locator(MultipleLocator(1)) 
    ax.yaxis.set_major_locator(MultipleLocator(1)) 
    axis.update_axisLimits(parameters1_temp, parameters2_temp)
    axis.rescale_axis()   
     
    # Set the ticks
    if not interpolate: 
        if 0.0 in parameters1 and 8.0 in parameters1:
            parameters1 = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5]
            parameters2 = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5]
        plt.xticks(parameters1[:-1])
        plt.yticks(parameters2[:-1])        
        
    # Show the figure
    if show_figure: plt.show()
    return cbar
 
  
#===============================================================================
#                              CALCULATE THE DATA                              #
#=============================================================================== 

def get_parameters(research, key, knob): 
    
    # Intitiate the parameters
    parameters = [] 
     
    # Get the requested parameter
    for experiment in research.plotted_experiments: 
        for simulation in experiment.plotted_simulations:
            if key!="explicit_option" and key!="nfield_periods":
                value = float(simulation.input.inputParameters[knob][key])
            if key=="nfield_periods": 
                if simulation.vmec:     value = simulation.input.nfield_periods
                if simulation.miller:   value = simulation.input.nperiod 
            if key=="explicit_option": 
                if (simulation.input.inputParameters[knob][key]=="rk2"): value = 2 
                if (simulation.input.inputParameters[knob][key]=="rk3"): value = 3
            parameters.append(value)
     
    # Sort the simulations by this parameter
    parameters = sorted(list(set(parameters)))
    return parameters

#------------------
def get_linearData(research, parameters1, parameters2, knob1, key1, knob2, key2):
    
    # Initiate the data
    gamma = np.empty((len(parameters1), len(parameters2)))*[np.NaN]
    omega = np.empty((len(parameters1), len(parameters2)))*[np.NaN]
    ky    = np.empty((len(parameters1), len(parameters2)))*[np.NaN] 

    # Add the (0,0) point since we can't simulate it
    if 0 in parameters1 and 0 in parameters2: 
        gamma[0,0] = 0
        omega[0,0] = 0
        ky[0,0] = 0 
            
    # Iterate over the simulations
    for experiment in research.plotted_experiments: 
        for simulation in experiment.plotted_simulations:
         
            # Update the progress bar of the GUI  
            displayProgressGUI.move_progressBar() 
            
            # Get the value of the parameters for this <simulation>
            parameter1 = simulation.inputParameters[plot.knob1][plot.key1]
            parameter2 = simulation.inputParameters[plot.knob2][plot.key2]
            
            # Get the index of the parameters  
            a = parameters1.index(parameter1)
            b = parameters2.index(parameter2)
            
            # Get the linear data
            iky = [mode.ky if mode.lineardata.unstable else 0 for mode in simulation.plotted_modes]
            iomega = [mode.lineardata.omega_avg if mode.lineardata.unstable else 0 for mode in simulation.plotted_modes]
            igamma = [mode.lineardata.gamma_avg if mode.lineardata.unstable else 0 for mode in simulation.plotted_modes]
              
            # Find the most unstable mode
            index = np.nanargmax(igamma) 
            gamma[a,b] = igamma[index]
            omega[a,b] = iomega[index]
            ky[a,b] = iky[index]
            
    return gamma, omega, ky  

#----------------
def remove_nansFromGrid(y_data, parameters1, parameters2): 
    # Nans represent missing simulation, use an interpolation to fill them
    if len(np.argwhere(np.isnan(y_data)))>0:
        y_data = np.ma.masked_invalid(y_data) 
        x, y = np.meshgrid(parameters1, parameters2)
        xx = x[~y_data.mask]
        yy = y[~y_data.mask]
        zz = y_data[~y_data.mask]
        y_data = interpolate.griddata((xx, yy), zz.ravel(), (x, y), method='cubic')  
    return y_data

#----------------
def interpolate_data(step, y_data, parameters1, parameters2, filter_gamma_zero): 
    # Remember the original arrays 
    parameters1_temp = [p for p in parameters1]
    parameters2_temp = [p for p in parameters2]  
    # Interpolate   
    function = interp2d(parameters1, parameters2, y_data, kind='linear')
    xnew = np.linspace(0, parameters1[-1], int(len(parameters1)+1)*step)
    ynew = np.linspace(0, parameters2[-1], int(len(parameters2)+1)*step)
    y_data = function(xnew,ynew)
    parameters1, parameters2 = np.meshgrid(xnew, ynew)
    parameters1 = parameters1[0,:]
    parameters2 = parameters2[:,0]
    # Get where we had zeros (stable modes)
    filter_gamma_zero_dummy = np.zeros(np.shape(y_data)) 
    for i in range(len(parameters1_temp)-1):
        for j in range(len(parameters1_temp)-1):
            indexes_k = [k for k in range(len(parameters1)) if parameters1[k]>=parameters1_temp[i] and parameters1[k]<parameters1_temp[i+1]]
            indexes_l = [k for k in range(len(parameters2)) if parameters2[k]>=parameters2_temp[j] and parameters2[k]<parameters2_temp[j+1]]
            for k in indexes_k:
                for l in indexes_l:
                    filter_gamma_zero_dummy[k,l] = filter_gamma_zero[i,j]
    filter_gamma_zero = filter_gamma_zero_dummy.astype(np.bool)
    return parameters1, parameters2, y_data, filter_gamma_zero
 
#-----------------------
def get_colorMaps():
        from matplotlib.colors import LinearSegmentedColormap
         
        # Jet map where the zero is set to black for stable modes
        cdict_jet = {'red': 
                            ((0., 0, 0),     # Black
                            (0.01, 0, 0),    # Black
                            (0.11, 0, 0),    # Blue
                            (0.66, 1, 1),    # Yellow
                            (0.89, 1, 1),    # Orange
                            (1, 0.5, 0.5)),  # Dark red
                    'green': 
                            ((0., 0, 0),    # Black
                            (0.01, 0, 0),    # Black
                            (0.11, 0, 0),    # Blue
                            (0.375, 1, 1),   # Greenish
                            (0.64, 1, 1),    # Greenish
                            (0.91, 0, 0),    # Orange
                            (1, 0, 0)),      # Dark red
                    'blue': 
                            ((0., 0, 0),     # Black
                            (0.01, 0, 1),    # Black
                            (0.11, 1, 1),    # Blue
                            (0.34, 1, 1),    # Light blue
                            (0.65, 0, 0),    # Yellow
                            (1, 0, 0))}      # Dark red
         
        cdict_jet = LinearSegmentedColormap('my_colormap1',cdict_jet,256)
        cdict_jet = copy.copy(plt.get_cmap("jet"))
         
        # Seismic map with red for positive and blue for negative values, 
        # Where zero is put to black for stable modes
        cdict_seismic = {'red': 
                            [[0.  , 0.  , 0.  ],
                            [0.25,  0.  , 0.  ],
                            [0.495, 1.  , 0.  ],     # White
                            [0.5,   0.  , 0.  ],     # Black
                            [0.505, 0.  , 1.  ],     # White
                            [0.75,  1.  , 1.  ],
                            [1.  ,  0.5 , 0.5 ]], 
                        'green': 
                            [[0.  , 0.  , 0.  ],
                            [0.25, 0.  , 0.  ],
                            [0.495, 1.  , 0.  ],     # White
                            [0.5,   0.  , 0.  ],     # Black
                            [0.505, 0.  , 1.  ],     # White
                            [0.75, 0.  , 0.  ],
                            [1.  , 0.  , 0.  ]], 
                        'blue': 
                            [[0. , 0.3 , 0.3 ],     # Darkbue
                            [0.25, 1.  , 1.  ],     # Blue
                            [0.495, 1.  , 0.  ],     # White
                            [0.5,   0.  , 0.  ],     # Black
                            [0.505, 0.  , 1.  ],     # White
                            [0.75, 0.  , 0.  ],     # Red
                            [1.  , 0.  , 0.  ]]}    # Darkred
         
        cdict_seismic = LinearSegmentedColormap('my_colormap2',cdict_seismic,256)
        cdict_seismic = copy.copy(plt.get_cmap("seismic")) 
         
        # Put nan values in purple
        cdict_jet.set_bad(color='black')
        cdict_seismic.set_bad(color='black')
        return cdict_jet, cdict_seismic

