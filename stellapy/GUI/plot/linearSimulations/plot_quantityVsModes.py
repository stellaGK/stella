
# Load the modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.interpolate import interp1d

# Personal modules
from stellapy.utils.decorators.verbose import noverbose
from stellapy.GUI.plot.utils.get_dataForPlots import get_dataForPlots

# Simulation modules
from stellapy.simulations.utils.get_modes import get_modes
from stellapy.simulations.utils.get_simulations import get_simulations
from stellapy.simulations.utils.get_experiments import get_experiments

# Plotting modules
from stellapy.GUI.plot.utils.Plot import Plot
from stellapy.GUI.plot.utils import Axis, Legend
from stellapy.GUI.plot.utils import create_figure
from stellapy.GUI.utils import DisplayProgressGUI
from stellapy.plot.utils.style import get_styleForLinesAndMarkers 

#===============================================================================
#                    PLOT QUANTITY(K) FOR LINEAR SIMULATIONS                   #
#===============================================================================

@noverbose
def plot_quantityVsModes(\
            # Specify which simulations to plot
                research=None,
                experiment_id="All experiments",
                simulation_id="All simulations",
            # Data on the (x,y) axis  
                x_quantity='ky',
                y_quantity='gamma_avg',
                units="normalized", 
                x_range=None,
                y_range=None,
            # Modes
                modes_id="all",\
                kx_range=[0,0],
                ky_range=[0,100],  
            # Labels
                x_label=None,
                y_label=None,
                title=None, 
            # Toggles 
                fontsize=18,
                show_error=False, 
                maxima=False, 
                biggestmaxima=False,
                interpolate=False,
            # Figure 
                ax=None,  
                mfc=None,
                color=None,
                Progress=None,
                removelegend=False,
                show_figure=False): 
    
    # Recognize x-quantity
    if kx_range[0]==kx_range[1]: x_quantity = 'ky'
    if ky_range[0]==ky_range[1]: x_quantity = 'kx'

    # Save the plotting details to a <plot> object
    plot = Plot()  
    plot.update_toggles(show_error=show_error, maxima=maxima, biggestmaxima=biggestmaxima, interpolate=interpolate)
    plot.update_simulations(experiment_id, simulation_id, species=[0]) 
    plot.update_xyData(x_quantity, y_quantity, x_range, y_range, units)
    plot.update_modes(modes_id, kx_range, ky_range)
    plot.update_figure(Progress, ax, show_figure)
    plot.update_labels(title, x_label, y_label)
    plot.update_legend(fontsize=fontsize, removelegend=removelegend)
    plot.process_plottingVariables(research)  
    
    # Update the progress bar of the GUI 
    displayProgressGUI = DisplayProgressGUI(research, "Plot "+plot.ydim+" versus "+plot.xdim, plot.xdim, plot)  
    
    # Create the figure, the legend and the axis objects
    ax = create_figure(ax, plot, research) 
    legend = Legend(ax, plot); axis = Axis(ax, plot)
    axis.set_limits(xbot=0, ybot_pos=0)
    
    # Iterate over the experiments, simulations and modes, which are plotted for each kx
    for experiment in get_experiments(research, experiment_id): 
        for simulation in get_simulations(experiment, simulation_id):  
            
            # Iterate over the kx-modes, update the plotted modes and plot the data
            if x_quantity=="ky":
                for kx in [ kx for kx in simulation.vec.kx if (kx >= kx_range[0] and kx <= kx_range[1])]:   
                    simulation.plotted_modes = get_modes(simulation, [kx,kx], ky_range, modes_id) 
                    plot_lineardata(ax, simulation, experiment, research, plot, legend, axis, displayProgressGUI, ky_range, color, mfc)
                    
            # Iterate over the ky-modes, update the plotted modes and plot the data
            if x_quantity=="kx":
                for ky in [ ky for ky in simulation.vec.ky if (ky >= ky_range[0] and ky <= ky_range[1])]:   
                    simulation.plotted_modes = get_modes(simulation, kx_range, [ky,ky], modes_id) 
                    plot_lineardata(ax, simulation, experiment, research, plot, legend, axis, displayProgressGUI, ky_range, color, mfc)
                                        
    # Add the legend and rescale the axis
    legend.add_legend()
    axis.rescale_axis()
    
    # Show the figure
    if show_figure: plt.show()
    return 

#===============================================================================
#                      PLOT QUANTITY(K) ALONG FIXED KX OR KY                   #
#===============================================================================

def plot_lineardata(ax, simulation, experiment, research, plot, legend, axis, displayProgressGUI, ky_range, color, mfc):        

    # Update the progress bar of the GUI
    displayProgressGUI.move_progressBar()
    
    # Calculate the linear data for these plotted modes and save them as e.g. simulation.lineardata.gamma_avg_vs_ky
    simulation.lineardata.get_linearDataPerSimulation(simulation.plotted_modes, plot) 
    
    # Get the lines and marker styles for the experiments and simulations
    style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation) 
    if "marker" not in style: style["marker"] = experiment.marker_style  
    if isinstance(color, (np.ndarray, np.generic)): style["color"] = color; style["mfc"] = color 
    if isinstance(mfc, (np.ndarray, np.generic)): style["mfc"] = mfc
    
    # Plot the (modes, quantity) data
    if not plot.show_error:
        vec_modes, vec_quantity = get_dataForPlots(plot, simulation, plot.quantity) 
        ax.plot(vec_modes, vec_quantity, **style)    

    # Plot the (modes, quantity, error) data
    if plot.show_error: 
        vec_modes, vec_quantity = get_dataForPlots(plot, simulation, plot.quantity); del style['marker']
        vec_modes, vec_error = get_dataForPlots(plot, simulation, plot.quantity+"_error")    
        ax.errorbar(vec_modes, vec_quantity, yerr=vec_error, fmt="o", capsize=2, **style)      

    # Keep track of the axis limits 
    axis.update_axisLimits(vec_modes, vec_quantity)    
    
    # Add the interpolation of (modes, quantity)
    if plot.interpolate: 
        vec_modes, vec_quantity = interpolateData(ax, vec_modes, vec_quantity, ky_range, plot.quantity, style['color'])
        color = style["mfc"] if "mfc" in style else style["color"] 
        ax.plot(vec_modes, vec_quantity, lw=2 ,linestyle='-', color=lighten_color(color,0.5)) 
    
    # Highlight the mode with the maximum growth rate with a star 
    if plot.maxima or plot.biggestmaxima:
        if "omega" in plot.quantity: 
            try:
                vec_ky = [getattr(mode, plot.xdim) for mode in simulation.plotted_modes]  
                vec_gamma = [mode.lineardata.gamma_avg if mode.lineardata.unstable else 0 for mode in simulation.plotted_modes]      
                if plot.interpolate: vec_ky, vec_gamma = interpolateData(ax, vec_ky, vec_gamma, ky_range, plot.quantity, style['color'])
                if plot.biggestmaxima: ky, _ = get_absoluteMaximumOfGraph(vec_ky, vec_gamma)
                elif plot.maxima: ky, _ = get_maximaOfGraph(vec_ky, vec_gamma, plot.interpolate) 
                vec_quantity = np.array(vec_quantity)[vec_modes==ky]; vec_modes = ky
            except: pass
        if "gamma" in plot.quantity: 
            try:
                if plot.biggestmaxima: vec_modes, vec_quantity = get_absoluteMaximumOfGraph(vec_modes, vec_quantity)
                elif plot.maxima: vec_modes, vec_quantity = get_maximaOfGraph(vec_modes, vec_quantity, plot.interpolate)
            except: pass
        ax.plot(vec_modes, vec_quantity, lw=0, marker='*',mec='black',mfc='gold',label="_nolegend_",ms=12,zorder=10)    
    return 

#===============================================================================
#                                   METHODS                                    #
#===============================================================================

def interpolateData(ax,x,y,xrange,quantity,color,step=0.5):
    
    # Can't interpolate a point
    if len(x)<=1: return x,y
    # Count the jumps in the spectrum
    xjumps = []
    # For omega interpolate each frequency branch separately
    if "omega" in quantity: 
        for i in range(len(y)-1):  
            if np.abs(y[i+1]-y[i])>step: 
                xjumps.append((x[i]+x[i+1])/2)
                ax.axvline(x=(x[i]+x[i+1])/2,color=color,label="_nolegend_") 
    # For gamma interpolate where gamma>=0 
    if "gamma" in quantity:
        print("Force interpolation to go below zero.")
        if len(y[y<=0])>0:
            x = x[y>0]; x = list(x); x.append(2.1)
            y = y[y>0]; y = list(y); y.append(-1)
#         for i in range(len(y)-1):  
#             if y[i]>0 and y[i+1]<=0: 
#                 xjumps.append((x[i]+x[i+1])/2)
#                 ax.axvline(x=(x[i]+x[i+1])/2,color=color,label="_nolegend_")  
    # Interpolate each branch seperatly
    if len(xjumps)>0:
        xnew = np.linspace(xrange[0], xrange[1], num=1000, endpoint=True)
        xjumps = [xrange[0]] + xjumps + [xrange[1]+1]; ynew = []
        for i in range(len(xjumps)-1): 
            xcut = x[(x>=xjumps[i]) & (x<xjumps[i+1])]
            ycut = y[(x>=xjumps[i]) & (x<xjumps[i+1])] 
            xnewcut = xnew[(xnew>=xjumps[i]) & (xnew<xjumps[i+1])] 
            if np.all(ycut==0): ynew += list([0]*len(xnewcut))
            elif len(xcut)==1: ynew += list([ycut]*len(xnewcut))
            else:
                if len(xcut)==2: xcut = [xcut[0]] + [(xcut[0]+xcut[1])/2] + [xcut[1]]
                if len(ycut)==2: ycut = [ycut[0]] + [(ycut[0]+ycut[1])/2] + [ycut[1]]
                if len(xcut)==3: xcut = [xcut[0]] + [(xcut[0]+xcut[1])/2] + [xcut[1],xcut[2]]
                if len(ycut)==3: ycut = [ycut[0]] + [(ycut[0]+ycut[1])/2] + [ycut[1],ycut[2]]
                function = interp1d(xcut, ycut, kind='cubic', fill_value="extrapolate")
                ynew += list(function(xnewcut))
    # Interpolate the entire data, and extrapolate if necessary
    if len(xjumps)==0:       
        xnew = np.linspace(xrange[0], xrange[1], num=1000, endpoint=True)
        function = interp1d(x, y, kind='cubic', fill_value="extrapolate")
        ynew = function(xnew)
    return xnew, ynew

def get_maximaOfGraph(x, y, interpolate):
    if len(x)<=1: return x,y
    if interpolate:
        peaks, _ = find_peaks(y) 
        return x[peaks], y[peaks]
    elif not interpolate: 
        x = list(x); 
        y = list(y)
        ymax = [iy for iy in y[1:-1] if (iy>y[y.index(iy)-1] and iy>y[y.index(iy)+1]) ] 
        xmax = [x[y.index(iy)] for iy in ymax] 
        return xmax, ymax
    
def get_absoluteMaximumOfGraph(x, y): 
    return x[np.nanargmax(y)], np.nanmax(y)

def lighten_color(color, amount=0.5): 
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

################################################################################
#                  USE THIS PLOTTING FUNCTION AS A MAIN SCRIPT                 #
################################################################################ 
    
if __name__ == "__main__":   
    
    import pathlib, time
    from stellapy.simulations.Research import create_research  
    start = time.time()
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI_PREVIOUS")    
    folder = pathlib.Path("/home/hanne/CIEMAT/RUNS/TEST_NEW_GUI_PREVIOUS/fprim2tprim6/")  
    research = create_research(folders=folder)#, knob1="species_parameters_1", key1="tprim", knob2="species_parameters_1", key2="fprim")
    research.print_research()  
    plot_quantityVsModes(research, show_error=True, maxima=True, interpolate=True)
    print("    ---> It took "+str(time.time() - start)+" seconds to plot a quantity versus modes.")
    plt.show()
    
