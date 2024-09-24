
import numpy as np

def get_timeFrame(simulation, normalize=False, y_quantity=None): 
    
    # Read the time range defined for this simulation
    t_start = simulation.time.tstart
    t_stop  = simulation.time.tend
    
    # Find the last time value of the time axis
    vec_time = simulation.fluxes.qflux_vs_ts.t
    t_last = vec_time[len(vec_time[~np.isnan(vec_time)]) - 1]
    
    # The time range ends at t_stop or the last time value
    t_stop = np.nanmin([t_last, t_stop])  

    # T_start can not be bigger than t_stop
    if t_start >= t_stop:
        t_start = t_stop-10

    # We can rescale the axis, where t=1 corresponds to the end of the linear growth of the heat flux
    if normalize and y_quantity in ['qflux', 'pflux', 'vflux']:
        t_start = t_start/simulation.get_peakTime()
        t_stop  = t_stop/simulation.get_peakTime()
             
    # Return the experiment and simulation label
    return t_start, t_stop
 
