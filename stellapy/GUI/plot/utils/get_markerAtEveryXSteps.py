
import numpy as np 

def get_markerAtEveryXSteps(x, y, time_of_peak, log=False, units="normalized", dx=None):
    ''' Instead of displaying markers at every point, only display 
    a marker after a fixed dx interval of the vector x.'''
    
    # Define the time step
    if dx==None:
        if log==False and units=="normalized":  dx = 25
        if log==False and units!="normalized":  dx = 0.05
        if log==True  and units=="normalized":  dx = 15
        if log==True  and units!="normalized":  dx = 0.05 
 
    # For normal scales, plot markers after the exponentional growth
    if log == False:    
        y = y[x >= time_of_peak]
        x = x[x >= time_of_peak]
        
    # for log scales, plot markers on the exponentional growth
    if log == True:  
        y = y[x <= time_of_peak]
        x = x[x <= time_of_peak]
        
    # Plot a marker at t=time_step, 2*time_step, ...  
    try:
        marker_times = np.arange(x[0],x[-1]+dx, dx)
        marker_times = list(marker_times) + [2*marker_times[-1]] 
        temp_x= []; temp_y = []; index = 0
        for ix in x:
            if ix >= marker_times[index]:
                temp_x.append(ix)
                temp_y.append(y[list(x).index(ix)])
                index = index + 1
    except:
        return x,y
 
    # Return the new (x,y) vectors with values every time_step in x
    return np.array(temp_x), np.array(temp_y)