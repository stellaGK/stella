
import numpy as np

def update_limits(x_range, y_range, xdata, ydata):
    if len(xdata)==0: return x_range, y_range 
    if np.all(np.isnan(ydata)): return x_range, y_range 
    x_range[0] = np.nanmin([x_range[0], np.nanmin(xdata)])
    x_range[1] = np.nanmax([x_range[1], np.nanmax(xdata)])
    y_range[0] = np.nanmin([y_range[0], np.nanmin(ydata)])
    y_range[1] = np.nanmax([y_range[1], np.nanmax(ydata)])  
    return x_range, y_range