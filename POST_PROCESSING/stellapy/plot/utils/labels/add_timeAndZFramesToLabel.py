
import numpy as np

def add_timeFrameToLabel(label, tstart, tend): 
    
    # Round the number to the closest multiple of 10 
    tstart = str(int(np.round(tstart/10,0)*10))
    tend = str(int(np.round(tend/10,0)*10)) if tend!=None else None
    
    # Return the string of the time frame
    if tend!=None: label = "$\\langle$"+label+"$\\rangle_{t=["+tstart+", "+tend+"]}$" 
    if tend==None: label =  label+" at $ t="+tstart+"$" 
    return label
    
def add_timeAndZFramesToLabel(label, tstart, tend, z): 
    
    # Round the number to the closest multiple of 10 
    tstart = str(int(round(tstart/10,0)*10))
    tend = str(int(np.round(tend/10,0)*10)) if tend!=None else None
    
    # Round z to one digit
    z = str(round(z,1)) if z!=None else None
    
    # Return the string of the time frame
    if tend!=None:
        if z==None: label = "$\\langle$"+label+"$\\rangle_{z,t=["+tstart+", "+tend+"]}$" 
        if z!=None: label = "$\\langle$"+label+"$\\rangle_{t=["+tstart+", "+tend+"]} $ at $ z = "+z+"$" 
    if tend==None:
        if z==None: return "$\\langle$"+label+"$\\rangle_{z} $ at $ t="+tstart+"$"
        if z!=None: return label+" at $t="+tstart+"; z="+z+"$"
    return label
    
    
    