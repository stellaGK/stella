from stellapy.plot.utils.devices.recognize_device import recognize_device


#----------------------
def get_colorsPerDevice(device): 
    """Give each device a fixed color."""
    
    # Initiate the color
    color = None
    
    # Recognize the exact device from a related string 
    device = recognize_device(device)
    
    # Define the colors for each device
    if device=="TOK":   color = "black"  
    if device=="CBC":   color = "black"  
    if device=="NCSX":  color = "green"
    if device=="LHD":   color = "red"
    if device=="TJII":  color = "orange" 
    if device=="W7X":   color = "navy"
    if device=="ASDEX": color = "magenta"
    if device=="AUG":   color = "magenta" 
    return color

#----------------------
def get_markerPerDevice(device): 
    """Give each device a fixed color."""
    
    # Initiate the color
    marker = None
    
    # Recognize the exact device from a related string 
    device = recognize_device(device)
    
    # Define the colors for each device
    if device=="TOK":   marker = "x"  
    if device=="CBC":   marker = "x"  
    if device=="NCSX":  marker = "X"
    if device=="LHD":   marker = "s"
    if device=="TJII":  marker = "^" 
    if device=="W7X":   marker = "o"
    if device=="ASDEX": marker = "+"
    if device=="AUG":   marker = "+" 
    return marker


#------------------------------
def replaceWithCustomLineColors(experiments, line_labels, line_colors):
    
    # Initiate the devices 
    devices = []  
    
    # First count whether we have multiple devices
    for experiment in experiments:  
        device = recognize_device(line_labels[experiment])
        if device!=None: devices.append(device)
        
    # If we have multiple devices, change the line color
    if len(list(set(devices)))>0:
        for experiment in experiments:  
            device = recognize_device(line_labels[experiment])
            if device!=None: line_colors[experiment] = get_colorsPerDevice(device)
            
    # Return the updates line_colors
    return line_colors
