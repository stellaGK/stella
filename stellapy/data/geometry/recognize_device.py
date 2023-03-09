
def recognize_device(name): 
    "Recognize the device (fusion reactor) from a string. "
    
    # Initatiate the device
    device = None
    name = str(name)
    
    # Recognize the device from a string (general)
    if "ncsx" in name.lower(): device = "NCSX"
    if "lhd" in  name.lower(): device = "LHD"
    if "w7x" in  name.lower(): device = "W7X"
    if "asd" in  name.lower(): device = "AUG"
    if "aug" in  name.lower(): device = "AUG"
    if "tok" in name.lower():  device = "CBC"
    if "cbc" in name.lower():  device = "CBC"
    if "tjii" in name.lower(): device = "TJII"
    if "tj2"  in name.lower(): device = "TJII"
    if name=='wout*.nc':       device = "CBC"
    
    # Recognize the device from a string (specific)
    if name=="W7X003":         device = "W7-X (EIM)"
    if name=="W7X027":         device = "W7-X (KJM)"
    if name=="LHDIS":          device = "LHD (IS)"
    if name=="LHD375":         device = "LHD (STD)" 
    
    # Return the device name 
    return device

#-----------------------------
def get_device_name(device):
    name = device 
    if device=="W7X": name = "W7-X" 
    if device=="TJII": name = "TJ-II" 
    return name
    