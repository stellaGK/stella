
import numpy as np
from stellapy.utils.files import initiate_nesteddict

def get_differenceInDictionaries(dict1, dict2):
    ''' Returns the difference in keys between two dictionaries. '''
    
    # Initiate the different dictionary
    dict_difference = initiate_nesteddict() 
    
    # The following two vmec files are the same (linearly)
    ignore = ["wout_tj20.nc", "wout_tjii_edi.nc"] 
    
    # Get the difference of a two layered dictionary
    if isinstance(dict1[list(dict1.keys())[0]], dict):
        for knob, sub_dict1 in dict1.items():
            for key, value in sub_dict1.items(): 
                
                # Compare two numbers
                if type(value)==float or type(value)==int or isinstance(value, np.floating):
                    if value != np.nan and not np.isnan(value) and value != None and dict2[knob][key] != None:  
                        if dict2[knob][key] != value:  
                            if value == 0.0 or dict2[knob][key] == 0.0: change_percent = 2
                            else: change_percent = ((float(dict2[knob][key])-value)/value)*100 
                            if abs(change_percent)>1: 
                                if knob not in dict_difference:
                                    dict_difference[knob] = [key]
                                else:
                                    dict_difference[knob].append(key)
                                    
                # Compare two strings
                else:
                    if dict2[knob][key] != value:  
                        if dict2[knob][key] not in ignore or value not in ignore:
                            if knob not in dict_difference:
                                dict_difference[knob] = [key]
                            else:
                                dict_difference[knob].append(key)
    
    # Some keys present the same thing  
    if "zgrid_parameters" in dict_difference.keys():
        if "nzgrid" in dict_difference["zgrid_parameters"]:
            if "nzed" in dict_difference["zgrid_parameters"]:
                dict_difference["zgrid_parameters"].remove("nzgrid")
    if "kt_grids_box_parameters" in dict_difference.keys():
        if "nx" in dict_difference["kt_grids_box_parameters"]:
            if "nakx" in dict_difference["kt_grids_box_parameters"]:
                dict_difference["kt_grids_box_parameters"].remove("nakx")
    if "kt_grids_box_parameters" in dict_difference.keys():
        if "ny" in dict_difference["kt_grids_box_parameters"]:
            if "naky" in dict_difference["kt_grids_box_parameters"]:
                dict_difference["kt_grids_box_parameters"].remove("naky")
    if "kt_grids_box_parameters" in dict_difference.keys():
        if "kx max" in dict_difference["kt_grids_box_parameters"]:
            if "nx" in dict_difference["kt_grids_box_parameters"]:
                dict_difference["kt_grids_box_parameters"].remove("nx")  
    if "kt_grids_box_parameters" in dict_difference.keys():
        if "ky max" in dict_difference["kt_grids_box_parameters"]:
            if "ny" in dict_difference["kt_grids_box_parameters"]:
                dict_difference["kt_grids_box_parameters"].remove("ny")  
    if "kt_grids_box_parameters" in dict_difference.keys():
        if "y0" in dict_difference["kt_grids_box_parameters"]:
            if "dkx" in dict_difference["kt_grids_box_parameters"]:
                    dict_difference["kt_grids_box_parameters"].remove("dkx")
            if "dky" in dict_difference["kt_grids_box_parameters"]:
                dict_difference["kt_grids_box_parameters"].remove("dky")
            if "Lx" in dict_difference["kt_grids_box_parameters"]:
                dict_difference["kt_grids_box_parameters"].remove("Lx")
            if "Ly" in dict_difference["kt_grids_box_parameters"]:
                dict_difference["kt_grids_box_parameters"].remove("Ly")
    if "vpamu_grids_parameters" in dict_difference.keys():
        if "nvgrid" in dict_difference["vpamu_grids_parameters"]:
            if "dvpa" in dict_difference["vpamu_grids_parameters"]:
                dict_difference["vpamu_grids_parameters"].remove("dvpa")
        if "vpa_max" in dict_difference["vpamu_grids_parameters"]:
            if "dvpa" in dict_difference["vpamu_grids_parameters"]:
                dict_difference["vpamu_grids_parameters"].remove("dvpa")
    if "vpamu_grids_parameters" in dict_difference.keys():
        if "nmu" in dict_difference["vpamu_grids_parameters"]:
            if "dmu" in dict_difference["vpamu_grids_parameters"]:
                dict_difference["vpamu_grids_parameters"].remove("dmu")
        if "vperp_max" in dict_difference["vpamu_grids_parameters"]:
            if "dmu" in dict_difference["vpamu_grids_parameters"]:
                dict_difference["vpamu_grids_parameters"].remove("dmu")  
                
    # Count the number of differences
    numberOfDifferences = 0
    for knob in dict_difference.keys():
        numberOfDifferences += len(dict_difference[knob])
    return dict_difference, numberOfDifferences
