
import numpy as np
from stellapy.utils.files import initiate_nesteddict

def get_differenceInDictionaries(dict1, dict2, ignore_resolution=False):
    ''' Returns the difference in keys between two dictionaries. '''

    # Initiate the different dictionary
    dict_difference = initiate_nesteddict()
    
    # The following two vmec files are the same (linearly)
    ignore = ['wout_tj20.nc', 'wout_tjii_edi.nc']
    
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
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'z_grid', 'nzed', 'nzgrid')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'z_grid', 'nzed', 'nz')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'kxky_grid_box', 'nx', 'nakx')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'kxky_grid_box', 'ny', 'naky')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'kxky_grid_box', 'kx max', 'nx')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'kxky_grid_box', 'ky max', 'ny')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'kxky_grid_box', 'y0', 'dkx')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'kxky_grid_box', 'y0', 'dky')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'kxky_grid_box', 'y0', 'Lx')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'kxky_grid_box', 'y0', 'Ly')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'velocity_grids', 'nvgrid', 'dvpa')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'velocity_grids', 'vpa_max', 'dvpa')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'velocity_grids', 'nmu', 'dmu')
    dict_difference = remove_repetitionsFromDictionary(dict_difference, 'velocity_grids', 'vperp_max', 'dmu')
                
    # For a linear simulations (nakx, naky) do not matter
    if dict1['gyrokinetic_terms']['include_nonlinear']==False: 
        dict_difference = remove_keyKnobFromDictionary(dict_difference, 'kxky_grid_box', 'nakx')
        dict_difference = remove_keyKnobFromDictionary(dict_difference, 'kxky_grid_box', 'naky')
    
    # Do not take into account resolution
    if ignore_resolution==True:
        dict_difference = remove_keyKnobFromDictionary(dict_difference, 'time_step', 'delt')
                
    # Count the number of differences
    numberOfDifferences = 0
    for knob in dict_difference.keys():
        numberOfDifferences += len(dict_difference[knob])
    return dict_difference, numberOfDifferences

#-------------------------------------
def remove_repetitionsFromDictionary(dictionary, knob, key1, key2):
    if knob in dictionary.keys():
        if key1 in dictionary[knob]:
            if key2 in dictionary[knob]:
                dictionary[knob].remove(key2)
    return dictionary
            
#-------------------------------------
def remove_keyKnobFromDictionary(dictionary, knob, key):
    if knob in dictionary.keys():
        if key in dictionary[knob]:
            dictionary[knob].remove(key)
    return dictionary
