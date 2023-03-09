""" 

#===============================================================================
#            Get boolean map for each area i defined in <areas_map>            #
#===============================================================================

A map versus (variables1, variables2) can be split into different areas
labeled {1, 2, 3, ...} where each area will be interpolated seperaptly. For
linear maps each area is a specific mode which changes continuously, while the 
mode structure, growth rate, and frequency jump when entering a different 
area, or a different mode (e.g. ITG to TEM to ETG).

areas_map
---------
Define the <areas_map> through a dictionary:
    areas_map[variables1] = { 0 :  [variables2[0], variables2[3]], ... }  

For example, the original map is divided in the following areas:
    areas_map = [0, 0, 0, 1
                 0, 1, 1, 2
                 0, 1, 2, 2
                 1, 2, 2, 2]
                     
Therefore, the <areas_map> can be defined as:   
    variables1 = [0.0, 1.0, 2.0, 3.0]
    variables2 = [0.0, 1.0, 2.0, 3.0]
    areas_ids = [0, 1, 2] 
    areas_map =  {
        0.0  : {0 : [0.0, 2.0], 1 : [3.0, 3.0]},
        1.0  : {0 : [0.0, 0.0], 1 : [1.0, 2.0], 2 : [3.0, 3.0]},
        2.0  : {0 : [0.0, 0.0], 1 : [1.0, 1.0], 2 : [2.0, 3.0]},
        3.0  : {1 : [0.0, 0.0], 2 : [1.0, 3.0]}}   

"""

import numpy as np

#===============================================================================
#            Get boolean map for each area i defined in <areas_map>            #
#===============================================================================
    
def get_booleanMapsPerAreaFromAreaMap(areas_map, variables1=None, variables2=None, flip_columns_rows=False):
    
    # If we only want the areas_map, we don't have the variables1 and variables2
    if variables1==None: return []
    if areas_map==None: return None
    
    # Check whether we processed each variable
    variable2_keys = [float(i) for i in areas_map.keys()]
    for variable2 in variables2:
        if variable2 not in variable2_keys:
            print("ABORT: variable2 = "+str(variable2)+" was missing in the areas_map map.")
            import sys; sys.exit()
            
    # Initiate
    ivariables2 = range(len(variables2))
    ivariables1 = range(len(variables1)) 
    boolean_maps = []
    area_ids = []
    
    # Get the mode ids {0, 1, 2, ...}
    for variable2 in areas_map.keys(): area_ids += [ int(i) for i in areas_map[variable2].keys()] 
    area_ids = list(set(area_ids))

    # For each <area_id>, create a map of 0s and 1s where the 1s belong to the mode
    for area_id in area_ids: 
        boolean_map = np.zeros((len(variables1),len(variables2))) 
        for variable2 in variables2:
            if area_id in areas_map[variable2].keys():
                ivariable2 = [i for i in ivariables2 if variables2[i]==variable2] 
                ivariable1 = [i for i in ivariables1 if variables1[i]>=areas_map[variable2][area_id][0] and variables1[i]<=areas_map[variable2][area_id][1]]
                boolean_map[ivariable1, ivariable2] = 1 
            if -area_id in areas_map[variable2].keys():
                ivariable2 = [i for i in ivariables2 if variables2[i]==variable2] 
                ivariable1 = [i for i in ivariables1 if variables1[i]>=areas_map[variable2][-area_id][0] and variables1[i]<=areas_map[variable2][-area_id][1]]
                boolean_map[ivariable1, ivariable2] = 1   
        boolean_maps.append(boolean_map.astype(bool))
        
    # If we accidently defined the maps as a function of variables2 then transpose them
    if flip_columns_rows:
        for i in range(len(boolean_maps)):
            boolean_maps[i] = boolean_maps[i].T
            
    # Return the boolean maps
    return boolean_maps
