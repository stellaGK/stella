
import numpy as np 

def get_inputFilesWhoOnlyDifferInResolution(input_files, input_parameters, ignore_resolution=True, include_knobs=[]): 

    # Nonlinear simulations can not be merged
    if input_parameters[0]["physics_flags"]["nonlinear"]==True or ignore_resolution==False:
        return input_files, input_parameters, input_files
        
    # Save the input parameters of the first <input_file> 
    saved_input_files = {1 : {"input_parameters" : input_parameters[0], "input_files" : [input_files[0]]}}
    
    # Knobs that can't differ since they would definitely represent a "new" simulation 
    knobs = ["physics_flags", "parameters", "species_knobs", "geo_knobs", "millergeo_parameters", "vmec_parameters"]
    knobs += ["species_parameters_"+str(specie) for specie in np.arange(1,5)]
    knobs += include_knobs
    
    # Identify the unique input files who differ beyond their resolution
    for i, input_file in enumerate(input_files[1:]):
        
        # Read the input parameters of this <input_file>
        input_parameters_new = input_parameters[i+1] 
        
        # Initiate
        similar_input_file_found = False
        
        # Iterate over the saved input files  
        for j in saved_input_files.keys():
            
            # Get the input parameters of the saved <input_file> 
            input_parameters_saved = saved_input_files[j]["input_parameters"] 
            
            # Check if <input_parameters_new> differs from <input_parameters_saved> in the relevant parameters 
            for knob in knobs:
                if input_parameters_new[knob]!=input_parameters_saved[knob]: 
                    break                
        
            # If we managed to check all the knobs, the input_file is similar
            else:
                similar_input_file_found = True
                saved_input_files[j]["input_files"].append(input_file)
                break
                                
        # If the input parameters if <input_file> didn't match any of the saved files, then save it
        if similar_input_file_found==False: 
            numbers = list(saved_input_files.keys())
            k = np.max([int(number) for number in numbers]) + 1
            saved_input_files[k] = {"input_parameters" : input_parameters_new, "input_files" : [input_file]} 
 
    # Reconstruct <input_files> and <input_parameters>, and create <similar_input_files>
    numbers = list(saved_input_files.keys())
    input_files = [saved_input_files[i]["input_files"][0] for i in numbers]
    input_parameters = [saved_input_files[i]["input_parameters"] for i in numbers]
    similar_input_files = [saved_input_files[i]["input_files"] for i in numbers]
    return input_files, input_parameters, similar_input_files

