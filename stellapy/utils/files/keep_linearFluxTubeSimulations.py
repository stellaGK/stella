
#===============================================================================
#            Keep only input files from linear flux tube simulations            
#=============================================================================== 

def keep_linearFluxTubeSimulations(input_files):
    
    # List of input files to keep
    input_files_to_keep = []
    
    # Iterate over the input files
    for input_file in input_files: 

        # Default values in stella 
        full_flux_surface = False
        nonlinear = False
        
        # Open each input files
        with open(input_file, 'r') as input_data:
            
            # Read the text if the input file
            input_text = input_data.read()
            
            # Read what <nonlinear> is set to
            if "nonlinear" in input_text:
                boolean = input_text.split("nonlinear")[1].split('=')[1].split('\n')[0].replace(' ', '')
                if boolean==".false." or boolean=="F": nonlinear = False
                if boolean==".true." or boolean=="T":  nonlinear = True
            
            # Read what <full_flux_surface> is set to
            if "full_flux_surface" in input_text:
                boolean = input_text.split("full_flux_surface")[1].split('=')[1].split('\n')[0].replace(' ', '')
                if boolean==".false." or boolean=="F": full_flux_surface = False
                if boolean==".true." or boolean=="T":  full_flux_surface = True
                
        # If the simulation is linear and not full flux surface, then keep it
        if (not nonlinear) and (not full_flux_surface): input_files_to_keep.append(input_file)
                
    # Return the input files which are kept
    return input_files_to_keep
