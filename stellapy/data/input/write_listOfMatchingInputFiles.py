"""

#===============================================================================
#                    WRITE LIST OF MATCHING INPUT FILES
#===============================================================================
 
For linear flux-tube simulations, we often launch 1 mode per file, since we 
want to put these modes into a single <simulation>, make a list of input files
that correspond to the same <simulation>.

Hanne Thienpondt
18/09/2022
 
"""

#!/usr/bin/python3
import copy 
import os, sys 
import pathlib 

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.utils.files.keep_linearFluxTubeSimulations import keep_linearFluxTubeSimulations
from stellapy.utils.files.sort_listByNumbers import sort_listByNumbers
from stellapy.utils.files.get_filesInFolder import get_filesInFolder  
from stellapy.data.input.read_inputFile import read_inFile  

#===============================================================================
#                    WRITE LIST OF MATCHING INPUT FILES
#===============================================================================
 
def write_listOfMatchingInputFiles(folder):  
    
    # Get the input_files inside <folder>
    # Only keep linear flux tube simulations that have output files
    input_files = get_filesInFolder(folder, end=".in") 
    input_files = [i for i in input_files if os.path.isfile(i.with_suffix(".out"))]
    input_files = keep_linearFluxTubeSimulations(input_files) 
    
    # Collect the parent directories since we will write a list for each directory
    # For linear restarts, the old files are put in an "OLD" folder, include these in the list of OLD.parent
    parent_directories = list(set([i.parent if "OLD" not in i.parent.name else i.parent.parent for i in input_files]))
    
    # Only look at the input files in each parent folder
    for directory in parent_directories:    

        # Initialize a dictionary simulation{ID}{list,input_parameters}
        simulations = {}; input_files_in_dummy_files = []
        
        # Check the lists that already exist
        dummy_input_files = get_filesInFolder(directory, end="_dummy.in") 
        if dummy_input_files!=None:
            for dummy_input_file in dummy_input_files:
                
                # Read the input parameters
                input_parameters, _ = read_relevantInputParameters(dummy_input_file)
                
                # Read the input files in the dummy input file
                input_files_in_dummy_file = read_inputFilesInDummyInputFile(dummy_input_file)
                        
                # Create the simulation
                simulations[str(dummy_input_file)] = {}
                simulations[str(dummy_input_file)]["input_parameters"] = input_parameters
                simulations[str(dummy_input_file)]["input_files"] = input_files_in_dummy_file 
        
        # Get the input files in <directory>
        input_files_directory = [i for i in input_files if (i.parent==directory) or (i.parent.parent==directory and "OLD" in i.parent.name)]
            
        # Get the input files that are already listed
        for simulation_id in simulations.keys(): 
            input_files_in_dummy_files += simulations[simulation_id]["input_files"]
        
        # Iterate through the input files 
        for i, input_file in enumerate(input_files_directory):   
            
            # Read the input parameters
            input_parameters, kt_grids_range_parameters = read_relevantInputParameters(input_file)
            nakxnaky = kt_grids_range_parameters["naky"]*kt_grids_range_parameters["nakx"]
            
            # Assume this is a new simulation
            new_simulation = True
            
            # Compare these input parameters with the saved input parameters 
            if nakxnaky==1: 
                for simulation_id in simulations.keys(): 
                    if simulations[simulation_id]["input_parameters"]==input_parameters:
                        
                        # The simulation was probably relaunched, the input file will be removed from this dummy list
                        if input_file in input_files_in_dummy_files and input_file not in simulations[simulation_id]["input_files"]:
                            input_files_in_dummy_files.remove(input_file)
                            for simulation_id in simulations.keys(): 
                                if input_file in simulations[simulation_id]["input_files"]:
                                    simulations[simulation_id]["input_files"].remove(input_file); break
                            
                        # Add a new input file to the list of dummy input files if it's not already there
                        if input_file not in input_files_in_dummy_files:
                            simulations[simulation_id]["input_files"].append(input_file)
                            new_simulation = False 
                        elif input_file in input_files_in_dummy_files:
                            new_simulation = False 
                    
            # If the input parameters didn't match those of a simulation already saved, create a new simulation
            if new_simulation==True:
                
                # The simulation was probably relaunched, the input file will be removed from this dummy list
                if input_file in input_files_in_dummy_files:  
                    input_files_in_dummy_files.remove(input_file)
                    for simulation_id in simulations.keys(): 
                        if input_file in simulations[simulation_id]["input_files"]:
                            simulations[simulation_id]["input_files"].remove(input_file); break
                    
                # Create a new dummy input file list
                if input_file not in input_files_in_dummy_files:
                    simulations[str(input_file)] = {}
                    simulations[str(input_file)]["input_parameters"] = input_parameters
                    simulations[str(input_file)]["input_files"] = [input_file] 
                
        # Get the names of the dummy input files, saved to simulations[simulation_id]["input_file_name"]
        get_namesForDummyInputFiles(simulations)
        
        # Save each dummy input file to a txt file
        save_dummyInputFiles(simulations, folder, directory)

    return
 
#--------------------------------------
def read_relevantInputParameters(input_file):
    """ Returns the <input_parameters> in <input_file> which define unique simulations. """
            
    # Read the input parameters 
    input_parameters = read_inFile(input_file)
    kt_grids_range_parameters = copy.deepcopy(input_parameters['kt_grids_range_parameters'])
    
    # The (kx,ky) values are allowed to be different in the "same" simulation
    del input_parameters['kt_grids_range_parameters']
    
    # Remove other parameters that are allowed to differ
    del input_parameters['knobs']['mat_gen']
    del input_parameters['knobs']['nstep']
    del input_parameters['knobs']['tend']
    del input_parameters['knobs']['delt_option']
    del input_parameters['knobs']['delt_max']
    del input_parameters['knobs']['delt_min']
    del input_parameters['init_g_knobs']['phiinit']
    del input_parameters['init_g_knobs']['restart_file']
    del input_parameters['init_g_knobs']['restart_dir']
    del input_parameters['init_g_knobs']['ginit_option'] 
    del input_parameters['stella_diagnostics_knobs'] 
    return input_parameters, kt_grids_range_parameters

#--------------------------------------
def read_inputFilesInDummyInputFile(dummy_input_file):  
    """ Returns the list of input files in 'input_dummy.in' or 'input_v1_dummy.in'. """
           
    # Read the input files
    with open(dummy_input_file, 'r') as input_data:
        input_text = input_data.read().replace(" ","")
        input_text = input_text.split("&input_files")[-1].split("/\n")[0]
        input_files = input_text.split("\n")
        input_files = [i for i in input_files if i!=""] 
        input_files = [dummy_input_file.parent/i for i in input_files]

    # Check whether the files exist
    for input_file in input_files:
        if not os.path.isfile(input_file):
            input_files.remove(input_file)
    
    # Return the list of input files
    return input_files

#--------------------------------------
def read_inputFilesInDummyInputFiles(dummy_input_files):
    """ Returns the list of input files in 'input_dummy.in' or 'input_v1_dummy.in'. """
    
    # Return an empty list
    if dummy_input_files==None: return []  
    
    # Initiate
    input_files = [] 
    
    # Read the input files
    for dummy_input_file in dummy_input_files:
        with open(dummy_input_file, 'r') as input_data:
            input_text = input_data.read().replace(" ","")
            input_text = input_text.split("&input_files")[-1].split("/\n")[0]
            input_files_temp = input_text.split("\n")
            input_files_temp = [i for i in input_files_temp if i!=""]  
            input_files_temp = [dummy_input_file.parent/i for i in input_files_temp]
        input_files += input_files_temp

    # Check whether the files exist
    for input_file in input_files:
        if not os.path.isfile(input_file): 
            input_files.remove(input_file)
    
    # Return the list of input files
    return input_files


#--------------------------------------
def get_namesForDummyInputFiles(simulations):
    """ For each <simulation_id> in <simulations>, add simulations[simulation_id]["input_file_name"]. """

    # Initialize  
    input_file_names = []
    non_unique_ids = {} 
    
    # Get the existing version numbers 
    existing_ids = get_existingVersionNumbers(simulations) 
    
    # Get the names of the dummy input files 
    for simulation_id in simulations.keys(): 
        
        # Assume the dummy file isn't written yet (then input_file_name is e.g. "input_ky0.125.in") 
        input_file_name = pathlib.Path(simulation_id).name
        version_number = existing_ids[simulation_id]
        dummy_file_is_already_written = True if (version_number>0) else False  
                    
        # Get an appropriate name for the dummy file
        if dummy_file_is_already_written==False:
            if "_dummy.in" in input_file_name: input_file_name = input_file_name.split("_dummy.in")[0]
            if "_kx" in input_file_name: input_file_name = input_file_name.split("_kx")[0]
            if "_ky" in input_file_name: input_file_name = input_file_name.split("_ky")[0]
            if ".in" in input_file_name: input_file_name = input_file_name.split(".in")[0]
            new_version_number = get_newVersionNumber(existing_ids) 
            input_file_name = input_file_name + "_v"+str(new_version_number)+"_dummy.in"
            
        # Save the name for the dummy input file
        simulations[simulation_id]["input_file_name"] = input_file_name
        input_file_names.append(input_file_name)
        
    # Check if there are double names
    if len(list(set(input_file_names)))<len(input_file_names): 
        
        # Save the non unique ID's
        for input_file_name in list(set(input_file_names)): 
            if input_file_names.count(input_file_name)>1:
                non_unique_ids[input_file_name] = 1
                
        # For the non unique ID's, add a "_v1_dummy.in", "_v2_dummy.in", ...
        for simulation_id in simulations.keys():  
            if simulations[simulation_id]["input_file_name"] in non_unique_ids.keys(): 
                input_file_name = simulations[simulation_id]["input_file_name"]  
                new_input_file_name = input_file_name.split("_dummy.in")[0]+"_v"+str(non_unique_ids[input_file_name])+"_dummy.in"
                simulations[simulation_id]["input_file_name"] = new_input_file_name 
                non_unique_ids[input_file_name] += 1   
    return


#--------------------------------------    
def get_existingVersionNumbers(simulations):
    existing_ids = {}
    for simulation_id in simulations.keys(): 
        input_file_name = pathlib.Path(simulation_id).name
        if "_v" in input_file_name and "_dummy.in" in input_file_name:
            for i in range(100):
                if "_v"+str(i)+"_dummy.in" in input_file_name:
                    existing_ids[simulation_id] = i
                    break
        if simulation_id not in existing_ids.keys(): 
            existing_ids[simulation_id] = -1 
    return existing_ids

#--------------------------------------    
def get_newVersionNumber(existing_ids):
    version_numbers = [existing_ids[i] for i in existing_ids.keys()]
    for i in range(1,100): 
        if i not in version_numbers:
            new_version_number = i
            existing_ids[i] = i
            break 
    return new_version_number
                
#--------------------------------------    
def save_dummyInputFiles(simulations, folder, directory):

    # In <directory>, save the list of input files that correspond to the same linear flux-tube simulation 
    for simulation_id in simulations.keys():
          
        # If we only have input file, there's no need for a dummy input file
        if len(simulations)==1:
            if len(simulations[simulation_id]["input_files"])==1:
                continue
          
        # Get the list of input files for <simulation_id> 
        try: input_files = sort_listByNumbers(simulations[simulation_id]["input_files"])
        except: pass 
          
        # Copy the text of one of the input files to <input_file_text>
        input_file = open(simulations[simulation_id]["input_files"][0], "r" )
        input_file_text = copy.deepcopy(input_file.read())
        input_file.close() 
          
        # Add a knob for the input files to <input_file_text>
        knob_input_files = "&input_files \n" 
        for input_file in input_files:
            path_input_file = " "+input_file.name if input_file.parent==directory else " "+input_file.parent.name+"/"+input_file.name # could be "OLD/input.in"
            knob_input_files += path_input_file + "\n"
      
        # Add the knob to the input file
        input_file_text = knob_input_files + "/ \n" + input_file_text
       
        # Write the new input file and close it
        new_file = open(directory / simulations[simulation_id]["input_file_name"], "w" )  
        new_file.write(input_file_text) 
        new_file.close()   
              
        # Message
        length = len(input_files); name = simulations[simulation_id]["input_file_name"]
        path = simulations[simulation_id]["input_files"][0].parent 
        path = path if "OLD" not in path.name else path.parent
        path = folder.name+str(path).split(folder.name)[-1]
        if length==1: print("   Combined "+str(length)+" input file inside "+path+"/"+name+".") 
        if length>1: print("   Combined "+str(length)+" input files inside "+path+"/"+name+".") 
