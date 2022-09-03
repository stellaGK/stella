 
import numpy as np 
import os, pathlib, configparser
from stellapy.utils.decorators.verbose import noverbose   
from stellapy.data.input.load_inputObject import save_inputFile
from stellapy.data.input.read_inputFile import read_inputFile, read_inFile 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder 
from stellapy.utils.files import remove_simulationsWithoutOutncFile

#===============================================================================
#                          WRITE THE GEOMETRY FILE
#===============================================================================
 
@noverbose
def write_iniFileForInputs(folder=None, input_files=None):  
    
    # Get the input_files
    if not input_files: input_files = get_filesInFolder(folder, end=".in")
    
    # Only look at the input files in each parent folder
    for folder in list(set([i.parent if "OLD" not in i.parent.name else i.parent.parent for i in input_files])):    
        
        # Initialize
        modes = []
        
        # Get the files in this folder
        inputfiles = [i for i in input_files if (i.parent==folder) or (i.parent.parent==folder and "OLD" in i.parent.name)]
        inputfiles = remove_simulationsWithoutOutncFile(inputfiles)   
        if inputfiles==[]: continue
        
        # Create a simulation for this folder 
        from stellapy.simulations.Simulation import create_simulations
        simulations = create_simulations(None, inputfiles, write_uniqueFiles=True)  
          
        # For each <mode> save the input to an *.unique.inputX file 
        if simulations[0].linear:   
            for simulation in simulations: 
                modes += simulation.modes
            unique_modes = get_uniqueInputFiles(modes)
            for mode in unique_modes:   
                if not os.path.isfile(mode.path.input):
                    save_inputFile(mode)    
             
        # For each <simulation> save the input to an *.input file  
        if simulations[0].nonlinear:  
            for simulation in simulations:  
                if not os.path.isfile(simulation.path.input): 
                    simulation.inputParameters = read_inFile(simulation.path.input_file)
                    save_inputFile(simulation)  

    return
    
#-------------------------------
def get_uniqueInputFiles(modes):

    # Save each mode in unique_inputs[i] = [mode1, ..., modeX] 
    unique_inputs = {}
    
    # Read the input parameters for each <mode>
    for mode in modes: 
         
        # Read the input parameters and remove the variables that allowed to differ between different <modes> 
        mode.inputParameters = read_inputFile(mode.input_file)
        del mode.inputParameters["kt_grids_range_parameters"]["akx_min"]
        del mode.inputParameters["kt_grids_range_parameters"]["akx_max"]
        del mode.inputParameters["kt_grids_range_parameters"]["aky_min"]
        del mode.inputParameters["kt_grids_range_parameters"]["aky_max"]
        del mode.inputParameters["stella_diagnostics_knobs"]["save_for_restart"]
        del mode.inputParameters["stella_diagnostics_knobs"]["nwrite"]
        del mode.inputParameters["stella_diagnostics_knobs"]["nsave"] 
        del mode.inputParameters["init_g_knobs"]["restart_file"] 
        del mode.inputParameters["init_g_knobs"]["restart_dir"] 
        del mode.inputParameters["knobs"]["tend"]
        
    # Find the unique input files
    for mode in modes: 
        for ireference in unique_inputs.keys():
            if mode.inputParameters==unique_inputs[ireference][0].inputParameters: 
                unique_inputs[ireference].append(mode)
                break 
        else: 
            unique_inputs[len(list(unique_inputs.keys()))+1] = [mode] 
    
    # Save the unique input files to a configuration file 
    unique_modes = save_uniqueInputFiles(modes, unique_inputs)

    # Return a reference input file for each geometry 
    return unique_modes

#-------------------------------
def save_uniqueInputFiles(modes, unique_inputs, difference={}, unique_modes=[]):
 
    # Define the configuration file
    ini_path = modes[0].path.folder / (modes[0].path.name +".list.inputs.ini") 
    ini_data = configparser.ConfigParser() 
    
    # Give some general information: the amount of unique inputs and the parent folder
    ini_data.add_section('Unique Input Files')
    ini_data['Unique Input Files']['folder'] = modes[0].path.folder.name
    ini_data['Unique Input Files']['amount'] = count = str(len(list(set(list(unique_inputs.keys())))))
    if count!="1": print("  We found "+count+" unique input files inside "+modes[0].path.folder.name+" containing "+str(len(modes))+" modes.") 
    if count=="1": print("  We found "+count+" unique input file inside "+modes[0].path.folder.name+" containing "+str(len(modes))+" modes.") 

    # Find the knobs and values that are varied between the input files
    keys = list(unique_inputs.keys())
    knobs = list(modes[0].inputParameters.keys())
    for i, j in ((i,j) for i in keys for j in keys): 
        iparameters = unique_inputs[i][0].inputParameters
        jparameters = unique_inputs[j][0].inputParameters
        for knob in knobs:
            if iparameters[knob]!=jparameters[knob]:
                if knob not in difference.keys():
                    difference[knob] = {}
                for parameter in iparameters[knob].keys(): 
                    if iparameters[knob][parameter]!=jparameters[knob][parameter]:
                        if parameter in difference[knob].keys():
                            difference[knob][parameter].append(iparameters[knob][parameter])
                            difference[knob][parameter].append(jparameters[knob][parameter])
                        if parameter not in difference[knob].keys():
                            difference[knob][parameter] = [iparameters[knob][parameter]]
                            difference[knob][parameter].append(jparameters[knob][parameter])
    
    # Save the information of the varied parameters to the configuration file
    ini_data.add_section('Varied parameters')    
    if len(list(difference.keys()))==0:
        ini_data['Varied parameters']['None'] = str([])
    for knob in difference.keys():
        for parameter in difference[knob].keys(): 
            ini_data['Varied parameters'][parameter] = str(sorted(list(set(difference[knob][parameter]))))
    
    # Save the input_files that belong to each reference input file
    ireferences = sorted(list(unique_inputs.keys()))
    for ireference in ireferences:
        
        # Make a section for each unique input file
        ini_data.add_section("Input "+str(ireference)) 
        
        # Save the unique mode and create its file path
        reference_mode = unique_inputs[ireference][0]
        unique_modes.append(reference_mode)
        reference_mode.path.input = reference_mode.path.folder / (reference_mode.path.name + ".unique.input"+str(ireference)+".ini")
        
        # Sort the modes
        modes = unique_inputs[ireference]
        numbers = [mode.ky for mode in modes]
        sorted_indexes = list(np.array(numbers).argsort()) 
        modes = [modes[i] for i in sorted_indexes] 
        
        # Add the modes that belong to this input file 
        for mode in modes:
            mode.path.input = reference_mode.path.input
            name = "v0 ("+str(mode.kx)+", "+str(mode.ky)+")"; i=0
            while name in ini_data["Input "+str(ireference)]: name = name.replace("v"+str(i),"v"+str(i+1)); i += 1 
            ini_data["Input "+str(ireference)][name] = str(mode.input_file)

    # Save the list of unique inputs as a configuration file 
    with open(ini_path, 'w') as ini_file:
        ini_data.write(ini_file) 
    
    # Return the unique input files  
    return unique_modes
 

################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    folder = pathlib.Path("/home/hanne/CIEMAT/PREVIOUSRUNS/LINEARMAPS/W7Xstandard_rho0.7_aLTe0/LinearMap/fprim4tprim4")
    write_iniFileForInputs(folder)

