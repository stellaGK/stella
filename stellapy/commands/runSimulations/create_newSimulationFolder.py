#!/usr/bin/python3
''' Run create_newSimulationFolder() as a bash command on the command prompt.'''

# Tell python where to find the personal modules since this script is run from the terminal
import sys, os, shutil
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.config import CONFIG
from stellapy.utils.decorators import exit_program  
from stellapy.data import write_profile
from stellapy.utils import get_filesInFolder
from stellapy.commands.utils.get_bashArguments import get_bashArguments
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.utils.decorators import print_functionInformationOnCommandPrompt

#################################################################
#                     MAIN PROGRAM
#################################################################

# TODO: Automatically add fprim and tprim if rho changes
# TODO: Automatically make restart folders for nonlinear simulations

def create_newSimulationFolder(ref_folder=None, folder_extension=None, research=None, specific_research=None):
    ''' Create a new folder to start a new simulations based on the folder "<ref_folder>+ref".
    
    The files in the reference folder will be copied to a new folder "<ref_folder> + research".
    This folder needs to be placed in the path "NEWRUNS" specified in the configuration file.
    
    In case only one input file is present, this file will be copied. In case multiple input
    files are present with the prefix "input", one of these files needs to be selected by 
    setting <research> to the string coming after "input_".
    
    In case <specific_research> is specified, the <research> string is assumed to refer
    to the variable which will be varied in multiple simulations. For example if <research>
    is "rho_*" the value of torflux will be varied and a folder will be created for each
    value of torflux that is specified in the list <specific_research>. The specific value
    will be added to the folder name and the value will be automatically filled in the input file.
    
    The reference folder <ref_folder> needs to contain the following standard files: 
        - "eu.slurm"        -->    Used by the supercomputer to launch the simulation
        - "wout*.nc"        -->    Magnetic field file
        - "*.in"            -->    Standard input file
        - "*profile*.dat"   -->    The raw profile data 
        
    Parameters
    ----------
    ref_folder : str
        Path of the reference folder which contains the standard files.
        
    folder_extension : str, optional
        In case research and specific_research are not specified, this will be added
        behind the <ref_folder> as the name of the new folder. Otherwise the folder
        extension is research+specific_research.
    
    research : str, optional
        String that comes after "input_" in "*.in" to know which input file is selected. 
        It can also be used to specify which variable is varied in the input parameters 
        of stella, in this case only the string before the next "_" is used.
        
    specific_research : float or list of floats, optional
        The specific values that are selected to scan the parameter defined in 
        the first part of the string "research.
        
        
    Example bash command
    --------------------
    create_newSimulationFolder -r "rho_kyscan_ad" -s "[0.5, 0.6]"
    '''
    
    # Get the absolute path of the reference folder 
    path_refFolder = ref_folder
    
    # It needs to be performed in a folder inside the $NEWRUNS directory
    newruns_directory = CONFIG['PATHS']['NEWRUNS']
    if (newruns_directory not in path_refFolder):
        exit_program("Please use this command on a folder inside CONFIG['PATHS']['NEWRUNS'].", create_newSimulationFolder, sys._getframe().f_lineno) 
    elif len(path_refFolder.split(newruns_directory)[-1])<2:
        exit_program("Please use this command on a folder inside CONFIG['PATHS']['NEWRUNS'].", create_newSimulationFolder, sys._getframe().f_lineno) 
        
    # If specific_research is specified: create a new folder for each value
    if specific_research is None:
        specific_research = [""]
    if type(specific_research) is not list:
        specific_research = [specific_research]
    
    # Count the number of input files
    input_files = [file for file in os.listdir(path_refFolder) if ((".in" in file) and not "~" in file)]
    numberOfInputfiles = len(input_files)
        
    # Iterate over the new folders that need to be created
    for specific_value in specific_research:
        
        # Define the new folder name
        if "_ref" in ref_folder: 
            ref_folder = ref_folder.split("_ref")[0]
        if not folder_extension and (research and specific_research):
            folder_extension_temp = "_" + research.split("_")[0] + str(specific_value)
            folder_name = ref_folder + str(folder_extension_temp)
        elif folder_extension:     
            folder_name = ref_folder + "_" + str(folder_extension)
        elif not folder_extension:
            exit_program("Please specify the name of the new folder.", create_newSimulationFolder, sys._getframe().f_lineno) 
        
        # Get the absolute path of the new folder
        path_newFolder = path_refFolder.split(ref_folder)[0] + folder_name
        
        # Create a new folder for the simulation
        if not os.path.exists(path_newFolder):
            os.makedirs(path_newFolder)
            
        # Make sure we have the profile.txt file if there is a profile.dat file
        profile_files = [file for file in os.listdir(path_refFolder) if "profile" in file]
        if len (profile_files) > 0: 
            raw_profile   = [file for file in profile_files if ".dat" in file][0]
            if not "profile.txt" in profile_files:
                write_profile(source_file=raw_profile, folder=path_refFolder)
            
        # Copy the data to the newfolder 
        for file in os.listdir(path_refFolder):  
            if ("~" not in file) and (".dat" not in file) and file.startswith(".")==False:
                if ("input") not in file \
                or numberOfInputfiles==1 \
                or (research is not None) and (research in file):
                    old_file = path_refFolder + "/" + file.split('/')[-1]
                    new_file = path_newFolder + "/" + file.split('/')[-1]
                    shutil.copyfile(old_file, new_file, follow_symlinks=False)
                    
        # Find the name of the magnetic field file
        for file in os.listdir(path_newFolder): 
            if (".in" in file) and (not "~" in file):
                input_file = open(path_newFolder + "/" + file, "r" )
                text_input = input_file.read()
                vmec_filename = text_input.split("vmec_filename")[-1].split("\n")[0].split("=")[-1]
                vmec_filename = vmec_filename.replace("'","").replace('"','').replace(' ','')
                input_file.close()

        # Add symbolic links to the stella code and the magnetic field file that will work on the supercomputer 
        source = CONFIG['PATHS']['RUNS_SUPERCOMPUTER'] 
        if not os.path.islink(path_newFolder+"/stella"):
            os.symlink(source+"/stella", path_newFolder+"/stella")
            os.symlink(source+"/"+vmec_filename, path_newFolder+"/"+vmec_filename)
        
        # Add the values of <specific_research> to the input file
        input_file = None
        if specific_research != ['']:
            for file in os.listdir(path_newFolder): 
                if (".in" in file) and (not "~" in file):
                    
                    # Remember the input file for the next task
                    input_file = file
                    
                    # Make a new file and save/read the old one
                    path_file = path_newFolder + "/" + file
                    shutil.move(path_file, path_file+"~" )
                    destination = open(path_file, "w" )
                    source = open(path_file+"~", "r" )
                    source_text = source.read()
                    
                    # Get the variable
                    variable = research.split("_")[0]
                    
                    # Sometimes its a synonym
                    if variable == "rho":    variable = "torflux"
                    
                    # Check if this variable is in the input file
                    if variable not in source_text:
                        exit_reason = "Could not find the paramater "+variable+' in '+file+'.'
                        exit_program(exit_reason, create_newSimulationFolder, sys._getframe().f_lineno) 
                    else:
                        text_before = source_text.split(variable)[0]+variable
                        text_after  = source_text.split(variable)[-1].split("\n")[1:]
                        text_after  = '\n'.join(text_after)
                        new_text    = text_before + " = " + str(specific_value) + "\n" + text_after
                        destination.write(new_text)  
                        source.close()
                        destination.close()
                    
        # If input_file wasn't defined assume we only have one: grab that one
        if input_file==None:
            input_file = get_filesInFolder(path_newFolder, end=".in")[0]
                    
        # Automatically add the input_file name to the eu.slurm file
        path_file = path_newFolder + "/eu.slurm" 
        shutil.move(path_file, path_file+"~" )
        destination = open(path_file, "w" )
        source = open(path_file+"~", "r" )
        source_text = source.read()
        source_text = source_text.split("$EXEFILE")[0] + "$EXEFILE " + input_file
        destination.write(source_text)  
        source.close()
        destination.close()
    
#################################################################
#                     BASH INTERFACE
#################################################################

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        function = "simulations.create_newSimulationFolder"
        default_arguments = {\
                'ref_folder': None, \
                'folder_extension': None, \
                'research': None, \
                'specific_research': None}
        arguments = get_arguments(function, default_arguments) 
        create_newSimulationFolder(**arguments)

# Explain the program
def print_help(function, default_arguments):
    print_decoratorOnCommandPrompt("begin")
    print_functionInformationOnCommandPrompt(function, default_arguments, None)
    print("\n The optional arguments are:")
    print("      -r         <research>" )
    print("      -s         <specific_research>" )
    print("      -e         <folder exntension>" )
    print("      --folder   <case folder>" )
    print_decoratorOnCommandPrompt("end")

# Get the arguments from the terminal
def get_arguments(function, default_arguments):

    # Make a new dictionary which will hold the new arguments
    arguments = default_arguments.copy()

    # Initiate mandetory variables
    arguments['ref_folder'] = os.getcwd()
    
    # Read the options and arguments from the terminal
    options = "hr:s:e:"
    long_options = ["folder=","help"]
    opts_args = get_bashArguments(options, long_options, print_help, function, default_arguments)

    # Asign the options and arguments to variables
    for opt, arg in opts_args:
        if opt in ("-h", "--help"):
            print_help(function, default_arguments)
            sys.exit()
        elif opt == '--folder':
            arguments['ref_folder'] = str(arg)
        if opt in ("-e"):
            arguments['folder_extension'] = arg
        if opt in ("-r"):
            arguments['research'] = arg
        if opt in ("-s"):
            arguments['specific_research'] = arg
            if "[" in arguments['specific_research']:
                arguments['specific_research'] = arguments['specific_research'].split('[')[-1]
                arguments['specific_research'] = arguments['specific_research'].split(']')[0]
                arguments['specific_research'] = [float(i) for i in arguments['specific_research'].split(",")]  

    # Show the parameters
    print_functionInformationOnCommandPrompt(function, default_arguments, arguments)
    return arguments

# Execute the python script
main_program()





    




