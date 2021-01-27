#!/usr/bin/python3
''' Run kyscan(identifier) as a bash command on the command prompt.'''

# Tell python where to find the personal modules since this script is run from the terminal
import sys, os, pathlib
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.utils import get_filesInFolder 
from stellapy.commands.utils.get_bashArguments import get_bashArguments
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.utils.decorators import print_functionInformationOnCommandPrompt

#################################################################
#                     MAIN PROGRAM
#################################################################


def restart(folder): 
    
    # Check whether the correct files are present
    input_file, default_text = check_presenceFiles(folder)
    if input_file==None: return 
    
    # If we have a nonlinear simulation, make sure there is a restart folder
    restart_folder = None; count=2
    if "restart1" not in os.listdir(folder) or "restart" not in os.listdir(folder):
        restart_folder = "restart1" 
        os.system("mkdir restart1")
    else:
        while (restart_folder==None):
            if "restart"+str(count) not in os.listdir(folder):
                restart_folder = "restart"+str(count)
                os.system("cp -r restart"+str(count-1) + " restart"+str(count))
            count += 1
    
    # Move the current files to a subfolder
    print("#####################################################################")
    print("                       MOVE FILES TO A SUBFOLDER                     ")
    print("#####################################################################")
    os.system("$STELLAPY/commands/runSimulations/move_simulationToSubfolder.sh")
    print("\n")
        
    # Change the input file
    print("#####################################################################")
    print("                EDIT INPUT FILE AND RELAUNCH SIMULATION              ")
    print("#####################################################################")
    create_inputFile(folder, input_file, default_text, restart_folder)
 
    # Now launch the simulation
    print("sbatch -D " + str(folder) +" eu.slurm    --->     "+input_file.stem)
    os.system("sbatch -D " + str(folder) +" eu.slurm")
    print()

    return

#################################################################
#                     CREATE THE INPUT FILE
#################################################################

def create_inputFile(folder, input_file, default_text, restart_folder):
    
    # If ginit_option==noise then change the restart options:
    if "noise" in default_text:
        
        # Create the new input file
        new_file = open(input_file, "w" )
        
        # Replace ginit_option
        text_before = default_text.split("ginit_option")[0]
        text_after  = default_text.split("ginit_option")[-1].split("\n")[1:]
        text_after  = '\n'.join(text_after)
        new_text    = text_before + 'ginit_option = "many"' + "\n" + text_after
        
        # Use a new restart directory
        text_before = new_text.split("restart_dir")[0]
        text_after  = new_text.split("restart_dir")[-1].split("\n")[1:]
        text_after  = '\n'.join(text_after)
        new_text    = text_before + 'restart_dir = "'+restart_folder+'"' + "\n" + text_after
        
        # Replace delt
        if "delt=" in new_text:   delt = "delt="
        if "delt =" in new_text:  delt = "delt ="
        if "delt  =" in new_text: delt = "delt  ="
        text_before = new_text.split(delt)[0]
        text_after  = new_text.split(delt)[-1].split("\n")[1:]
        text_after  = '\n'.join(text_after)
        new_text    = text_before + 'delt_option = "check_restart"' + "\n" + text_after
        
        # Make sure write moments is off
        if "write_moments" in new_text:
            text_before = new_text.split("write_moments")[0]
            text_after  = new_text.split("write_moments")[-1].split("\n")[1:]
            text_after  = '\n'.join(text_after)
            new_text    = text_before + 'write_moments = .false.' + "\n" + text_after

        # Write the new input file
        new_file.write(new_text) 
        new_file.close() 

    return
                    
#################################################################
#          CHECK WHETHER THE CORRECT FILES ARE PRESENT
#################################################################

def check_presenceFiles(folder):
    
    # Select the default input file
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [f for f in input_files if "run0" not in str(f)]
    input_files = [f for f in input_files if "run1" not in str(f)]
    
    # If there are multiple default input files: exit
    if len(input_files) > 1:
        print(" Multiple default input files were found inside "+str(folder)+":")
        for f in input_files:
            print("     ", f)
        print(" Please make sure there is only one default input file.\n")
        return None, None, None
    
    # If there are no input_files: exit
    if len(input_files) == 0:
        print(" No input file found in "+str(folder)+".\n")
        return None, None, None
    
    # Read the magnetic field file name and make sure the file is in the folder
    input_data = open(input_files[0], 'r')
    input_text = input_data.read() 
    vmec_file  = str(input_text.split("vmec_filename")[-1].split("\n")[0])
    vmec_file  = vmec_file.replace(" ", "").replace("=", "").replace("'", "")
    if not os.path.isfile(folder / vmec_file) and not os.path.islink(folder / vmec_file):
        print(" The vmec file "+vmec_file+" was not found in "+str(folder)+".\n")
        return None, None, None
    
    # Make sure we have a stella executable
    if not os.path.isfile(folder / "stella") and not os.path.islink(folder / "stella"):
        print(" No stella executable found in "+str(folder)+".\n")
        return None, None, None
    
    # Read the default input file
    input_file   = input_files[0]
    default_file = open(input_file, "r" )
    default_text = default_file.read()
    default_file.close()
    
    # If we made it true, return the input file text
    return input_file, default_text
    

#################################################################
#                     BASH INTERFACE
#################################################################

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        function = "runSimulations.restart"
        default_arguments = {\
                'folder': None}
        arguments = get_arguments(function, default_arguments) 
        restart(**arguments) 

# Explain the program
def print_help(function, default_arguments):
    print_decoratorOnCommandPrompt("begin")
    print_functionInformationOnCommandPrompt(function, default_arguments, None)
    print("\n The arguments are:") 
    #print("      -i         <identifier>" )  
    print_decoratorOnCommandPrompt("end")

# Get the arguments from the terminal
def get_arguments(function, default_arguments):

    # Make a new dictionary which will hold the new arguments
    arguments = default_arguments.copy()

    # Initiate mandetory variables
    arguments['folder'] = pathlib.Path(os.getcwd())
    
    # Read the options and arguments from the terminal
    options = "h"
    long_options = ["folder=","help"]
    opts_args = get_bashArguments(options, long_options, print_help, function, default_arguments)

    # Asign the options and arguments to variables
    for opt, arg in opts_args:
        if opt in ("-h", "--help"):
            print_help(function, default_arguments)
            sys.exit()
        elif opt == '--folder':
            arguments['folder'] = str(arg)
            
    # Show the parameters
    print_functionInformationOnCommandPrompt(function, default_arguments, arguments)
    return arguments

# Execute the python script
main_program()
