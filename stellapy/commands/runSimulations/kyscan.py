#!/usr/bin/python3
''' Run kyscan(identifier) as a bash command on the command prompt.'''

# Tell python where to find the personal modules since this script is run from the terminal
import sys, os, shutil, configparser, pathlib
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.utils import get_filesInFolder 
from stellapy.commands.utils.get_bashArguments import get_bashArguments
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.utils.decorators import print_functionInformationOnCommandPrompt

#################################################################
#                     MAIN PROGRAM
#################################################################


def kyscan(folder, identifier="ten", confirmation=False):
    ''' Simulate one mode per file for different values of ky. 
    
    Make sure we have the required files: wout.nc, eu.slurm, input.in and stella. '''
    
    # Select the default input file
    input_files = get_filesInFolder(folder, end=".in")
    input_files = [f for f in input_files if "_ky" not in str(f)]
    
    # If there are multiple default input files: exit
    if len(input_files) > 1:
        print(" Multiple default input files were found inside "+str(folder)+":")
        for f in input_files:
            print("     ", f)
        print(" Please make sure there is only one default input file.\n")
        return
    
    # If there are no input_files: exit
    if len(input_files) == 0:
        print(" No input file found in "+str(folder)+".\n")
        return
    
    # Read the magnetic field file name and make sure the file is in the folder
    input_data = open(input_files[0], 'r')
    input_text = input_data.read() 
    vmec_file  = str(input_text.split("vmec_filename")[-1].split("\n")[0])
    vmec_file  = vmec_file.replace(" ", "").replace("=", "").replace("'", "")
    if not os.path.isfile(folder / vmec_file) and not os.path.islink(folder / vmec_file):
        print(" The vmec file "+vmec_file+" was not found in "+str(folder)+".\n")
        return
    
    # Make sure we have the required files: eu.slurm and stella
    if not get_filesInFolder(folder, end=".slurm"):
        print(" No *.pbs or *.slurm file found in "+str(folder)+".\n")
        return
    if not os.path.isfile(folder / "stella") and not os.path.islink(folder / "stella"):
        print(" No stella executable found in "+str(folder)+".\n")
        return
        
    # Read the list of ky modes based on the "stellapy/config/kyscan.ini" file
    ky_file = configparser.ConfigParser() 
    ky_path = pathlib.Path(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
    ky_path = ky_path / "stellapy/config/kyscan.ini"
    ky_file.read(ky_path)
    if identifier in ky_file["SCANS"].keys():
        ky_modes = ky_file["SCANS"][identifier]
        ky_modes = ky_modes.replace("\n", "").replace("\t", "").replace(" ", "")
        ky_modes = ky_modes.split("[")[-1].split("]")[0].split(",")  
    else:
        print(' The identifier "'+identifier+'" is not present in stellapy/config/kyscan.ini, please add it')
        print(' or choose one of the existing identifiers:\n')
        for key in ky_file["SCANS"].keys():
            ky_modes = ky_file["SCANS"][key]
            ky_modes = ky_modes.replace("\n", "").replace("\t", "").replace(" ", "")
            ky_modes = ky_modes.split("[")[-1].split("]")[0].split(",") 
            ky_modes = "[" + ", \n\t\t".join([", ".join(ky_modes[i:i+10]) for i in range(0,len(ky_modes),10)]) + "]"
            print("     ", key)
            print("\t\t"+ky_modes)
        print()
        return
    
    # Ask for confirmation to proceed: so we can check the modes before launching
    if not confirmation:
        modes = len(ky_modes)
        ky_modes = "[" + ", \n\t".join([", ".join(ky_modes[i:i+10]) for i in range(0,len(ky_modes),10)]) + "]"
        print(" Are you sure you want to launch "+str(modes)+" simulations for the following ky values?")
        print("\t"+ky_modes)
        print(' If yes, confirm the selection with "kyscan -c" or "kyscan -ci <identifier>".')
        print(' In order to list the different ranges of ky values use "kyscan -i help".\n')
        return
        
    # Read the default input file
    input_file   = input_files[0]
    default_file = open(input_file, "r" )
    default_text = default_file.read()
    
    # Make sure n_aky and n_akx are set to 1
    nakx = int(default_text.split("nakx")[-1].split("\n")[0].replace("=", ""))
    naky = int(default_text.split("naky")[-1].split("\n")[0].replace("=", ""))
    if nakx!=1 or naky!=1:
        print(" Please set nakx and naky to 1 if you want to run one mode per input file.\n")
        return

    # Replace the ky value in ".in", replace the EXEFILE in ".slurm" and launch the simulation
    for ky in ky_modes: 
        
        # Replace the ky value in each input file
        file_name = str(input_file.stem)+"_ky"+ky+".in"
        new_file = open(input_file.parent / file_name, "w" )
        text_before = default_text.split("aky_min")[0]
        text_after  = default_text.split("aky_min")[-1].split("\n")[1:]
        text_after  = '\n'.join(text_after)
        new_text    = text_before + "aky_min = " + ky + "\n" + text_after
        text_before = new_text.split("aky_max")[0]
        text_after  = new_text.split("aky_max")[-1].split("\n")[1:]
        text_after  = '\n'.join(text_after)
        new_text    = text_before + "aky_max = " + ky + "\n" + text_after
        new_file.write(new_text) 
        new_file.close() 
        print("Created "+file_name)
        
        # Automatically add the input_file name to the eu.slurm file 
        shutil.move(input_file.parent / "eu.slurm", input_file.parent / "eu.slurm~" )
        new_file = open(input_file.parent / "eu.slurm", "w" )
        old_file = open(input_file.parent / "eu.slurm~", "r" )
        old_text = old_file.read()
        new_text = old_text.split("$EXEFILE")[0] + "$EXEFILE " + file_name
        new_file.write(new_text)  
        new_file.close()
        old_file.close()
        
        # Now launch the simulation
        os.system("sbatch eu.slurm") # qsub is an alias which isn't loaded in the python3 shell
        print()
        
    print()
    default_file.close()
    
    print("Finished launching the "+str(len(ky_modes))+" ky modes.\n")


#################################################################
#                     BASH INTERFACE
#################################################################

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        function = "runSimulations.kyscan"
        default_arguments = {\
                'folder': None, \
                'identifier': "ten", \
                'confirmation': False}
        arguments = get_arguments(function, default_arguments) 
        kyscan(**arguments) 

# Explain the program
def print_help(function, default_arguments):
    print_decoratorOnCommandPrompt("begin")
    print_functionInformationOnCommandPrompt(function, default_arguments, None)
    print("\n The arguments are:") 
    print("      -i         <identifier>" )
    print("      -c         Sets confirmation=True" )  
    print_decoratorOnCommandPrompt("end")

# Get the arguments from the terminal
def get_arguments(function, default_arguments):

    # Make a new dictionary which will hold the new arguments
    arguments = default_arguments.copy()

    # Initiate mandetory variables
    arguments['folder'] = pathlib.Path(os.getcwd())
    
    # Read the options and arguments from the terminal
    options = "hci:"
    long_options = ["folder=", "identifier=","help"]
    opts_args = get_bashArguments(options, long_options, print_help, function, default_arguments)

    # Asign the options and arguments to variables
    for opt, arg in opts_args:
        if opt in ("-h", "--help"):
            print_help(function, default_arguments)
            sys.exit()
        elif opt == '--folder':
            arguments['folder'] = str(arg)
        if opt in ("-i", "--identifier"):
            arguments['identifier'] = arg
        if opt in ("-c"):
            arguments['confirmation'] = True 
            
    # Show the parameters
    print_functionInformationOnCommandPrompt(function, default_arguments, arguments)
    return arguments

# Execute the python script
main_program()
