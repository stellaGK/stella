#!/usr/bin/python3
''' Run reduce_sizeNetcdf() as a bash command on the command prompt.'''

# Tell python where to find the personal modules
import sys, os, pathlib

# Tell python where to find the personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.commands.utils.get_bashArguments import get_bashArguments
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.utils.decorators import print_functionInformationOnCommandPrompt
from stellapy.config import CONFIG

# Get the python function we want to execute here
from stellapy.data import reduce_sizeNetcdf 

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        function = "data.utils.reduce_sizeNetcdf"
        default_arguments = {\
                'folder' : "MANDETORY ARGUMENT", \
                'write_again' : True}
        arguments = get_arguments(function, default_arguments) 
        reduce_sizeNetcdf(**arguments)

# Explain the program
def print_help(function, default_arguments):
    print_decoratorOnCommandPrompt("begin")
    print_functionInformationOnCommandPrompt(function, default_arguments, None)
    print_decoratorOnCommandPrompt("end")

# Get the arguments from the terminal
def get_arguments(function, default_arguments):

    # Make a new dictionary which will hold the new arguments
    arguments = default_arguments.copy()

    # Initiate mandetory variables
    arguments['folder'] = pathlib.Path(os.getcwd())
    
    # Read the options and arguments from the terminal
    options = "hr:s:"
    long_options = ["help","folder=","dontwriteagain"]
    opts_args = get_bashArguments(options, long_options, print_help, function, default_arguments)

    # Asign the options and arguments to variables
    for opt, arg in opts_args:
        if opt in ("-h", "--help"):
            print_help(function, default_arguments)
            sys.exit()
        elif opt == '--folder':
            arguments['folder'] = pathlib.Path(arg)  
        elif opt == '--dontwriteagain':
            arguments['write_again'] = False

    # Don't perform this command in the following folders
    if "not_important" in str(arguments['folder'])\
    or "failed" in str(arguments['folder'])\
    or "figure" in str(arguments['folder'])\
    or "restart" in str(arguments['folder'])\
    or "oldddd" in str(arguments['folder']):
        arguments['folder'] = None

    # Show the parameters
    if CONFIG['DEFAULT']['use_verbosewrapper'] != "False":
        print_functionInformationOnCommandPrompt(function, default_arguments, arguments)
    return arguments

# Execute the python script
main_program()






    




