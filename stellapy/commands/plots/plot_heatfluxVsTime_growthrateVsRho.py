#!/usr/bin/python3
''' Run plot_heatfluxVsTime_growthrateVsRho() as a bash command on the command prompt.'''

# Load the modules
import sys, os

# Tell python where to find the personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.bashscripts.get_bashArguments import get_bashArguments
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.utils.decorators import print_functionInformationOnCommandPrompt

# Get the python function we want to execute here
from stellapy.plot.fluxes import plot_heatfluxVsTime_growthrateVsRho

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        function = "plot.fluxes.plot_heatfluxVsTime_growthrateVsRho"
        default_arguments = {\
              'folder' : "MANDETORY ARGUMENT", \
              'save' : True, \
              'log' : True, \
              'versus_time' : True}
        arguments = get_arguments(function, default_arguments)
        plot_heatfluxVsTime_growthrateVsRho(**arguments)

# Explain the program
def print_help(function, default_arguments):
    print_decoratorOnCommandPrompt("begin")
    print_functionInformationOnCommandPrompt(function, default_arguments, None)
    print("\n The optional arguments are:")
    print("      --folder      <folder>" )
    print_decoratorOnCommandPrompt("end")

# Get the arguments from the terminal
def get_arguments(function, default_arguments):
    
    # Make a new dictionary which will hold the new arguments
    arguments = default_arguments.copy()

    # Initiate mandetory variables
    arguments['folder'] = os.getcwd().split('RUNS/')[-1]

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





    




