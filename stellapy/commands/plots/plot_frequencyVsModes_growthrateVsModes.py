#!/usr/bin/python3
''' Run plot_frequencyVsModes_growthrateVsModes() as a bash command on the command prompt.'''

# Load the modules
import sys, os

# Tell python where to find the personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.bashscripts.get_bashArguments import get_bashArguments
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.utils.decorators import print_functionInformationOnCommandPrompt

# Get the python function we want to execute here
from stellapy.plot.frequencyAndGrowthrate import plot_frequencyVsModes_growthrateVsModes

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        function = "plot.frequencyAndGrowthrate.plot_frequencyVsModes_growthrateVsModes"
        default_arguments = {\
            'folder' : 'MANDETORY ARGUMENT', \
            'linear_data' : 'MANDETORY ARGUMENT', \
            'input_values' : 'MANDETORY ARGUMENT', \
            'y_quantity' : 'omega', \
            'ax' : None, \
            'kx' : -0.0, \
            'ki' : 'ky', \
            'k_value' : 9.0, \
            'k_delta' : 1.0, \
            'specific_rho' : '0.5', \
            'ymin' : None, \
            'ymax' : None, \
            'x_range' : None, \
            'SI' : False, \
            'verbose' : True, \
            'save' : False}
        arguments = get_arguments(function, default_arguments)
        plot_frequencyVsModes_growthrateVsModes(**arguments)
        
# Explain the program
def print_help(function, default_arguments):
    print_decoratorOnCommandPrompt("begin")
    print_functionInformationOnCommandPrompt(function, default_arguments, None)
    print("\n The optional arguments are:")
    print("      --folder      <case folder>" )
    print("      --ki          <ki>                 {'ky', 'kx'}")
    print("      --kx          <plot at this kx>    {-0.0}")
    print("      --kvalue      <k reflectometry>    {9.0}")
    print("      --kdelta      <k reflectometry>    {1.0}")
    print("      --xrange      <xmax>               {[0, xmax]}")
    print("\n The options are:")
    print("      --SI          --noSI           Plot gamma and ki in SI units")
    print("      --save         --nosave        Toggle savong the figure (for interactive mode)")
    print_decoratorOnCommandPrompt("end")

# Get the arguments from the terminal
def get_arguments(function, default_arguments):
    
    # Make a new dictionary which will hold the new arguments
    arguments = default_arguments.copy()

    # Initiate mandetory variables
    arguments['folder'] = os.getcwd().split('RUNS/')[-1]

    # Read the options and arguments from the terminal
    options = "h"
    long_options = ["folder=","ki=","kx=","kvalue=","kdelta=","xrange=","ymin=","ymax=","save","nosave","SI","noSI","help"]
    opts_args = get_bashArguments(options, long_options, print_help, function, default_arguments)

    # Asign the options and arguments to variables
    for opt, arg in opts_args:
        if opt in ("-h", "--help"):
            print_help(function, default_arguments)
            sys.exit()
        elif opt == '--folder':
            arguments['folder'] = str(arg)
        elif opt == '--ki':
            arguments['ki'] = str(arg)
        elif opt == '--kx':
            arguments['kx'] = str(arg)
        elif opt == '--kvalue':
            arguments['k_value'] = str(arg)
        elif opt == '--kdelta':
            arguments['k_delta'] = str(arg)
        elif opt == '--xrange':
            arguments['x_range'] = [0, int(arg)]
        elif opt == '--ymin':
            arguments['ymin'] = float(arg)
        elif opt == '--ymax':
            arguments['ymax'] = float(arg)
        elif opt == '--SI':
            arguments['SI'] = True
        elif opt == '--noSI':
            arguments['SI'] = False
        elif opt == '--save':
            arguments['save'] = True
        elif opt == '--nosave':
            arguments['save'] = False

    # Show the parameters
    print_functionInformationOnCommandPrompt(function, default_arguments, arguments)
    return arguments
    
# Execute the python script
main_program()





    




