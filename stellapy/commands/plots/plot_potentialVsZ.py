#!/usr/bin/python3
''' Run plot_potentialVsZ() as a bash command on the command prompt.'''

# Load the modules
import sys, os

# Tell python where to find the personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.bashscripts.get_bashArguments import get_bashArguments
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.utils.decorators import print_functionInformationOnCommandPrompt

# Get the python function we want to execute here
from stellapy.plot.potential import plot_potentialVsZ

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        function = "plot.potential.plot_potentialVsZ"
        default_arguments = {\
            # Specify which simulations to plot
                'folders' : 'MANDETORY ARGUMENT', \
                'input_files' : None,\
            # Other options
                'convergence_studies' : False,\
                'verbose' : True, \
            # Specify data range
                'y_quantity' : "phi_real",\
                'z_quantity' : "z",\
                'units' : "zeta",\
                'kx' : -0.0,\
                'ky' : None,\
                'kmin' : 0,\
                'kmax' : 100,\
                'plotted_modes' : "unstable",\
                'y_range_real' : None,\
                'y_range_imag' : None,\
                'y_range_squared' : None,\
                'z_range' : None,\
                'converged_modes' : None,\
                'unconverged_modes' : None,\
            # For the GUI the figure object and axes already exist
                'show_figure' : True,\
                'figure' : None,\
            # Appearance options
                'fontsize' : 20,\
                'handlelength' : 1}
        arguments = get_arguments(function, default_arguments)
        plot_potentialVsZ(**arguments)
        
# Explain the program
def print_help(function, default_arguments):
    print_decoratorOnCommandPrompt("begin")
    print_functionInformationOnCommandPrompt(function, default_arguments, None)
    print("\n The optional arguments are:")
    print("      --folder      <folder>" )
    print("      --kx          <plot at this kx>    {-0.0}")
    print("      --kmin        <plot ky > kmin>     {0.0}")
    print("      --kmax        <plot ky < kmax>     {100}")
    print("      --ymax        <max of y_range>     {auto}")
    print("      --zman        <max of z_range>     {auto}")
    print_decoratorOnCommandPrompt("end")

# Get the arguments from the terminal
def get_arguments(function, default_arguments):
    
    # Make a new dictionary which will hold the new arguments
    arguments = default_arguments.copy()

    # Initiate mandetory variables
    arguments['folders'] = os.getcwd().split('RUNS/')[-1]

    # Read the options and arguments from the terminal
    options = "h"
    long_options = ["folder=","kx=","kmin=","kmax=","zmax=","ymax=","stable","movie","help"]
    opts_args = get_bashArguments(options, long_options, print_help, function, default_arguments)

    # Asign the options and arguments to variables
    for opt, arg in opts_args:
        if opt in ("-h", "--help"):
            print_help(function, default_arguments)
            sys.exit()
        elif opt == '--folder':
            arguments['folders'] = str(arg)
        elif opt == '--kx':
            arguments['kx'] = str(arg)
        elif opt == '--kmax':
            arguments['kmax'] = str(arg)
        elif opt == '--kmin':
            arguments['kmin'] = str(arg)
        elif opt == '--zmax':
            arguments['z_range'] = [0,1]
            arguments['z_range'][1] = float(arg)

    # Show the parameters
    print_functionInformationOnCommandPrompt(function, default_arguments, arguments)
    return arguments
    
# Execute the python script
main_program()





    




