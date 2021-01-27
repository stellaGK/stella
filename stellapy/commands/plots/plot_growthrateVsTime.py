#!/usr/bin/python3
''' Run plot_growthrateVsTime() as a bash command on the command prompt.'''

# Load the modules
import sys, os

# Tell python where to find the personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.bashscripts.get_bashArguments import get_bashArguments
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.utils.decorators import print_functionInformationOnCommandPrompt

# Get the python function we want to execute here
from stellapy.plot.frequencyAndGrowthrate import plot_growthrateVsTime

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        function = "plot.frequencyAndGrowthrate.plot_growthrateVsTime"
        default_arguments = {\
              # Specify which simulations to plot
                'case_folders' : None,\
                'input_files' : None,\
            # Other options
                'convergence_studies' : False,\
                'verbose' : True, \
            # Specify data range
                'kx' : -0.0,
                'kmin' : 0,\
                'kmax' : 100,\
                'plotted_modes' : "unstable",\
                'y_range_omega' : None,\
                'y_range_gamma' : None,\
                'x_range' : None,\
            # For the GUI the figure object and axes already exist
                'show_figure' : True,\
                'figure' : None,\
            # Appearance options
                'fontsize' : 20,\
                'handlelength' : 1}
        arguments = get_arguments(function, default_arguments)
        plot_growthrateVsTime(**arguments)

# Explain the program
def print_help(function, default_arguments):
    print_decoratorOnCommandPrompt("begin")
    print_functionInformationOnCommandPrompt(function, default_arguments, None)
    print("\n The optional arguments are:")
    print("      --folder      <case folder>" )
    print("      --kx          <plot at this kx>    {-0.0}")
    print("      --kmin        <plot ky > kmin>     {0.0}")
    print("      --kmax        <plot ky < kmax>     {100}")
    print("      --tmin        <min range time>     {auto}")
    print("      --tmax        <max range time>     {auto}")
    print("      --omax        <max range omega>    {[0,1]}")
    print("      --gmax        <max range gamma>    {[0,0.2]}")
    print("\n The options are:")
    print("      --stable       {False, True}       Plot the stable modes for phi2(t_last) < 100")
    print_decoratorOnCommandPrompt("end")

# Get the arguments from the terminal
def get_arguments(function, default_arguments):
    
    # Make a new dictionary which will hold the new arguments
    arguments = default_arguments.copy()

    # Initiate mandetory variables
    arguments['folder'] = os.getcwd().split('RUNS/')[-1]

    # Read the options and arguments from the terminal
    options = "h"
    long_options = ["folder=","kx=","kmin=","kmax=","tmin=","tmax=","omax=","gmax=","stable","help","i="]
    opts_args = get_bashArguments(options, long_options, print_help, function, default_arguments)

    # Asign the options and arguments to variables
    for opt, arg in opts_args:
        if opt in ("-h", "--help"):
            print_help(function, default_arguments)
            sys.exit()
        elif opt == '--folder':
            arguments['folder'] = str(arg)
        elif opt == '--kx':
            arguments['kx'] = str(arg)
        elif opt == '--kmax':
            arguments['kmax'] = float(arg)
        elif opt == '--kmin':
            arguments['kmin'] = float(arg)
        elif opt == '--tmin':
            arguments['x_range'] = [0,800]
            arguments['x_range'][0] = float(arg)
        elif opt == '--tmax':
            if arguments['x_range'] is None:
                arguments['x_range'] = [0,800]
            arguments['x_range'][1] = float(arg)
        elif opt == '--omax':
            arguments['y_range_omega'] = [0,1]
            arguments['y_range_omega'][1] = float(arg)
        elif opt == '--gmax':
            arguments['y_range_gamma'] = [0,1]
            arguments['y_range_gamma'][1] = float(arg)
        elif opt == '--stable':
            arguments['plotted_modes'] = 'stable'
        elif opt == '--i':
            arguments['input_files'] = str(arg)

    # Show the parameters
    print_functionInformationOnCommandPrompt(function, default_arguments, arguments)
    return arguments
    
# Execute the python script
main_program()





    




