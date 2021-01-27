#!/usr/bin/python3
''' Run reduce_sizeNetcdf() as a bash command on the command prompt.'''

# Tell python where to find the personal modules
import sys, os

# Tell python where to find the personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.commands.utils.get_bashArguments import get_bashArguments
from stellapy.config import turnOffVerboseWrapper_configurationFile, turnOnVerboseWrapper_configurationFile

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        state = get_arguments() 
        if state == "on":
            turnOnVerboseWrapper_configurationFile()
        if state == "off":
            turnOffVerboseWrapper_configurationFile()
            
# Get the arguments from the terminal
def get_arguments():
    
    # Read the options and arguments from the terminal
    options = "of"
    long_options = ["on","off"]
    opts_args = get_bashArguments(options, long_options, None, None, None)

    # Asign the options and arguments to variables
    for opt in opts_args[0]:
        if opt == '--on':
            state = 'on'
        elif opt == '--off':
            state = 'off'

    return state

# Execute the python script
main_program()






    




