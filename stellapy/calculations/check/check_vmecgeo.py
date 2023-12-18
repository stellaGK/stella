#!/usr/bin/python3
''' Run read_vmecgeo() as a bash command on the command prompt.'''

# Load the modules
import sys, os

# Tell python where to find the personal modules
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.utils.decorators import print_decoratorOnCommandPrompt 
from stellapy.simulations.Simulation import create_simulations 
from stellapy.data import read_vmecgeo
from stellapy.utils.config import turnOnVerboseWrapper_configurationFile, turnOffVerboseWrapper_configurationFile


# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        
        # Create simulation objects
        simulations = create_simulations(folders=os.getcwd())
        
        # Check the inputs
        for self in simulations:
            for input_file in self.input_files:
                geo_data = read_vmecgeo(input_file)
                turnOnVerboseWrapper_configurationFile()
                check_wout(geo_data)
                turnOffVerboseWrapper_configurationFile()

def check_wout(geo_data):
    from stellapy.utils.decorators.verbose import indent
    print("\n", indent, 'Check the parameters in the ".vmec_geo" file.')
    for key, value in geo_data.items():
        print(indent, '    ', '{:10}'.format(key+":"), " ", value)

# Execute the python script
turnOffVerboseWrapper_configurationFile()
print_decoratorOnCommandPrompt("begin")
main_program()
print_decoratorOnCommandPrompt("end") 

