#!/usr/bin/python3
import sys, os
import numpy as np

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep) 
from stellapy.simulations.Simulation import create_simulations
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                            CHECK INPUT VARIABLES                             #
#===============================================================================

def check_inputs(): 
    ''' Print the most important input parameters of the stella simulation to the command prompt. '''
    
    # Create simulation objects
    simulations = create_simulations(folders=os.getcwd()) 
    
    # Get an input_file in folder 
    for simulation in simulations:
        
        # Get the input parameters
        inputParameters = simulation.input.inputParameters
        nonlinear = "True" if inputParameters["physics_flags"]["nonlinear"] else "False"
        
        # Gather the data
        indent = "    "; length = 70
        d1 = "{0:<15}".format("nonlinear:")+"{0:<15}".format(nonlinear)
        d2 = "{0:<15}".format("rho:")+"{0:<15}".format(inputParameters["vmec_parameters"]["rho"])
        d3 = "{0:<15}".format("nspec:")+"{0:<15}".format(inputParameters['species_knobs']['nspec'])
        d4 = "{0:<15}".format("tite:")+"{0:<15}".format(inputParameters['parameters']['tite'])
        d5 = "{0:<15}".format("ref B:")+"{0:<15}".format(simulation.geometry.bref)
        d6 = "{0:<15}".format("ref a:")+"{0:<15}".format(simulation.geometry.aref)
        s1 = "{0:<15}".format("s1_dens:")+"{0:<15}".format(inputParameters['species_parameters_1']['dens'])
        s2 = "{0:<15}".format("s1_temp:")+"{0:<15}".format(inputParameters['species_parameters_1']['temp'])
        s3 = "{0:<15}".format("s1_fprim:")+"{0:<15}".format(inputParameters['species_parameters_1']['fprim'])
        s4 = "{0:<15}".format("s1_tprim:")+"{0:<15}".format(inputParameters['species_parameters_1']['tprim'])
        s5 = "{0:<15}".format("sgn(b0):")+"{0:<15}".format(simulation.geometry.sign_B)
        s6 = "{0:<15}".format("iota:")+"{0:<15}".format(np.round(simulation.geometry.iota,3))

        # Print the relevant data
        print()
        print(" "+"".center(length,"-")) 
        print('    Most important input parameters for the simulation:')
        print('       ', simulation.path.input_file.parent.parent.name+"/"+simulation.path.input_file.parent.name+"/"+simulation.path.input_file.name, '\n')
        print(indent+d1+s1)     
        print(indent+d2+s2)     
        print(indent+d3+s3)     
        print(indent+d4+s4)     
        print(indent+d5+s5)     
        print(indent+d6+s6)     
        print(" "+"".center(length,"-")) 
        print()

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    bash = Bash(check_inputs)   
    args = bash.get_arguments()
    del args['folder']
    check_inputs(**args)   
    




