#!/usr/bin/python3
''' Run check_inputs() as a bash command on the command prompt.'''

# Tell python where to find the personal modules
import sys, os

# Tell python where to find the personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.simulations.Simulation import create_simulations 
from stellapy.utils.decorators import print_decoratorOnCommandPrompt
from stellapy.config import turnOffVerboseWrapper_configurationFile

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        
        # Create simulation objects
        simulations = create_simulations(folders=os.getcwd())
        
        # Check the inputs
        check_inputs(simulations)


def check_inputs(simulations): 
    ''' Print the most important input parameters of the stella simulation to the command prompt. '''

    # Load the modules
    from stellapy.utils.decorators.print_inTwoColumns import print_inTwoColumns3 as print_inTwoColumns
    from stellapy.utils.decorators.verbose_wrapper import indent

    # Get an input_file in folder 
    for simulation in simulations:
        
        # Get the input parameters
        inputParameters = simulation.inputParameters

        # Print the relevant data
        print()
        print(indent, 'Most important input parameters for the simulation:', simulation.id)
        print_inTwoColumns("non-linear:",simulation.inputParameters["physics_flags"]["nonlinear"], \
                           "s1_dens:", inputParameters['species_parameters_1']['dens'])  
        print_inTwoColumns("rho:",inputParameters["vmec_parameters"]["rho"], \
                           "s1_temp:", inputParameters['species_parameters_1']['temp'])   
        print_inTwoColumns("nspec:", inputParameters['species_knobs']['nspec'], \
                           "s1_fprim:", inputParameters['species_parameters_1']['fprim'])   
        print_inTwoColumns("tite:",inputParameters['parameters']['tite'], \
                           "s1_tprim:", inputParameters['species_parameters_1']['tprim'])   
        print_inTwoColumns("ref B:",simulation.ref_B, \
                           "sgn(b0):", simulation.sign_B)   
        print_inTwoColumns("ref a:",simulation.ref_a, \
                           "iota:", simulation.iota)   
        print()

# Execute the python script
turnOffVerboseWrapper_configurationFile()
print_decoratorOnCommandPrompt("begin")
main_program()
print_decoratorOnCommandPrompt("end") 




    




