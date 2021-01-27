#!/usr/bin/python3
''' Run read_wout() as a bash command on the command prompt.'''

# Load the modules
import sys, os

# Tell python where to find the personal modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)).split("stellapy/")[0])
from stellapy.utils.decorators import print_decoratorOnCommandPrompt 
from stellapy.simulations.Simulation import create_simulations 
from stellapy.config import turnOnVerboseWrapper_configurationFile, turnOffVerboseWrapper_configurationFile

# Get the python function we want to execute here
from stellapy.data import read_referenceUnits

# If the program is run from the terminal, then collect the arguments and execute write_profiles()
def main_program():
    if __name__ == "__main__":
        
        # Create simulation objects
        simulations = create_simulations(folders=os.getcwd())
        
        # Check the inputs
        for self in simulations:
            ref, input_values = read_referenceUnits(self.inputParameters, self.ref_a, self.ref_B, self.prof_n, self.prof_T, verbose=True)
        
        # Now print the parameters
        turnOnVerboseWrapper_configurationFile()
        check_referenceUnits(ref, input_values)
        
def check_referenceUnits(ref, input_values):

    from stellapy.utils.decorators import print_inTwoColumns 

    print()
    print_inTwoColumns("Input parameters:", "Reference units:")
    print_inTwoColumns("rho=r/a:",         input_values['rho'], "",              "rho_i:", ref['rho_i']*100, "cm")    
    print_inTwoColumns("s1_temp:",         input_values['s1_tempI'], "keV",       "v_th:", ref['vthermal'], "m/s")  
    print_inTwoColumns("s1_dens:",         input_values['s1_densI'], "10^19 m^-3",       "rho_i:", ref['rho_i']*100, "cm")    
    print_inTwoColumns("s1_temp_prof:",    input_values['s1_temp'], "keV",       "v_th:", ref['vthermal'], "m/s")  
    print_inTwoColumns("s1_dens_prof:",    input_values['s1_dens'], "10^19 m^-3",       "gyro_w:",  ref['omega'], "1/s")  
    print_inTwoColumns("s1_mass:",         ref['mass']*10**27, "10^27 kg",       "gyrobohm_p:",  ref['flux_part'], "") 
    print_inTwoColumns("B:",               ref['B'], "T",                        "gyrobohm_v:",  ref['flux_mom'], "") 
    print_inTwoColumns("a:",               ref['length'], "m",                   "gyrobohm_q:",  ref['flux_heat'], "") 
    print_inTwoColumns("Z:",               ref['Z'], " ", " ", " ", " ") 
    print()
    
# Execute the python script
turnOffVerboseWrapper_configurationFile()
print_decoratorOnCommandPrompt("begin")
main_program()
print_decoratorOnCommandPrompt("end")

