#!/usr/bin/python3
import numpy as np

#===============================================================================
#                              GET SCANNED PARAMETERS                          #
#===============================================================================

def get_parameters(research_or_experiment, key, knob):
    
    # Initiate the simulations and the parameters
    simulations = []
    parameters = [] 
    
    # Get all the simulations
    if hasattr(research_or_experiment, "experiments"):
        for experiment in research_or_experiment.experiments:
            simulations += experiment.simulations
    if hasattr(research_or_experiment, "simulations"):
        simulations += research_or_experiment.simulations
        
    # Get the requested parameter
    for simulation in simulations:
        simulation.input.inputParameters["parameters"]["teti"] = 1/simulation.input.inputParameters["parameters"]["tite"]
        if key!="explicit_option" and key!="nfield_periods" and key!="boundary_option":
            value = float(simulation.input.inputParameters[knob][key])
        if key=="nfield_periods": 
            if simulation.input.vmec:     value = simulation.input.nfield_periods
            if simulation.input.miller:   value = simulation.input.nperiod  
        if key=="poloidal_turns": 
            if simulation.input.vmec:     value = float(simulation.input.inputParameters[knob][key])
            if simulation.input.miller:   value = 2*(simulation.input.nperiod-1)+1 
        if key=="explicit_option": 
            if (simulation.input.inputParameters[knob][key]=="rk2"): value = 2 
            if (simulation.input.inputParameters[knob][key]=="rk3"): value = 3
        if key=="boundary_option": 
            value = len(simulation.input.inputParameters[knob][key])
        parameters.append(value)
     
    # Sort the simulations by this parameter
    sorted_indexes = list(np.array(parameters).argsort(kind="heapsort"))
    simulations = [simulations[i] for i in sorted_indexes]
    parameters  = [parameters[i] for i in sorted_indexes]
    return parameters, simulations



