
#===============================================================================
#                              GET UNIQUE PARAMETERS                           #
#===============================================================================

def get_uniqueParameters(research, key, knob): 
    
    # Intitiate the parameters
    parameters = [] 
     
    # Get the requested parameter
    for experiment in research.experiments: 
        for simulation in experiment.simulations:
            if key!="explicit_option" and key!="nfield_periods":
                value = float(simulation.input.inputParameters[knob][key])
            if key=="nfield_periods": 
                if simulation.vmec:     value = simulation.input.nfield_periods
                if simulation.miller:   value = simulation.input.nperiod 
            if key=="explicit_option": 
                if (simulation.input.inputParameters[knob][key]=="rk2"): value = 2 
                if (simulation.input.inputParameters[knob][key]=="rk3"): value = 3
            parameters.append(value)
     
    # Sort the simulations by this parameter
    parameters = sorted(list(set(parameters)))
    return parameters
