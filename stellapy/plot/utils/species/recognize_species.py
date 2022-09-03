from stellapy.utils.decorators.printv import printv
 
def recognize_species(research, specie_id): 
    
    # Assume we haven't identified the species
    specie = None
    
    # Iterate through the experiments and find the species 
    for experiment in research.experiments: 
        for simulation in experiment.simulations:
            
            # Get the mass and charge of the reference species
            z1 = simulation.inputParameters['species_parameters_1']['z']
            m1 = simulation.inputParameters['species_parameters_1']['mass']
            
            # Get the mass and charge of the second species
            if 'species_parameters_2' in simulation.inputParameters.keys():
                z2 = simulation.inputParameters['species_parameters_2']['z']
                m2 = simulation.inputParameters['species_parameters_2']['mass']
            else: z2 = None; m2 = None
            
            # Get the mass and charge of the species
            knob = 'species_parameters_'+str(int(specie_id+1))
            if knob in simulation.inputParameters.keys():
                z = simulation.inputParameters[knob]['z']
                m = simulation.inputParameters[knob]['mass']
                
                # Identify the specie
                specie_temp = identify_specie(z,m,z1,m1,z2,m2,specie_id)
                
                # Overwrite the first species
                if specie==None: 
                    specie = specie_temp
                elif specie!=None:
                    if specie!=specie_temp:
                        print('WARNING: Creating a label for the axis when we have different species!')
    
    # Standard label
    s = "_s"
    
    # Based on the specie, return the index
    if specie=="ions":      s = "i"
    if specie=="electrons": s = "e" 
    if specie=="impurity":  s = "z" 
    if specie=="carbon 6":  s = "{\\text{C}^6}" 
    if specie=="iron 16":   s = "{\\text{Fe}^{16}}" 
    return s


def identify_specie(z,m,z1,m1,z2,m2,specie_id):
    ''' m = 0.000543867  for electrons with adiabatic hydrogen
        m = 0.0002719335 for electrons with adiabatic deuterium '''

    # Assume we haven't found the species
    specie = None
    
    # Identify the ions
    if z==1:
        if int(specie_id)==0:
            if z2!=None:
                if z2==-1:
                    if 0.00026 < m2 < 0.00028: specie = 'ions'; printv("We're simulating Deuterium and kinetic electrons")  
                    if 0.00053 < m2 < 0.00055: specie = 'ions'; printv("We're simulating Hydrogen and kinetic electrons")   
            if z2==None:
                specie = 'ions'; print("We're simulating ions with adiabatic electrons")  
    
    # Identify the electrons      
    elif z==-1:
        if int(specie_id)==0:
            if 0.00026 < m < 0.00028: specie = 'electrons'; printv("We're simulating electrons with adiabatic Deuterium")
            if 0.00053 < m < 0.00055: specie = 'electrons'; printv("We're simulating electrons with adiabatic Hydrogen")   
        if int(specie_id)==1:
            if 0.00026 < m < 0.00028: specie = 'electrons'; printv("We're simulating Deuterium and kinetic electrons")   
            if 0.00053 < m < 0.00055: specie = 'electrons'; printv("We're simulating Hydrogen and kinetic electrons")
            
    # Identify the impurities
    elif z>z or m>1: 
        specie = 'impurity'; printv("We've added impurities to the simulation.")
        if z==6: specie = 'carbon 6'
        if z==16: specie = 'iron 16'
            
    # If we didn't identify the species, give a notification so it can be added
    if specie==None:
        print("WARNING: the species could not be identified automatically in stellapy/simulations/utils/get_species.py")
        
    # Return the species
    return specie
            
            
            
