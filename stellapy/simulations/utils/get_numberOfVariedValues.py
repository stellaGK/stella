
import numpy as np

def get_numberOfVariedValues(experiments):
    ''' Get the total number of varied values, to determine whether we need to use markers
    # to differentiate between different simulations or whether the line colors are sufficient. '''
    number_ofVariedValues = 0
    for experiment in experiments: 
        number_ofVariedValues = np.max([number_ofVariedValues,len(experiment.variedValues)])
    return number_ofVariedValues