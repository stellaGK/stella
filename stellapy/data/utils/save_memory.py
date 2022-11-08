
import numpy as np

def save_memory(memory, quant, array):
    
    # Calculate the size
    size=1; shape = np.shape(array)
    for i in range(len(shape)):
        size = size*shape[i]
    
    # Save the name and size
    memory['names'].append(quant)
    memory['sizes'].append(size)
    
    # Return this data
    return memory