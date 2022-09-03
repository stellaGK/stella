
import numpy as np

def print_memory(memory):
    
    # Find the total size
    total = np.sum(memory['sizes'])  
    print()
    print("     ===================================================")
    print("                          MEMORY USAGE")
    print("     ====================================================")
    print('    ', '{0:20}'.format("Quantity"), '{0:15}'.format("           Size"), '{0:15}'.format("    Percentage"))
    print('    ', '{0:20}'.format("--------"), '{0:15}'.format("           ----"), '{0:15}'.format("    ----------"))
    print()
    for i in range(len(memory['names'])): 
        percentage = round(memory['sizes'][i]/total*100,2)
        print('    ', '{0:20}'.format(memory['names'][i]), '{0:15}'.format(memory['sizes'][i]), '{0:15}'.format(percentage))
    
    