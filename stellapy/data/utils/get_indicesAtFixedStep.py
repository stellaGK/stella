 
import numpy as np
 
def get_indicesAtFixedStep(vector, dx, xmax=None):
    '''From an array only keep the values [O, dx, 2*dx, 3*dx, ...]. '''
    
    # We could ask for the last time step
    if dx==-1:
        return [-1]
    
    # Only reduce the array if <dx> is defined
    if dx!=None:
        
        # Initiate the indices and start looking at x=0
        indices = []
        ideal_x = 0 
        
        # Don't go past xmax
        if not xmax: xmax = np.max(vector)
        
        # Iterate through the vector
        for i in range(len(vector)):
            
            # When the value in the vector is bigger than the ideal x, 
            # save the index of this value and increase the ideal x by dx 
            if (vector[i] >= ideal_x) and (vector[i]<=xmax):
                indices.append(i)
                ideal_x += dx   
                
        # At the last time value if it is missing
        if i not in indices:
            indices.append(i) 
            
        
    # If no step size is defined, return all indices
    else: indices = range(len(vector))
        
    # Return the indices for the array
    return indices





