
import numpy as np

def calculate_volumeAverage(self, quantity, axis_t=0, axis_z=1, axis_kx=2, axis_ky=3):
    """ Returns the volume average of the quantity and return |quantity|**2.
    The quantity has dimensions (t,tubes,z,x,y,re/im) or (t,tubes,z,kx,ky,re/im).
    Remove the dimensions (tubes,z,x,y) by taking the volume average. Which is 
    equivalent to summing over the (tubes,x,y) dimensions and taking the field
    line average to remove the z-dimension. When doing volume averages, note the following:
        int dxdy g(x,y)^2 = sum_ky |g(ky=0,kx)|^2 + 2 * sum_{kx,ky} |g(ky>0,kx)|^2
    Therefore, a factor of 2 is accounted for by means of extra_factor(iky). 
    
    Assume the dimensions are those of phi(t,z,kx,ky)."""

    # Get the extra factor to account for the double presence of ky modes. 
    extra_factor = 2*np.ones((self.dim_ky))
    if (self.vec_ky[0]==0): extra_factor[0] = 1.0
    
    # Get the extra factor "l_over_B" (=dz/B=dz*J) in the field line average 
    dl_over_B = self.dl_over_B
    
    # Determine the number of dimensions  
    dimensions = list(np.shape(quantity)) 

    # The time axis can be removed already to save computation time
    if len(dimensions)==3:
        
        # Redefine the basic axis
        if (axis_z==1 and axis_kx==2 and axis_ky==3):
            axis_z  = 0
            axis_kx = 1
            axis_ky = 2

        # Initiate the volume average
        quantity_volumeAveraged = 0
    
        # Iterate over (z,kx,ky) to remove these dimensions 
        for ikx in range(self.dim_kx):
            for iky in range(self.dim_ky):
                factor = extra_factor[iky]*dl_over_B[:]
                quantity_squared = np.real(quantity[:,ikx,iky]*np.conj(quantity[:,ikx,iky]))
                quantity_volumeAveraged = quantity_volumeAveraged + np.nansum(quantity_squared*factor)
                
    # The time axis can be removed already to save computation time
    if len(dimensions)==4:

        # Initiate the volume average as a function of time
        quantity_volumeAveraged = np.zeros((dimensions[axis_t]))
    
        # Iterate over (z,kx,ky) for each time to remove these dimensions 
        for it in range(dimensions[axis_t]):
            for ikx in range(self.dim_kx):
                for iky in range(self.dim_ky):
                    factor = extra_factor[iky]*dl_over_B[:]  
                    quantity_squared = np.real(quantity[it,:,ikx,iky]*np.conj(quantity[it,:,ikx,iky]))
                    quantity_volumeAveraged[it] += np.nansum(quantity_squared*factor)

    return quantity_volumeAveraged

