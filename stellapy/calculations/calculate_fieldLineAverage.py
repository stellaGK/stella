
import numpy as np
 
def calculate_fieldLineAverage(self, quantity, z_axis=0, z_specific=None, i_kx=None, i_ky=None):
    """Calculates the field line average over a given quantity along z, 
    which removes the 'z' dimension. This is achieved by multiplying each
    z-component of the quantity with dz/B or dz*Jacobian and summing this.
    Note that the factor dl_over_B needs to go with |phi|**2."""
    
    # Get the extra factor "dl_over_B" (=dz/B=dz*J) in the field line average 
    if self.nonlinear==True:  dl_over_B = self.dl_over_B 
    if self.nonlinear==False: dl_over_B = self.dl_over_B[:,i_kx,i_ky] 

    # Get the field line average along z
    if not z_specific: 
        
        # Determine the number of dimensions and the dimensions over which to iterate
        dimensions = list(np.shape(quantity)) 
        iterate_dimensions = list(range(len(dimensions)))
        iterate_dimensions.pop(z_axis) 
        
        # Field line average over one dimension (which is z)
        if len(iterate_dimensions)==0: 
            quantity_fieldLineAveraged = 0 
            quantity_fieldLineAveraged = np.nansum(quantity[:]*dl_over_B[:])   
            
        # Field line average over two dimensions (z + one extra dimension)
        if len(iterate_dimensions)==1: 
            quantity_fieldLineAveraged = np.zeros((dimensions[iterate_dimensions[0]]), dtype=complex)
            for i in range(dimensions[iterate_dimensions[0]]): 
                if z_axis==0: quantity_fieldLineAveraged[i] = np.nansum(quantity[:,i]*dl_over_B[:]) 
                if z_axis==1: quantity_fieldLineAveraged[i] = np.nansum(quantity[i,:]*dl_over_B[:])   
                
        # Field line average over three dimensions (z + two extra dimension)
        if len(iterate_dimensions)==2: 
            dim0 = dimensions[iterate_dimensions[0]]
            dim1 = dimensions[iterate_dimensions[1]]
            quantity_fieldLineAveraged = np.zeros((dim0, dim1), dtype=complex)
            for i in range(dim0): 
                for j in range(dim1): 
                    if z_axis==0: quantity_fieldLineAveraged[i,j] = np.nansum(quantity[:,i,j]*dl_over_B[:]) 
                    if z_axis==1: quantity_fieldLineAveraged[i,j] = np.nansum(quantity[i,:,j]*dl_over_B[:])   
                    if z_axis==2: quantity_fieldLineAveraged[i,j] = np.nansum(quantity[i,j,:]*dl_over_B[:])   
                    
        # Field line average over four dimensions (z + three extra dimensions)
        if len(iterate_dimensions)==3: 
            dim0 = dimensions[iterate_dimensions[0]]
            dim1 = dimensions[iterate_dimensions[1]]
            dim2 = dimensions[iterate_dimensions[2]]
            quantity_fieldLineAveraged = np.zeros((dim0, dim1, dim2), dtype=complex)
            for i in range(dim0): 
                for j in range(dim1): 
                    for k in range(dim2):  
                        if z_axis==0: quantity_fieldLineAveraged[i,j,k] = np.nansum(quantity[:,i,j,k]*dl_over_B[:]) 
                        if z_axis==1: quantity_fieldLineAveraged[i,j,k] = np.nansum(quantity[i,:,j,k]*dl_over_B[:])   
                        if z_axis==2: quantity_fieldLineAveraged[i,j,k] = np.nansum(quantity[i,j,:,k]*dl_over_B[:])
                        if z_axis==3: quantity_fieldLineAveraged[i,j,k] = np.nansum(quantity[i,j,k,:]*dl_over_B[:])
                        
        # Field line average over five dimension (z + four extra dimensions)
        if len(iterate_dimensions)==4: 
            dim0 = dimensions[iterate_dimensions[0]]
            dim1 = dimensions[iterate_dimensions[1]]
            dim2 = dimensions[iterate_dimensions[2]]
            dim3 = dimensions[iterate_dimensions[3]]
            quantity_fieldLineAveraged = np.zeros((dim0, dim1, dim2, dim3), dtype=complex)
            for i in range(dim0): 
                for j in range(dim1): 
                    for k in range(dim2): 
                        for l in range(dim3): 
                            if z_axis==0: quantity_fieldLineAveraged[i,j,k,l] = np.nansum(quantity[:,i,j,k,l]*dl_over_B[:]) 
                            if z_axis==1: quantity_fieldLineAveraged[i,j,k,l] = np.nansum(quantity[i,:,j,k,l]*dl_over_B[:])   
                            if z_axis==2: quantity_fieldLineAveraged[i,j,k,l] = np.nansum(quantity[i,j,:,k,l]*dl_over_B[:])
                            if z_axis==3: quantity_fieldLineAveraged[i,j,k,l] = np.nansum(quantity[i,j,k,:,l]*dl_over_B[:])
                            if z_axis==4: quantity_fieldLineAveraged[i,j,k,l] = np.nansum(quantity[i,j,k,l,:]*dl_over_B[:])
          
        # Return the field line average of the quantity
        return quantity_fieldLineAveraged

    # Get the quantity at a specific z
    if z_specific:
                        
        # Get the index along the time axis where we encounter <t_specific>
        z_index = np.argwhere(self.vec_z > z_specific)[0][0]

        # Quantity at a specific time when there is one dimension (which is t)
        if len(list(np.shape(quantity)))==1:  
            quantity_atSpecificZ = quantity[z_index]
            
        # Quantity at a specific time when there are two dimensions (t + one extra dimension)
        if len(list(np.shape(quantity)))==2:   
            if z_axis==0: quantity_atSpecificZ = quantity[z_index,:]
            if z_axis==1: quantity_atSpecificZ = quantity[:,z_index]
                
        # Quantity at a specific time when there are three dimensions (t + two extra dimension)
        if len(list(np.shape(quantity)))==3:  
            if z_axis==0: quantity_atSpecificZ = quantity[z_index,:,:]
            if z_axis==1: quantity_atSpecificZ = quantity[:,z_index,:]
            if z_axis==2: quantity_atSpecificZ = quantity[:,:,z_index]
                
        # Quantity at a specific time when there are four dimensions (t + three extra dimension)
        if len(list(np.shape(quantity)))==4:  
            if z_axis==0: quantity_atSpecificZ = quantity[z_index,:,:,:]
            if z_axis==1: quantity_atSpecificZ = quantity[:,z_index,:,:]
            if z_axis==2: quantity_atSpecificZ = quantity[:,:,z_index,:]
            if z_axis==3: quantity_atSpecificZ = quantity[:,:,:,z_index]
             
        # Quantity at a specific time when there are five dimensions (t + four extra dimension)
        if len(list(np.shape(quantity)))==5:     
            if z_axis==0: quantity_atSpecificZ = quantity[z_index,:,:,:,:]
            if z_axis==1: quantity_atSpecificZ = quantity[:,z_index,:,:,:]
            if z_axis==2: quantity_atSpecificZ = quantity[:,:,z_index,:,:]
            if z_axis==3: quantity_atSpecificZ = quantity[:,:,:,z_index,:]
            if z_axis==4: quantity_atSpecificZ = quantity[:,:,:,:,z_index]
 
        # Return the field line average of the quantity
        return quantity_atSpecificZ
    
