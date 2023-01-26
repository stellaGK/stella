
import numpy as np
import scipy.fft as fftpack 

def calculate_inverseFourierTransform(self, quantity_kxky, axis_kx=2, axis_ky=3):
    """ Return the discrete inverse Fast Fourier transform of the quantity_kxky in (kx,ky)."""

    # Determine the number of dimensions and the dimensions over which to iterate 
    dimensions = list(np.shape(quantity_kxky)) 
    iterate_dimensions = list(range(len(dimensions)))
    
    # The (nakx,naky) dimensions become (nx,ny)
    dimensions[axis_kx] = self.nx
    dimensions[axis_ky] = self.ny  
    
    # Don't iterate over (kx,ky)
    if axis_ky < axis_kx:
        iterate_dimensions.pop(axis_kx)
        iterate_dimensions.pop(axis_ky)
    if axis_kx < axis_ky:
        iterate_dimensions.pop(axis_ky)
        iterate_dimensions.pop(axis_kx)
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # And we have already put all quantities in ascending order with respect to kx.
    # However for the Fourier transform we want again the ordering of vec_kx_stella.
    sorted_indexes = [ self.vec.kx.index(kx) for kx in self.vec.kx_stella]
    
    # Change the ordering of the kx-axis
    if len(dimensions)==2: 
        if axis_kx==0: quantity_kxky = quantity_kxky[sorted_indexes, :] 
        if axis_kx==1: quantity_kxky = quantity_kxky[:, sorted_indexes] 
    if len(dimensions)==3: 
        if axis_kx==0: quantity_kxky = quantity_kxky[sorted_indexes, :, :] 
        if axis_kx==1: quantity_kxky = quantity_kxky[:, sorted_indexes, :] 
        if axis_kx==2: quantity_kxky = quantity_kxky[:, :, sorted_indexes] 
    if len(dimensions)==4: 
        if axis_kx==0: quantity_kxky = quantity_kxky[sorted_indexes, :, :, :] 
        if axis_kx==1: quantity_kxky = quantity_kxky[:, sorted_indexes, :, :] 
        if axis_kx==2: quantity_kxky = quantity_kxky[:, :, sorted_indexes, :] 
        if axis_kx==3: quantity_kxky = quantity_kxky[:, :, :, sorted_indexes] 
    if len(dimensions)==5: 
        if axis_kx==0: quantity_kxky = quantity_kxky[sorted_indexes, :, :, :, :] 
        if axis_kx==1: quantity_kxky = quantity_kxky[:, sorted_indexes, :, :, :] 
        if axis_kx==2: quantity_kxky = quantity_kxky[:, :, sorted_indexes, :, :] 
        if axis_kx==3: quantity_kxky = quantity_kxky[:, :, :, sorted_indexes, :] 
        if axis_kx==4: quantity_kxky = quantity_kxky[:, :, :, :, sorted_indexes]  

    # Get the number of modes in real space and Fourier space
    nx = self.nx;  nakx = self.nakx
    ny = self.ny;  naky = self.naky
    
    # Initiate the quantity in real space (x,y)
    quantity_xy = np.ones(dimensions)*np.NaN
    
    # Calculate the inverse Fourier in two dimensions (kx,ky)
    if len(dimensions)==2:
        
        # Make sure we have quantity_kxky(kx,ky)
        if axis_kx==0: pass
        if axis_kx==1: quantity_kxky = np.swapaxes(quantity_kxky, 0, 1)
        
        # Return the inverse Fourier
        quantity_xy = calculate_IFFT(self, quantity_kxky, nx, ny, nakx, naky)
        
        # Make sure we return quantity_xy with the original axis order
        if axis_kx==0: pass
        if axis_kx==1: quantity_xy = np.swapaxes(quantity_xy, 0, 1)
        
    # Calculate the inverse Fourier in three dimensions (i,kx,ky)
    if len(dimensions)==3:
        
        # Make sure we have quantity_kxky(i,kx,ky)
        if axis_kx==0: quantity_kxky = np.swapaxes(quantity_kxky, 0, 1)
        if axis_kx==1: pass
        if axis_kx==2: quantity_kxky = np.swapaxes(quantity_kxky, 2, 1)
        if axis_ky==0: quantity_kxky = np.swapaxes(quantity_kxky, 0, 2)
        if axis_ky==1: quantity_kxky = np.swapaxes(quantity_kxky, 1, 2)
        if axis_ky==2: pass
        
        # Iterate over the extra dimension and calculate the inverse Fourier
        for i in range(dimensions[iterate_dimensions[0]]): 
            quantity_xy[i,:,:] = calculate_IFFT(self, quantity_kxky[i,:,:], nx, ny, nakx, naky)

        # Make sure we return quantity_xy with the original axis order
        if axis_kx==0: quantity_xy = np.swapaxes(quantity_xy, 0, 1)
        if axis_kx==1: pass
        if axis_kx==2: quantity_xy = np.swapaxes(quantity_xy, 2, 1)
        if axis_ky==0: quantity_xy = np.swapaxes(quantity_xy, 0, 2)
        if axis_ky==1: quantity_xy = np.swapaxes(quantity_xy, 1, 2)
        if axis_ky==2: pass
        
    # Calculate the inverse Fourier in four dimensions (i,j,kx,ky)
    if len(dimensions)==4:
        
        # Make sure we have quantity_kxky(i,j,kx,ky)
        if axis_kx==0: quantity_kxky = np.swapaxes(quantity_kxky, 0, 2)
        if axis_kx==1: quantity_kxky = np.swapaxes(quantity_kxky, 1, 2)
        if axis_kx==2: pass
        if axis_kx==3: quantity_kxky = np.swapaxes(quantity_kxky, 3, 2)
        if axis_ky==0: quantity_kxky = np.swapaxes(quantity_kxky, 0, 3)
        if axis_ky==1: quantity_kxky = np.swapaxes(quantity_kxky, 1, 3)
        if axis_ky==2: quantity_kxky = np.swapaxes(quantity_kxky, 2, 3)
        if axis_ky==3: pass
        
        # Iterate over the extra dimensions and calculate the inverse Fourier
        for i in range(dimensions[iterate_dimensions[0]]):
            for j in range(dimensions[iterate_dimensions[1]]):
                quantity_xy[i,j,:,:] = calculate_IFFT(self, quantity_kxky[i,j,:,:], nx, ny, nakx, naky)

        # Make sure we return quantity_xy with the original axis order
        if axis_kx==0: quantity_xy = np.swapaxes(quantity_kxky, 0, 2)
        if axis_kx==1: quantity_xy = np.swapaxes(quantity_kxky, 1, 2)
        if axis_kx==2: pass
        if axis_kx==3: quantity_xy = np.swapaxes(quantity_kxky, 3, 2)
        if axis_ky==0: quantity_xy = np.swapaxes(quantity_kxky, 0, 3)
        if axis_ky==1: quantity_xy = np.swapaxes(quantity_kxky, 1, 3)
        if axis_ky==2: quantity_xy = np.swapaxes(quantity_kxky, 2, 3)
        if axis_ky==3: pass
        

    # Return the quantity in real space (x,y)
    return quantity_xy

def calculate_IFFT(self, quantity_kxky, nx, ny, nakx, naky):

    # In stella we have less modes along kx than along x so add zero padding in 
    # the middle (= highest and lowest wavenumbers) to increase the resolution
    ikx_max  = int(nakx/2+1)                      
    ipad_upx = int(ikx_max + nx - nakx)-1    
    quantity_kxky_temp = np.zeros((nx, ny),dtype = 'complex_') 
    quantity_kxky_temp[:ikx_max, :naky] = quantity_kxky[:ikx_max,:naky]         
    quantity_kxky_temp[ipad_upx+1:, :naky] = quantity_kxky[ikx_max:, :naky] 
        
    # In stella we only have the positive modes along ky, therefore, add their
    # complex conjugates to the end of the matrix, the middle is again zero
    # padded since we have less modes in Fourier space than in real space
    # Thanks to the complex cunjugates, the imaginary part of the ifft is zero.
    quantity_kxky_temp[1:, ny-naky+1:] = np.flip(np.flip(np.conj(quantity_kxky_temp[1:, 1:naky]), axis=1), axis=0)  
    quantity_kxky_temp[0, ny-naky+1:]  = np.flip(np.conj(quantity_kxky_temp[0, 1:naky]), axis=0)
      
    # Increase the dimensions of quantity(kx,ky) from (nakx,naky) to (nx, ny)
    quantity_kxky = quantity_kxky_temp
          
    # There is some ordered noise on the imaginary part of the order of E-13
    quantity_kxky = np.real(quantity_kxky.round(decimals=9)) + 1j * np.imag(quantity_kxky.round(decimals=9))
         
    # Take the inverse Fourier transform along (kx,ky) of the given quantity_kxky
    quantity_xy = np.real(fftpack.ifft2(quantity_kxky))
    
    # Check parsevals theorem 
    #print("------ PLANCHERAL -----------")  
    #print("SUM IN FOURIER:  ", np.sum(np.abs(quantity_kxky)**2.0))
    #print("SUM IN REAL:     ", np.sum(np.abs(quantity_xy)**2.0)*nx*ny)
    
    # Return the quantity in real space (x,y)
    return quantity_xy

