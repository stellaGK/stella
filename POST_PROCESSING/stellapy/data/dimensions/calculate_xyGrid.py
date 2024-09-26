""" 

#===============================================================================
#                     Calculate the real space grid (x,y)                      #
#===============================================================================

Based on the Fourier grid (kx,ky) we can define the real space grid (x,y).
    Lx = 2*pi/dkx        dkx = kx[1]-kx[0]      dx = Lx/(nx-1)
    Ly = 2*pi/dky        dky = ky[1]-ky[0]      dy = Ly/(ny-1)

Hanne Thienpondt
20/01/2023

"""

#!/usr/bin/python3
import os, sys
import pathlib
import numpy as np 

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_nx
from stellapy.data.geometry.calculate_gridDivisionsAndSize import calculate_ny

#===============================================================================
#                     Calculate the real space grid (x,y)                      #
#===============================================================================

def calculate_xyGrid(vec): 
    
    # Read the Fourier grid divisions
    nakx = len(vec.kx)
    naky = len(vec.ky)
    
    # Get the Fourier grid step size
    dkx = vec.kx[1] - vec.kx[0]
    dky = vec.ky[1] - vec.ky[0]
    
    # Calculate the real space grid divisions
    nx = calculate_nx(nakx)  
    ny = calculate_ny(naky)   
    
    # Calculate the real space box size
    Lx = 2*np.pi/dkx
    Ly = 2*np.pi/dky
    
    # Calculate the real space grid
    vec_x = np.linspace(-Lx/2,Lx/2,nx) 
    vec_y = np.linspace(-Ly/2,Ly/2,ny) 
    return nx, ny, vec_x, vec_y

#===============================================================================
#             Attach the real space grid (x,y) to the <vec> object             #
#===============================================================================

def get_xyGrid(self):
    
    # Calculate the real space grid
    nx, ny, vec_x, vec_y = calculate_xyGrid(self.vec)
    
    # Attach the real space grid
    self.dim.x = nx
    self.dim.y = ny 
    self.vec.x = vec_x 
    self.vec.y = vec_y 
    return 
    
    