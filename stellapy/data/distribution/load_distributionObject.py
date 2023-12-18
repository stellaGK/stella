""" 

#===============================================================================
#                       Create the <distribution> object                       #
#===============================================================================   

Attach the data related to the distribution function "g" to a <distribution> object.

Attributes
----------
    g2_vs_tsz = simulation.distribution.g2_vs_tsz.g2
    g2_vs_tsmu = simulation.distribution.g2_vs_tsmu.g2
    g2_vs_tsvpa = simulation.distribution.g2_vs_tsvpa.g2 
    g2_vs_tsvpaz = simulation.distribution.g2_vs_tsvpaz.g2 
    g2_vs_tsmuvpa = simulation.distribution.g2_vs_tsmuvpa.g2 
    vec_time = simulation.distribution.g2_vs_tsz.t
    vec_z = simulation.distribution.g2_vs_tsz.z 
    
For linear simulations, simulating one mode (kx,ky) per simulation, we also have,
    g2_vs_tskxky = simulation.distribution.g2_vs_tskxky.g2
    g2_vs_szkxky = simulation.distribution.g2_vs_szkxky.g2
    g2_vs_smukxky = simulation.distribution.g2_vs_smukxky.g2
    g2_vs_svpakxky = simulation.distribution.g2_vs_svpakxky.g2
        
Average distribution squared in real space
------------------------------------------
In stella the Fourier components hat{g}(kx,ky) are related to the convential
Fourier components through hat{g}(kx,ky) = hat{G}(kx,ky)/NxNy. In other words,
the inverse Fourier transformation is defined as
    
    g_{xy} = 1/NxNy sum_{kx,ky} hat{G}(kx,ky) exp(i*kx*x+i*ky*y)
           = sum_{kx,ky} hat{g}(kx,ky) exp(i*kx*x+i*ky*y)
             
As a result, Parseval's theorem is given by

    sum_{xy} |g_{xy}|^2 = 1/NxNy sum_{kx,ky} |hat{G}(kx,ky)|^2
                        = NxNy sum_{kx,ky} |hat{g}(kx,ky)|^2  
                          
Therefore, the average distribution squared in real space, is simply given by the sum
over the (kx,ky) Fourier components as calculated by stella. 

    |g|^2 = 1/NxNy sum_{xy} |phi(x,y)|^2 
          = 1/(NxNy)^2 sum_{kxky} |hat{G}(kx,ky)|^2
          = sum_{kxky} |hat{g}(kx,ky)|^2
          
Taking into account the mirror condition, the calculation in stella is,
    |g(t,s,z,mu,vpa)|^2  = 1/(NxNy)^2 [ sum_{kx} |hat{G}(ky=0,kx)|^2 + 2 * sum_{kx,ky} |hat{G}(ky>0,kx)|^2 ]
                         = sum_{kx} |hat{g}(ky=0,kx)|^2 + 2 * sum_{kx,ky} |hat{g}(ky>0,kx)|^2
            
Thanks to the "calculate_attributeWhenReadFirstTime" wrapper, the <distribution> 
attributes will not be loaded/read/calculated until they have been requested.
    
Hanne Thienpondt
20/01/2023

"""
 
#!/usr/bin/python3   
import pathlib
import os, sys 
 
# Stellapy package  
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)    
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime 
from stellapy.data.distribution.read_distributionVsMuOrVpaOrZ import get_distributionVsMuOrVpaOrZ
from stellapy.data.distribution.read_distributionVsTimeLinear import get_distributionVsTimeLinear
from stellapy.data.distribution.read_distributionVsTime import get_distributionVsTime
from stellapy.data.distribution.read_distribution3D import get_distribution3D
from stellapy.data.distribution.read_distribution4D import get_distribution4D

#===============================================================================
#                       Create the <distribution> object                       #
#=============================================================================== 
  
class Distribution:
    
    # Copy the data from <simulation> that is needed to construct <distribution>
    def __init__(self, simulation):
        
        # Remember the file paths and the progress
        self.path = simulation.path 
        self.Progress = simulation.Progress
        
        # Remember the dimensions and vectors
        self.nakxnaky = simulation.nakxnaky
        self.dim = simulation.dim 
        self.vec = simulation.vec    
        return
    
    # Get the 2D distribution data
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_ts(self):          get_distributionVsTime(self);          return self.g2_vs_ts
    
    # Get the 3D distribution data
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_tsz(self):         get_distribution3D(self);              return self.g2_vs_tsz
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_tsmu(self):        get_distribution3D(self);              return self.g2_vs_tsmu 
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_tsvpa(self):       get_distribution3D(self);              return self.g2_vs_tsvpa  
    
    # Get the 4D distribution data
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_tsvpaz(self):      get_distribution4D(self);              return self.g2_vs_tsvpaz
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_tsmuvpa(self):     get_distribution4D(self);              return self.g2_vs_tsmuvpa  
    
    #=========================
    #   LINEAR SIMULATIONS   #
    #=========================
    
    # Get the distribution versus time 
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_tskxky(self):      get_distributionVsTimeLinear(self);    return self.g2_vs_tskxky 
    
    # Get the distribution along (vpa, mu) at the final time step for each species 
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_szkxky(self):      get_distributionVsMuOrVpaOrZ(self);    return self.g2_vs_szkxky
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_smukxky(self):     get_distributionVsMuOrVpaOrZ(self);    return self.g2_vs_smukxky
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_svpakxky(self):    get_distributionVsMuOrVpaOrZ(self);    return self.g2_vs_svpakxky
    
#-------------------------
def load_distributionObject(self): 
    self.distribution = Distribution(self) 
