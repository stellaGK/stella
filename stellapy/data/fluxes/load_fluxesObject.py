#!/usr/bin/python3  
import sys, os   
import pathlib 

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime
from stellapy.data.fluxes.read_quasilinearfluxes import get_quasilinearfluxes
from stellapy.data.fluxes.read_fluxesFile import get_fluxes 
from stellapy.data.fluxes.read_fluxes3D import get_fluxes3D
from stellapy.data.fluxes.read_fluxes4D import get_fluxes4D
 
#===============================================================================
#                        CREATE THE FLUXES OBJECT
#===============================================================================

class Fluxes:
    
    # Copy the data from <simulation> that is needed to construct <fluxes>
    def __init__(self, simulation):
        
        # Remember the file paths and the progress
        self.path = simulation.path 
        self.Progress = simulation.Progress
        self.simulation = simulation  
        
        # Remember the dimensions and vectors
        self.dim = simulation.dim 
        self.vec = simulation.vec      
        
        # Remember the sign of the magnetic field
        self.sign_B = simulation.geometry.sign_B 
        
        # Toggles
        self.includeFluxNorm = False
        return

    # Read the fluxes data
    @calculate_attributeWhenReadFirstTime 
    def date(self):                 get_fluxes(self);      return self.date
    @calculate_attributeWhenReadFirstTime 
    def dim_time(self):             get_fluxes(self);      return self.dim_time 
    @calculate_attributeWhenReadFirstTime 
    def pflux_vs_ts(self):          get_fluxes(self);      return self.pflux_vs_ts
    @calculate_attributeWhenReadFirstTime 
    def qflux_vs_ts(self):          get_fluxes(self);      return self.qflux_vs_ts
    @calculate_attributeWhenReadFirstTime  
    def vflux_vs_ts(self):          get_fluxes(self);      return self.vflux_vs_ts  
    
    # Get the 3D fluxes data 
    @calculate_attributeWhenReadFirstTime 
    def qflux_vs_tskx(self):        get_fluxes3D(self);          return self.qflux_vs_tskx
    @calculate_attributeWhenReadFirstTime 
    def pflux_vs_tskx(self):        get_fluxes3D(self);          return self.pflux_vs_tskx
    @calculate_attributeWhenReadFirstTime 
    def vflux_vs_tskx(self):        get_fluxes3D(self);          return self.vflux_vs_tskx
    @calculate_attributeWhenReadFirstTime 
    def qflux_vs_tsky(self):        get_fluxes3D(self);          return self.qflux_vs_tsky
    @calculate_attributeWhenReadFirstTime 
    def pflux_vs_tsky(self):        get_fluxes3D(self);          return self.pflux_vs_tsky
    @calculate_attributeWhenReadFirstTime 
    def vflux_vs_tsky(self):        get_fluxes3D(self);          return self.vflux_vs_tsky
    @calculate_attributeWhenReadFirstTime 
    def qflux_vs_tsz(self):         get_fluxes3D(self);          return self.qflux_vs_tsz
    @calculate_attributeWhenReadFirstTime 
    def pflux_vs_tsz(self):         get_fluxes3D(self);          return self.pflux_vs_tsz
    @calculate_attributeWhenReadFirstTime 
    def vflux_vs_tsz(self):         get_fluxes3D(self);          return self.vflux_vs_tsz
    
    # Get the 4D fluxes data 
    @calculate_attributeWhenReadFirstTime 
    def qflux_vs_tskxky(self):      get_fluxes4D(self);          return self.qflux_vs_tskxky
    @calculate_attributeWhenReadFirstTime 
    def pflux_vs_tskxky(self):      get_fluxes4D(self);          return self.pflux_vs_tskxky
    @calculate_attributeWhenReadFirstTime 
    def vflux_vs_tskxky(self):      get_fluxes4D(self);          return self.vflux_vs_tskxky
    
    # Get the 4D linear fluxes data (at last time step)
    @calculate_attributeWhenReadFirstTime 
    def phi2_vs_kxky(self):         get_quasilinearfluxes(self);    return self.phi2_vs_kxky 
    @calculate_attributeWhenReadFirstTime 
    def qflux_vs_skxky(self):       get_quasilinearfluxes(self);    return self.qflux_vs_skxky 
    @calculate_attributeWhenReadFirstTime 
    def pflux_vs_skxky(self):       get_quasilinearfluxes(self);    return self.pflux_vs_skxky 
    @calculate_attributeWhenReadFirstTime 
    def vflux_vs_skxky(self):       get_quasilinearfluxes(self);    return self.vflux_vs_skxky 
    
#-------------------------
def load_fluxesObject(self):
    self.fluxes = Fluxes(self)
    return
