 
import os
import numpy as np
from stellapy.data.geometry.read_wout import get_woutData, get_woutDataExtra
from stellapy.data.geometry.read_output import get_outputData
from stellapy.data.geometry.read_geometry import get_geometryData
from stellapy.data.geometry.read_geometry import get_geometryDataExtra
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime
from stellapy.data.input.read_inputFile import read_inputFile
 
#===============================================================================
#                        CREATE THE GEOMETRY OBJECT
#===============================================================================
# Combine the data from the *.geometry, *.vmec_geo, *.out.nc and the VMEC file,
# together with some calculations to save all geometry related quantities.
#===============================================================================

class Geometry:
    
    # Initialize the geometry data
    def __init__(self, simulation): 

        # Remember the paths to the data files
        self.path = simulation.path 
        
        # Remember the quantities needed to calculate [iota, shat, ...]
        if os.path.isfile(self.path.input): inputParameters = simulation.input.inputParameters
        if not os.path.isfile(self.path.input): inputParameters = read_inputFile(simulation.input_file)
        self.y0 = inputParameters['kxky_grid_box']['y0']
        self.svalue = inputParameters['geometry_vmec']['torflux']
        self.vmec_filename = inputParameters['geometry_vmec']['vmec_filename']
        self.nfield_periods = inputParameters['geometry_vmec']['nfield_periods']
        if self.vmec_filename=='wout*.nc': self.nfield_periods = np.nan
        if self.vmec_filename=='wout*.nc': self.svalue = inputParameters['geometry_miller']['rhoc']*inputParameters['geometry_miller']['rhoc']
        return
    
    # Load all the geometry data
    def load_data(self):
        get_woutData(self)
        get_geometryData(self)
        get_outputData(self)
        get_woutDataExtra(self)
        get_geometryDataExtra(self)
        
    # Basic geometry scalars in *.geometry
    @calculate_attributeWhenReadFirstTime
    def aref(self):         get_geometryData(self);          return self.aref
    @calculate_attributeWhenReadFirstTime
    def bref(self):         get_geometryData(self);          return self.bref
    @calculate_attributeWhenReadFirstTime
    def qinp(self):         get_geometryData(self);          return self.qinp
    @calculate_attributeWhenReadFirstTime
    def rhotor(self):       get_geometryData(self);          return self.rhotor
    
    # Extra geometry scalars in *.geometry
    @calculate_attributeWhenReadFirstTime
    def alpha(self):        get_geometryDataExtra(self);    return self.alpha
    @calculate_attributeWhenReadFirstTime
    def zed(self):          get_geometryDataExtra(self);    return self.zed
    @calculate_attributeWhenReadFirstTime
    def gradpar(self):      get_geometryDataExtra(self);    return self.gradpar
    @calculate_attributeWhenReadFirstTime
    def gds2(self):         get_geometryDataExtra(self);    return self.gds2
    @calculate_attributeWhenReadFirstTime
    def gds21(self):        get_geometryDataExtra(self);    return self.gds21
    @calculate_attributeWhenReadFirstTime
    def gds22(self):        get_geometryDataExtra(self);    return self.gds22
    @calculate_attributeWhenReadFirstTime
    def gds23(self):        get_geometryDataExtra(self);    return self.gds23
    @calculate_attributeWhenReadFirstTime
    def gds24(self):        get_geometryDataExtra(self);    return self.gds24
    @calculate_attributeWhenReadFirstTime
    def cvdrift(self):      get_geometryDataExtra(self);    return self.cvdrift
    @calculate_attributeWhenReadFirstTime
    def gbdrift0(self):     get_geometryDataExtra(self);    return self.gbdrift0
    @calculate_attributeWhenReadFirstTime
    def bmag_psi0(self):    get_geometryDataExtra(self);    return self.bmag_psi0
    
    # Basic data in the VMEC file
    @calculate_attributeWhenReadFirstTime
    def nfp(self):          get_woutData(self);             return self.nfp
    @calculate_attributeWhenReadFirstTime
    def sign_B(self):       get_woutData(self);             return self.sign_B
      
    # Extra data in the VMEC file
    @calculate_attributeWhenReadFirstTime
    def b0(self):           get_woutDataExtra(self);        return self.b0
    @calculate_attributeWhenReadFirstTime
    def aminor(self):       get_woutDataExtra(self);        return self.aminor
    @calculate_attributeWhenReadFirstTime
    def rmajor(self):       get_woutDataExtra(self);        return self.rmajor
    @calculate_attributeWhenReadFirstTime
    def volume(self):       get_woutDataExtra(self);        return self.volume
    @calculate_attributeWhenReadFirstTime
    def vmec_path(self):    get_woutDataExtra(self);        return self.vmec_path
    
    # Variables used in stella that can be derived from the VMEC file
    @calculate_attributeWhenReadFirstTime
    def P(self):            get_woutDataExtra(self);        return self.P
    @calculate_attributeWhenReadFirstTime
    def Lx(self):           get_woutDataExtra(self);        return self.Lx
    @calculate_attributeWhenReadFirstTime
    def Ly(self):           get_woutDataExtra(self);        return self.Ly
    @calculate_attributeWhenReadFirstTime
    def dkx(self):          get_woutDataExtra(self);        return self.dkx
    @calculate_attributeWhenReadFirstTime
    def dky(self):          get_woutDataExtra(self);        return self.dky
    @calculate_attributeWhenReadFirstTime
    def shat(self):         get_woutDataExtra(self);        return self.shat
    @calculate_attributeWhenReadFirstTime
    def iota(self):         get_woutDataExtra(self);        return self.iota
    @calculate_attributeWhenReadFirstTime
    def jtwist(self):       get_woutDataExtra(self);        return self.jtwist
    @calculate_attributeWhenReadFirstTime
    def diotaOverds(self):  get_woutDataExtra(self);        return self.diotaOverds
    @calculate_attributeWhenReadFirstTime
    def twist_and_shift_geo_fac(self):  get_woutDataExtra(self); return self.twist_and_shift_geo_fac
    
    # Variables related to the geometry that are saved in the output file
    @calculate_attributeWhenReadFirstTime
    def bmag(self):         get_outputData(self);           return self.bmag
    @calculate_attributeWhenReadFirstTime
    def jacob(self):        get_outputData(self);           return self.jacob
    @calculate_attributeWhenReadFirstTime
    def vec_z(self):        get_outputData(self);           return self.vec_z
    @calculate_attributeWhenReadFirstTime
    def dl_over_B(self):    get_outputData(self);           return self.dl_over_B
    @calculate_attributeWhenReadFirstTime
    def mu_weights(self):   get_outputData(self);           return self.mu_weights
    @calculate_attributeWhenReadFirstTime
    def vpa_weights(self):  get_outputData(self);           return self.vpa_weights

#===============================================================================
#                           LOAD THE GEOMETRY OBJECT                           #
#===============================================================================

def load_geometryObject(self):
    self.geometry = Geometry(self)
    return


