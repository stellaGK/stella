 
import numpy as np
import configparser
import os, h5py, pathlib
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
        self.y0 = inputParameters["kt_grids_box_parameters"]["y0"]
        self.svalue = inputParameters["vmec_parameters"]["torflux"]
        self.vmec_filename = inputParameters["vmec_parameters"]["vmec_filename"] 
        self.nfield_periods = inputParameters["vmec_parameters"]["nfield_periods"]    
        if self.vmec_filename=='wout*.nc': self.nfield_periods = np.nan
        if self.vmec_filename=='wout*.nc': self.svalue = inputParameters["millergeo_parameters"]["rhoc"]*inputParameters["millergeo_parameters"]["rhoc"]
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
    
    # Basic geometry arrays in *.geometry
    @calculate_attributeWhenReadFirstTime 
    def zeta(self):         get_geometryData(self);         return self.zeta
    @calculate_attributeWhenReadFirstTime 
    def gbdrift(self):      get_geometryData(self);         return self.gbdrift
    
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
    
    # Add the geometry for a simulation
    if self.object=="Simulation":
        self.geometry = Geometry(self)
        return
        
    # For a linear simulation, add the geometry per mode. Since most geometries are
    # identical, only read each unique geometry once (self==mode)
    if self.object=="Mode":
        load_onlyReferenceGeometries(self)
        return

def load_geometryObjectDirectly(self):  
    self.geometry = Geometry(self)  
     
#===============================================================================
#                          LOAD REFERENCE GEOMETRIES                           #
#===============================================================================

def load_onlyReferenceGeometries(self):
    """ Load a geometry object for each unique geometry. """
    
    # Initialize
    reference_geometries = {}
    
    # Only read one geometry per unique geometry
    for imode, mode in enumerate(self.simulation.modes): 
        
        # If the geometry is already read, create a reference to the existing geometry
        if mode.path.geometry in self.simulation.path.loaded_geometries: 
            mode.geometry = self.simulation.modes[reference_geometries[mode.path.geometry]].geometry
            mode.geometry.path = mode.path
            
        # If the geometry isn't read, read it now
        if mode.path.geometry not in self.simulation.path.loaded_geometries:
            self.simulation.path.loaded_geometries.append(mode.path.geometry) 
            reference_geometries[mode.path.geometry] = imode
            mode.geometry = Geometry(mode)  
    return  

#-----------------------------------
def read_uniqueGeometries(path):  

    # Make sure we have the overview of unique geometries 
    path_of_reference_geometries = path.folder / (path.name + ".list.geometries.ini")
    if not os.path.isfile(path_of_reference_geometries):
        from stellapy.data.geometry.write_h5FileForGeometry import write_h5FileForGeometry 
        write_h5FileForGeometry(path.folder)
        
    # Read the overview of unique geometries
    list_of_reference_geometries = configparser.ConfigParser()
    list_of_reference_geometries.read(path_of_reference_geometries)
    
    # Get the identifier of each reference geometry: "Geometry X"
    ids_of_reference_geometries = sorted([ g for g in list_of_reference_geometries.keys() if "Geometry " in g]) 
    return list_of_reference_geometries, ids_of_reference_geometries

#-----------------------------------
def find_referenceGeometry(path, mode, list_of_reference_geometries, ids_of_reference_geometries):
    
    # Initiate 
    reference_geometry_path = None
    
    # If it is written on marconi, translate to local paths
    for i in ids_of_reference_geometries:
        for key in list_of_reference_geometries[i].keys():
            if 'marconi' in list_of_reference_geometries[i][key]:   
                list_of_reference_geometries[i][key] = str(path.folder)+"/"+(list_of_reference_geometries[i][key].split("/"+path.folder.name+"/")[-1])
            elif '/mnt/lustre/' in list_of_reference_geometries[i][key]:   
                list_of_reference_geometries[i][key] = str(path.folder)+"/"+(list_of_reference_geometries[i][key].split("/"+path.folder.name+"/")[-1])
            else:
                break  
            
    # Iterate through the reference geometies to find <input_file>
    for reference_geometry_id in ids_of_reference_geometries:    
        if str(mode.input_file) in list_of_reference_geometries[reference_geometry_id].values():   
            reference_geometry_path = path.folder / (path.name + ".unique.geometry" + reference_geometry_id.split("Geometry ")[-1])
            break
        
    # Return the reference geometry and the path to the file
    if reference_geometry_path==None:
        print("\n      Couldn't find the reference geometry for:")
        print("             ", mode.input_file, "\n")  
        import sys; sys.exit()
    return reference_geometry_path
    
#===============================================================================
#                           SAVE THE GEOMETRY OBJECT                           #
#===============================================================================
# Replace the "*.geometry", "*.vmec_geo" and "wout*.nc" files with a "*.geo.h5" file.

def save_geometryObject(simulation):
    
    # Load all the geometry data  
    simulation.geometry = Geometry(simulation)  
    simulation.geometry.load_data()
    
    # Save the relevant input parameters
    for knob in ["geo_knobs", "vmec_parameters", "zgrid_parameters", "millergeo_parameters", "vpamu_grids_parameters"]: 
        setattr(simulation.geometry, knob, type('knob', (object,), {})) 
        for key, value in simulation.input.inputParameters[knob].items():
            setattr(getattr(simulation.geometry, knob), key, value) 

    # Check whether we succeeded in loading the data
    if np.isnan(simulation.geometry.jacob[0]):
        print("    ---> The geometry data could not be found for:")
        print("             ", simulation.input_file)
        import sys; sys.exit()
    
    # Save the geometry data to the "*.geo.h5" file  
    with h5py.File(simulation.path.geometry, 'w') as f:
        for key, value in vars(simulation.geometry).items():  
            if key not in ["path", "geo_knobs", "vmec_parameters", "zgrid_parameters", "millergeo_parameters", "vpamu_grids_parameters"]:  
                f.create_dataset(key, data=value) 
            elif key in ["geo_knobs", "vmec_parameters", "zgrid_parameters", "millergeo_parameters", "vpamu_grids_parameters"]:  
                knob = key; knob_group = f.create_group(knob) 
                for key, value in vars(getattr(simulation.geometry, knob)).items(): 
                    if "__" not in key:
                        if value==None: value = "None"
                        knob_group.create_dataset(key, data=value) 
                       
    # Track progress
    print("    ----> Saved the geometry file as " + simulation.path.geometry.name)  
    return 

#----------------------------
def read_geometryObject(input_file):
    geometry_file = input_file.with_suffix(".geo.h5") 
    geometry = Geometry(type('obj', (object,), {'input_files' : [input_file]}) ) 
    with h5py.File(geometry_file, 'r') as f:
        for key in f.keys():
            if key not in ["geo_knobs", "vmec_parameters", "zgrid_parameters", "millergeo_parameters", "vpamu_grids_parameters"]:  
                setattr(geometry, key, f[key].value) 
            elif key in ["geo_knobs", "vmec_parameters", "zgrid_parameters", "millergeo_parameters", "vpamu_grids_parameters"]:  
                knob = key
                for key in getattr(f, knob).keys(): 
                    value = f[key].value if (f[key].value!="None") else None
                    setattr(getattr(geometry, knob), key, value) 
    
################################################################################
#                     USE THESE FUNCTIONS AS A MAIN SCRIPT                     #
################################################################################
if __name__ == "__main__":  
    
    from stellapy.simulations.Simulation import create_simulations
    import copy, timeit; start = timeit.timeit()

    folder = pathlib.Path("/home/hanne/CIEMAT/PREVIOUSRUNS/LINEARMAPS/W7Xstandard_rho0.7_aLTe0/ResolutionScan/fprim4tprim4_ky1.5/nmu") 
    folder = pathlib.Path("/home/hanne/CIEMAT/PREVIOUSRUNS/LINEARMAPS/W7Xstandard_rho0.7_aLTe0/LinearMap/fprim4tprim4")      
    folder = pathlib.Path("/home/hanne/CIEMAT/PREVIOUSRUNS/LINEARMAPS/W7Xstandard_rho0.7_aLTe0/ResolutionScan/fprim4tprim4_ky1.5/nzed") 
    print("CREATE SIMULATIONS")
    simulations = create_simulations(folders=folder, input_files=None, ignore_resolution=True, number_variedVariables=5)
    print("We have "+str(len(simulations))+" simulations. Test whether we have") 
    print("references to the data instead of copies.")
    for simulation in simulations:
        print("\nSimulation:", simulation.id)  
        b0_test = copy.deepcopy(simulation.modes[0].geometry.b0) 
        zeta_test = copy.deepcopy(simulation.modes[0].geometry.zeta) 
        jacob_test = copy.deepcopy(simulation.modes[0].geometry.jacob) 
        jtwist_test = copy.deepcopy(simulation.modes[0].geometry.jtwist)  
        for mode in simulation.modes: 
            b0 = mode.geometry.b0
            zeta = mode.geometry.zeta
            jacob = mode.geometry.jacob 
            jtwist = mode.geometry.jtwist 
            reference_test1f = b0 is b0_test
            reference_test2f = zeta is zeta_test
            reference_test3f = jacob is jacob_test
            reference_test4f = jtwist is jtwist_test
            reference_test1 = b0 is simulation.modes[0].geometry.b0
            reference_test2 = zeta is simulation.modes[0].geometry.zeta
            reference_test3 = jacob is simulation.modes[0].geometry.jacob
            reference_test4 = jtwist is simulation.modes[0].geometry.jtwist 
            name = "("+str(mode.kx)+", "+str(mode.ky)+")" 
            print("{:<15}".format(name), len(mode.geometry.zeta), mode.path.geometry.name, " References:", reference_test1, reference_test2, reference_test3, reference_test4,\
                                             "   TEST:", reference_test1f, reference_test2f, reference_test3f, reference_test4f) 

        
        