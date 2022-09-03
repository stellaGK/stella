 
import numpy as np 
import os, pathlib, h5py
from scipy.io import netcdf as scnetcdf   
from stellapy.utils.config import CONFIG
from stellapy.utils.files.get_filesInFolder import get_filesInFolder 
from stellapy.data.input.load_inputObject import read_uniqueInputs, find_referenceInput
from stellapy.data.geometry.load_geometryObject import read_uniqueGeometries, find_referenceGeometry
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime  
from stellapy.data.input.read_inputFile import read_svalueFromInputFile

#===============================================================================
#                           CREATE THE PATHS OBJECT                            #
#===============================================================================  
#
# We want to read from the following files for nonlinear sims (<simulation>):
#
#    *.in                      input               txt 
#    *.slurm                   slurm               txt 
#    *.dimensions              slurm               h5 
#    *.rho70.geometry          dimensions          h5
# 
#    *.dt10.moments            output              h5    (upar_vs_tszkxkyri, dens_vs_tszkxkyri, temp_vs_tszkxkyri)
#    *.dt10.potential3D        potential3D         h5    (phi_vs_tz_zonal, phi_vs_tz_nozonal, phi_vs_tz, phi2_vs_tz)
#    *.dt10.potential4D        potential4D         h5    (phi2_vs_tkxky, phi_vs_tkxky)
#    *.dt10.potential5D        potential5D         h5    (phi_tzkxky)
#    *.dt10.distribution3D     distribution3D      h5    (g_vs_tsz, g_vs_tsmu, g_vs_tsvpa)
#    *.dt10.distribution4D     distribution4D      h5    (g_vs_tsmuvpa, g_vs_tsvpaz)
#    *.dt10.fluxes_vs_t        fluxes              txt   (pflux_vs_ts, ...) 
#    *.dt10.phi_vs_t           phi_vs_t            txt   (phi_vs_t, phi2_vs_t)
#
# We want to read from the following files for linear sims (<mode>):  
#
#    *.list.inputs.ini        
#    *.list.inputs.geometries        
#    *.unique.inputX.ini       input               ini
#    *.unique.geometryX        geometry            h5
#
#    *.slurm                   slurm               txt 
#    *.dimensions              slurm               h5 
#
#    *.dt10.omega_vs_t         omega_vs_t          txt    (omega_vs_t, gamma_vs_t)
#    *.lineardata              lineardata          ini    (omega_last, gamma_last)
#
#    *.phi_vs_z                phi_vs_z            txt    (phi_vs_z, phi2_vs_z)
#    *.dt10.phi_vs_t           phi_vs_t            txt    (phi_vs_t, phi2_vs_t)
#    *.dt10.dphiz_vs_t         dphiz_vs_t          txt    (dphiz_vs_t)
#    *.dt10.potential3D        potential3D         h5     (phi_vs_tz, phi2_vs_tz)
#
#    *.g_vs_z                  g_vs_z              txt    (g_vs_sz)
#    *.g_vs_mu                 g_vs_mu             txt    (g_vs_smu)
#    *.g_vs_vpa                g_vs_vpa            txt    (g_vs_svpa) 
#    *.dt10.g_vs_t             g_vs_t              txt    (g_vs_st)  
#    *.dt10.distribution3D     distribution        h5     (g_vs_tsmu, g_vs_tsvpa)
#    *.dt10.distribution4D     distribution        h5     (g_vs_tsmuvpa)
#
#
# Stella outputs the following files:
#    *.in                      input_file           
#    *.final_fields            finalphi_stella
#    *.geometry                geometry_stella
#    *.vmec.geo                vmecgeo_stella
#    *.vmec_geo                vmecgeo_stella
#    *.fluxes                  fluxes_stella 
#    *.out.nc                  output_stella 
#    wout*.nc                  vmec    
#    XXXXXX_*.out              slurmout_stella  
#    slurm_*                   slurmin_stella
#    slurm-*.out               slurmstatus_stella 
#
# Useless files
#    *species.input   
#    *.error
#    *.vmec_geo
#=============================================================================== 


class Paths:
    
    # Save the paths of all the data files corresponding to <simulation>
    def __init__(self, simulation): 
        
        # Remember the parent simulation and its input file  
        self.simulation = simulation if simulation.object=="Simulation" else simulation.simulation
        
        # Get the input file
        if simulation.object=="Simulation":
            if simulation.linear:       self.input_file = simulation.modes[0].input_file
            if simulation.nonlinear:    self.input_file = simulation.input_file 
        if simulation.object=="Mode":
            if simulation.linear:       self.input_file = simulation.input_file 

        # Remember whether we have a <simulation> or a <mode>
        self.object = simulation.object
        self.linear = simulation.linear
        self.nonlinear = simulation.nonlinear
        
        # Link itself for more logical access to the variables
        self.path = self   
        return 
    
    # Basic information of the simulation: folder and name
    @calculate_attributeWhenReadFirstTime 
    def name(self):                 get_nameSimulation(self);       return self.name 
    @calculate_attributeWhenReadFirstTime 
    def folder(self):               get_folderPath(self);           return self.folder 
        
    # Set the paths when we ask for them
    @calculate_attributeWhenReadFirstTime 
    def vmec(self):                 get_vmecPath(self);             return self.vmec 
    @calculate_attributeWhenReadFirstTime 
    def input(self):                get_inputPath(self);            return self.input 
    @calculate_attributeWhenReadFirstTime 
    def output(self):               get_outputPath(self);           return self.output 
    @calculate_attributeWhenReadFirstTime 
    def omega(self):                get_omegaPath(self);            return self.omega 
    @calculate_attributeWhenReadFirstTime 
    def fluxes(self):               get_fluxesPath(self);           return self.fluxes 
    @calculate_attributeWhenReadFirstTime 
    def geometry(self):             get_geometryPath(self);         return self.geometry 
    @calculate_attributeWhenReadFirstTime 
    def moments3D(self):            get_moments3DPath(self);        return self.moments3D 
    @calculate_attributeWhenReadFirstTime 
    def moments4D(self):            get_moments4DPath(self);        return self.moments4D 
    @calculate_attributeWhenReadFirstTime 
    def moments5D(self):            get_moments5DPath(self);        return self.moments5D 
    @calculate_attributeWhenReadFirstTime 
    def dimensions(self):           get_dimensionsPath(self);       return self.dimensions 
    @calculate_attributeWhenReadFirstTime 
    def saturated(self):            get_saturatedPath(self);        return self.saturated
    @calculate_attributeWhenReadFirstTime 
    def fluxes3D(self):             get_fluxes3DPath(self);         return self.fluxes3D
    @calculate_attributeWhenReadFirstTime 
    def fluxes4D(self):             get_fluxes4DPath(self);         return self.fluxes4D
    @calculate_attributeWhenReadFirstTime 
    def potential3D(self):          get_potential3DPath(self);      return self.potential3D
    @calculate_attributeWhenReadFirstTime 
    def potential4D(self):          get_potential4DPath(self);      return self.potential4D
    @calculate_attributeWhenReadFirstTime 
    def potential5D(self):          get_potential5DPath(self);      return self.potential5D
    @calculate_attributeWhenReadFirstTime 
    def phaseshifts(self):          get_phaseshiftsPath(self);      return self.phaseshifts 
    @calculate_attributeWhenReadFirstTime 
    def distribution3D(self):       get_distribution3DPath(self);   return self.distribution3D
    @calculate_attributeWhenReadFirstTime 
    def distribution4D(self):       get_distribution4DPath(self);   return self.distribution4D
    @calculate_attributeWhenReadFirstTime 
    def trappedWeights5D(self):     get_trappedWeights5DPath(self); return self.trappedWeights5D 
    
    # Check linear convergence to identify jumpers
    @calculate_attributeWhenReadFirstTime 
    def dphiz_vs_t(self):           get_dphizvstPath(self);         return self.dphiz_vs_t 
    
    # Data as a funtion of time
    @calculate_attributeWhenReadFirstTime 
    def g_vs_t(self):               get_gvstPath(self);             return self.g_vs_t 
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_t(self):             get_phivstPath(self);           return self.phi_vs_t 
    
    # Data at the last time step (linear)
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_z(self):             get_phivszPath(self);           return self.phi_vs_z
    @calculate_attributeWhenReadFirstTime 
    def g_vs_z(self):               get_finalgPath(self);           return self.g_vs_z 
    @calculate_attributeWhenReadFirstTime 
    def g_vs_mu(self):              get_finalgPath(self);           return self.g_vs_mu 
    @calculate_attributeWhenReadFirstTime 
    def g_vs_vpa(self):             get_finalgPath(self);           return self.g_vs_vpa
    
    # For linear simulations we remember the unique inputs and geometries 
    @calculate_attributeWhenReadFirstTime 
    def loaded_inputs(self):        get_inputPath(self);            return self.loaded_inputs 
    @calculate_attributeWhenReadFirstTime 
    def loaded_geometries(self):    get_geometryPath(self);         return self.loaded_geometries 
    
    # Remember the data in the output file
    @calculate_attributeWhenReadFirstTime 
    def output_keys(self):          get_outputKeys(self);           return self.output_keys 
    @calculate_attributeWhenReadFirstTime 
    def output_stella_keys(self):   get_outputKeys(self);           return self.output_stella_keys 
    
    # Get the stella files
    @calculate_attributeWhenReadFirstTime 
    def omega_stella(self):         get_stellaPaths(self);          return self.omega_stella 
    @calculate_attributeWhenReadFirstTime 
    def fluxes_stella(self):        get_stellaPaths(self);          return self.fluxes_stella 
    @calculate_attributeWhenReadFirstTime 
    def output_stella(self):        get_stellaPaths(self);          return self.output_stella 
    @calculate_attributeWhenReadFirstTime 
    def geometry_stella(self):      get_stellaPaths(self);          return self.geometry_stella 
    @calculate_attributeWhenReadFirstTime 
    def finalphi_stella(self):      get_stellaPaths(self);          return self.finalphi_stella 
    
#----------------------    
def load_pathObject(self): 
    self.path = Paths(self) 
    return
     
#===============================================================================
#                           CREATE DUMMY PATH OBJECT                           #
#=============================================================================== 

def create_dummyPathObject(input_file, vmec_filename, nonlinear):
    """ Needed to access all the paths of a simulation, without actually going
    through the <create_simulations> routine. Return a working <path> object."""
    
    # Create a dummy path and input object
    path_object  = type('Paths', (object,), {})
    input_object = type('input', (object,), {'vmec_filename' : vmec_filename})

    # Create a dummy mode and simulation object
    attributes   = {'input' : input_object, 'input_file' : input_file, 'linear' : not nonlinear, 'nonlinear' : nonlinear}
    mode_object  = type('Mode',  (object,), {'object' : 'Mode', 'path' : path_object, **attributes}) 
    simulation_object = type('Simulation', (object,), {'object' : 'Simulation', 'modes' : [mode_object], **attributes})
    mode_object.simulation = simulation_object 
    
    # Load the paths for the mode and simulation
    load_pathObject(simulation_object) 
    load_pathObject(simulation_object.modes[0])
    
    # Return the finished <path> object of the mode or simulation
    if not nonlinear:   return simulation_object.modes[0].path
    elif nonlinear:     return simulation_object.path

#===============================================================================
#                      DEFINE THE DATA IN THE OUTPUT FILE                      #
#=============================================================================== 

def get_outputKeys(self):
    
    # Initiate the keys
    self.path.output_keys = []
    self.path.output_stella_keys = []
    
    # Read the keys of the *.out.h5 file
    if os.path.isfile(self.path.output):
        with h5py.File(self.path.output, 'r') as h5file:
            self.path.output_keys = list(h5file.keys()) + ["dimensions"]
            
    # Read the keys of the *.out.nc file
    if os.path.isfile(self.path.output_stella):
        try:
            with scnetcdf.netcdf_file(self.path.output_stella,'r') as ncfile:
                self.path.output_stella_keys += list(ncfile.variables.keys()) + ["dimensions"]
        except:
            import netCDF4 as nc4 
            with nc4.Dataset(self.path.output_stella) as ncfile:
                self.path.output_stella_keys = list(ncfile.variables.keys()) + ["dimensions"]
    return
 
#===============================================================================
#                               DEFINE THE PATHS                               #
#=============================================================================== 

def get_folderPath(self): 
    """" Get the folder which is the parent, except when the parent is "OLD". """
    if "OLD" not in self.input_file.parent.name: self.path.folder = self.input_file.parent
    elif "OLD" in self.input_file.parent.name:   self.path.folder = self.input_file.parent.parent

#------------------------
def get_nameSimulation(self): 
    """ Get the name of the input file, which is without "kyX" for a linear simulation. """
    self.path.name = self.input_file.name.split("_ky")[0] if "_ky" in self.input_file.name else self.input_file.stem
    return

#------------------------
def get_saturatedPath(self):
    self.path.saturated = get_pathWithLargestTend(self, "out.nc.t")
    return

#------------------------
def get_outputPath(self):
    self.path.output = get_pathWithSmallestTimeStep(self, "output", ".out.h5")
    return

def get_omegaPath(self):  
    self.path.omega = get_pathWithSmallestTimeStep(self, "omega_vs_t", ".omega_reduced") 
    return

def get_fluxesPath(self): 
    self.path.fluxes = get_pathWithSmallestTimeStep(self, "fluxes_vs_t", ".fluxes_reduced")
    return

def get_fluxes3DPath(self): 
    self.path.fluxes3D = get_pathWithSmallestTimeStep(self, "fluxes3D")
    return 

def get_fluxes4DPath(self): 
    self.path.fluxes4D = get_pathWithSmallestTimeStep(self, "fluxes4D")
    return 

def get_potential3DPath(self): 
    self.path.potential3D = get_pathWithSmallestTimeStep(self, "potential3D")
    return

def get_potential4DPath(self): 
    self.path.potential4D = get_pathWithSmallestTimeStep(self, "potential4D")
    return

def get_potential5DPath(self): 
    self.path.potential5D = get_pathWithSmallestTimeStep(self, "potential5D")
    return

def get_distribution3DPath(self): 
    self.path.distribution3D = get_pathWithSmallestTimeStep(self, "distribution3D")
    return

def get_distribution4DPath(self): 
    self.path.distribution4D = get_pathWithSmallestTimeStep(self, "distribution4D")
    return

def get_trappedWeights5DPath(self): 
    self.path.trappedWeights5D = get_pathWithSmallestTimeStep(self, "trappedWeights5D") 
    return

def get_phaseshiftsPath(self): 
    self.path.phaseshifts = get_pathWithSmallestTimeStep(self, "phaseshifts")
    return

def get_moments3DPath(self): 
    self.path.moments3D = get_pathWithSmallestTimeStep(self, "moments3D")
    return

def get_moments4DPath(self): 
    self.path.moments4D = get_pathWithSmallestTimeStep(self, "moments4D")
    return

def get_moments5DPath(self): 
    self.path.moments5D = get_pathWithSmallestTimeStep(self, "moments5D")
    return

def get_dphizvstPath(self): 
    self.path.dphiz_vs_t = get_pathWithSmallestTimeStep(self, "dphiz_vs_t")
    return

def get_phivstPath(self): 
    self.path.phi_vs_t = get_pathWithSmallestTimeStep(self, "phi_vs_t")
    return

def get_gvstPath(self): 
    self.path.g_vs_t = get_pathWithSmallestTimeStep(self, "g_vs_t")
    return

#------------------------
def get_dimensionsPath(self):  
    self.path.dimensions = self.input_file.with_suffix(".dimensions")
    return

def get_phivszPath(self):  
    self.path.phi_vs_z = self.input_file.with_suffix(".phi_vs_z")
    return

def get_finalgPath(self):  
    self.path.g_vs_z = self.input_file.with_suffix(".g_vs_z")
    self.path.g_vs_mu = self.input_file.with_suffix(".g_vs_mu")
    self.path.g_vs_vpa = self.input_file.with_suffix(".g_vs_vpa")
    return
    
#------------------------
def get_stellaPaths(self): 
    """ Get the paths to the stella output files. """
    
    # Get the standard output files for stella 
    self.path.finalphi_stella = self.input_file.with_suffix(".final_fields")
    self.path.geometry_stella = self.input_file.with_suffix(".geometry")
    self.path.fluxes_stella   = self.input_file.with_suffix(".fluxes")
    self.path.omega_stella    = self.input_file.with_suffix(".omega")
    self.path.output_stella   = self.input_file.with_suffix(".out.nc")
    return

#------------------------
def get_vmecPath(self): 
    
    # Get the path of the VMEC file
    self.vmec_filename = vmec_filename = self.simulation.input.vmec_filename
    self.path.vmec = self.path.input.parent / vmec_filename if (vmec_filename != 'wout*.nc') else "pathNotFound"
    if (not os.path.isfile(self.path.vmec)) and (self.path.vmec!=None): 
        
        # Decide whether we're on the local pc or the supercomputer 
        vmec_folder = pathlib.Path(CONFIG["PATHS"]['VMEC']) 
         
        # See if the VMEC file is present in the main VMEC folder  
        if vmec_filename!='wout*.nc':   
            vmec_files = get_filesInFolder(vmec_folder, start=vmec_filename, end=".nc")
            self.path.vmec = vmec_files[0] if (vmec_files!=None) else self.path.vmec  
        if vmec_filename=='wout*.nc':    
            self.path.vmec = 'MillerDoesntHaveAVmec'  
    return

#------------------------
def get_inputPath(self):
    """" Get the input file of the <simulation> or <mode> """
    
    # Read the input parameters from the ini file
    if self.linear:     self.path.input = self.input_file 
    if self.nonlinear:  self.path.input = self.input_file.with_suffix(".ini")
    
    # For linear simulations, read the *.ini file instead
    if self.object=="Simulation" and self.linear:
            
            # Read the unique inputs
            list_of_reference_inputs, ids_of_reference_inputs = read_uniqueInputs(self.path)

            # Find the reference input for each <mode> 
            for mode in self.simulation.modes:   
                mode.path.input = find_referenceInput(self.path, mode, list_of_reference_inputs, ids_of_reference_inputs)

            # Remember which inputs have already been read
            self.loaded_inputs = [] 
            
            # Load one of the files for the simulation
            self.path.input = mode.path.input
          
    # If we are directly asking it for the mode, load it for the entire simulation  
    if self.object=="Mode": self.simulation.path.input  
    return

#------------------------
def get_geometryPath(self): 
    """ Get the path to the *.geo.h5 file. """ 

    # Set the geometry path of a simulation
    if self.object=="Simulation":
        
        # Link the *.geo.h5 file for a nonlinear simulation
        if self.nonlinear:
            
            # Read the geometry from the *.geo.h5 file 
            rho = np.sqrt(read_svalueFromInputFile(self.input_file))*100 
            rho = int(rho) if (int(rho)==rho) else round(rho,2)  
            self.geometry = self.input_file.with_suffix(".rho"+str(rho)+".geometry") 
            
            # Make sure the file is written
            if not os.path.isfile(self.geometry) and os.path.isfile(self.output_stella):
                from stellapy.data.geometry.write_h5FileForGeometry import write_h5FileForGeometry 
                write_h5FileForGeometry(self.folder, self.geometry)
                    
        # Link the *unique.geometryX for a linear simulation
        if self.linear:

#             # Read the geometry from the *.geo.h5 file 
#             for mode in self.simulation.modes: 
#                 rho = np.sqrt(read_svalueFromInputFile(self.input_file))*100 
#                 rho = int(rho) if (int(rho)==rho) else round(rho,2)  
#                 mode.path.geometry = mode.input_file.with_suffix(".rho"+str(rho)+".geometry")  
#                 if not os.path.isfile(mode.path.geometry) and os.path.isfile(mode.path.output_stella):
#                     from stellapy.data.geometry.write_h5FileForGeometry import write_h5FileForGeometry 
#                     write_h5FileForGeometry(self.folder, mode.path.geometry, mode)
#             if not os.path.isfile(mode.path.geometry):
                
            # Read the unique geometries if we cant read from the specific geometry file
            list_of_reference_geometries, ids_of_reference_geometries = read_uniqueGeometries(self.path)

            # Find the reference geometry for each <mode>
            for mode in self.simulation.modes: 
                mode.path.geometry = find_referenceGeometry(self.path, mode, list_of_reference_geometries, ids_of_reference_geometries)
                    
            # Remember which geometries have already been read
            self.loaded_geometries = [] 
            
            # Make sure <simulation> also has a path.geometry variable
            self.simulation.path.geometry = "Don't use the geometry of the <simulation>, use the geometry of the <mode>."
          
    # If we are directly asking it for the mode, load it for the entire simulation  
    if self.object=="Mode": self.simulation.path.geometry
    return 

#------------------------
def get_pathWithSmallestTimeStep(self, identifier, oldidentifier=None, path="/not/found/"):  
    
    # Get the *.dt10.identifier file  
    files = get_filesInFolder(self.path.folder, start=self.input_file.stem+".", end=identifier) 
    if files!= None:  
        files = [f for f in files if "dt" in f.name]  
        files = [f for f in files if "OLD" not in f.parent.name] + [f for f in files if "OLD" in f.parent.name]
        if len(files)!=0:
            times = [float(f.name.split("dt")[-1].split("."+identifier)[0]) for f in files] 
            min_dt = np.min(np.array(times))
            files = [files[i] for i in range(len(files)) if times[i]==min_dt]
            path = files[0] 
        
    # If it doesn't exist find the old file
    if path and oldidentifier:
        if not os.path.isfile(path): 
            path = self.input_file.with_suffix(oldidentifier)
             
    return path 

#------------------------
def get_pathWithLargestTend(self, identifier, path="/not/found/"):  
    
    # Get the *.identifier.t500-1000 file  
    files = get_filesInFolder(self.path.folder, start=self.input_file.stem+".") 
    if files!=None:  
        if self.linear: 
            files = [f for f in files if "out.nc.tlast" in f.name] 
            if len(files)!=0: return files[0]
        files = [f for f in files if identifier in f.name] 
        if len(files)!=0:   
            files_no_skip = [f for f in files if "skip" not in f.name] 
            if len(files_no_skip)>0: files = files_no_skip
            if len(files_no_skip)==0:
                skips = [float(f.name.split("-")[-1]) for f in files] 
                min_skip = np.min(np.array(skips))  
                files = [files[i] for i in range(len(files)) if skips[i]==min_skip]
            tends = [float(f.name.split("-")[-1]) if ("skip" not in f.name) else float(f.name.split("-")[-2]) for f in files]   
            max_tend = np.max(np.array(tends))
            files = [files[i] for i in range(len(files)) if tends[i]==max_tend]
            path = files[0]   
    return path 
    
