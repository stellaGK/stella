 
#!/usr/bin/python3   
import sys, os
import numpy as np 
import pathlib, h5py
from scipy.io import netcdf as scnetcdf   
from stellapy.data.input.write_listOfMatchingInputFiles import read_inputFilesInDummyInputFile

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)    
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime   
from stellapy.data.input.read_inputFile import read_svalueFromInputFile,\
    read_numberOfModesFromInputFile, read_modeFromInputFile
from stellapy.utils.files.get_filesInFolder import get_filesInFolder   
from stellapy.utils.config import CONFIG

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
#    *.dt10.distribution3D     distribution3D      h5    (g2_vs_tsz, g2_vs_tsmu, g2_vs_tsvpa)
#    *.dt10.distribution4D     distribution4D      h5    (g2_vs_tsmuvpa, g2_vs_tsvpaz)
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
#    *.g2_vs_z                  g2_vs_z              txt    (g2_vs_sz)
#    *.g2_vs_mu                 g2_vs_mu             txt    (g2_vs_smu)
#    *.g2_vs_vpa                g2_vs_vpa            txt    (g2_vs_svpa) 
#    *.dt10.g2_vs_t             g2_vs_t              txt    (g2_vs_st)  
#    *.dt10.distribution3D     distribution        h5     (g2_vs_tsmu, g2_vs_tsvpa)
#    *.dt10.distribution4D     distribution        h5     (g2_vs_tsmuvpa)
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
        self.simulation = simulation  
        self.input_file = simulation.input_file   
        self.nonlinear = simulation.nonlinear
        self.linear = simulation.linear
        
        # Link itself for more logical access to the variables
        self.path = self   
        return 
    
    # Basic information of the simulation: folder and name
    @calculate_attributeWhenReadFirstTime 
    def name(self):                 get_nameSimulation(self);       return self.name 
    @calculate_attributeWhenReadFirstTime 
    def folder(self):               get_folderPath(self);           return self.folder 
    @calculate_attributeWhenReadFirstTime 
    def multiple_input_files(self): get_multipleInputFiles(self);   return self.multiple_input_files 
    @calculate_attributeWhenReadFirstTime 
    def input_files(self):          get_multipleInputFiles(self);   return self.input_files 
    @calculate_attributeWhenReadFirstTime 
    def dummy_input_file(self):     get_multipleInputFiles(self);   return self.dummy_input_file 
        
    # Set the paths when we ask for them
    @calculate_attributeWhenReadFirstTime 
    def vmec_filename(self):        get_vmecPath(self);             return self.vmec_filename 
    @calculate_attributeWhenReadFirstTime 
    def vmec(self):                 get_vmecPath(self);             return self.vmec 
    @calculate_attributeWhenReadFirstTime 
    def input(self):                get_inputPath(self);            return self.input 
    @calculate_attributeWhenReadFirstTime 
    def output(self):               get_outputPath(self);           return self.output 
    @calculate_attributeWhenReadFirstTime 
    def omega(self):                get_omegaPath(self);            return self.omega 
    @calculate_attributeWhenReadFirstTime 
    def geometry(self):             get_geometryPath(self);         return self.geometry 
    @calculate_attributeWhenReadFirstTime 
    def moments2D(self):            get_moments2DPath(self);        return self.moments2D 
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
    def QLfluxes(self):             get_QLfluxesPath(self);         return self.QLfluxes 
    @calculate_attributeWhenReadFirstTime 
    def fluxes(self):               get_fluxesPath(self);           return self.fluxes 
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
    def g2_vs_t(self):              get_g2vstPath(self);            return self.g2_vs_t 
    @calculate_attributeWhenReadFirstTime # Old name
    def g_vs_t(self):               get_gvstPath(self);             return self.g_vs_t 
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_t(self):             get_phivstPath(self);           return self.phi_vs_t 
    
    # Data at the last time step (linear)
    @calculate_attributeWhenReadFirstTime 
    def phi_vs_z(self):             get_phivszPath(self);           return self.phi_vs_z
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_z(self):               get_finalgPath(self);           return self.g2_vs_z 
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_mu(self):              get_finalgPath(self);           return self.g2_vs_mu 
    @calculate_attributeWhenReadFirstTime 
    def g2_vs_vpa(self):             get_finalgPath(self);           return self.g2_vs_vpa
    
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
#                            GET DUMMY INPUT FILES                            #
#=============================================================================== 
     
def get_multipleInputFiles(self):
    """ For 1 mode per simulation, we created a dummy input file. """
    
    # Attach whether we have a dummy input file 
    self.multiple_input_files = self.simulation.multiple_input_files
    self.input_files = None
    self.dummy_input_file = None
    
    # If we have a multiple input files, attach the input files to the <path> object 
    if self.multiple_input_files==True:
        
        # Get the input files corresponding to a similar simulation with different (kx,ky)
        self.input_files = self.simulation.input_files; self.paths = []
        self.input_files = [i for i in self.input_files if pathlib.Path(str(i).replace(".in", "_kx0.0.in")) not in self.input_files]

        # Create dummy path objects for each input file 
        for input_file in self.input_files:  
            self.paths.append(create_dummyPathObject(input_file, "/not/used"))
        
        # For each input file, remember the modes inside
        for path in self.paths: 
            path.dummy_input_file = None
            nakx, naky = read_numberOfModesFromInputFile(path.input_file)
            kx, ky = read_modeFromInputFile(path.input_file)
            path.nakxnaky = nakx*naky
            path.kx = kx 
            path.ky = ky
            if path.nakxnaky==1:
                path.dim_kx = 1
                path.dim_ky = 1
                path.vec_kx = [kx]
                path.vec_ky = [ky]
            if path.nakxnaky>1 or "_dummy.in" in str(path.input_file):
                with h5py.File(path.dimensions, 'r') as f:  
                    path.dim_kx = f["dim_kx"][()] 
                    path.dim_ky = f["dim_ky"][()] 
                    path.vec_kx = f["vec_kx"][()] 
                    path.vec_ky = f["vec_ky"][()]   
                    
        # For each input file, remember if it is part of a dummy input file
        for input_file in self.input_files: 
            if "_dummy.in" in str(input_file):
                dummy_input_files = read_inputFilesInDummyInputFile(input_file)   
                for path in self.paths: 
                    if path.input_file in dummy_input_files: path.dummy_input_file = input_file         
    return 

#===============================================================================
#                           CREATE DUMMY PATH OBJECT                           #
#=============================================================================== 

def create_dummyPathObject(path_input_file, vmec_filename):
    """ Needed to access all the paths of a simulation, without actually going
    through the <create_simulations> routine. Returns a working <path> object."""
    
    # Create a dummy input and simulation object 
    input_object = type('input', (object,), {'vmec_filename' : vmec_filename}) 
    dummy_object = type('Simulation', (object,), {'input' : input_object, 'input_file' : path_input_file, 'linear' : 'dummy', 'nonlinear' : 'dummy'})
    
    # Load the paths for the mode and simulation
    load_pathObject(dummy_object)  
    
    # Return the finished <path> object of the mode or simulation
    return dummy_object.path

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
#                        MAKE COMPATIBLE WITH OLD FILES                        #
#=============================================================================== 

def get_geometryPath(self): 
    self.path.geometry = self.input_file.with_suffix(".geo")       
    if not os.path.isfile(self.path.geometry):
        rho = np.sqrt(read_svalueFromInputFile(self.input_file))*100 
        rho = int(rho) if (int(rho)==rho) else round(rho,2)  
        self.geometry = self.input_file.with_suffix(".rho"+str(rho)+".geometry") 
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
def get_inputPath(self): 
    self.path.input = self.input_file.with_suffix(".ini")  
    return
 
def get_saturatedPath(self):
    self.path.saturated = get_pathWithLargestTend(self, "out.nc.t")
    return

def get_outputPath(self):
    self.path.output = get_pathWithSmallestTimeStep(self, "output", ".out.h5")
    return

def get_omegaPath(self):  
    self.path.omega = get_pathWithSmallestTimeStep(self, "omega_vs_t", ".omega_reduced") 
    return

def get_fluxesPath(self): 
    self.path.fluxes = get_pathWithSmallestTimeStep(self, "fluxes_vs_t", ".fluxes_reduced")
    return

def get_QLfluxesPath(self): 
    self.path.QLfluxes = self.input_file.with_suffix(".QLfluxes")  
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

def get_moments2DPath(self): 
    self.path.moments2D = get_pathWithSmallestTimeStep(self, "moments2D")
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

def get_g2vstPath(self): 
    self.path.g2_vs_t = get_pathWithSmallestTimeStep(self, "g2_vs_t")
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
    self.path.g2_vs_z = self.input_file.with_suffix(".g_vs_z")
    self.path.g2_vs_mu = self.input_file.with_suffix(".g_vs_mu")
    self.path.g2_vs_vpa = self.input_file.with_suffix(".g_vs_vpa")
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
    self.vmec_filename = self.simulation.input.vmec_filename
    vmec_filename = self.vmec_filename 
    if os.path.isfile(self.vmec_filename):
        self.path.vmec = pathlib.Path(self.vmec_filename)
    else:
        self.path.vmec = self.path.input.parent / vmec_filename if (vmec_filename != 'wout*.nc') else "pathNotFound"
    if (not os.path.isfile(self.path.vmec)) and (self.path.vmec!=None): 
        
        # Decide whether we're on the local pc or the supercomputer 
        vmec_folder = pathlib.Path(CONFIG["PATHS"]['VMEC']) 
         
        # See if the VMEC file is present in the main VMEC folder  
        if vmec_filename!='wout*.nc':   
            if os.path.isfile(vmec_folder / vmec_filename):
                self.path.vmec = vmec_folder / vmec_filename
            else:
                vmec_files = get_filesInFolder(vmec_folder, start=vmec_filename, end=".nc")
                self.path.vmec = vmec_files[0] if (vmec_files!=[]) else self.path.vmec  
        if vmec_filename=='wout*.nc':    
            self.path.vmec = 'MillerDoesntHaveAVmec'  
    return

#------------------------
def get_pathWithSmallestTimeStep(self, identifier, oldidentifier=None, path="/not/found/"):  
    
    # Get the *.dt10.identifier file  
    files = get_filesInFolder(self.path.folder, start=self.input_file.stem+".", end=identifier, search_in_subfolders=False) 
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
    files = get_filesInFolder(self.path.folder, start=self.input_file.stem+".", search_in_subfolders=False) 
    if files!=None:  
        if self.simulation.linear: 
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
    
