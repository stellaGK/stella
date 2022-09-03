
import numpy as np
import os, sys, h5py  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program 
from stellapy.data.paths.load_pathObject import get_moments3DPath
from stellapy.data.moments.write_h5FileForMoments3D import write_h5FileForMoments3D

#===============================================================================
#                        ATTACH THE DIMENSIONS DATA
#===============================================================================

def read_phaseshifts(path):
    ''' Read the moments from *.phaseshifts or *.out.nc ''' 
               
    # Make sure that the *.phaseshifts file exists
    if not os.path.isfile(path.phaseshifts) and os.path.isfile(path.output_stella): 
        write_h5FileForMoments3D(path.folder)
        get_moments3DPath(path)   
                           
    # Read from the *.phaseshifts file 
    if os.path.isfile(path.phaseshifts):
        return read_fromMomentsFile(path.phaseshifts)  

    # Critical error if we didn't find any data
    exit_reason = "The phase shifts data could not be found."
    exit_program(exit_reason, read_phaseshifts, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromMomentsFile(path):
    moments={} 
    try:
        with h5py.File(path, 'r') as f:  
            if "Finished" not in f.keys():   
                print("Remove this phase shift file: ")
                print(path)
                sys.exit()
            for key in ["vec_time", "tstarts", "tends"]: 
                if key in f.keys():   
                    moments[key] = f[key][()]  
            for key in ["phi_and_n", "phi_and_T"]: 
                for extra in ["0_vs_tangleky", "0_vs_tanglekx", "1_vs_tangleky", "1_vs_tanglekx", "0_weights", "1_weights"]: 
                    if key+extra in f.keys():   
                        moments[key+extra] = f[key+extra][()] 
            for key in ["n0_and_T0", "n1_and_T1", "n0_and_n1", "T0_and_T1"]: 
                for extra in ["_vs_tangleky", "_vs_tanglekx", "_weights"]: 
                    if key+extra in f.keys():   
                        moments[key+extra] = f[key+extra][()] 
    except: 
        print("Something went wrong when reading:")
        print("    ", str(path))
        sys.exit()
    return moments 
            
#===============================================================================
#                  ATTACH THE MOMENTS TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_phaseshifts(self): 
    
    # Read the moments
    moments = read_phaseshifts(self.path) 
    angles = np.linspace(-180, 180, 361)
    
    # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
    # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
    sorted_indexes = np.argsort(self.vec.kx_stella)  
    moments["phi_and_n0_vs_tanglekx"] = moments["phi_and_n0_vs_tanglekx"][:,:,sorted_indexes]
    moments["phi_and_T0_vs_tanglekx"] = moments["phi_and_T0_vs_tanglekx"][:,:,sorted_indexes] 
    moments["n0_and_T0_vs_tanglekx"] = moments["n0_and_T0_vs_tanglekx"][:,:,sorted_indexes] 
    if "phi_and_n1_vs_tanglekx" in moments.keys(): 
        moments["phi_and_n1_vs_tanglekx"] = moments["phi_and_n1_vs_tanglekx"][:,:,sorted_indexes] 
        moments["phi_and_T1_vs_tanglekx"] = moments["phi_and_T1_vs_tanglekx"][:,:,sorted_indexes] 
        moments["n1_and_T1_vs_tanglekx"] = moments["n1_and_T1_vs_tanglekx"][:,:,sorted_indexes] 
        moments["n0_and_n1_vs_tanglekx"] = moments["n0_and_n1_vs_tanglekx"][:,:,sorted_indexes] 
        moments["T0_and_T1_vs_tanglekx"] = moments["T0_and_T1_vs_tanglekx"][:,:,sorted_indexes] 
    
    # Save the moments
    self.tends = moments["tends"]
    self.tstarts = moments["tstarts"]
    self.n0_and_T0_weights = moments["n0_and_T0_weights"]
    self.phi_and_n0_weights = moments["phi_and_n0_weights"]
    self.phi_and_T0_weights = moments["phi_and_T0_weights"]
    self.phase_shifts_phi_and_n0_vs_tanglekx = Data(["phaseshift","t","angle","kx"], moments["phi_and_n0_vs_tanglekx"], moments["vec_time"], angles, self.vec.kx) 
    self.phase_shifts_phi_and_T0_vs_tanglekx = Data(["phaseshift","t","angle","kx"], moments["phi_and_T0_vs_tanglekx"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.kx)
    self.phase_shifts_n0_and_T0_vs_tanglekx = Data(["phaseshift","t","angle","kx"], moments["n0_and_T0_vs_tanglekx"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.kx)
    self.phase_shifts_phi_and_n0_vs_tangleky = Data(["phaseshift","t","angle","ky"], moments["phi_and_n0_vs_tangleky"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.ky) 
    self.phase_shifts_phi_and_T0_vs_tangleky = Data(["phaseshift","t","angle","ky"], moments["phi_and_T0_vs_tangleky"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.ky)
    self.phase_shifts_n0_and_T0_vs_tangleky = Data(["phaseshift","t","angle","ky"], moments["n0_and_T0_vs_tangleky"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.ky)
    if "phi_and_n1_vs_tanglekx" in moments.keys(): 
        self.n0_and_n1_weights = moments["n0_and_n1_weights"]
        self.T0_and_T1_weights = moments["T0_and_T1_weights"]
        self.n1_and_T1_weights = moments["n1_and_T1_weights"]
        self.phi_and_n1_weights = moments["phi_and_n1_weights"]
        self.phi_and_T1_weights = moments["phi_and_T1_weights"]
        self.phase_shifts_phi_and_n1_vs_tanglekx = Data(["phaseshift","t","angle","kx"], moments["phi_and_n1_vs_tanglekx"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.kx) 
        self.phase_shifts_phi_and_T1_vs_tanglekx = Data(["phaseshift","t","angle","kx"], moments["phi_and_T1_vs_tanglekx"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.kx)
        self.phase_shifts_n1_and_T1_vs_tanglekx = Data(["phaseshift","t","angle","kx"], moments["n1_and_T1_vs_tanglekx"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.kx)
        self.phase_shifts_n0_and_n1_vs_tanglekx = Data(["phaseshift","t","angle","kx"], moments["n0_and_n1_vs_tanglekx"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.kx)
        self.phase_shifts_T0_and_T1_vs_tanglekx = Data(["phaseshift","t","angle","kx"], moments["T0_and_T1_vs_tanglekx"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.kx)
        self.phase_shifts_phi_and_n1_vs_tangleky = Data(["phaseshift","t","angle","ky"], moments["phi_and_n1_vs_tangleky"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.ky) 
        self.phase_shifts_phi_and_T1_vs_tangleky = Data(["phaseshift","t","angle","ky"], moments["phi_and_T1_vs_tangleky"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.ky)
        self.phase_shifts_n1_and_T1_vs_tangleky = Data(["phaseshift","t","angle","ky"], moments["n1_and_T1_vs_tangleky"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.ky)
        self.phase_shifts_n0_and_n1_vs_tangleky = Data(["phaseshift","t","angle","ky"], moments["n0_and_n1_vs_tangleky"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.ky)
        self.phase_shifts_T0_and_T1_vs_tangleky = Data(["phaseshift","t","angle","ky"], moments["T0_and_T1_vs_tangleky"], self.phase_shifts_phi_and_n0_vs_tanglekx.t, angles, self.vec.ky)
    else:
        self.phase_shifts_phi_and_n1_vs_tanglekx = None
        self.phase_shifts_phi_and_T1_vs_tanglekx = None
        self.phase_shifts_n1_and_T1_vs_tanglekx = None
        self.phase_shifts_n0_and_n1_vs_tanglekx = None
        self.phase_shifts_T0_and_T1_vs_tanglekx = None
        self.phase_shifts_phi_and_n1_vs_tangleky = None
        self.phase_shifts_phi_and_T1_vs_tangleky = None
        self.phase_shifts_n1_and_T1_vs_tangleky = None
        self.phase_shifts_n0_and_n1_vs_tangleky = None
        self.phase_shifts_T0_and_T1_vs_tangleky = None
        self.n0_and_n1_weights = None
        self.T0_and_T1_weights = None
        self.n1_and_T1_weights = None
        self.phi_and_n1_weights = None
        self.phi_and_T1_weights = None
    return 

