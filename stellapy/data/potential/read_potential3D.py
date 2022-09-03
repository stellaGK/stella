 
import numpy as np
import os, sys, h5py  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program  
from stellapy.data.paths.load_pathObject import get_potential3DPath
from stellapy.data.potential.write_h5FileForPotential3D import write_h5FileForPotential3D

#===============================================================================
#                        ATTACH THE POTENTIAL DATA
#===============================================================================

def read_potential3D(path): 
    ''' Read the potential from *.potential3D or *.out.nc ''' 
               
    # Make sure that the *.dt10.write_h5FileForPotential3D file exists
    if not os.path.isfile(path.potential3D) and os.path.isfile(path.output_stella): 
        write_h5FileForPotential3D(path.folder)
        get_potential3DPath(path)  
                              
    # Read from the *.dt10.potential3D file 
    if os.path.isfile(path.potential3D):
        return read_fromPotentialFile(path.potential3D)    
    
    # This file doesn't exist for old simulations
    else:  
        return read_fromH5File(path)
    
    # Critical error if we didn't find any data
    exit_reason = "The potential data could not be found."
    exit_program(exit_reason, read_potential3D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromPotentialFile(path): 
    potential={} 
    with h5py.File(path, 'r') as f: 
        for key in ["vec_time", "phi_vs_tz", "phi2_vs_tz", "phi2_vs_tkx", "phi2_vs_tky", "phi_vs_tkx", "phi_vs_tky"]:
            for extra in ["", "_zonal", "_nozonal"]: 
                if key+extra in f.keys():   
                    potential[key+extra] = f[key+extra][()]  
    return potential 

#===============================================================================
#                            READ FROM OLD H5 FILE                             #
#===============================================================================

def read_fromH5File(path):  
    potential_data={} 
    with h5py.File(path.output, 'r') as f:    
        potential_data["vec_time"] = f["vec_time"][()]
        potential_data["phi_vs_tz"] = f["phi_vs_tzri"][()][:,:,0]+1j*f["phi_vs_tzri"][()][:,:,1]
        potential_data["phi_vs_tz_zonal"] = f["phi_vs_tzri_zonal"][()][:,:,0]+1j*f["phi_vs_tzri_zonal"][()][:,:,1]
        potential_data["phi_vs_tz_nozonal"] = f["phi_vs_tzri_nozonal"][()][:,:,0]+1j*f["phi_vs_tzri_nozonal"][()][:,:,1]  
        potential_data["phi2_vs_tz"] = np.ones((np.shape(potential_data["phi_vs_tz"])))*np.nan
        potential_data["phi2_vs_tz_zonal"] = np.ones((np.shape(potential_data["phi_vs_tz"])))*np.nan
        potential_data["phi2_vs_tz_nozonal"] = np.ones((np.shape(potential_data["phi_vs_tz"])))*np.nan
        phi2_vs_tkxky = f["phi2_vs_kxky"][()] if ("phi2_vs_kxky" in f.keys()) else f["phi2_vs_tkxky"][()]  
        potential_data["phi2_vs_tkx"] = phi2_vs_tkxky[:,:,0] + 2*np.sum(phi2_vs_tkxky[:,:,1:],axis=2)  
        potential_data["phi2_vs_tky"] = np.sum(phi2_vs_tkxky[:,:,:],axis=1)  
    return potential_data

#===============================================================================
#                  ATTACH THE POTENTIAL TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_potential3D(self): 

    # Read the potential
    potential = read_potential3D(self.path)   
    
    # Save the potential
    self.phi_vs_tz = Data(["phi","t","z"], potential["phi_vs_tz"], potential["vec_time"], self.vec.z) 
    self.phi2_vs_tz = Data(["phi2","t","z"], potential["phi2_vs_tz"], self.phi_vs_tz.t, self.phi_vs_tz.z)   
    
    # Save the potential along (kx) or (ky)
    if "phi2_vs_tkx" in potential: 
        self.phi2_vs_tkx = Data(["phi2","t","kx"], potential["phi2_vs_tkx"], self.phi_vs_tz.t, self.vec.kx)   
        self.phi2_vs_tky = Data(["phi2","t","ky"], potential["phi2_vs_tky"], self.phi_vs_tz.t, self.vec.ky)  
        
        # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
        # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
        sorted_indexes = np.argsort(self.vec.kx_stella)   
        self.phi2_vs_tkx.phi2 = self.phi2_vs_tkx.phi2[:,sorted_indexes] 
    else:
        self.phi2_vs_tkx = np.nan
        self.phi2_vs_tky = np.nan
    
    # Save the potential along (kx) or (ky)
    if "phi_vs_tkx" in potential: 
        self.phi_vs_tkx = Data(["phi","t","kx"], potential["phi_vs_tkx"], self.phi_vs_tz.t, self.vec.kx)   
        self.phi_vs_tky = Data(["phi","t","ky"], potential["phi_vs_tky"], self.phi_vs_tz.t, self.vec.ky)  
        
        # Stella has vec_kx_stella = [0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx]
        # So put pflx_kxky(t, kx, ky, s) in ascending order with respect to kx.
        sorted_indexes = np.argsort(self.vec.kx_stella)   
        self.phi_vs_tkx.phi2 = self.phi_vs_tkx.phi[:,sorted_indexes] 
    else:
        self.phi_vs_tkx = np.nan
        self.phi_vs_tky = np.nan
    
    # Save the potential split in zonal/nonzonal contributions
    if "phi_vs_tz_zonal" in potential: 
        self.phi_vs_tz_zonal    = Data(["phi_zonal","t","z"],   potential["phi_vs_tz_zonal"],    self.phi_vs_tz.t, self.phi_vs_tz.z)
        self.phi_vs_tz_nozonal  = Data(["phi_nozonal","t","z"], potential["phi_vs_tz_nozonal"],  self.phi_vs_tz.t, self.phi_vs_tz.z)
        self.phi2_vs_tz_zonal   = Data(["phi2_zonal","t","z"],  potential["phi2_vs_tz_zonal"],   self.phi_vs_tz.t, self.phi_vs_tz.z)
        self.phi2_vs_tz_nozonal = Data(["phi2_nozonal","t","z"],potential["phi2_vs_tz_nozonal"], self.phi_vs_tz.t, self.phi_vs_tz.z)
    return 
