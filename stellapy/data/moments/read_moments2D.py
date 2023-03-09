 
import os, sys, h5py  
from stellapy.data.utils import Data
from stellapy.utils.decorators.exit_program import exit_program 
from stellapy.data.paths.load_pathObject import get_moments2DPath
from stellapy.data.moments.write_h5FileForMoments2D import write_h5FileForMoments2D

#===============================================================================
#                        ATTACH THE DIMENSIONS DATA
#===============================================================================

def read_moments2D(path):
    ''' Read the moments from *.moments or *.out.nc ''' 
               
    # Make sure that the *.moments2D file exists
    if not os.path.isfile(path.moments2D) and os.path.isfile(path.output_stella): 
        write_h5FileForMoments2D(path.folder)
        get_moments2DPath(path)   
                           
    # Read from the *.moments2D file 
    if os.path.isfile(path.moments2D):
        return read_fromMomentsFile(path.moments2D)  

    # Critical error if we didn't find any data
    exit_reason = "The 2D moments data could not be found."
    exit_program(exit_reason, read_moments2D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromMomentsFile(path):
    moments={}
    with h5py.File(path, 'r') as f: 
        for key in ["vec_time", "upar_vs_ts", "dens_vs_ts", "temp_vs_ts", "upar2_vs_ts", "dens2_vs_ts", "temp2_vs_ts"]: 
                moments[key] = f[key][()]   
    return moments 
            
#===============================================================================
#                  ATTACH THE MOMENTS TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_moments2D(self): 
    
    # Read the moments
    moments = read_moments2D(self.path) 
    
    # Save the moments
    self.upar_vs_ts = Data(["upar","t","s"], moments["upar_vs_ts"], moments["vec_time"], self.vec.species) 
    self.dens_vs_ts = Data(["dens","t","s"], moments["dens_vs_ts"], moments["vec_time"], self.vec.species) 
    self.temp_vs_ts = Data(["temp","t","s"], moments["temp_vs_ts"], moments["vec_time"], self.vec.species) 
    self.upar2_vs_ts = Data(["upar2","t","s"], moments["upar2_vs_ts"], moments["vec_time"], self.vec.species) 
    self.dens2_vs_ts = Data(["dens2","t","s"], moments["dens2_vs_ts"], moments["vec_time"], self.vec.species) 
    self.temp2_vs_ts = Data(["temp2","t","s"], moments["temp2_vs_ts"], moments["vec_time"], self.vec.species) 

    return 

