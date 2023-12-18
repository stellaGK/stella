""" 

#===============================================================================
#                        Read the 3D distribution data                         #
#===============================================================================
 
Read the distribution function "g" from the *.dt10.distribution3D or *.out.nc files.

Returns
-------
    {vec_time, g2_vs_tsz, g2_vs_tsmu, g2_vs_tsvpa}

Hanne Thienpondt
27/01/2023

"""

#!/usr/bin/python3  
import h5py
import os, sys 
import pathlib

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.distribution.write_h5FileForDistribution3D import write_h5FileForDistribution3D 
from stellapy.data.paths.load_pathObject import get_distribution3DPath 
from stellapy.utils.decorators.exit_program import exit_program  
from stellapy.data.utils import Data

#===============================================================================
#                        Read the 3D distribution data                         #
#===============================================================================

def read_distribution3D(path): 
    ''' Read the distribution from *.dt10.distribution3D or *.out.nc ''' 
               
    # Make sure that the *.dt10.distribution3D file exists
    if not os.path.isfile(path.distribution3D) and os.path.isfile(path.output_stella): 
        write_h5FileForDistribution3D(path.folder)
        get_distribution3DPath(path)  
                              
    # Read from the *.dt10.distribution3D file 
    if os.path.isfile(path.distribution3D):
        return read_distributionFromH5File(path.distribution3D)    
    
    # Critical error if we didn't find any data
    exit_reason = "The 3D distribution data could not be found:\n" 
    exit_reason += "   "+path.distribution3D+"\n" 
    exit_program(exit_reason, read_distribution3D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_distributionFromH5File(path_distribution3D):  
    
    # Initiate
    distribution={}
    
    # Read the *.dt10.distribution3D file
    with h5py.File(path_distribution3D, 'r') as f: 
        for key in ["vec_time", "g2_vs_tsz", "g2_vs_tsmu", "g2_vs_tsvpa"]: 
            if key in f.keys():   
                distribution[key] = f[key][()]  
        for old_key in ["g_vs_tsz", "g_vs_tsmu", "g_vs_tsvpa"]: 
            if old_key in f.keys():   
                distribution[old_key.replace("g_","g2_")] = f[old_key][()]  
    return distribution 

#===============================================================================
#         Attach the distribution squared to the <distribution> object         #
#===============================================================================

def get_distribution3D(self): 
    
    # Read the distribution squared
    distribution = read_distribution3D(self.path)  
    
    # Save the distribution squared
    self.g2_vs_tsz = Data(["g2","t","z"], distribution["g2_vs_tsz"], distribution["vec_time"], self.vec.species, self.vec.z) 
    self.g2_vs_tsmu = Data(["g2","t","mu"], distribution["g2_vs_tsmu"], self.g2_vs_tsz.t, self.vec.species, self.vec.mu) 
    self.g2_vs_tsvpa = Data(["g2","t","vpa"], distribution["g2_vs_tsvpa"], self.g2_vs_tsz.t, self.vec.species, self.vec.vpa)  
    return 
