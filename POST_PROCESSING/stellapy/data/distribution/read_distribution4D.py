""" 

#===============================================================================
#                        Read the 4D distribution data                         #
#===============================================================================
 
Read the distribution function "g" from the *.dt10.distribution4D or *.out.nc files.

Returns
-------
    {vec_time, g2_vs_tsvpaz, g2_vs_tsmuvpa}

Hanne Thienpondt
27/01/2023

"""

#!/usr/bin/python3  
import h5py   
import sys, os
import pathlib

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)  
from stellapy.data.distribution.write_h5FileForDistribution4D import write_h5FileForDistribution4D
from stellapy.data.paths.load_pathObject import get_distribution4DPath
from stellapy.utils.decorators.exit_program import exit_program  
from stellapy.data.utils import Data

#===============================================================================
#                        Read the 4D distribution data                         #
#===============================================================================

def read_distribution4D(path): 
    ''' Read the distribution from *.dt10.distribution4D or *.out.nc ''' 
               
    # Make sure that the *.dt10.distribution4D file exists
    if not os.path.isfile(path.distribution4D) and os.path.isfile(path.output_stella): 
        write_h5FileForDistribution4D(path.folder)
        get_distribution4DPath(path)   
                              
    # Read from the *.dt10.distribution4D file 
    if os.path.isfile(path.distribution4D):
        return read_distributionFromH5File(path.distribution4D)    

    # Critical error if we didn't find any data
    exit_reason = "The 4D distribution data could not be found:\n" 
    exit_reason += "   "+path.distribution4D+"\n" 
    exit_program(exit_reason, read_distribution4D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_distributionFromH5File(path_distribution4D): 
    
    # Initiate
    distribution={}
    
    # Read the *.dt10.distribution4D file  
    with h5py.File(path_distribution4D, 'r') as f: 
        for key in ["vec_time", "g2_vs_tsvpaz", "g2_vs_tsmuvpa"]: 
            if key in f.keys():  
                distribution[key] = f[key][()]  
    return distribution 
            
#===============================================================================
#         Attach the distribution squared to the <distribution> object         #
#===============================================================================

def get_distribution4D(self): 
    
    # Read the distribution squared
    distribution = read_distribution4D(self.path)  
    
    # Save the distribution squared
    self.g2_vs_tsvpaz = Data(["g2","t","s","vpa","z"], distribution["g2_vs_tsvpaz"], distribution["vec_time"], self.vec.species, self.vec.vpa, self.vec.z) 
    self.g2_vs_tsmuvpa = Data(["g2","t","s","mu","vpa"], distribution["g2_vs_tsmuvpa"], self.g2_vs_tsvpaz.t, self.vec.species, self.vec.mu, self.vec.vpa) 
    return 
