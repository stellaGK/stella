 
import os, sys, h5py  
from stellapy.utils.decorators.exit_program import exit_program  
from stellapy.data.paths.load_pathObject import get_distribution4DPath
from stellapy.data.distribution.write_h5FileForDistribution4D import write_h5FileForDistribution4D

#===============================================================================
#                        ATTACH THE DISTRIBUTION DATA
#===============================================================================

def read_distribution4D(path): 
    ''' Read the distribution from *.distribution4D or *.out.nc ''' 
               
    # Make sure that the *.dt10.distribution4D file exists
    if not os.path.isfile(path.distribution4D) and os.path.isfile(path.output_stella): 
        write_h5FileForDistribution4D(path.folder)
        get_distribution4DPath(path)   
                              
    # Read from the *.dt10.distribution4D file 
    if os.path.isfile(path.distribution4D):
        return read_fromdistributionFile(path.distribution4D)    

    # Critical error if we didn't find any data
    exit_reason = "The distribution data could not be found."
    exit_program(exit_reason, read_distribution4D, sys._getframe().f_lineno)   
    return

#-------------------------------------
def read_fromdistributionFile(path, distribution={}):   
    with h5py.File(path, 'r') as f: 
        for key in ["vec_time", "g_vs_tsvpaz", "g_vs_tsmuvpa"]: 
            if key in f.keys():  
                distribution[key] = f[key][()]  
    return distribution 
            
#===============================================================================
#                  ATTACH THE distribution TO THE SIMULATION OBJECT                 #
#===============================================================================

def get_distribution4D(self): 
    
    # Read the distribution
    distribution = read_distribution4D(self.path)  
    
    # Save the distribution
    self.vec_time_4D = distribution["vec_time"]
    self.g_vs_tsvpaz = distribution["g_vs_tsvpaz"]
    self.g_vs_tsmuvpa = distribution["g_vs_tsmuvpa"] 
    return 
