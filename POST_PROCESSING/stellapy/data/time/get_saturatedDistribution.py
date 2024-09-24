
import numpy as np
from datetime import datetime  

#===============================================================================
#                ATTACH THE SATURATED DISTRIBUTION TO <FLUXES>                 #
#=============================================================================== 
 
def get_saturatedDistribution(self):
     
    # Check whether the fluxes file has been changed since calculating the data  
    rewrite_distribution = False  
    if "gi" not in self.section: rewrite_distribution = True
    if rewrite_distribution==False: 
        if self.section["gi"]=="/": rewrite_distribution = True
    if self.fluxes.date > datetime.strptime(self.section["fluxdate"], '%Y-%m-%d %H:%M:%S.%f') or rewrite_distribution:
        print(" -> REWRITE SATURATED DISTRIBUTION BECAUSE THE SATURATED FILE HAS BEEN TOUCHED") 
        self.section["fluxdate"] = str(datetime.now()) 
        self.saturatedDistribution, self.satDistStdErrors, self.satDistMinimum, self.satDistMaximum = calculate_saturatedDistribution(self) 
        write_saturatedDistribution(self.section, self.saturatedDistribution)
        write_saturatedDistribution(self.section, self.satDistMinimum, "min")
        write_saturatedDistribution(self.section, self.satDistMaximum, "max")
        write_saturatedDistribution(self.section, self.satDistStdErrors, "std")
        self.file.write(open(self.path.folder/"timeFrames.ini", 'w'))  
    
    # Otherwise read the saturatedDistribution data
    else:
        self.saturatedDistribution = read_saturatedDistribution(self.section)
        self.satDistMinimum  = read_saturatedDistribution(self.section, "min")
        self.satDistMaximum  = read_saturatedDistribution(self.section, "max")
        self.satDistStdErrors = read_saturatedDistribution(self.section, "std") 
    return   
            
#----------------------------------
def write_saturatedDistribution(section, satpot, extra=""):
    for pot in ['ge', 'gi']: 
        section[pot+extra] = str(satpot[pot])
    return
    
#----------------------------------
def read_saturatedDistribution(section, extra=""):
    satpot = {"ge" : [], "gi" : []}
    for pot in ['ge', 'gi']:
        satpot[pot] = float(section[pot+extra]) 
    return satpot
    
#===============================================================================
#              AVERAGE THE POTENTIAL OVER THE SELECTED TIME FRAME              #
#=============================================================================== 
 
def calculate_saturatedDistribution(self):   
    """ The object is a <time> object. """
          
    # Get the time vector and the time filter
    vec_time = self.distribution.g_vs_ts.t
    tfilter = (vec_time >= self.tstart) & (vec_time <= self.tend)
    
    # Initiate the saturated fluxes dictionary 
    saturatedDistribution = {}
    satDistStdErrors = {}
    satDistMinimum = {}
    satDistMaximum = {}
  
    # Iterate over the (ge, gi):
    for pot in ['ge', 'gi']:
          
        # Get the distribution
        if pot=="gi": vec_g = self.distribution.g_vs_ts.g[tfilter,0] 
        if pot=="ge": vec_g = self.distribution.g_vs_ts.g[tfilter,1] 

        # Calculate the saturated distribution 
        saturatedDistribution[pot] = np.nanmean(vec_g)
        satDistStdErrors[pot] = np.std(vec_g)
        satDistMinimum[pot] = - np.min(vec_g) + saturatedDistribution[pot]
        satDistMaximum[pot] = np.max(vec_g) - saturatedDistribution[pot]
 
    return saturatedDistribution, satDistStdErrors, satDistMinimum, satDistMaximum
