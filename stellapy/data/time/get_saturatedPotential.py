  
import numpy as np
from datetime import datetime  

#===============================================================================
#                    ATTACH THE SATURATED FLUXES TO <FLUXES>                   #
#=============================================================================== 
 
def get_saturatedPotential(self):
     
    # Check whether the fluxes file has been changed since calculating the data  
    rewrite_potential = False  
    if "phi2" not in self.section: rewrite_potential = True
    if rewrite_potential==False: 
        if self.section["phi2"]=="/": rewrite_potential = True
    if self.fluxes.date > datetime.strptime(self.section["fluxdate"], '%Y-%m-%d %H:%M:%S.%f') or rewrite_potential:
        print(" -> REWRITE SATURATED POTENTIAL BECAUSE THE SATURATED FILE HAS BEEN TOUCHED") 
        self.section["fluxdate"] = str(datetime.now()) 
        self.saturatedPotential, self.satPotStdErrors, self.satPotMinimum, self.satPotMaximum = calculate_saturatedPotential(self) 
        write_saturatedPotential(self.section, self.saturatedPotential)
        write_saturatedPotential(self.section, self.satPotMinimum, "min")
        write_saturatedPotential(self.section, self.satPotMaximum, "max")
        write_saturatedPotential(self.section, self.satPotStdErrors, "std")
        self.file.write(open(self.path.folder/"timeFrames.ini", 'w'))  
    
    # Otherwise read the saturatedPotential data
    else:
        self.saturatedPotential = read_saturatedPotential(self.section)
        self.satPotMinimum  = read_saturatedPotential(self.section, "min")
        self.satPotMaximum  = read_saturatedPotential(self.section, "max")
        self.satPotStdErrors = read_saturatedPotential(self.section, "std") 
    return   
            
#----------------------------------
def write_saturatedPotential(section, satpot, extra=""):
    for pot in ['phi2', 'phi2zonal', 'phi2nozonal']: 
        section[pot+extra] = str(satpot[pot])
    return
    
#----------------------------------
def read_saturatedPotential(section, extra=""):
    satpot = {"phi2" : [], "phi2zonal" : [], "phi2nozonal" : []}
    for pot in ['phi2', 'phi2zonal', 'phi2nozonal']:
        satpot[pot] = float(section[pot+extra]) 
    return satpot
    
#===============================================================================
#              AVERAGE THE POTENTIAL OVER THE SELECTED TIME FRAME              #
#=============================================================================== 
 
def calculate_saturatedPotential(self):   
    """ The object is a <time> object. """
          
    # Get the time vector and the time filter
    vec_time = self.potential.phi_vs_t.t
    tfilter = (vec_time >= self.tstart) & (vec_time <= self.tend)
    
    # Initiate the saturated fluxes dictionary 
    saturatedPotential = {}
    satPotStdErrors = {}
    satPotMinimum = {}
    satPotMaximum = {}
  
    # Iterate over the (phi2, phi2zonal, phi2nozonal):
    for pot in ['phi2', 'phi2zonal', 'phi2nozonal']:
          
        # Get the potential
        if pot=="phi2":         vec_phi2 = self.potential.phi2_vs_t.phi2[tfilter] 
        if pot=="phi2zonal":    vec_phi2 = self.potential.phi2_vs_t_zonal.phi2_zonal[tfilter] 
        if pot=="phi2nozonal":  vec_phi2 = self.potential.phi2_vs_t_nozonal.phi2_nozonal[tfilter]  

        # Calculate the saturated potetial 
        saturatedPotential[pot] = np.nanmean(vec_phi2)
        satPotStdErrors[pot] = np.std(vec_phi2)
        satPotMinimum[pot] = - np.min(vec_phi2) + saturatedPotential[pot]
        satPotMaximum[pot] = np.max(vec_phi2) - saturatedPotential[pot]
 
    return saturatedPotential, satPotStdErrors, satPotMinimum, satPotMaximum

