 
import numpy as np
from datetime import datetime  

#===============================================================================
#                    ATTACH THE SATURATED FLUXES TO <FLUXES>                   #
#=============================================================================== 
 
def get_saturatedFluxes(self):
     
    # Check whether the fluxes file has been changed since calculating the data    
    if self.fluxes.date > datetime.strptime(self.section["fluxdate"], '%Y-%m-%d %H:%M:%S.%f') or self.section["qflux"]=="/":
        print(" -> REWRITE SATURATED FLUXES BECAUSE THE FLUXES FILE HAS BEEN TOUCHED") 
        self.section["fluxdate"] = str(datetime.now()) 
        self.saturatedFluxes, self.satFluxStdErrors, self.satFluxMinimum, self.satFluxMaximum = calculate_saturatedFluxes(self) 
        write_saturatedFluxes(self.section, self.saturatedFluxes)
        write_saturatedFluxes(self.section, self.satFluxMinimum, "min")
        write_saturatedFluxes(self.section, self.satFluxMaximum, "max")
        write_saturatedFluxes(self.section, self.satFluxStdErrors, "std")
        self.file.write(open(self.path.folder/"timeFrames.ini", 'w'))  
    
    # Otherwise read the saturatedFluxes data
    else:
        self.saturatedFluxes = read_saturatedFluxes(self.section)
        self.satFluxMinimum  = read_saturatedFluxes(self.section, "min")
        self.satFluxMaximum  = read_saturatedFluxes(self.section, "max")
        self.satFluxStdErrors = read_saturatedFluxes(self.section, "std") 
    return   
            
#----------------------------------
def write_saturatedFluxes(section, satflux, extra=""):
    for flux in ['qflux', 'pflux', 'vflux']:
        satflux[flux] = [str(i) for i in satflux[flux]] 
        section[flux+extra] = "[" + ", ".join(satflux[flux]) + "]" 
        satflux[flux] = [float(i) for i in satflux[flux]] 
    return
    
#----------------------------------
def read_saturatedFluxes(section, extra=""):
    satflux = {"qflux" : [], "pflux" : [], "vflux" : []}
    for flux in ['qflux', 'pflux', 'vflux']:
        satflux[flux] = section[flux+extra].split("[")[-1].split("]")[0].split(", ")
        satflux[flux] = [float(f) for f in satflux[flux]]
    return satflux
    
#===============================================================================
#                AVERAGE THE FLUX OVER THE SELECTED TIME FRAME                 #
#=============================================================================== 
 
def calculate_saturatedFluxes(self):   
    """ The object is a <time> object. """
          
    # Get the time vector and the time filter
    vec_time = self.fluxes.qflux_vs_ts.t
    tfilter = (vec_time >= self.tstart) & (vec_time <= self.tend)
    
    # Initiate the saturated fluxes dictionary
    satFluxMinimum = {"qflux" : [], "pflux" : [], "vflux" : []}
    satFluxMaximum = {"qflux" : [], "pflux" : [], "vflux" : []}
    saturatedFluxes = {"qflux" : [], "pflux" : [], "vflux" : []}
    satFluxStdErrors = {"qflux" : [], "pflux" : [], "vflux" : []}
  
    # Iterate over the species and fluxes
    for s in range(self.dim.species): 
        for flux in ['qflux', 'pflux', 'vflux']:
              
            # Get the fluxes
            if flux=="qflux": vec_flux = self.fluxes.qflux_vs_ts.qflux[tfilter,s] 
            if flux=="pflux": vec_flux = self.fluxes.pflux_vs_ts.pflux[tfilter,s] 
            if flux=="vflux": vec_flux = self.fluxes.vflux_vs_ts.vflux[tfilter,s] 
                  
            # For the fluxes we can look at Q or Q/A   
            if flux=="qflux" and self.fluxes.includeFluxNorm:    
                from stellapy.data.fluxes.calculate_fluxnorm import calculate_fluxnorm  
                vec_flux = vec_flux*calculate_fluxnorm(self.input.vmec_filename, self.input.poloidal_turns) 

            # Calculate the saturated flux 
            saturatedFluxes[flux].append(np.nanmean(vec_flux))
            satFluxStdErrors[flux].append(np.std(vec_flux))
            satFluxMinimum[flux].append(- np.min(vec_flux) + saturatedFluxes[flux][s])
            satFluxMaximum[flux].append(np.max(vec_flux) - saturatedFluxes[flux][s]) 
 
    return saturatedFluxes, satFluxStdErrors, satFluxMinimum, satFluxMaximum
