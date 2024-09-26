  
import numpy as np 
from datetime import datetime   

#===============================================================================
#                      READ THE TIME OF THE OVERSHOOT/PEAK                     #
#=============================================================================== 

def read_peakTime(self): 
    
    # Check whether the fluxes file has been changed since calculating the data    
    if self.fluxes.date > datetime.strptime(self.section["peakdate"], '%Y-%m-%d %H:%M:%S.%f') or self.section["peak"]=="/":
        print(" -> REWRITE PEAK TIME BECAUSE THE FLUXES FILE HAS BEEN TOUCHED") 
        self.section["peakdate"] = str(datetime.now())
        self.section["peak"] = str(calculate_peakTime(self.fluxes.qflux_vs_ts.t, self.fluxes.qflux_vs_ts.qflux[:,0]))
        self.file.write(open(self.path.folder/"timeFrames.ini", 'w'))   
        
    # Return the peak time
    return float(self.section["peak"])

#===============================================================================
#  CALCULATE THE TIME OF THE OVERSHOOT OF THE QFLUX OF A NONLINEAR SIMULATION  #
#=============================================================================== 

def calculate_peakTime(vec_time, vec_flux, box=4):
    ''' Calculate the time of the peak of the ramp up. '''  
    time = vec_time[vec_flux>0.001] 
    fluxes = vec_flux[vec_flux>0.001]
    fluxes = [np.mean(fluxes[i:i+box]) for i in range(len(fluxes)-box-1)]
    fluxes = [ fluxes[i+1]-fluxes[i] for i in range(len(fluxes)-1)]
    index_peak = [ i for i in range(len(fluxes)) if fluxes[i]<0] 
    fluxes = vec_flux[vec_flux>0.001][index_peak[0]:index_peak[0]+box]
    index_peak = index_peak + np.argmax(fluxes)
    time_peak = time[index_peak[0]]
    return time_peak 
        
#===============================================================================
#           ATTACH THE TIME OF THE OVERSHOOT/PEAK TO THE SIMULATION            #
#=============================================================================== 

def get_peakTime(self):
    self.peakTime = read_peakTime(self)
    return
