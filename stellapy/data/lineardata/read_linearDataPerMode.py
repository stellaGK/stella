
import os
import numpy as np
import configparser 
from stellapy.GUI.widgets.Progress import show_progressWhenReadingFiles

#===============================================================================
#                              GET THE LINEAR DATA                             #
#=============================================================================== 

def get_linearDataPerMode(self):
    ''' Read the linear data from the simulation file or calculate it.  
    Grab omega(t_last); gamma(t_last) and fluxes(t_last). Here t_last is the last
    time value of the simulation, assuming the simulation has converged at this point. '''

    # Check whether the file exists 
    if not os.path.isfile(self.input_file.with_suffix(".lineardata.ini")):
        calculate_linearData(self, self.input_file)  
    
    # If the linear data file is older than omega(t), then rewrite it
    try: time_omegaFile = self.path.omega.stat().st_mtime
    except: time_omegaFile = self.path.omega_stella.stat().st_mtime
    time_linearDataFile = self.input_file.with_suffix(".lineardata.ini").stat().st_mtime
    if time_omegaFile>time_linearDataFile:
        show_progressWhenReadingFiles(self, "Calculating linear data")  
        calculate_linearData(self, self.input_file)  
        return
                        
    # If the linear data isn't written to "*.lineardata" then calculate and write it.
    if not os.path.isfile(self.input_file.with_suffix(".lineardata.ini")):  
        show_progressWhenReadingFiles(self, "Calculating linear data")  
        calculate_linearData(self, self.input_file)  
        return
        
    # If the linear data is written to "*.lineardata" then read it
    elif os.path.isfile(self.input_file.with_suffix(".lineardata.ini")):  
        show_progressWhenReadingFiles(self, "Reading linear data")  
        read_linearDataPerMode(self, self.input_file)  
        return
        
#===============================================================================
#                             READ THE LINEAR DATA                             #
#===============================================================================  

def read_linearDataPerMode(self, input_file):
    
    # Read the linear data from "*.lineardata"  
    file = configparser.ConfigParser() 
    path = input_file.with_suffix(".lineardata.ini")
    file.read(path) 
     
    # Read the linear data of the (kx,ky) mode
    mode_identifier = '('+str(self.kx)+', '+str(self.ky)+')'
    if mode_identifier not in file: calculate_linearData(self, self.input_file); return 
    data = file[mode_identifier] 
    
    # Check whether we have the data
    if "stable" not in data:
        calculate_linearData(self, input_file)  
        return
    
    # Save this data in the (kx,ky) matrixes 
    self.stable     = True if (data["stable"]=="True") else False
    self.unstable   = True if (data["unstable"]=="True") else False  
    self.converged  = True if (data["converged"]=="True") else False  
    self.omega_last = float(data["omega_last"])
    self.omega_avg  = float(data["omega_avg"])
    self.omega_min  = float(data["omega_min"])
    self.omega_max  = float(data["omega_max"])
    self.gamma_last = float(data["gamma_last"])
    self.gamma_avg  = float(data["gamma_avg"])
    self.gamma_min  = float(data["gamma_min"])
    self.gamma_max  = float(data["gamma_max"])    
    return 
            
#===============================================================================
#                      CALCULATE AND WRITE THE LINEAR DATA                     #
#===============================================================================  

def calculate_linearData(self, input_file, verbose=False):
    ''' Calculate the linear data for a single mode. ''' 

    # Write the linear data to "*.linear_data" so open the file 
    file = configparser.ConfigParser() 
    path = input_file.with_suffix(".lineardata.ini")
    file.read(path) 
    
    # Gamma and omega can be nans 
    if len(self.omega.omega_vs_t.omega[:][np.isfinite(self.omega.omega_vs_t.omega)])==0:       
        file['('+str(self.kx)+', '+str(self.ky)+')'] = {
            'stable'     : np.nan,\
            'unstable'   : np.nan,\
            'converged'  : np.nan,\
            'gamma_max'  : np.nan,\
            'omega_last' : np.nan,\
            'gamma_last' : np.nan,\
            'omega_avg'  : np.nan,\
            'gamma_avg'  : np.nan,\
            'omega_min'  : np.nan,\
            'omega_max'  : np.nan,\
            'gamma_min'  : np.nan,\
            'gamma_max'  : np.nan} 
        
    # If no wave is examined then omega=0 and gamma=0 and we divide by 0
    elif not (np.isclose(self.kx,0) and np.isclose(self.ky,0)):  
        
        # Read (omega; gamma) to calculate the linear data                                  
        omega = self.omega.omega_vs_t.omega[:]
        gamma = self.omega.gamma_vs_t.gamma[:]  
        
        # Remove the infinite and NaN values 
        gamma = gamma[np.isfinite(omega)]
        omega = omega[np.isfinite(omega)]  
        
        # Grab omega and gamma at the last time value where omega is finite. 
        self.omega_last = omega[-1] 
        self.gamma_last = gamma[-1]
        
        # Reduce (omega; gamma) to the last x% of their data to get the stable state 
        omega = omega[int(np.size(omega)*self.percentage):]
        gamma = gamma[int(np.size(gamma)*self.percentage):] 
        
        # Average over the last x% since before it converges it osicillates around its convergent point.
        self.omega_avg = np.average(omega)
        self.gamma_avg = np.average(gamma)
        
        # Calculate the error bars over the last x%:
        self.omega_min = np.abs(np.nanmin(omega)-self.omega_avg)
        self.omega_max = np.abs(np.nanmax(omega)-self.omega_avg)
        self.gamma_min = np.abs(np.nanmin(gamma)-self.gamma_avg)
        self.gamma_max = np.abs(np.nanmax(gamma)-self.gamma_avg) 
        
        # Intiate each mode to be stable and not converged
        self.stable = True
        self.unstable = False
        self.converged = False
             
        # A mode is unstable when gamma does not change sign in the last 20% of the time trace 
        gamma = self.omega.gamma_vs_t.gamma[np.isfinite(self.omega.gamma_vs_t.gamma)] 
        last_20percent  = int(np.size(gamma)*8/10) 
        gamma_20percent = gamma[last_20percent:]
        sign_changes    = np.where(np.diff(np.sign(gamma_20percent)))[0] 
        if len(sign_changes)==0 and gamma[-1]>0:  
            self.stable = False
            self.unstable = True 
            
        # For the unstable modes, check whether the simulation has converged (constant gamma)
        # Consider the mode converged if gamma/omega doesn't change more than 2% over the last 10% of time
        omega = self.omega.omega_vs_t.omega[np.isfinite(self.omega.omega_vs_t.omega)] 
        gamma = self.omega.gamma_vs_t.gamma[np.isfinite(self.omega.gamma_vs_t.gamma)]  
        self.converged = np.all( np.isclose(omega[int(len(omega)*0.90):], omega[-1], rtol=0.02) )
        self.converged = np.all( np.isclose(gamma[int(len(gamma)*0.90):], gamma[-1], rtol=0.02) ) & self.converged
        
        # Add some information to the command prompt if the mode is not converged
        if self.converged == False and verbose: 
            fluctuation_o = (max(omega[int(len(omega)*0.90):])-min(omega[int(len(omega)*0.90):]))/omega[-1]*100
            fluctuation_g = (max(gamma[int(len(gamma)*0.90):])-min(gamma[int(len(gamma)*0.90):]))/gamma[-1]*100
            print("\nThe simulation for (kx,kx) = ("+str(self.kx)+", "+str(self.ky)+") has not converged yet:")
            print("    The relative fluctuation of omega between (0.9*t, t) is:  "+str(fluctuation_o)+"%.")
            print("    The relative fluctuation of gamma between (0.9*t, t) is:  "+str(fluctuation_g)+"%.")

        # Write the data to the "*.linear_data" file 
        file['('+str(self.kx)+', '+str(self.ky)+')'] = {
            'stable'     : self.stable,\
            'unstable'   : self.unstable,\
            'converged'  : self.converged,\
            'gamma_max'  : self.gamma_max,\
            'omega_last' : self.omega_last,\
            'gamma_last' : self.gamma_last,\
            'omega_avg'  : self.omega_avg,\
            'gamma_avg'  : self.gamma_avg,\
            'omega_min'  : self.omega_min,\
            'omega_max'  : self.omega_max,\
            'gamma_min'  : self.gamma_min,\
            'gamma_max'  : self.gamma_max} 
        
    # Write the "*.linear_data" file
    file.write(open(path, 'w'))
    return    


