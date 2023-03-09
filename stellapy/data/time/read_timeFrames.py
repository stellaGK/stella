
#!/usr/bin/python3  
import pathlib
import sys, os
import numpy as np
import configparser
from datetime import datetime   

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.utils.decorators.exit_program import exit_program

#===============================================================================
#                              READ THE TIME FRAMES                            #
#=============================================================================== 

def read_timeFrames(self): 

    # Write the file if it doesn't exist
    if not os.path.isfile(self.path.folder / "timeFrames.ini"):  
        self.file = write_txtFimeForTimeFrames(self.path)
    
    # Read the data from the "timeFrames" file 
    elif os.path.isfile(self.path.folder / "timeFrames.ini"):   
        self.file = configparser.ConfigParser()  
        self.file.read(self.path.folder / "timeFrames.ini")
        
    # Read the section related to the selected time frame
    self.timeframe = self.file["General"]["timeframe"].replace("%","%%")  
    self.section = self.file[self.timeframe]
    
    # Get tstart and tend
    try: self.tstart, self.tend = get_tstartAndTend(self.timeframe, self.fluxes.qflux_vs_ts.t)
    except: 
        exit_reason = "The time frame in the timeFrames.ini file doesn't work, \n"
        exit_reason += "correct it or remove the timeFrames.ini file inside: \n"
        exit_reason += str(self.path.folder)
        exit_program(exit_reason, read_timeFrames, sys._getframe().f_lineno)
    
#===============================================================================
#                             GET TSTART AND TEND                              #
#=============================================================================== 

def get_tstartAndTend(timeframe, vec_time):
    
    # Get [tstart, tend] from the "timeframe" file  
    tstart, tend = timeframe.replace(" ","").replace("%%","%").split("(")[-1].split(")")[0].split(",")
    
    # Make sure the times represent an actual time in vec_time
    tstart = get_timePoint(tstart, vec_time)
    tend = get_timePoint(tend, vec_time)
    
    # The end time can not be bigger than tlast and tstart can not be bigger than tend
    tend = np.nanmin([tend, vec_time[-1]])  
    if tstart >= tend: tstart = tend-10 
    return tstart, tend

#-----------------------------
def get_timePoint(time, vec_time):
    if time=="-1":     time = vec_time[-1]
    elif time=="100%": time = vec_time[-1]
    elif time=="0%":   time = vec_time[0] 
    elif "%" in time:  time = vec_time[-1]*float(time.replace("%",""))/100
    time = float(time) 
    return time 

#===============================================================================
#                             UPDATE THE TIME FRAME                            #
#=============================================================================== 

def update_timeFrame(self, time_range): 
    
    # Read the timeframes file and add a new section for the selected time frame 
    read_timeFrames(self)
    timeframe = "("+str(time_range[0])+", "+str(time_range[1])+")" 
    self.file["General"]["timeframe"] = timeframe
    if timeframe not in self.file:
        self.file[timeframe] = {"peakdate" : datetime.now(),\
                           "peak" : "/",\
                           "fluxdate" : datetime.now(),\
                           # Fluxes
                           "qflux" : "/", "pflux" : "/", "vflux" : "/",\
                           "qfluxstd" : "/", "pfluxstd" : "/", "vfluxstd" : "/",\
                           "qfluxmin" : "/", "pfluxmin" : "/", "vfluxmin" : "/",\
                           "qfluxmax" : "/", "pfluxmax" : "/", "vfluxmax" : "/",\
                           # Potential
                           "phi2" : "/", "phi2zonal" : "/", "phi2nozonal" : "/",\
                           "phi2std" : "/", "phi2zonalstd" : "/", "phi2nozonalstd" : "/",\
                           "phi2max" : "/", "phi2zonalmax" : "/", "phi2nozonalmax" : "/",\
                           "phi2min" : "/", "phi2zonalmin" : "/", "phi2nozonalmin" : "/",\
                           # Distribution
                           "ge" : "/", "gi" : "/",\
                           "gestd" : "/", "gistd" : "/",\
                           "gemax" : "/", "gimax" : "/",\
                           "gemin" : "/", "gimin" : "/"} 
    self.file.write(open(self.path.folder/"timeFrames.ini", 'w'))  
    
    # Re-initiate the time class to make sure we use the new time frame
    from stellapy.data.time.load_timeObject import load_timeObject
    load_timeObject(self.simulation)
    return

#===============================================================================
#                     WRITE A TEXT FILE FOR THE TIME FRAMES                    #
#=============================================================================== 

def write_txtFimeForTimeFrames(path): 
    
    # Initiate the "timeFrames" file 
    file = configparser.ConfigParser()  
    file.read(path.folder / "timeFrames.ini") 
    
    # Put some basic information inside 
    file["General"] = {"timeframe" : "(50%%, 100%%)"}
    file["(50%%, 100%%)"] = {"peakdate" : datetime.now(), "fluxdate" : datetime.now(), "peak" : "/", "qflux" : "/", "pflux" : "/", "vflux" : "/"} 
        
    # Write the the "timeFrames" file 
    file.write(open(path.folder / "timeFrames.ini", 'w'))
    return file