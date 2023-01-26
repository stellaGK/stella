
def get_timeFrameString(tstart, tend):
    ''' Input tstart = [min, max] and tend = [min, max] and return a string
    "[tstart:tend]" to indicate over which time frame the quantitiy was averaged. 
    Alternatively "[tstartmin:tstartmax, tendmin:tendmax]" is given when multiple 
    simulations with different time frames are plotted. Here tstart and tend are 
    rounded to the closest multiple of 10. '''
    
    # Round the number to the closest multiple of 10
    tstart[0] = int(round(tstart[0]/10,0)*10)
    tstart[1] = int(round(tstart[1]/10,0)*10)
    tend[0] = int(round(tend[0]/10,0)*10)
    tend[1] = int(round(tend[1]/10,0)*10)
    
    # If tendmin==tendmax then str=tend otherwise str=tendmin:tendmax 
    tend = str(tend[0]) if tend[0]==tend[1] else str(tend[0])+":"+str(tend[1])
    tstart = str(tstart[0]) if tstart[0]==tstart[1] else str(tstart[0])+":"+str(tstart[1])
    
    # Return the string of the time frame
    return "["+tstart+", "+tend+"]"
    