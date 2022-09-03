
import numpy as np
from stellapy.data.time.get_timeFrame import get_timeFrame
 
def calculate_timeAverage(simulation, quantity, specie, t_axis=0, vec_time=None, t_range=None, t_specific=None, quant="phi"): 
    
    # Get the quantity as an average over a time frame
    if not t_specific: 
        
        # Get the correct time axis
        if vec_time is None:
            if quant=="g":  vec_time = simulation.g_time
            else:           vec_time = simulation.vec_time
        
        # Get the time frame to average over the corresponding filter
        if t_range==None: t_start, t_stop = get_timeFrame(simulation, specie)  
        if t_range!=None: t_start, t_stop = t_range[0], t_range[1]
        time_filter = (vec_time > t_start) & (vec_time < t_stop)

        # Determine the number of dimensions and the dimensions over which to iterate
        dimensions = list(np.shape(quantity)) 
        iterate_dimensions = list(range(len(dimensions)))
        iterate_dimensions.pop(t_axis) 
        
        # Time average over one dimension (which is t)
        if len(iterate_dimensions)==0:  
            quantity_timeAveraged = np.nanmean(quantity[time_filter])
            
        # Time average over two dimensions (t + one extra dimension)
        if len(iterate_dimensions)==1:    
            if t_axis==0: quantity_timeAveraged = np.nanmean(quantity[time_filter,:], axis=0)
            if t_axis==1: quantity_timeAveraged = np.nanmean(quantity[:,time_filter], axis=1)
                
        # Time average over three dimensions (t + two extra dimensions)
        if len(iterate_dimensions)==2:  
            if t_axis==0: quantity_timeAveraged = np.nanmean(quantity[time_filter,:,:], axis=0)
            if t_axis==1: quantity_timeAveraged = np.nanmean(quantity[:,time_filter,:], axis=1)
            if t_axis==2: quantity_timeAveraged = np.nanmean(quantity[:,:,time_filter], axis=2) 
                
        # Time average over four dimensions (t + three extra dimensions)
        if len(iterate_dimensions)==3:    
            if t_axis==0: quantity_timeAveraged = np.nanmean(quantity[time_filter,:,:,:], axis=0)
            if t_axis==1: quantity_timeAveraged = np.nanmean(quantity[:,time_filter,:,:], axis=1)
            if t_axis==2: quantity_timeAveraged = np.nanmean(quantity[:,:,time_filter,:], axis=2)
            if t_axis==3: quantity_timeAveraged = np.nanmean(quantity[:,:,:,time_filter], axis=3)
             
        # Time average over five dimensions (t + four extra dimensions)
        if len(iterate_dimensions)==4:    
            if t_axis==0: quantity_timeAveraged = np.nanmean(quantity[time_filter,:,:,:,:], axis=0)
            if t_axis==1: quantity_timeAveraged = np.nanmean(quantity[:,time_filter,:,:,:], axis=1)
            if t_axis==2: quantity_timeAveraged = np.nanmean(quantity[:,:,time_filter,:,:], axis=2)
            if t_axis==3: quantity_timeAveraged = np.nanmean(quantity[:,:,:,time_filter,:], axis=3)
            if t_axis==4: quantity_timeAveraged = np.nanmean(quantity[:,:,:,:,time_filter], axis=4)
            
        # Return the time average of the quantity
        return quantity_timeAveraged
            
    # Get the quantity at a specific time
    if t_specific:
        
        # Last time point
        if "last" in str(t_specific):
            index_time = len(simulation.vec_time)-1
        
        # Get the index along the time axis where we encounter <t_specific> 
        else:
            index_time = np.argwhere(simulation.vec_time > t_specific)[0][0]

        # Quantity at a specific time when there is one dimension (which is t)
        if len(list(np.shape(quantity)))==1:  
            quantity_atSpecificTime = quantity[index_time]
            
        # Quantity at a specific time when there are two dimensions (t + one extra dimension)
        if len(list(np.shape(quantity)))==2:   
            if t_axis==0: quantity_atSpecificTime = quantity[index_time,:]
            if t_axis==1: quantity_atSpecificTime = quantity[:,index_time]
                
        # Quantity at a specific time when there are three dimensions (t + two extra dimension)
        if len(list(np.shape(quantity)))==3:  
            if t_axis==0: quantity_atSpecificTime = quantity[index_time,:,:]
            if t_axis==1: quantity_atSpecificTime = quantity[:,index_time,:]
            if t_axis==2: quantity_atSpecificTime = quantity[:,:,index_time]
                
        # Quantity at a specific time when there are four dimensions (t + three extra dimension)
        if len(list(np.shape(quantity)))==4:  
            if t_axis==0: quantity_atSpecificTime = quantity[index_time,:,:,:]
            if t_axis==1: quantity_atSpecificTime = quantity[:,index_time,:,:]
            if t_axis==2: quantity_atSpecificTime = quantity[:,:,index_time,:]
            if t_axis==3: quantity_atSpecificTime = quantity[:,:,:,index_time]
             
        # Quantity at a specific time when there are five dimensions (t + four extra dimension)
        if len(list(np.shape(quantity)))==5:     
            if t_axis==0: quantity_atSpecificTime = quantity[index_time,:,:,:,:]
            if t_axis==1: quantity_atSpecificTime = quantity[:,index_time,:,:,:]
            if t_axis==2: quantity_atSpecificTime = quantity[:,:,index_time,:,:]
            if t_axis==3: quantity_atSpecificTime = quantity[:,:,:,index_time,:]
            if t_axis==4: quantity_atSpecificTime = quantity[:,:,:,:,index_time]
 
        # Return the field line average of the quantity
        return quantity_atSpecificTime

