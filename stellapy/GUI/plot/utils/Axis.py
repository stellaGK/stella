
import numpy as np

class Axis:
    
    def __init__(self, ax, plot, overshoot_x=1.0, overshoot_y=1.1,
            xbot=None, xtop=None, ybot=None, ytop=None, xbot_pos=None, ybot_pos=None, xtop_neg=None, ytop_neg=None, percentage=None):
        
        # Save the data
        self.ax = ax
        self.logx = plot.logx
        self.logy = plot.logy
        self.x_range = plot.x_range
        self.y_range = plot.y_range
        self.tick_style = plot.tick_style
        self.overshoot_x = overshoot_x
        self.overshoot_y = overshoot_y
        
        # Set the limits
        self.set_limits()        
        self.percentage = percentage
        self.xbot = xbot
        self.xtop = xtop
        self.ybot = ybot
        self.ytop = ytop
        self.xbot_pos = xbot_pos
        self.ybot_pos = ybot_pos
        self.xtop_neg = xtop_neg
        self.ytop_neg = ytop_neg
        
        # Initiate the calculated limits
        self.xlims = [np.nan, np.nan]
        self.ylims = [np.nan, np.nan] 
        return
    
    #--------------------------     
    def set_limits(self, xbot=None, xtop=None, ybot=None, ytop=None, xbot_pos=None, ybot_pos=None, xtop_neg=None, ytop_neg=None, percentage=None):
        self.percentage = percentage
        self.xbot = xbot
        self.xtop = xtop
        self.ybot = ybot
        self.ytop = ytop
        self.xbot_pos = xbot_pos
        self.ybot_pos = ybot_pos
        self.xtop_neg = xtop_neg
        self.ytop_neg = ytop_neg
        return
    
    #--------------------------     
    def update_axisLimits(self, x, y): 
        
        # Grab the arrays
        x = np.array(x)
        y = np.array(y)   
        
        # Grab the limit of the y-axis based on the last part
        if self.percentage!=None:
            y = y[int(np.size(y)*self.percentage):] 
        
        # Calculate the limits
        if not np.isnan(y).all():    
            if self.logx==True: x = abs(x[x != 0])
            if self.logy==True: y = abs(y[y != 0])
            xmin = self.overshoot_x*np.nanmin(x) if np.nanmin(x)<0 else np.nanmin(x)
            ymin = self.overshoot_y*np.nanmin(y) if np.nanmin(y)<0 else np.nanmin(y)
            xmax = self.overshoot_x*np.nanmax(x) if np.nanmax(x)>0 else np.nanmax(x)
            ymax = self.overshoot_y*np.nanmax(y) if np.nanmax(y)>0 else np.nanmax(y)
            self.xlims[0] = np.nanmin([xmin, self.xlims[0]])
            self.xlims[1] = np.nanmax([xmax, self.xlims[1]]) 
            self.ylims[0] = np.nanmin([ymin, self.ylims[0]])
            self.ylims[1] = np.nanmax([ymax, self.ylims[1]])   
    
    #--------------------------       
    def rescale_axis(self):
        
        # Round the time axis to a round number
        if 50 < self.xlims[1] < 1000: self.xlims[1] = np.round(self.xlims[1],0)
        
        # First add an autoscale
        self.ax.autoscale() 
        
        # Depending on the tick style, we change the tick format
        if self.tick_style=='sci': 
            try: self.ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            except: pass
        if self.tick_style==None:
            try: self.ax.ticklabel_format(useOffset=False)
            except: pass
        
        # Change the axis to logarithmic scales 
        if self.logx == True: self.ax.set_xscale('log')
        if self.logy == True: self.ax.set_yscale('log')
        
        # If no range was manually set, use the calculated ranges for x and y
        if self.x_range==None: self.x_range = self.xlims
        if self.y_range==None: self.y_range = self.ylims
        
        # Since we use indexing, make sure we have a list and not a tuple
        self.x_range = list(self.x_range)
        self.y_range = list(self.y_range)
        
        # Make sure the lower and upper limits are not equal
        if self.x_range[0]==self.x_range[1]: self.x_range[1]=self.x_range[0]+1
        if self.y_range[0]==self.y_range[1]: self.y_range[1]=self.y_range[0]+1
        
        # Manually set the ranges  
        if self.xbot!=None: self.x_range[0] = self.xbot
        if self.xtop!=None: self.x_range[1] = self.xtop
        if self.ybot!=None: self.y_range[0] = self.ybot
        if self.ytop!=None: self.y_range[1] = self.ytop
        if self.xbot_pos!=None: self.x_range[0] = self.xbot_pos if self.x_range[0]>0 else self.x_range[0]
        if self.ybot_pos!=None: self.y_range[0] = self.ybot_pos if self.y_range[0]>0 else self.y_range[0]
        if self.xtop_neg!=None: self.x_range[1] = self.xtop_neg if self.x_range[1]<0 else self.x_range[1]
        if self.ytop_neg!=None: self.y_range[1] = self.ytop_neg if self.y_range[1]<0 else self.y_range[1]

        # Change the limits of the axis 
        if not np.isnan(self.x_range[0]) and not np.isnan(self.x_range[1]):
            self.ax.set_xlim(self.x_range)  
        if not np.isnan(self.y_range[0]) and not np.isnan(self.y_range[1]):
            self.ax.set_ylim(self.y_range)  
        return

 
