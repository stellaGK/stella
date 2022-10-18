

from matplotlib.lines import Line2D
from stellapy.plot.utils.labels.standardLabels import standardLabels   

class Legend:
    
    def __init__(self, ax, plot): 
        self.labels_points = []
        self.labels_lines = []
        self.plot = plot
        self.loc = plot.loc
        self.shadow = plot.shadow
        self.ax = ax 
        return
    
    #---------------------------------
    def add_legend(self):
        ''' Add the legend to the plot. ''' 
        if not self.plot.removelegend: 
            if len(self.plot.y_quantities)==1: self.add_legendForOneYAxis()
            if len(self.plot.y_quantities)>1:  self.add_legendForMultipleYAxis()
    
    #---------------------------------
    def add_legendForOneYAxis(self): 
        
        # First sort the legend labels
        labels, handles, emptyLegend = self.sort_legend()
        
        # Legend location 
        if self.loc==None:
            if self.plot.splitInKineticAdiabatic==False: self.loc = 'best'
            if self.plot.splitInKineticAdiabatic==True:  self.loc = 'upper left'
        
        # Only add a legend if it's usefull
        if not emptyLegend:  
            legend = self.ax.legend(handles, labels, \
                       labelspacing=0.0,\
                       shadow=self.shadow,\
                       loc=self.loc,\
                       prop={'size':self.plot.fontsize},\
                       handlelength=self.plot.handlelength,\
                       ncol=self.plot.ncol)

        # Add the difference between kinetic/adiabatic simulations
        if self.plot.splitInKineticAdiabatic: 
            loc = 'upper right'; loc = None
            if self.plot.adiabaticIons==True and self.plot.adiabaticElectrons==True:
                lines = [Line2D([0], [0], color="black", linewidth=3, linestyle=s) for s in ["-", ":", "--"]]
                labels = ['Kinetic electrons', 'Adiabatic electrons', 'Adiabatic ions']
            if self.plot.adiabaticIons==True and self.plot.adiabaticElectrons==False:
                lines = [Line2D([0], [0], color="black", linewidth=3, linestyle=s) for s in ["-", "--"]]
                labels = ['Kinetic electrons', 'Adiabatic ions']
            if self.plot.adiabaticIons==False and self.plot.adiabaticElectrons==True:
                lines = [Line2D([0], [0], color="black", linewidth=3, linestyle=s) for s in ["-", ":"]]
                labels = ['Kinetic electrons', 'Adiabatic electrons']
            _ = self.ax.legend(
               lines, labels,
               labelspacing=0.0,\
               shadow=True,\
               loc=loc,\
               prop={'size':self.plot.fontsize},\
               handlelength=self.plot.handlelength)  
        
        # Add the legend
        if not emptyLegend: self.ax.add_artist(legend) 
        return
        
    #---------------------------------
    def add_legendForMultipleYAxis(self): 

        # Define the colors and labels of the plot
        colors = ["black", "red", "green"]
        labels = standardLabels[self.plot.units]
        labels = [labels[key] for key in self.plot.y_quantities]             
        
        # Change the labels and colors
        lines = self.ax.get_lines()
        for i in range(len(self.plot.y_quantities)):   
            lines[i].set_c(colors[i])    
            lines[i].set_label(labels[i]) 
            lines[i].set_marker(None)
            
        # Make and add the legend 
        self.ax.legend(
           labelspacing=0.0,\
           shadow=True,\
           prop={'size':self.plot.fontsize},\
           handlelength=self.plot.handlelength,\
           ncol=self.plot.ncol)    
        return 
        
    #---------------------------------
    def sort_legend(self):
        
        # Grab the handles and legend of the default legend
        handles, labels = self.ax.get_legend_handles_labels()    
        
        # Sort the labels of the markers by the numbers in them
        from stellapy.utils.files.sort_listByNumbers import sort_listByNumbers 
        self.labels_points = sort_listByNumbers(self.labels_points)  
        
        # Remove the labels of points if they are already a label of the line
        self.labels_points = [label for label in self.labels_points if label not in self.labels_lines]
        
        # We want to first show the labels of the lines, and then those of the markers
        labels_order = self.labels_lines + self.labels_points   
        
        # Make sure we only look at labels that are actually plotted 
        labels_order = [label for label in labels_order if label in labels]
        
        # Sort <labels> and <handles> to match the order in <labels_order>
        labels_index = [labels.index(label) for label in labels_order] 
        labels  = [labels[i] for i in labels_index]
        handles = [handles[i] for i in labels_index] 
        
        # Check whether the legend would be empty or pointless
        emptyLegend = (len(labels_order)) <= 1 or (len(self.labels_lines)==0 and len(self.labels_points)==0) 
        return labels, handles, emptyLegend
 
