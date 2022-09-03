
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#################################################################
#                   CLASS FOR THE AX OBJECT
#################################################################

# Make the axis instance and the required attributes for the options window
# Either the plotting function is a class or we manually give it some class attributes like here
class Graph:
    
    def __init__(self, tab, figure, i):
        
        # Select the current figure
        plt.figure(figure.number)
        
        # Change the margins
        self.grid_specifications = gridspec.GridSpec(1, 1)
        if i==0: self.grid_specifications.update(top=0.88, left=0.15, right=0.95, bottom=0.15)  
        if i==1: self.grid_specifications.update(top=0.88, left=0.15, right=0.95, bottom=0.15)  
        
        # Set the names for the axis variables shown in the OptionsWindow
        if i==0: self.x_name = "Time";          self.y_name = "Electrostatic potential"
        if i==1: self.x_name = "Z";             self.y_name = "Electrostatic potential" 
        
        # Create the axis object
        self.ax = plt.subplot(self.grid_specifications[0])
        
        # Set the variables for the OptionsWindow
        self.range = {"units" : "normalized", "x_scale" : "linear", "y_scale" : "linear"}    
        self.label = {"x" : None, "y" : None, "title" : None}
        self.layout = {"fontsize" : "N.A.", 'handlelength' : "N.A."}
        
        # Set the names for the axis variables shown in the OptionsWindow
        if i==0: self.x_name = "Time";          self.y_name = "Fluxes"
        if i==1: self.x_name = "Parameter";     self.y_name = "Saturated fluxes" 
        
        # Set the ranges and labels to None when the graph is reset, this will triger the defaults
        self.load_defaults()
        
        # Save the graph and id
        self.id = i
        self.tab = tab
        return 
    
    #---------------
    def load_defaults(self):       
        self.range["x"] = None
        self.range["y"] = None
        for key in self.label.keys(): 
            if key != "title":
                self.label[key] = None
        return

    #---------------
    def update_quantitiesAndKeys(self):
#         if self.id==0:
#             # TODO: Change this
#             y_quantity = self.tab.PlotLinearSpectrum.y_quantity
#             if y_quantity=="omega":         self.y_name = "Frequency" 
#             if y_quantity=="gamma":         self.y_name = "Growth rate" 
#             if y_quantity=="gamma/ky**2":   self.y_name = "Growthrate/ky**2 versus modes"
#         if self.id==1:
#             # TODO: Change this
#             self.x_name = self.tab.PlotParameterInfluence.key 
#             y_quantity  = self.tab.PlotParameterInfluence.y_quantity
#             if y_quantity=="omega":         self.y_name = "Frequency" 
#             if y_quantity=="gamma":         self.y_name = "Growth rate" 
#             if y_quantity=="gamma/ky**2":   self.y_name = "Growthrate/ky**2 versus modes"
#         if self.id==2:
#             # TODO: Change this
#             self.x_name = self.tab.PlotLinearMap.key1
#             self.y_name = self.tab.PlotLinearMap.key2 
#             z_quantity  = self.tab.PlotLinearMap.z_quantity
#             if z_quantity=="omega":         self.z_name = "Frequency" 
#             if z_quantity=="gamma":         self.z_name = "Growth rate" 
#             if z_quantity=="ky":            self.z_name = "Wavenumber"
        return
    
    #---------------
    def update_rangesAndLabels(self):

        # Get the ranges, labels and titles
        self.range["x"]     = self.ax.get_xlim()
        self.range["y"]     = self.ax.get_ylim()
        self.label["x"]     = self.ax.get_xlabel()
        self.label["y"]     = self.ax.get_ylabel()
        self.label["title"] = self.ax.get_title()

        # The titles of the option window sections
        # TODO: Change this
#         if self.var_plot.get()==1: self.y_name = "Heat flux"
#         if self.var_plot.get()==2: self.y_name = "Particle flux"
#         if self.var_plot.get()==3: self.y_name = "Momentum flux"




