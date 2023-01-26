 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib; matplotlib.use("TkAgg")   
from stellapy.plot.utils.labels.standardNames import standardNames 
from stellapy.plot.utils.labels.standardHeaders import standardHeaders
 
################################################################################
#                          CLASS TO CREATE CANVASSES
################################################################################
    
class Axis: 
    ''' The class <Axis> is attached to the class <Canvas> and holds the <ax>
    object of the matplotlib <figure> that is displayed on the <Canvas>. All
    data related to the axis and figure is attached to this object. 
    
    
    Attributes
    ----------
    figure : matplotlib.pyplot instance
        The figure object that is displayed on the <Canvas>
    ax : Matplotlib.axes.Axes instance
        The axis object of the <figure>
    '''
    
    def __init__(self, Canvas, root, identifier):
        ''' Create the <ax> object, and define its (top, bottom, left, right) 
        margins through <grid_specifications>. Set default values for the (x,y) 
        labels, the (x,y) ranges and the (x,y) variables, which are needed for 
        the <OptionsWindow>. '''
        
        # Remember the root, Canvas and the figure identifier
        self.root = root
        self.Canvas = Canvas
        self.identifier = identifier

        # Create a matplotlib figure and style it 
        self.initiate_figure()
        
        # Create the axis object
        self.initiate_axis()
        
        # Set the standard ranges/labels for this graph and the options window
        self.initiate_axisInformation()  
        
        # Initiate the (kx,ky) ranges
        self.kx_range = 0.0             
        self.ky_range = [0,100]     
        
        # Initiate the empty dictionaries
        self.axisQuantities = {}; 
        self.plot_id = 0
        return 
    

    #----------------------------
    def initiate_figure(self):
        
        # Create the matplotlib figure 
        self.figure = plt.figure(self.identifier)   

        # Style the figure  
        self.figure.set_tight_layout(False)  
        self.figure.patch.set_facecolor(self.root.color['canvas']) 
        
        # Set the color of the axis
        plt.rcParams['text.color']       = self.root.color['fg']
        plt.rcParams['axes.edgecolor']   = self.root.color['fg']
        plt.rcParams['axes.labelcolor']  = self.root.color['fg']
        plt.rcParams['xtick.color']      = self.root.color['fg']
        plt.rcParams['ytick.color']      = self.root.color['fg']
        plt.rcParams['axes.facecolor']   = self.root.color['canvas']
        plt.rcParams['figure.facecolor'] = self.root.color['canvas']
        plt.rcParams['savefig.facecolor']= "white"
        return

    #----------------------------
    def initiate_axis(self):
        
        # Create a <grid_specifications> object to set the margins easily
        self.grid_specifications = gridspec.GridSpec(1, 1, figure=self.figure)
        self.grid_specifications.update(top=0.95, left=0.08, right=0.95, bottom=0.1)
        
        # Create the axis object
        self.ax = plt.subplot(self.grid_specifications[0])
        
    #----------------------------
    def initiate_axisInformation(self):
        ''' Initiate all the information that is related to the axis and plot.
        After plotting a figure, the label and ranges information is attached to
        the <Axis> object, which allows one to re-use these settings when the 
        graph gets plotted again. Since one can manually set the ranges/labels
        through the optionswindow and we don't want these changes to be removed
        when the graph gets plotted again. When we want to forget the ranges and
        labels that were set manually, execute this function again. ''' 
        # Data ranges
        self.x_range  = None
        self.y_range  = None
        self.z_range  = None 
        # Axis labels      
        self.x_label = None
        self.y_label = None  
        self.z_label = None  
        self.title   = None
        # Other
        self.x_scale = "linear"    
        self.y_scale = "linear"  
        self.z_scale = "linear"
        self.units = "normalized"
        # Layout
        self.fontsize = 20
        self.handlelength = 1
        # Plotting and research arguments
        self.plotting_arguments = None
        self.research_arguments = None
        # Clear the axis
        self.ax.clear()
        self.root.Progress.show_loadingCursor("nothing") 

################################################################################
#                   METHODS FOR THE AXIS AND OPTIONSWINDOW
################################################################################

    def reset_axis(self): 
        self.initiate_axisInformation() 
        return

    #----------------------------
    def update_gridSpeficications(self, top=0.95, left=0.08, right=0.95, bottom=0.1):
        self.grid_specifications.update(top=top, left=left, right=right, bottom=bottom)
        return
    
    #---------------------------
    def set_axisQuantities(self, plot_id, x_quantity, y_quantity, z_quantity=None, name=None):
        ''' For <plot_id> define which quantities are plotted on the (x,y) axis. '''
        self.axisQuantities[plot_id] = {'x' : x_quantity, 'y' : y_quantity, 'z' : z_quantity, 'name' : name}
        return

    #---------------------------
    def get_axisQuantities(self):
        ''' For <plot_id> get the defined quantities that are plotted on the (x,y) axis. '''
        return list(self.axisQuantities.keys()), self.axisQuantities
        return
    
    #---------------------------
    def get_plotid(self):
        ''' When the plot is changed, make sure the Axis object knows. '''
        return self.plot_id
    
    #---------------------------
    def get_standardHeader(self):
        ''' For a popped out window or for a labelFrame we can set a custom 
        title to clarify what is being plotted. '''
        return standardHeaders[self.y_quantity][self.x_quantity] 
        
    #---------------------------
    def set_plot(self, plot_id):
        ''' When the plot is changed, make sure the Axis object knows. '''
        self.plot_id = plot_id
        self.update_axisNames()
        self.update_axisQuantities()
        return 
    
    #---------------------------
    def update_labelFrame(self, frame_canvas):
        ''' The frame of a <Canvas> can be a label frame, then we can 
        change the label of it according to the graph. '''   
        if self.z_quantity==None: y = self.y_quantity 
        if self.z_quantity!=None: y = self.z_quantity
        try: text = standardHeaders[y][self.x_quantity] 
        except: print("NO STANDARD HEADER SET FOR:", y, " vs ", self.x_quantity)
        try: frame_canvas.config(text="   "+text+"  ")
        except: pass
        if True: return 

        
##########################
# METHODS BEFORE PLOTTING
##########################

    def start_plotting(self, research_arguments, input_files): 
        self.check_whetherItIsTheSameResearch(research_arguments, input_files)
        self.root.Progress.show_loadingCursor("start")
        self.ax.clear()
        return
    
    #---------------------------
    def check_whetherItIsTheSameResearch(self, research_arguments, input_files):
        ''' If the research object has changed, load the default axis data. '''
        if self.research_arguments != research_arguments: self.reset_axis()
        elif self.input_files != input_files: self.reset_axis()
        if True: return 

        
##########################
# METHODS AFTER PLOTTING
##########################

    def finish_plotting(self, research_arguments, input_files):
        self.update_axisInformation()
        self.save_currentPlottedData(research_arguments, input_files)
        self.Canvas.draw_idle()
        self.root.Progress.show_loadingCursor("finished")
        return

    #---------------------------
    def save_currentPlottedData(self, research_arguments, input_files):
        ''' Save the plotting and research arguments in order to know whether 
        during a next plotting call, the requested plot is already displayed, 
        this can prevent having to replot. Copying a dictionary will make a 
        shallow copy, meaning that changes to the old dictionary will not affect
        the new copy, however, changes to lists/classes/objects in the old 
        dictionary will also change in the new dictionary since the dictionary
        holds a reference to lists/classes/objects, not the actual data. Thus, 
        save the input file with a deep copy since it is a list.'''   
        self.research_arguments = research_arguments.copy()  
        self.input_files = input_files.copy()  
        return
    
    #---------------------------
    def save_plottingArguments(self, plotting_arguments): 
        ''' See previous function. '''
        self.plotting_arguments = plotting_arguments.copy()  
        return
        
    #---------------------------
    def update_axisQuantities(self):
        ''' A canvas/figure/axis can be reused to plot different (x,y) variables,
        therefore, we need to update the quantities for the plot. '''
        self.x_quantity = self.axisQuantities[self.plot_id]['x']
        self.y_quantity = self.axisQuantities[self.plot_id]['y']  
        self.z_quantity = self.axisQuantities[self.plot_id]['z']  
        return

    #---------------------------
    def update_axisNames(self):
        ''' A canvas/figure/axis can be reused to plot different (x,y) variables,
        therefore, we need to update the names for the OptionsWindow. '''
        self.x_name = standardNames[self.axisQuantities[self.plot_id]['x']]
        self.y_name = standardNames[self.axisQuantities[self.plot_id]['y']]  
        self.z_name = standardNames[self.axisQuantities[self.plot_id]['z']]  
    
    #---------------------------
    def update_axisInformation(self):
        ''' After plotting a figure, attach the label and ranges information 
        to the <Axis> object, this allows us to re-use these settings when the 
        graph gets plotted again. Since one can manually set the ranges/labels
        through the optionswindow and we don't want these changes to be removed
        when the graph gets plotted again. '''
        self.x_range  = self.ax.get_xlim()
        self.y_range  = self.ax.get_ylim()
        self.x_label  = self.ax.get_xlabel()
        self.y_label  = self.ax.get_ylabel()
        self.title    = self.ax.get_title()            
        return
    
    #---------------------------
    def copy_dataFromOtherAxis(self, Axis):
        ''' When creating a popped out window, we want an exact copy of the 
        original figure, therefore copy all the axis and optionswindow data.'''
        self.kx_range       = Axis.kx_range
        self.ky_range       = Axis.ky_range
        self.x_range        = Axis.x_range
        self.y_range        = Axis.y_range
        self.z_range        = Axis.z_range
        self.x_label        = Axis.x_label
        self.y_label        = Axis.y_label
        self.z_label        = Axis.z_label
        self.title          = Axis.title 
        self.x_scale        = Axis.x_scale    
        self.y_scale        = Axis.y_scale  
        self.z_scale        = Axis.z_scale
        self.units          = Axis.units
        self.fontsize       = Axis.fontsize
        self.handlelength   = Axis.handlelength
        self.plot_id        = Axis.plot_id
        self.axisQuantities = Axis.axisQuantities 
        self.update_axisQuantities()
        self.update_axisNames()
        return
        
        
        
        
        
        
        
        
        
        
    

