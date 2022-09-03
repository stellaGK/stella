
import tkinter as tk
from tkinter import ttk  
from stellapy.GUI.widgets import Canvas, PAD_TITLE2, PAD_LABEL2 
from stellapy.GUI.plot.linearSimulations.plot_velocityDistribution import plot_velocityDistribution

################################################################################
#        CLASS TO PLOT THE VELOCITY DISTRIBUTION OF LINEAR SIMULATIONS
################################################################################
   
class PlotVelocityDistribution: 
    ''' Holds and manipulates two Canvas and Axis objects for the following plots:
            (1) ... '''
    
    
    def __init__(self, tab, notebook_tab):

        # Set an identifier for this plotting class, used as the figure name
        self.identifier = "TabLinear:PlotVelocityDistribution"
        
        # Make the parent objects available
        self.tab = tab
        self.root = tab.root  
        
        # Remember whether we have a Canvas yet
        self.initiated_canvas = False 

        # Create the frames and variables for the plotting arguments
        self.initiate_optionVariablesAndFrames(notebook_tab)
        
        # When the tab is visible, make sure its <frame_canvas> is displayed
        notebook_tab.bind("<Visibility>", lambda _, i=self.identifier: self.tab.display_canvas(i)) 
        return 

################################################################################
#                            CLASSES FOR THE PLOTS
################################################################################

    def initiate_canvas(self, frame1, frame2):
        ''' The <Canvas> is not generated until its parent frame is visible, this
        saves GUI loading time, since only the visisble Canvasses are created. 
        Create the Canvas of this plot in the given frame, configure and save it. '''
        
        # Initiate the <Canvas>
        if self.initiated_canvas==False: 
            
            # Create the <Canvas> and remember also the <Axis>
            self.Canvas1 = Canvas(self.root, frame1, self.identifier+"_vpa")
            self.Canvas2 = Canvas(self.root, frame2, self.identifier+"_mu") 
            self.Axis1 = self.Canvas1.Axis; self.Canvas1.CanvasII = self.Canvas2 
            self.Axis2 = self.Canvas2.Axis; self.Canvas2.CanvasII = self.Canvas1
    
            # Attach this plotting class to the canvas
            self.Canvas1.add_plottingClass(self)
            self.Canvas2.add_plottingClass(self)
            
            # Resize the <matplotlib figure>
            self.Axis1.update_gridSpeficications(top=0.94, left=0.1, right=0.95, bottom=0.12)
            self.Axis2.update_gridSpeficications(top=0.94, left=0.1, right=0.95, bottom=0.12)
            
            # Set the (x,y) quantities and titles for the canvas
            self.Axis1.set_axisQuantities(1, "vpa", "g"); self.Axis1.set_plot(1)
            self.Axis2.set_axisQuantities(1, "mu", "g");  self.Axis2.set_plot(1) 
            self.Axis1.update_labelFrame(self.Canvas1.frame)
            self.Axis2.update_labelFrame(self.Canvas2.frame)
    
            # Remember that the <Canvas> is initiaited
            self.initiated_canvas = True 
            
            # Plot vpa left and mu right
            self.Canvas1.x_quantity = "vpa"
            self.Canvas2.x_quantity = "mu"
        return     

    #---------------------- 
    def get_argumentsForPlot(self, Canvas):
        ''' Gather the arguments for <plot_quantityVsParameter>. '''
        Canvas.Axis.update_axisQuantities() 
        plotting_arguments = {
            # Specify which simulations to plot
                "research"      : self.root.Research,\
                "experiment_id" : self.tab.Simulations.experiment_id,\
                "simulation_id" : self.tab.Simulations.simulation_id,\
                "species"       : self.species,\
            # Specify data range
                "x_quantity"    : Canvas.x_quantity,\
                "x_range"       : Canvas.Axis.x_range,\
                "y_range"       : Canvas.Axis.y_range,\
                "units"         : Canvas.Axis.units,\
            # Labels
                "x_label"       : Canvas.Axis.x_label,\
                "y_label"       : Canvas.Axis.y_label,\
                "title"         : Canvas.Axis.title,\
            # Figure 
                "ax"            : Canvas.Axis.ax,\
                "show_figure"   : False,\
                "Progress"      : self.tab.Progress}   
        return plotting_arguments

    #---------------------- 
    def plotOnTheGUI(self):
        ''' When we plot on the GUI, plot on the attached <Canvas>. '''
        if self.initiated_canvas: self.plot(self.Canvas1)
        if self.initiated_canvas: self.plot(self.Canvas2)
        return
    
    #----------------------  
    def plot(self, Canvas):
        ''' Plot on the given <Canvas>. '''
        
        # Get the research arguments, plotting arguments and input_files
        research_arguments = self.root.research_arguments
        plotting_arguments = self.get_argumentsForPlot(Canvas)
        input_files = self.root.Research.input_files
        
        # Only plot if there are input_files
        if input_files != []:

            # Plot if the plotting arguments do not match those of the <Canvas> 
            if plotting_arguments != Canvas.Axis.plotting_arguments:
                    
                # Plot the parallel mode structure
                Canvas.Axis.start_plotting(research_arguments, input_files)
                plot_velocityDistribution(**plotting_arguments)
                Canvas.Axis.finish_plotting(research_arguments, input_files)
                Canvas.Axis.save_plottingArguments(self.get_argumentsForPlot(Canvas))
    
        # If no simulations have been selected, clear the current figure
        if input_files == []: Canvas.Axis.reset_axis() 
        return   
    

################################################################################
#                      METHODS AND CLASSES FOR THE OPTIONS
################################################################################
 
    def initiate_optionVariablesAndFrames(self, frame):
        ''' Enable the user to switch between the different plotting options
        of <plot_quantityVsParameter>. Create the variables, the tkinter 
        variables and the frames to display the options. These widgets are added
        to the <frame> which is a tab of a <ttk.Notebook> widget on <tabLinear>.'''
        
        # Attributes for the plot
        self.species = [0]
        
        # Create a subframe for each set of plotting options 
        self.subframe_options    = ttk.Frame(frame)
        
        # Configure the rows/columns of the <frame> to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1)      
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=1)       
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_options

        # Add the subframes to the main frame 
        self.subframe_options.grid(row=1, column=0, padx=25, pady=(0,0), stick='NSEW')
        
        # Fill the frames with widgets 
        self.Options = self.Frame_Options(self, self.subframe_options)
        return  
    

    #---------------------- 
    class Frame_Options:
        
        def __init__(self, PlotVelocityDistribution, frame):

            # Make the parents available
            self.PlotVelocityDistribution = PlotVelocityDistribution  
            
            # Variables 
            self.var_species  = tk.IntVar(value=1) 
              
            # Create the widgets             
            self.lbl_species = ttk.Label(frame, text="Species:", style='prefTitle.TLabel')
            self.rbn_ions = ttk.Radiobutton(frame, text=' Ions') 
            self.rbn_electrons = ttk.Radiobutton(frame, text=' Electrons') 
            
            # Add the commands
            self.rbn_ions.config(value=1, variable=self.var_species, command=self.change_species)
            self.rbn_electrons.config(value=2, variable=self.var_species, command=self.change_species) 
            
            # Configure the frame
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Add the options to the frame
            self.lbl_species.grid(  row=0, column=0, **PAD_TITLE2)
            self.rbn_ions.grid(     row=1, column=0, **PAD_LABEL2)
            self.rbn_electrons.grid(row=2, column=0, **PAD_LABEL2) 
        
        #---------------------- 
        def change_species(self, *_): 
            if self.var_species.get()==1: self.PlotVelocityDistribution.species = [0]
            if self.var_species.get()==2: self.PlotVelocityDistribution.species = [1]
            self.PlotVelocityDistribution.Axis1.reset_axis()
            self.PlotVelocityDistribution.Axis2.reset_axis()
            self.PlotVelocityDistribution.plotOnTheGUI() 
    
    
    
    
    
    
    
    
    