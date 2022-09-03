
import tkinter as tk
from tkinter import ttk 
from stellapy.GUI.widgets import Canvas, PAD_TITLE2, PAD_LABEL2
from stellapy.GUI.plot.nonlinearSimulations.plot_quantityVsTime import plot_quantityVsTime 

################################################################################
#           CLASS TO PLOT THE TIME EVOLUTION OF NONLINEAR QUANTITIES
################################################################################
 
class PlotTimeEvolution: 
    ''' Holds and manipulates the Canvas and Axis object for the following plots:
            (1) ... '''

    def __init__(self, tab, notebook_tab): 

        # Set an identifier for this plotting class, used as the figure name
        self.identifier = "TabNonlinearSpatial:PlotTimeEvolution"
        
        # Make the parent objects available
        self.tab = tab
        self.root = tab.root  
        
        # Remember whether we have a <Canvas> yet
        self.initiated_canvas = False
        
        # When the tab is visible, make sure its <frame_canvas> is displayed
        notebook_tab.bind("<Visibility>", lambda _, i=self.identifier: self.tab.display_canvas(i))  
        return 

################################################################################
#                            CLASSES FOR THE PLOTS
################################################################################

    def initiate_canvas(self, frame, notebook_tab):
        ''' The <Canvas> is not generated until its parent frame is visible, this
        saves GUI loading time, since only the visisble Canvasses are created. 
        Create the Canvas of this plot in the given frame, configure and save it. '''
        
        # Initiate the <Canvas>
        if self.initiated_canvas==False: 
            
            # Create the <Canvas> and remember also the <Axis>
            self.Canvas = Canvas(self.root, frame, self.identifier)
            self.Axis = self.Canvas.Axis
    
            # Attach this plotting class to the canvas
            self.Canvas.add_plottingClass(self)

            # Resize the <matplotlib figure>
            self.Axis.update_gridSpeficications(top=0.94, left=0.1, right=0.95, bottom=0.15)
            
            # Set the (x,y) quantities and titles for the canvas
            for i in range(len(self.tab.Simulations.options_yquantityKeys)): 
                self.Axis.set_axisQuantities(i, "t", self.tab.Simulations.options_yquantityKeys[i], name = self.tab.Simulations.options_yquantityNames[i])

            # Create the frames and variables for the plotting arguments
            self.initiate_optionVariablesAndFrames(notebook_tab)
             
            # Tell the <Axis> which quantity to plot, and update the label frame
            self.Axis.set_plot(self.tab.Simulations.options_yquantityNames.index(self.tab.Simulations.var_quantity.get()))
            self.Axis.update_labelFrame(self.Canvas.frame)
        
            # Remember that the <Canvas> is initiaited
            self.initiated_canvas = True
            
        # Bind <ctrl+s> to replot the graphs
        self.root.bind('<Control-s>', self.Canvas.reset_graph)
        return
    
    #---------------------- 
    def get_argumentsForPlot(self, Canvas):
        ''' Gather the arguments for <plot_quantityVsZ>. '''
        self.Axis.update_axisQuantities()
        plotting_arguments = {
            # Specify which simulations to plot
                "research"      : self.root.Research,\
                "experiment_id" : self.tab.Simulations.experiment_id,\
                "simulation_id" : self.tab.Simulations.simulation_id,\
                "species"       : self.tab.Simulations.species,\
            # Speficy which data to plot 
                "x_quantity"    : Canvas.Axis.x_quantity,\
                "y_quantity"    : Canvas.Axis.y_quantity,\
            # Specify data range
                "x_range"       : Canvas.Axis.x_range,\
                "y_range"       : Canvas.Axis.y_range,\
                "units"         : Canvas.Axis.units,\
            # Labels
                "x_label"       : Canvas.Axis.x_label,\
                "y_label"       : Canvas.Axis.y_label,\
                "title"         : Canvas.Axis.title,\
            # Options 
                "normalize"     : self.normalizeBySaturation,\
                "show_timeFrame": self.show_timeFrame,\
                "show_figure"   : False,\
            # For the GUI the figure object already exists 
                "ax"            : Canvas.Axis.ax,\
                "Progress"      : self.tab.Progress}
        return plotting_arguments
    
    #---------------------- 
    def plotOnTheGUI(self):
        ''' When we plot on the GUI, plot on the attached <Canvas>. '''
        if self.initiated_canvas: self.plot(self.Canvas)
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
                plot_quantityVsTime(**plotting_arguments)
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
        of <plot_quantityVsTime>. Create the variables, the tkinter variables
        and the frames to display the options. These widgets are added to the 
        <frame> which is a tab of a <ttk.Notebook> widget on <TabNonlinearTime>.'''
        
        # Attributes for the plot
        self.show_timeFrame = False  
        self.normalizeBySaturation = False
         
        # Create a subframe for each set of plotting options
        self.subframe_yQuantity = ttk.Frame(frame)
        self.subframe_options = ttk.Frame(frame)
         
        # Configure the main frame to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1)      
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=1)      
        tk.Grid.rowconfigure(   frame, 3, weight=1)       
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_options
 
        # Add the subframes to the main frame
        self.subframe_options.grid(row=1, column=0, padx=25, pady=(0,0), stick='NSEW')
         
        # Fill the frames with their options
        self.Options   = self.Frame_Options(  self, self.subframe_options)
        return    
    
    #----------------------  
    class Frame_Options:
         
        def __init__(self, PlotTimeEvolution, frame):
 
            # Make the parents available
            self.PlotTimeEvolution = PlotTimeEvolution 
            self.tab = PlotTimeEvolution.tab
            
            # Variables
            self.var_timeFrame = tk.IntVar(value=self.PlotTimeEvolution.show_timeFrame) 
            self.var_norm = tk.IntVar(value=self.PlotTimeEvolution.normalizeBySaturation)
               
            # Create the widgets
            self.lbl_options = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_vspan   = ttk.Checkbutton(frame, text=" Show time frames") 
            self.chk_norm    = ttk.Checkbutton(frame, text=" Normalize by saturation")
             
            # Add the commands
            self.chk_vspan.config(variable=self.var_timeFrame, command=self.update_vspan)
            self.chk_norm.config(variable=self.var_norm, command=self.update_norm) 
             
            # Configure the frame
            tk.Grid.columnconfigure(frame, 0, weight=1) 
             
            # Add the options to the frame
            self.lbl_options.grid( row=0, column=0, **PAD_TITLE2)
            self.chk_vspan.grid(   row=1, column=0, **PAD_LABEL2)
            self.chk_norm.grid(    row=2, column=0, **PAD_LABEL2) 
            return 
         
        #-----------------------
        def update_vspan(self, *_): 
            self.PlotTimeEvolution.show_timeFrame = self.var_timeFrame.get() 
            self.PlotTimeEvolution.Axis.reset_axis()
            self.PlotTimeEvolution.plotOnTheGUI()
            return 
         
        #-----------------------
        def update_norm(self, *_): 
            self.PlotTimeEvolution.normalizeBySaturation = self.var_norm.get() 
            self.PlotTimeEvolution.Axis.reset_axis()
            self.PlotTimeEvolution.plotOnTheGUI()
            self.tab.DisplayTimeFrames.update_timeFramesAndSaturatedFluxes()
            if True: return 
