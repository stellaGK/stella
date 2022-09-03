
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
        self.identifier = "TabNonlinearTime:PlotTimeEvolution"
        
        # Make the parent objects available
        self.tab = tab
        self.root = tab.root  
        
        # Remember whether we have a Canvas yet
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
            self.Axis.update_gridSpeficications(top=0.94, left=0.12, right=0.95, bottom=0.12)
            
            # Set the (x,y) quantities and titles for the canvas
            self.Axis.set_axisQuantities(1,  "t", "qflux",       name = "Heat flux")
            self.Axis.set_axisQuantities(2,  "t", "vflux",       name = "Momentum flux")
            self.Axis.set_axisQuantities(3,  "t", "pflux",       name = "Particle flux") 
            self.Axis.set_axisQuantities(4,  "t", "g",           name = "Distribution") 
            self.Axis.set_axisQuantities(5,  "t", "phi2",        name = "Potential squared") 
            self.Axis.set_axisQuantities(6,  "t", "phi_real",    name = "Real part of the potential") 
            self.Axis.set_axisQuantities(7,  "t", "phi_imag",    name = "Imaginary part of the potential") 
            self.Axis.set_axisQuantities(8,  "t", "phi2_nozonal",name = "Potential squared without zonal modes") 
            self.Axis.set_axisQuantities(9,  "t", "phiRealImag", name = "Real and imaginary part of phi") 
            self.Axis.set_axisQuantities(10, "t", "phi2_split",  name = "Phi² split in zonal modes") 
            self.Axis.set_axisQuantities(11, "t", "phi2AndQFlux",name = "Phi² without zonal modes and heat flux") 

            # Create the frames and variables for the plotting arguments
            self.initiate_optionVariablesAndFrames(notebook_tab)
             
            # Tell the <Axis> which quantity to plot, and update the label frame
            self.Axis.set_plot(self.YQuantity.var_yQuantity.get())
            self.Axis.update_labelFrame(self.Canvas.frame)
        
            # Remember that the <Canvas> is initiaited
            self.initiated_canvas = True
        return

    #---------------------- 
    def get_argumentsForPlot(self, Canvas):
        ''' Gather the arguments for <plot_quantityVsZ>. '''
        self.Axis.update_axisQuantities()
        plotting_arguments = {
            # Simulations
                "research"      : self.root.Research,\
                "experiment_id" : self.tab.Simulations.experiment_id,\
                "simulation_id" : self.tab.Simulations.simulation_id,\
                "species"       : self.tab.Simulations.species,\
            # Data
                "x_quantity"    : Canvas.Axis.x_quantity,\
                "y_quantity"    : Canvas.Axis.y_quantity,\
                "x_range"       : Canvas.Axis.x_range,\
                "y_range"       : Canvas.Axis.y_range,\
                "units"         : Canvas.Axis.units,\
            # Specify data range 
                "x_label"       : Canvas.Axis.x_label,\
                "y_label"       : Canvas.Axis.y_label,\
                "title"         : Canvas.Axis.title,\
            # Toggles 
                "normalize"     : self.normalizeBySaturation,\
                "show_timeFrame": self.show_timeFrame,\
                #"includeFluxNorm" : self.includeFluxNorm,\
            # Figure
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
        
        # Adjustable plotting arguments 
        self.y_quantity = "qflux" 
        self.show_timeFrame = False 
        self.normalizeBySaturation = False
        self.units = "normalized"
        self.includeFluxNorm = False
        
        # Create a subframe for each set of plotting options
        self.subframe_yQuantity = ttk.Frame(frame)
        self.subframe_options = ttk.Frame(frame)
         
        # Configure the main frame to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1)      
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=1)      
        tk.Grid.rowconfigure(   frame, 3, weight=1)       
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_yQuantity
        tk.Grid.columnconfigure(frame, 1, weight=0) # subframe_options
 
        # Add the subframes to the main frame
        self.subframe_yQuantity.grid(row=1, column=0, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_options.grid(  row=1, column=1, padx=25, pady=(0,0), stick='NSEW')
         
        # Fill the frames with their options
        self.YQuantity = self.Frame_YQuantity(self, self.subframe_yQuantity)
        self.Options   = self.Frame_Options(  self, self.subframe_options)
        return     
    
    #----------------------  
    class Frame_YQuantity:
         
        def __init__(self, PlotTimeEvolution, frame):
             
            # Make the tabs available
            self.PlotTimeEvolution = PlotTimeEvolution 
            self.tab = PlotTimeEvolution.tab
            
            # Variables
            self.var_yQuantity = tk.IntVar(value=1)
            
            # Get the possible y-quantities that were set in <initiate_canvas>
            plot_ids, axis_quantities = self.PlotTimeEvolution.Axis.get_axisQuantities()
            self.radiobuttons = [None]*len(plot_ids)
        
            # Configure the frame <self.subframe_yQuantity>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            tk.Grid.columnconfigure(frame, 1, weight=1) 
            tk.Grid.columnconfigure(frame, 2, weight=1) 
             
            # Create the widgets
            self.lbl_yQuantity = ttk.Label(frame, text="Quantity on the y-axis", style='prefTitle.TLabel')
            for i, plot_id in enumerate(plot_ids): self.radiobuttons[i] = ttk.Radiobutton(frame, text='  '+axis_quantities[plot_id]['name']) 
            for i, plot_id in enumerate(plot_ids): self.radiobuttons[i].config(value=plot_id, variable=self.var_yQuantity, command=self.change_yQuantity)
             
            # Add the radio buttons to the frame
            self.lbl_yQuantity.grid(   row=1, column=0, **PAD_TITLE2)
            self.radiobuttons[0].grid( row=2, column=0, **PAD_LABEL2)
            self.radiobuttons[1].grid( row=3, column=0, **PAD_LABEL2)
            self.radiobuttons[2].grid( row=4, column=0, **PAD_LABEL2)
            self.radiobuttons[3].grid( row=5, column=0, **PAD_LABEL2)
            self.radiobuttons[4].grid( row=2, column=1, **PAD_LABEL2)
            self.radiobuttons[5].grid( row=3, column=1, **PAD_LABEL2)
            self.radiobuttons[6].grid( row=4, column=1, **PAD_LABEL2)
            self.radiobuttons[7].grid( row=5, column=1, **PAD_LABEL2)
            self.radiobuttons[8].grid( row=2, column=2, **PAD_LABEL2)
            self.radiobuttons[9].grid( row=3, column=2, **PAD_LABEL2)
            self.radiobuttons[10].grid(row=4, column=2, **PAD_LABEL2)
            return 

        #---------------------
        def change_yQuantity(self, *_):
            ''' Update the (x,y) quantities, reset the axis since the data has
            changed, and change the title of the label frame of the canvas. ''' 
            self.PlotTimeEvolution.Axis.set_plot(self.var_yQuantity.get())
            self.PlotTimeEvolution.Axis.reset_axis()
            self.PlotTimeEvolution.Axis.update_labelFrame(self.PlotTimeEvolution.Canvas.frame)
            self.PlotTimeEvolution.plotOnTheGUI()
            return
        
    #---------------------- 
    class Frame_Options:
         
        def __init__(self, PlotTimeEvolution, frame):
 
            # Make the tabs available
            self.PlotTimeEvolution = PlotTimeEvolution  
            self.tab = PlotTimeEvolution.tab
            
            # Variables
            self.var_timeFrame = tk.IntVar(value=self.PlotTimeEvolution.show_timeFrame) 
            self.var_norm = tk.IntVar(value=self.PlotTimeEvolution.normalizeBySaturation)
            self.var_SI = tk.IntVar(value=0)
            self.var_includefluxnorm = tk.IntVar(value=0)
               
            # Create the widgets
            self.lbl_options   = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_timeFrame = ttk.Checkbutton(frame, text=" Show time frames") 
            self.chk_norm      = ttk.Checkbutton(frame, text=" Normalize by saturation")
            self.chk_SI        = ttk.Checkbutton(frame, text=" Plot in SI Units")
            self.chk_includefluxnorm = ttk.Checkbutton(frame, text=" Include 1/|<nabla r>|")
             
            # Add the commands
            self.chk_timeFrame.config(variable=self.var_timeFrame, command=self.update_timeFrame)
            self.chk_norm.config(     variable=self.var_norm,      command=self.update_norm) 
            self.chk_SI.config(       variable=self.var_SI,        command=self.update_SI) 
            self.chk_includefluxnorm.config(variable=self.var_includefluxnorm, command=self.update_includefluxnorm)
             
            # Configure the frame
            tk.Grid.columnconfigure(frame, 0, weight=1) 
             
            # Add the options to the frame
            self.lbl_options.grid(   row=0, column=0, **PAD_TITLE2)
            self.chk_timeFrame.grid( row=1, column=0, **PAD_LABEL2)
            self.chk_SI.grid(        row=2, column=0, **PAD_LABEL2) 
            self.chk_norm.grid(      row=3, column=0, **PAD_LABEL2) 
            self.chk_includefluxnorm.grid(row=4, column=0, **PAD_LABEL2) 
            return 
         
        #-----------------------
        def update_timeFrame(self, *_): 
            if self.var_timeFrame.get()==0: self.PlotTimeEvolution.show_timeFrame = False
            if self.var_timeFrame.get()==1: self.PlotTimeEvolution.show_timeFrame = True
            self.PlotTimeEvolution.plotOnTheGUI()
            return

        #-----------------------
        def update_SI(self, *_): 
            if self.var_SI.get()==0: self.PlotTimeEvolution.units = "normalized"
            if self.var_SI.get()==1: self.PlotTimeEvolution.units = "SI"
            self.PlotTimeEvolution.Axis.reset_axis()
            self.PlotTimeEvolution.plotOnTheGUI()
            return
            
        #-----------------------
        def update_norm(self, *_): 
            if self.var_norm.get()==0: self.PlotTimeEvolution.normalizeBySaturation = False
            if self.var_norm.get()==1: self.PlotTimeEvolution.normalizeBySaturation = True
            self.PlotTimeEvolution.Axis.reset_axis()
            self.PlotTimeEvolution.plotOnTheGUI()
            return        
        
        #----------------------  
        def update_includefluxnorm(self, *_): 
            if self.var_includefluxnorm.get()==0: self.PlotTimeEvolution.includeFluxNorm = False
            if self.var_includefluxnorm.get()==1: self.PlotTimeEvolution.includeFluxNorm = True 
            if self.var_includefluxnorm.get()==0: self.tab.PlotSaturatedFluxVersusParameter.includeFluxNorm = False
            if self.var_includefluxnorm.get()==1: self.tab.PlotSaturatedFluxVersusParameter.includeFluxNorm = True 
            self.PlotTimeEvolution.plotOnTheGUI()
               
