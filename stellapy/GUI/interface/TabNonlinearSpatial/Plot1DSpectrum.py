
import tkinter as tk
from tkinter import ttk  
from stellapy.GUI.widgets import Canvas, PAD_TITLE2, PAD_LABEL2
from stellapy.GUI.plot.nonlinearSimulations.plot_oneDimensionalSpectrum import plot_oneDimensionalSpectrum 

#################################################################
#                 CLASS TO PLOT THE LINEAR MAP
#################################################################

class Plot1DSpectrum: 
    ''' Holds and manipulates the Canvas and Axis object for the following plots:
            (1) ... '''

    def __init__(self, tab, notebook_tab): 

        # Set an identifier for this plotting class, used as the figure name
        self.identifier = "TabNonlinearSpatial:Plot1DSpectrum"
        
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

    def initiate_canvas(self, frame1, frame2, notebook_tab):
        ''' The <Canvas> is not generated until its parent frame is visible, this
        saves GUI loading time, since only the visisble Canvasses are created. 
        Create the Canvas of this plot in the given frame, configure and save it. '''
        
        # Initiate the <Canvas>
        if self.initiated_canvas==False: 

            # Create a subclass of <Canvas>
            class CanvasForTwoPlots(Canvas):
                def reset_graph(self): 
                    ''' Reset the axis and then plot again. '''  
                    if "poppepout" not in self.identifier:
                        self.Axis.reset_axis() 
                        self.CanvasII.Axis.reset_axis() 
                        self.plottingClass.plot(self, self.CanvasII)
                    if "poppepout" in self.identifier:
                        self.Axis.reset_axis() 
                        self.plottingClass.plot(self)
                    return
                
            # Create the <Canvas> and remember also the <Axis>
            self.Canvas1 = CanvasForTwoPlots(self.root, frame1, self.identifier+"_kx")
            self.Canvas2 = CanvasForTwoPlots(self.root, frame2, self.identifier+"_ky")
            self.Axis1 = self.Canvas1.Axis; self.Canvas1.CanvasII = self.Canvas2 
            self.Axis2 = self.Canvas2.Axis; self.Canvas2.CanvasII = self.Canvas1
    
            # Attach this plotting class to the canvas           
            self.Canvas1.add_plottingClass(self)
            self.Canvas2.add_plottingClass(self)

            # Resize the <matplotlib figure>
            self.Axis1.update_gridSpeficications(top=0.94, left=0.15, right=0.95, bottom=0.15)
            self.Axis2.update_gridSpeficications(top=0.94, left=0.15, right=0.95, bottom=0.15)
            
            # Set the (x,y) quantities and titles for the canvas
            for i in range(len(self.tab.Simulations.options_yquantityKeys)):  
                x1_quantity = "kx" if (self.tab.Simulations.options_yquantityKeys[i] not in ["gzvs","gvmus"]) else "vpa"
                x2_quantity = "ky" if (self.tab.Simulations.options_yquantityKeys[i] not in ["gzvs","gvmus"]) else "mu" 
                self.Axis1.set_axisQuantities(i, x1_quantity, self.tab.Simulations.options_yquantityKeys[i], name = self.tab.Simulations.options_yquantityNames[i])
                self.Axis2.set_axisQuantities(i, x2_quantity, self.tab.Simulations.options_yquantityKeys[i], name = self.tab.Simulations.options_yquantityNames[i])

            # Create the frames and variables for the plotting arguments
            self.initiate_optionVariablesAndFrames(notebook_tab)
             
            # Tell the <Axis> which quantity to plot, and update the label frame
            self.Axis1.set_plot(self.tab.Simulations.options_yquantityNames.index(self.tab.Simulations.var_quantity.get()))
            self.Axis2.set_plot(self.tab.Simulations.options_yquantityNames.index(self.tab.Simulations.var_quantity.get()))
            self.Axis1.update_labelFrame(self.Canvas1.frame)
            self.Axis2.update_labelFrame(self.Canvas2.frame)
        
            # Remember that the <Canvas> is initiaited
            self.initiated_canvas = True
            
        # Bind <ctrl+s> to replot the graphs
        self.root.bind('<Control-s>', self.Canvas1.reset_graph)
        return
    
    #---------------------- 
    def get_argumentsForPlot(self, Canvas1, Canvas2):
        ''' Gather the arguments for <plot_quantityVsZ>. ''' 
        if Canvas2!=None:
            Canvas1.Axis.update_axisQuantities()
            Canvas2.Axis.update_axisQuantities()
            plotting_arguments = {
                # Specify which simulations to plo
                    "research"      : self.root.Research,\
                    "experiment_id" : self.tab.Simulations.experiment_id,\
                    "simulation_id" : self.tab.Simulations.simulation_id,\
                    "species"       : self.tab.Simulations.species,\
                # Speficy which data to plot
                    "t_specific"    : self.tab.DataSelection.t_specific,\
                    "z_specific"    : self.tab.DataSelection.z_specific,\
                    "x1_quantity"   : Canvas1.Axis.x_quantity,\
                    "x2_quantity"   : Canvas2.Axis.x_quantity,\
                    "y_quantity"    : Canvas1.Axis.y_quantity,\
                # Specify data range
                    "x1_range"      : Canvas1.Axis.x_range,\
                    "y1_range"      : Canvas1.Axis.y_range,\
                    "x2_range"      : Canvas2.Axis.x_range,\
                    "y2_range"      : Canvas2.Axis.y_range,\
                    "units"         : Canvas1.Axis.units,\
                # Labels
                    "x1_label"      : Canvas1.Axis.x_label,\
                    "y1_label"      : Canvas1.Axis.y_label,\
                    "x2_label"      : Canvas2.Axis.x_label,\
                    "y2_label"      : Canvas2.Axis.y_label,\
                    "title1"        : Canvas1.Axis.title,\
                    "title2"        : Canvas2.Axis.title,\
                # Options
                    "loglog"        : self.loglog,\
                    "loglin"        : self.loglin,\
                    "splitZonal"    : self.splitZonal,\
                    "removeZonal"   : self.removeZonal,\
                    "cascade"       : self.cascade,\
                    "scaleToOne"    : self.scaleToOne,\
                # For the GUI the figure object already exists 
                    "show_figure"   : False,\
                    "ax1"           : Canvas1.Axis.ax,\
                    "ax2"           : Canvas2.Axis.ax,\
                    "Progress"      : self.tab.Progress}
        if Canvas2==None:
            Canvas1.Axis.update_axisQuantities()
            plotting_arguments = {
                # Specify which simulations to plot
                    "research"      : self.root.Research,\
                    "experiment_id" : self.tab.Simulations.experiment_id,\
                    "simulation_id" : self.tab.Simulations.simulation_id,\
                # Specify data range
                    "t_specific"    : self.tab.DataSelection.t_specific,\
                    "z_specific"    : self.tab.DataSelection.z_specific,\
                    "x1_quantity"   : Canvas1.Axis.x_quantity,\
                    "x2_quantity"   : Canvas1.Axis.x_quantity,\
                    "y_quantity"    : Canvas1.Axis.y_quantity,\
                # Specify data range
                    "x1_range"      : Canvas1.Axis.x_range,\
                    "y1_range"      : Canvas1.Axis.y_range,\
                    "units"         : Canvas1.Axis.units,\
                # Labels
                    "x1_label"      : Canvas1.Axis.x_label,\
                    "y1_label"      : Canvas1.Axis.y_label,\
                    "title1"        : Canvas1.Axis.title,\
                # Options
                    "loglog"        : self.loglog,\
                    "loglin"        : self.loglin,\
                    "splitZonal"    : self.splitZonal,\
                    "removeZonal"   : self.removeZonal,\
                    "cascade"       : self.cascade,\
                    "scaleToOne"    : self.scaleToOne,\
                # For the GUI the figure object already exists 
                    "show_figure"   : False,\
                    "ax1"           : Canvas1.Axis.ax,\
                    "Progress"      : self.tab.Progress} 
        return plotting_arguments
    
    #---------------------- 
    def plotOnTheGUI(self):
        ''' When we plot on the GUI, plot on the attached <Canvas>. '''
        if self.initiated_canvas: self.plot(self.Canvas1, self.Canvas2)
        return
    
    #----------------------  
    def plot(self, Canvas1, Canvas2=None):
        ''' Plot on the given <Canvas>. '''
         
        # Get the research arguments, plotting arguments and input_files
        research_arguments = self.root.research_arguments
        plotting_arguments = self.get_argumentsForPlot(Canvas1, Canvas2)
        input_files = self.root.Research.input_files

        # Only plot if there are input_files
        if input_files != []:

            # Plot if the plotting arguments do not match those of the <Canvas> 
            if plotting_arguments != Canvas1.Axis.plotting_arguments:
                
                # Plot the 1D spectra along (kx) or (ky)
                Canvas1.Axis.start_plotting(research_arguments, input_files)
                if Canvas2!=None: Canvas2.Axis.start_plotting(research_arguments, input_files)
                plot_oneDimensionalSpectrum(**plotting_arguments)
                Canvas1.Axis.finish_plotting(research_arguments, input_files)
                Canvas1.Axis.save_plottingArguments(self.get_argumentsForPlot(Canvas1, Canvas2)) 
                if Canvas2!=None: Canvas2.Axis.finish_plotting(research_arguments, input_files)
                if Canvas2!=None: Canvas2.Axis.save_plottingArguments(self.get_argumentsForPlot(Canvas1, Canvas2)) 

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
        self.splitZonal = False
        self.removeZonal = False
        self.cascade = False
        self.loglog = False
        self.loglin = False
        self.scaleToOne = False
         
        # Create a subframe for each set of plotting options  
        self.subframe_options = ttk.Frame(frame)
         
        # Configure the main frame to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1,  uniform="yes")
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=2,  uniform="yes")  
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_options    
 
        # Add the subframes to the main frame  
        self.subframe_options.grid(row=1, column=0, padx=25, pady=(0,0), stick='NSEW')
        
        # Fill the frames with their options  
        self.Options = self.Frame_Options(    self, self.subframe_options)
        return 
    
    #----------------------  
    class Frame_Options:
         
        def __init__(self, Plot1DSpectrum, frame):
         
            # Make the parents available
            self.Plot1DSpectrum = Plot1DSpectrum  
                            
            # Variables 
            self.var_splitZonal  = tk.IntVar(value=self.Plot1DSpectrum.splitZonal) 
            self.var_removeZonal = tk.IntVar(value=self.Plot1DSpectrum.removeZonal) 
            self.var_loglog      = tk.IntVar(value=self.Plot1DSpectrum.loglog) 
            self.var_loglin      = tk.IntVar(value=self.Plot1DSpectrum.loglin) 
            self.var_cascade     = tk.IntVar(value=self.Plot1DSpectrum.cascade)
            self.var_scaleToOne  = tk.IntVar(value=self.Plot1DSpectrum.scaleToOne)
        
            # Create the widgets
            self.lbl_options     = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_splitZonal  = ttk.Checkbutton(frame, text=" Split in zonal contributions")  
            self.chk_removeZonal = ttk.Checkbutton(frame, text=" Remove the zonal mode")  
            self.chk_scaleToOne  = ttk.Checkbutton(frame, text=" Rescale to one")   
            self.chk_loglog      = ttk.Checkbutton(frame, text=" Use (log,log) scales")  
            self.chk_loglin      = ttk.Checkbutton(frame, text=" Use (log,lin) scales")  
            self.chk_cascade     = ttk.Checkbutton(frame, text=" Add the 7/3 power law (cascade)")  
     
            # Add the commands
            self.chk_splitZonal.config( variable=self.var_splitZonal,  command=self.update_splitZonal)  
            self.chk_removeZonal.config(variable=self.var_removeZonal, command=self.update_removeZonal) 
            self.chk_scaleToOne.config( variable=self.var_scaleToOne,  command=self.update_scaleToOne) 
            self.chk_loglog.config(     variable=self.var_loglog,      command=self.update_loglog) 
            self.chk_loglin.config(     variable=self.var_loglin,      command=self.update_loglin) 
            self.chk_cascade.config(    variable=self.var_cascade,     command=self.update_cascade) 
             
            # Configure the frame
            tk.Grid.rowconfigure(frame, 0, weight=0) 
            tk.Grid.rowconfigure(frame, 1, weight=0) 
            tk.Grid.rowconfigure(frame, 2, weight=0) 
            tk.Grid.columnconfigure(frame, 0, weight=0) 
            tk.Grid.columnconfigure(frame, 1, weight=0) 
            tk.Grid.columnconfigure(frame, 2, weight=0)  
            tk.Grid.columnconfigure(frame, 3, weight=1) 
             
            # Add the options to the frame
            self.lbl_options.grid(     row=0, column=0, **PAD_TITLE2)
            self.chk_splitZonal.grid(  row=1, column=0, **PAD_LABEL2)
            self.chk_removeZonal.grid( row=2, column=0, **PAD_LABEL2)
            self.chk_scaleToOne.grid(  row=3, column=0, **PAD_LABEL2)
            self.chk_loglog.grid(      row=1, column=1, **PAD_LABEL2)
            self.chk_loglin.grid(      row=2, column=1, **PAD_LABEL2)
            self.chk_cascade.grid(     row=3, column=1, **PAD_LABEL2)
            return 
                            
        #-------------------
        def update_splitZonal(self, *_): 
            self.Plot1DSpectrum.splitZonal = self.var_splitZonal.get() 
            self.Plot1DSpectrum.plotOnTheGUI()
            return 
               
        #-------------------
        def update_removeZonal(self, *_): 
            self.Plot1DSpectrum.removeZonal = self.var_removeZonal.get() 
            self.Plot1DSpectrum.plotOnTheGUI()
            return 

        #-------------------
        def update_scaleToOne(self, *_): 
            self.Plot1DSpectrum.scaleToOne = self.var_scaleToOne.get() 
            self.Plot1DSpectrum.plotOnTheGUI()
            return 
        
        #-------------------
        def update_loglog(self, *_): 
            self.Plot1DSpectrum.loglog = self.var_loglog.get() 
            self.Plot1DSpectrum.Axis1.reset_axis()
            self.Plot1DSpectrum.Axis2.reset_axis()
            self.Plot1DSpectrum.plotOnTheGUI()
            return 
        
        #-------------------
        def update_loglin(self, *_): 
            self.Plot1DSpectrum.loglin = self.var_loglin.get() 
            self.Plot1DSpectrum.Axis1.reset_axis()
            self.Plot1DSpectrum.Axis2.reset_axis()
            self.Plot1DSpectrum.plotOnTheGUI()
            return 
        
        #-------------------
        def update_cascade(self, *_): 
            self.Plot1DSpectrum.cascade = self.var_cascade.get() 
            self.Plot1DSpectrum.plotOnTheGUI()
            return 
        
 
  