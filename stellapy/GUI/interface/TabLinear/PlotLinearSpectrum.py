
import tkinter as tk
from tkinter import ttk   
from stellapy.GUI.widgets import Canvas, PAD_TITLE2, PAD_LABEL2
from stellapy.GUI.plot.linearSimulations.plot_quantityVsModes import plot_quantityVsModes
 
################################################################################
#            CLASS TO PLOT THE GROWTH RATE OR FREQUENCY SPECTRA
################################################################################
                    
class PlotLinearSpectrum: 
    ''' Holds and manipulates the Canvas and Axis object for the following plots:
            (1) Omega versus k
            (2) Gamma versus k
            (3) Gamma/k**2 versus k '''
    
    def __init__(self, tab, notebook_tab):

        # Set an identifier for this plotting class, used as the figure name
        self.identifier = "TabLinear:PlotLinearSpectrums"
        
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

    def initiate_canvas(self, frame):
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
            self.Axis.update_gridSpeficications(top=0.94, left=0.15, right=0.95, bottom=0.12)
            
            # Set the (x,y) quantities and titles for the canvas
            self.Axis.set_axisQuantities(1, "ky", "omega_avg")
            self.Axis.set_axisQuantities(2, "ky", "gamma_avg")
            self.Axis.set_axisQuantities(3, "ky", "gamma_avg/ky**2") 
            
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
            # Data
                "x_quantity"    : Canvas.Axis.x_quantity,\
                "y_quantity"    : Canvas.Axis.y_quantity,\
                "x_range"       : Canvas.Axis.x_range,\
                "y_range"       : Canvas.Axis.y_range,\
                "units"         : Canvas.Axis.units,\
            # Modes
                "modes_id"      : "unstable",\
                "kx_range"      : Canvas.Axis.kx_range,\
                "ky_range"      : Canvas.Axis.ky_range,\
            # Labels
                "x_label"       : Canvas.Axis.x_label,\
                "y_label"       : Canvas.Axis.y_label,\
                "title"         : Canvas.Axis.title,\
            # Toggles
                "show_error"    : self.show_error,\
                "maxima"        : self.maxima,\
                "biggestmaxima" : self.biggestmaxima,\
                "interpolate"   : self.interpolate,\
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
                plot_quantityVsModes(**plotting_arguments)
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
        of <plot_quantityVsModes>. Create the variables, the tkinter 
        variables and the frames to display the options. These widgets are added
        to the <frame> which is a tab of a <ttk.Notebook> widget on <tabLinear>.'''
    
        # Adjustable plotting arguments 
        self.show_error = True
        self.interpolate = False
        self.maxima = False
        self.biggestmaxima = True 
        
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
        
        def __init__(self, PlotLinearSpectrum, frame):
            
            # Make the tabs available
            self.PlotLinearSpectrum = PlotLinearSpectrum 
            
            # Variables 
            self.var_yQuantity = tk.IntVar(value=2)
        
            # Configure the frame <self.subframe_yQuantity>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Create the widgets
            self.lbl_spectrum  = ttk.Label(frame, text="Quantity on the y-axis", style='prefTitle.TLabel')
            self.rbn_freq      = ttk.Radiobutton(frame, text='  Frequency')
            self.rbn_growth    = ttk.Radiobutton(frame, text='  Growth rate')
            self.rbn_growthky2 = ttk.Radiobutton(frame, text='  Growth rate on ky**2')
    
            # Add the values and commands to the radiobuttons
            self.rbn_freq.config(      value=1, variable=self.var_yQuantity, command=self.change_yQuantity)
            self.rbn_growth.config(    value=2, variable=self.var_yQuantity, command=self.change_yQuantity)
            self.rbn_growthky2.config( value=3, variable=self.var_yQuantity, command=self.change_yQuantity)
            
            # Add the options to the frame
            self.lbl_spectrum.grid( row=1,  column=0, **PAD_TITLE2)
            self.rbn_freq.grid(     row=2,  column=0, **PAD_LABEL2)
            self.rbn_growth.grid(   row=3,  column=0, **PAD_LABEL2)
            self.rbn_growthky2.grid(row=4,  column=0, **PAD_LABEL2)
            return 
        
        #---------------------
        def change_yQuantity(self, *_):
            ''' Update the (x,y) quantities, reset the axis since the data has
            changed, and change the title of the label frame of the canvas. ''' 
            self.PlotLinearSpectrum.Axis.set_plot(self.var_yQuantity.get())
            self.PlotLinearSpectrum.Axis.reset_axis()
            self.PlotLinearSpectrum.Axis.update_labelFrame(self.PlotLinearSpectrum.Canvas.frame)
            self.PlotLinearSpectrum.plotOnTheGUI()
            return

    #----------------------  
    class Frame_Options:
        
        def __init__(self, PlotLinearSpectrum, frame):

            # Make the tabs available
            self.PlotLinearSpectrum = PlotLinearSpectrum  
            
            # Variables 
            self.var_maxima = tk.IntVar(value=self.PlotLinearSpectrum.maxima)
            self.var_maxima2 = tk.IntVar(value=self.PlotLinearSpectrum.biggestmaxima)
            self.var_interpolate = tk.IntVar(value=self.PlotLinearSpectrum.interpolate)
            self.var_plotErrorBars = tk.IntVar(value=self.PlotLinearSpectrum.show_error) 
              
            # Create the widgets
            self.lbl_options = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_error   = ttk.Checkbutton(frame, text=" Include error bars")
            self.chk_maxima  = ttk.Checkbutton(frame, text=" Find all the relative maxima")
            self.chk_maxima2 = ttk.Checkbutton(frame, text=" Find the absolute maximum")
            self.chk_interp  = ttk.Checkbutton(frame, text=" Interpolate and find maxima") 
            
            # Add the commands
            self.chk_error.config(  variable=self.var_plotErrorBars, command=self.update_error)
            self.chk_maxima.config( variable=self.var_maxima,        command=self.update_maxima)
            self.chk_maxima2.config(variable=self.var_maxima2,       command=self.update_maxima)
            self.chk_interp.config( variable=self.var_interpolate,   command=self.update_interp) 
            
            # Configure the frame
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Add the options to the frame
            self.lbl_options.grid( row=0, column=0, **PAD_TITLE2)
            self.chk_error.grid(   row=1, column=0, **PAD_LABEL2)
            self.chk_maxima.grid(  row=2, column=0, **PAD_LABEL2)
            self.chk_maxima2.grid( row=3, column=0, **PAD_LABEL2)
            self.chk_interp.grid(  row=4, column=0, **PAD_LABEL2) 
            return 
        
        #-----------------------
        def update_error(self): 
            if self.var_plotErrorBars.get()==0: self.PlotLinearSpectrum.show_error = False
            if self.var_plotErrorBars.get()==1: self.PlotLinearSpectrum.show_error = True
            self.PlotLinearSpectrum.plotOnTheGUI()
            return 
        
        #-----------------------
        def update_maxima(self): 
            if self.var_maxima.get()==0:  self.PlotLinearSpectrum.maxima = False
            if self.var_maxima.get()==1:  self.PlotLinearSpectrum.maxima = True
            if self.var_maxima2.get()==0: self.PlotLinearSpectrum.biggestmaxima = False
            if self.var_maxima2.get()==1: self.PlotLinearSpectrum.biggestmaxima = True
            self.PlotLinearSpectrum.plotOnTheGUI()
            return 
                   
        #-----------------------
        def update_interp(self): 
            if self.var_interpolate.get()==0: self.PlotLinearSpectrum.interpolate = False
            if self.var_interpolate.get()==1: self.PlotLinearSpectrum.interpolate = True
            self.PlotLinearSpectrum.plotOnTheGUI()
            if True: return 
                    
 