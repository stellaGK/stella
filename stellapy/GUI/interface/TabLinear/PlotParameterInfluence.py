 
import tkinter as tk
from tkinter import ttk  
from stellapy.GUI.widgets import Canvas, PAD_TITLE2, PAD_LABEL2 
from stellapy.GUI.plot.linearSimulations.plot_quantityVsParameter import plot_quantityVsParameter
from stellapy.plot.utils.labels.standardParameters import standardParameters
from stellapy.plot.utils.labels.standardParameters import standardParametersInOrder 
################################################################################
#            CLASS TO PLOT THE MOST UNSTABLE MODE VERSUS PARAMETER
################################################################################
   
class PlotParameterInfluence: 
    ''' Holds and manipulates the Canvas and Axis object for the following plots:
            (1) Omega versus k
            (2) Gamma versus k
            (3) Gamma/k**2 versus k '''
    
    def __init__(self, tab, notebook_tab):

        # Set an identifier for this plotting class, used as the figure name
        self.identifier = "TabLinear:PlotParameterInfluence"
        
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
            self.Axis.update_gridSpeficications(top=0.9, left=0.15, right=0.95, bottom=0.12)
            
            # Set the (x,y) quantities and titles for the canvas
            self.Axis.set_axisQuantities(1, "parameter", "omega_avg")
            self.Axis.set_axisQuantities(2, "parameter", "gamma_avg")
            self.Axis.set_axisQuantities(3, "parameter", "gamma_avg/ky**2")  
            
            # Tell the <Axis> which quantity to plot
            self.Axis.set_plot(self.YQuantity.var_yQuantity.get()) 
            self.Axis.update_labelFrame(self.Canvas.frame)
            
            # Remember that the <Canvas> is initiaited
            self.initiated_canvas = True
        return     

    #---------------------- 
    def get_argumentsForPlot(self, Canvas):
        ''' Gather the arguments for <plot_quantityVsParameter>. '''
        self.Axis.update_axisQuantities()
        plotting_arguments = {
            # Simulations
                "research"      : self.root.Research,\
                "experiment_id" : self.tab.Simulations.experiment_id,\
                "simulation_id" : self.tab.Simulations.simulation_id,\
            # Parameter
                "knob"          : self.knob,\
                "key"           : self.key,\
            # Data
                "x_quantity"    : Canvas.Axis.x_quantity,\
                "y_quantity"    : Canvas.Axis.y_quantity,\
                "units"         : Canvas.Axis.units,\
                "x_range"       : Canvas.Axis.x_range,\
                "y_range"       : Canvas.Axis.y_range,\
            # Modes
                "kx_range"      : Canvas.Axis.kx_range,\
                "ky_range"      : Canvas.Axis.ky_range,\
            # Labels
                "x_label"       : Canvas.Axis.x_label,\
                "y_label"       : Canvas.Axis.y_label,\
                "title"         : Canvas.Axis.title,\
            # Toggles
                "add_dottedLine": self.add_dottedLine,\
                "show_error"    : self.show_error,\
                "splitInKineticAdiabatic" : self.splitInKineticAdiabatic,\
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
                plot_quantityVsParameter(**plotting_arguments)
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
        self.k_value = "max"
        self.show_error = False
        self.add_dottedLine = False
        self.knob = "-"
        self.key  = "-" 
        self.splitInKineticAdiabatic = False
        
        # Create a subframe for each set of plotting options
        self.subframe_yQuantity  = ttk.Frame(frame)
        self.subframe_parameter  = ttk.Frame(frame)
        self.subframe_extraction = ttk.Frame(frame)
        self.subframe_options    = ttk.Frame(frame)
        
        # Configure the rows/columns of the <frame> to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1)      
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=1)      
        tk.Grid.rowconfigure(   frame, 3, weight=1)      
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_yQuantity
        tk.Grid.columnconfigure(frame, 1, weight=0) # subframe_parameter
        tk.Grid.columnconfigure(frame, 2, weight=0) # subframe_extraction
        tk.Grid.columnconfigure(frame, 3, weight=0) # subframe_options

        # Add the subframes to the main frame
        self.subframe_yQuantity.grid( row=1, column=0, padx=25, pady=(0,0), stick='NSEW', rowspan=2)
        self.subframe_parameter.grid( row=1, column=1, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_extraction.grid(row=1, column=2, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_options.grid(   row=1, column=3, padx=25, pady=(0,0), stick='NSEW', rowspan=2)
        
        # Fill the frames with widgets
        self.YQuantity  = self.Frame_YQuantity( self, self.subframe_yQuantity)
        self.Parameter  = self.Frame_Parameter( self, self.subframe_parameter)
        self.Extraction = self.Frame_Extraction(self, self.subframe_extraction)
        self.Options    = self.Frame_Options(   self, self.subframe_options)
        return  

    #----------------------  
    class Frame_YQuantity:
        
        def __init__(self, PlotParameterInfluence, frame):

            # Make the parents available  
            self.PlotParameterInfluence = PlotParameterInfluence 
            
            # Variables
            self.var_yQuantity = tk.IntVar(value=2)
        
            # Configure the frame <self.subframe_yQuantity>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Create the widgets
            self.lbl_spectrum  = ttk.Label(frame, text="Quantity on the Y-axis", style='prefTitle.TLabel')
            self.rbn_freq      = ttk.Radiobutton(frame, text='  Frequency')
            self.rbn_growth    = ttk.Radiobutton(frame, text='  Growth rate')
            self.rbn_growthky2 = ttk.Radiobutton(frame, text='  Growth rate on ky**2')
            self.rbn_growthInt = ttk.Radiobutton(frame, text='  Integrated growth rate')
    
            # Add the values and commands to the radiobuttons
            self.rbn_freq.config(      value=1, variable=self.var_yQuantity, command=self.change_yQuantity)
            self.rbn_growth.config(    value=2, variable=self.var_yQuantity, command=self.change_yQuantity)
            self.rbn_growthky2.config( value=3, variable=self.var_yQuantity, command=self.change_yQuantity)
            self.rbn_growthInt.config( value=4, variable=self.var_yQuantity, command=self.change_yQuantity)
            
            # Add the options to the frame
            self.lbl_spectrum.grid( row=1,  column=0, **PAD_TITLE2)
            self.rbn_freq.grid(     row=2,  column=0, **PAD_LABEL2)
            self.rbn_growth.grid(   row=3,  column=0, **PAD_LABEL2)
            self.rbn_growthky2.grid(row=4,  column=0, **PAD_LABEL2)
            self.rbn_growthInt.grid(row=5,  column=0, **PAD_LABEL2)
            return 
        
        #---------------------
        def change_yQuantity(self, *_):
            ''' Update the (x,y) quantities, reset the axis since the data has
            changed, and change the header of the label frame of the canvas. ''' 
            self.PlotParameterInfluence.Axis.set_plot(self.var_yQuantity.get())
            self.PlotParameterInfluence.Axis.reset_axis()
            self.PlotParameterInfluence.Axis.update_labelFrame(self.PlotParameterInfluence.Canvas.frame)
            self.PlotParameterInfluence.plotOnTheGUI()
            return
    
    #---------------------- 
    class Frame_Parameter:
        
        def __init__(self, PlotParameterInfluence, frame):

            # Make the parents available 
            self.root = PlotParameterInfluence.tab.root 
            self.PlotParameterInfluence = PlotParameterInfluence 
            
            # Variables 
            self.var_parameter = tk.StringVar(value=standardParametersInOrder[0])
        
            # Configure the frame <self.subframe_parameter>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Choose the parameter for the x-axis  
            self.lbl_prm = ttk.Label(frame, text="Parameters", style='prefTitle.TLabel'); width=10
            self.mnu_par = ttk.OptionMenu(frame, self.var_parameter, standardParametersInOrder[0], *standardParametersInOrder, style='option.TMenubutton')
            self.mnu_par["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
            self.var_parameter.trace_id = self.var_parameter.trace('w', self.update_parameter); self.mnu_par.config(width=width)

            # Add the widgets to the frame
            self.lbl_prm.grid(row=0, column=0, **PAD_TITLE2)
            self.mnu_par.grid(row=1, column=0, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)

        #---------------------- 
        def update_parameter(self, calledFromGui=True, *_):
            ''' The parameter plotted on the x-axis can be chose from a dropdown
            menu. When a parameter is selected, read the stella <knob> and <key> 
            from the <standardParameters> dictionary, update the x-quantity of
            the 3 possibly plots and reset the axis of the maplotlib figure. '''
            
            # Read which parameter is selected
            parameter = self.var_parameter.get()
            
            # Select the stella <knob> and <key> for the selected <parameter>
            self.PlotParameterInfluence.knob = standardParameters[parameter]["knob"]
            self.PlotParameterInfluence.key  = standardParameters[parameter]["key"]

            # Set the (x,y) quantities for the axis and the canvas
            self.PlotParameterInfluence.Axis.set_axisQuantities(1, parameter, "omega_avg")
            self.PlotParameterInfluence.Axis.set_axisQuantities(2, parameter, "gamma_avg")
            self.PlotParameterInfluence.Axis.set_axisQuantities(3, parameter, "gamma_avg/ky**2")  
            self.PlotParameterInfluence.Axis.update_axisQuantities()
            
            # Update the label in the label frame
            self.PlotParameterInfluence.Axis.update_labelFrame(self.PlotParameterInfluence.Canvas.frame)
            
            # Reset the axis when we change the parameter plotted on the x-axis
            if calledFromGui: self.PlotParameterInfluence.Axis.reset_axis()  
            if calledFromGui: self.PlotParameterInfluence.plotOnTheGUI()
            return

        #----------------------------
        def update_parameterWithoutPlotting(self, *_):
            self.update_parameter(calledFromGui=False)
            return
    
    #----------------------     
    class Frame_Extraction:
        
        def __init__(self, PlotParameterInfluence, frame):

            # Make the parents available 
            self.root = PlotParameterInfluence.tab.root 
            self.PlotParameterInfluence = PlotParameterInfluence 
            
            # Variables
            self.options_k = ["-"]
            self.var_extraction     = tk.IntVar(value=1)
            self.var_specificKvalue = tk.StringVar(value=self.options_k[0])
            
            # Configure the frame <self.subframe_extraction>
            tk.Grid.columnconfigure(frame, 0, weight=0) 
            tk.Grid.columnconfigure(frame, 1, weight=0) 
            tk.Grid.columnconfigure(frame, 2, weight=1) 
            
            # Choose which growth rate we extract: the maximum one or at a specific wavenumber
            self.lbl_extraction = ttk.Label(frame, text="Extraction", style='prefTitle.TLabel')
            self.rbn_gammamax = ttk.Radiobutton(frame, text='  Plot the most unstable mode')
            self.rbn_gammak   = ttk.Radiobutton(frame, text='  Plot the mode at ky =')
      
            # Add the commands to the radio buttons
            self.rbn_gammamax.config( value=1, variable=self.var_extraction, command=self.change_extraction)
            self.rbn_gammak.config(   value=2, variable=self.var_extraction, command=self.change_extraction)
      
            # Drop down menu to choose the k-value
            self.mnu_kvalue = ttk.OptionMenu(frame, self.var_specificKvalue, self.options_k[0], *self.options_k, style='option.TMenubutton'); width=4
            self.mnu_kvalue["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
            self.mnu_kvalue.config(width=width)
            self.var_specificKvalue.trace_id = self.var_specificKvalue.trace('w', self.update_selectedMode) # link function to a change of the dropdown options
        
            # Add the widgets to the frame
            self.lbl_extraction.grid(row=0,  column=0, **PAD_TITLE2, columnspan=3)
            self.rbn_gammamax.grid(  row=1,  column=0, **PAD_LABEL2, columnspan=3)
            self.rbn_gammak.grid(    row=2,  column=0, **PAD_LABEL2)
            self.mnu_kvalue.grid(    row=2,  column=1, **PAD_LABEL2) 
            return 
        
        #---------------------- 
        def update_selectedMode(self, calledFromGui=True, *_):
            if self.var_extraction.get()==1: self.PlotParameterInfluence.k_value = "max"
            if self.var_extraction.get()==2: self.PlotParameterInfluence.k_value = float(self.var_specificKvalue.get())
            if calledFromGui: self.PlotParameterInfluence.plotOnTheGUI() 

        #----------------------------
        def update_selectedModeWithoutPlotting(self, *_):
            self.update_selectedMode(calledFromGui=False)
        
        #---------------------- 
        def change_extraction(self):
            if self.var_extraction.get()==1: self.PlotParameterInfluence.k_value = "max"
            if self.var_extraction.get()==2: self.PlotParameterInfluence.k_value = float(self.var_specificKvalue.get())
            self.PlotParameterInfluence.plotOnTheGUI()
            return

    #---------------------- 
    class Frame_Options:
        
        def __init__(self, PlotParameterInfluence, frame):

            # Make the parents available
            self.PlotParameterInfluence = PlotParameterInfluence  
            
            # Variables 
            self.var_plotErrorBars  = tk.IntVar(value=0)
            self.var_plotDottedLine = tk.IntVar(value=0)
            self.var_splitInKineticAdiabatic = tk.IntVar(value=0)
              
            # Create the widgets
            self.lbl_options = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_error = ttk.Checkbutton(frame, text=" Include error bars")
            self.chk_addLine = ttk.Checkbutton(frame, text=" Add dotted line")
            self.chk_split = ttk.Checkbutton(frame, text=" Split in kinetic/adiabatic")
            
            # Add the commands
            self.chk_error.config(variable=self.var_plotErrorBars, command=self.update_error)
            self.chk_addLine.config(variable=self.var_plotDottedLine, command=self.update_dottedLine)
            self.chk_split.config(variable=self.var_splitInKineticAdiabatic, command=self.update_split)
            
            # Configure the frame
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Add the options to the frame
            self.lbl_options.grid(row=0, column=0, **PAD_TITLE2)
            self.chk_error.grid(  row=1, column=0, **PAD_LABEL2)
            self.chk_addLine.grid(row=2, column=0, **PAD_LABEL2)
            self.chk_split.grid(  row=3, column=0, **PAD_LABEL2)
        
        #---------------------- 
        def update_error(self, *_): 
            if self.var_plotErrorBars.get()==0: self.PlotParameterInfluence.show_error = False
            if self.var_plotErrorBars.get()==1: self.PlotParameterInfluence.show_error = True
            self.PlotParameterInfluence.plotOnTheGUI()
        
        #---------------------- 
        def update_dottedLine(self, *_): 
            if self.var_plotDottedLine.get()==0: self.PlotParameterInfluence.add_dottedLine = False
            if self.var_plotDottedLine.get()==1: self.PlotParameterInfluence.add_dottedLine = True
            self.PlotParameterInfluence.plotOnTheGUI()
            
        #---------------------- 
        def update_split(self, *_): 
            if self.var_splitInKineticAdiabatic.get()==0: self.PlotParameterInfluence.splitInKineticAdiabatic = False
            if self.var_splitInKineticAdiabatic.get()==1: self.PlotParameterInfluence.splitInKineticAdiabatic = True
            self.PlotParameterInfluence.plotOnTheGUI()
    
   