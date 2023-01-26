
# Load modules
import tkinter as tk
from tkinter import ttk    
from stellapy.plot.utils.labels.standardParameters import standardParameters, standardParametersInOrder 
from stellapy.GUI.widgets import Canvas, PAD_TITLE2, PAD_LABEL2 
from stellapy.GUI.plot.nonlinearSimulations.plot_saturatedfluxVsParameter import plot_saturatedfluxVsParameter

#######################################################################
#   CLASS TO PLOT THE INFLUENCE OF A PARAMETER ON THE SATURATED FLUX
#######################################################################
   
class PlotSaturatedFluxVersusParameter: 
    ''' Holds and manipulates the Canvas and Axis object for the following plots:
            (1) ... '''

    def __init__(self, tab, notebook_tab):

        # Set an identifier for this plotting class, used as the figure name
        self.identifier = "TabNonlinearTime:PlotSaturatedFluxVersusParameter"
        
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
            self.Axis.update_gridSpeficications(top=0.9, left=0.15, right=0.95, bottom=0.12)
            
            # Set the (x,y) quantities and titles for the canvas  
            self.Axis.set_axisQuantities(1, "parameter", "qflux", name = "Heat flux")
            self.Axis.set_axisQuantities(2, "parameter", "vflux", name = "Momentum flux")
            self.Axis.set_axisQuantities(3, "parameter", "pflux", name = "Particle flux") 

            # Create the frames and variables for the plotting arguments
            self.initiate_optionVariablesAndFrames(notebook_tab)
        
            # Tell the <Axis> which quantity to plot
            self.Axis.set_plot(1) 
            self.Axis.update_labelFrame(self.Canvas.frame)
            
            # Remember that the <Canvas> is initiaited
            self.initiated_canvas = True
        return     

    #---------------------- 
    def get_argumentsForPlot(self, Canvas):
        ''' Gather the arguments for <plot_quantityVsParameter>. '''
        self.Axis.update_axisQuantities()
        plotting_arguments = {
            # Specify which simulations to plot
                "research"      : self.root.Research,\
                "experiment_id" : self.tab.Simulations.experiment_id,\
                "simulation_id" : self.tab.Simulations.simulation_id,\
                "species"       : self.tab.Simulations.species,\
            # Specify which data to plot
                "knob"          : self.knob,\
                "key"           : self.key,\
                "x_quantity"    : Canvas.Axis.x_quantity,\
                "y_quantity"    : Canvas.Axis.y_quantity,\
            # Specify data range
                "x_range"       : Canvas.Axis.x_range,\
                "y_range"       : Canvas.Axis.y_range,
                "units"         : Canvas.Axis.units,\
            # Labels
                "x_label"       : Canvas.Axis.x_label,\
                "y_label"       : Canvas.Axis.y_label,\
                "title"         : Canvas.Axis.title,\
            # Options
                "normalize"     : self.normalize,\
                "show_errorStd" : self.show_errorStd,\
                "show_errorMinMax" : self.show_errorMinMax,\
                "splitInKineticAdiabatic" : self.splitInKineticAdiabatic,\
                "splitInKineticAdiabaticTEM" : self.splitInKineticAdiabaticTEM,\
            # For the GUI the figure object already exists 
                "show_figure"   : False,\
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
                plot_saturatedfluxVsParameter(**plotting_arguments)
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
        self.knob = "-"
        self.key  = "-"
        self.show_errorStd = True
        self.show_errorMinMax = False
        self.normalize = False
        self.splitInKineticAdiabatic = False
        self.splitInKineticAdiabaticTEM = False
        self.includeFluxNorm = False
       
        # Create a subframe for each set of plotting options
        self.subframe_yQuantity = ttk.Frame(frame)
        self.subframe_parameter = ttk.Frame(frame)
        self.subframe_options = ttk.Frame(frame)
         
        # Configure the rows/columns of the <frame> to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1,  uniform="yes")
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=3,  uniform="yes")
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_yQuantity
        tk.Grid.columnconfigure(frame, 1, weight=0) # subframe_parameter
        tk.Grid.columnconfigure(frame, 2, weight=0) # subframe_options 
 
        # Add the subframes to the main frame
        self.subframe_yQuantity.grid(  row=1, column=0, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_parameter.grid(  row=1, column=1, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_options.grid(    row=1, column=2, padx=25, pady=(0,0), stick='NSEW')
         
        # Fill the frames with their options
        self.YQuantity = self.Frame_YQuantity(self, self.subframe_yQuantity)
        self.Parameter = self.Frame_Parameter(self, self.subframe_parameter)
        self.Options   = self.Frame_Options(self, self.subframe_options)
        return     

    #----------------------  
    class Frame_YQuantity:
         
        def __init__(self, PlotSaturatedFluxVersusParameter, frame):
             
            # Make the tabs available
            self.PlotSaturatedFluxVersusParameter = PlotSaturatedFluxVersusParameter 
            self.tab = PlotSaturatedFluxVersusParameter.tab
            
            # Get the possible y-quantities that were set in <initiate_canvas>
            plot_ids, axis_quantities = self.PlotSaturatedFluxVersusParameter.Axis.get_axisQuantities()
            self.radiobuttons = [None]*len(plot_ids)

            # Variables
            self.var_yQuantity = tk.IntVar(value=self.PlotSaturatedFluxVersusParameter.Axis.get_plotid())
            
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
            return 

        #---------------------
        def change_yQuantity(self, *_):
            ''' Update the (x,y) quantities, reset the axis since the data has
            changed, and change the title of the label frame of the canvas. ''' 
            self.PlotSaturatedFluxVersusParameter.Axis.set_plot(self.var_yQuantity.get()) 
            self.PlotSaturatedFluxVersusParameter.Axis.reset_axis()
            self.PlotSaturatedFluxVersusParameter.Axis.update_labelFrame(self.PlotSaturatedFluxVersusParameter.Canvas.frame)
            self.PlotSaturatedFluxVersusParameter.plotOnTheGUI()
            return
        
    #----------------------
    class Frame_Parameter:
         
        def __init__(self, PlotSaturatedFluxVersusParameter, frame):
 
            # Make the parents available
            self.PlotSaturatedFluxVersusParameter = PlotSaturatedFluxVersusParameter  
            self.root = PlotSaturatedFluxVersusParameter.tab.root 
            self.tab = PlotSaturatedFluxVersusParameter.tab
            
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
 
            # Initiate the stella <knob> and <key> for the selected <parameter>
            self.PlotSaturatedFluxVersusParameter.knob = standardParameters[self.var_parameter.get()]["knob"]
            self.PlotSaturatedFluxVersusParameter.key  = standardParameters[self.var_parameter.get()]["key"]
            return 
        
        #---------------------- 
        def update_parameter(self, calledFromGui=True, *_):
             
            # Read which parameter is selected
            parameter = self.var_parameter.get()
            if parameter=='nfield_periods': parameter = "nfield"
             
            # Select the stella <knob> and <key> for the selected <parameter>
            self.PlotSaturatedFluxVersusParameter.knob = standardParameters[parameter]["knob"]
            self.PlotSaturatedFluxVersusParameter.key  = standardParameters[parameter]["key"]

            # Set the (x,y) quantities for the axis and the canvas
            plot_ids, axis_quantities = self.tab.PlotTimeEvolution.Axis.get_axisQuantities()
            for plot_id in plot_ids: self.PlotSaturatedFluxVersusParameter.Axis.set_axisQuantities(plot_id, parameter, axis_quantities[plot_id]['y']) 
            self.PlotSaturatedFluxVersusParameter.Axis.update_axisQuantities()
            
            # Update the label in the label frame
            self.PlotSaturatedFluxVersusParameter.Axis.update_labelFrame(self.PlotSaturatedFluxVersusParameter.Canvas.frame)
            
            # Reset the axis when we change the parameter plotted on the x-axis
            if calledFromGui: self.PlotSaturatedFluxVersusParameter.Axis.reset_axis() 
            if calledFromGui: self.PlotSaturatedFluxVersusParameter.plotOnTheGUI()
            return

        #----------------------------
        def update_parameterWithoutPlotting(self, *_):
            self.update_parameter(calledFromGui=False)
            return
    
    #----------------------    
    class Frame_Options:
         
        def __init__(self, PlotSaturatedFluxVersusParameter, frame):
 
            # Make the parents available
            self.PlotSaturatedFluxVersusParameter = PlotSaturatedFluxVersusParameter  
            self.tab = PlotSaturatedFluxVersusParameter.tab
            
            # Variables  
            self.var_plotErrorBars1 = tk.IntVar(value=self.PlotSaturatedFluxVersusParameter.show_errorStd)
            self.var_plotErrorBars2 = tk.IntVar(value=self.PlotSaturatedFluxVersusParameter.show_errorMinMax)
            self.var_normalize = tk.IntVar(value=self.PlotSaturatedFluxVersusParameter.normalize)
            self.var_splitInKineticAdiabatic = tk.IntVar(value=0)
            self.var_splitInKineticAdiabaticTEM = tk.IntVar(value=0)
            self.var_includefluxnorm = tk.IntVar(value=0)
               
            # Create the widgets
            self.lbl_options = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_error1 = ttk.Checkbutton(frame, text=" Include error bars (std)")
            self.chk_error2 = ttk.Checkbutton(frame, text=" Include error bars (min/max)")
            self.chk_norm   = ttk.Checkbutton(frame, text=" Normalize by saturation")
            self.chk_split = ttk.Checkbutton(frame, text=" Split in kinetic/adiabatic")
            self.chk_splitTEM = ttk.Checkbutton(frame, text=" Split in kinetic/adiabatic/TEM")
            self.chk_includefluxnorm = ttk.Checkbutton(frame, text=" Include 1/|<nabla r>|")
             
            # Add the commands
            self.chk_error1.config(variable=self.var_plotErrorBars1, command=self.update_error1)
            self.chk_error2.config(variable=self.var_plotErrorBars2, command=self.update_error2)
            self.chk_norm.config(variable=self.var_normalize, command=self.update_norm)
            self.chk_split.config(variable=self.var_splitInKineticAdiabatic, command=self.update_split)
            self.chk_splitTEM.config(variable=self.var_splitInKineticAdiabaticTEM, command=self.update_splitTEM)
            self.chk_includefluxnorm.config(variable=self.var_includefluxnorm, command=self.update_includefluxnorm)
             
            # Configure the frame
            tk.Grid.columnconfigure(frame, 0, weight=1) 
             
            # Add the options to the frame
            self.lbl_options.grid(row=0, column=0, **PAD_TITLE2)
            self.chk_error1.grid( row=1, column=0, **PAD_LABEL2)
            self.chk_error2.grid( row=2, column=0, **PAD_LABEL2)
            self.chk_norm.grid(   row=3, column=0, **PAD_LABEL2)
            self.chk_split.grid(  row=1, column=1, **PAD_LABEL2)
            self.chk_splitTEM.grid(row=2, column=1, **PAD_LABEL2)
            self.chk_includefluxnorm.grid(row=3, column=1, **PAD_LABEL2)
 
        #----------------------  
        def update_error1(self, *_): 
            if self.var_plotErrorBars1.get()==0: self.PlotSaturatedFluxVersusParameter.show_errorStd = False
            if self.var_plotErrorBars1.get()==1: self.PlotSaturatedFluxVersusParameter.show_errorStd = True
            self.PlotSaturatedFluxVersusParameter.plotOnTheGUI()
     
        #----------------------  
        def update_error2(self, *_): 
            if self.var_plotErrorBars2.get()==0: self.PlotSaturatedFluxVersusParameter.show_errorMinMax = False
            if self.var_plotErrorBars2.get()==1: self.PlotSaturatedFluxVersusParameter.show_errorMinMax = True
            self.PlotSaturatedFluxVersusParameter.plotOnTheGUI()
        
        #----------------------  
        def update_norm(self, *_): 
            if self.var_normalize.get()==0: self.PlotSaturatedFluxVersusParameter.normalize = False
            if self.var_normalize.get()==1: self.PlotSaturatedFluxVersusParameter.normalize = True
            self.PlotSaturatedFluxVersusParameter.Axis.reset_axis() 
            self.PlotSaturatedFluxVersusParameter.plotOnTheGUI()
 
        #----------------------  
        def update_split(self, *_): 
            if self.var_splitInKineticAdiabatic.get()==0: self.PlotSaturatedFluxVersusParameter.splitInKineticAdiabatic = False
            if self.var_splitInKineticAdiabatic.get()==1: self.PlotSaturatedFluxVersusParameter.splitInKineticAdiabatic = True
            self.PlotSaturatedFluxVersusParameter.plotOnTheGUI()
            
        #----------------------  
        def update_splitTEM(self, *_): 
            if self.var_splitInKineticAdiabaticTEM.get()==0: self.PlotSaturatedFluxVersusParameter.splitInKineticAdiabaticTEM = False
            if self.var_splitInKineticAdiabaticTEM.get()==1: self.PlotSaturatedFluxVersusParameter.splitInKineticAdiabaticTEM = True
            self.PlotSaturatedFluxVersusParameter.plotOnTheGUI()
               
        #----------------------  
        def update_includefluxnorm(self, *_): 
            if self.var_includefluxnorm.get()==0: self.PlotSaturatedFluxVersusParameter.includeFluxNorm = False
            if self.var_includefluxnorm.get()==1: self.PlotSaturatedFluxVersusParameter.includeFluxNorm = True 
            if self.var_includefluxnorm.get()==0: self.tab.PlotTimeEvolution.includeFluxNorm = False
            if self.var_includefluxnorm.get()==1: self.tab.PlotTimeEvolution.includeFluxNorm = True 
            self.PlotSaturatedFluxVersusParameter.plotOnTheGUI()
               
            