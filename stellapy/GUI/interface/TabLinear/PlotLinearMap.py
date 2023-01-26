
import tkinter as tk
from tkinter import ttk  
from stellapy.GUI.widgets import Canvas, PAD_TITLE2, PAD_LABEL2
from stellapy.plot.utils.labels.standardParameters import standardParameters
from stellapy.plot.utils.labels.standardParameters import standardParametersInOrder 
#from stellapy.GUI.plot.linearSimulations.plot_quantityVsParameterVsParameter import plot_quantityVsParameterVsParameter

#################################################################
#                 CLASS TO PLOT THE LINEAR MAP
#################################################################

class PlotLinearMap: 
    ''' Holds and manipulates the Canvas and Axis object for the following plots:
            (1) ... '''
    
    def __init__(self, tab, notebook_tab):

        # Set an identifier for this plotting class, used as the figure name
        self.identifier = "TabLinear:PlotLinearMap"
        
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
            self.Axis.update_gridSpeficications(top=0.94, left=0.1, right=0.95, bottom=0.12)
            
            # Set the (x,y) quantities and titles for the canvas
            self.Axis.set_axisQuantities(1, "parameter", "parameter", "gamma")
            self.Axis.set_axisQuantities(2, "parameter", "parameter", "omega")
            self.Axis.set_axisQuantities(3, "parameter", "parameter", "ky") 
            
            # Tell the <Axis> which quantity to plot
            self.Axis.set_plot(self.ZQuantity.var_zQuantity.get())
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
                "research"       : self.root.Research,\
                "experiment_id"  : self.tab.Simulations.experiment_id,\
                "simulation_id"  : self.tab.Simulations.simulation_id,\
            # Speficy which data to plot
                "parameter_knob1": self.knob1,\
                "parameter_key1" : self.key1,\
                "parameter_knob2": self.knob2,\
                "parameter_key2" : self.key2,\
                "k_value"        : self.k_value,\
                "z_quantity"     : Canvas.Axis.z_quantity,\
            # Specify data range
                "x_range"        : Canvas.Axis.x_range,\
                "y_range"        : Canvas.Axis.y_range,\
                "kx_range"       : Canvas.Axis.kx_range,\
                "ky_range"       : Canvas.Axis.ky_range,\
                "units"          : Canvas.Axis.units,\
                "lineardata"     : "average",\
                "k_value"        : self.k_value,\
            # Labels
                "x_label"        : Canvas.Axis.x_label,\
                "y_label"        : Canvas.Axis.y_label,\
                "title"          : Canvas.Axis.title,\
            # Options
                "interpolate"    : self.interpolate,\
                "step"           : self.step,\
                "log"            : self.log,\
                "dynamickyrange"  : self.tab.PlotLinearSpectrum.dynamickyrange,\
            # Save time by reading the data elsewhere
                "parameters1"    : self.root.Research.data['linear_map'][self.key1],\
                "parameters2"    : self.root.Research.data['linear_map'][self.key2],\
                "gamma_mostUnstableMode" : self.root.Research.data['linear_map']['gamma'],\
                "omega_mostUnstableMode" : self.root.Research.data['linear_map']['omega'],\
                "ky_mostUnstableMode" : self.root.Research.data['linear_map']['ky'],\
            # For the GUI the figure object already exists 
                "ax"             : Canvas.Axis.ax,\
                "Progress"       : self.tab.Progress,\
                "root"           : self.root}
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
                
                # Remove the colorbar
                try: self.colorbar.remove()
                except: pass
                
                # Make sure the saved linear map has the correct keys
                keys = self.root.Research.data['linear_map'].keys()
                for key in ['key1', 'key2', 'gamma', 'omega', 'ky', self.key1, self.key2]: 
                    if key not in keys:
                        self.root.Research.data['linear_map'][key] = None
                
                # Be able to flip the parameter axis if the data is saved
                if  self.key1==self.root.Research.data['linear_map']['key2']\
                and self.key2==self.root.Research.data['linear_map']['key1']:
                    self.root.Research.data['linear_map']['key1']  = self.key1
                    self.root.Research.data['linear_map']['key2']  = self.key2
                    self.root.Research.data['linear_map']['gamma'] = self.root.Research.data['linear_map']['gamma'].transpose()
                    self.root.Research.data['linear_map']['omega'] = self.root.Research.data['linear_map']['omega'].transpose()
                    self.root.Research.data['linear_map']['ky']    = self.root.Research.data['linear_map']['ky'].transpose()
                
                # Plot the linear map
                Canvas.Axis.start_plotting(research_arguments, input_files)
                self.colorbar, data = plot_quantityVsParameterVsParameter(**plotting_arguments)
                Canvas.Axis.finish_plotting(research_arguments, input_files)
                Canvas.Axis.save_plottingArguments(self.get_argumentsForPlot(Canvas)) 
        
                # Save the data to the research object
                if data is not None: self.root.Research.data['linear_map'] = data
            
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
        self.z_quantity = "gamma"
        self.k_value = "max"
        self.interpolate = False
        self.step  = 4
        self.knob1 = "-"
        self.key1  = "-"
        self.knob2 = "-"
        self.key2  = "-"
        self.colorbar = None
        self.log = False
        
        # Make sure there is a linear_map object 
        if 'linear_map' not in self.root.Research.data.keys():
            self.root.Research.data['linear_map'] = {
                'keys' : ['-', '-'],\
                'gamma' : None,\
                'omega' : None,\
                'ky' : None} 
        
        # Create a subframe for each set of plotting options
        self.subframe_zQuantity = ttk.Frame(frame)
        self.subframe_parameter = ttk.Frame(frame)
        self.subframe_extraction = ttk.Frame(frame)
        self.subframe_options = ttk.Frame(frame)
        
        # Configure the main frame to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1)      
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=1)     
        tk.Grid.rowconfigure(   frame, 3, weight=1)     
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_yQuantity
        tk.Grid.columnconfigure(frame, 1, weight=0) # subframe_parameter
        tk.Grid.columnconfigure(frame, 2, weight=0) # subframe_extraction
        tk.Grid.columnconfigure(frame, 3, weight=0) # subframe_options

        # Add the subframes to the main frame
        self.subframe_zQuantity.grid( row=1, column=0, padx=25, pady=(0,0), stick='NSEW', rowspan=2)
        self.subframe_parameter.grid( row=1, column=1, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_extraction.grid(row=1, column=2, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_options.grid(   row=1, column=3, padx=25, pady=(0,0), stick='NSEW', rowspan=2)
        
        # Fill the frames with their options
        self.ZQuantity  = self.Frame_ZQuantity( self, self.subframe_zQuantity)
        self.Parameter  = self.Frame_Parameter( self, self.subframe_parameter)
        self.Extraction = self.Frame_Extraction(self, self.subframe_extraction)
        self.Options    = self.Frame_Options(   self, self.subframe_options)
        return   
    
    class Frame_ZQuantity:
        
        def __init__(self, PlotLinearMap, frame):
            
            # Make the parents available
            self.PlotLinearMap = PlotLinearMap  
            
            # Variables
            self.var_zQuantity = tk.IntVar(value=1)
        
            # Configure the frame <self.subframe_yQuantity>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Create the widgets
            self.lbl_zdata = ttk.Label(frame, text="Z-data", style='prefTitle.TLabel')
            self.rbn_gamma = ttk.Radiobutton(frame, text='  Growth rate')
            self.rbn_omega = ttk.Radiobutton(frame, text='  Frequency')
            self.rbn_ky    = ttk.Radiobutton(frame, text='  Wavenumber')
    
            # Add the values and commands to the radiobuttons
            self.rbn_gamma.config(value=1, variable=self.var_zQuantity, command=self.change_zQuantity)
            self.rbn_omega.config(value=2, variable=self.var_zQuantity, command=self.change_zQuantity)
            self.rbn_ky.config(   value=3, variable=self.var_zQuantity, command=self.change_zQuantity)
            
            # Add the options to the frame
            self.lbl_zdata.grid(row=0, column=0, **PAD_TITLE2)
            self.rbn_gamma.grid(row=1, column=0, **PAD_LABEL2)
            self.rbn_omega.grid(row=2, column=0, **PAD_LABEL2)
            self.rbn_ky.grid(   row=3, column=0, **PAD_LABEL2)
            return 
        
        #---------------------
        def change_zQuantity(self, *_):
            ''' Update the (z) quantity, reset the axis since the data has
            changed, and change the header of the label frame of the canvas. ''' 
            self.PlotLinearMap.Axis.set_plot(self.var_zQuantity.get())
            self.PlotLinearMap.Axis.reset_axis()
            self.PlotLinearMap.Axis.update_labelFrame(self.PlotLinearMap.Canvas.frame)
            self.PlotLinearMap.plotOnTheGUI()
            return
    
    #---------------------
    class Frame_Parameter:
        
        def __init__(self, PlotLinearMap, frame):

            # Make the parents available
            self.root = PlotLinearMap.tab.root
            self.PlotLinearMap = PlotLinearMap  
            
            # Variables 
            self.var_parameter1 = tk.StringVar(value=standardParametersInOrder[0])
            self.var_parameter2 = tk.StringVar(value=standardParametersInOrder[0])
        
            # Configure the frame <self.subframe_parameter>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Choose the parameter for the x-axis and y-axis
            self.lbl_par  = ttk.Label(frame, text="Parameters", style='prefTitle.TLabel'); width=10
            self.mnu_par1 = ttk.OptionMenu(frame, self.var_parameter1, standardParametersInOrder[0], *standardParametersInOrder, style='option.TMenubutton')
            self.mnu_par1["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
            self.mnu_par1.config(width=width)
            self.mnu_par2 = ttk.OptionMenu(frame, self.var_parameter2, standardParametersInOrder[0], *standardParametersInOrder, style='option.TMenubutton')
            self.mnu_par2["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
            self.mnu_par2.config(width=width) 
            self.var_parameter1.trace_id = self.var_parameter1.trace('w', self.update_parameter) 
            self.var_parameter2.trace_id = self.var_parameter2.trace('w', self.update_parameter) 
        
            # Add the widgets to the frame
            self.lbl_par.grid( row=0, column=0, **PAD_TITLE2)
            self.mnu_par1.grid(row=1, column=0, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
            self.mnu_par2.grid(row=2, column=0, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
            return

        #---------------------- 
        def update_parameter(self, calledFromGui=True, *_):
            ''' The parameter plotted on the x-axis can be chose from a dropdown
            menu. When a parameter is selected, read the stella <knob> and <key> 
            from the <standardParameters> dictionary, update the x-quantity of
            the 3 possibly plots and reset the axis of the maplotlib figure. '''
            
            # Read which parameter is selected
            parameter1 = self.var_parameter1.get()
            parameter2 = self.var_parameter2.get()
            
            # Select the stella <knob> and <key> for the selected <parameter>
            self.PlotLinearMap.knob1 = standardParameters[parameter1]["knob"]
            self.PlotLinearMap.key1  = standardParameters[parameter1]["key"]
            self.PlotLinearMap.knob2 = standardParameters[parameter2]["knob"]
            self.PlotLinearMap.key2  = standardParameters[parameter2]["key"]

            # Set the (x,y) quantities for the axis and the canvas 
            self.PlotLinearMap.Axis.set_axisQuantities(1, parameter1, parameter2, "gamma")
            self.PlotLinearMap.Axis.set_axisQuantities(2, parameter1, parameter2, "omega")
            self.PlotLinearMap.Axis.set_axisQuantities(3, parameter1, parameter2, "ky") 
            self.PlotParameterInfluence.Axis.update_axisQuantities()
            
            # Update the label in the label frame
            self.PlotParameterInfluence.Axis.update_labelFrame(self.PlotParameterInfluence.Canvas.frame)
            
            # Reset the axis when we change the parameter plotted on the x-axis
            if calledFromGui: self.PlotLinearMap.Axis.reset_axis() 
            return

        #----------------------------
        def update_parameterWithoutPlotting(self, *_):
            if self.PlotLinearMap.initiated_canvas: self.update_parameter(calledFromGui=False)
            return
        
    #----------------------     
    class Frame_Extraction:
        
        def __init__(self, PlotLinearMap, frame):

            # Make the parents available 
            self.root = PlotLinearMap.tab.root 
            self.PlotLinearMap = PlotLinearMap 
            
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
            if self.var_extraction.get()==1: self.PlotLinearMap.k_value = "max"
            if self.var_extraction.get()==2: self.PlotLinearMap.k_value = float(self.var_specificKvalue.get())
            if calledFromGui: self.PlotLinearMap.plotOnTheGUI() 

        #----------------------------
        def update_selectedModeWithoutPlotting(self, *_):
            self.update_selectedMode(calledFromGui=False)
        
        #---------------------- 
        def change_extraction(self):
            if self.var_extraction.get()==1: self.PlotLinearMap.k_value = "max"
            if self.var_extraction.get()==2: self.PlotLinearMap.k_value = float(self.var_specificKvalue.get())
            self.PlotLinearMap.plotOnTheGUI()
            return

    #---------------------
    class Frame_Options:
        
        def __init__(self, PlotLinearMap, frame):
        
            # Make the parents available
            self.root = PlotLinearMap.tab.root
            self.PlotLinearMap = PlotLinearMap
        
            # Variables
            self.var_interpolate = tk.IntVar(value=self.PlotLinearMap.interpolate)
            self.var_step = tk.StringVar(value=self.PlotLinearMap.step)
            self.var_log = tk.IntVar(value=self.PlotLinearMap.log)
        
            # Create the widgets
            self.lbl_options = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_interp = ttk.Checkbutton(frame, text=" Interpolate z-data with step = ")
            self.chk_log = ttk.Checkbutton(frame, text=" Logaritmic z-axis")
            self.btn_rewrite = ttk.Button(frame, text="Reread linear map", width=20)
    
            # Add the commands
            self.chk_interp.config(variable=self.var_interpolate, command=self.update_interp)
            self.chk_log.config(variable=self.var_log, command=self.update_log)
            self.btn_rewrite.config(command=lambda: self.rewrite_linearMap())
            
            # Drop down menu to choose the interpolation step 
            self.ent_interp = ttk.Entry(frame, textvariable=self.var_step, font=("Courier New", 11), style='opt_valueCBold.TEntry', width=5)
            self.ent_interp.bind('<Return>', self.update_step) 
            
            # Configure the frame
            tk.Grid.rowconfigure(frame, 0, weight=0) 
            tk.Grid.rowconfigure(frame, 1, weight=0) 
            tk.Grid.columnconfigure(frame, 0, weight=0) 
            tk.Grid.columnconfigure(frame, 1, weight=0) 
            tk.Grid.columnconfigure(frame, 2, weight=1) 
            
            # Add the options to the frame
            self.lbl_options.grid(row=0, column=0, **PAD_TITLE2, columnspan=2)
            self.chk_log.grid(    row=1, column=0, **PAD_LABEL2, columnspan=2)
            self.chk_interp.grid( row=2, column=0, **PAD_LABEL2)
            self.ent_interp.grid( row=2, column=1)
            self.btn_rewrite.grid(row=3, column=0, **PAD_LABEL2, columnspan=2) 
            return 

        #-------------------
        def update_log(self, *args): 
            if self.var_log.get()==0: self.PlotLinearMap.log = False
            if self.var_log.get()==1: self.PlotLinearMap.log = True
            self.PlotLinearMap.plotOnTheGUI()
            return 
                   
        #-------------------
        def update_interp(self, *args): 
            if self.var_interpolate.get()==0: self.PlotLinearMap.interpolate = False
            if self.var_interpolate.get()==1: self.PlotLinearMap.interpolate = True
            self.PlotLinearMap.plotOnTheGUI()
            return 
                           
        #-------------------  
        def update_step(self, *args):
            self.PlotLinearMap.step = int(self.var_step.get()) 
            self.PlotLinearMap.plotOnTheGUI()
            return 
        
        #----------------------
        def rewrite_linearMap(self):
            self.root.Research.data['linear_map'] = {
                'parameters1' : None,\
                'parameters2' : None,\
                'gamma' : None,\
                'omega' : None,\
                'ky' : None} 
            self.PlotLinearMap.plotOnTheGUI()
        