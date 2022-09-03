
import tkinter as tk
from tkinter import ttk  
from stellapy.plot.utils.labels import standardParameters, standardParametersInOrder 
from stellapy.GUI.widgets import Canvas, PAD_TITLE2, PAD_LABEL2 
from stellapy.GUI.plot.nonlinearSimulations.plot_saturatedFluxVsParameterVsParameter import plot_saturatedFluxVsParameterVsParameter

#################################################################
#                 CLASS TO PLOT THE LINEAR MAP
#################################################################

class PlotNonlinearMap: 
    ''' Holds and manipulates the Canvas and Axis object for the following plots:
            (1) ... '''
    
    def __init__(self, tab, notebook_tab):

        # Set an identifier for this plotting class, used as the figure name
        self.identifier = "TabNonlinearTime:PlotNonlinearMap"
        
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
            self.Axis.update_gridSpeficications(top=0.9, left=0.1, right=0.95, bottom=0.12)
            
            # Set the (x,y) quantities and titles for the canvas (copy from <PlotTimeEvolution>)
            self.Axis.set_axisQuantities(1,  "parameter", "parameter", "qflux", " Saturated heat flux")
            self.Axis.set_axisQuantities(2,  "parameter", "parameter", "qflux_range", " max(Q(kx or ky))-min(Q(kx or ky))")
            self.Axis.set_axisQuantities(3,  "parameter", "parameter", "phi2_range", " max(Phi²(kx or ky))-min(Phi²(kx or ky))")

            # Create the frames and variables for the plotting arguments
            self.initiate_optionVariablesAndFrames(notebook_tab)
            
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
                "research"        : self.root.Research,\
                "experiment_id"   : self.tab.Simulations.experiment_id,\
                "simulation_id"   : self.tab.Simulations.simulation_id,\
                "species"         : self.tab.Simulations.species,\
            # Speficy which data to plot
                "parameter_knob1" : self.knob1,\
                "parameter_key1"  : self.key1,\
                "parameter_knob2" : self.knob2,\
                "parameter_key2"  : self.key2,\
                "x_quantity"      : Canvas.Axis.x_quantity,\
                "y_quantity"      : Canvas.Axis.y_quantity,\
                "z_quantity"      : Canvas.Axis.z_quantity,\
            # Specify data range
                "x_range"       : Canvas.Axis.x_range,\
                "y_range"       : Canvas.Axis.y_range,
                "units"         : Canvas.Axis.units,\
            # Labels
                "x_label"       : Canvas.Axis.x_label,\
                "y_label"       : Canvas.Axis.y_label,\
                "title"         : Canvas.Axis.title,\
            # Options
                "interpolate"   : self.interpolate,\
                "step"          : self.step,\
                "log"           : self.log,\
                "removeZonal"   : self.removeZonal,\
                "spectrum_along_axis" : self.spectrum_along_axis,\
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

                # Try to remove the colorbar since each plot will add an additional one
                try: self.colorbar.remove()
                except: pass
                
                # Plot the parallel mode structure 
                Canvas.Axis.start_plotting(research_arguments, input_files)
                self.colorbar = plot_saturatedFluxVsParameterVsParameter(**plotting_arguments)
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
        self.knob1 = "-"
        self.key1  = "-"
        self.knob2 = "-"
        self.key2  = "-"
        self.spectrum_along_axis = "kx"
        self.removeZonal = True
        self.z_quantity = "qflux"
        
        # Extra options
        self.interpolate = False
        self.step  = 4
        self.colorbar = None
        self.log = False
         
        # Create a subframe for each set of plotting options 
        self.subframe_zQuantity = ttk.Frame(frame) 
        self.subframe_parameter = ttk.Frame(frame) 
        self.subframe_options = ttk.Frame(frame)
         
        # Configure the rows/columns of the <frame> to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1,  uniform="yes")
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=2,  uniform="yes")  
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_zQuantity 
        tk.Grid.columnconfigure(frame, 1, weight=0) # subframe_parameter 
        tk.Grid.columnconfigure(frame, 2, weight=0) # subframe_options
 
        # Add the subframes to the main frame 
        self.subframe_zQuantity.grid( row=1, column=0, padx=25, pady=(0,0), stick='NSEW') 
        self.subframe_parameter.grid( row=1, column=1, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_options.grid(   row=1, column=2, padx=25, pady=(0,0), stick='NSEW')
         
        # Fill the frames with their options 
        self.ZQuantity  = self.Frame_ZQuantity( self, self.subframe_zQuantity) 
        self.Parameter  = self.Frame_Parameter( self, self.subframe_parameter) 
        self.Options    = self.Frame_Options(   self, self.subframe_options)
        return     
 
    #----------------------
    class Frame_ZQuantity:
         
        def __init__(self, PlotNonlinearMap, frame):
             
            # Make the parents available
            self.PlotNonlinearMap = PlotNonlinearMap  
         
            # Variables
            self.var_zQuantity = tk.IntVar(value=1)

            # Get the possible z-quantities that were set in <initiate_canvas>
            plot_ids, axis_quantities = self.PlotNonlinearMap.Axis.get_axisQuantities()
            self.radiobuttons = [None]*len(plot_ids)
            
            # Configure the frame <self.subframe_zQuantity>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
             
            # Create the widgets
            self.lbl_zQuantity = ttk.Label(frame, text="Quantity on the z-axis", style='prefTitle.TLabel')
            for i, plot_id in enumerate(plot_ids): self.radiobuttons[i] = ttk.Radiobutton(frame, text=axis_quantities[plot_id]['name']) 
            for i, plot_id in enumerate(plot_ids): self.radiobuttons[i].config(value=plot_id, variable=self.var_zQuantity, command=self.change_zQuantity)
             
            # Add the options to the frame
            self.lbl_zQuantity.grid( row=1, column=0, **PAD_TITLE2)
            for i, plot_id in enumerate(plot_ids): self.radiobuttons[i].grid(row=2+i, column=0, **PAD_LABEL2) 
            return 
         
        #---------------------
        def change_zQuantity(self, *args):
            ''' Update the (x,y) quantities, reset the axis since the data has
            changed, and change the title of the label frame of the canvas. ''' 
            self.PlotNonlinearMap.Axis.set_plot(self.var_zQuantity.get())
            self.PlotNonlinearMap.Axis.reset_axis()
            self.PlotNonlinearMap.Axis.update_labelFrame(self.PlotNonlinearMap.Canvas.frame)
            self.PlotNonlinearMap.plotOnTheGUI()
            return
    
    #----------------------
    class Frame_Parameter:
         
        def __init__(self, PlotNonlinearMap, frame):
 
            # Make the parents available
            self.PlotNonlinearMap = PlotNonlinearMap  
            self.root = PlotNonlinearMap.root
            self.tab = PlotNonlinearMap.tab
            
            # Variables
            self.var_parameter1 = tk.StringVar(value=standardParametersInOrder[0])
            self.var_parameter2 = tk.StringVar(value=standardParametersInOrder[0])
            self.var_parameter3 = tk.StringVar(value=["kx", "ky"][0])
        
            # Configure the frame <self.subframe_parameter>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
             
            # Choose the parameter for the x-axis and y-axis
            self.lbl_par  = ttk.Label(frame, text="Parameters", style='prefTitle.TLabel'); width=10
            self.mnu_par1 = ttk.OptionMenu(frame, self.var_parameter1, standardParametersInOrder[0], *standardParametersInOrder, style='option.TMenubutton')
            self.mnu_par2 = ttk.OptionMenu(frame, self.var_parameter2, standardParametersInOrder[0], *standardParametersInOrder, style='option.TMenubutton')
            self.mnu_par3 = ttk.OptionMenu(frame, self.var_parameter3, ["kx", "ky"][0], *["kx", "ky"], style='option.TMenubutton')
            self.mnu_par1["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
            self.mnu_par2["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
            self.mnu_par3["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])           
            self.mnu_par1.config(width=width); self.mnu_par2.config(width=width); self.mnu_par3.config(width=width)
            def update_parameter1(*_): self.update_parameter(i=1); return
            def update_parameter2(*_): self.update_parameter(i=2); return 
            def update_parameter3(*_): self.update_parameter(i=3); return 
            self.var_parameter1.trace('w', update_parameter1) 
            self.var_parameter2.trace('w', update_parameter2) 
            self.var_parameter3.trace('w', update_parameter3) 
         
            # Add the widgets to the frame
            self.lbl_par.grid( row=0, column=0, **PAD_TITLE2)
            self.mnu_par1.grid(row=1, column=0, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
            self.mnu_par2.grid(row=2, column=0, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
            self.mnu_par3.grid(row=3, column=0, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
            return
 
        #-----------------------------
        def update_parameter(self, i):
            
            # Assign the correct parameter
            if i in [1,2]:
                
                # Read which parameter is selected
                parameter1 = self.var_parameter1.get()
                parameter2 = self.var_parameter2.get()
             
                # Select the stella <knob> and <key> for the selected <parameter>
                self.PlotNonlinearMap.knob1 = standardParameters[parameter1]["knob"]
                self.PlotNonlinearMap.key1  = standardParameters[parameter1]["key"]
                self.PlotNonlinearMap.knob2 = standardParameters[parameter2]["knob"]
                self.PlotNonlinearMap.key2  = standardParameters[parameter2]["key"]

                # Set the (x,y) quantities for the axis and the canvas
                plot_ids, axis_quantities = self.tab.PlotNonlinearMap.Axis.get_axisQuantities()
                for plot_id in plot_ids: self.PlotNonlinearMap.Axis.set_axisQuantities(plot_id, parameter1, parameter2, axis_quantities[plot_id]['z']) 
                self.PlotNonlinearMap.Axis.update_axisQuantities()
                
            # Change the sum alongst kx or ky
            if i==3: self.PlotNonlinearMap.spectrum_along_axis = self.var_parameter3.get() 
                
            # Update the label in the label frame, reset the axis and replot
            self.PlotNonlinearMap.Axis.update_labelFrame(self.PlotNonlinearMap.Canvas.frame) 
            self.PlotNonlinearMap.Axis.reset_axis() 
            self.PlotNonlinearMap.plotOnTheGUI()
            return
 
    #-----------------------------
    class Frame_Options:
         
        def __init__(self, PlotNonlinearMap, frame):
         
            # Make the parents available
            self.PlotNonlinearMap = PlotNonlinearMap  
            
            # Variables
            self.var_interpolate = tk.IntVar(value=0)
            self.var_step = tk.StringVar(value="4")
            self.var_log = tk.IntVar(value=self.PlotNonlinearMap.log)
            self.var_removeZonal = tk.IntVar(value=self.PlotNonlinearMap.removeZonal)
        
            # Create the widgets
            self.lbl_options = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_interp = ttk.Checkbutton(frame, text=" Interpolate z-data with step = ") 
            self.chk_log = ttk.Checkbutton(frame, text=" Logaritmic z-axis") 
            self.chk_removeZonal = ttk.Checkbutton(frame, text=" Remove the zonal modes") 
     
            # Add the commands
            self.chk_interp.config(variable=self.var_interpolate, command=self.update_interp) 
            self.chk_log.config(variable=self.var_log, command=self.update_log) 
            self.chk_removeZonal.config(variable=self.var_removeZonal, command=self.update_removeZonal) 
             
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
            self.chk_removeZonal.grid(row=3, column=0, **PAD_LABEL2, columnspan=2)
            return 

        #-------------------
        def update_log(self, *_): 
            self.PlotNonlinearMap.log = self.var_log.get() 
            self.PlotNonlinearMap.plotOnTheGUI()
            return 
                           
        #-------------------
        def update_interp(self, *_): 
            self.PlotNonlinearMap.interpolate = self.var_interpolate.get() 
            self.PlotNonlinearMap.plotOnTheGUI()
            return 
                            
        #-------------------  
        def update_step(self, *_):
            self.PlotNonlinearMap.step = int(self.var_step.get()) 
            self.PlotNonlinearMap.plotOnTheGUI()
            return 
        
        #-------------------
        def update_removeZonal(self, *_): 
            self.PlotNonlinearMap.removeZonal = self.var_removeZonal.get()         
            self.PlotNonlinearMap.Axis.update_labelFrame(self.PlotNonlinearMap.Canvas.frame)
            self.PlotNonlinearMap.plotOnTheGUI()
            return 
             
 
 
  