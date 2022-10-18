
import tkinter as tk
from tkinter import ttk  
from stellapy.GUI.widgets import Canvas, PAD_TITLE2, PAD_LABEL2
from stellapy.GUI.plot.nonlinearSimulations.plot_potentialVsKxVsKy import plot_potentialVsKxVsKy

#################################################################
#                 CLASS TO PLOT THE LINEAR MAP
#################################################################

class PlotSpectrumKxKy: 
    ''' Holds and manipulates the Canvas and Axis object for the following plots:
            (1) ... '''

    def __init__(self, tab, notebook_tab): 

        # Set an identifier for this plotting class, used as the figure name
        self.identifier = "TabNonlinearSpatial:PlotSpectrumKxKy"
        
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
            self.Axis.update_gridSpeficications(top=0.94, left=0.1, right=0.85, bottom=0.15)
            
            # Set the (x,y) quantities and titles for the canvas
            for i in range(len(self.tab.Simulations.options_yquantityKeys)): 
                self.Axis.set_axisQuantities(i, "kx", "ky",  self.tab.Simulations.options_yquantityKeys[i])

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
                "z_specific"    : self.tab.DataSelection.z_specific,\
                "t_specific"    : self.tab.DataSelection.t_specific,\
                "x_quantity"    : Canvas.Axis.x_quantity,\
                "y_quantity"    : Canvas.Axis.y_quantity,\
                "z_quantity"    : Canvas.Axis.z_quantity,\
            # Specify data range
                "x_range"       : Canvas.Axis.x_range,\
                "y_range"       : Canvas.Axis.y_range,\
                "c_range"       : None,\
                "units"         : Canvas.Axis.units,\
            # Labels
                "x_label"       : Canvas.Axis.x_label,\
                "y_label"       : Canvas.Axis.y_label,\
                "title"         : Canvas.Axis.title,\
            # Options
                "log"           : self.log,\
                "interpolate"   : self.interpolate,\
                "step"          : self.step,\
                "ordersOfMagnitude" : self.ordersOfMagnitude,\
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

                # Remove the colorbar
                try: self.colorbar.remove()
                except: pass
                
                # Plot the parallel mode structure
                Canvas.Axis.start_plotting(research_arguments, input_files)
                self.colorbar = plot_potentialVsKxVsKy(**plotting_arguments)
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
        self.interpolate = False
        self.log = False
        self.step  = 20
        self.colorbar = None
        self.ordersOfMagnitude = None
         
        # Create a subframe for each set of plotting options  
        self.subframe_xyQuantity = ttk.Frame(frame) 
        self.subframe_buttons    = ttk.Frame(frame) 
        self.subframe_options = ttk.Frame(frame)
         
        # Configure the main frame to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1,  uniform="yes")
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=0)      
        tk.Grid.rowconfigure(   frame, 3, weight=2,  uniform="yes")  
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_xyQuantity   
        tk.Grid.columnconfigure(frame, 1, weight=0) # subframe_options
 
        # Add the subframes to the main frame 
        self.subframe_xyQuantity.grid(  row=1, column=0, padx=25, pady=(0,0), stick='NSEW') 
        self.subframe_buttons.grid(     row=2, column=0, padx=25, pady=(0,0), stick='NSEW')   
        self.subframe_options.grid(     row=1, column=1, padx=25, pady=(0,0), stick='NSEW', rowspan=2)
        
        # Fill the frames with their options 
        self.XYQuantity = self.Frame_XYQuantity(self, self.subframe_xyQuantity) 
        self.Buttons    = self.Frame_Buttons(   self, self.subframe_buttons) 
        self.Options    = self.Frame_Options(   self, self.subframe_options)
        return     

    #----------------------  
    class Frame_XYQuantity:
         
        def __init__(self, PlotSpectrumKxKy, frame):
             
            # Make the parents available
            self.PlotSpectrumKxKy = PlotSpectrumKxKy 
            self.tab = PlotSpectrumKxKy.tab
         
            # Configure the frame <self.subframe_xyQuantity>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Variables
            self.var_xyQuantity = tk.IntVar(value=2)
             
            # Create the widgets
            self.lbl_xyQuantity = ttk.Label(frame, text="Quantity on the xy-axis", style='prefTitle.TLabel')
            self.rbn_xy         = ttk.Radiobutton(frame, text='  Real space (x,y)')
            self.rbn_kxky       = ttk.Radiobutton(frame, text='  Reciprocal (kx,ky)') 
     
            # Add the values and commands to the radiobuttons
            self.rbn_xy.config(  value=1, variable=self.var_xyQuantity, command=self.change_xyQuantity)
            self.rbn_kxky.config(value=2, variable=self.var_xyQuantity, command=self.change_xyQuantity) 
             
            # Add the options to the frame
            self.lbl_xyQuantity.grid(row=1, column=0, **PAD_TITLE2)
            self.rbn_xy.grid(        row=2, column=0, **PAD_LABEL2)
            self.rbn_kxky.grid(      row=3, column=0, **PAD_LABEL2) 
            return 
         
        #---------------------
        def change_xyQuantity(self, *_):
             
            # Extract what needs to be plotted on the xy-axis 
            if self.var_xyQuantity.get()==1: x_quantity = "x";  y_quantity = "y"
            if self.var_xyQuantity.get()==2: x_quantity = "kx"; y_quantity = "ky"
             
            # Set the (x,y,z) quantities for the axis and the canvas
            plot_ids, axis_quantities = self.PlotSpectrumKxKy.Axis.get_axisQuantities()
            for plot_id in plot_ids: self.PlotSpectrumKxKy.Axis.set_axisQuantities(plot_id, x_quantity, y_quantity, axis_quantities[plot_id]['z']) 
            self.PlotSpectrumKxKy.Axis.update_axisQuantities()
            
            # Plot the graph
            self.PlotSpectrumKxKy.Axis.reset_axis()
            self.PlotSpectrumKxKy.Axis.update_labelFrame(self.PlotSpectrumKxKy.Canvas.frame)
            self.PlotSpectrumKxKy.plotOnTheGUI() 
            return

    #----------------------  
    class Frame_Buttons:
         
        def __init__(self, PlotSpectrumKxKy, frame):
             
            # Make the parents available
            self.PlotSpectrumKxKy = PlotSpectrumKxKy 
            self.tab = PlotSpectrumKxKy.tab
            self.root = PlotSpectrumKxKy.tab.root
         
            # Configure the frame <self.subframe_buttons>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            tk.Grid.columnconfigure(frame, 1, weight=0) 
            tk.Grid.columnconfigure(frame, 2, weight=0) 
            tk.Grid.columnconfigure(frame, 0, weight=1) 
             
            # Create the widgets
            self.btn_next = ttk.Button(frame, text="Next >", width=10)
            self.btn_next.config(command=lambda: self.next_simulation())
            self.btn_prev = ttk.Button(frame, text="< Previous", width=10)
            self.btn_prev.config(command=lambda: self.prev_simulation())

            # Add the options to the frame
            self.btn_prev.grid( row=0, column=1, padx=(0,2), pady=(2,2), ipadx=1, ipady=1) 
            self.btn_next.grid( row=0, column=2, padx=(2,0), pady=(2,2), ipadx=2, ipady=1)
         
        #---------------------
        def next_simulation(self, *_):

            # Select the next simulation
            simulations = self.tab.Simulations.options_simulations
            simulation_id = self.tab.Simulations.var_simulation.get()
            simulation_index = simulations.index(simulation_id)
            if simulation_index == len(simulations)-1: simulation_index=-1
            self.tab.Simulations.var_simulation.set(simulations[simulation_index+1])

            # Since the simulation changes, perhaps the z where |phi**2|(z) is maximum changed
            if self.tab.DataSelection.ZSelection.var_z.get()==2:
                self.tab.DataSelection.ZSelection.calculate_zIndexWherePhi2IsMaximum()
                
            # Plot the graph
            self.tab.reset_Axis(); self.tab.plot()
            return
        
        #---------------------
        def prev_simulation(self, *_):
             
            # When changing the plots, reset the axis of the graph
            self.tab.Graph[self.plot.graph_id].load_defaults()
            self.tab.Graph[self.plot.graph_id].update_quantitiesAndKeys()

            # Select the previous simulation
            simulations = self.tab.Simulations.options_simulations
            simulation_id = self.tab.Simulations.var_simulation.get()
            simulation_index = simulations.index(simulation_id)
            if simulation_index == 0: len(simulations)
            self.tab.Simulations.var_simulation.set(simulations[simulation_index-1])

            # Since the simulation changes, perhaps the z where |phi**2|(z) is maximum changed
            if self.tab.DataSelection.ZSelection.var_z.get()==2:
                self.tab.DataSelection.ZSelection.calculate_zIndexWherePhi2IsMaximum()
                
            # Plot the graph
            self.tab.reset_Axis(); self.tab.plot()
            return
        
    #----------------------  
    class Frame_Options:
         
        def __init__(self, PlotSpectrumKxKy, frame):
         
            # Make the parents available
            self.PlotSpectrumKxKy = PlotSpectrumKxKy 
            self.tab = PlotSpectrumKxKy.tab
            self.root = PlotSpectrumKxKy.tab.root
            
            # Variables
            self.var_interpolate = tk.IntVar(value=self.PlotSpectrumKxKy.interpolate)
            self.var_log = tk.IntVar(value=self.PlotSpectrumKxKy.log)
            self.var_orders = tk.IntVar(value=self.PlotSpectrumKxKy.ordersOfMagnitude)
            self.var_step = tk.StringVar(value="20")
            self.var_ordersOfMagnitude = tk.StringVar(value="2")
        
            # Create the widgets
            self.lbl_options = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_interp = ttk.Checkbutton(frame, text=" Interpolate with step = ") 
            self.chk_log = ttk.Checkbutton(frame, text=" Logarithmic z-axis") 
            self.chk_orders= ttk.Checkbutton(frame, text=" Orders of magnitude = ") 
     
            # Add the commands
            self.chk_interp.config(variable=self.var_interpolate, command=self.update_interp)
            self.chk_log.config(variable=self.var_log, command=self.update_log)  
            self.chk_orders.config(variable=self.var_orders, command=self.update_orders)  
             
            # Entry widget to choose the interpolation step 
            self.ent_interp = ttk.Entry(frame, textvariable=self.var_step, font=("Courier New", 11), style='opt_valueCBold.TEntry', width=5)
            self.ent_interp.bind('<Return>', self.update_step) 
            self.ent_orders = ttk.Entry(frame, textvariable=self.var_ordersOfMagnitude, font=("Courier New", 11), style='opt_valueCBold.TEntry', width=5)
            self.ent_orders.bind('<Return>', self.update_ordersOfMagnitude) 
             
            # Configure the frame
            tk.Grid.rowconfigure(frame, 0, weight=0) 
            tk.Grid.rowconfigure(frame, 1, weight=0) 
            tk.Grid.columnconfigure(frame, 0, weight=0) 
            tk.Grid.columnconfigure(frame, 1, weight=0) 
            tk.Grid.columnconfigure(frame, 2, weight=1) 
             
            # Add the options to the frame
            self.lbl_options.grid(row=0, column=0, **PAD_TITLE2, columnspan=2)
            self.chk_interp.grid( row=1, column=0, **PAD_LABEL2)
            self.ent_interp.grid( row=1, column=1) 
            self.chk_log.grid(    row=2, column=0, **PAD_LABEL2, columnspan=2)
            self.chk_orders.grid( row=3, column=0, **PAD_LABEL2)
            self.ent_orders.grid( row=3, column=1) 
            return 
                            
        #-------------------
        def update_interp(self, *_): 
            self.PlotSpectrumKxKy.interpolate = self.var_interpolate.get() 
            self.PlotSpectrumKxKy.Axis.reset_axis() 
            self.PlotSpectrumKxKy.plotOnTheGUI()
            return 

        #-------------------
        def update_log(self, *_): 
            self.PlotSpectrumKxKy.log = self.var_log.get() 
            self.PlotSpectrumKxKy.plotOnTheGUI()
            return 
                            
        #-------------------  
        def update_step(self, *_):
            self.PlotSpectrumKxKy.step = int(self.var_step.get()) 
            self.PlotSpectrumKxKy.plotOnTheGUI()
            return 
    
        #-------------------
        def update_orders(self, *_):  
            if self.var_interpolate.get()==True:  self.PlotSpectrumKxKy.ordersOfMagnitude = float(self.var_ordersOfMagnitude.get()) 
            if self.var_interpolate.get()==False: self.PlotSpectrumKxKy.ordersOfMagnitude = None 
            self.PlotSpectrumKxKy.Axis.reset_axis() 
            self.PlotSpectrumKxKy.plotOnTheGUI()
            return 
 
        #-------------------  
        def update_ordersOfMagnitude(self, *_):
            self.PlotSpectrumKxKy.ordersOfMagnitude = float(self.var_ordersOfMagnitude.get()) 
            self.PlotSpectrumKxKy.plotOnTheGUI()
            return 
 
  