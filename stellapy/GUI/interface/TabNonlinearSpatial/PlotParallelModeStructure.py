

import tkinter as tk
from tkinter import ttk   
from stellapy.GUI.widgets import Canvas, PAD_TITLE2, PAD_LABEL2 
from stellapy.GUI.plot.nonlinearSimulations.plot_quantityVsZ import plot_quantityVsZ 

################################################################################
#          CLASS TO PLOT THE PARALLEL MODE STRUCTURE OF THE MODES
################################################################################
   
class PlotParallelModeStructure: 
    ''' Holds and manipulates the Canvas and Axis object for the following plots:
            (1) ... '''

    def __init__(self, tab, notebook_tab): 

        # Set an identifier for this plotting class, used as the figure name
        self.identifier = "TabNonlinearSpatial:PlotParallelModeStructure"
        
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
            self.Axis.update_gridSpeficications(top=0.94, left=0.12, right=0.95, bottom=0.15)
            
            # Set the (x,y) quantities and titles for the canvas
            for i in range(len(self.tab.Simulations.options_yquantityKeys)): 
                self.Axis.set_axisQuantities(i, "z", self.tab.Simulations.options_yquantityKeys[i], name=self.tab.Simulations.options_yquantityNames[i])

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
                "specie"        : self.tab.Simulations.species[0],\
            # Speficy which data to plot 
                "x_quantity"    : Canvas.Axis.x_quantity,\
                "y_quantity"    : Canvas.Axis.y_quantity,\
            # Specify data range
                "x_range"       : Canvas.Axis.x_range,\
                "y_range"       : Canvas.Axis.y_range,\
                "units"         : Canvas.Axis.units,\
            # Eliminate time dimension
                "percentage"    : 50,\
                "time"          : "averaged",\
            # Labels
                "x_label"       : Canvas.Axis.x_label,\
                "y_label"       : Canvas.Axis.y_label,\
                "title"         : Canvas.Axis.title,\
            # Options 
                "rescaleToOne"  : self.normalizeToOne,\
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
                plot_quantityVsZ(**plotting_arguments)
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
        self.splitZonal = False
        self.show_zFrame = True   
        self.normalizeToOne = False 
       
        # Create a subframe for each set of plotting options
        self.subframe_xQuantity = ttk.Frame(frame)
        self.subframe_yQuantity = ttk.Frame(frame)
        self.subframe_options = ttk.Frame(frame)
         
        # Configure the main frame to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1,  uniform="yes")
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=3,  uniform="yes")
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_xQuantity
        tk.Grid.columnconfigure(frame, 1, weight=0) # subframe_yQuantity
        tk.Grid.columnconfigure(frame, 2, weight=0) # subframe_options 
 
        # Add the subframes to the main frame
        self.subframe_xQuantity.grid(  row=1, column=0, padx=25, pady=(0,0), stick='NSEW') 
        self.subframe_yQuantity.grid(  row=1, column=1, padx=25, pady=(0,0), stick='NSEW') 
        self.subframe_options.grid(    row=1, column=2, padx=25, pady=(0,0), stick='NSEW')
         
        # Fill the frames with their options
        self.XQuantity = self.Frame_XQuantity(self, self.subframe_xQuantity)
        self.YQuantity = self.Frame_YQuantity(self, self.subframe_yQuantity) 
        self.Options   = self.Frame_Options(self, self.subframe_options)
        return    

    #----------------------  
    class Frame_XQuantity:
         
        def __init__(self, PlotParallelModeStructure, frame):
             
            # Make the parents available
            self.PlotParallelModeStructure = PlotParallelModeStructure  
         
            # Configure the frame <self.subframe_xQuantity>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Variables
            self.var_xQuantity = tk.IntVar(value=1)
             
            # Create the widgets
            self.lbl_xQuantity  = ttk.Label(frame, text="Quantity on the x-axis", style='prefTitle.TLabel')
            self.rbn_z          = ttk.Radiobutton(frame, text='  z/a [-pi, pi]')
            self.rbn_zeta       = ttk.Radiobutton(frame, text='  Zeta')
            self.rbn_pol        = ttk.Radiobutton(frame, text='  Number of poloidal turns')
            self.rbn_tor        = ttk.Radiobutton(frame, text='  Number of toroidal turns')
     
            # Add the values and commands to the radiobuttons
            self.rbn_z.config(   value=1, variable=self.var_xQuantity, command=self.change_xQuantity)
            self.rbn_zeta.config(value=2, variable=self.var_xQuantity, command=self.change_xQuantity)
            self.rbn_pol.config( value=3, variable=self.var_xQuantity, command=self.change_xQuantity)
            self.rbn_tor.config( value=4, variable=self.var_xQuantity, command=self.change_xQuantity)
             
            # Add the options to the frame
            self.lbl_xQuantity.grid(row=1, column=0, **PAD_TITLE2)
            self.rbn_z.grid(        row=2, column=0, **PAD_LABEL2)
            self.rbn_zeta.grid(     row=3, column=0, **PAD_LABEL2)
            self.rbn_pol.grid(      row=4, column=0, **PAD_LABEL2)
            self.rbn_tor.grid(      row=5, column=0, **PAD_LABEL2)
            return 
         
        #---------------------- 
        def change_xQuantity(self, *_):
             
            # Extract what needs to be plotted on the y-axis 
            if self.var_xQuantity.get()==1: x_quantity = "z"  
            if self.var_xQuantity.get()==2: x_quantity = "zeta" 
            if self.var_xQuantity.get()==3: x_quantity = "pol"
            if self.var_xQuantity.get()==4: x_quantity = "tor"

            # Set the (x,y) quantities for the axis and the canvas
            plot_ids, axis_quantities = self.PlotParallelModeStructure.Axis.get_axisQuantities()
            for plot_id in plot_ids: self.PlotParallelModeStructure.Axis.set_axisQuantities(plot_id, x_quantity, axis_quantities[plot_id]['y']) 
            self.PlotParallelModeStructure.Axis.update_axisQuantities()
            
            # Reset the axis and replot the graph
            self.PlotParallelModeStructure.Axis.reset_axis() 
            self.PlotParallelModeStructure.plotOnTheGUI()
            return

    #----------------------  
    class Frame_YQuantity:
         
        def __init__(self, PlotParallelModeStructure, frame):
             
            # Make the parents available
            self.PlotParallelModeStructure = PlotParallelModeStructure  
         
            # Configure the frame <self.subframe_yQuantity>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Variables
            self.var_yQuantity = tk.IntVar(value=1)
             
            # Create the widgets
            self.lbl_yQuantity  = ttk.Label(frame, text="Quantity on the y-axis", style='prefTitle.TLabel')
            self.rbn_phi2       = ttk.Radiobutton(frame, text='  Potential squared')
            self.rbn_phi2_nozonal = ttk.Radiobutton(frame, text='  Phi² no zonal')
            self.rbn_phi2_zonal = ttk.Radiobutton(frame, text='  Phi² only zonal')
            self.rbn_phi_real   = ttk.Radiobutton(frame, text='  Real(phi)')
     
            # Add the values and commands to the radiobuttons
            self.rbn_phi2.config(        value=1, variable=self.var_yQuantity, command=self.change_yQuantity)
            self.rbn_phi2_nozonal.config(value=2, variable=self.var_yQuantity, command=self.change_yQuantity)
            self.rbn_phi2_zonal.config(  value=3, variable=self.var_yQuantity, command=self.change_yQuantity)
            self.rbn_phi_real.config(    value=4, variable=self.var_yQuantity, command=self.change_yQuantity)
             
            # Add the options to the frame
            self.lbl_yQuantity.grid(   row=1, column=0, **PAD_TITLE2)
            self.rbn_phi2.grid(        row=2, column=0, **PAD_LABEL2)
            self.rbn_phi2_nozonal.grid(row=3, column=0, **PAD_LABEL2)
            self.rbn_phi2_zonal.grid(  row=4, column=0, **PAD_LABEL2)
            self.rbn_phi_real.grid(    row=5, column=0, **PAD_LABEL2)
            return 
         
        #---------------------- 
        def change_yQuantity(self, *_):
             
            # Extract what needs to be plotted on the y-axis 
            if self.var_yQuantity.get()==1: y_quantity = "phi2"  
            if self.var_yQuantity.get()==2: y_quantity = "phi2_nozonal" 
            if self.var_yQuantity.get()==3: y_quantity = "phi2_zonal"
            if self.var_yQuantity.get()==4: y_quantity = "phi_real"

            # Set the (x,y) quantities for the axis and the canvas
            plot_ids, axis_quantities = self.PlotParallelModeStructure.Axis.get_axisQuantities()
            for plot_id in plot_ids: self.PlotParallelModeStructure.Axis.set_axisQuantities(plot_id, axis_quantities[plot_id]['x'], y_quantity) 
            self.PlotParallelModeStructure.Axis.update_axisQuantities()
            
            # Reset the axis and replot the graph
            self.PlotParallelModeStructure.Axis.reset_axis() 
            self.PlotParallelModeStructure.plotOnTheGUI()
            return
        
    #---------------------- 
    class Frame_Options:
         
        def __init__(self, PlotParallelModeStructure, frame):
 
            # Make the parents available
            self.PlotParallelModeStructure = PlotParallelModeStructure  
            
            # Variables
            self.var_splitZonal = tk.IntVar(value=self.PlotParallelModeStructure.splitZonal) 
            self.var_zFrame = tk.IntVar(value=self.PlotParallelModeStructure.show_zFrame) 
            self.var_normalizeToOne = tk.IntVar(value=self.PlotParallelModeStructure.normalizeToOne)
               
            # Create the widgets
            self.lbl_options = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_splitZonal = ttk.Checkbutton(frame, text=" Split in zonal contributions")  
            self.chk_normalizeToOne = ttk.Checkbutton(frame, text=" Normalize to one")  
            self.chk_vspan = ttk.Checkbutton(frame, text=" Show z frames")  
             
            # Add the commands
            self.chk_splitZonal.config( variable=self.var_splitZonal,  command=self.update_splitZonal) 
            self.chk_normalizeToOne.config(variable=self.var_normalizeToOne, command=self.update_normalizeToOne) 
            self.chk_vspan.config(variable=self.var_zFrame, command=self.update_vspan) 
             
            # Configure the frame
            tk.Grid.columnconfigure(frame, 0, weight=1) 
             
            # Add the options to the frame
            self.lbl_options.grid(       row=0, column=0, **PAD_TITLE2)
            self.chk_splitZonal.grid(    row=1, column=0, **PAD_LABEL2)
            self.chk_normalizeToOne.grid(row=2, column=0, **PAD_LABEL2)
            self.chk_vspan.grid(         row=3, column=0, **PAD_LABEL2) 
            return 
         
        #---------------------- 
        def update_vspan(self, *_): 
            self.PlotParallelModeStructure.show_zFrame = self.var_zFrame.get() 
            self.PlotParallelModeStructure.plotOnTheGUI()
 
        #---------------------- 
        def update_normalizeToOne(self, *_): 
            self.PlotParallelModeStructure.normalizeToOne = self.var_normalizeToOne.get() 
            self.PlotParallelModeStructure.plotOnTheGUI()
            return 
        
        #---------------------- 
        def update_splitZonal(self, *_): 
            self.PlotParallelModeStructure.splitZonal = self.var_splitZonal.get() 
            self.PlotParallelModeStructure.plotOnTheGUI()
            return 
            