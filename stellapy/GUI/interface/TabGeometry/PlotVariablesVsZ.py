
# Load modules
import tkinter as tk
from tkinter import ttk   
from stellapy.GUI.widgets import PAD_TITLE2, PAD_LABEL2  
from stellapy.plot.geometry.plot_geometryVsZ import plot_geometryVsZ 

#######################################################################
#   CLASS TO PLOT THE INFLUENCE OF A PARAMETER ON THE SATURATED FLUX
#######################################################################
   
class PlotVariablesVsZ: 
    
    def __init__(self, parent, frame, graph_id):
        
        # Make the parents available
        self.tab = parent               # The tab "Geometry"
        self.root = parent.root         # The root of the GUI
        self.frame = frame              # The frame "frame_graph0" 
        self.graph_id = graph_id        # The id for self.Graph and self.Canvas
        
        # Keep track of what is plotted
        self.replot = False
        self.plotted_inputFiles = None
        
        # Data to plot
        self.x_quantity = "B_zeta"  
        
        # Toggles
        self.showOnlyMarkers = False 
        self.normalize = False
       
        # Create a subframe for each set of plotting options
        self.subframe_xQuantity = ttk.Frame(frame) 
        self.subframe_options = ttk.Frame(frame)
         
        # Configure the main frame to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1,  uniform="yes")
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=3,  uniform="yes")
        tk.Grid.columnconfigure(frame, 0, weight=0) 
        tk.Grid.columnconfigure(frame, 1, weight=0) 
 
        # Add the subframes to the main frame
        self.subframe_xQuantity.grid(  row=1, column=0, padx=25, pady=(0,0), stick='NSEW') 
        self.subframe_options.grid(    row=1, column=1, padx=25, pady=(0,0), stick='NSEW')
         
        # Fill the frames with their options
        self.XQuantity  = self.Frame_XQuantity(self, self.subframe_xQuantity) 
        self.Options    = self.Frame_Options(self, self.subframe_options)
         
        # When the tab is visible, make sure its canvas and frame are loaded
        self.frame.bind("<Visibility>", self.load_figure) 
        return     
    
######################################## 
#        LOAD AND PLOT FIGURE
########################################

    def load_figure(self, event):
        
        # Bind <ctrl+s> to replot the graphs
        self.root.bind('<Control-s>', self.replot_graph)
        return 
        
    #----------------------
    def replot_graph(self):
        self.replot = True
        self.plot_graph(None)
        return
            
    #-------------------------------------
    def plot_graph(self, poppedout_id=None):
         
        # Perhaps the button was clicked on the popped out window
        if poppedout_id == None: Graph = self.tab.Graph[self.graph_id]
        if poppedout_id != None: Graph = self.root.graph_poppedOut[poppedout_id]
         
        # Only proceed if there are input_files
        if self.input_files != []:
             
            # Get the input_files that need to be plotted
            input_files = self.root.Research.input_files
             
            # If the selection of input_files changed, then reset the axis
            if (self.plotted_inputFiles != input_files):
                Graph.load_defaults()
                 
            # If there are files plot them, if they changed plot them, if its a popout plot them
            if (input_files != [] and self.plotted_inputFiles != input_files) or poppedout_id!=None or self.replot==True:
                 
                # Clear the axis and remember which input_files are plotted
                Graph.ax.clear()
                self.plotted_inputFiles = input_files.copy() 
                self.replot = False
                 
                # While plotting show a progress bar turn the cursor into a waiting cursor
                self.tab.show_progress("start", poppedout_id)   
                    
                # Plot the graph   
                plot_geometryVsZ(\
                    # Specify which simulations to plot
                        research=self.root.Research,\
                        experiment_id=self.tab.Simulations.experiment_id,\
                        simulation_id=self.tab.Simulations.simulation_id,\
                    # Specify data range 
                        x_quantity=self.x_quantity,\
                        y_quantity=self.tab.Simulations.y_quantity,\
                        x_range=Graph.range["x"],\
                        y_range=Graph.range["y"],\
                    # Labels
                        x_label=Graph.label["x"],\
                        y_label=Graph.label["y"],\
                        title=Graph.label["title"],\
                    # For the GUI the figure object already exists 
                        show_figure = False,\
                        ax=Graph.ax, \
                        Progress=self.tab.Progress,\
                    # Toggles
                        tick_style='sci',\
                        normalize=self.normalize,\
                        showOnlyMarkers=self.showOnlyMarkers,\
                        log=False,\
                        units="normalized") 
                     
                # Update the <graph> class so the option window can use them
                Graph.update_rangesAndLabels()
             
                # When finished show the canvas and change the cursor back.
                self.tab.show_progress("finished", poppedout_id) 
     
        # If no simulations have been selected, clear the current figure
        if self.input_files == [] or self.root.Research.input_files == []:
            Graph.ax.clear()
            Graph.load_defaults()
            self.tab.show_progress("nothing", poppedout_id) 
         
        # Update screen
        if poppedout_id==None: self.tab.Canvas[self.graph_id].draw_idle()
        if poppedout_id!=None: self.root.canvasPoppedOut[poppedout_id].draw_idle()
        if True: return

######################################## 
#        CLASSES FOR THE OPTIONS
########################################

    class Frame_XQuantity:
         
        def __init__(self, parent, frame):
             
            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
         
            # Configure the frame <self.subframe_xQuantity>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Variables
            self.var_xQuantity = tk.IntVar(value=2)
             
            # Create the widgets
            self.lbl_xQuantity  = ttk.Label(frame, text="Quantity on the X-axis", style='prefTitle.TLabel')
            self.rbn_z          = ttk.Radiobutton(frame, text='  z/a [-pi, pi]')
            self.rbn_zeta       = ttk.Radiobutton(frame, text='  Zeta') 
     
            # Add the values and commands to the radiobuttons
            self.rbn_z.config(   value=1, variable=self.var_xQuantity, command=self.change_xQuantity)
            self.rbn_zeta.config(value=2, variable=self.var_xQuantity, command=self.change_xQuantity) 
             
            # Add the options to the frame
            self.lbl_xQuantity.grid(row=1, column=0, **PAD_TITLE2)
            self.rbn_z.grid(        row=2, column=0, **PAD_LABEL2)
            self.rbn_zeta.grid(     row=3, column=0, **PAD_LABEL2) 
            return 
         
        #---------------------
        def change_xQuantity(self, *args):
             
            # When changing the plots, reset the axis of the graph
            self.tab.Graph[self.plot.graph_id].load_defaults()
            self.tab.Graph[self.plot.graph_id].update_quantitiesAndKeys()
             
            # Extract what needs to be plotted on the y-axis
            x_quantity = self.var_xQuantity.get()
            if x_quantity==1: self.plot.x_quantity = "B_zed";    
            if x_quantity==2: self.plot.x_quantity = "B_zeta";    
             
            # Change the header of the frame and plot the graph 
            self.plot.replot_graph()
            return

    class Frame_Options:
         
        def __init__(self, parent, frame):
 
            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
            
            # Variables
            self.var_showMarkers = tk.IntVar(value=0) 
            self.var_normalize = tk.IntVar(value=0)  
        
            # Create the widgets
            self.lbl_options = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_showMarker   = ttk.Checkbutton(frame, text=" Only show markers")  
            self.chk_normalize   = ttk.Checkbutton(frame, text=" Normalize data")  
             
            # Add the commands
            self.chk_showMarker.config(variable=self.var_showMarkers, command=self.update_showMarker) 
            self.chk_normalize.config(variable=self.var_normalize, command=self.update_normalize) 
             
            # Configure the frame
            tk.Grid.columnconfigure(frame, 0, weight=1) 
             
            # Add the options to the frame
            self.lbl_options.grid(    row=0, column=0, **PAD_TITLE2)
            self.chk_showMarker.grid( row=1, column=0, **PAD_LABEL2) 
            self.chk_normalize.grid(  row=2, column=0, **PAD_LABEL2) 
            return 
         
        #-----------------------
        def update_showMarker(self, *args): 
            if self.var_showMarkers.get()==0: self.plot.showOnlyMarkers = False
            if self.var_showMarkers.get()==1: self.plot.showOnlyMarkers = True
            self.plot.replot_graph()
 
        #-----------------------
        def update_normalize(self, *args): 
            if self.var_normalize.get()==0: self.plot.normalize = False
            if self.var_normalize.get()==1: self.plot.normalize = True
            self.plot.replot_graph()
            