
import tkinter as tk
from tkinter import ttk 
from stellapy.GUI.interface.TabConvergence.Simulations import Simulations
from stellapy.GUI.interface.TabConvergence.PlotTimeEvolution import PlotTimeEvolution 
from stellapy.GUI.interface.TabConvergence.PlotParallelModeStructure import PlotParallelModeStructure

################################################################################
#     CLASS TO STUDY THE TIME AND SPACE CONVERGENCE OF LINEAR SIMULATIONS
################################################################################

class TabConvergence:
    ''' Initiate the frames on the tab "Convergence" to investigate the 
    convergence in time and space of linear modes simulated with stella.
    The tab can plot the following graphs: 
            (1) Omega versus time
            (2) Gamma versus time
            (3) Phi**2 versus time
            (4) Re(Phi) versus zeta
            (5) Im(Phi) versus zeta
            (6) Phi**2 versus zeta  ''' 

    def __init__(self, tab): 

        # Make the parent objects available      
        self.dict_awthemes = tab.root.awthemes             
        self.root = tab.root               

        # Create the frames for the two canvasses and options on <tabConvergence>
        self.frame_canvas1 = ttk.LabelFrame(tab, text="   Time evolution of the modes  ",          **self.dict_awthemes['labelframe2'])
        self.frame_canvas2 = ttk.LabelFrame(tab, text="   Parallel mode structure of the modes  ", **self.dict_awthemes['labelframe2'])
        self.frame_options = ttk.Frame(tab)    
        
        # Configure the rows/columns of the grid on <tabConvergence>
        tk.Grid.rowconfigure(   tab, 0, weight=1) 
        tk.Grid.columnconfigure(tab, 0, weight=10, uniform="columns") 
        tk.Grid.columnconfigure(tab, 1, weight=3,  uniform="columns")
            
        # Add the two canvasses and options frames to <tabConvergence>
        self.frame_canvas1.grid(in_=tab, row=0, column=0, padx=(20,5), pady=(20,20), stick='NSEW')
        self.frame_canvas2.grid(in_=tab, row=0, column=0, padx=(20,5), pady=(20,20), stick='NSEW')
        self.frame_options.grid(in_=tab, row=0, column=1, padx=(5,20), pady=(20,20), stick='NSEW')        
        self.frame_canvas2.grid_remove()      
         
        # Bind focus on every click on a widget: this unfoces the entry widgets 
        self.frame_canvas1.bind("<1>", lambda _: self.frame_canvas1.focus_set())
        self.frame_canvas2.bind("<1>", lambda _: self.frame_canvas2.focus_set())
        self.frame_options.bind("<1>", lambda _: self.frame_options.focus_set())  
        
        # Fill the subframes with widgets
        self.Simulations = Simulations(self, self.frame_options)
        self.PlotTimeEvolution = PlotTimeEvolution(self)
        self.PlotParallelModeStructure = PlotParallelModeStructure(self)
        
        #========================
        # LOAD TAB AND CANVASSES
        #========================
         
        # When the tab is visible, make sure the tab loaded the correct widgets, 
        # and that the canvasses are initiated.
        self.frame_options.bind("<Visibility>", self.load_tab)
        self.frame_canvas1.bind("<Visibility>", self.load_canvas1)
        self.frame_canvas2.bind("<Visibility>", self.load_canvas2)
        if True: return
 
################################################################################
#                                     METHODS
################################################################################

    def plot(self):
        ''' From <Simulations> we can ask to plot the figure. ''' 
        if self.Simulations.var_plot.get() in [1,2,3]:   
            self.PlotTimeEvolution.plotOnTheGUI()
        if self.Simulations.var_plot.get() in [4,5,6]:  
            self.PlotParallelModeStructure.plotOnTheGUI()
        return
    
    #-------------------------------------  
    def reset_axes(self):
        ''' From <Simulations> we can ask to reset the figures. '''
        if self.PlotTimeEvolution.initiated_canvas: self.PlotTimeEvolution.Canvas.Axis.reset_axis()
        if self.PlotParallelModeStructure.initiated_canvas: self.PlotParallelModeStructure.Canvas.Axis.reset_axis()

################################################################################
#                          LOAD TAB AND CANVASSES
################################################################################
    
    def load_tab(self, *_):
        '''  When the tab is visible, make sure the dropdown menus for the 
        experiment and simulations are updated. Connect the progress widget. '''
        
        # Connect the Progress widget of this tab to the research object
        self.root.Progress = self.Progress
           
        # Make sure the drpdown menus for the experiments/simulations are updated
        if self.root.Research.input_files != []: 
            self.Simulations.display_plottedModesAndExperiments()
        return
    
    #-------------------------------------  
    def load_canvas1(self, *_):
        '''  When the frame of the canvas is visible, make sure the canvas is
        initiated, so that we have the figure and axis object. ''' 
        self.PlotTimeEvolution.initiate_canvas(self.frame_canvas1)
        return
 
    #-------------------------------------  
    def load_canvas2(self, *_):
        '''  When the frame of the canvas is visible, make sure the canvas is
        initiated, so that we have the figure and axis object. ''' 
        self.PlotParallelModeStructure.initiate_canvas(self.frame_canvas2)
        return 
    
    