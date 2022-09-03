 
# TODO: Change "impurity x" to something that recognizes the impurity
 
#################################################################
#                   CLASS FOR THE FIRST TAB
#################################################################
''' 
'''
 
# Load modules
import tkinter as tk
from tkinter import ttk 
from stellapy.GUI.widgets import Progress 
from stellapy.GUI.widgets.DisplayTimeFrames import DisplayTimeFrames
from stellapy.GUI.interface.TabNonlinearSpatial.Simulations import Simulations
from stellapy.GUI.interface.TabNonlinearSpatial.DataSelection import DataSelection
from stellapy.GUI.interface.TabNonlinearSpatial.PlotTimeEvolution import PlotTimeEvolution
from stellapy.GUI.interface.TabNonlinearSpatial.PlotParallelModeStructure import PlotParallelModeStructure    
from stellapy.GUI.interface.TabNonlinearSpatial.PlotSpectrumKxKy import PlotSpectrumKxKy   
from stellapy.GUI.interface.TabNonlinearSpatial.Plot1DSpectrum import Plot1DSpectrum    
 
#################################################################
#                   CLASS FOR THE FIRST TAB
#################################################################
class TabNonlinearSpatial:
    ''' Initiate the frames on the tab "TabNonlinearSpatial" to investigate the results 
    of nonlinear simulations with stella. The tab can plot the following graphs: 
            (1) ...  '''
  
    def __init__(self, tab): 
 
        # Make the parent objects available      
        self.dict_awthemes = tab.root.awthemes             
        self.root = tab.root     
 
        # Create a frame for each graph and one for the options and add them to the tab "Nonlinear time traces" 
        self.frame_canvas1 = ttk.LabelFrame(tab, text="   Time evolution  ", **self.dict_awthemes['labelframe2'])
        self.frame_canvas2 = ttk.LabelFrame(tab, text="   Parallel mode structure  ", **self.dict_awthemes['labelframe2'])
        self.frame_canvas3 = ttk.LabelFrame(tab, text="   Surface plot versus (kx,ky)", **self.dict_awthemes['labelframe2'])
        self.frame_canvas4 = ttk.LabelFrame(tab, text="   1D Spectrum  ", **self.dict_awthemes['labelframe2'])
        self.frame_canvas5 = ttk.LabelFrame(tab, text="   1D Spectrum  ", **self.dict_awthemes['labelframe2'])
        self.frame_options = ttk.Frame(tab)
 
        # Configure the rows/columns of the grid on <TabNonlinearSpatial>
        tk.Grid.rowconfigure(   tab, 0, weight=0) 
        tk.Grid.rowconfigure(   tab, 1, weight=1)  
        tk.Grid.columnconfigure(tab, 0, weight=1,  uniform="for sure") 
        tk.Grid.columnconfigure(tab, 1, weight=1,  uniform="for sure") 
  
        # Add the frames to the main window (= tab "Nonlinear spectra and parallel mode structures" )
        self.frame_options.grid(in_=tab, row=0, column=0, padx=(20,20), pady=(5,10), stick='NSEW', columnspan=3)
        self.frame_canvas1.grid(in_=tab, row=1, column=0, padx=(20,10), pady=(10,20),stick='NSEW')
        self.frame_canvas2.grid(in_=tab, row=1, column=1, padx=(10,10), pady=(10,20),stick='NSEW')
        self.frame_canvas3.grid(in_=tab, row=1, column=1, padx=(10,20), pady=(10,20),stick='NSEW')
        self.frame_canvas4.grid(in_=tab, row=1, column=0, padx=(10,10), pady=(10,20),stick='NSEW')
        self.frame_canvas5.grid(in_=tab, row=1, column=1, padx=(10,20), pady=(10,20),stick='NSEW')
        self.frame_canvas4.grid_remove() # graph1 and graph4 overlap so hide one 
        self.frame_canvas5.grid_remove() # graph2 and graph5 overlap so hide one 
         
        # Bind focus on every click on a widget: this makes sure the entry widgets become unfocused
        self.frame_canvas1.bind("<1>", lambda _: self.frame_canvas1.focus_set())
        self.frame_canvas2.bind("<1>", lambda _: self.frame_canvas2.focus_set())
        self.frame_canvas3.bind("<1>", lambda _: self.frame_canvas3.focus_set()) 
        self.frame_canvas4.bind("<1>", lambda _: self.frame_canvas4.focus_set()) 
        self.frame_canvas5.bind("<1>", lambda _: self.frame_canvas5.focus_set()) 
        self.frame_options.bind("<1>", lambda _: self.frame_options.focus_set()) 

        # Fill the options frame with subframes, and create tabs for the graph options
        self.fill_optionsFrame(self.frame_options)
        
        # When <tabLinear> is visible, make sure the Progress object is linked,
        # that the dropdown menus are updated, and that the canvasses are initiated.
        self.frame_options.bind("<Visibility>", self.load_tab)
        self.frame_canvas1.bind("<Visibility>", self.load_canvas1)
        self.frame_canvas2.bind("<Visibility>", self.load_canvas2)
        self.frame_canvas3.bind("<Visibility>", self.load_canvas3) 
        self.frame_canvas4.bind("<Visibility>", self.load_canvas4)  
        self.frame_canvas5.bind("<Visibility>", self.load_canvas5) 
        
    #-----------------------------
    def fill_optionsFrame(self, frame):
        ''' Create a subframe for the dropdown menus of the experiments and 
        simulations. Create a subframe for the Progress object. Create a tabbed
        widget, with a tab for each graph to hold the various options. '''
     
        # Create the subframes in the options frame
        self.subframe_simulations   = ttk.LabelFrame(frame, text="   Simulations  ", **self.dict_awthemes['labelframe2'])
        self.subframe_dataSelection = ttk.LabelFrame(frame, text="   Data Selection  ", **self.dict_awthemes['labelframe2'])
        self.subframe_tabs          = ttk.Frame(frame)
        self.frame_progress         = ttk.Frame(frame)
 
        # Add the progress object
        self.Progress = Progress(self.root, self.frame_progress, anchor="fill", length=300) 
        self.Progress.move(0,"Please select simulations.")
         
        # Configure the frame
        tk.Grid.rowconfigure(   frame, 0, weight=1)  
        tk.Grid.rowconfigure(   frame, 1, weight=0)     
        tk.Grid.columnconfigure(frame, 0, weight=0)   
        tk.Grid.columnconfigure(frame, 1, weight=0)   
        tk.Grid.columnconfigure(frame, 2, weight=1)    
         
        # Add the labelframes to <frame_options>
        self.subframe_simulations.grid(  row=0, column=0, padx=(0,5), pady=(0,0), stick='NSEW', rowspan=2)
        self.subframe_dataSelection.grid(row=0, column=1, padx=(5,5), pady=(0,0), stick='NSEW')
        self.frame_progress.grid(        row=1, column=1, padx=(5,5), pady=(5,0), stick='NSEW')
        self.subframe_tabs.grid(         row=0, column=2, padx=(5,0), pady=(0,0), stick='NSEW', rowspan=2)
         
        # Fill the subframes with widgets
        self.Simulations = Simulations(self, self.subframe_simulations)
        self.DataSelection = DataSelection(self, self.subframe_dataSelection) 
         
        # Tabbed display for "Plot_linearSpectrum", "Plot_parameterInfluence" and "Plot_linearMap"
        self.notebook = ttk.Notebook(self.subframe_tabs, style='species.TNotebook')  
        self.notebook.grid(row=0, column=0, sticky="NSEW") 
        tk.Grid.rowconfigure(self.subframe_tabs, 0, weight=1)  
        tk.Grid.columnconfigure(self.subframe_tabs, 0, weight=1)  
 
        # Create frames for each tab
        self.notebook_tab1 = ttk.Frame(self.notebook, padding=(0,0,0,0))  # Time evolution
        self.notebook_tab2 = ttk.Frame(self.notebook, padding=(0,0,0,0))  # Parallel mode structure
        self.notebook_tab3 = ttk.Frame(self.notebook, padding=(0,0,0,0))  # Spectrum (kx,ky) or (x,y)
        self.notebook_tab4 = ttk.Frame(self.notebook, padding=(0,0,0,0))  # 1D spectra along kx and ky
        self.notebook.add(self.notebook_tab1, text="Time evolution")  
        self.notebook.add(self.notebook_tab2, text="Parallel mode structure")  
        self.notebook.add(self.notebook_tab3, text="2D Spectrum (kx,ky) or (x,y)")  
        self.notebook.add(self.notebook_tab4, text="1D Spectrum") 
        self.notebook.select(self.notebook_tab1)
                
        # Fill the subframes with widgets
        self.PlotTimeEvolution = PlotTimeEvolution(self, self.notebook_tab1)
        self.PlotParallelModeStructure = PlotParallelModeStructure(self, self.notebook_tab2)
        self.PlotSpectrumKxKy = PlotSpectrumKxKy(self, self.notebook_tab3) 
        self.Plot1DSpectrum = Plot1DSpectrum(self, self.notebook_tab4) 
        self.DisplayTimeFrames = DisplayTimeFrames(self)
        
        # Save the plots
        self.Plots = [self.PlotTimeEvolution, self.PlotParallelModeStructure, self.PlotSpectrumKxKy, self.Plot1DSpectrum]
        if True: return
    
################################################################################
#                                     METHODS
################################################################################

    def plot(self):
        ''' From <Simulations> we can ask to plot the figure. ''' 
        if self.frame_canvas1.grid_info()!={}:   
            self.PlotTimeEvolution.plotOnTheGUI()
        if self.frame_canvas2.grid_info()!={}:  
            self.PlotParallelModeStructure.plotOnTheGUI()
        if self.frame_canvas3.grid_info()!={}:  
            self.PlotSpectrumKxKy.plotOnTheGUI() 
        if self.frame_canvas4.grid_info()!={}:  
            self.Plot1DSpectrum.plotOnTheGUI() 
        if self.frame_canvas5.grid_info()!={}:  
            self.Plot1DSpectrum.plotOnTheGUI() 
        return
        
    #-------------------------------------  
    def reset_axes(self):
        ''' From <Simulations> we can ask to reset the figures. '''
        if self.PlotTimeEvolution.initiated_canvas: self.PlotTimeEvolution.Canvas.Axis.reset_axis()
        if self.PlotParallelModeStructure.initiated_canvas: self.PlotParallelModeStructure.Canvas.Axis.reset_axis()
        if self.PlotSpectrumKxKy.initiated_canvas: self.PlotSpectrumKxKy.Canvas.Axis.reset_axis()
        if self.Plot1DSpectrum.initiated_canvas: self.Plot1DSpectrum.Canvas1.Axis.reset_axis()
        if self.Plot1DSpectrum.initiated_canvas: self.Plot1DSpectrum.Canvas2.Axis.reset_axis()
        return 
    
    #-------------------------------------  
    def display_canvas(self, identifier):
        ''' Attach the correct frame with the appropriate <Canvas> to the GUI. '''
        if identifier=="TabNonlinearSpatial:PlotTimeEvolution" and self.frame_canvas1.grid_info()=={}:
            self.frame_canvas1.grid() 
            self.frame_canvas4.grid_remove()  
            self.frame_canvas5.grid_remove()  
            if self.frame_canvas2.grid_info()=={} and self.frame_canvas3.grid_info()=={}: self.frame_canvas2.grid() 
        if identifier=="TabNonlinearSpatial:PlotParallelModeStructure" and self.frame_canvas2.grid_info()=={} or self.frame_canvas3.grid_info()!={}:
            self.frame_canvas1.grid() 
            self.frame_canvas2.grid() 
            self.frame_canvas3.grid_remove() 
            self.frame_canvas4.grid_remove()  
            self.frame_canvas5.grid_remove()  
        if identifier=="TabNonlinearSpatial:PlotSpectrumKxKy" and self.frame_canvas3.grid_info()=={}:
            self.frame_canvas1.grid() 
            self.frame_canvas3.grid()  
            self.frame_canvas2.grid_remove() 
            self.frame_canvas4.grid_remove() 
            self.frame_canvas5.grid_remove() 
        if identifier=="TabNonlinearSpatial:Plot1DSpectrum" and self.frame_canvas4.grid_info()=={}:
            self.frame_canvas4.grid() 
            self.frame_canvas5.grid() 
            self.frame_canvas1.grid_remove()  
            self.frame_canvas2.grid_remove()  
            self.frame_canvas3.grid_remove()  
        if identifier=="TabNonlinearSpatial:Plot1DSpectrum" and self.frame_canvas5.grid_info()=={}:
            self.frame_canvas4.grid() 
            self.frame_canvas5.grid() 
            self.frame_canvas1.grid_remove()  
            self.frame_canvas2.grid_remove()  
            self.frame_canvas3.grid_remove()  
        return 
    
################################################################################
#                          LOAD TAB AND CANVASSES
################################################################################

    def load_tab(self, *_):
        '''  When the tab is visible, make sure the correct figure is loaded. 
        When the frame of the canvas is visible, make sure the canvas is initiated, 
        so that we have the figure and axis object. '''       
         
        # Connect the Progress widget of this tab to the research object
        self.root.Progress = self.Progress
             
        # Make sure the dropdown menus for the experiments/simulations are updated
        if self.root.Research.input_files != []: 
            self.Simulations.display_plottedModesAndExperiments()

        # Make sure the time frame is filled in 
        self.DataSelection.Frame_TSelection.set_defaultTimeFrame(self.DataSelection.TSelection) 
        return
         
    #-------------------------------------  
    # When the frame of the canvas is visible, make sure the canvas is
    # initiated, so that we have the figure and axis object.  
    def load_canvas1(self, *_): self.PlotTimeEvolution.initiate_canvas(self.frame_canvas1, self.notebook_tab1)  
    def load_canvas2(self, *_): self.PlotParallelModeStructure.initiate_canvas(self.frame_canvas2, self.notebook_tab2) 
    def load_canvas3(self, *_): self.PlotSpectrumKxKy.initiate_canvas(self.frame_canvas3, self.notebook_tab3) 
    def load_canvas4(self, *_): self.Plot1DSpectrum.initiate_canvas(self.frame_canvas4, self.frame_canvas5, self.notebook_tab4) 
    def load_canvas5(self, *_): self.Plot1DSpectrum.initiate_canvas(self.frame_canvas4, self.frame_canvas5, self.notebook_tab4) 

