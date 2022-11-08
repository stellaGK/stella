 
 
import tkinter as tk
from tkinter import ttk 
from stellapy.plot.utils.labels import standardParameters
from stellapy.GUI.widgets import Progress 
from stellapy.GUI.widgets.DisplayTimeFrames import DisplayTimeFrames
from stellapy.GUI.interface.TabNonlinearTime.Simulations import Simulations
from stellapy.GUI.interface.TabNonlinearTime.PlotTimeEvolution import PlotTimeEvolution
from stellapy.GUI.interface.TabNonlinearTime.PlotSaturatedFluxVersusParameter import PlotSaturatedFluxVersusParameter  
from stellapy.GUI.interface.TabNonlinearTime.PlotNonlinearMap import PlotNonlinearMap
 
################################################################################
#            CLASS TO PLOT NONLINEAR QUANTITIES IN FUNCTION OF TIME
################################################################################

class TabNonlinearTime:
    ''' Initiate the frames on the tab "NonlinearTime" to investigate the results 
    of nonlinear simulations with stella. The tab can plot the following graphs: 
            (1) Time evolution
            (3) Satured flux versus parameter
            (4) Nonlinear map ''' 

    def __init__(self, tab): 
 
        # Make the parent objects available      
        self.dict_awthemes = tab.root.awthemes             
        self.root = tab.root    
 
        # Create a frame for each canvas and a frame for the options  
        self.frame_canvas1 = ttk.LabelFrame(tab, text="   Time evolution  ", **self.dict_awthemes['labelframe2'])
        self.frame_canvas2 = ttk.LabelFrame(tab, text="   Satured flux versus parameter  ", **self.dict_awthemes['labelframe2'])
        self.frame_canvas3 = ttk.LabelFrame(tab, text="   Nonlinear map  ", **self.dict_awthemes['labelframe2'])
        self.frame_options = ttk.Frame(tab) 
 
        # Configure the rows/columns of the grid on <tabNonlinearTime>
        tk.Grid.rowconfigure(   tab, 0, weight=0) 
        tk.Grid.rowconfigure(   tab, 1, weight=1) 
        tk.Grid.columnconfigure(tab, 0, weight=1,  uniform="for sure") 
        tk.Grid.columnconfigure(tab, 1, weight=1,  uniform="for sure") 
  
        # Add the frames to the main window (= tab "Nonlinear time traces" )
        self.frame_options.grid(in_=tab, row=0, column=0, padx=(20,20), pady=(5,10), stick='NSEW', columnspan=2)
        self.frame_canvas1.grid(in_=tab, row=1, column=0, padx=(20,10), pady=(10,20),stick='NSEW')
        self.frame_canvas2.grid(in_=tab, row=1, column=1, padx=(10,20), pady=(10,20),stick='NSEW')
        self.frame_canvas3.grid(in_=tab, row=1, column=1, padx=(10,20), pady=(10,20),stick='NSEW') 
        self.frame_canvas3.grid_remove()  
         
        # Bind focus on every click on a widget: this makes sure the entry widgets become unfocused
        self.frame_canvas1.bind("<1>", lambda _: self.frame_canvas1.focus_set())
        self.frame_canvas2.bind("<1>", lambda _: self.frame_canvas2.focus_set())
        self.frame_canvas3.bind("<1>", lambda _: self.frame_canvas3.focus_set()) 
        self.frame_options.bind("<1>", lambda _: self.frame_options.focus_set())  

        # Fill the options frame with subframes, and create tabs for the graph options
        self.fill_optionsFrame(self.frame_options)
         
        # When <TabNonlinearTime> is visible, make sure the Progress object is linked,
        # that the dropdown menus are updated, and that the canvasses are initiated.
        self.frame_options.bind("<Visibility>", self.load_tab) 
        self.frame_canvas1.bind("<Visibility>", self.load_canvas1)
        self.frame_canvas2.bind("<Visibility>", self.load_canvas2)
        self.frame_canvas3.bind("<Visibility>", self.load_canvas3) 
        return 

    #-----------------------------
    def fill_optionsFrame(self, frame):
        ''' Create a subframe for the dropdown menus of the experiments and 
        simulations. Create a subframe for the Progress object. Create a tabbed
        widget, with a tab for each graph to hold the various options. '''
        
        # Create the subframes in the options frame
        self.subframe_simulations = ttk.LabelFrame(frame, text="   Simulations  ", **self.dict_awthemes['labelframe2'])
        self.subframe_tabs        = ttk.Frame(frame)
        self.frame_progress       = ttk.Frame(frame)
 
        # Add the progress object
        self.Progress = Progress(self.root, self.frame_progress, anchor="fill", length=300) 
        self.Progress.move(0,"Please select simulations.")
         
        # Configure the frame
        tk.Grid.rowconfigure(   frame, 0, weight=1)  
        tk.Grid.rowconfigure(   frame, 1, weight=0)     
        tk.Grid.columnconfigure(frame, 0, weight=0)   
        tk.Grid.columnconfigure(frame, 1, weight=1)  
         
        # Add the labelframes to <frame_options>
        self.subframe_simulations.grid(row=0, column=0, padx=(0,5), pady=(0,0), stick='NSEW')
        self.frame_progress.grid(      row=1, column=0, padx=(0,5), pady=(5,0), stick='NSEW')
        self.subframe_tabs.grid(       row=0, column=1, padx=(5,0), pady=(0,0), stick='NSEW', rowspan=2)
         
        # Tabbed display for "PlotTimeEvolution", "PlotSaturatedFluxVersusParameter" and "PlotNonlinearMap"
        self.notebook = ttk.Notebook(self.subframe_tabs, style='species.TNotebook')  
        self.notebook.grid(row=0, column=0, sticky="NSEW") 
        tk.Grid.rowconfigure(self.subframe_tabs, 0, weight=1)  
        tk.Grid.columnconfigure(self.subframe_tabs, 0, weight=1)  
 
        # Create frames for each tab
        self.notebook_tab1 = ttk.Frame(self.notebook, padding=(0,0,0,0)) 
        self.notebook_tab2 = ttk.Frame(self.notebook, padding=(0,0,0,0)) 
        self.notebook_tab3 = ttk.Frame(self.notebook, padding=(0,0,0,0))  
        self.notebook.add(self.notebook_tab1, text="Time evolution")   
        self.notebook.add(self.notebook_tab2, text="Saturated flux versus parameter") 
        self.notebook.add(self.notebook_tab3, text="Nonlinear map")  

        # Fill <subframe_simulations> with the dropdown menus for the experiment/simulations
        self.Simulations = Simulations(self, self.subframe_simulations)
        
        # Fill <notebook> with the options for each graph
        self.PlotTimeEvolution = PlotTimeEvolution(self, self.notebook_tab1)
        self.PlotSaturatedFluxVersusParameter = PlotSaturatedFluxVersusParameter(self, self.notebook_tab2)
        self.PlotNonlinearMap = PlotNonlinearMap(self, self.notebook_tab3)
        self.DisplayTimeFrames = DisplayTimeFrames(self)
        if True: return
        
################################################################################
#                                     METHODS
################################################################################

    def plot(self):
        ''' From <Simulations> we can ask to plot the figure. ''' 
        if self.frame_canvas1.grid_info()!={}:   
            self.PlotTimeEvolution.plotOnTheGUI()
        if self.frame_canvas2.grid_info()!={}:  
            self.PlotSaturatedFluxVersusParameter.plotOnTheGUI()
        if self.frame_canvas3.grid_info()!={}:  
            self.PlotNonlinearMap.plotOnTheGUI() 
        return
        
    #-------------------------------------  
    def reset_axes(self):
        ''' From <Simulations> we can ask to reset the figures. '''
        if self.PlotTimeEvolution.initiated_canvas: self.PlotTimeEvolution.Canvas.Axis.reset_axis()
        if self.PlotSaturatedFluxVersusParameter.initiated_canvas: self.PlotSaturatedFluxVersusParameter.Canvas.Axis.reset_axis()
        if self.PlotNonlinearMap.initiated_canvas: self.PlotNonlinearMap.Canvas.Axis.reset_axis()
        return 
    
    #-------------------------------------  
    def display_canvas(self, identifier):
        ''' Attach the correct frame with the appropriate <Canvas> to the GUI. '''
        if identifier=="TabNonlinearTime:PlotSaturatedFluxVersusParameter" and self.frame_canvas2.grid_info()=={}:
            self.frame_canvas2.grid() 
            self.frame_canvas3.grid_remove()  
        if identifier=="TabNonlinearTime:PlotNonlinearMap" and self.frame_canvas3.grid_info()=={}:
            self.frame_canvas3.grid() 
            self.frame_canvas2.grid_remove()  
        return 
    
    #-------------------------------------  
    def detect_theVariedParameter(self):
        ''' Try to automatically detect the stella <knob> and <key> that is varied. ''' 
#         try:
        # Detect the varied parameter for the following plotting class
        Plot = self.PlotSaturatedFluxVersusParameter
        var_parameter = Plot.Parameter.var_parameter
        if Plot.knob=="-" and len(self.root.Research.experiments)!=0:
            # Get the varied variables per simulation
            if len(self.Simulations.experiment.variedVariables)!=0: 
                if len(self.Simulations.experiment.variedVariables)>=3: variable = self.Simulations.experiment.variedVariables[2]
                if len(self.Simulations.experiment.variedVariables)<3: variable = self.Simulations.experiment.variedVariables[0]
                Plot.knob = variable.split(":")[0]
                Plot.key = variable.split(": ")[-1]
            # Get the varied variables per experiment
            if len(self.Simulations.experiment.variedVariables)==0:
                Plot.key = self.root.Research.experiment_key
                Plot.knob = standardParameters[Plot.key]["knob"]
            # Replace the key with the label in the dropdown menu
            if Plot.key=="tprim": Plot.key = "tiprim"
            if Plot.key=="delt":  Plot.key = "delta t"
            # Update the parameter of <PlotSaturatedFluxVersusParameter> to the detected parameter
            var_parameter.trace_vdelete("w", var_parameter.trace_id) 
            var_parameter.trace_id = var_parameter.trace('w', Plot.Parameter.update_parameterWithoutPlotting) 
            var_parameter.set(Plot.key)
            var_parameter.trace_vdelete("w", var_parameter.trace_id) 
            var_parameter.trace_id = var_parameter.trace('w', Plot.Parameter.update_parameter) 
#         except: pass
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

        # Detect the varied variables
        self.detect_theVariedParameter()
        return

    #-------------------------------------  
    # When the frame of the canvas is visible, make sure the canvas is
    # initiated, so that we have the figure and axis object.  
    def load_canvas1(self, *_): self.PlotTimeEvolution.initiate_canvas(self.frame_canvas1, self.notebook_tab1)  
    def load_canvas2(self, *_): self.PlotSaturatedFluxVersusParameter.initiate_canvas(self.frame_canvas2, self.notebook_tab2) 
    def load_canvas3(self, *_): self.PlotNonlinearMap.initiate_canvas(self.frame_canvas3, self.notebook_tab3) 

