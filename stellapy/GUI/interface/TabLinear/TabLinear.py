 
import tkinter as tk
from tkinter import ttk  
from stellapy.GUI.widgets.Progress import Progress
from stellapy.GUI.interface.TabLinear.Simulations import Simulations
from stellapy.GUI.interface.TabLinear.PlotLinearMap import PlotLinearMap
from stellapy.GUI.interface.TabLinear.PlotLinearSpectrum import PlotLinearSpectrum
from stellapy.GUI.interface.TabLinear.PlotParameterInfluence import PlotParameterInfluence    
from stellapy.GUI.interface.TabLinear.PlotVelocityDistribution import PlotVelocityDistribution
from stellapy.plot.utils.labels.standardParameters import standardParameters

################################################################################
#              CLASS TO STUDY THE RESULTS OF LINEAR SIMULATIONS
################################################################################

class TabLinear:
    ''' Initiate the frames on the tab "Linear" to investigate the 
    results of linear simulations with stella.
    The tab can plot the following graphs: 
            (1) ... ''' 
 
    def __init__(self, tab): 
 
        # Make the parent objects available      
        self.dict_awthemes = tab.root.awthemes             
        self.root = tab.root    
 
        # Create a frame for each canvas and a frame for the options  
        self.frame_canvas1 = ttk.LabelFrame(tab, text="   Linear spectrum  ",          **self.dict_awthemes['labelframe2'])
        self.frame_canvas2 = ttk.LabelFrame(tab, text="   Parameter influence  ",      **self.dict_awthemes['labelframe2'])
        self.frame_canvas3 = ttk.LabelFrame(tab, text="   Linear map  ",               **self.dict_awthemes['labelframe2'])
        self.frame_canvas4 = ttk.LabelFrame(tab, text="   Parallel velocity  ",        **self.dict_awthemes['labelframe2'])
        self.frame_canvas5 = ttk.LabelFrame(tab, text="   Perpendicular velocity  ",   **self.dict_awthemes['labelframe2'])
        self.frame_options = ttk.Frame(tab) 
 
        # Configure the rows/columns of the grid on <tabLinear>
        tk.Grid.rowconfigure(   tab, 0, weight=0) 
        tk.Grid.rowconfigure(   tab, 1, weight=1) 
        tk.Grid.columnconfigure(tab, 0, weight=1,  uniform="for sure") 
        tk.Grid.columnconfigure(tab, 1, weight=1,  uniform="for sure") 
  
        # Add the canvasses and frames to <tabLinear>
        self.frame_options.grid(in_=tab, row=0, column=0, padx=(20,20), pady=(20,10),stick='NSEW', columnspan=2)
        self.frame_canvas1.grid(in_=tab, row=1, column=0, padx=(20,10), pady=(10,20),stick='NSEW')
        self.frame_canvas2.grid(in_=tab, row=1, column=1, padx=(10,20), pady=(10,20),stick='NSEW')
        self.frame_canvas3.grid(in_=tab, row=1, column=1, padx=(10,20), pady=(10,20),stick='NSEW')
        self.frame_canvas4.grid(in_=tab, row=1, column=0, padx=(20,10), pady=(10,20),stick='NSEW')
        self.frame_canvas5.grid(in_=tab, row=1, column=1, padx=(10,20), pady=(10,20),stick='NSEW')
        self.frame_canvas3.grid_remove() 
        self.frame_canvas4.grid_remove()  
        self.frame_canvas5.grid_remove()  
         
        # Bind focus on every click on a widget: this unfoces the entry widgets 
        self.frame_canvas1.bind("<1>", lambda _: self.frame_canvas1.focus_set())
        self.frame_canvas2.bind("<1>", lambda _: self.frame_canvas2.focus_set())
        self.frame_canvas3.bind("<1>", lambda _: self.frame_canvas3.focus_set())
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
         
        # Tabbed display for "PlotLinearSpectrum", "PlotParameterInfluence", "PlotLinearMap" and "PlotVelocityDistribution"
        self.notebook = ttk.Notebook(self.subframe_tabs, style='species.TNotebook')  
        self.notebook.grid(row=0, column=0, sticky="NSEW") 
        tk.Grid.rowconfigure(self.subframe_tabs, 0, weight=1)  
        tk.Grid.columnconfigure(self.subframe_tabs, 0, weight=1)  
 
        # Create frames for each tab
        self.notebook_tab1 = ttk.Frame(self.notebook, padding=(0,0,0,0)) 
        self.notebook_tab2 = ttk.Frame(self.notebook, padding=(0,0,0,0)) 
        self.notebook_tab3 = ttk.Frame(self.notebook, padding=(0,0,0,0)) 
        self.notebook_tab4 = ttk.Frame(self.notebook, padding=(0,0,0,0)) 
        self.notebook.add(self.notebook_tab1, text="Linear spectrum") 
        self.notebook.add(self.notebook_tab2, text="Parameter influence") 
        self.notebook.add(self.notebook_tab3, text="Parameter map") 
        self.notebook.add(self.notebook_tab4, text="Velocity distribution") 

        # Fill <subframe_simulations> with the dropdown menus for the experiment/simulations
        self.Simulations = Simulations(self, self.subframe_simulations)
                
        # Fill <notebook> with the options for each graph
        self.PlotLinearSpectrum = PlotLinearSpectrum(self, self.notebook_tab1)
        self.PlotParameterInfluence = PlotParameterInfluence(self, self.notebook_tab2)
        self.PlotLinearMap = PlotLinearMap(self, self.notebook_tab3)
        self.PlotVelocityDistribution = PlotVelocityDistribution(self, self.notebook_tab4)
        
        # Save the plots
        self.Plots = [self.PlotLinearSpectrum, self.PlotParameterInfluence, self.PlotLinearMap, self.PlotVelocityDistribution]
        if True: return
     
################################################################################
#                                     METHODS
################################################################################
 
    def plot(self):
        ''' From <Simulations> we can ask to plot the figure. ''' 
        if self.frame_canvas1.grid_info()!={}:   
            self.PlotLinearSpectrum.plotOnTheGUI()
        if self.frame_canvas2.grid_info()!={}:  
            self.PlotParameterInfluence.plotOnTheGUI()
        if self.frame_canvas3.grid_info()!={}:  
            self.PlotLinearMap.plotOnTheGUI()
        if self.frame_canvas4.grid_info()!={}:  
            self.PlotVelocityDistribution.plotOnTheGUI()
        return
        
    #-------------------------------------  
    def reset_axes(self):
        ''' From <Simulations> we can ask to reset the figures. '''
        if self.PlotLinearMap.initiated_canvas: self.PlotLinearMap.Canvas.Axis.reset_axis()
        if self.PlotLinearSpectrum.initiated_canvas: self.PlotLinearSpectrum.Canvas.Axis.reset_axis()
        if self.PlotParameterInfluence.initiated_canvas: self.PlotParameterInfluence.Canvas.Axis.reset_axis()
        if self.PlotVelocityDistribution.initiated_canvas: self.PlotVelocityDistribution.Canvas1.Axis.reset_axis()
        if self.PlotVelocityDistribution.initiated_canvas: self.PlotVelocityDistribution.Canvas2.Axis.reset_axis()

    #-------------------------------------  
    def display_canvas(self, identifier):
        ''' Attach the correct frame with the appropriate <Canvas> to the GUI. '''
        if identifier=="TabLinear:PlotLinearSpectrums" and self.frame_canvas1.grid_info()=={}:
            self.frame_canvas1.grid() 
            self.frame_canvas4.grid_remove() 
        if identifier=="TabLinear:PlotParameterInfluence" and self.frame_canvas2.grid_info()=={}:
            self.frame_canvas2.grid()
            self.frame_canvas3.grid_remove()
            self.frame_canvas5.grid_remove()
        if identifier=="TabLinear:PlotLinearMap" and self.frame_canvas3.grid_info()=={}:
            self.frame_canvas3.grid()
            self.frame_canvas2.grid_remove()
            self.frame_canvas5.grid_remove()
        if identifier=="TabLinear:PlotVelocityDistribution" and (self.frame_canvas4.grid_info()=={} or self.frame_canvas5.grid_info()=={}):
            self.frame_canvas4.grid()
            self.frame_canvas5.grid()
            self.frame_canvas1.grid_remove()
            self.frame_canvas2.grid_remove()
            self.frame_canvas3.grid_remove()

    #-------------------------------------  
    def detect_theVariedParameter(self):
        ''' Try to automatically detect the stella <knob> and <key> that is varied. '''
        # Detect the varied parameter for the following plotting class 
        Plot = self.PlotParameterInfluence
        var_parameter = Plot.Parameter.var_parameter
        if Plot.knob=="-" and len(self.root.Research.experiments)!=0:
            # Get the varied variables per simulation
            if len(self.Simulations.experiment.variedVariables)!=0: 
                variable = self.Simulations.experiment.variedVariables[0]
                Plot.knob = variable.split(":")[0]
                Plot.key = variable.split(": ")[-1]
            # Get the varied variables per experiment
            if len(self.Simulations.experiment.variedVariables)==0:
                Plot.key = self.root.Research.creationDetails.key1 
                Plot.key = Plot.key if (Plot.key in standardParameters) else "rho"
                Plot.knob = standardParameters[Plot.key]["knob"]
            # Replace the key with the label in the dropdown menu
            if Plot.key=="tprim": Plot.key = "tiprim"
            if Plot.key=="delt":  Plot.key = "delta t"
            # Update the parameter of <PlotParameterInfluence> to the detected parameter
            var_parameter.trace_vdelete("w", var_parameter.trace_id) 
            var_parameter.trace_id = var_parameter.trace('w', Plot.Parameter.update_parameterWithoutPlotting) 
            var_parameter.set(Plot.key)
            var_parameter.trace_vdelete("w", var_parameter.trace_id) 
            var_parameter.trace_id = var_parameter.trace('w', Plot.Parameter.update_parameter) 
        return  
                
    #-------------------------------------  
    def detect_theVariedParameters(self):
        ''' Try to automatically detect the stella <knob> and <key> that is varied. ''' 
        # Detect the varied parameter for the following plotting class 
        Plot = self.PlotLinearMap
        var_parameter1 = Plot.Parameter.var_parameter1
        var_parameter2 = Plot.Parameter.var_parameter2
        if self.PlotLinearMap.knob1 == "-" and len(self.root.Research.experiments)!=0:
            # Get the varied variables per simulation
            if len(self.Simulations.experiment.variedVariables)>1:
                variable1 = self.Simulations.experiment.variedVariables[0]
                variable2 = self.Simulations.experiment.variedVariables[1]
                Plot.knob1 = variable1.split(":")[0];    Plot.knob2 = variable2.split(":")[0]
                Plot.key1  = variable1.split(": ")[-1];  Plot.key2  = variable2.split(": ")[-1]
            # Get the varied variables per experiment
            if len(self.Simulations.experiment.variedVariables)<=1:
                Plot.key1 = self.root.Research.creationDetails.key1; Plot.key1 = "rho" if Plot.key1=="vmec_filename" else Plot.key1 
                Plot.key2 = self.root.Research.creationDetails.key2; Plot.key2 = "rho" if Plot.key2=="vmec_filename" else Plot.key2
                Plot.knob1 = standardParameters[Plot.key1]["knob"]
                Plot.knob2 = standardParameters[Plot.key2]["knob"] 
            # Replace the key with the label in the dropdown menu
            if Plot.key1=="tprim": Plot.key1 = "tiprim"
            if Plot.key1=="delt":  Plot.key1 = "delta t"
            if Plot.key2=="tprim": Plot.key2 = "tiprim"
            if Plot.key2=="delt":  Plot.key2 = "delta t"
            # Update the parameter of <PlotLinearMap> to the detected parameter
            var_parameter1.trace_vdelete("w", var_parameter1.trace_id) 
            var_parameter2.trace_vdelete("w", var_parameter2.trace_id) 
            var_parameter1.trace_id = var_parameter1.trace('w', Plot.Parameter.update_parameterWithoutPlotting) 
            var_parameter2.trace_id = var_parameter1.trace('w', Plot.Parameter.update_parameterWithoutPlotting) 
            var_parameter1.set(Plot.key1) 
            var_parameter2.set(Plot.key2)  
            var_parameter1.trace_vdelete("w", var_parameter1.trace_id) 
            var_parameter2.trace_vdelete("w", var_parameter2.trace_id) 
            var_parameter1.trace_id = var_parameter1.trace('w', Plot.Parameter.update_parameter) 
            var_parameter2.trace_id = var_parameter2.trace('w', Plot.Parameter.update_parameter) 
        return 
                
################################################################################
#                          LOAD TAB AND CANVASSES
################################################################################

    def load_tab(self, *_):
        '''  When the tab is visible, make sure the correct figure is loaded. 
        If the simulations were changed, replot the graph. When the frame of 
        the canvas is visible, make sure the canvas is initiated, so that we 
        have the figure and axis object. '''
 
        # Connect the Progress widget of this tab to the research object
        self.root.Progress = self.Progress
             
        # Make sure the dropdown menus for the experiments/simulations are updated
        if self.root.Research.input_files != []: 
            self.Simulations.display_plottedModesAndExperiments()
     
        # Make sure that the loaded research object had a linear_map object     
        if 'linear_map' not in self.root.Research.data.keys():
            self.root.Research.data['linear_map'] = {
                'parameters1' : None,\
                'parameters2' : None,\
                'gamma' : None,\
                'omega' : None,\
                'ky' : None} 
            
        # Detect the varied variables
        self.detect_theVariedParameter()
        self.detect_theVariedParameters()
        return
    
    #-------------------------------------  
    # When the frame of the canvas is visible, make sure the canvas is
    # initiated, so that we have the figure and axis object.  
    def load_canvas1(self, *_): self.PlotLinearSpectrum.initiate_canvas(self.frame_canvas1)  
    def load_canvas2(self, *_): self.PlotParameterInfluence.initiate_canvas(self.frame_canvas2) 
    def load_canvas3(self, *_): self.PlotLinearMap.initiate_canvas(self.frame_canvas3) 
    def load_canvas4(self, *_): self.PlotVelocityDistribution.initiate_canvas(self.frame_canvas4, self.frame_canvas5) 
    def load_canvas5(self, *_): self.PlotVelocityDistribution.initiate_canvas(self.frame_canvas4, self.frame_canvas5) 

 