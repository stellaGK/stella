
import numpy as np
import tkinter as tk
from tkinter import ttk    
from stellapy.simulations.utils.get_simulations import get_simulations
from stellapy.simulations.utils.get_experiments import get_experiments 
from stellapy.GUI.widgets import PAD_LABEL2, PAD_ENTRY2
        
#################################################################
#                  CLASS TO CHOOSE SIMULATIONS
#################################################################
   
class DataSelection: 
    
    def __init__(self, parent, frame):
        ''' This options frame controls which time and z frame are chosen. '''
        
        # Make the parents available
        self.tab = parent
        self.root = parent.root
        self.frame = frame 
        
        # Attributes 
        self.t_specific = None
        self.z_specific = None
        
        # Create a subframe for each set of plotting options  
        self.subframe_zSelection = ttk.Frame(frame)
        self.subframe_tSelection = ttk.Frame(frame)
         
        # Configure the main frame to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1) # subframe_tSelection
        tk.Grid.rowconfigure(   frame, 1, weight=1) # subframe_zSelection
        tk.Grid.columnconfigure(frame, 0, weight=1) 
 
        # Add the subframes to the main frame  
        self.subframe_zSelection.grid(row=0, column=0, padx=(2,2), pady=(0,2), stick='NSEW') 
        self.subframe_tSelection.grid(row=1, column=0, padx=(2,2), pady=(2,0), stick='NSEW')  
        
        # Fill the frames with their options  
        self.ZSelection = self.Frame_ZSelection(self, self.subframe_zSelection) 
        self.TSelection = self.Frame_TSelection(self, self.subframe_tSelection)  
        if True: return


#################################################################
#                          METHODS
#################################################################

    class Frame_ZSelection:
         
        def __init__(self, parent, frame):
             
            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
            self.root = parent.tab.root
         
            # Configure the frame <self.subframe_zSelection>
            tk.Grid.columnconfigure(frame, 0, weight=0) 
            tk.Grid.columnconfigure(frame, 1, weight=0) 
            tk.Grid.columnconfigure(frame, 2, weight=1) 
            
            # Variables
            font = ("Courier New", 11); width=5
            self.var_z = tk.IntVar(value=1)
            self.var_zSpecific = tk.StringVar(value="0")
             
            # Create the widgets
            self.rbn_zFrame       = ttk.Radiobutton(frame, text='  Average over the z-axis')
            self.rbn_zSpecific    = ttk.Radiobutton(frame, text='  At z = ')
            self.ent_zSpecific    = ttk.Entry(frame, textvariable=self.var_zSpecific, font=font, style='opt_valueCBold.TEntry', width=width)
            self.ent_zSpecific.bind('<Return>', self.change_zSelection)
     
            # Add the values and commands to the radiobuttons 
            self.rbn_zFrame.config(    value=1, variable=self.var_z, command=self.change_zSelection)
            self.rbn_zSpecific.config( value=3, variable=self.var_z, command=self.change_zSelection) 
             
            # Add the options to the frame
            self.rbn_zFrame.grid(   row=1, column=0, **PAD_LABEL2, columnspan=3)
            self.rbn_zSpecific.grid(row=2, column=0, **PAD_LABEL2)
            self.ent_zSpecific.grid(row=2, column=1, **PAD_ENTRY2)
            return 
        
        #---------------------
        def change_zSelection(self, *_):
             
            # When changing the plots, reset the axis of the graph
            self.tab.Graph[self.tab.PlotParallelModeStructure.graph_id].load_defaults()
            self.tab.Graph[self.tab.PlotParallelModeStructure.graph_id].update_quantitiesAndKeys()
            self.tab.Graph[self.tab.PlotSpectrumKxKy.graph_id].load_defaults()
            self.tab.Graph[self.tab.PlotSpectrumKxKy.graph_id].update_quantitiesAndKeys()
             
            # Extract what needs to be plotted on the y-axis
            z_quantity = self.var_z.get()
            if z_quantity==1: self.plot.z_specific = None
            if z_quantity==2: self.plot.z_specific = "max"
            if z_quantity==3: self.plot.z_specific = float(self.var_zSpecific.get())
            
            # Calculate where |phi**2|(z) is maximum along z
            if z_quantity==2:
                self.calculate_zIndexWherePhi2IsMaximum()
                
            # Change the header of the frame and plot the graph
            self.tab.plot()
            return
    
    #---------------------- 
    class Frame_TSelection:
         
        def __init__(self, parent, frame):
             
            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
            self.root = parent.tab.root
         
            # Configure the frame <self.subframe_zSelection>
            tk.Grid.columnconfigure(frame, 0, weight=0) 
            tk.Grid.columnconfigure(frame, 1, weight=1)   
            
            # Variables
            font = ("Courier New", 11)
            self.var_time = tk.IntVar(value=1)
            self.var_timeSpecific = tk.StringVar(value="-") 
            
            # Get the default time frame: the time frame that was chosen
            self.set_defaultTimeFrame()
            
            # Create the widgets
            self.rbn_timeFrame    = ttk.Radiobutton(frame, text='  Average over the selected time frames')
            self.rbn_timeSpecific = ttk.Radiobutton(frame, text='  At specific t = '); width=8
            self.ent_timeSpecific = ttk.Entry(frame, textvariable=self.var_timeSpecific, font=font, style='opt_valueCBold.TEntry', width=width)
            self.ent_timeSpecific.bind('<Return>', self.change_timeSelection)

            # Add the values and commands to the radiobuttons 
            self.rbn_timeFrame.config(    value=1, variable=self.var_time, command=self.change_timeSelection)
            self.rbn_timeSpecific.config( value=2, variable=self.var_time, command=self.change_timeSelection) 
             
            # Add the options to the frame
            self.rbn_timeFrame.grid(   row=2, column=0, **PAD_LABEL2, columnspan=2)
            self.rbn_timeSpecific.grid(row=3, column=0, **PAD_LABEL2) 
            self.ent_timeSpecific.grid(row=3, column=1, **PAD_ENTRY2)
            return 
        
        #--------------------
        def set_defaultTimeFrame(self):
            if self.root.TabSelectedFiles.input_files != []:
                experiment = get_experiments(self.root.Research, self.tab.Simulations.experiment_id)[0]
                simulation = get_simulations(experiment, self.tab.Simulations.simulation_id)[0]
                t_range = [simulation.time.tstart, simulation.time.tend]
                self.var_timeSpecific.set(str(t_range[0]))
                if self.var_time.get()==1: 
                    self.plot.t_specific = None
                if self.var_time.get()==2: 
                    self.plot.t_specific = t_range[0]
            return 
         
        #---------------------
        def change_timeSelection(self, *_):
             
            # When changing the plots, reset the axis of the graph
            self.tab.reset_axes()
             
            # Extract what needs to be plotted on the y-axis
            t_button   = self.var_time.get()
            t_specific = float(self.var_timeSpecific.get())
            if t_button==1: self.plot.t_specific = None 
            if t_button==2: self.plot.t_specific = t_specific 
             
            # Replot the graphs
            self.tab.plot()
            return
    

