#  
# #################################################################
# #                 CLASS FOR THE GEOMETRY TAB
# #################################################################
# ''' 
# '''
# 
# # Load modules
# import tkinter as tk
# from tkinter import ttk 
# from stellapy.GUI.widgets import Progress
# from stellapy.GUI.interface.TabGeometry.Canvasses import Canvasses 
# from stellapy.GUI.interface.TabGeometry.Simulations import Simulations  
# from stellapy.GUI.interface.TabGeometry.PlotVariablesVsZ import PlotVariablesVsZ
# 
# #################################################################
# #                   CLASS FOR THE FIRST TAB
# #################################################################
# class TabGeometry:
#     ''' This class governs the tab "Geometry".
#     
#     There are four frames attached to display the following graphs:
#         graph_id = 0: Variable versus z
#         graph_id = 1: ???
#     '''
#  
# #################################################################
# #                          WIDGETS
# #################################################################
# 
#     # Initate the tab for "Nonlinear time traces"
#     def __init__(self, tab): 
# 
#         #======================================================
#         # Safe info from the tab and root that can be passed on
#         #======================================================
# 
#         self.window = tab                      # Direct accesses to the window of tab "Nonlinear time traces" 
#         self.root = tab.root                   # Needed to center windows based on the root screen
#         self.dict_awthemes = tab.root.awthemes # To make the tk widgets look like the ttk widgets
# 
#         #===========
#         # VARIABLES
#         #===========
#         
#         # Only initiate the canvasses when the tab is used
#         self.initiated_canvas = None
#         
#         # Store the Graph and Canvas objects in a list
#         self.Graph = [None, None]     
#         self.Canvas = [None, None]   
#         
#         #===========
#         # FRAMES
#         #===========
# 
#         # Create a frame for each graph and one for the options and add them to the tab "Nonlinear time traces" 
#         self.frame_graph0   = ttk.Frame(self.window)
#         self.frame_graph1   = ttk.Frame(self.window) 
#         self.frame_options  = ttk.Frame(self.window) 
# 
#         # Configure the columns and rows of the frames     
#         tk.Grid.rowconfigure(   self.window, 0, weight=0) 
#         tk.Grid.rowconfigure(   self.window, 1, weight=1,  uniform="please") 
#         tk.Grid.rowconfigure(   self.window, 2, weight=9,  uniform="please") 
#         tk.Grid.rowconfigure(   self.window, 3, weight=1,  uniform="please") 
#         tk.Grid.columnconfigure(self.window, 0, weight=1,  uniform="for sure") 
#         tk.Grid.columnconfigure(self.window, 1, weight=1,  uniform="for sure")   
#  
#         # Add the frames to the main window (= tab "Nonlinear time traces" )
#         self.frame_options.grid(in_=self.window, row=0, column=0, padx=(20,20), pady=(5,10), stick='NSEW', columnspan=2)
#         self.frame_graph0.grid( in_=self.window, row=2, column=0, padx=(20,10), pady=(10,20),stick='NSEW')
#         self.frame_graph1.grid( in_=self.window, row=2, column=1, padx=(10,20), pady=(10,20),stick='NSEW')  
#         
#         # Bind focus on every click on a widget: this makes sure the entry widgets become unfocused
#         self.frame_graph0.bind("<1>", lambda event: self.frame_graph1.focus_set())
#         self.frame_graph1.bind("<1>", lambda event: self.frame_graph2.focus_set()) 
#         self.frame_options.bind("<1>", lambda event: self.frame_options.focus_set()) 
#         
#         #==============================================================================
#         # FRAMES FOR THE OPTIONS: SIMULATIONS FRAME, PROGRESS FRAME AND TABS PER GRAPH
#         #==============================================================================
#     
#         # Create the subframes in the options frame
#         self.subframe_simulations   = ttk.LabelFrame(self.frame_options, text="   Simulations  ", **self.dict_awthemes['labelframe2'])
#         self.subframe_tabs          = ttk.Frame(self.frame_options)
#         self.frame_progress         = ttk.Frame(self.frame_options)
# 
#         # Add the progress object
#         self.Progress = Progress(self, anchor="fill", length=300)
#         self.Progress.move(0,"Please select simulations.")
#         
#         # Configure the frame
#         tk.Grid.rowconfigure(   self.frame_options, 0, weight=1)  
#         tk.Grid.rowconfigure(   self.frame_options, 1, weight=0)     
#         tk.Grid.columnconfigure(self.frame_options, 0, weight=0)   
#         tk.Grid.columnconfigure(self.frame_options, 1, weight=1)     
#         
#         # Add the labelframes to <frame_options>
#         self.subframe_simulations.grid(  row=0, column=0, padx=(0,5), pady=(0,0), stick='NSEW')
#         self.frame_progress.grid(        row=1, column=0, padx=(0,5), pady=(5,0), stick='NSEW') 
#         self.subframe_tabs.grid(         row=0, column=1, padx=(5,0), pady=(0,0), stick='NSEW', rowspan=2)
#         
#         # Fill the subframes with widgets
#         self.Simulations = Simulations(self, self.subframe_simulations) 
#         
#         #===========================================
#         # EACH GRAPH HAS A TAB IN THE OPTIONS FRAME
#         #===========================================
#         
#         # Tabbed display for ....
#         self.tab_plots = ttk.Notebook(self.subframe_tabs, style='species.TNotebook')  
#         self.tab_plots.grid(row=0, column=0, sticky="NSEW") 
#         tk.Grid.rowconfigure(self.subframe_tabs, 0, weight=1)  
#         tk.Grid.columnconfigure(self.subframe_tabs, 0, weight=1)  
# 
#         # Create frames for each tab
#         self.tab_graph0 = ttk.Frame(self.tab_plots, padding=(0,0,0,0))   
#         self.tab_graph1 = ttk.Frame(self.tab_plots, padding=(0,0,0,0))   
#         self.tab_plots.add(self.tab_graph0, text="Variables along the flux tube")  
#         self.tab_plots.add(self.tab_graph1, text="...")   
#         self.tab_plots.select(self.tab_graph0)
#                
#         # Fill the subframes with widgets
#         self.PlotVariablesVsZ = PlotVariablesVsZ(self, self.tab_graph0, 0)  
#         
#         #=============
#         # LOAD FIGURE
#         #=============
#         
#         # When the tab is visible, make sure the correct figure is loaded
#         self.frame_options.bind("<Visibility>", self.load_figure)
#         if True: return
#     
# #################################################################
# #                          METHODS
# #################################################################
# 
#     def load_figure(self, *args):
#         '''  When the tab is visible, make sure the correct figure is loaded. 
#         If the simulations were changed, replot the graph. '''        
#         
#         # Load the canvasses only when the tab is visible, to save time.
#         if self.initiated_canvas==None:
#             self.CanvasClass = Canvasses(self)
#             self.initiated_canvas = True
#         
#         # Connect the Progress class to the research object
#         self.root.Progress = self.Progress
#         for experiment in self.root.Research.experiments:
#             experiment.Progress = self.Progress
#             for simulation in experiment.simulations:
#                 simulation.Progress = self.Progress
#         
#         # Make sure we display the correct experiments
#         if self.input_files != []:
#             self.Simulations.display_plottedModesAndExperiments() 
#         self.show_progress("finished", None)
#         
#         # Make sure the titles are set
#         self.Simulations.change_yQuantity(plot=False)
#         return
#     
#     #-------------------------      
#     def show_progress(self, status="start", poppedout_id=None):
#         # While plotting show a progress bar turn the cursor into a waiting cursor.
#         if status=="start" and poppedout_id == None: 
#             self.Progress.move(0,"Start plotting.")
#             self.root.config(cursor="watch")
#             self.root.update_idletasks() 
#         # When finished show the canvas and change the cursor back.
#         if status=="finished" and poppedout_id == None: 
#             self.Progress.move(100, "No tasks are running.")
#             self.root.config(cursor="")
#         # When there are no simulations revert to the loading bar
#         if status=="nothing" and poppedout_id == None: 
#             self.Progress.move(0,"Please select simulations.")
#         return
# 
#     #-------------------------------------  
#     def reset_graph(self, axis_id, poppedout_id=None):
#         self.CanvasClass.reset_graph(axis_id, poppedout_id)
#         
#     #-------------------------------------  
#     def popout_window(self, axis_id):
#         self.CanvasClass.popout_window(axis_id)
#         
# 


