# 
# 
# import matplotlib.pyplot as plt
# from stellapy.GUI.widgets import CanvasForGraphs, PoppedOutWindow
# from stellapy.GUI.interface.TabGeometry.Graph import Graph
# 
# #################################################################
# #                    CLASS FOR THE CANVAS
# #################################################################
#    
# class Canvasses: 
#     
#     def __init__(self, parent):
#         ''' Initiate the canvas for the matplotlib plots.
#         
#         The toolbar requires the class to have a reset_graph and popout_window function
#         The optionswindow requires the Graph objects to have attributes: ax, x_name, x_key, y_name, y_key, range, label
#         '''
# 
#         # Make the parent available
#         self.tab = parent
#         self.root = parent.root
#         
#         # Create the figures for each canvas
#         self.figure1 = plt.figure("geometry1") 
#         self.figure2 = plt.figure("geometry2")        
#         self.figure1.set_tight_layout(False) 
#         self.figure2.set_tight_layout(False) 
#         self.figure1.patch.set_facecolor(self.root.color['canvas']) 
#         self.figure2.patch.set_facecolor(self.root.color['canvas'])   
#         
#         # Put the canvas on the screen
#         CanvasForGraphs(self.root, self.tab.frame_graph0, self.root.tab_Geometry, self.figure1, axis_id=0)
#         CanvasForGraphs(self.root, self.tab.frame_graph1, self.root.tab_Geometry, self.figure2, axis_id=1) 
# 
#         # Initiate the plotting class for the main window
#         self.initiate_GraphClass()
# 
#     #--------------------------------------------
#     def initiate_GraphClass(self, figure=None):
#         
#         # Set the color of the axis
#         plt.rcParams['text.color']       = self.root.color['fg']
#         plt.rcParams['axes.edgecolor']   = self.root.color['fg']
#         plt.rcParams['axes.labelcolor']  = self.root.color['fg']
#         plt.rcParams['xtick.color']      = self.root.color['fg']
#         plt.rcParams['ytick.color']      = self.root.color['fg']
#         plt.rcParams['axes.facecolor']   = self.root.color['canvas']
#         plt.rcParams['figure.facecolor'] = self.root.color['canvas'] 
#         
#         # Initiate the plotting class for the main window
#         if figure is None:
#             # Make the axis instance and the required attributes for the options window through the class <graph>
#             self.tab.Graph[0] = Graph(self.tab, self.figure1, 0)
#             self.tab.Graph[1] = Graph(self.tab, self.figure2, 1) 
#             # Put the empty figure on the GUI
#             self.tab.Canvas[0].draw_idle(); self.root.update_idletasks()
#             self.tab.Canvas[1].draw_idle(); self.root.update_idletasks() 
#         
#         # Or initiate the plotting class for a popped out window
#         if figure is not None:
#             # Make the axis instance and the required attributes for the options window through the class <graph>
#             self.root.graph_poppedOut.append(Graph(self.tab, figure, 0))
#             
#         if True: return 
#         
# #################################################################
# #                          METHODS
# #################################################################
# 
#     def popout_window(self, axis_id):
#         '''Replot the figure in a seperate window when the button on the toolbar is clicked.'''
#         
#         # Create a popped out window
#         poppedout_window = PoppedOutWindow(self.root)
#         poppedout_id = poppedout_window.poppedout_id
#         
#         # Add a fitting title to the popped out window
#         if axis_id==0:
#             y_quantity = 'Heat flux'
#             if y_quantity=='Heat flux':       poppedout_window.set_title("Stellapy: Heat flux")
#             if y_quantity=='Momentum flux':   poppedout_window.set_title("Stellapy: Momentum flux")
#             if y_quantity=='Particle flux':   poppedout_window.set_title("Stellapy: Particle flux")
#             # ...
# 
#         # Initiate the plotting class
#         self.initiate_GraphClass(figure=poppedout_window.figure)
#         
#         # Initiate the canvas        
#         CanvasForGraphs(self.root, poppedout_window.frame, self.root.tab_Linear, poppedout_window.figure, poppedout_id=poppedout_id)
#         
#         # Parse the data from the main canvas to the popped out canvas
#         self.root.graph_poppedOut[poppedout_id].range = self.tab.Graph[axis_id].range
#         self.root.graph_poppedOut[poppedout_id].label = self.tab.Graph[axis_id].label
#         
#         # Now plot the figure on the popped out canvas
#         if axis_id==0: Plot = self.tab.PlotVariablesVsZ 
#         
#         # Replot
#         if axis_id in [0,1,2]: Plot.replot = True; Plot.plot_graph(poppedout_id)
#         if axis_id == 3:       Plot.replot = True; Plot.plot_graph(poppedout_id, "kx")
#         if axis_id == 4:       Plot.replot = True; Plot.plot_graph(poppedout_id, "ky")
#         
#         # Update screen
#         self.root.canvasPoppedOut[poppedout_id].draw_idle()
#         return 
#     
#     #------------------------------------------------
#     def reset_graph(self, axis_id, poppedout_id=None):
#         ''' Method called by the "Reset graph" button on the canvas. '''
#         
#         # Perhaps the button was clicked on the popped out window
#         if poppedout_id == None: Graph = self.tab.Graph[axis_id]
#         if poppedout_id != None: Graph = self.root.graph_poppedOut[poppedout_id]
#         
#         # Reset the axis
#         Graph.load_defaults()
#         
#         # Now plot the graph 
#         if axis_id==0: Plot = self.tab.PlotVariablesVsZ 
#         Plot.replot = True; Plot.plot_graph(None)
#         return 

