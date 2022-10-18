
from stellapy.GUI.widgets import Canvas 
from stellapy.GUI.plot.linearSimulations.plot_quantityVsTime import plot_quantityVsTime 

################################################################################
#       CLASS TO PLOT THE TIME EVOLUTION OF MODES OF LINEAR SIMULATIONS
################################################################################
            
class PlotTimeEvolution: 
    ''' Holds and manipulates the Canvas and Axis object for the following plots:
            (1) Omega versus time
            (2) Gamma versus time
            (3) Phi**2 versus time  '''

    def __init__(self, tab):

        # Set an identifier for this plotting class, used as the figure name
        self.identifier = "TabConvergence:PlotTimeEvolution"
        
        # Make the parent objects available
        self.tab = tab
        self.root = tab.root  
        
        # Remember whether we have a Canvas yet
        self.initiated_canvas = False
        
    #---------------------- 
    def initiate_canvas(self, frame):
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
            
            # Set the (x,y) quantities and titles for the canvas
            self.Axis.set_axisQuantities(1, "t", "omega")
            self.Axis.set_axisQuantities(2, "t", "gamma")
            self.Axis.set_axisQuantities(3, "t", "phi2") 
        
            # Initiate the plot that is currently set
            self.Axis.set_plot(self.tab.Simulations.var_plot.get())
            
            # Remember that the <Canvas> is initiaited 
            self.initiated_canvas = True
            
            # Update the kx/ky modes
            self.tab.Simulations.display_plottedModes()

        return

    #---------------------- 
    def get_argumentsForPlot(self, Canvas):
        ''' Gather the arguments for <plot_quantityVsTime>. '''
        plotting_arguments = {
            # Specify which simulations to plot
                "research"      : self.root.Research,\
                "experiment_id" : self.tab.Simulations.experiment_id,\
                "simulation_id" : self.tab.Simulations.simulation_id,\
            # Data on the (x,y) axis  
                "y_quantity"    : Canvas.Axis.y_quantity,\
                "x_range"       : Canvas.Axis.x_range,\
                "y_range"       : Canvas.Axis.y_range,\
                "units"         : Canvas.Axis.units,\
            # Modes
                "modes_id" : self.tab.Simulations.modes_id,\
                "kx_range"      : Canvas.Axis.kx_range,\
                "ky_range"      : Canvas.Axis.ky_range,\
            # Labels
                "x_label"       : Canvas.Axis.x_label,\
                "y_label"       : Canvas.Axis.y_label,\
                "title"         : Canvas.Axis.title,\
            # Figue-re
                "ax"            : Canvas.Axis.ax,\
                "Progress"      : self.tab.Progress,\
            # Legend
                "fontsize"     : Canvas.Axis.fontsize,\
                "handlelength" : Canvas.Axis.handlelength}
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
                plot_quantityVsTime(**plotting_arguments)
                Canvas.Axis.finish_plotting(research_arguments, input_files)
                Canvas.Axis.save_plottingArguments(self.get_argumentsForPlot(Canvas))
        
                # Update the non-converged simulations
                self.tab.Simulations.write_informationConvergence()           
            
        # If no simulations have been selected, clear the current figure
        if self.root.Research.input_files == []:
            Canvas.Axis.reset_axis()
            self.tab.Simulations.txt_convergence2.set("Please select simulations.")
        return

    



 





