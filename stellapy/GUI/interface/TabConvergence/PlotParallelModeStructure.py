
from stellapy.GUI.widgets import Canvas 
from stellapy.GUI.plot.linearSimulations.plot_quantityVsZ import plot_quantityVsZ

################################################################################
#       CLASS TO PLOT THE PARALLEL MODE STRUCTURE OF LINEAR SIMULATIONS
################################################################################
            
class PlotParallelModeStructure: 
    ''' Holds and manipulates the Canvas and Axis object for the following plots:
            (1) Re(Phi) versus zeta
            (2) Im(Phi) versus zeta
            (3) Phi**2 versus zeta  '''

    def __init__(self, tab):
        ''' Initiate the class that plots the parallel mode structure. '''
        
        # Set an identifier for this plotting class, used as the figure name
        self.identifier = "TabConvergence:PlotParallelModeStructure"
        
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

            # Resize the <matplotlib figure>
            self.Axis.update_gridSpeficications(top=0.95, left=0.1, right=0.95, bottom=0.1)
            
            # Set the (x,y) quantities and titles for the canvas
            self.Axis.set_axisQuantities(1, "zeta", "phi_real")
            self.Axis.set_axisQuantities(2, "zeta", "phi_imag")
            self.Axis.set_axisQuantities(3, "zeta", "phi2") 
            
            # Initiate the plot that is currently set
            self.Canvas.Axis.set_plot(self.tab.Simulations.var_plot.get()-3)

            # Since this canvas is initiated after that of PlotTimeEvolution, 
            # copy the (kx,ky) selection from the PlotTimeEvolution Canvas
            self.Canvas.Axis.kx_range = self.tab.PlotTimeEvolution.Canvas.Axis.kx_range
            self.Canvas.Axis.ky_range = self.tab.PlotTimeEvolution.Canvas.Axis.ky_range 
            
            # Remember that the <Canvas> is initiaited
            self.initiated_canvas = True
            
        return
    
    #---------------------- 
    def get_argumentsForPlot(self, Canvas):
        ''' Gather the arguments for <plot_quantityVsZ>. '''
        plotting_arguments = {
            # Specify which simulations to plot
                "research"      : self.root.Research,\
                "experiment_id" : self.tab.Simulations.experiment_id,\
                "simulation_id" : self.tab.Simulations.simulation_id,\
            # Data on the (x,y) axis  
                "x_quantity"    : Canvas.Axis.x_quantity,\
                "y_quantity"    : Canvas.Axis.y_quantity,\
                "x_range"       : Canvas.Axis.x_range,\
                "y_range"       : Canvas.Axis.y_range,\
                "units"         : Canvas.Axis.units,\
            # Modes
                "modes_id"      : self.tab.Simulations.modes_id,\
                "kx_range"      : Canvas.Axis.kx_range,\
                "ky_range"      : Canvas.Axis.ky_range,\
            # Labels
                "x_label"       : Canvas.Axis.x_label,\
                "y_label"       : Canvas.Axis.y_label,\
                "title"         : Canvas.Axis.title,\
            # Figure
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
                
                # Update the non-converged simulations
                self.tab.Simulations.write_informationConvergence()           
            
        # If no simulations have been selected, clear the current figure
        if input_files == []:
            Canvas.Axis.reset_axis()
            self.tab.Simulations.txt_convergence2.set("Please select simulations.")
        return

    



 





