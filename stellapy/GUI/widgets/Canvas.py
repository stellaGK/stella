
import tkinter as tk 
import matplotlib; matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvas
from stellapy.GUI.widgets.CustomToolbar import CustomToolbar 
from stellapy.GUI.widgets.PoppedOutWindow import PoppedOutWindow
from stellapy.GUI.widgets.Axis import Axis
from stellapy.GUI.widgets.OptionsWindow import OptionsWindow

################################################################################
#                          CLASS TO CREATE CANVASSES
################################################################################
    
class Canvas: 
    ''' The class <Canvas> will create a ttk.frame that holds a <FigureCanvas>
    and a modified version of <NavigationToolbar2Tk> from the matplotlib package. 
    Together they allow one to display plots in the GUI. '''
    
    def __init__(self, root, frame, identifier, direction="x"):
        ''' For each canvas+toolbar embedded in the GUI this class is called. The
        <root> of the GUI is required and the <frame> where the canvas+toolbar will
        be embedded. The <direction> determines the orientation of the toolbar. 
        The <graph> class holds all information needed for the custom toolbar and
        the options window, and to set the titles of the frames. '''
        
        # Save the variables
        self.root = root
        self.frame = frame  
        self.identifier = identifier 
        self.direction = direction
    
        # Create an <Axis> object that holds the matplotlib <figure> that will 
        # be displayed in the <Canvas>, and create an <ax> object for the figure
        self.Axis = Axis(self, root, identifier)
        
        # Configure the frame and add a canvas and toolbar to it
        self.configure_frameCanvasAndToolbar()
        self.add_canvas()
        self.add_toolbar()
        
        # Put the empty figure on the canvas so its visible in the GUI
        self.draw_idle()
        self.root.update_idletasks()        
                                                           
    #----------- Frame ----------
    def configure_frameCanvasAndToolbar(self):
        ''' Create the frame.'''
        
        # Configure the frame that holds the canvas and toolbar
        if self.direction=="x":
            tk.Grid.rowconfigure(   self.frame, 0, weight=1) 
            tk.Grid.rowconfigure(   self.frame, 1, weight=0) 
            tk.Grid.columnconfigure(self.frame, 0, weight=1)  
        if self.direction=="y":
            tk.Grid.columnconfigure(self.frame, 0, weight=0) 
            tk.Grid.columnconfigure(self.frame, 1, weight=1) 
            tk.Grid.rowconfigure(   self.frame, 0, weight=1)  
        
    #----------- Canvas -----------
    def add_canvas(self):
        ''' Add a canvas. '''

        # Add the figure to a Canvas
        self.FigureCanvas = FigureCanvas(self.Axis.figure, self.frame) 
        self.FigureCanvas.get_tk_widget().configure(background=self.root.color['canvas'])
        self.FigureCanvas.draw()
        self.widget = self.FigureCanvas.get_tk_widget() 
        if self.direction=="x": self.widget.grid(row=0, column=0, stick='NSEW')
        if self.direction=="y": self.widget.grid(row=0, column=1, stick='NSEW')

    #----------- Toolbar -----------
    def add_toolbar(self):
        ''' Add a toolbar. '''
        
        # Add a toolbar to the canvas
        self.toolbar = CustomToolbar(self.root, self, self.FigureCanvas, self.frame, self.direction)
        self.toolbar.config(background=self.root.color['bg'])
        self.toolbar.update()

################################################################################
#                      METHODS FOR THE TOOLBAR BUTTONS
################################################################################

    def add_plottingClass(self, plottingClass):
        ''' To each canvas, a plotting function needs to be linked. '''
        self.plottingClass = plottingClass
        return
    
    #----------------------------------- 
    def draw_idle(self):
        ''' Draw idle on the actual <FigureCanvas> object. '''
        self.FigureCanvas.draw_idle()
    
    #----------------------------------- 
    def plot(self):
        ''' Plot the plotting function on this canvas. '''
        self.plottingClass.plot(Canvas)
    
    #----------------------------------- 
    def popout_window(self):
        '''Replot the figure in a seperate window when the "popout" button on
        the toolbar is clicked. Attach this function to the <Canvas> object of
        the GUI, since it is called from the button on this specific <Canvas>.'''
        
        # Create a popped out window and its id in <self.root.canvasPoppedOut>
        poppedout_window = PoppedOutWindow(self.root)
        identifier = "poppepoutwindow_"+str(poppedout_window.poppedout_id)

        # Create a <Canvas>, <figure> and <ax> in the <frame>        
        PoppedOutCanvas = Canvas(self.root, poppedout_window.frame, identifier, self.direction)
        self.root.canvasPoppedOut.append(PoppedOutCanvas)
        
        # When creating a popped out window, we want an exact copy of the 
        # original figure, therefore copy all the axis and optionswindow data.
        PoppedOutCanvas.Axis.copy_dataFromOtherAxis(self.Axis)
        PoppedOutCanvas.plottingClass = self.plottingClass
        
        # Change the margins to better fit a popped out window
        PoppedOutCanvas.Axis.update_gridSpeficications(top=0.92, left=0.1, right=0.95, bottom=0.15)
        
        # Add a fitting title to the popped out window
        poppedout_window.set_title("Stellapy: "+self.Axis.get_standardHeader()) 
        
        # Now plot the figure on the popped out canvas
        PoppedOutCanvas.plottingClass.plot(PoppedOutCanvas)
                
    #----------------------------------- 
    def reset_graph(self): 
        ''' Reset the axis and then plot again. '''
        self.Axis.reset_axis() 
        self.plottingClass.plot(self)
        return
        
    #----------------------------------- 
    def open_optionsWindow(self): 
        ''' Open a window with options to modify the axis/ranges/labels. '''
        self.OptionsWindow = OptionsWindow(self.root, self)  
        self.OptionsWindow.open_optionsWindow()
        return
        
