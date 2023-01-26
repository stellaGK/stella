

import tkinter as tk
from tkinter import ttk
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk  
matplotlib.rcParams["toolbar"] = "toolbar2"
from stellapy.utils.config import CONFIG   
from stellapy.GUI.widgets.ToolTip import ToolTip  

#################################################################
#                       CUSTOM TOOLBAR
#################################################################

class CustomToolbar(NavigationToolbar2Tk): 
    """ Customized Navigator object (Toolbar on matplotlib canvas).
    - Removed mouse_move event.x and event.y display
    - Changed buttons layout
    """ 

    def __init__(self, root, Canvas, FigureCanvas, frame_canvasAndToolbar, direction="x"):

        # Safe the root, canvas and direction of the toolbar
        self.root = root
        self.Canvas = Canvas
        self.direction = direction
        
        # Safe the custom buttons
        self.customButtons = []
        
        # Initiate the standard toolbar (which is empty) and manually add it to the canvas so we can use grid()   
        self.toolitems = ()
        NavigationToolbar2Tk.__init__(self,FigureCanvas,frame_canvasAndToolbar,pack_toolbar=False)
        if self.direction=="x": self.grid(row=1, column=0, stick='NSEW')
        if self.direction=="y": self.grid(row=0, column=0, stick='NSEW')

        # There is a tk.Label object at the end of the toolbar which has the wrong color so remove it
        # See /usr/local/lib/python3.6/dist-packages/matplotlib/backends/_backend_tk.py line 521
        list_widgets = self.pack_slaves()
        for i in range(len(list_widgets)):
            list_widgets[i].destroy()
        
        # Define our own custom buttons
        self.tool_buttons = []; num = 0
        self.extratoolitems = (
            ('Save', 'Save the figure', 'icon', 'custom_save_figure'),
            ('Home', 'Reset original view', 'icon', 'home'),
            ('Back', 'Back to previous view', 'icon', 'back'),
            ('Forward', 'Forward to next view', 'icon', 'forward'),
            ('Pan', 'Pan axes with left mouse, zoom with right', 'icon', 'pan'),
            ('Zoom', 'Zoom to rectangle', 'icon', 'zoom'),
            (None, None, None, None),
            ('Reset', 'Reset the figure', 'icon', 'reset'),
            ('Options', 'Figure options', 'icon', 'optionsWindow'),
            ('PopOut', 'Open the figure in a new window', 'icon', 'popout_window'))
        
        # If the time is dark, use white icons instead of black icons
        extension = "_inv" if self.root.theme=="awdark" else ""
        self.toolbar_icons = ["filesave","home","back","forward","move","zoom_to_rect","reload","subplots","qt4_editor_options"]
        self.toolbar_icons = [ CONFIG['CODE']['Stellapy']+"GUI/images/"+icon+extension+".png" for icon in self.toolbar_icons ]
        self.toolbar_icons.insert(6, None)
        self.toolbar_icons.append(None)
        
        # Add the icons to the toolbar
        for text, tooltip_text, image_file, callback in self.extratoolitems: #@unusedvariable
            if text is None:
                self.add_customLabelSpacing(num)
            else:
                try:
                    button = self.add_customButtonWithIcon(text=text, file=self.toolbar_icons[num], command=getattr(self, callback), i=num)
                    if tooltip_text is not None:
                        # Custom tooltip class with a half second delay
                        ToolTip(self.root, button, tooltip_text)
                except IndexError:
                    pass
            num+=1
    
        # Configure the grid: put fillers before and after the buttons so its centered
        if self.direction=="x": 
            tk.Grid.rowconfigure(self, 0, weight=0) 
            tk.Grid.columnconfigure(self, 0, weight=1)
            tk.Grid.columnconfigure(self, num+1, weight=1)
        if self.direction=="y": 
            tk.Grid.columnconfigure(self, 0, weight=0) 
            tk.Grid.rowconfigure(self, 0, weight=1)
            tk.Grid.rowconfigure(self, num+1, weight=1)
            tk.Grid.rowconfigure(self, num+2, weight=1)

################################################################################
#                            METHODS FOR THE TOOLBAR 
################################################################################

    def add_customButtonWithIcon(self, text, file, command, i):
        im = tk.PhotoImage(file=file)
        self.customButtons.append(ttk.Button(master=self, text=text, image=im, command=command, style='toolbar.TButton'))
        self.customButtons[-1]._ntimage = im # Attach image to button, otherwise it disappears
        if self.direction=="x": self.customButtons[-1].grid(row=0, column=i+1)
        if self.direction=="y": self.customButtons[-1].grid(row=i+1, column=0)
        return self.customButtons[-1]
    
    #----------------------
    def add_customLabelSpacing(self, i):
        label = ttk.Label(master=self, text="   ")
        label.grid(row=0, column=i+1)
        self.tool_buttons.append(label)
        return label

    #----------------------     
    def add_customButton(self, text, command, **kwargs):
        button = ttk.Button(master=self, text=text, command=command, style='toolbar.TButton', **kwargs)
        button.pack(side=tk.LEFT,fill="y", padx=5)
        return button
    
################################################################################
#                               BUTTON METHODS 
################################################################################

    def custom_save_figure(self):
        ''' Before calling save_figure, set the background to white. '''
        # Set the figure, axis and legend to white
        self.canvas.figure.patch.set_facecolor('white')  
        self.Canvas.Axis.ax.patch.set_facecolor('white')  
        self.makeLegendWhite()
        self.save_figure()
        # Return the original colors
        self.canvas.figure.patch.set_facecolor(self.root.color['canvas'])  
        self.Canvas.Axis.ax.patch.set_facecolor(self.root.color['canvas'])  
#         if legend!=None: 
#             frame = legend.get_frame(); 
#             frame.set_color(self.root.color['canvas'])
#             frame.set_edgecolor(self.root.color['fg'])
        return
    
    #----------------------  
    def reset(self):
        ''' Replot the figure with its default ranges and labels. '''
        self.Canvas.reset_graph()
        return

    #----------------------  
    def optionsWindow(self):
        ''' Open a window with options to modify the axis/ranges/labels. '''
        self.Canvas.open_optionsWindow()
        return 
    
    #----------------------  
    def popout_window(self):
        ''' Replot the figure in a new window. '''
        self.Canvas.popout_window() 
        return 

    def makeLegendWhite(self, **kwargs):
        
        l = self.Canvas.Axis.ax.legend_
        if l!=None:
            
            import matplotlib as mpl
            defaults = dict(
                loc = l._loc,
                numpoints = l.numpoints,
                markerscale = l.markerscale,
                scatterpoints = l.scatterpoints,
                scatteryoffsets = l._scatteryoffsets,
                prop = l.prop,
                borderpad = l.borderpad,
                labelspacing = l.labelspacing,
                handlelength = l.handlelength,
                handleheight = l.handleheight,
                handletextpad = l.handletextpad,
                borderaxespad = l.borderaxespad,
                columnspacing = l.columnspacing,
                ncol = l._ncol, 
                mode = l._mode,
                fancybox = type(l.legendPatch.get_boxstyle())==mpl.patches.BoxStyle.Round,
                shadow = l.shadow,
                title = l.get_title().get_text() if l._legend_title_box.get_visible() else None,
                framealpha = l.get_frame().get_alpha(),
                bbox_to_anchor = l.get_bbox_to_anchor()._bbox, 
                bbox_transform = l.get_bbox_to_anchor()._transform, 
                handler_map = l._custom_handler_map)
     
            changed_values = {"facecolor" : "white"}
            self.Canvas.Axis.ax.legend(**dict(list(defaults.items()) + list(changed_values.items())))
