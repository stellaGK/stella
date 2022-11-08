

import tkinter as tk
from tkinter import ttk
from stellapy.GUI.widgets import PAD_TITLE, PAD_LABEL, PAD_ENTRY  

################################################################################
#                OPEN A WINDOW WITH OPTIONS TO EDIT THE AXIS
################################################################################

class OptionsWindow:
    ''' Open a window with options to modify the axis/ranges/labels. '''

    def __init__(self, root, Canvas):
        
        # Save the root, Canvas and Axis object
        self.root = root
        self.Canvas = Canvas 
        
        # Get the width and height of the window
        self.height = root.winfo_height() 
        self.width =  root.winfo_width()  
        return
        
    #-------------------------
    def open_optionsWindow(self):
        
        # Create the options window, keep it on top and give it a good size
        self.window_options = tk.Toplevel(self.root)
        self.window_options.title("Modify the figure") 
        self.window_options.attributes('-topmost', 'true')
        self.rescale_window()
        
        # Create a tabbed view with the possible settings
        self.tab_header = ttk.Notebook(self.window_options, style='header.TNotebook')
        
        # Add frames to the tab_header which are the tab windows
        self.tab_axis   = ttk.Frame(self.tab_header) 
        self.tab_curves = ttk.Frame(self.tab_header)  
        self.tab_legend = ttk.Frame(self.tab_header)   
        
        # Add the tabs to the tab header
        self.tab_header.add(self.tab_axis,   text='Axis')
        self.tab_header.add(self.tab_curves, text='Curves') 
        self.tab_header.add(self.tab_legend, text='Legend') 
        self.tab_header.pack(expand=1, fill='both')
        
        # Fill the tabs with widgets through classes
        self.tabAxes    = tabAxis(  self.root, self.Canvas, self.tab_axis)
        self.tabCurves  = tabCurves(self.root, self.Canvas, self.tab_curves)
        self.tabLegend  = tabLegend(self.root, self.Canvas, self.tab_legend)
        return

    #-------------------------
    def rescale_window (self):    
        ''' Make sure the size of the options window is decent. '''
        
        # Standard size of the screen
        winx = 200
        winy = 600
        
        # Adjust the height if we have a third axis
        if self.Canvas.Axis.z_quantity!=None:
            winy = 750
                
        # Center the new window in the screen
        x = self.width/2  - winx/2
        y = self.height/2 - winy/2
        self.window_options.geometry("+%d+%d" % (x, y))        

################################################################################
#                                AXIS OPTIONS
################################################################################
   
class tabAxis:
    
    def __init__(self, root, Canvas, tab):
        
        # Get data from the GUI
        self.root    = root 
        self.Canvas  = Canvas
        self.Axis    = Canvas.Axis
        self.tab     = tab  
        self.row     = 0

        # Add some variables
        self.options_scale = sorted(("Linear", "Logaritmic"))
        self.options_units = sorted(("Normalized", "SI units"))
        
        # Create the frame 
        self.frame = ttk.Frame(tab)
        self.frame.pack(expand=1, fill=tk.BOTH)
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.frame, 0, weight=0) # x title
        tk.Grid.rowconfigure(   self.frame, 1, weight=0) # x min
        tk.Grid.rowconfigure(   self.frame, 2, weight=0) # x max
        tk.Grid.rowconfigure(   self.frame, 3, weight=0) # x label
        tk.Grid.rowconfigure(   self.frame, 4, weight=0) # x scale
        tk.Grid.rowconfigure(   self.frame, 5, weight=0) # y title
        tk.Grid.rowconfigure(   self.frame, 6, weight=0) # y min
        tk.Grid.rowconfigure(   self.frame, 7, weight=0) # y max
        tk.Grid.rowconfigure(   self.frame, 8, weight=0) # y label
        tk.Grid.rowconfigure(   self.frame, 9, weight=0) # y scale 
        tk.Grid.columnconfigure(self.frame, 0, weight=1, uniform="options") 
        tk.Grid.columnconfigure(self.frame, 0, weight=1, uniform="options") 

        # Add elements to the frame
        self.init_options()
        self.init_xAxis()
        self.init_yAxis() 
        if self.Axis.z_quantity!=None:
            self.init_zAxis()
        return 
    
    #-------------------------
    def init_options(self):
        ''' Change the title and units of the Axis. '''
        
        def update_title(event):
            ''' Change the title of the figure. '''
            try: fontsize = int(self.var_font.get())
            except: fontsize = str(self.var_font.get())
            self.Axis.title = self.var_title.get() 
            self.Axis.ax.set_title(self.var_title.get(), fontsize=fontsize) 
            self.Canvas.draw_idle()
        
        def update_units(*args):
            ''' Change the units of the figure, replot the figure immediatly, to
            set the (x,y) ranges correctly, since they will change. '''
            if self.var_units.get()=="Normalized":  
                if self.Axis.units!="normalized":
                    self.Axis.units = "normalized"
                    self.Canvas.plot()
                    number_format = "{:.2}" 
            if self.var_units.get()=="SI units":    
                if self.Axis.units!="SI units":
                    self.Axis.units = "SI units"
                    self.Canvas.plot()
                    number_format = "{:.2e}"
            self.var_xMin.set(number_format.format(self.Axis.x_range[0],2))
            self.var_xMax.set(number_format.format(self.Axis.x_range[1],2))
            self.var_yMin.set(number_format.format(self.Axis.y_range[0],2))
            self.var_yMax.set(number_format.format(self.Axis.y_range[1],2))
            self.var_xlabel.set(self.Axis.x_label)
            self.var_ylabel.set(self.Axis.y_label)
            
        def update_fontsize(*_):
            try: fontsize = int(self.var_font.get())
            except: fontsize = str(self.var_font.get()) 
            self.Axis.ax.set_title(self.var_title.get(), fontsize=fontsize) 
            self.Canvas.draw_idle()
                    
        # Change the Axis title
        self.lbl_Axis = ttk.Label(self.frame, text="Title", style='prefTitle.TLabel')
        self.var_title = tk.StringVar(value=self.Axis.title)
        self.lbl_title = ttk.Label(self.frame, text="Title")
        self.ent_title = ttk.Entry(self.frame, textvariable=self.var_title, width=20, style='opt_valueR.TEntry')
        self.ent_title.bind('<Return>', update_title)

        # Choose the font size of the title
        self.var_font = tk.StringVar(value=str(20))
        self.lbl_font = ttk.Label(self.frame, text="Font size: ")
        self.ent_font = ttk.Entry(self.frame, textvariable=self.var_font, width=20, style='opt_valueR.TEntry')
        self.ent_font.bind('<Return>', update_fontsize) 
        
        # Change the units
        self.var_units = tk.StringVar(value=self.options_units[0])
        self.lbl_units = ttk.Label(self.frame, text="Units")
        self.mnu_units = ttk.OptionMenu(self.frame, self.var_units, self.options_units[0], *self.options_units, style='option.TMenubutton')
        self.mnu_units["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'])
        if self.Axis.units=="normalized": self.var_units.set(self.options_units[0])
        if self.Axis.units=="SI units":   self.var_units.set(self.options_units[1])
        self.var_units.trace('w', update_units) # link function to a change of the dropdown options
    
        # Add the widgets to the frame
        self.lbl_Axis.grid( row=self.row, column=0, columnspan=2, **PAD_TITLE); self.row += 1
        self.lbl_title.grid(row=self.row, column=0, **PAD_LABEL)
        self.ent_title.grid(row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_font.grid( row=self.row, column=0, **PAD_LABEL)
        self.ent_font.grid( row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_units.grid(row=self.row, column=0, **PAD_LABEL)
        self.mnu_units.grid(row=self.row, column=1, **PAD_ENTRY); self.row += 1
    
    #-------------------------
    def init_xAxis(self):
        ''' Modify the x-axis of the Axis. '''
        
        def update_xAxis(event):
            new_range = [float(self.var_xMin.get()), float(self.var_xMax.get())]
            self.Axis.x_range = new_range
            self.Axis.ax.set_xlim(new_range)
            self.Canvas.draw_idle()
            
        def update_xLabel(event):
            self.Axis.x_label = self.var_xlabel.get()
            self.Axis.ax.set_xlabel(self.Axis.x_label)
            self.Canvas.draw_idle()
            
        def update_xscale(*args):
            if self.var_xscale.get()=="Linear":     scale='linear'
            if self.var_xscale.get()=="Logaritmic": scale='log'
            self.Axis.x_scale = scale
            self.Axis.ax.set_xscale(scale)
            # Logaritmic axis needs a positive start
            if float(self.var_xMin.get()) <= 0: 
                if float(self.var_xMin.get()) > 1:       self.var_xMin.set(1) 
                elif float(self.var_xMin.get()) > 0.1:   self.var_xMin.set(0.1) 
                elif float(self.var_xMin.get()) > 0.01:  self.var_xMin.set(0.01) 
                else:                                    self.var_xMin.set(0.001) 
                new_range = [float(self.var_xMin.get()), float(self.var_xMax.get())]
                self.Axis.x_range = new_range
                self.Axis.ax.set_xlim(new_range)
            self.Canvas.draw_idle()
        
        # Minimum and maximum of the x-axis
        x_range = self.Axis.x_range if (self.Axis.x_range!=None) else [0,0] 
        self.lbl_xTitle = ttk.Label(self.frame, text="x-axis: "+self.Axis.x_name, style='prefTitle.TLabel')
        self.lbl_xMin = ttk.Label(self.frame, text="Minimum")
        self.lbl_xMax = ttk.Label(self.frame, text="Maximum")
        self.var_xMin = tk.StringVar(value=round(x_range[0],2))
        self.var_xMax = tk.StringVar(value=round(x_range[1],2))
        self.ent_xMin = ttk.Entry(self.frame, textvariable=self.var_xMin, width=5, style='opt_valueR.TEntry')
        self.ent_xMax = ttk.Entry(self.frame, textvariable=self.var_xMax, width=5, style='opt_valueR.TEntry')
        self.ent_xMin.bind('<Return>', update_xAxis)
        self.ent_xMax.bind('<Return>', update_xAxis)
        
        # Label for the x-axis
        self.var_xlabel = tk.StringVar(value=self.Axis.x_label)
        self.lbl_xlabel = ttk.Label(self.frame, text="Label")
        self.ent_xlabel = ttk.Entry(self.frame, textvariable=self.var_xlabel, width=20, style='opt_valueR.TEntry')
        self.ent_xlabel.bind('<Return>', update_xLabel)
                
        # Choice between linear and log scales for the x-axis 
        self.var_xscale = tk.StringVar(value=self.options_scale[0])
        self.lbl_xscale = ttk.Label(self.frame, text="Scale")
        self.mnu_xscale = ttk.OptionMenu(self.frame, self.var_xscale, self.options_scale[0], *self.options_scale, style='option.TMenubutton')
        self.mnu_xscale["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'])
        if self.Axis.x_scale=="linear": self.var_xscale.set(self.options_scale[0])
        if self.Axis.x_scale=="log":    self.var_xscale.set(self.options_scale[1])
        self.var_xscale.trace('w', update_xscale) # link function to a change of the dropdown options
    
        # Add the labels to the frame 
        self.lbl_xTitle.grid( row=self.row, column=0, columnspan=2, **PAD_TITLE); self.row += 1
        self.lbl_xMin.grid(   row=self.row, column=0, **PAD_LABEL)
        self.ent_xMin.grid(   row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_xMax.grid(   row=self.row, column=0, **PAD_LABEL)
        self.ent_xMax.grid(   row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_xlabel.grid( row=self.row, column=0, **PAD_LABEL)
        self.ent_xlabel.grid( row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_xscale.grid( row=self.row, column=0, **PAD_LABEL)
        self.mnu_xscale.grid( row=self.row, column=1, **PAD_ENTRY); self.row += 1

    #-------------------------
    def init_yAxis(self):
        ''' Modify the y-axis of the Axis. '''

        def update_yAxis(event):
            new_range = [float(self.var_yMin.get()), float(self.var_yMax.get())]
            self.Axis.y_range = new_range
            self.Axis.ax.set_ylim(new_range)
            self.Canvas.draw_idle()
            
        def update_yLabel(event):
            self.Axis.y_label = self.var_ylabel.get()
            self.Axis.ax.set_ylabel(self.Axis.y_label)
            self.Canvas.draw_idle()
            
        def update_yscale(*args):
            if self.var_yscale.get()=="Linear":     scale='linear'
            if self.var_yscale.get()=="Logaritmic": scale='log'
            self.Axis.y_scale = scale
            self.Axis.ax.set_yscale(scale)
            # Logaritmic axis needs a positive start
            if float(self.var_yMin.get()) <= 0: 
                if float(self.var_yMin.get()) > 1:      self.var_yMin.set(1) 
                elif float(self.var_yMin.get()) > 0.1:  self.var_yMin.set(0.1) 
                elif float(self.var_yMin.get()) > 0.01: self.var_yMin.set(0.01)
                else:                                   self.var_yMin.set(0.001)
                new_range = [float(self.var_yMin.get()), float(self.var_yMax.get())]
                self.Axis.y_range = new_range
                self.Axis.ax.set_ylim(new_range)
            self.Canvas.draw_idle()             
        
        # Minimum and maximum of the x-axis 
        y_range = self.Axis.y_range if (self.Axis.y_range!=None) else [0,0] 
        self.lbl_yTitle = ttk.Label(self.frame, text="y-axis: "+self.Axis.y_name, style='prefTitle.TLabel')
        self.lbl_space = ttk.Label(self.frame, text="    ", style='prefTitle.TLabel')
        self.lbl_yMin = ttk.Label(self.frame, text="Minimum")
        self.lbl_yMax = ttk.Label(self.frame, text="Maximum")
        self.var_yMin = tk.StringVar(value=round(y_range[0],2))
        self.var_yMax = tk.StringVar(value=round(y_range[1],2))
        self.ent_yMin = ttk.Entry(self.frame, textvariable=self.var_yMin, width=5, style='opt_valueR.TEntry')
        self.ent_yMax = ttk.Entry(self.frame, textvariable=self.var_yMax, width=5, style='opt_valueR.TEntry')
        self.ent_yMin.bind('<Return>', update_yAxis)
        self.ent_yMax.bind('<Return>', update_yAxis)
        
        # Label for the y-axis
        self.var_ylabel = tk.StringVar(value=self.Axis.y_label)
        self.lbl_ylabel = ttk.Label(self.frame, text="Label")
        self.ent_ylabel = ttk.Entry(self.frame, textvariable=self.var_ylabel, width=20, style='opt_valueR.TEntry')
        self.ent_ylabel.bind('<Return>', update_yLabel)
        
        # Choice between linear and log scales for the x-axis 
        self.var_yscale = tk.StringVar(value=self.options_scale[0])
        self.lbl_yscale = ttk.Label(self.frame, text="Scale")
        self.mnu_yscale = ttk.OptionMenu(self.frame, self.var_yscale, self.options_scale[0], *self.options_scale, style='option.TMenubutton')
        self.mnu_yscale["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'])
        if self.Axis.y_scale=="linear": self.var_yscale.set(self.options_scale[0])
        if self.Axis.y_scale=="log":    self.var_yscale.set(self.options_scale[1])
        self.var_yscale.trace('w', update_yscale) # link function to a change of the dropdown options
    
        # Add the labels to the frame 
        self.lbl_yTitle.grid( row=self.row, column=0, columnspan=2, **PAD_TITLE); self.row += 1
        self.lbl_yMin.grid(   row=self.row, column=0, **PAD_LABEL)
        self.ent_yMin.grid(   row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_yMax.grid(   row=self.row, column=0, **PAD_LABEL)
        self.ent_yMax.grid(   row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_ylabel.grid( row=self.row, column=0, **PAD_LABEL)
        self.ent_ylabel.grid( row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_yscale.grid( row=self.row, column=0, **PAD_LABEL)
        self.mnu_yscale.grid( row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_space.grid(  row=self.row, column=0, columnspan=2, **PAD_TITLE); self.row += 1
        return
    
    #------------------------- 
    def init_zAxis(self):
        ''' Modify the z-axis of the Axis. '''

        def update_zAxis(event):
            new_range = [float(self.var_zMin.get()), float(self.var_zMax.get())]
            self.Axis.z_range = new_range
            self.Axis.ax.clim(new_range)
            self.Canvas.draw_idle()
            
        def update_zLabel(event):
            self.Axis.label["z"] = self.var_zlabel.get()
            self.Canvas.draw_idle() 
            
        def update_zscale(*args):
            if self.var_zscale.get()=="Linear":     scale='linear'
            if self.var_zscale.get()=="Logaritmic": scale='log'
            self.Axis.range["z_scale"] = scale
            self.Axis.ax_twin.set_yscale(scale)
            # Logaritmic axis needs a positive start
            if float(self.var_zMin.get()) <= 0: 
                if float(self.var_zMax.get()) > 10:    self.var_zMin.set(1) 
                elif float(self.var_zMax.get()) > 1:   self.var_zMin.set(0.1) 
                elif float(self.var_zMax.get()) > 0.1: self.var_zMin.set(0.01)
                range_ = [float(self.var_zMin.get()), float(self.var_zMax.get())]
                self.Axis.range["z"] = range_
                self.Axis.ax_twin.set_ylim(range_)
            self.Canvas.draw_idle()             
        
        # Minimum and maximum of the twinned y-axis
        self.lbl_zTitle = ttk.Label(self.frame, text="y-axis: "+self.z_name, style='prefTitle.TLabel')
        self.lbl_space = ttk.Label(self.frame, text="    ", style='prefTitle.TLabel')
        self.lbl_zMin = ttk.Label(self.frame, text="Minimum")
        self.lbl_zMax = ttk.Label(self.frame, text="Maximum")
        self.var_zMin = tk.StringVar(value=round(self.Axis.range["z"][0],2))
        self.var_zMax = tk.StringVar(value=round(self.Axis.range["z"][1],2))
        self.ent_zMin = ttk.Entry(self.frame, textvariable=self.var_zMin, width=5, style='opt_valueR.TEntry')
        self.ent_zMax = ttk.Entry(self.frame, textvariable=self.var_zMax, width=5, style='opt_valueR.TEntry')
        self.ent_zMin.bind('<Return>', update_zAxis)
        self.ent_zMax.bind('<Return>', update_zAxis)
        
        # Label for the z-axis
        self.var_zlabel = tk.StringVar(value=self.Axis.z_label)
        self.lbl_zlabel = ttk.Label(self.frame, text="Label")
        self.ent_zlabel = ttk.Entry(self.frame, textvariable=self.var_zlabel, width=20, style='opt_valueR.TEntry')
        self.ent_zlabel.bind('<Return>', update_zLabel)
        
        # Choice between linear and log scales for the x-axis 
        self.var_zscale = tk.StringVar(value=self.options_scale[0])
        self.lbl_zscale = ttk.Label(self.frame, text="Scale")
        self.mnu_zscale = ttk.OptionMenu(self.frame, self.var_zscale, self.options_scale[0], *self.options_scale, style='option.TMenubutton')
        self.mnu_zscale["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'])
        if self.Axis.z_scale=="linear": self.var_zscale.set(self.options_scale[0])
        if self.Axis.z_scale=="log":    self.var_zscale.set(self.options_scale[1])
        self.var_zscale.trace('w', update_zscale) # link function to a change of the dropdown options
    
        # Add the labels to the frame 
        self.lbl_zTitle.grid( row=self.row, column=0, columnspan=2, **PAD_TITLE); self.row += 1
        self.lbl_zMin.grid(   row=self.row, column=0, **PAD_LABEL)
        self.ent_zMin.grid(   row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_zMax.grid(   row=self.row, column=0, **PAD_LABEL)
        self.ent_zMax.grid(   row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_zlabel.grid( row=self.row, column=0, **PAD_LABEL)
        self.ent_zlabel.grid( row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_zscale.grid( row=self.row, column=0, **PAD_LABEL)
        self.mnu_zscale.grid( row=self.row, column=1, **PAD_ENTRY); self.row += 1
        self.lbl_space.grid(  row=self.row, column=0, columnspan=2, **PAD_TITLE); self.row += 1
        if True: return

################################################################################
#                                CURVE OPTIONS
################################################################################
       
class tabCurves:
    
    def __init__(self, root, Canvas, tab):
        
        # Get data from the GUI
        self.tab    = tab
        self.root   = root
        self.Canvas = Canvas
        self.Axis   = Canvas.Axis
        
        # Create the frame 
        self.frame = ttk.Frame(tab)
        self.frame.pack(expand=1, fill=tk.BOTH)
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.frame, 0, weight=0) # Appearanc title
        tk.Grid.rowconfigure(   self.frame, 1, weight=0) # Background
        tk.Grid.rowconfigure(   self.frame, 2, weight=0) # Font size
        tk.Grid.rowconfigure(   self.frame, 3, weight=0) # Handle length
        tk.Grid.columnconfigure(self.frame, 0, weight=1, uniform="options") 
        tk.Grid.columnconfigure(self.frame, 0, weight=1, uniform="options") 
        
        # Add elements to the frame
        self.init_appearance()
        return
    
    #-------------------------
    def init_appearance(self):
        
        # Width of the entry widget
        width=5
                
        # Change background color, font size and handle length
        self.lbl_aTitle = ttk.Label(self.frame, text="Appearance", style='prefTitle.TLabel')
        self.var_bg = tk.StringVar(value=self.root.color['bg'])
        self.var_fs = tk.StringVar(value=self.Axis.fontsize)
        self.var_hl = tk.StringVar(value=self.Axis.handlelength)
        self.lbl_bg = ttk.Label(self.frame, text="Background color")
        self.lbl_fs = ttk.Label(self.frame, text="Font size")
        self.lbl_hl = ttk.Label(self.frame, text="Handle length")
        self.ent_bg = ttk.Entry(self.frame, textvariable=self.var_bg, width=width, style='opt_valueR.TEntry')
        self.ent_fs = ttk.Entry(self.frame, textvariable=self.var_fs, width=width, style='opt_valueR.TEntry')
        self.ent_hl = ttk.Entry(self.frame, textvariable=self.var_hl, width=width, style='opt_valueR.TEntry')
    
        # Add the labels to the frame
        self.lbl_aTitle.grid( row=1, column=0, columnspan=2, **PAD_TITLE)
        self.lbl_bg.grid( row=2, column=0, **PAD_LABEL)
        self.ent_bg.grid( row=2, column=1, **PAD_ENTRY)
        self.lbl_fs.grid( row=3, column=0, **PAD_LABEL)
        self.ent_fs.grid( row=3, column=1, **PAD_ENTRY)
        self.lbl_hl.grid( row=4, column=0, **PAD_LABEL)
        self.ent_hl.grid( row=4, column=1, **PAD_ENTRY)
        if True: return
  
################################################################################
#                                LEGEND OPTIONS
################################################################################
    
class tabLegend:
    
    def __init__(self, root, Canvas, tab):
        
        # Get data from the GUI
        self.tab    = tab
        self.root   = root
        self.Canvas = Canvas
        self.Axis   = Canvas.Axis

        # Add some variables
        self.options_loc = ['best', 'upper right', 'upper left', 'lower left', 'lower right', 'right',\
                            'center left', 'center right', 'lower center', 'upper center', 'center'] 
        self.options_bbox = ['Use bbox', 'Disabled bbox']
        
        # Create the frame 
        self.frame = ttk.Frame(tab)
        self.frame.pack(expand=1, fill=tk.BOTH)
        
        # Configure the frame
        tk.Grid.rowconfigure(   self.frame, 0, weight=0) # x title
        tk.Grid.rowconfigure(   self.frame, 1, weight=0) # x min
        tk.Grid.rowconfigure(   self.frame, 2, weight=0) # x max
        tk.Grid.rowconfigure(   self.frame, 3, weight=0) # x label
        tk.Grid.rowconfigure(   self.frame, 4, weight=0) # x scale
        tk.Grid.rowconfigure(   self.frame, 5, weight=0) # y title
        tk.Grid.rowconfigure(   self.frame, 6, weight=0) # y min
        tk.Grid.rowconfigure(   self.frame, 7, weight=0) # y max
        tk.Grid.rowconfigure(   self.frame, 8, weight=0) # y label
        tk.Grid.rowconfigure(   self.frame, 9, weight=0) # y scale 
        tk.Grid.columnconfigure(self.frame, 0, weight=1, uniform="options") 
        tk.Grid.columnconfigure(self.frame, 0, weight=1, uniform="options") 

        # Add elements to the frame
        #self.init_legend() 
    
    #-------------------------   
    def init_legend(self):
        ''' Modify the legend of the Axis. '''

        def modify_legend(**kwargs):
            import matplotlib as mpl
        
            l = self.Axis.ax.legend_
            if l==None: self.Axis.ax.legend()
            l = self.Axis.ax.legend_
        
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
                handler_map = l._custom_handler_map,
            )
        
            if "fontsize" in kwargs and "prop" not in kwargs:
                defaults["prop"].set_size(kwargs["fontsize"])
        
            self.Axis.ax.legend(**dict(list(defaults.items()) + list(kwargs.items())))
            return defaults
    
        def update_location(*args): 
            index = self.options_bbox.index(self.var_bbox.get())
            if index==0: 
                x0 = float(self.var_xbbox.get())
                y0 = float(self.var_ybbox.get())
                modify_legend(loc = self.var_loc.get(), bbox_to_anchor = (x0, y0)) 
            if index==1: 
                modify_legend(loc = self.var_loc.get(), bbox_to_anchor = None)
            self.Canvas.draw_idle()  

        def update_bbox(event):   
            index = self.options_bbox.index(self.var_bbox.get())
            if index==0: 
                x0 = float(self.var_xbbox.get())
                y0 = float(self.var_ybbox.get())
                modify_legend(bbox_to_anchor = (x0, y0)) 
            if index==1: 
                modify_legend(bbox_to_anchor = None)
            self.Canvas.draw_idle()  
            
        # Get the current parameters of the legend
        defaults = modify_legend() 
            
        # Title of the position section
        self.lbl_posTitle = ttk.Label(self.frame, text="Position", style='prefTitle.TLabel') 
        
        # Change the location 
        index = defaults["loc"]
        self.var_loc = tk.StringVar(value=self.options_loc[index])
        self.lbl_loc = ttk.Label(self.frame, text="Location: ")
        self.mnu_loc = ttk.OptionMenu(self.frame, self.var_loc, self.options_loc[index], *self.options_loc, style='option.TMenubutton')
        self.mnu_loc["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg']) 
        self.var_loc.trace('w', update_location) 
        
        # Use bbox  
        self.var_bbox = tk.StringVar(value=self.options_bbox[0])
        self.lbl_bbox = ttk.Label(self.frame, text="Bbox to anchor: ")
        self.mnu_bbox = ttk.OptionMenu(self.frame, self.var_bbox, self.options_bbox[0], *self.options_bbox, style='option.TMenubutton')
        self.mnu_bbox["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg']) 
        self.var_bbox.trace('w', update_location) 
        
        # Change the bbox 
        self.var_xbbox = tk.StringVar(value=defaults["bbox_to_anchor"].x0)
        self.lbl_xbbox = ttk.Label(self.frame, text="   bbox x-value: ")
        self.ent_xbbox = ttk.Entry(self.frame, textvariable=self.var_xbbox, width=20, style='opt_valueR.TEntry')
        self.ent_xbbox.bind('<Return>', update_bbox)
        self.var_ybbox = tk.StringVar(value=defaults["bbox_to_anchor"].y0)
        self.lbl_ybbox = ttk.Label(self.frame, text="   bbox x-value: ")
        self.ent_ybbox = ttk.Entry(self.frame, textvariable=self.var_ybbox, width=20, style='opt_valueR.TEntry')
        self.ent_ybbox.bind('<Return>', update_bbox)
    
        # Add the widgets to the frame
        self.lbl_posTitle.grid( row=0, column=0, columnspan=2, **PAD_TITLE); i=1
        self.lbl_loc.grid(   row=i, column=0, **PAD_LABEL)
        self.mnu_loc.grid(   row=i, column=1, **PAD_ENTRY); i+= 1
        self.lbl_bbox.grid(  row=i, column=0, **PAD_LABEL)
        self.mnu_bbox.grid(  row=i, column=1, **PAD_ENTRY); i+= 1
        self.lbl_xbbox.grid( row=i, column=0, **PAD_LABEL)
        self.ent_xbbox.grid( row=i, column=1, **PAD_ENTRY); i+= 1
        self.lbl_ybbox.grid( row=i, column=0, **PAD_LABEL)
        self.ent_ybbox.grid( row=i, column=1, **PAD_ENTRY); i+= 1
