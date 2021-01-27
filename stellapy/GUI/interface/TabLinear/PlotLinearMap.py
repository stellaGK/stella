
import tkinter as tk
from tkinter import ttk 
from stellapy.GUI.graph_tools import PAD_TITLE2, PAD_LABEL2, PAD_ENTRY2, PAD_ENTRYR
from stellapy.plot.plot_frequencyVsParameter import plot_frequencyVsParameter
from stellapy.plot.plot_frequencyVsParameterVsParameter import plot_frequencyVsParameterVsParameter

#################################################################
#                 CLASS TO PLOT THE LINEAR MAP
#################################################################

class PlotLinearMap: 
    
    def __init__(self, parent, frame):
        
        # Make the parents available
        self.tab = parent
        self.root = parent.root
        self.frame = frame
        
        # Keep track of what is plotted
        self.replot = False
        self.plotted_inputFiles = None
        
        # Attributes for the plot
        self.z_quantity = "gamma"
        self.k_value = "max"
        self.interpolate = False
        self.step  = 4
        self.knob1 = "-"
        self.key1  = "-"
        self.knob2 = "-"
        self.key2  = "-"
        self.colorbar = None
        
        # Make sure there is a linear_map object 
        if 'linear_map' not in self.root.Research.data.keys():
            self.root.Research.data['linear_map'] = {
                'keys' : ['-', '-'],\
                'gamma' : None,\
                'omega' : None,\
                'ky' : None} 
        
        # Option lists
        self.options_par1 = ["rho", "tiprim", "teprim", "fprim", "delta t", "nmu", "nvgrid"]
        self.options_par2 = ["-", "rho", "tiprim", "teprim", "fprim", "delta t", "nmu", "nvgrid"]
        self.options_k    = ["-"]
        
        # Variables to switch the plotting options
        self.var_zQuantity = tk.IntVar(value=1)
        self.var_parameter1 = tk.StringVar(value=self.options_par1[0])
        self.var_parameter2 = tk.StringVar(value=self.options_par2[0])
        self.var_specificKvalue = tk.StringVar(value=self.options_k[0])
        self.var_extraction = tk.IntVar(value=1)
        self.var_interpolate = tk.IntVar(value=0)
        self.var_step = tk.StringVar(value="4")
        
        # Create a subframe for each set of plotting options
        self.subframe_zQuantity = ttk.Frame(frame)
        self.subframe_parameter = ttk.Frame(frame)
        self.subframe_extraction = ttk.Frame(frame)
        self.subframe_options = ttk.Frame(frame)
        
        # Configure the main frame to hold the subframes
        tk.Grid.rowconfigure(   frame, 0, weight=1)      
        tk.Grid.rowconfigure(   frame, 1, weight=0)      
        tk.Grid.rowconfigure(   frame, 2, weight=1)     
        tk.Grid.rowconfigure(   frame, 3, weight=1)     
        tk.Grid.columnconfigure(frame, 0, weight=0) # subframe_yQuantity
        tk.Grid.columnconfigure(frame, 1, weight=0) # subframe_parameter
        tk.Grid.columnconfigure(frame, 2, weight=0) # subframe_extraction
        tk.Grid.columnconfigure(frame, 3, weight=0) # subframe_options

        # Add the subframes to the main frame
        self.subframe_zQuantity.grid( row=1, column=0, padx=25, pady=(0,0), stick='NSEW', rowspan=2)
        self.subframe_parameter.grid( row=1, column=1, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_extraction.grid(row=1, column=2, padx=25, pady=(0,0), stick='NSEW')
        self.subframe_options.grid(   row=1, column=3, padx=25, pady=(0,0), stick='NSEW', rowspan=2)
        
        # Fill the frames with their options
        self.ZQuantity  = self.Frame_ZQuantity( self, self.subframe_zQuantity)
        self.Parameter  = self.Frame_Parameter( self, self.subframe_parameter)
        self.Extraction = self.Frame_Extraction(self, self.subframe_extraction)
        self.Options    = self.Frame_Options(   self, self.subframe_options)
        
        # When the tab is visible, make sure its canvas and frame are loaded
        self.frame.bind("<Visibility>", self.load_figure) 
        return     
    
    #------------------- 
    def load_figure(self, event):
        
        # Attach the correct frame to the GUI
        self.tab.frame_graph3.grid()
        self.tab.frame_graph2.grid_remove()
        
        # Try to automatically detect the correct knob and key
        if self.knob1 == "-" and len(self.root.Research.experiments)!=0:
            if len(self.tab.Simulations.experiment.variedVariables)!=0:
                variable1 = self.tab.Simulations.experiment.variedVariables[0]
                variable2 = self.tab.Simulations.experiment.variedVariables[1]
                self.knob1 = variable1.split(":")[0]
                self.key1  = variable1.split(": ")[-1]
                self.knob2 = variable2.split(":")[0]
                self.key2  = variable2.split(": ")[-1]
                self.var_parameter1.set(self.key1) 
                self.var_parameter2.set(self.key2) 
    
    class Frame_ZQuantity:
        
        def __init__(self, parent, frame):
            
            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
        
            # Configure the frame <self.subframe_yQuantity>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Create the widgets
            self.lbl_zdata    = ttk.Label(frame, text="Z-data", style='prefTitle.TLabel')
            self.rbn_gamma    = ttk.Radiobutton(frame, text='  Growth rate')
            self.rbn_omega    = ttk.Radiobutton(frame, text='  Frequency')
            self.rbn_ky       = ttk.Radiobutton(frame, text='  Wavenumber')
    
            # Add the values and commands to the radiobuttons
            self.rbn_gamma.config(  value=1, variable=self.plot.var_zQuantity, command=self.change_zQuantity)
            self.rbn_omega.config(  value=2, variable=self.plot.var_zQuantity, command=self.change_zQuantity)
            self.rbn_ky.config(     value=3, variable=self.plot.var_zQuantity, command=self.change_zQuantity)
            
            # Add the options to the frame
            self.lbl_zdata.grid(row=0, column=0, **PAD_TITLE2)
            self.rbn_gamma.grid(row=1, column=0, **PAD_LABEL2)
            self.rbn_omega.grid(row=2, column=0, **PAD_LABEL2)
            self.rbn_ky.grid(   row=3, column=0, **PAD_LABEL2)
            return 
        
        #---------------------
        def change_zQuantity(self, *args):
            
            # When changing the plots, reset the axis of the graph
            self.tab.Graph[2].load_defaults()
            self.tab.Graph[2].update_quantitiesAndKeys()
            
            # Extract what needs to be plotted on the y-axis
            z_quantity = self.plot.var_zQuantity.get()
            if z_quantity==1: self.plot.z_quantity = "gamma"; header = "Growthrate versus parameters" 
            if z_quantity==2: self.plot.z_quantity = "omega"; header = "Frequency versus parameters"
            if z_quantity==3: self.plot.z_quantity = "ky";    header = "Wavenumber versus parameters"
            
            # Change the header of the frame and plot the graph
            self.tab.frame_graph3.config(text="   "+header+"  ")
            self.plot.replot = True; self.plot.plot_graph(None)
            return
    
        
    class Frame_Parameter:
        
        def __init__(self, parent, frame):

            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
            self.root = parent.tab.root
        
            # Configure the frame <self.subframe_parameter>
            tk.Grid.columnconfigure(frame, 0, weight=1) 
            
            # Choose the parameter for the x-axis and y-axis
            self.lbl_par  = ttk.Label(frame, text="Parameters", style='prefTitle.TLabel'); width=10
            self.mnu_par1 = ttk.OptionMenu(frame, self.plot.var_parameter1, self.plot.options_par1[0], *self.plot.options_par1, style='option.TMenubutton')
            self.mnu_par1["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
            self.mnu_par1.config(width=width)
            self.mnu_par2 = ttk.OptionMenu(frame, self.plot.var_parameter2, self.plot.options_par2[0], *self.plot.options_par2, style='option.TMenubutton')
            self.mnu_par2["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
            self.mnu_par2.config(width=width)
            def update_parameter1(*args): self.update_parameter(i=1); return
            def update_parameter2(*args): self.update_parameter(i=2); return 
            self.plot.var_parameter1.trace('w', update_parameter1) # link function to a change of the dropdown options
            self.plot.var_parameter2.trace('w', update_parameter2) # link function to a change of the dropdown options
        
            # Add the widgets to the frame
            self.lbl_par.grid( row=0, column=0, **PAD_TITLE2)
            self.mnu_par1.grid(row=1, column=0, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
            self.mnu_par2.grid(row=2, column=0, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
            return

        #-----------------------------
        def update_parameter(self, i):
            
            # Load the defaults when we touch this
            self.tab.Graph[2].load_defaults()
            
            # Get the parameter
            if i==1: var_par = self.plot.var_parameter1.get()
            if i==2: var_par = self.plot.var_parameter2.get()
            
            # Assign the correct parameter
            if var_par=="rho":      knob = "vmec_parameters";        key="rho"
            if var_par=="tprim":    knob = "species_parameters_1";   key="tprim"
            if var_par=="tiprim":   knob = "species_parameters_1";   key="tprim"
            if var_par=="teprim":   knob = "species_parameters_2";   key="tprim"
            if var_par=="fprim":    knob = "species_parameters_1";   key="fprim"
            if var_par=="delta t":  knob = "knobs";                  key="delt"
            if var_par=="delt":     knob = "knobs";                  key="delt"
            if var_par=="nmu":      knob = "vpamu_grids_parameters"; key="rho"
            if var_par=="nvgrid":   knob = "vpamu_grids_parameters"; key="rho"
            if var_par=="-":        knob = "-";                      key="-"
            
            # Save the knobs and keys
            if i==1: self.plot.knob1 = knob;     self.plot.key1 = key
            if i==2: self.plot.knob2 = knob;     self.plot.key2 = key
            return

    class Frame_Extraction:
        
        def __init__(self, parent, frame):

            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
            self.root = parent.tab.root
            
            # Configure the frame <self.subframe_extraction>
            tk.Grid.columnconfigure(frame, 0, weight=0) 
            tk.Grid.columnconfigure(frame, 1, weight=0) 
            tk.Grid.columnconfigure(frame, 2, weight=1) 
            
            # Choose which growth rate we extract: the maximum one or at a specific wavenumber
            self.lbl_extraction = ttk.Label(frame, text="Extraction", style='prefTitle.TLabel')
            self.rbn_gammamax = ttk.Radiobutton(frame, text='  Plot the most unstable mode')
            self.rbn_gammak   = ttk.Radiobutton(frame, text='  Plot the mode at ky =')
      
            # Add the commands to the radio buttons
            self.rbn_gammamax.config( value=1, variable=self.plot.var_extraction, command=self.change_extraction)
            self.rbn_gammak.config(   value=2, variable=self.plot.var_extraction, command=self.change_extraction)
      
            # Drop down menu to choose the k-value
            self.mnu_kvalue = ttk.OptionMenu(frame, self.plot.var_specificKvalue, self.plot.options_k[0], *self.plot.options_k, style='option.TMenubutton'); width=4
            self.mnu_kvalue["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
            self.mnu_kvalue.config(width=width)
            self.plot.var_specificKvalue.trace('w', self.update_selectedMode) # link function to a change of the dropdown options
        
            # Add the widgets to the frame
            self.lbl_extraction.grid(row=0,  column=0, **PAD_TITLE2, columnspan=3)
            self.rbn_gammamax.grid(  row=1,  column=0, **PAD_LABEL2, columnspan=3)
            self.rbn_gammak.grid(    row=2,  column=0, **PAD_LABEL2)
            self.mnu_kvalue.grid(    row=2,  column=1, **PAD_LABEL2) 
            return 
        
        #-----------------------------------------
        def update_selectedMode(self, *args):
            self.plot.k_value = float(self.plot.var_specificKvalue.get())
            self.plot.replot=True; self.plot.plot_graph(None) 
                
        #-----------------------------------------
        def change_extraction(self):
            if self.plot.var_extraction.get()==1: self.plot.k_value = "max"
            if self.plot.var_extraction.get()==2: self.plot.k_value = float(self.plot.var_specificKvalue.get())
            self.plot.replot=True; self.plot.plot_graph(None) 
            return
        
    class Frame_Options:
        
        def __init__(self, parent, frame):
        
            # Make the parents available
            self.plot = parent 
            self.tab = parent.tab
            self.root = parent.tab.root
              
            # Create the widgets
            self.lbl_options = ttk.Label(frame, text="Options", style='prefTitle.TLabel')
            self.chk_interp = ttk.Checkbutton(frame, text=" Interpolate z-data with step = ")
            self.btn_rewrite = ttk.Button(frame, text="Reread linear map", width=20)
    
            # Add the commands
            self.chk_interp.config(variable=self.plot.var_interpolate, command=self.update_interp)
            self.btn_rewrite.config(command=lambda: self.rewrite_linearMap())
            
            # Drop down menu to choose the interpolation step 
            self.ent_interp = ttk.Entry(frame, textvariable=self.plot.var_step, font=("Courier New", 11), style='opt_valueCBold.TEntry', width=5)
            self.ent_interp.bind('<Return>', self.update_step) 
            
            # Configure the frame
            tk.Grid.rowconfigure(frame, 0, weight=0) 
            tk.Grid.rowconfigure(frame, 1, weight=0) 
            tk.Grid.columnconfigure(frame, 0, weight=0) 
            tk.Grid.columnconfigure(frame, 1, weight=0) 
            tk.Grid.columnconfigure(frame, 2, weight=1) 
            
            # Add the options to the frame
            self.lbl_options.grid(row=0, column=0, **PAD_TITLE2, columnspan=2)
            self.chk_interp.grid( row=1, column=0, **PAD_LABEL2)
            self.ent_interp.grid( row=1, column=1)
            self.btn_rewrite.grid(row=2, column=0, **PAD_LABEL2, columnspan=2) 
            return 
                           
        #-------------------
        def update_interp(self, *args): 
            if self.plot.var_interpolate.get()==0: self.plot.interpolate = False
            if self.plot.var_interpolate.get()==1: self.plot.interpolate = True
            self.plot.replot=True; self.plot.plot_graph(None) 
            return 
                           
        #-------------------  
        def update_step(self, *args):
            self.plot.step = int(self.plot.var_step.get()) 
            self.plot.replot=True; self.plot.plot_graph(None) 
            return 
        
        #----------------------
        def rewrite_linearMap(self):
            self.root.Research.data['linear_map'] = {
                'parameters1' : None,\
                'parameters2' : None,\
                'gamma' : None,\
                'omega' : None,\
                'ky' : None} 
            self.plot.replot=True; self.plot.plot_graph(None) 
            


    def plot_graph(self, poppedout_id=None):
        
        # Perhaps the button was clicked on the popped out window
        if poppedout_id == None: Graph = self.tab.Graph[2]
        if poppedout_id != None: Graph = self.root.graph_poppedOut[poppedout_id]
        
        # Only proceed if there are input_files
        if self.root.input_files != []:
            
            # Get the input_files that need to be plotted
            input_files = self.root.Research.input_files
            
            # If the selection of input_files changed, then reset the axis
            if (self.plotted_inputFiles != input_files):
                Graph.load_defaults()
                
            # If there are files plot them, if they changed plot them, if its a popout plot them
            if (input_files != [] and self.plotted_inputFiles != input_files) or poppedout_id!=None or self.replot==True:
                
                # Clear the axis and remember which input_files are plotted
                Graph.ax.clear()
                self.plotted_inputFiles = input_files.copy() 
                self.replot = False
                
                # While plotting show a progress bar turn the cursor into a waiting cursor
                self.tab.show_progress("start", poppedout_id) 
                
                # Remove the colorbar
                try: self.colorbar.remove()
                except: pass
                
                # Make sure the saved linear map has the correct keys
                keys = self.root.Research.data['linear_map'].keys()
                for key in ['key1', 'key2', 'gamma', 'omega', 'ky', self.key1, self.key2]: 
                    if key not in keys:
                        self.root.Research.data['linear_map'][key] = None
                
                # Be able to flip the parameter axis if the data is saved
                if  self.key1==self.root.Research.data['linear_map']['key2']\
                and self.key2==self.root.Research.data['linear_map']['key1']:
                    self.root.Research.data['linear_map']['key1']  = self.key1
                    self.root.Research.data['linear_map']['key2']  = self.key2
                    self.root.Research.data['linear_map']['gamma'] = self.root.Research.data['linear_map']['gamma'].transpose()
                    self.root.Research.data['linear_map']['omega'] = self.root.Research.data['linear_map']['omega'].transpose()
                    self.root.Research.data['linear_map']['ky']    = self.root.Research.data['linear_map']['ky'].transpose()
                
                # Plot the graph
                self.colorbar, data = plot_frequencyVsParameterVsParameter(\
                    # Specify which simulations to plot
                        research=self.root.Research,\
                        experiment_id="All experiments",#self.tab.Simulations.experiment_id,\
                        simulation_id="All simulations",#self.tab.Simulations.simulation_id,\
                        parameter_knob1=self.knob1,\
                        parameter_key1=self.key1,\
                        parameter_knob2=self.knob2,\
                        parameter_key2=self.key2,\
                    # Specify data range
                        z_quantity=self.z_quantity,\
                        x_range=Graph.range["x"],\
                        y_range=Graph.range["y"],\
                    # Details of the modes
                        kx_range=Graph.kx, \
                        ky_range=Graph.ky, \
                        k_value=self.k_value, \
                        lineardata="average",\
                    # Save time by reading the data elsewhere
                        parameters1=self.root.Research.data['linear_map'][self.key1],\
                        parameters2=self.root.Research.data['linear_map'][self.key2],\
                        gamma_mostUnstableMode=self.root.Research.data['linear_map']['gamma'],\
                        omega_mostUnstableMode=self.root.Research.data['linear_map']['omega'],\
                        ky_mostUnstableMode=self.root.Research.data['linear_map']['ky'],\
                    # Labels
                        x_label=Graph.label["x"],\
                        y_label=Graph.label["y"],\
                        title=None,\
                    # For the GUI the figure object already exists  
                        show_figure = False,\
                        ax=Graph.ax, \
                        Progress=self.tab.Progress,\
                        root=self.root,\
                    # Extra options
                        interpolate=self.interpolate,\
                        step=self.step)
                    
                # Save the t_range to the research object
                if data is not None: self.root.Research.data['linear_map'] = data
                
                # Update the <graph> class so the option window can use them
                Graph.update_rangesAndLabels()
            
                # When finished show the canvas and change the cursor back.
                self.tab.show_progress("finished", poppedout_id) 
    
        # If no simulations have been selected, clear the current figure
        if self.root.input_files == [] or self.root.Research.input_files == []:
            Graph.ax.clear()
            Graph.load_defaults()
            self.tab.show_progress("nothing", poppedout_id) 
        
        # Update screen
        if poppedout_id==None: self.tab.Canvas[2].draw_idle()
        if poppedout_id!=None: self.root.canvasPoppedOut[poppedout_id].draw_idle()
        if True: return