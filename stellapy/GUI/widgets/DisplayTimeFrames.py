
import numpy as np
import tkinter as tk
from tkinter import ttk  
from bisect import bisect
from stellapy.GUI.utils.get_centerCoordinatesOfMonitor import get_centerCoordinatesOfMonitor
from stellapy.data.time.get_timeFrame import get_timeFrame

################################################################################
#   CREATE A TREEVIEW TO DISPLAY AND EDIT THE TIME FRAMES FOR EACH SIMULATION
################################################################################

class DisplayTimeFrames:

    def __init__(self, parent):

        # Make the parents available
        self.tab = parent
        self.root = parent.root
        
        # Attributes for the plot 
        self.tree_experiments = {}
        self.frame_species = {} 
        self.tree = {}  
        
    #--------------------------
    def select_timeFrames(self):
        
        # Create the preferences window and keep it on top
        self.window_timeFrames = tk.Toplevel(self.root, bg=self.root.color['bg'])
        self.window_timeFrames.title("Select time frames") 
        self.window_timeFrames.attributes('-topmost', 'true')
        self.rescale_window()
        
        # Add the tabbed view to the popout window   
        self.initiate_frameOptions_timeFrame()   
        self.tab_species.grid(row=0, column=0, sticky="NSEW")
        self.update_speciesFrame()
        self.update_timeFramesAndSaturatedFluxes()
        
        def on_closing():
            self.tree_experiments = {}
            self.frame_species = {} 
            self.tree = {} 
            self.window_timeFrames.destroy()
        self.window_timeFrames.protocol("WM_DELETE_WINDOW", on_closing)
        
    #-------------------------
    def rescale_window (self):    
        ''' Make sure the size of the options window is decent. '''
    
        # Make sure the new window is centered
        self.window_timeFrames.withdraw()
        self.window_timeFrames.update_idletasks()   
        x, y = get_centerCoordinatesOfMonitor(self.root, self.window_timeFrames)  
        
        # Standard size of the screen
        self.winx = 1000
        self.winy = 600
                       
        # Center the new window in the screen
        x = x  - self.winx/2
        y = y - self.winy/2
        self.window_timeFrames.geometry("+%d+%d" % (x, y))  
        self.window_timeFrames.deiconify() 
        
    #-------------------------------------------       
    def initiate_frameOptions_timeFrame(self):
        ''' Choose the time frame to measure the saturated fluxes and
        display some information about the saturation of the fluxes. '''
         
        # Tabbed display for species 1, 2, 3, ....
        self.tab_species = ttk.Notebook(self.window_timeFrames, style='species.TNotebook', width=self.winx, height=self.winy)
        self.initiate_speciesTab(0,"  Ions  ")
        tk.Grid.rowconfigure(self.window_timeFrames, 0, weight=1) 
        tk.Grid.columnconfigure(self.window_timeFrames, 0, weight=1)  
        tk.Grid.rowconfigure(self.tab_species, 0, weight=1) 
        tk.Grid.columnconfigure(self.tab_species, 0, weight=1)  
        return 

    #-------------------------------------------    
    def initiate_speciesTab(self, specie_id, title):
         
        # For each specie we make a tabbed view
        self.frame_species[specie_id] = ttk.Frame(self.tab_species, padding=(10,10,10,10))
        self.tab_species.add(self.frame_species[specie_id], text=title) 
        tk.Grid.rowconfigure(self.frame_species[specie_id], 0, weight=1)
        tk.Grid.columnconfigure(self.frame_species[specie_id], 0, weight=1)
        self.tree_experiments[specie_id] = []
         
        # The list of time frames will be displayed inside a treeview widget
        self.tree[specie_id] = ttk.Treeview(self.frame_species[specie_id], show="tree",columns=("#0","#1","#2","#3"))
        self.tree[specie_id].grid(row=0, column=0, padx=0, pady=0, sticky='NESW')
        tk.Grid.rowconfigure(self.tree[specie_id], 0, weight=1)
        tk.Grid.columnconfigure(self.tree[specie_id], 0, weight=1)

        # Create a treeview widget and bind keypresses to the widget
        self.tree[specie_id].bind("<Double-1>", self.onDoubleClick)
         
        # Add columns to the tree view 
        self.tree[specie_id].heading("#0", text="Simulation",     anchor=tk.W) 
        self.tree[specie_id].heading("#1", text="Available time", anchor=tk.W) 
        self.tree[specie_id].heading("#2", text="Time frame",     anchor=tk.W) 
        self.tree[specie_id].heading("#3", text="Saturated flux", anchor=tk.W) 
        self.tree[specie_id]["displaycolumns"] = ("#0", "#1", "#2", "#3")
 
        # When the tab is visible, make sure the columns have a good size
        def resize_columns(*_): 
            self.tree[specie_id].column("#0", width=int(self.tree[specie_id].winfo_width()*4/10), stretch=tk.YES)
            self.tree[specie_id].column("#1", width=int(self.tree[specie_id].winfo_width()*2/10), stretch=tk.YES) 
            self.tree[specie_id].column("#2", width=int(self.tree[specie_id].winfo_width()*2/10), stretch=tk.YES)  
            self.tree[specie_id].column("#3", width=int(self.tree[specie_id].winfo_width()*2/10), stretch=tk.YES)   
        self.frame_species[specie_id].bind("<Visibility>", resize_columns)
        self.tab_species.bind("<Visibility>", resize_columns)
        if True: return

    def update_speciesFrame(self):

        # Find the maximum number of species
        dim_species = 0
        for experiment in self.root.Research.experiments:
            dim_species = max(dim_species, experiment.simulations[0].dim.species)
            
        # If there are more species, add more tabs 
        frame_ids = list(self.frame_species.keys()) 
        
        # When all simulations are removed, remove the additional species tabs
        try: unique_folders = self.root.Research.unique_folders
        except: unique_folders = []
        if len(unique_folders)==0: dim_species=1      
        
        # Add more tabs if there aren't enough
        if len(frame_ids) < dim_species: 
            for i in range(len(frame_ids), dim_species):
                self.initiate_speciesTab(i, "Species "+str(i+1))
                
        # Remove tabs if there are too many
        if dim_species < len(frame_ids): 
            for i in range(dim_species, len(frame_ids)):
                self.tab_species.forget(self.frame_species["tab "+str(i+1)])
                self.frame_species["tab "+str(i+1)].destroy()                                       
        
        # Change the name if the tabs are named wrongly
        if not dim_species==1:
            if dim_species == 2: # Assume the second species is always kinetic electrons 
                try: self.tab_species.tab(self.frame_species[1], text = 'Kinetic electrons')
                except: pass
            elif dim_species== 3: # Assume second species is always kinetic electrons and third impurities
                self.tab_species.tab(self.frame_species[1], text = 'Kinetic electrons')
                self.tab_species.tab(self.frame_species[2], text = 'Impurities')
        return 
    
    #---------------------------------------              
    def update_timeFramesAndSaturatedFluxes(self):
        ''' Update the time frames to calculate the saturated flux. '''

        # Only update when there are input files 
        if len(self.root.TabSelectedFiles.input_files) != 0:
            
            # Remove all the items from the tree view
            for specie_id in range(len(self.tree_experiments)):
                self.tree[specie_id].delete(*self.tree[specie_id].get_children())
                self.tree_experiments[specie_id] = []

            # Iterate over the experiments
            for experiment in self.root.Research.experiments:
                
                # Iterate over the species to fill each tab
                for specie_id in range(experiment.simulations[0].dim.species):

                    # Check whether the experiment is already in the treeview
                    if not self.tree[specie_id].exists(experiment.id):
                        contents = [self.tree[specie_id].item(child)["text"] for child in self.tree[specie_id].get_children("")]
                        shown_numbers = [ int("0"+"".join([s for s in name if s.isdigit()])) for name in contents ]
                        new_number = int("0"+"".join([s for s in experiment.id if s.isdigit()]))
                        index = bisect(shown_numbers, new_number)
                        text = experiment.line_label
                        self.tree_experiments[specie_id].append(self.tree[specie_id].insert(parent="", index=index, iid=experiment.id, text=text, values=("", "")))
                        self.tree[specie_id].item(experiment.id, open=True)
            
                    # Add the values to the experiment: marker label; available time; time frame; Q_sat
                    for tree_experiment in self.tree_experiments[specie_id]:
                        self.tree[specie_id].selection_set(tree_experiment)  
                        experiment_iid = self.tree[specie_id].selection()[0]
                        if experiment_iid==experiment.id:
                            for simulation in experiment.simulations:
                                vec_time = simulation.fluxes.qflux_vs_ts.t
                                m = simulation.marker_label.replace("$", "").replace("\,", " ")
                                q = "Q = " + "{:.2e}".format(simulation.time.saturatedFluxes["qflux"][specie_id]) 
                                t = int(vec_time[len(vec_time[~np.isnan(vec_time)]) - 1])
                                t = "    t_end = " + str(t)
                                try: 
                                    dt = get_timeFrame(simulation, normalize=False, y_quantity=None)
                                    if self.tab.PlotTimeEvolution.normalizeBySaturation:
                                        time_peak = simulation.get_peakTime()
                                        dt = [round(dt[0]/time_peak,2), round(dt[1]/time_peak,2)]   
                                    if not self.tab.PlotTimeEvolution.normalizeBySaturation:
                                        dt = [int(dt[0]), int(dt[1])]
                                    dt = "dt = [" + str(dt[0]) + ", " + str(dt[1]) + "]"
                                except: dt = "[ - ]"
                                contents = [self.tree[specie_id].item(child)["text"] for child in self.tree[specie_id].get_children(experiment.id)]
                                shown_numbers = [ int("0"+"".join([s for s in name if s.isdigit()])) for name in contents ]
                                new_number = int("0"+"".join([s for s in m if s.isdigit()]))
                                if self.tree[specie_id].exists(simulation.id):
                                    self.tree[specie_id].item(simulation.id, text=m, values=(t, dt, q)) 
                                if not self.tree[specie_id].exists(simulation.id):
                                    self.tree[specie_id].insert(tree_experiment, bisect(shown_numbers, new_number), iid=simulation.id, text=m, values=(t, dt, q))           
                                
                                
        # Remove focus from items
        for specie_id in range(len(self.tree_experiments)):
            for item in self.tree[specie_id].selection():
                self.tree[specie_id].selection_remove(item)
        return 
        
    #----------------------------------
    def update_treeViewAfterChange(self, values, simulation_iid, experiment_iid, specie_id):
        ''' After manually changing the time frames, update the figure. '''
        
        # Identify the simulation
        simulation = None
        for experiment_temp in self.root.Research.experiments:
            for simulation_temp in experiment_temp.simulations: 
                if simulation_iid == simulation_temp.id:
                    simulation = simulation_temp
                    
        # Add the normalization factor
        norm = simulation.get_peakTime() if self.tab.PlotTimeEvolution.normalizeBySaturation else 1
                        
        # Update the t_range object  
        tstart = values[1].split("[")[-1].split(",")[0]
        tend = values[1].split(",")[-1].split("]")[0]
        simulation.time.update_timeFrame([tstart, tend]*norm)  
        
        # Replot the figures of the tab 
        self.tab.plot()
        return 
    
    #----------------------------------
    def onDoubleClick(self, event):
        ''' Executed, when a row is double-clicked. Opens EntryPopup above the item's column, 
        so it is possible to change the text. '''
    
        # Close previous popup if there is one
        try: self.entryPopup.on_return()
        except: pass
    
        # Get the treeview that was selected
        specie_id = self.tab_species.index(self.tab_species.select())
    
        # Select row and column that was clicked on
        rowid = self.tree[specie_id].identify_row(event.y)
        columnid = self.tree[specie_id].identify_column(event.x)
    
        # Get column position info
        x,y,_,height = self.tree[specie_id].bbox(rowid, columnid)
    
        # Place Entry popup properly 
        if columnid=="#2":
            text = self.tree[specie_id].item(rowid, 'text')
            values = self.tree[specie_id].item(rowid, 'values')
            self.entryPopup = EntryPopup(self.root, self.tree[specie_id], specie_id, rowid, columnid, text, values, self.update_treeViewAfterChange)
            self.entryPopup.place(x=x, y=y+height//2, anchor=tk.W)
        if True: return 
     
     
########################################################
#                   ENTRY POPUP
########################################################

class EntryPopup(ttk.Entry):

    def __init__(self, root, tree, specie_id, iid, columnid, text, values, update_treeViewAfterChange, **kw):
        super().__init__(tree, style='tree.TEntry', **kw)
        
        # Save the information about the treeview
        self.specie_id = specie_id
        self.tree = tree 
        self.update_tree = update_treeViewAfterChange
        self.iid = iid
        self.columnid = columnid
        self.text = text
        self.values = list(values)
        self['exportselection'] = False
        
        # Insert the text of the cell in the entry widget: either its text or one of the values
        if "0" in self.columnid: self.insert(0, text) 
        for i in range(1,len(values)+1):
            if str(i) in self.columnid: self.insert(0, values[i-1]) 

        # Bind the key presses to the entry widget
        self.focus_force()
        self.bind("<Return>", self.on_return)
        self.bind("<Control-a>", self.select_all)
        self.bind("<Escape>", lambda *ignore: self.destroy())
        return 
    
    #-------------------------
    def on_return(self, event=None):
        
        # Get the text/value that was changed  
        if "0" in self.columnid: self.text = self.get()
        for i in range(1,len(self.values)+1):
            if str(i) in self.columnid: 
                self.values[i-1] = self.get()
                
        # Replace the text/value in the treeview with this new text/value
        self.tree.item(self.iid, text=self.text, values = self.values) 

        # Also update the t_range
        simulation_iid = self.tree.selection()[0]
        experiment_iid = self.tree.parent(simulation_iid)
        self.update_tree(self.values, simulation_iid, experiment_iid, self.specie_id)
        self.destroy()
        return 
    
    #-------------------------
    def select_all(self, *ignore):
        ''' Set selection on the whole text '''
        self.selection_range(0, 'end')

        # returns 'break' to interrupt default key-bindings
        return 'break'


