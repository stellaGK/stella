####################################################################
#              CLASS FOR THE FRAME "RESEARCH"
####################################################################
''' SUBFRAME RESEARCH ON TAB 1 '''

import pickle
import pathlib
import tkinter as tk
import matplotlib as mpl 
from tkinter import ttk
from bisect import bisect 
from stellapy.utils.config import CONFIG 
from stellapy.GUI.interface import update_GUI 
from stellapy.simulations.Research import create_research
from stellapy.GUI.widgets import PAD_TITLE2, PAD_LABEL2, PAD_ENTRY2  #@UnusedImport 
from stellapy.data.stella.load_stellaKnobsAndKeys import load_stellaKnobsAndKeys
stellaKnobsAndKeys = load_stellaKnobsAndKeys() 

####################################################################
# INITIALIZE THE FRAME WITH INFORMATION ON THE SELECTED SIMULATIONS
####################################################################

class Research():
    '''
    Calling this class initiates the "Research" LabelFrame on the middle left side 
    of the tab "Simulations" attached to the root.
        

    Parent widgets
    --------------
    root:               Tk()
    tabheader:          ttk.Notebook(root)
    tab:               ttk.Frame(tabheader)          
    frame_left:         ttk.Frame(tab)
    frame_research:     ttk.LabelFrame(frame_left)    --> frame

    '''

#################################################################
#                          WIDGETS
#################################################################

    def __init__(self,tab):
        '''
        Initialize the widgets in the "Research" LabelFrame of the tab "Simulations".
            - A scrollable frame to show the experiments and simulations
            - A checkbox to ignore the resolution
            - A entry box to choose the name of the research (for reloading)
            - A button to load previous researches
            - A button to save the current research (will be quite automatic tho)

        '''

        #====================================================================
        # Safe info from the root so that it can be passed on or used locally
        #====================================================================
        
        self.tab = tab            
        self.root = tab.root 
        self.frame = tab.frame_research
                    
        #===========
        # VARIABLES
        #===========
        
        # Basic variables
        self.research_name = "DEFAULT"
        self.ignore_resolution = False
        self.number_variedVariables = 1
        self.folderIsExperiment = False
        self.resolutionScan = False
        self.experiment_knob = "vmec_parameters"
        self.experiment_key = "vmec_filename"
        self.experiment_knob2 = "-"
        self.experiment_key2 = "-"
        self.tree_experiments = []
        self.initialdir = CONFIG['PATHS']['GUI_Pickles']

        # Dropdown menus
        self._options_resolution = ["   Ignore resolution", "   Include resolution"]
        self._options_variedParam = ["  0 varied parameters",\
                                     "  1 varied parameter", "  2 varied parameters", "  3 varied parameters",
                                     "  4 varied parameters", "  5 varied parameters"]
        
        # Knob and key options 
        self.options_knobs = list(stellaKnobsAndKeys.keys()) 
        self.options_knobs.sort()
        
        #============================
        # WIDGETS CREATION FOR FRAME
        #============================
        
        # The list of researches will be displayed inside a treeview widget
        self.tree = ttk.Treeview(self.frame, show="tree", columns=("#0"))
        style = ttk.Style()
        style.configure("Treeview.Heading", font=(None, 100))

        # "Save research" button which allows to set the name of the research.
        self.button_saveResearch = ttk.Button(self.frame, text="Save research", width=14)
        self.button_saveResearch.config(command=lambda: self.save_research())

        # "Load research" button to choose one of the previous researches
        self.button_loadResearch = ttk.Button(self.frame, text="Load research", width=14)
        self.button_loadResearch.config(command=lambda: self.load_research())

        # "Ignore resolution" and "Include resolution" option 
        self.button_knobkey = ttk.Button(self.frame, text="Sort by knob/key", width=16) 
        self.button_knobkey.config(command=lambda: self.openMenuToChooseKnobAndKey())    
        
        # "X varied parameters" option        
        self.popup_variedParam = tk.Menu(self.frame, tearoff=0)
        self.popup_variedParam.add_command(label="0 varied parameters", command= lambda: update_btn2("0 varied parameters")) 
        self.popup_variedParam.add_command(label="1 varied parameter",  command= lambda: update_btn2("1 varied parameter")) 
        self.popup_variedParam.add_command(label="2 varied parameters", command= lambda: update_btn2("2 varied parameters")) 
        self.popup_variedParam.add_command(label="3 varied parameters", command= lambda: update_btn2("3 varied parameters")) 
        self.popup_variedParam.add_command(label="4 varied parameters", command= lambda: update_btn2("4 varied parameters")) 
        self.popup_variedParam.add_command(label="5 varied parameters", command= lambda: update_btn2("5 varied parameters"))
        self.popup_variedParam.add_command(label="10 varied parameters", command= lambda: update_btn2("10 varied parameters"))
        self.button_variedParam = ttk.Button(self.frame, text="1 varied parameter", width=16)  
        def update_btn2(x):
            # When an option is chosen from the popup menu, change the text on the button
            self.button_variedParam.config(text=x)
            # Get the setting for ignoring the resolution
            if x=="0 varied parameters":        self.number_variedVariables = 0
            if x=="1 varied parameter":         self.number_variedVariables = 1
            if x=="2 varied parameters":        self.number_variedVariables = 2
            if x=="3 varied parameters":        self.number_variedVariables = 3
            if x=="4 varied parameters":        self.number_variedVariables = 4
            if x=="5 varied parameters":        self.number_variedVariables = 5
            if x=="10 varied parameters":       self.number_variedVariables = 10
            # Remake the simulation object with the new setting
            self.create_researchObject()
            # Update the GUI
            update_GUI(self.root)
        def btn_popup2(event):
            try: self.popup_variedParam.tk_popup(event.x_root, event.y_root)
            finally: self.popup_variedParam.grab_release()
        self.button_variedParam.bind("<Button-1>", btn_popup2)
        
        # A dot dot dot button with extra options
        self.popup_dotdotdot = tk.Menu(self.frame, tearoff=0)
        self.popup_dotdotdot.add_command(label='Add "1 folder=1 experiment" option', command = lambda: self.update_menu(1))
        self.popup_dotdotdot.add_command(label='Remove "1 folder=1 experiment" option', command = lambda: self.update_menu(2))
        self.popup_dotdotdot.add_command(label='Include resolution', command = lambda: self.update_menu(3))
        self.popup_dotdotdot.add_command(label='Ignore resolution', command = lambda: self.update_menu(4)) 
        self.popup_dotdotdot.add_command(label='Resolution scan (yes)', command = lambda: self.update_menu(5)) 
        self.popup_dotdotdot.add_command(label='Resolution scan (no)', command = lambda: self.update_menu(6)) 
        self.button_dotdotdot = ttk.Button(self.frame, text="...", width=4) 
        def popup_after_leftclick(event):
            try: self.popup_dotdotdot.tk_popup(event.x_root, event.y_root)
            finally: self.popup_dotdotdot.grab_release()
        self.button_dotdotdot.bind("<Button-1>", popup_after_leftclick)  
        
        
        #=================
        # CONFIGURE FRAME
        #=================
        tk.Grid.rowconfigure(   self.frame, 0, weight=1) # Scrollable canvas
        tk.Grid.rowconfigure(   self.frame, 1, weight=0) # 3 buttons to edit the simulation selection
 
        tk.Grid.columnconfigure(self.frame, 0, weight=1) # Column for button_saveResearch
        tk.Grid.columnconfigure(self.frame, 1, weight=1) # Column for button_loadResearch
        tk.Grid.columnconfigure(self.frame, 2, weight=1) # Column for button_resolution
        tk.Grid.columnconfigure(self.frame, 3, weight=1) # Column for button_variedParam
        tk.Grid.columnconfigure(self.frame, 4, weight=0) # Column for button_dotdotdot

        #======================
        # WIDGETS ARRANGEMENT
        #=====================
        self.tree.grid(                in_=self.frame, row=0, column=0, padx=8, pady=4, sticky='nesw', columnspan=4)
        self.button_saveResearch.grid( in_=self.frame, row=1, column=0, padx=8, pady=4, sticky='nw', ipady=7)
        self.button_loadResearch.grid( in_=self.frame, row=1, column=1, padx=8, pady=4, sticky='nw', ipady=7) 
        self.button_knobkey.grid(      in_=self.frame, row=1, column=2, padx=8, pady=4, sticky='nw', ipady=7) 
        self.button_variedParam.grid(  in_=self.frame, row=1, column=3, padx=8, pady=4, sticky='nw', ipady=7) 
        self.button_dotdotdot.grid(    in_=self.frame, row=1, column=4, padx=8, pady=4, sticky='nw', ipady=7) 

        #============
        # TREE VIEW
        #============
        
        # Add columns  
        self.tree.heading("#0", text="Simulation",anchor=tk.W) 
        self.tree["displaycolumns"] = ("#0")

        # When the tab is visible, make sure the columns have a good size
        def resize_columns(*args): 
            self.tree.column("#0", width=int(self.tree.winfo_width()), stretch=tk.YES)
        self.frame.bind("<Visibility>", resize_columns)
    
        # Prevent folding of header comment on next line
        if True: return

#==================================================
# For the popup menu on x varied parameters button
#===================================================

    def popup_menu(self, event):
        ''' Rightclick event linked to the scrollable canvas. '''
        try:        self.rightclick_menu.tk_popup(event.x_root, event.y_root)
        finally:    self.rightclick_menu.grab_release()
        return
        
    #------------------------
    def update_menu(self, i):
        if i==1: self.folderIsExperiment = True
        if i==2: self.folderIsExperiment = False
        if i==3: self.ignore_resolution = False 
        if i==4: self.ignore_resolution = True 
        if i==5: self.resolutionScan = True 
        if i==6: self.resolutionScan = False  
        if len(self.tab.input_files) != 0 and i in [1,2,3,4,5,6]: 
            self.create_researchObject()
            update_GUI(self.root)
        return 
    
    #---------------------
    def openMenuToChooseKnobAndKey(self):
        
        # Create a top window to ask for the parameter knob and key
        research_window = tk.Toplevel(bg=self.root.color['bg'])
        research_window.title("Sort experiments by stella knob and key")
        research_window.withdraw()
        self.root.eval(f'tk::PlaceWindow {str(research_window)} center')
        
        # When updating the stella knob, set the options list for the stella key
        def update_key(*_):
            knob = self.var_knob.get()
            key_list = stellaKnobsAndKeys[knob]; key_list.sort()
            self.mnu_key['menu'].delete(0, 'end')
            for knob in key_list:
                self.mnu_key['menu'].add_command(label=knob, command=tk._setit(self.var_key, knob))
            knob2 = self.var_knob2.get()
            key_list = stellaKnobsAndKeys[knob2]; key_list.sort()
            self.mnu_key2['menu'].delete(0, 'end')
            for knob2 in key_list:
                self.mnu_key2['menu'].add_command(label=knob2, command=tk._setit(self.var_key2, knob2))
    
        # Choose the stella knob
        self.lbl_title  = ttk.Label(research_window, text="Choose the stella knob and key:") 
        self.var_knob = tk.StringVar(value=self.options_knobs[22]); width=30; 
        self.lbl_knob = ttk.Label(research_window, text="    Knob: ", style='prefTitle.TLabel')
        self.mnu_knob = ttk.OptionMenu(research_window, self.var_knob, self.options_knobs[22], *self.options_knobs, style='option.TMenubutton')
        self.mnu_knob["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_knob.config(width=width)
        self.var_knob.trace('w', update_key) # link function to a change of the dropdown options
  
        # Choose the second stella knob
        self.lbl_title2  = ttk.Label(research_window, text="Choose a second stella knob and key:")
        self.var_knob2 = tk.StringVar(value=self.options_knobs[0]); width=30; 
        self.lbl_knob2 = ttk.Label(research_window, text="    Knob: ", style='prefTitle.TLabel')
        self.mnu_knob2 = ttk.OptionMenu(research_window, self.var_knob2, self.options_knobs[0], *self.options_knobs, style='option.TMenubutton')
        self.mnu_knob2["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_knob2.config(width=width)
        self.var_knob2.trace('w', update_key) # link function to a change of the dropdown options
        
        # Choose the stella key
        knob = self.var_knob.get()
        key_list = stellaKnobsAndKeys[knob]; key_list.sort()
        self.var_key = tk.StringVar(value=key_list[5])
        self.lbl_key = ttk.Label(research_window, text="    Key: ", style='prefTitle.TLabel') 
        self.mnu_key = ttk.OptionMenu(research_window, self.var_key, key_list[5], *key_list, style='option.TMenubutton')
        self.mnu_key["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_key.config(width=width)
        
        # Choose the second stella key
        knob2 = self.var_knob2.get()
        key_list = stellaKnobsAndKeys[knob2]; key_list.sort()
        self.var_key2 = tk.StringVar(value=key_list[0])
        self.lbl_key2 = ttk.Label(research_window, text="    Key: ", style='prefTitle.TLabel') 
        self.mnu_key2 = ttk.OptionMenu(research_window, self.var_key2, key_list[0], *key_list, style='option.TMenubutton')
        self.mnu_key2["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_key2.config(width=width)

        # Configure the frame
        tk.Grid.rowconfigure(   research_window, 0, weight=0) 
        tk.Grid.rowconfigure(   research_window, 1, weight=0) 
        tk.Grid.rowconfigure(   research_window, 2, weight=0) 
        tk.Grid.rowconfigure(   research_window, 3, weight=0) 
        tk.Grid.rowconfigure(   research_window, 4, weight=0) 
        tk.Grid.rowconfigure(   research_window, 5, weight=0) 
        tk.Grid.columnconfigure(research_window, 0, weight=0)  
        tk.Grid.columnconfigure(research_window, 0, weight=0)  
        
        # Place the widgets in the frame
        self.lbl_title.grid( row=0, column=0, padx=50,     pady=(50,5), sticky='nesw', ipady=2, ipadx=5, columnspan=2)
        self.lbl_knob.grid(  row=1, column=0, padx=(50,5), pady=(5,5),  sticky='nesw', ipady=2, ipadx=5)
        self.mnu_knob.grid(  row=1, column=1, padx=(5,50), pady=(5,5),  sticky='nesw', ipady=2, ipadx=5)
        self.lbl_key.grid(   row=2, column=0, padx=(50,5), pady=(5,50),  sticky='nesw', ipady=2, ipadx=5)
        self.mnu_key.grid(   row=2, column=1, padx=(5,50), pady=(5,50), sticky='nesw', ipady=2, ipadx=5)
        self.lbl_title2.grid(row=3, column=0, padx=50,     pady=(5,5), sticky='nesw', ipady=2, ipadx=5, columnspan=2)
        self.lbl_knob2.grid( row=4, column=0, padx=(50,5), pady=(5,5),  sticky='nesw', ipady=2, ipadx=5)
        self.mnu_knob2.grid( row=4, column=1, padx=(5,50), pady=(5,5),  sticky='nesw', ipady=2, ipadx=5)
        self.lbl_key2.grid(  row=5, column=0, padx=(50,5), pady=(5,50),  sticky='nesw', ipady=2, ipadx=5)
        self.mnu_key2.grid(  row=5, column=1, padx=(5,50), pady=(5,50), sticky='nesw', ipady=2, ipadx=5)

        # Closing event
        def on_closing(*_):
            # Save the stella knob and key
            self.experiment_knob  = self.var_knob.get()
            self.experiment_key   = self.var_key.get() 
            self.experiment_knob2 = self.var_knob2.get()
            self.experiment_key2  = self.var_key2.get() 
            # Create a new research object
            if len(self.tab.input_files) != 0: 
                self.create_researchObject()
                update_GUI(self.root)
            # Close the window
            research_window.destroy()
        research_window.protocol("WM_DELETE_WINDOW", on_closing)
        research_window.deiconify()
        if True: return 
    
#====================
# Research object
#====================

    def create_researchObject(self):

        # Gather the arguments of the research object
        self.root.research_arguments = {
            # To create the simulations we need their location  
            "folders" : None, "input_files" : self.tab.input_files, 
            # To group by experiments and simulations we need the following data
            "variables" : self.number_variedVariables,
            "knob1" : self.experiment_knob, "key1" : self.experiment_key, 
            "knob2" : self.experiment_knob2, "key2" : self.experiment_key2,
            "ignoreResolution" : self.ignore_resolution, "folderIsExperiment" : self.folderIsExperiment, 
            "resolutionScan" : self.resolutionScan }
                
        # For the selected files, create a Research object
        self.root.Research = create_research(**self.root.research_arguments)
        
        # Save some information from the plots    
        try: data = self.root.Research.data
        except: data = {}
        self.root.Research.data = data
        
        # Print some information to the command prompt 
        self.root.Research.print_research() 
        
    #---------------------------
    def save_research(self):
        
        # Choose files and start selection in the standard run folder. 
        title           = "Save pickle file"
        filetypes       = (("in files","*.pickle"),("all files","*.*"))  
        selected_file   = tk.filedialog.asksaveasfile(initialdir=self.initialdir,title=title,filetypes=filetypes)

        # Only continue if a file was selected
        if selected_file!=None and selected_file != "DEFAULT" and self.tab.input_files != []:

            # Create a new research object 
            self.research_name = selected_file.name.replace(" ", "_").replace(".pickle", "")
            self.research_path = "/".join(self.research_name.split("/")[:-1])
            self.research_name = self.research_name.split("/")[-1]
            self.create_researchObject()   
            
            # Save it with pickle
            pickle_path = self.research_path+"/"+str(self.research_name+".pickle")
            pickly_file = open(pickle_path, "wb")
            pickle.dump(self.root.Research, pickly_file)            
            pickly_file.close()
            
            self.research_name = "DEFAULT"
            self.initialdir = self.research_path  
            if pathlib.Path(self.initialdir)!=pathlib.Path(CONFIG['PATHS']['GUI_Pickles']):
                mpl.rcParams["savefig.directory"] = self.initialdir
        
    #-----------------------
    def load_research(self):
        
        # Choose files and start selection in the standard run folder.
        title           = "Select pickle file"
        filetypes       = (("in files","*.pickle"),("all files","*.*"))  
        selected_files  = tk.filedialog.askopenfilenames(initialdir=self.initialdir,title=title,filetypes=filetypes)

        # Only continue if a file was selected
        if len(selected_files) > 0 and selected_files!=None:
            
            # Remember which folder we were loading 
            self.initialdir = "/".join(selected_files[0].split("/")[:-1]) 
            mpl.rcParams["savefig.directory"] = self.initialdir
            
            # Remove the input_files that are displayed now
            self.root.TabSelectedFiles.class_simulations.clear_simulations()
        
            # Open the file with pickle
            pickly_file = open(selected_files[0], 'rb')
            self.root.Research = pickle.load(pickly_file)            
            pickly_file.close()
            
            # Display the current input_files
            for experiment in self.root.Research.experiments:
                for simulation in experiment.simulations:
                    if simulation.nonlinear: input_files = [simulation.input_file]
                    if simulation.linear:    input_files = [mode.input_file for mode in simulation.modes]
                    self.tab.input_files += input_files 
            self.root.TabSelectedFiles.class_simulations.update_treeView()
            
            # Update the resolution parameter
            self.ignore_resolution = self.root.Research.creationDetails.ignoreResolution 
            try: self.folderIsExperiment = self.root.Research.creationDetails.folderIsExperiment 
            except: self.folderIsExperiment = False
            try: self.resolutionScan = self.root.Research.creationDetails.resolutionScan 
            except: self.resolutionScan = False
            
            # Update the varied parameters button 
            self.number_variedVariables = self.root.Research.creationDetails.variables 
            if self.number_variedVariables==-1: x = "No grouping"      
            if self.number_variedVariables==0:  x = "0 varied parameters" 
            if self.number_variedVariables==1:  x = "1 varied parameter"  
            if self.number_variedVariables==2:  x = "2 varied parameters"
            if self.number_variedVariables==3:  x = "3 varied parameters" 
            if self.number_variedVariables==4:  x = "4 varied parameters"
            if self.number_variedVariables==5:  x = "5 varied parameters"
            if self.number_variedVariables==10: x = "10 varied parameters"
            self.button_variedParam.config(text=x)
            
            # Remember the stella kob and key
            self.experiment_knob = self.root.Research.creationDetails.knob1
            self.experiment_key  = self.root.Research.creationDetails.key1
            self.experiment_knob2 = self.root.Research.creationDetails.knob2
            self.experiment_key2  = self.root.Research.creationDetails.key2
            
            # Update the GUI   
            update_GUI(self.root)
            return 
        
        # Prevent the collapsing of the header comment on the next line
        if True: return
        
#====================
# Update canvas
#====================

    def clear_simulations(self):
        ''' Remove all simulations from the research frame. '''
        
        # Remove all the items from the tree view
        self.tree.delete(*self.tree.get_children())
        self.tree_experiments = []
        
        # Remove the Research object
        self.root.Research = type('Dummy', (object,), {'content':{}})()
        self.root.Research.experiments = []
        return

    #----------------------------------
    def update_treeView(self):
        '''
        Update the scrollable canvas which displays the selected input files, 
        the selection is sorted by the parents folders. First make a frame
        and title label for each folder, next add the label, cross button and 
        checkmark widget for each input file or overwrite the existing labels
        with the correct data, overwriting reduces flickering of the GUI.
        
        Attributes
        ----------
        input_files: dict[folder][input_file] = tk.IntVar()
            Go thorough the simulations to display them on the GUI
            
        widgets : dict['folder']['title', 'label', 'button', 'check'] = list of widgets; dict['folder']['var_lbl'] = list of tk.StringVar
            Stores the newly created widgets for the simulations 
            
        frame_folders: dict['folder'] = ttk.Frame(master=self.scrollableCanvas.scrollable_frame)
            Stores the newly created widgets for the folders
        '''

        # Only update when there are input files
        if len(self.tab.input_files) != 0:
            
            # Remove all the items from the tree view
            self.tree.delete(*self.tree.get_children())
            self.tree_experiments = []
            
            # Create the research object
            self.create_researchObject()
    
            # Iterate over the folders
            for experiment in self.root.Research.experiments:

                # Check whether the experiment is already in the treeview
                # If it is not, add the experiment to the treeview
                if not self.tree.exists(experiment.id):
                    contents = [self.tree.item(child)["text"] for child in self.tree.get_children("")]
                    self.tree_experiments.append(self.tree.insert(parent="", index=bisect(contents, experiment.id), iid=experiment.id, text="Experiment: "+experiment.id, values=("","")))
                    self.tree.item(experiment.id, open=True, tags=['bold'])
                    self.tree.tag_configure("bold", font=(None, 11, 'bold'))
                    
                # Add the simulations to the experiment
                for tree_experiment in self.tree_experiments:
                    self.tree.selection_set(tree_experiment)  
                    experiment_iid = self.tree.selection()[0]
                    if experiment_iid==experiment.id:
                        for simulation in experiment.simulations:
                            # If the simulation is not in the treeview, add it
                            if not self.tree.exists(simulation.id): 
                                i = experiment.simulations.index(simulation)
                                text = experiment.variedValues[i].replace("\\","").replace("$","").replace(",", " ")
                                contents = [self.tree.item(child)["text"] for child in self.tree.get_children(experiment.id)]
                                self.tree.insert(tree_experiment, bisect(contents, simulation.id), iid=simulation.id, text=text, values=(" "))           
                            # If the input file is already in the treeview, make sure it has the correct columns
                            if self.tree.exists(simulation.id):
                                i = experiment.simulations.index(simulation)
                                text = experiment.variedValues[i].replace("\\","").replace("$","").replace(",", " ")
                                self.tree.item(simulation.id, text=text, values=(" "))

            # Remove focus from items
            for item in self.tree.selection():
                self.tree.selection_remove(item)

