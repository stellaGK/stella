
import tkinter as tk
from tkinter import ttk

#===============================================================================
#                   SHOW THE PROGRESS WHEN READING FILES
#===============================================================================

def show_progressWhenReadingFiles(self, message): 
    if self.Progress==None:         return
    if self.object=="Mode":         show_progressWhenReadingFilesForModes(self, message)
    if self.object=="Simulation":   show_progressWhenReadingFilesForSimulations(self, message)
    return

def show_progressWhenReadingFilesForModes(self, message):
    current_mode = self.simulation.modes.index(self.mode)
    total_modes = len(self.simulation.modes)
    position = current_mode/total_modes*100
    message = message + "("+str(current_mode)+"/"+str(total_modes)+")"
    self.Progress.move(position, message)
    return
    
def show_progressWhenReadingFilesForSimulations(self, message): 
    current_simulation = self.simulation.simulations.index(self.simulation)
    total_simulations = len(self.simulation.simulations)
    position = current_simulation/total_simulations*100
    message = message + "("+str(current_simulation)+"/"+str(total_simulations)+")"
    self.Progress.move(position, message) 
    return
        
#===============================================================================
#           SHOW LIVE PROGRESS ON THE BACKGROUND PROCESSES OF THE GUI
#===============================================================================

class Progress:
 
    def __init__(self, root, frame, anchor=None, length=None):

        # Safe the root and frame
        self.root = root  
        self.frame = frame
        
        # Message to show progress  
        self.txt_information = tk.StringVar()               
    
        # Only create if we have a frame, we can make a Progress bar without
        # a frame just to avoid having to initialize it first in every routine
        if self.frame!=None:

            # Load the current style
            self.style = ttk.Style(self.root)
            
            # Create progress bar
            self.bar = ttk.Progressbar(frame, orient = tk.HORIZONTAL, style="Labeledself", length=length) 
            self.bar.config(mode = 'determinate', maximum=100, value = 0)
            self.style.configure("Labeledself", text="No tasks are running.")
    
            # Add the widgets
            if anchor == None:
                tk.Grid.rowconfigure(   frame, 0, weight=0) # Progress bar, make as small as possible
                tk.Grid.columnconfigure(frame, 0, weight=1) # Elongate along y
                self.bar.grid(in_=frame, row=0, column=0, padx=2, pady=0, sticky='nesw', ipady=5)
            if anchor == "center":
                tk.Grid.rowconfigure(   frame, 0, weight=1) # Spacing above
                tk.Grid.rowconfigure(   frame, 1, weight=0) # Progress bar, make as small as possible
                tk.Grid.rowconfigure(   frame, 2, weight=1) # Spacing below
                tk.Grid.columnconfigure(frame, 0, weight=1) # Spacing left
                tk.Grid.columnconfigure(frame, 1, weight=0) # Progress bar, make as small as possible
                tk.Grid.columnconfigure(frame, 2, weight=1) # Spacing right
                self.bar.grid(in_=frame, row=1, column=1, padx=2, pady=0, sticky='nesw', ipady=5)
            if anchor == "fill":
                tk.Grid.rowconfigure(   frame, 0, weight=0) # Progress bar, make as small as possible
                tk.Grid.columnconfigure(frame, 0, weight=1) # Elongate along y
                self.bar.grid(in_=frame, row=0, column=0, padx=0, pady=0, sticky='nesw', ipady=5)

    #------------------------------
    def start(self, message):
        ''' Start the progress bar. '''
        if self.frame!=None:
            self.root.update_idletasks()
            self.bar['value'] = 0
            self.style.configure("Labeledself", text=message)
            self.root.update_idletasks()
        return

    #------------------------------
    def move(self, position, message):
        ''' Move the dark bar inside the progress bar. '''
        if self.frame!=None:
            self.root.update_idletasks()
            self.bar['value'] = position
            self.style.configure("Labeledself", text=message)
            self.root.update_idletasks()
        return
 
    #------------------------------
    def finish(self):
        ''' Finish the progress bar. '''
        if self.frame!=None:
            self.root.update_idletasks()
            self.bar['value'] = 0
            self.style.configure("Labeledself", text="No tasks are running.")
            self.root.update_idletasks()
        return

    #------------------------------
    def show_loadingCursor(self, status="start"):
        ''' Show a loading cursor when plotting and move the progress bar. '''
        if self.frame!=None:
            # While plotting show a progress bar turn the cursor into a waiting cursor.
            if status=="start": 
                self.start("Start plotting.")
                self.root.config(cursor="watch")
                self.root.update_idletasks() 
            # When finished show the canvas and change the cursor back.
            if status=="finished": 
                self.finish()
                self.root.config(cursor="")
            # When there are no simulations revert to the loading bar
            if status=="nothing": 
                self.start("Please select simulations.")
        return














