
import tkinter as tk
from tkinter import ttk   
        
#################################################################
#                  CLASS TO CHOOSE SIMULATIONS
#################################################################
   
class Simulations: 
    
    def __init__(self, parent, frame):
        ''' This options frame controls which experiments/simulations are plotted. '''
        
        # Make the parents available
        self.tab = parent
        self.root = parent.root
        self.frame = frame
        
        # Attributes
        self.experiment = None
        self.simulation = None
        self.experiment_id = "All experiments"
        self.simulation_id = "All simulations"
        self.y_quantity = "bmag"
        
        # Option lists
        self.options_experiment = ["First plot data"]
        self.options_simulation = ["First plot data"] 
        self.options_simulationsids = []
        self.options_quantity = ["alpha","zed", "zeta", "bmag", "gradpar", "gds2", "gds21", "gds22", "gds23", "gds24", "gbdrift", "cvdrift", "gbdrift0", "bmag_psi0"] 
        
        # Variables
        self.var_experiment = tk.StringVar(value=self.options_experiment[0])
        self.var_simulation = tk.StringVar(value=self.options_simulation[0]) 
        self.var_quantity   = tk.StringVar(value=self.options_quantity[0]) 
        
        # Configure the frame
        tk.Grid.columnconfigure(frame, 0, weight=1) 
        
        # Choice between the experiments
        self.lbl_experiment = ttk.Label(self.frame, text="Experiment: "); width=20
        self.mnu_experiment = ttk.OptionMenu(self.frame, self.var_experiment, self.options_experiment[0], *self.options_experiment, style='option.TMenubutton')
        self.mnu_experiment["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_experiment.config(width=width)
        self.var_experiment.trace('w', self.change_plottedExperiment) # link function to a change of the dropdown options
        
        # Choice between the simulations
        self.lbl_simulation = ttk.Label(self.frame, text="Simulation: ")
        self.mnu_simulation = ttk.OptionMenu(self.frame, self.var_simulation, self.options_simulation[0], *self.options_simulation, style='option.TMenubutton')
        self.mnu_simulation["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_simulation.config(width=width)
        self.var_simulation.trace('w', self.change_plottedSimulation) # link function to a change of the dropdown options

        # Choice between the potential; fluxes and moments
        self.lbl_quantity = ttk.Label(self.frame, text="Quantity: ")
        self.mnu_quantity = ttk.OptionMenu(frame, self.var_quantity, self.options_quantity[3], *self.options_quantity, style='option.TMenubutton')
        self.mnu_quantity["menu"].config(bg=self.root.color['bbg'], fg=self.root.color['fg'], activebackground=self.root.color['bg'], activeforeground=self.root.color['fg'])
        self.mnu_quantity.config(width=width)
        self.var_quantity.trace('w', self.change_yQuantity) 
        
        # Configure the subframe 
        tk.Grid.columnconfigure(frame, 1, weight=1)  
        
        # Place the widgets in the frame
        self.lbl_experiment.grid(row=0, column=0, stick='NSEW', padx=(0,0), pady=(2,2))
        self.mnu_experiment.grid(row=0, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=1)
        self.lbl_simulation.grid(row=1, column=0, stick='NSEW', padx=(0,0), pady=(2,2))
        self.mnu_simulation.grid(row=1, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=1) 
        self.lbl_quantity.grid(  row=3, column=0, stick='NSEW', padx=(0,0), pady=(2,2))
        self.mnu_quantity.grid(  row=3, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=1) 
        if True: return


#################################################################
#                          METHODS
#################################################################

    def change_yQuantity(self, plot=True, *args):
         
        # When changing the plots, reset the axis of the graph
        self.tab.Graph[self.tab.PlotVariablesVsZ.graph_id].load_defaults()
        self.tab.Graph[self.tab.PlotVariablesVsZ.graph_id].update_quantitiesAndKeys() 
         
        # Extract what needs to be plotted on the y-axis 
        self.y_quantity = self.var_quantity.get() 
        
        # Plot the graph
        if plot:
            if self.tab.frame_graph0.grid_info() != {}: self.tab.PlotVariablesVsZ.replot_graph()
        return     
        
    def change_plottedExperiment(self, *args):
        
        # Change the experiment that is plotted
        self.experiment_id = self.var_experiment.get()
    
        # Get a reference to the experiment
        self.experiment = None
        for experiment in self.root.Research.experiments:
            if experiment.id == self.experiment_id:
                self.experiment = experiment
        if self.experiment==None and self.experiment_id!="All experiments":
            self.experiment = self.root.Research.experiments[0]
            self.experiment_id = self.experiment.id
        if self.experiment==None and self.experiment_id=="All experiments":
            self.experiment = self.root.Research.experiments[0]
            
        # Reset the graph classes
        self.tab.Graph[0].load_defaults()
        self.tab.Graph[1].load_defaults() 
 
        # Replot the graphs
        #if self.tab.frame_graph0.grid_info() != {}: self.tab.PlotTimeEvolution.plot_graph()
        #if self.tab.frame_graph1.grid_info() != {}: self.tab.PlotParallelModeStructure.plot_graph()
        #if self.tab.frame_graph2.grid_info() != {}: self.tab.PlotSpectrumKxKy.plot_graph()
            
    #------------------------------------
    def change_plottedSimulation(self, *args):
        
        # Change the simulation that is plotted
        self.simulation_id = self.var_simulation.get()
        self.simulation_id = self.options_simulationsids[self.options_simulations.index(self.simulation_id)]
        
        # Get a reference to the simulation
        self.simulation = None
        for simulation in self.experiment.simulations:
            if simulation.id == self.simulation_id:
                self.simulation = simulation
        if self.simulation==None and self.simulation_id!="All simulations":
            self.simulation = self.experiment.simulations[0]
            self.simulation_id = self.simulation.id
        if self.simulation==None and self.simulation_id=="All simulations":
            self.simulation = self.experiment.simulations[0] 
        self.replot=True

        # Replot the graphs
        if self.tab.frame_graph0.grid_info() != {}: self.tab.PlotVariablesVsZ.plot_graph()
        return 

    #-------------------------------------------
    def display_plottedModesAndExperiments(self):
        ''' When the graph is plotted, adjust the menus on the GUI of the modes and experiments. '''
        
        # Get the experiment plotted by the function
        self.options_experiments = ["All experiments"] + [ e.id for e in self.root.Research.experiments ]
        if self.experiment==None:
            self.experiment = self.root.Research.experiments[0]
            self.var_experiment.set(self.options_experiments[0])        
        if self.experiment.id not in self.options_experiments:
            self.experiment = self.root.Research.experiments[0]
            self.var_experiment.set(self.options_experiments[0])      
        if self.experiment_id not in self.options_experiments:
            self.experiment = self.root.Research.experiments[0]
            self.var_experiment.set(self.options_experiments[0])  
        if self.experiment_id == "All experiments":
            self.experiment = self.root.Research.experiments[0]
        for experiment in self.root.Research.experiments:
            if experiment.id == self.experiment_id:
                self.experiment = experiment
                self.var_experiment.set(self.experiment.id)
            
        # Reset the options in the menu
        self.mnu_experiment['menu'].delete(0, 'end')
        for experiment in self.options_experiments:
            self.mnu_experiment['menu'].add_command(label=experiment, command=tk._setit(self.var_experiment, experiment))

        # Get the simulation plotted by the function
        self.options_simulationsids = ["All simulations"] + [ s.id for s in self.experiment.simulations ]  
        self.options_simulations = ["All simulations"] + [ s.id.split("__")[-1] for s in self.experiment.simulations ]  
        self.options_simulations = ["All simulations"] + [ v for v in self.experiment.variedValues ]  
        self.options_simulations = [ s.replace('\\', '').replace(',', '').replace('$', '') for s in self.options_simulations ] 
        if self.simulation==None:
            self.simulation = self.experiment.simulations[0]
            self.var_simulation.set(self.options_simulations[0])
        if self.simulation.id not in self.options_simulationsids:
            self.simulation = self.experiment.simulations[0]
            self.var_simulation.set(self.options_simulations[0])
        if self.simulation_id not in self.options_simulationsids:
            self.simulation = self.experiment.simulations[0]
            self.var_simulation.set(self.options_simulations[0])
        for simulation in self.experiment.simulations:
            if simulation.id == self.simulation_id:
                self.simulation = simulation
                self.var_simulation.set(self.options_simulations[1+self.experiment.simulations.index(simulation)])
        
        # Reset the options in the menu
        self.mnu_simulation['menu'].delete(0, 'end')
        for simulation in self.options_simulations:
            self.mnu_simulation['menu'].add_command(label=simulation, command=tk._setit(self.var_simulation, simulation))
        return
    

