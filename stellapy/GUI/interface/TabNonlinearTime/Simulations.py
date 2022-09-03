
import tkinter as tk
from tkinter import ttk   
        
################################################################################
#         CLASS TO SELECT THE EXPERIMENTS/SIMULATIONS TO BE PLOTTED
################################################################################
# TODO: Also update the species option list dynamically

class Simulations: 
    ''' This frame allows one to choose which experiments/simulations are plotted. '''
    
    def __init__(self, parent, frame): 
        
        # Make the parents available
        self.tab = parent
        self.root = parent.root
        self.style = self.tab.root.awthemes  
        
        # Attributes
        self.species = [0]
        self.experiment = None
        self.simulation = None
        self.experiment_id = "All experiments"
        self.simulation_id = "All simulations"
        
        # Option lists
        self.options_exp = ["First plot data"]
        self.options_sim = ["First plot data"]
        self.options_species = ["Ions", "Electrons", "Impurity 1", "Impurity 2", "All species"]
        self.options_simsids = []
        
        # Variables for the experiment/simulation/species
        self.var_experiment = tk.StringVar(value=self.options_exp[0])
        self.var_simulation = tk.StringVar(value=self.options_sim[0])
        self.var_specie     = tk.StringVar(value=self.options_species[0])

        # Put a label in front of the menus to know what the dropdown menu is
        self.lbl_exp    = ttk.Label(frame, text="Experiment:")
        self.lbl_sim    = ttk.Label(frame, text="Simulation:")
        self.lbl_specie = ttk.Label(frame, text="Specie: ")

        # Dropdown menus for the experiments, simulations and species
        self.mnu_exp    = ttk.OptionMenu(frame, self.var_experiment, self.options_exp[0], *self.options_exp , style='option.TMenubutton')
        self.mnu_sim    = ttk.OptionMenu(frame, self.var_simulation, self.options_sim[0], *self.options_sim, style='option.TMenubutton')
        self.mnu_specie = ttk.OptionMenu(frame, self.var_specie, self.options_species[0], *self.options_species, style='option.TMenubutton')
        
        # Link a function to a change of the dropdown options
        self.var_experiment.trace_id = self.var_experiment.trace('w', self.change_plottedExperiment) 
        self.var_simulation.trace_id = self.var_simulation.trace('w', self.change_plottedSimulation) 
        self.var_specie.trace_id     = self.var_specie.trace('w', self.change_species) 
        
        # Make the menus pretty
        self.mnu_exp["menu"].config(**self.style['menu']); self.mnu_exp.config(width=20)
        self.mnu_sim["menu"].config(**self.style['menu']); self.mnu_sim.config(width=20) 
        self.mnu_specie["menu"].config(**self.style['menu']); self.mnu_sim.config(width=20)      

        # Button to display/modify the time frames for each simulation
        self.btn_timeFrame = ttk.Button(frame, text="Select time frames", width=30)
        self.btn_timeFrame.config(command=lambda: self.tab.DisplayTimeFrames.select_timeFrames())

        # Configure the rows/columns of the <frame>
        tk.Grid.columnconfigure(frame, 0, weight=1) 
        tk.Grid.columnconfigure(frame, 1, weight=1)  
        
        # Place the widgets in the <frame>
        self.lbl_exp.grid(      row=0, column=0, stick='NSEW', padx=(0,0), pady=(2,2))
        self.mnu_exp.grid(      row=0, column=1, stick='NSEW', padx=(0,0), pady=(2,2),  ipadx=2, ipady=1)
        self.lbl_sim.grid(      row=1, column=0, stick='NSEW', padx=(0,0), pady=(2,2))
        self.mnu_sim.grid(      row=1, column=1, stick='NSEW', padx=(0,0), pady=(2,2),  ipadx=2, ipady=1)
        self.lbl_specie.grid(   row=2, column=0, stick='NSEW', padx=(0,0), pady=(2,2))
        self.mnu_specie.grid(   row=2, column=1, stick='NSEW', padx=(0,0), pady=(2,2),  ipadx=2, ipady=1)
        self.btn_timeFrame.grid(row=3, column=0, stick='NSEW', padx=(0,0), pady=(10,2), ipadx=2, ipady=1, columnspan=2)
        if True: return

################################################################################
#                                     METHODS
################################################################################

    def display_plottedModesAndExperiments(self):
        ''' When the research is changed, we need to display the correct 
        experiments/simulations/modes. '''
        self.display_plottedExperiments() 
        self.display_plottedSimulations()  
        return

    #----------------------------
    def display_plottedExperiments(self):
        ''' The displayed experiment is choosen from a list of experiment id's 
        <self.options_exp>, this method puts the ids into a dropdown menu. When
        setting <self.var_experiment>, the method <change_plottedExperiment> is 
        called, which sets <self.experiment>. '''

        # Don't plot when changing <var_experiment> which executes <change_plottedExperiment>
        self.var_experiment.trace_vdelete("w", self.var_experiment.trace_id)
        self.var_experiment.trace_id = self.var_experiment.trace('w', self.change_plottedExperimentWithoutPlotting) 
        
        # Get the experiments that are included in the research
        self.options_exp = ["All experiments"] + [ e.id for e in self.root.Research.experiments ]
        
        # Make sure the dropdown menu displays the correct experiments
        self.mnu_exp['menu'].delete(0, 'end')
        for experiment_id in self.options_exp:
            self.mnu_exp['menu'].add_command(label=experiment_id, command=tk._setit(self.var_experiment, experiment_id))
        
        # If we want to show all experiments, put the first experiment in self.experiment. 
        if self.experiment_id == "All experiments":
            self.var_experiment.set("All experiments") 

        # If an experiment or experiment_id had already been set, make sure it is still a part 
        # of this research. If it is not, take again the first experiment in the research, 
        elif self.experiment.id not in self.options_exp: 
            self.var_experiment.set(self.options_exp[0])      
        elif self.experiment_id not in self.options_exp: 
            self.var_experiment.set(self.options_exp[0])  
            
        # If the experiment_id had already been set and is a part of the current research, 
        # make sure self.experiment is linked correctly. If we switched to a research with the same
        # experiments, the experiment id would still be valid but the experiment link would be wrong. 
        else: self.var_experiment.set(self.experiment_id)
        
        # Turn plotting back on when changing <var_experiment> in the GUI
        self.var_experiment.trace_vdelete("w", self.var_experiment.trace_id)
        self.var_experiment.trace_id = self.var_experiment.trace('w', self.change_plottedExperiment) 
        return
    
    #----------------------------
    def display_plottedSimulations(self):
        ''' The displayed simulation is choosen from a list of simulation id's 
        <self.options_sim>, this method puts the ids into a dropdown menu. When
        setting <self.var_simulation>, the method <change_plottedSimulation> is 
        called, which sets <self.simulation>. '''

        # Don't plot when changing <var_simulation> which executes <change_plottedsimulation>
        self.var_simulation.trace_vdelete("w", self.var_simulation.trace_id)
        self.var_simulation.trace_id = self.var_simulation.trace('w', self.change_plottedSimulationWithoutPlotting) 
        
        # Get the simulations that are included in the research, get a list of the simulation ids
        # to know which simulation to choose, but display the simulation marker labels in the dropdown
        # menu since these are usually easier to recognize than the simulaion ids
        self.options_simids = ["All simulations"] + [ s.id for s in self.experiment.simulations ]   
        self.options_sim = ["All simulations"] + [ v for v in self.experiment.variedValues ]  
        self.options_sim = [ s.replace('\\', '').replace(',', '').replace('$', '') for s in self.options_sim ] 
        
        # Make sure the dropdown menu displays the correct simulations
        self.mnu_sim['menu'].delete(0, 'end')
        for simulation in self.options_sim:
            self.mnu_sim['menu'].add_command(label=simulation, command=tk._setit(self.var_simulation, simulation))
        
        # If we want to show all simulations, put the first simulation in self.simulation. 
        if self.simulation_id == "All simulations":
            self.var_simulation.set(self.options_sim[0])
            
        # If a simulation or simulation_id had already been set, make sure it is still a part 
        # of this research. If it is not, take again the first simulation in the experiment.
        elif self.simulation.id not in self.options_simids:
            self.var_simulation.set(self.options_sim[0])
        elif self.simulation_id not in self.options_simids:
            self.var_simulation.set(self.options_sim[0])

        # If the simulation_id had already been set and is a part of the current experiment, 
        # make sure self.simulation is linked correctly. If we switched to a experiment with the same
        # simulations, the simulation id would still be valid but the simulation link would be wrong. 
        else: self.var_simulation.set(self.options_sim[self.options_simids.index(self.simulation_id)-1])
        
        # Turn plotting back on when changing <var_simulation> in the GUI
        self.var_simulation.trace_vdelete("w", self.var_simulation.trace_id)
        self.var_simulation.trace_id = self.var_simulation.trace('w', self.change_plottedSimulation) 
        return 
    
    #----------------------------
    def change_plottedExperiment(self, calledFromGui=True, *_):
        ''' When an experiment id is selected from the dropdown menu, link the 
        correct experiment and update the possible simulations and (kx,ky). '''
        
        # Get the experiment id that is selected from the dropdown menu
        self.experiment_id = self.var_experiment.get()
    
        # Based on the experiment id, set <self.experiment>
        if self.experiment_id=="All experiments":
            self.experiment = self.root.Research.experiments[0]
        if self.experiment_id!="All experiments":
            for experiment in self.root.Research.experiments:
                if experiment.id == self.experiment_id:
                    self.experiment = experiment; break
        
        # Update the dropdown menus of simulations and the possible kx_min/kx_max/ky_min/ky_max
        if calledFromGui: self.display_plottedSimulations() 
        
        # Reset the graph classes and plot the graph
        if calledFromGui: self.tab.reset_axes(); self.tab.plot()
        return 
        
    #----------------------------
    def change_plottedSimulation(self, calledFromGui=True, *_):
        ''' When a simulation id is selected from the dropdown menu, link the 
        correct experiment and update the possible (kx,ky). ''' 
        
        # Get the simulation id that is selected from the dropdown menu
        self.simulation_id = self.var_simulation.get()
        self.simulation_id = self.options_simids[self.options_sim.index(self.simulation_id)]
            
        # Based on the simulation id, set <self.simulation> 
        if self.simulation_id=="All simulations":
            self.simulation = self.experiment.simulations[0] 
        if self.simulation_id!="All simulations": 
            for simulation in self.experiment.simulations:
                if simulation.id == self.simulation_id:
                    self.simulation = simulation; break

        # Reset the graph classes and plot the graph
        if calledFromGui: self.tab.reset_axes(); self.tab.plot()
        return 
    
    #------------------------------------
    def change_species(self, *_):
        ''' Change the species that is plotted. '''
        if self.var_specie.get()=="All species": 
            dim_species = self.root.Research.experiments[0].simulations[0].dim_species
            self.species = [ i for i in range(dim_species)]
        if self.var_specie.get()!="All species":
            self.species = [int(self.options_species.index(self.var_specie.get()))]
        return
        
    #----------------------------
    def change_plottedExperimentWithoutPlotting(self, *_):
        self.change_plottedExperiment(calledFromGui=False)
        return

    #----------------------------
    def change_plottedSimulationWithoutPlotting(self, *_):
        self.change_plottedSimulation(calledFromGui=False)
        return
         

