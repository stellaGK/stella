
import tkinter as tk
from tkinter import ttk  
from stellapy.GUI.utils.display_information import display_information
        
################################################################################
#         CLASS TO SELECT THE EXPERIMENTS/SIMULATIONS TO BE PLOTTED
################################################################################
   
class Simulations: 
    ''' This frame allows one to choose which experiments and simulations 
    will be plotted, as well as which selection of (kx,ky) modes. '''
    
    def __init__(self, parent, frame):
        
        # Make the parents available
        self.tab = parent
        self.root = parent.root
        self.style = self.tab.root.awthemes  
        
        # Attributes
        self.experiment = None
        self.simulation = None
        self.experiment_id = "All experiments"
        self.simulation_id = "All simulations"
        
        # Option lists
        self.options_kx = [0.0]
        self.options_ky = [0.0, 100]
        self.options_exp = ["First plot data"]
        self.options_sim = ["First plot data"]
        self.options_simsids = []
        
        # Variables for the experiment/simulation/kxmin/kxmax/kymin/kymax
        self.var_experiment = tk.StringVar(value=self.options_exp[0])
        self.var_simulation = tk.StringVar(value=self.options_sim[0])
        self.var_kxmin = tk.StringVar(value=self.options_kx[0])
        self.var_kxmax = tk.StringVar(value=self.options_kx[-1])
        self.var_kymin = tk.StringVar(value=self.options_ky[0])
        self.var_kymax = tk.StringVar(value=self.options_ky[-1])
        
        # Put a label in front of the menus to know what the dropdown menu is
        self.lbl_exp   = ttk.Label(frame, text="Experiment:")
        self.lbl_sim   = ttk.Label(frame, text="Simulation:")
        self.lbl_kxmin = ttk.Label(frame, text="kx min:", style='opt_sign.TLabel') 
        self.lbl_kxmax = ttk.Label(frame, text=" kx max:", style='opt_sign.TLabel')
        self.lbl_kymin = ttk.Label(frame, text="ky min:", style='opt_sign.TLabel')
        self.lbl_kymax = ttk.Label(frame, text=" ky max:", style='opt_sign.TLabel')
        
        # Dropwon menus for the experiments, simulations, kxmin, kxmax, kymin and kymax
        self.mnu_exp   = ttk.OptionMenu(frame, self.var_experiment, self.options_exp[0], *self.options_exp , style='option.TMenubutton')
        self.mnu_sim   = ttk.OptionMenu(frame, self.var_simulation, self.options_sim[0], *self.options_sim, style='option.TMenubutton')
        self.mnu_kxmin = ttk.OptionMenu(frame, self.var_kxmin, self.options_kx[0], *self.options_kx, style='option.TMenubutton')
        self.mnu_kxmax = ttk.OptionMenu(frame, self.var_kxmax, self.options_kx[0], *self.options_kx, style='option.TMenubutton')
        self.mnu_kymin = ttk.OptionMenu(frame, self.var_kymin, self.options_ky[0], *self.options_ky, style='option.TMenubutton')
        self.mnu_kymax = ttk.OptionMenu(frame, self.var_kymax, self.options_ky[1], *self.options_ky, style='option.TMenubutton')

        # Link a function to a change of the dropdown options
        self.var_experiment.trace_id = self.var_experiment.trace('w', self.change_plottedExperiment) 
        self.var_simulation.trace_id = self.var_simulation.trace('w', self.change_plottedSimulation) 
        self.var_kxmin.trace_id = self.var_kxmin.trace('w', self.change_plottedModes)  
        self.var_kxmax.trace_id = self.var_kxmax.trace('w', self.change_plottedModes)  
        self.var_kymin.trace_id = self.var_kymin.trace('w', self.change_plottedModes)  
        self.var_kymax.trace_id = self.var_kymax.trace('w', self.change_plottedModes)  
        
        # Make the menus pretty
        self.mnu_exp["menu"].config(**self.style['menu']); self.mnu_exp.config(width=20)
        self.mnu_sim["menu"].config(**self.style['menu']); self.mnu_sim.config(width=20)       
        self.mnu_kxmin["menu"].config(**self.style['menu']); self.mnu_kxmin.config(width=4)
        self.mnu_kxmax["menu"].config(**self.style['menu']); self.mnu_kxmax.config(width=4)
        self.mnu_kymin["menu"].config(**self.style['menu']); self.mnu_kymin.config(width=4)
        self.mnu_kymax["menu"].config(**self.style['menu']); self.mnu_kymax.config(width=4)

        # Configure the rows/columns of the <frame>
        tk.Grid.columnconfigure(frame, 0, weight=0) 
        tk.Grid.columnconfigure(frame, 1, weight=1, uniform="2") 
        tk.Grid.columnconfigure(frame, 2, weight=0) 
        tk.Grid.columnconfigure(frame, 3, weight=1, uniform="2") 
        
        # Place the widgets in the <frame>
        self.lbl_exp.grid(  row=0, column=0, stick='NSEW', padx=(0,0), pady=(2,2), columnspan=2)
        self.mnu_exp.grid(  row=0, column=2, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=1, columnspan=2)
        self.lbl_sim.grid(  row=1, column=0, stick='NSEW', padx=(0,0), pady=(2,2), columnspan=2)
        self.mnu_sim.grid(  row=1, column=2, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=2, ipady=0, columnspan=2)
        self.lbl_kxmin.grid(row=2, column=0, stick='W',    padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_kxmin.grid(row=2, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.lbl_kxmax.grid(row=2, column=2, stick='W',    padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_kxmax.grid(row=2, column=3, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.lbl_kymin.grid(row=3, column=0, stick='W',    padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_kymin.grid(row=3, column=1, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.lbl_kymax.grid(row=3, column=2, stick='W',    padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        self.mnu_kymax.grid(row=3, column=3, stick='NSEW', padx=(0,0), pady=(2,2), ipadx=1, ipady=0)
        return

################################################################################
#                                     METHODS
################################################################################
 
    def display_plottedModesAndExperiments(self):
        ''' When the research is changed, we need to display the correct 
        experiments/simulations/modes. '''
        self.display_plottedExperiments() 
        self.display_plottedSimulations() 
        self.display_plottedModes()  
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
        for experiment in self.root.Research.experiments: sims = [ s.id for s in experiment.simulations ]  
        self.options_simids = ["All simulations"] + sims
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
    def display_plottedModes(self):
        ''' For the current simulation, show the range of kx/ky modes that can be plotted, this allows
        the user to set a specific k_min and k_max to display only a selection of modes. '''
        # Don't plot when changing <var_kmin> which executes <change_plottedModes>
        self.var_kxmin.trace_vdelete("w", self.var_kxmin.trace_id)
        self.var_kxmax.trace_vdelete("w", self.var_kxmax.trace_id)
        self.var_kymin.trace_vdelete("w", self.var_kymin.trace_id)
        self.var_kymax.trace_vdelete("w", self.var_kymax.trace_id)
        
        # Find the possible kx/ky values in the simulation/experiment
        if self.simulation_id == "All simulations":
            self.options_kx = [0] + self.experiment.vec.kx + [100]
            self.options_ky = [0] + self.experiment.vec.ky + [100]
        if self.simulation_id != "All simulations":
            self.options_kx = [0] + self.simulation.vec.kx + [100]
            self.options_ky = [0] + self.simulation.vec.ky + [100]
            
        # Limit the choices to 25
        while(len(self.options_kx)>25): self.options_kx = self.options_kx[::2]
        while(len(self.options_ky)>25): self.options_ky = self.options_ky[::2]
            
        # Make sure the dropdown menus display the correct kx/ky values    
        self.mnu_kxmin['menu'].delete(0, 'end');  self.mnu_kxmax['menu'].delete(0, 'end')
        self.mnu_kymin['menu'].delete(0, 'end');  self.mnu_kymax['menu'].delete(0, 'end')
        for kx in self.options_kx:
            self.mnu_kxmin['menu'].add_command(label=round(kx,2), command=tk._setit(self.var_kxmin, kx))
            self.mnu_kxmax['menu'].add_command(label=round(kx,2), command=tk._setit(self.var_kxmax, kx))
        for ky in self.options_ky:
            self.mnu_kymin['menu'].add_command(label=round(ky,2), command=tk._setit(self.var_kymin, ky))
            self.mnu_kymax['menu'].add_command(label=round(ky,2), command=tk._setit(self.var_kymax, ky))

        # Set the selected kx_min/kx_max/ky_min/ky_max to the minimum and maxim
        if self.tab.PlotLinearSpectrum.initiated_canvas: 
            kx_min = 0;     self.var_kxmin.set(kx_min)
            kx_max = 100;   self.var_kxmax.set(kx_max)  
            ky_min = 0;     self.var_kymin.set(ky_min)
            ky_max = 100;   self.var_kymax.set(ky_max) 
        
        # Attach the selected modes to the graph objects for the plotting
        self.change_plottedModes(calledFromGui=False)
        
        # Turn plotting back on when changing <var_simulation> in the GUI
        self.var_kxmin.trace_id = self.var_kxmin.trace('w', self.change_plottedModes) 
        self.var_kxmax.trace_id = self.var_kxmax.trace('w', self.change_plottedModes) 
        self.var_kymin.trace_id = self.var_kymin.trace('w', self.change_plottedModes) 
        self.var_kymax.trace_id = self.var_kymax.trace('w', self.change_plottedModes)
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
        if calledFromGui: self.display_plottedModes() 
        
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

        # Update the possible kx_min/kx_max/ky_min/ky_max  
        if calledFromGui: self.display_plottedModes()

        # Reset the graph classes and plot the graph
        if calledFromGui: self.tab.reset_axes(); self.tab.plot()
        return 
    
    #----------------------------
    def change_plottedModes(self, calledFromGui=True, *_): 
        ''' Tell the graph class which modes (kx,ky) to plot.''' 
        
        # Get the kx and ky ranges set in the GUI 
        kxmin = float(self.var_kxmin.get())
        kxmax = float(self.var_kxmax.get())
        kymin = float(self.var_kymin.get())
        kymax = float(self.var_kymax.get()) 

        # Tell the <Plots> which modes to plot
        for Plot in self.tab.Plots:
            if Plot.initiated_canvas: 
                if Plot.identifier!="TabLinear:PlotVelocityDistribution":
                    Plot.Axis.kx_range = [ kxmin, kxmax ]   
                    Plot.Axis.ky_range = [ kymin, kymax ] 
                if Plot.identifier=="TabLinear:PlotVelocityDistribution":  
                    Plot.Axis1.kx_range = [ kxmin, kxmax ]   
                    Plot.Axis1.ky_range = [ kymin, kymax ]   
                    Plot.Axis2.kx_range = [ kxmin, kxmax ]   
                    Plot.Axis2.ky_range = [ kymin, kymax ]    
            
        # Replot when we changed the (kx,ky) selection 
        if calledFromGui: self.tab.reset_axes(); self.tab.plot()
        return

    #----------------------------
    def change_plottedExperimentWithoutPlotting(self, *_):
        self.change_plottedExperiment(calledFromGui=False)
        return

    #----------------------------
    def change_plottedSimulationWithoutPlotting(self, *_):
        self.change_plottedSimulation(calledFromGui=False)
        return
         

