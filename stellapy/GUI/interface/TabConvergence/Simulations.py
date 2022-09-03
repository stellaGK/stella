 
import tkinter as tk
from tkinter import ttk  
from stellapy.GUI.widgets import Progress
from stellapy.GUI.utils.display_information import display_information
from stellapy.GUI.widgets import PAD_TITLE2, PAD_LABEL2 

################################################################################
#         CLASS TO SELECT THE EXPERIMENTS/SIMULATIONS TO BE PLOTTED
################################################################################
   
class Simulations: 
    ''' This frame allows one to choose which experiments and simulations 
    will be plotted, as well as which selection of (kx,ky) modes. Additionaly, 
    it controls which graph is currently plotted. The posibilities are:
        (1) Omega versus time
        (2) Gamma versus time
        (3) Phi**2 versus time
        (4) Re(Phi) versus zeta
        (5) Im(Phi) versus zeta
        (6) Phi**2 versus zeta  '''
    
    def __init__(self, parent, frame):
        ''' Create the frame <subframe_simulations> with the dropdown menus for
        the experiments, simulations and (kx,ky) ranges. Create the frame 
        <subframe_plots> that switches between the possible graphs. Create the
        frame <subframe_modes> that allows to choose to plot the stable, the 
        unstable, the unstable converged and the unstabled unconverged modes.
        Create the frame that gives some information about the seperate modes. '''

        # Make the parents available
        self.tab = parent
        self.root = parent.root
        self.style = self.tab.root.awthemes  

        # Attributes
        self.experiment = None
        self.simulation = None
        self.experiment_id = "All experiments"
        self.simulation_id = "All simulations"
        self.options_kx = [0.0, 100]
        self.options_ky = [0.0, 100]
        self.options_exp = ["First plot data"]
        self.options_sim = ["First plot data"]
        self.options_simids = []
        self.modes_id = "unstable"
        
        # Configure the rows/columns of the grid on <frame>
        tk.Grid.rowconfigure(   frame, 0, weight=0) 
        tk.Grid.rowconfigure(   frame, 1, weight=0) 
        tk.Grid.rowconfigure(   frame, 2, weight=0) 
        tk.Grid.rowconfigure(   frame, 3, weight=1)   
        tk.Grid.columnconfigure(frame, 0, weight=1) 
        
        # Create the <LabelFrames> on <frame>
        self.subframe_simulations = ttk.LabelFrame(frame, text="   Simulations  ", **self.style['labelframe2'])
        self.subframe_plots       = ttk.LabelFrame(frame, text="   Plot  ",        **self.style['labelframe2'])
        self.subframe_modes       = ttk.LabelFrame(frame, text="   Modes   ",      **self.style['labelframe2'])
        self.subframe_convergence = ttk.LabelFrame(frame, text="   Convergence   ",**self.style['labelframe2'])
        self.frame_progress       = ttk.Frame(frame)
        
        # Add the <LabelFrames> to the <frame>
        self.subframe_simulations.grid( row=0, column=0, padx=(5,20), pady=(0,10),  stick='NSEW')
        self.subframe_plots.grid(       row=1, column=0, padx=(5,20), pady=(10,10), stick='NSEW')
        self.subframe_modes.grid(       row=2, column=0, padx=(5,20), pady=(10,10), stick='NSEW')
        self.subframe_convergence.grid( row=3, column=0, padx=(5,20), pady=(10,10), stick='NSEW')
        self.frame_progress.grid(       row=4, column=0, padx=(5,20), pady=(10,0),  stick='NSEW')

        # Add the <Progress> object to <frame_progress>
        self.tab.Progress = Progress(self.root, self.frame_progress, anchor="fill", length=300)
        self.tab.Progress.move(0,"Please select simulations.")      
        
        # Fill the four subframes with widgets
        self.initiate_simulationsFrame(self.subframe_simulations)
        self.initiate_plotFrame(self.subframe_plots)
        self.initiate_modesFrame(self.subframe_modes)
        self.initiate_convergenceFrame(self.subframe_convergence)
        
    #----------- Simulations -------------
    def initiate_simulationsFrame(self, frame):
        ''' Fill the frame with the dropdown menus for the experiments, 
        simulations and (kx,ky) ranges. '''
        
        # Variables for the experiment/simulation/kxmin/kxmax/kymin/kymax
        self.var_experiment = tk.StringVar(value=self.options_exp [0])
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
        self.mnu_exp["menu"].config(**self.style['menu']); self.mnu_exp.config(width=15)
        self.mnu_sim["menu"].config(**self.style['menu']); self.mnu_sim.config(width=15)       
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
    
    #------------- Plots -------------    
    def initiate_plotFrame(self, frame):
        ''' Fill the frame <subframe_plots> that switches between the graphs: 
                (1) Omega versus time
                (2) Gamma versus time
                (3) Phi**2 versus time
                (4) Re(Phi) versus zeta
                (5) Im(Phi) versus zeta
                (6) Phi**2 versus zeta '''
        
        # Set the default value of the plot to "gamma versus time"
        self.var_plot  = tk.IntVar(value=2) 
        
        # Add the possible plots in two sections: "Time convergence" and "Space convergence"
        self.lbl_time  = ttk.Label(frame, text="Convergence in time", style='prefTitle.TLabel')
        self.rbn_omega = ttk.Radiobutton(frame, text='  Frequency versus time')
        self.rbn_gamma = ttk.Radiobutton(frame, text='  Growth rate versus time')
        self.rbn_phi2t = ttk.Radiobutton(frame, text='  Potential squared versus time')
        self.lbl_space = ttk.Label(frame, text="Confinement along the fluxtube", style='prefTitle.TLabel')
        self.rbn_real  = ttk.Radiobutton(frame, text='  Real potential versus z')
        self.rbn_imag  = ttk.Radiobutton(frame, text='  Imaginary potential versus z',)
        self.rbn_phi2  = ttk.Radiobutton(frame, text='  Potential squared versus z',)
        
        # Add the values and commands to the radiobuttons
        self.rbn_omega.config(value=1, variable=self.var_plot, command=lambda: self.change_plot())
        self.rbn_gamma.config(value=2, variable=self.var_plot, command=lambda: self.change_plot())
        self.rbn_phi2t.config(value=3, variable=self.var_plot, command=lambda: self.change_plot())
        self.rbn_real.config( value=4, variable=self.var_plot, command=lambda: self.change_plot())
        self.rbn_imag.config( value=5, variable=self.var_plot, command=lambda: self.change_plot())
        self.rbn_phi2.config( value=6, variable=self.var_plot, command=lambda: self.change_plot())
        
        # Configure the rows/columns of the <frame>
        tk.Grid.columnconfigure(frame, 0, weight=1) 
        
        # Add the options to the frame
        self.lbl_time.grid( row=0, column=0, **PAD_TITLE2)
        self.rbn_omega.grid(row=1, column=0, **PAD_LABEL2)
        self.rbn_gamma.grid(row=2, column=0, **PAD_LABEL2)
        self.rbn_phi2t.grid(row=3, column=0, **PAD_LABEL2)
        
        # Add the options to the frame
        self.lbl_space.grid(row=4, column=0, **PAD_TITLE2)
        self.rbn_real.grid( row=5, column=0, **PAD_LABEL2)
        self.rbn_imag.grid( row=6, column=0, **PAD_LABEL2)
        self.rbn_phi2.grid( row=7, column=0, **PAD_LABEL2)
        

    #-------------- Modes -------------
    def initiate_modesFrame(self, frame):
        ''' Create the frame <subframe_modes> that allows to choose to plot the stable,
        the unstable, the unstable converged and the unstabled unconverged modes. '''
        
        # Set the default value of the plotted modes 
        self.var_stab = tk.IntVar(value=2)

        # Stability and convergence
        self.rbn_stable      = ttk.Radiobutton(frame, text='  Stable modes')
        self.rbn_unstable    = ttk.Radiobutton(frame, text='  Unstable modes')
        self.rbn_converged   = ttk.Radiobutton(frame, text='  Unstable modes (converged)')
        self.rbn_unconverged = ttk.Radiobutton(frame, text='  Unstable modes (unconverged)')
        
        # Add the values and commands to the radiobuttons
        self.rbn_stable.config(     value=1, variable=self.var_stab, command=self.change_plottedSelectionModes)
        self.rbn_unstable.config(   value=2, variable=self.var_stab, command=self.change_plottedSelectionModes)
        self.rbn_converged.config(  value=3, variable=self.var_stab, command=self.change_plottedSelectionModes)
        self.rbn_unconverged.config(value=4, variable=self.var_stab, command=self.change_plottedSelectionModes)
        
        # Configure the rows/columns of the <frame>
        tk.Grid.columnconfigure(frame, 0, weight=1) 
        tk.Grid.columnconfigure(frame, 1, weight=0) 
        
        # Add the mode options to the frame
        self.rbn_stable.grid(     row=4, column=0, **PAD_LABEL2, columnspan=2)
        self.rbn_unstable.grid(   row=5, column=0, **PAD_LABEL2, columnspan=2)
        self.rbn_unconverged.grid(row=6, column=0, **PAD_LABEL2, columnspan=2)
        self.rbn_converged.grid(  row=7, column=0, **PAD_LABEL2, columnspan=2)


    #---------- Convergence ---------------        
    def initiate_convergenceFrame(self, frame):
        ''' Create the frame that gives some information about the seperate modes. '''
        
        # Some information about the convergence of the modes will be displayed:
        self.txt_convergence1 = tk.StringVar(value="")
        self.txt_convergence2 = tk.StringVar(value="No simulations are selected.") 
        self.txt_convergence3 = tk.StringVar(value="")      
        self.txt_convergence4 = tk.StringVar(value="")
        
        # Create 4 label widgets with the linked StringVar variables shown above
        self.lbl_convergence1 = ttk.Label(frame, textvariable=self.txt_convergence1, style='opt_paraC.TLabel')        
        self.lbl_convergence2 = ttk.Label(frame, textvariable=self.txt_convergence2, style='opt_valueC.TLabel') 
        self.lbl_convergence3 = ttk.Label(frame, textvariable=self.txt_convergence3, style='opt_paraC.TLabel')  
        self.lbl_convergence4 = ttk.Label(frame, textvariable=self.txt_convergence4, style='opt_valueC.TLabel')

        # Configure the rows/columns of the <frame>
        tk.Grid.rowconfigure(   frame, 0, weight=0) 
        tk.Grid.rowconfigure(   frame, 1, weight=0) 
        tk.Grid.rowconfigure(   frame, 2, weight=0)
        tk.Grid.rowconfigure(   frame, 3, weight=0)
        tk.Grid.columnconfigure(frame, 0, weight=1) 
        
        # Add the label widgets to the frame
        self.lbl_convergence1.grid(in_=frame, row=0, column=0, padx=(5,5), pady=(5,2),stick='NSEW')
        self.lbl_convergence2.grid(in_=frame, row=1, column=0, padx=(5,5), pady=(2,2),stick='NSEW')
        self.lbl_convergence3.grid(in_=frame, row=2, column=0, padx=(5,5), pady=(2,2),stick='NSEW')
        self.lbl_convergence4.grid(in_=frame, row=3, column=0, padx=(5,5), pady=(2,2),stick='NSEW')
        if True: return
  
################################################################################
#                                   METHODS
################################################################################

    def change_plot(self):
        ''' Make sure the correct canvas is visisble on the GUI, update the plot
        that is requested for the canvas, reset the axis and plot. '''
        if self.var_plot.get() in [1,2,3]:   
            self.tab.frame_canvas1.grid()
            self.tab.frame_canvas2.grid_remove() 
            self.tab.PlotTimeEvolution.Canvas.Axis.set_plot(self.var_plot.get())
        if self.var_plot.get() in [4,5,6]:   
            self.tab.frame_canvas2.grid()
            self.tab.frame_canvas1.grid_remove()
            self.tab.load_canvas2()
            self.tab.PlotParallelModeStructure.Canvas.Axis.set_plot(self.var_plot.get()-3)
        self.tab.reset_axes()
        self.tab.plot()
        return

    #----------------------------
    def change_plottedSelectionModes(self, *_):
        ''' Change the modes that will be plotted. '''
        if self.var_stab.get()==1: self.modes_id = "stable"
        if self.var_stab.get()==2: self.modes_id = "unstable"
        if self.var_stab.get()==3: self.modes_id = "converged"
        if self.var_stab.get()==4: self.modes_id = "unconverged" 
        self.tab.plot()
        return
    
    #----------------------------
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
            self.options_kx = self.experiment.vec.kx
            self.options_ky = self.experiment.vec.ky
        if self.simulation_id != "All simulations":
            self.options_kx = self.simulation.vec.kx
            self.options_ky = self.simulation.vec.ky
            
        # Make sure the dropdown menus display the correct kx/ky values    
        self.mnu_kxmin['menu'].delete(0, 'end');  self.mnu_kxmax['menu'].delete(0, 'end')
        self.mnu_kymin['menu'].delete(0, 'end');  self.mnu_kymax['menu'].delete(0, 'end')
        for kx in self.options_kx:
            self.mnu_kxmin['menu'].add_command(label=round(kx,2), command=tk._setit(self.var_kxmin, kx))
            self.mnu_kxmax['menu'].add_command(label=round(kx,2), command=tk._setit(self.var_kxmax, kx))
        for ky in self.options_ky:
            self.mnu_kymin['menu'].add_command(label=round(ky,2), command=tk._setit(self.var_kymin, ky))
            self.mnu_kymax['menu'].add_command(label=round(ky,2), command=tk._setit(self.var_kymax, ky)) 
            
        # Set the selected kx_min/kx_max/ky_min/ky_max to the minimum and maximum 
        if self.tab.PlotTimeEvolution.initiated_canvas:    
            kx_min = round(self.options_kx[0],2);     self.var_kxmin.set(kx_min)
            kx_max = round(self.options_kx[-1],2);    self.var_kxmax.set(kx_max) 
            ky_min = round(self.options_ky[0],2);     self.var_kymin.set(ky_min)
            ky_max = round(self.options_ky[-1],2);    self.var_kymax.set(ky_max)   
                    
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

        # Tell the graph class which modes to plot
        if self.tab.PlotTimeEvolution.initiated_canvas:
            self.tab.PlotTimeEvolution.Canvas.Axis.kx_range = [ kxmin, kxmax ] 
            self.tab.PlotTimeEvolution.Canvas.Axis.ky_range = [ kymin, kymax ] 
        if self.tab.PlotParallelModeStructure.initiated_canvas:
            self.tab.PlotParallelModeStructure.Canvas.Axis.kx_range = [ kxmin, kxmax ] 
            self.tab.PlotParallelModeStructure.Canvas.Axis.ky_range = [ kymin, kymax ] 
            
        # Replot when we changed the (kx,ky) selection 
        if calledFromGui: self.tab.reset_axes(); self.tab.plot()
        return
        
    #----------------------------
    def write_informationConvergence(self):
        ''' Write how many modes are unstable, and how many converged. ''' 
        
        # Make sure the elements are strings that are rounded to 2 digits
        vec_k            = [ str(round(mode.ky,2)) for mode in self.simulation.lineardata.allModes ]
        vec_kUnstable    = [ str(round(mode.ky,2)) for mode in self.simulation.lineardata.unstableModes ]  
        vec_kNotConverge = [ str(round(mode.ky,2)) for mode in self.simulation.lineardata.unconvergedModes ]
        vec_kStable      = [ str(round(mode.ky,2)) for mode in self.simulation.lineardata.stableModes ]
        
        # Print information if no modes are plotted
        if (self.var_stab.get()==1 and len(vec_kStable)==0) or (self.var_stab.get()!=1 and len(vec_kUnstable)==0):
            if self.var_stab.get()==1: # We want to see the stable modes, but there are none!
                self.txt_convergence1.set("From the "+str(len(vec_k))+" modes, all are unstable.")
            if self.var_stab.get()!=1: # We want to see the unstable modes, but there are none!
                self.txt_convergence1.set("From the "+str(len(vec_k))+" modes, all are stable.")
            self.txt_convergence2.set("There is nothing to plot.")
        
        # Overwrite the default message with information about the number of stable or unstable modes 
        else:  
            if len(vec_kStable) > 0 and self.var_stab.get()==1: # We want to see the stable modes
                txt_notconverge1 = " are stable:" if len(vec_kStable)>1 else " is stable:"
                txt_notconverge1 = "From the "+str(len(vec_k))+" modes, "+str(len(vec_kStable))+txt_notconverge1
                txt_notconverge2 = "[" + vec_kStable[0] + ", ..., " + vec_kStable[-1] + "]"
            if len(vec_kUnstable) > 0 and self.var_stab.get()!=1: # We want to see the unstable modes   
                txt_notconverge1 = " are unstable:" if len(vec_kUnstable)>1 else " is unstable:"
                txt_notconverge1 = "From the "+str(len(vec_k))+" modes, "+str(len(vec_kUnstable))+txt_notconverge1
                txt_notconverge2 = "[" + vec_kUnstable[0] + ", ..., " + vec_kUnstable[-1] + "]"
            self.txt_convergence1.set(txt_notconverge1)
            self.txt_convergence2.set(txt_notconverge2)
            
        # Print information about the number of unconverged modes
        if len(vec_kNotConverge) > 0:
            txt_notconverge3 = " modes have not converged:" if len(vec_kNotConverge)>1 else " mode has not converged:"
            txt_notconverge3 = "And " + str(len(vec_kNotConverge)) + txt_notconverge3
            txt_notconverge4 = "[" + "\n".join([", ".join(vec_kNotConverge[i:i+6]) for i in range(0,len(vec_kNotConverge),6)]) + "]"
            self.txt_convergence3.set(txt_notconverge3)
            self.txt_convergence4.set(txt_notconverge4)
        else:  
            self.txt_convergence3.set("All simulations have converged.")
            self.txt_convergence4.set("")
        return

    #----------------------------
    def change_plottedExperimentWithoutPlotting(self, *_):
        self.change_plottedExperiment(calledFromGui=False)

    #----------------------------
    def change_plottedSimulationWithoutPlotting(self, *_):
        self.change_plottedSimulation(calledFromGui=False)

 





