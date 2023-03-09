  
import numpy as np  
from stellapy.plot.utils.labels.standardLabels import standardLabels
from stellapy.plot.utils.labels.standardTitles import standardTitles  
from stellapy.plot.utils.species.recognize_species import recognize_species
from stellapy.plot.utils.style.load_styleFigures import load_styleFigures
 
################################################################################
#                            PLOTTING DETAILS                                  #
################################################################################

class Plot():
    ''' Remember all the plotting toggles for easy access. '''
    
    def __init__(self, quantity=None):
        self.update_simulations()
        self.update_parameter()
        self.update_toggles()
        self.update_xyzData() 
        self.update_labels()
        self.update_figure()
        self.update_legend()
        self.update_modes() 
        self.update_axis()
        self.quantity = quantity
        return

################################################################################
#                                 METHODS                                      #
################################################################################

    def process_plottingVariables(self, research):
        load_styleFigures()
        self.get_plottedExperimentsSimulationsAndModes(research)
        self.get_xyzDimensions()
        self.get_xyzLabels(research)  
        self.get_title()
        return

    #------------------------------------
    def get_plottedExperimentsSimulationsAndModes(self, research): 
        
        # Determine the plotted experiments/simulation/modes
        from stellapy.simulations.utils.get_simulations import get_simulations
        from stellapy.simulations.utils.get_experiments import get_experiments 
        research.plotted_experiments = get_experiments(research, self.experiment_id)
        research.numberOfPlottedExperiments = len(research.plotted_experiments)
        research.numberOfPlottedSimulations = 0
        research.numberOfPlottedModes = 0
        for experiment in research.plotted_experiments:
            experiment.plotted_simulations = get_simulations(experiment, self.simulation_id)
            research.numberOfPlottedSimulations += len(experiment.plotted_simulations) 
        
        # Make the <plot> class avaliable from the <research> class
        research.plot = self
        return
    
    #------------------------------------
    def get_xyzDimensions(self): 
        
        # We don't always use the full plot package
        if self.x_quantity==None: return
            
        # Remove certain words for the labels 
        x_quantity1 = self.x_quantity.replace("_last","").replace("_avg","") if self.x_quantity!=None else self.x_quantity
        y_quantity1 = self.y_quantity.replace("_last","").replace("_avg","") if self.y_quantity!=None else self.y_quantity 
        z_quantity1 = self.z_quantity.replace("_last","").replace("_avg","") if self.z_quantity!=None else self.z_quantity
        
        # Define what is on the (x,y) or (x,y,z) axis for the labels
        if self.z_quantity==None: 
            self.xname = x_quantity1
            self.yname = y_quantity1
            self.zname = None
        if self.z_quantity!=None: 
            self.xname = x_quantity1
            self.yname = y_quantity1
            self.zname = z_quantity1
        
        # Remove certain words for the name of the data
        x_quantity2 = self.x_quantity.replace("/ky**2","").replace("_real","").replace("_imag","").replace("_zonal","").replace("_nozonal","") if self.x_quantity!=None else self.x_quantity
        y_quantity2 = self.y_quantity.replace("/ky**2","").replace("_real","").replace("_imag","").replace("_zonal","").replace("_nozonal","") if self.y_quantity!=None else self.y_quantity
        z_quantity2 = self.z_quantity.replace("/ky**2","").replace("_real","").replace("_imag","").replace("_zonal","").replace("_nozonal","") if self.z_quantity!=None else self.z_quantity 
        
        # Define what is on the (x,y) or (x,y,z) axis for getting the data
        if self.z_quantity==None: 
            self.xdim = x_quantity2
            self.ydim = y_quantity2
            self.zdim = None
        if self.z_quantity!=None: 
            self.xdim = x_quantity2
            self.ydim = y_quantity2
            self.zdim = z_quantity2
        
        # For e.g. z we can plot {pol; tor; zeta} 
        xdim = "z" if self.xdim in ["pol", "tor", "zeta"] else self.xdim
        
        # Get the quantity like it is saved as an object
        self.quantity = self.ydim+"_vs_"+xdim
        self.quantity = "g_vs_ts" if self.quantity=="g_vs_t" else self.quantity
        self.quantity = "g_vs_smu" if self.quantity=="g_vs_mu" else self.quantity
        self.quantity = "g_vs_svpa" if self.quantity=="g_vs_vpa" else self.quantity 
        self.quantity = "qflux_vs_ts" if self.quantity=="qflux_vs_t" else self.quantity
        self.quantity = "pflux_vs_ts" if self.quantity=="pflux_vs_t" else self.quantity
        self.quantity = "vflux_vs_ts" if self.quantity=="vflux_vs_t" else self.quantity 
        return

    #------------------------------------
    def get_xyzLabels(self, research):  
        
        # Get the units
        units = self.units if self.rescaleToOne!=True else "rescaled"

        # Read the labels from the labels dictionary  
        if hasattr(self, "xname"):
            if self.x_label==None: self.x_label = standardLabels[units][self.xname]
            if self.y_label==None: self.y_label = standardLabels[units][self.yname]
            if self.z_label==None: self.z_label = standardLabels[units][self.zname]
        
        # If "_{s}" is present, replace it with the correct species 
        if self.species: 
            
            # Determine the species based on the input parameters
            if len(self.species)==1:
                self.s = recognize_species(research, self.species[0]) 
            
            # If we plot multiple species, then keep the label general
            else:  self.s = "_s"  
             
            # Replace "_{s}" ny the species
            if self.x_label!=None:
                self.x_label = self.x_label.replace("_{s}", "_"+self.s).replace(",{s}", self.s)
                self.y_label = self.y_label.replace("_{s}", "_"+self.s).replace(",{s}", self.s)
                self.z_label = self.z_label.replace("_{s}", "_"+self.s).replace(",{s}", self.s) 
 
        # Change the label when we plot a logarithmic quantity
        if self.takelogx: self.x_label =  "log("+self.x_label+")"  
        if self.takelogy: self.y_label =  "log("+self.y_label+")"   
        return

    #------------------------------------
    def get_title(self, verbose=False):  

        # We don't always use the full plot package
        if not hasattr(self, "ydim"): return   

        # Get a standard title 
        if self.title==None: 
            
            # Get the title 
            if self.ydim in standardTitles.keys(): 
                if self.xdim in standardTitles[self.ydim].keys(): 
                    self.title = standardTitles[self.ydim][self.xdim]
                else: 
                    self.title = ""; 
                    if verbose: print("NO TITLE SET FOR: "+str(self.ydim)+", "+self.xdim) 
            else: 
                self.title = ""; 
                if verbose: print("NO TITLE SET FOR:"+str(self.ydim)+", "+self.xdim) 
            
            # Put a \n in the title if title_length is small
            if self.title_length=='small' and '\n' not in self.title: 
                self.title = split_titleInTwoSentences(self.title)
                
            # Specify that we extracted the y-quantity of the most unstable mode or at a specific ky_value.  
            if self.specific_k!=None and self.title!="" and self.title!=None and '\n' not in self.title:
                y_quantity = self.title.split(" versus ")[0]
                x_quantity = self.title.split(" versus ")[-1]
                if y_quantity!="gamma_integrated": addition = " between ["+str(self.ky_range[0])+", "+str(self.ky_range[1])+"]"
                if self.specific_k==None: addition = " of the most unstable mode"
                if self.specific_k!=None: addition = " at $k_y =$ "+(str(self.specific_k))
                if len(x_quantity)>=33: self.title = y_quantity+addition+" versus \nthe "+x_quantity
                if len(x_quantity)<=33: self.title = y_quantity+addition+"\nversus the "+x_quantity
        
        # Replace the species
        if self.species and self.title!=None: self.title.replace("_{s}", self.s)
        
        # Font size
        if not self.fontsize_title:
            if self.title_size==None:      self.fontsize_title = 18
            if self.title_size=='small':   self.fontsize_title = 15 
        return
        
################################################################################
#                             INITIALIZATIONS                                  #
################################################################################
 
    def update_simulations(self, experiment_id="All experiments", simulation_id="All simulations", species=[0]):
        self.experiment_id = experiment_id
        self.simulation_id = simulation_id
        self.species = species
        return     
    
    #------------------------------------
    def update_modes(self, modes_id="all", kx_range=[0,0], ky_range=[0,100], specific_k=None):
        self.specific_k = specific_k
        self.modes_id = modes_id
        self.kx_range = kx_range
        self.ky_range = ky_range 
        return    
    
    #------------------------------------
    def update_parameter(self, knob="vmec_parameters", key="rho", knob1=None, key1=None, knob2=None, key2=None):
        self.knob = knob
        self.key = key 
        self.knob1 = knob1
        self.key1 = key1 
        self.knob2 = knob2
        self.key2 = key2 
        return

    #------------------------------------
    def update_xyData(self, x_quantity=None, y_quantity=None, x_range=None, y_range=None, \
                       units="normalized", y_quantities=[None], error=False, dimensions=2, logx=False, logy=False):
        self.y_quantities = y_quantities
        self.dimensions = dimensions
        self.x_quantity = x_quantity 
        self.y_quantity = y_quantity 
        self.x_range = x_range
        self.y_range = y_range 
        self.takelogx = logx
        self.takelogy = logy
        self.error = error
        self.units = units  
        return   
    
    #------------------------------------
    def update_xyzData(self, x_quantity=None, y_quantity=None, z_quantity=None, x_range=None, y_range=None, z_range=None, \
                       units="normalized", y_quantities=[None], error=False, dimensions=2, logx=False, logy=False):
        self.y_quantities = y_quantities
        self.dimensions = dimensions
        self.x_quantity = x_quantity 
        self.y_quantity = y_quantity 
        self.z_quantity = z_quantity 
        self.x_range = x_range
        self.y_range = y_range
        self.z_range = z_range
        self.takelogx = logx
        self.takelogy = logy
        self.error = error
        self.units = units 
        return   
    
    #------------------------------------
    def update_labels(self, title=None, x_label=None, y_label=None, z_label=None, 
                      title_size=None, fontsize_title=None, title_length="normal"):
        self.title = title
        self.x_label = x_label
        self.y_label = y_label
        self.z_label = z_label
        self.title_size = title_size
        self.title_length = title_length
        self.fontsize_title = fontsize_title
        return 
    
    #------------------------------------
    def update_figure(self, Progress=None, ax=None, show_figure=True, fig_size=(18, 9), logx=False, logy=False):
        self.show_figure = show_figure
        self.Progress = Progress
        self.fig_size = fig_size
        self.logx = logx 
        self.logy = logy 
        self.ax = ax 
        return 
    
    #------------------------------------
    def update_legend(self, fontsize=18, handlelength=1, ncol=1, loc=None, removelegend=False,\
                      shadow=True, splitInKineticAdiabatic=False, adiabaticIons=False, \
                      adiabaticElectrons=True, splitInKineticAdiabaticTEM=False):
        self.splitInKineticAdiabaticTEM = splitInKineticAdiabaticTEM
        self.splitInKineticAdiabatic = splitInKineticAdiabatic
        self.adiabaticElectrons = adiabaticElectrons
        self.adiabaticIons = adiabaticIons
        self.removelegend = removelegend
        self.handlelength = handlelength
        self.fontsize = fontsize
        self.shadow = shadow
        self.ncol = ncol
        self.loc = loc 
        return 
    
    #------------------------------------
    def update_axis(self, tick_style='sci'):
        self.tick_style = tick_style
        return
    
    #------------------------------------
    def update_toggles(self, rescaleToOne=False, normalize=False, show_error=False,  step=20,\
                       maxima=False, biggestmaxima=False, interpolate=False, add_dottedLine=False, \
                       show_timeFrame=False, show_errorStd=False, show_errorMinMax=False, cascade=False):
        self.add_dottedLine = add_dottedLine
        self.show_timeFrame = show_timeFrame
        self.biggestmaxima = biggestmaxima
        self.rescaleToOne = rescaleToOne
        self.interpolate = interpolate
        self.normalize = normalize
        self.cascade = cascade
        self.maxima = maxima
        self.step = step
        
        self.show_errorMinMax = show_errorMinMax
        self.show_errorStd = show_errorStd
        self.show_error = show_error
        self.show_errorMinMax = False if (self.show_errorStd==True) else self.show_errorMinMax
        return
    
#--------------------------------
def split_titleInTwoSentences(title):
    titleSplitInWords = title.split(" ")
    differenceInLength = [None]*len(titleSplitInWords)
    for i in range(len(titleSplitInWords)):
        chunck1 = ' '.join(titleSplitInWords[0:i]) 
        chunck2 = ' '.join(titleSplitInWords[i:]) 
        differenceInLength[i] = abs(len(chunck2)-len(chunck1))
    optimalSplitting = np.argmin(np.array(differenceInLength))
    title = ' '.join(titleSplitInWords[0:optimalSplitting])+"\n "+' '.join(titleSplitInWords[optimalSplitting:]) 
    return title
    
    