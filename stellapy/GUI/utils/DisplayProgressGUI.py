
################################################################################
#                 DISPLAY READING/PLOTTING PROGRESS IN THE GUI                 #
################################################################################

class DisplayProgressGUI:
    
    def __init__(self, research, message, count="simulations", plot=None):
        
        # Only do this for the GUI
        if not hasattr(research, 'plot'):
            self.progressBar = None
            return 
        
        # Unpack the needed data
        quantities = research.plot.y_quantities
        Progress = research.plot.Progress
        species = research.plot.species
        length = 0

        # Only do this for the GUI
        if not Progress: 
            self.progressBar = None
            return 
        
        # Count the experiments
        if count=="experiments":
            length = len(research.plotted_experiments)

        # Count the simulations 
        if count=="simulations":
            for experiment in research.plotted_experiments:   
                length += len(experiment.plotted_simulations)*len(quantities)*len(species)  

        # Count the modes 
        if count=="modes":
            for experiment in research.plotted_experiments:    
                for simulation in experiment.plotted_simulations: 
                    length += len(simulation.plotted_modes)*len(quantities)*len(species)  

        # Count the kx modes 
        if count=="kx":
            for experiment in research.plotted_experiments:    
                for simulation in experiment.plotted_simulations:  
                    vec_kx = [ kx for kx in simulation.vec.kx if (kx >= plot.kx_range[0] and kx <= plot.kx_range[1])]
                    length += len(vec_kx)*len(quantities)*len(species)  
                    
        # Count the ky modes 
        if count=="ky":
            for experiment in research.plotted_experiments:    
                for simulation in experiment.plotted_simulations:  
                    vec_ky = [ ky for ky in simulation.vec.ky if (ky >= plot.ky_range[0] and ky <= plot.ky_range[1])]
                    length += len(vec_ky)*len(quantities)*len(species)  
            
        # Save the length and which plot we have
        self.iplot = 1
        self.length = length
        self.message = message
        self.progressBar = Progress
        
        # Adjust the message on the GUI
        if "Plot" in message:
            self.message1 = "Start plotting the "+self.message+"."
            self.message2 = self.message+": plotting data ("
        if "Read" in message:
            self.message1 = "Start reading "+self.message+"."
            self.message2 = self.message+": reading data ("
        
        # Now initiate the progress bar
        self.initiate_progressBar()
        
    
    def initiate_progressBar(self):
        if self.progressBar:
            self.progressBar.move(0,self.message1)  
    
    def move_progressBar(self):
        if self.progressBar:
            text = self.message2+str(self.iplot)+"/"+str(self.length)+")"
            self.progressBar.move((self.iplot+1)/(self.length+1)*100,text); self.iplot +=1     

    def finilize_progressBar(self):
        if self.progressBar:
            self.progressBar.move(100,"Finish plotting the "+self.message+".")  
        