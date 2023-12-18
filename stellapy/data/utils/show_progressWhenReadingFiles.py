
#===============================================================================
#                   SHOW THE PROGRESS WHEN READING FILES
#===============================================================================

def show_progressWhenReadingFiles(self, message): 
    if self.Progress!=None: 
        current_simulation = self.simulation.simulations.index(self.simulation)
        total_simulations = len(self.simulation.simulations)
        position = current_simulation/total_simulations*100
        message = message + "("+str(current_simulation)+"/"+str(total_simulations)+")"
        self.Progress.move(position, message) 
    return 
