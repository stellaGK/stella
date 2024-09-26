
# Plot the electron and ion labels, and perhaps impurity labels
def display_specieLineLabels(ax, species):
    if len(species)>1:
        ax.plot(-1, -1, lw=3, linestyle='-',  color='black', label='Ions') 
        ax.plot(-1, -1, lw=3, linestyle='--',  color='black', label='Electrons') 
    if len(species)>2:
        ax.plot(-1, -1, lw=3, linestyle=':', color='black', label='Impurities') 
    if len(species)>3:
        ax.plot(-1, -1, lw=3, linestyle='-.',  color='black', label='Impurities 2') 

# Make sure we use the corresponding line style in the plots
def get_speciesLineStyle(species, specie):
    if   (len(species)>1) and (str(specie)=='0'): style = '-'  
    elif (len(species)>1) and (str(specie)=='1'): style = '--'  
    elif (len(species)>1) and (str(specie)=='2'): style = ':'
    elif (len(species)>1) and (str(specie)=='3'): style = '-.' 
    else: style = '-' 
    return style