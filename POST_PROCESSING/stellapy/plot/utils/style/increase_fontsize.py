
import matplotlib as matplotlib

def increase_fontsize(ax, fontsize, cbar=None, label=None, labelpad=20):
    
    # Increase the font size of the ticks
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize)
    
    # Increase the font size of the title and labels
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(fontsize) 
        
    # Increase the font size of the color bar
    if cbar!=None:
        cbar.ax.tick_params(labelsize=fontsize)  
        if label: cbar.set_label(label, labelpad=labelpad, y=0.52)
        text = cbar.ax.yaxis.label
        font = matplotlib.font_manager.FontProperties(size=fontsize)
        text.set_font_properties(font)
        label = text.get_text()  
        cbar.ax.yaxis.get_offset_text().set_fontsize(fontsize)
    return