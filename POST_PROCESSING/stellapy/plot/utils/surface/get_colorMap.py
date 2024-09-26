"""

#===============================================================================
#                             PERSONALIZED COLOR MAPS                          #
#===============================================================================

Get a standard color map (with nan=black), or one of the personalized maps:
    {"black-jet", "red-white-blue"}

There is a manually defined color map named "red-white-blue" which is blue for 
negative values, white for zero, and red for positive values. The "black-jet"
color map is black for zero, to black out the stable modes, followed by a jet
map for the other values.

Hanne Thienpondt
01/09/2022

"""

#!/usr/bin/python3
import copy
import matplotlib.pyplot as plt 
from matplotlib.colors import LinearSegmentedColormap

#===============================================================================
#                             PERSONALIZED COLOR MAPS                          #
#===============================================================================

def get_colorMap(color_map="jet"): 
        
    # Jet map where the zero is set to black for stable modes
    if color_map=="black-jet":
        cdict_jet = {'red': 
                            ((0., 0, 0),     # Black
                            (0.01, 0, 0),    # Black
                            (0.11, 0, 0),    # Blue
                            (0.66, 1, 1),    # Yellow
                            (0.89, 1, 1),    # Orange
                            (1, 0.5, 0.5)),  # Dark red
                    'green': 
                            ((0., 0, 0),     # Black
                            (0.01, 0, 0),    # Black
                            (0.11, 0, 0),    # Blue
                            (0.375, 1, 1),   # Greenish
                            (0.64, 1, 1),    # Greenish
                            (0.91, 0, 0),    # Orange
                            (1, 0, 0)),      # Dark red
                    'blue': 
                            ((0., 0, 0),     # Black
                            (0.01, 0, 1),    # Black
                            (0.11, 1, 1),    # Blue
                            (0.34, 1, 1),    # Light blue
                            (0.65, 0, 0),    # Yellow
                            (1, 0, 0))}      # Dark red
         
        cdict_jet = LinearSegmentedColormap('my_colormap1',cdict_jet,256)
        cdict_jet = copy.copy(plt.get_cmap("jet"))
        cdict_jet.set_bad(color='black')
        return cdict_jet
     
    # Seismic map with red for positive and blue for negative values, 
    # Where zero is put to black for stable modes
    if color_map=="red-white-blue":
        cdict_seismic = {'red': 
                            [[0.  , 0.  , 0.  ],
                            [0.25,  0.  , 0.  ],
                            [0.495, 1.  , 0.  ],     # White
                            [0.5,   0.  , 0.  ],     # Black
                            [0.505, 0.  , 1.  ],     # White
                            [0.75,  1.  , 1.  ],
                            [1.  ,  0.5 , 0.5 ]], 
                        'green': 
                            [[0.  , 0.  , 0.  ],
                            [0.25, 0.  , 0.  ],
                            [0.495, 1.  , 0.  ],     # White
                            [0.5,   0.  , 0.  ],     # Black
                            [0.505, 0.  , 1.  ],     # White
                            [0.75, 0.  , 0.  ],
                            [1.  , 0.  , 0.  ]], 
                        'blue': 
                            [[0. , 0.3 , 0.3 ],     # Darkbue
                            [0.25, 1.  , 1.  ],     # Blue
                            [0.495, 1.  , 0.  ],     # White
                            [0.5,   0.  , 0.  ],     # Black
                            [0.505, 0.  , 1.  ],     # White
                            [0.75, 0.  , 0.  ],     # Red
                            [1.  , 0.  , 0.  ]]}    # Darkred
         
        cdict_seismic = LinearSegmentedColormap('my_colormap2',cdict_seismic,256)
        cdict_seismic = copy.copy(plt.get_cmap("seismic")) 
        cdict_seismic.set_bad(color='black')
        return cdict_seismic 
    
    # Load any other map
    cmap = copy.copy(plt.get_cmap(color_map))
    cmap.set_bad(color='black')
    return cmap
