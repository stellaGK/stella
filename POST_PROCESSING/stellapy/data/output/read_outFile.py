
import os   
import numpy as np
      
def read_outFile(input_file):
    ''' Read the *.out file.'''

    # Find the final fields file corresponding to the input file
    out_file = input_file.with_suffix('.out')
    if not os.path.isfile(out_file):
            
        # Read the *.phi_vs_t file
        out_file = input_file.with_suffix('.dt1.phi_vs_t')
        if os.path.isfile(out_file):  
            data = np.loadtxt(out_file,skiprows=1,dtype='float').reshape(-1, 6)
            time_data = data[:,0]
            phi2_data = data[:,1] 
            return {"phi2" : phi2_data,\
                    "time" : time_data} 
        return False

    # Read the *.out_file file: 
    out_data   = open(out_file, 'r')
    out_text   = out_data.read().replace(' ', '')
    
    # Read the time axis
    time_data  = out_text.split("time=")[1:]
    time_data  = [ float(text.split("|phi|^2=")[0]) for text in time_data ]
    
    # Read the potential squared
    phi2_data  = out_text.split("|phi|^2=")[1:]
    phi2_data  = [ text.split("|apar|^2")[0] for text in phi2_data ]
    phi2_data  = [ float(text.replace("+","e+").replace("-","e-").replace("Ee","e")) for text in phi2_data ]
    
    # Return the data
    return {"phi2" : phi2_data,\
            "time" : time_data} 

    