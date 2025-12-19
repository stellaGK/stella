
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
            return {'phi2' : phi2_data,\
                    'time' : time_data}
        return False

    # Read the *.out_file file:
    out_data = open(out_file, 'r')
    data = np.loadtxt(out_file,skiprows=2,dtype='float').reshape(-1, 5)
    time_data = data[:,1]
    phi2_data = data[:,2]
    return {'phi2' : phi2_data,\
            'time' : time_data}

