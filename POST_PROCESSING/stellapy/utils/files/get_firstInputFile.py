
import os    
import pathlib

def get_firstInputFile(folder):
    ''' Find an input file in <folder> '''   
    for file_name in os.listdir(folder):
        if not file_name.startswith('.'):
            if file_name.endswith('.in') : 
                return folder / file_name 
    for directory in os.listdir(folder):
        if os.path.isdir(directory): 
            for file_name in os.listdir(directory):
                if not file_name.startswith('.'):
                    if file_name.endswith('.in') : 
                        return pathlib.Path(directory) / file_name 
    for directory in os.listdir(folder):
        if os.path.isdir(directory): 
            for sub_directory in os.listdir(folder):
                if os.path.isdir(sub_directory): 
                    for file_name in os.listdir(sub_directory):
                        if not file_name.startswith('.'):
                            if file_name.endswith('.in') : 
                                return pathlib.Path(sub_directory) / file_name 
    





