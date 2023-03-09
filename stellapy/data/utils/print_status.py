
import os 
import pathlib
import numpy as np

def print_status(file_path, file_name, input_file, input_files, nesting_depth):
    status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   " 
    message = "   ---> The "+file_name+" file is saved as " +  get_filePath(file_path, nesting_depth)
    print(status+message) 
    return
    
def print_statusSkipped(file_path, message, input_file, input_files, nesting_depth):
    status = "    ("+str(input_files.index(input_file)+1)+"/"+str(len(input_files))+")  " if len(input_files)>1 else "   "
    message = "   ===> " + message + "("+ get_filePath(file_path, nesting_depth)+")"
    print(status+message) 
    return
    
def get_nestingDepth(files):
    if files==[]: return -1
    common_parent_directory = pathlib.Path(os.path.commonprefix(files))  
    common_parent_directory = common_parent_directory if os.path.isdir(common_parent_directory) else common_parent_directory.parent
    nesting_depth = np.max([len(input_file.parts) for input_file in files]) - len(common_parent_directory.parts) 
    return nesting_depth

def get_filePath(file_path, nesting_depth):
    if nesting_depth==1:   file_path = file_path.name
    elif nesting_depth==2: file_path = file_path.parent.name+"/"+file_path.name
    elif nesting_depth==3: file_path = file_path.parent.parent.name+"/"+file_path.parent.name+"/"+file_path.name
    elif nesting_depth==4: file_path = file_path.parent.parent.parent.name+"/"+file_path.parent.parent.name+"/"+file_path.parent.name+"/"+file_path.name
    elif nesting_depth==5: file_path = file_path.parent.parent.parent.parent.name+"/"+file_path.parent.parent.parent.name+"/"+file_path.parent.parent.name+"/"+file_path.parent.name+"/"+file_path.name
    elif nesting_depth==6: file_path = file_path.parent.parent.parent.parent.parent.name+"/"+file_path.parent.parent.parent.name+"/"+file_path.parent.parent.name+"/"+file_path.parent.name+"/"+file_path.name
    else: file_path = file_path.name
    return file_path  
    