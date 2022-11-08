
import copy
import numpy as np
from pathlib import PosixPath
 
def sort_listByLabel(experiments): 
    line_label = [ e.line_label for e in experiments]   
    line_label_sorted = sorted(copy.deepcopy(line_label)) 
    line_label_sorted = line_label_sorted[::-1]
    sorted_indexes = [line_label.index(i) for i in line_label_sorted]
    experiments = [experiments[i] for i in sorted_indexes]    
    experiments.insert(0, experiments.pop()) 
    return experiments

def sort_listByNumbers(input_list, source=None):

    if isinstance(input_list, list): 
        if len(input_list)==0:
            return input_list
    
    if isinstance(input_list, list): 
        if isinstance(input_list[0], PosixPath): 
            if "_ky" in str(input_list[0]):
                input_list = [str(i) for i in input_list]
                numbers = [ i.replace(".in","").split("_ky")[-1] for i in input_list]  
                numbers = [ "".join([s for s in value if s.isdigit() or s=="." or s=="-"]) for value in numbers ]   
                numbers = [ number[:-1] if number.endswith(".") else number for number in numbers ]  
                numbers = [float(n) for n in numbers] 
                sorted_indexes = list(np.array(numbers).argsort())  
                input_list = [input_list[i] for i in sorted_indexes] 
                input_list = [PosixPath(i) for i in input_list]
                return input_list
            input_list = [str(i) for i in input_list]
            numbers = [ "".join([s for s in value if s.isdigit() or s=="."]) for value in input_list ] 
            numbers = [ number[:-1] if number.endswith(".") else number for number in numbers ] 
            numbers = [ float(number.split(".")[-1]) if "." in number else float(number) for number in numbers ] 
            sorted_indexes = list(np.array(numbers).argsort()) 
            input_list = [input_list[i] for i in sorted_indexes] 
            input_list = [PosixPath(i) for i in input_list]
            return input_list
    
    if source==None: 
        if "k_y\\rho_i =" in input_list[0] and ":" in input_list[0]:
            try:  
                input_list_temp = [f.split(":")[-1] for f in input_list]
                numbers = [ float("".join([s for s in value if s.isdigit() or s=="." or s=="-"])) for value in input_list_temp ]
                sorted_indexes = list(np.array(numbers).argsort()) 
                input_list = [input_list[i] for i in sorted_indexes]  
            except: pass
        else: 
            numbers = [ "".join([s for s in value if s.isdigit() or s=="." or s=="-"]) for value in input_list ]
            numbers = [ number[:-1] if number.endswith(".") else number for number in numbers ]
            numbers = [ number if number!='' else "0" for number in numbers] 
            numbers = [ number.split(".")[-2]+"."+number.split(".")[-1] if "." in number else number for number in numbers ]
            numbers = [ float("-"+number.split("-")[-1]) if "-" in number else float(number) for number in numbers ]
            sorted_indexes = list(np.array(numbers).argsort()) 
            input_list = [input_list[i] for i in sorted_indexes]   
    
    if source=="line_label":
        try:
            line_label = [ e.line_label for e in input_list]  
            numbers = [ float("".join([s for s in value if s.isdigit() or s=="." or s=="-"])) for value in line_label ]
            sorted_indexes = list(np.array(numbers).argsort())
            input_list = [input_list[i] for i in sorted_indexes]   
        except: pass
        
    if source=="id":
        try: 
            line_label = [ e.id for e in input_list]  
            numbers = [ "0"+"".join([s for s in value if s.isdigit() or s=="." or s=="-"])+".0.0" for value in line_label ]
            numbers = [ "-"+number.split("-")[-1] if "-" in number else number for number in numbers]
            numbers = [ float(number.split(".")[0]+"."+number.split(".")[1]) for number in numbers ]
            sorted_indexes = list(np.array(numbers).argsort())
            numbers = [numbers[i] for i in sorted_indexes]
            input_list = [input_list[i] for i in sorted_indexes] 
            line_label = [ e.id for e in input_list]   
        except: pass
        
    if source=="variedValues":
        try: 
            numbers = [ float("".join([s for s in value if s.isdigit() or s=="." or s=="-"])) for value in input_list.variedValues ]
            sorted_indexes = list(np.array(numbers).argsort())
            input_list.variedValues = [input_list.variedValues[i] for i in sorted_indexes]
            input_list.simulations  = [input_list.simulations[i] for i in sorted_indexes]   
        except: pass
    
    return input_list





