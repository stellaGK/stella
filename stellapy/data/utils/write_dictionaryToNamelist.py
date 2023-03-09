
import os
import f90nml
 
def write_dictionaryToNamelist(path, dictionary, indent="  ", sort_knobs=False):
    
    # Turn the input file into an empty namelist and preserve the comments
    inline_comments = write_emptyNamelist(path)
    
    # Add the dictionary to the empty namelist
    f90nml.patch(path, dictionary, str(path)+"_temp")
    
    # Remove the large indents
    with open(str(path)+"_temp", 'r') as file :
        filedata = file.read() 
    filedata = filedata.replace("    ", indent) 
    
    # Add the inline comments again  
    for knob_key, comment in inline_comments.items():
        knob = knob_key.split(": ")[0]
        key = knob_key.split(": ")[-1]
        before = filedata.split("&"+knob)[0]+"&"+knob
        after = "/"+"/".join(filedata.split("&"+knob)[-1].split("/")[1:])
        section = filedata.split("&"+knob)[-1].split("/")[0]  
        line = [line for line in section.split("\n") if (key+" = ") in line][0] 
        section = section.replace(line, line+comment)
        filedata = before+section+after 
        
    # Write the text file 
    with open(str(path), 'w') as file:
        file.write(filedata) 
    os.system("rm "+str(path)+"_temp")  
    
    # Sort the knobs in the input file
    if sort_knobs: sort_knobsInInputFile(path)
    return 

#--------------------------------------------
def write_emptyNamelist(path):

    # Read the input file
    with open(path, 'r') as file:
        filedata = file.read() 
    
    # Save the inline comments
    inline_comments = {}
    filedata_comments = [line for line in filedata.split("\n") if (("!" in line) and "=" in line)]
    for comment in filedata_comments: 
        knob = filedata.split(comment)[0].split("&")[-1].split("\n")[0] 
        key = comment.split("=")[0].replace(" ","")
        white_space = comment.split("!")[0][::-1]
        white_space = " "*(len(white_space) - len(white_space.lstrip()))
        inline_comments[knob+": "+key] = white_space+"!"+comment.split("!")[-1]
    
    # Remove all the data within the knobs, only keep knobs and comments
    filedata = filedata.split("&")
    for i in range(1,len(filedata)):
        filedata[i] = "&"+filedata[i].split("\n")[0]+"\n" + "/"+filedata[i].split("/")[-1]
    filedata = "".join(filedata)
    if filedata[-1]=="\n": filedata = filedata[:-1]
   
    # Write the "empty" namelist 
    with open(path, 'w') as file:
        file.write(filedata) 
    return inline_comments 

#-----------------------
def sort_knobsInInputFile(path):
    
    # Read the input file
    with open(path, 'r') as file :
        input_file_text = file.read()
        
    # Replace e.g. restart_file = 'restart/restart.nc' by restart_file = 'restart|restart.nc'
    lines = input_file_text.split('\n')
    lines = [l.replace("/","|") if (("'" in l or '"' in l) and "/" in l) else l for l in lines]
    input_file_text = "\n".join(lines)
    
    # Split the input file into its knobs
    input_file_sections = input_file_text.split('/')
    input_file_sections[0] = '\n'+input_file_sections[0] if input_file_sections[0][0]!="&" else input_file_sections[0]
    input_file_knobs = [section.split('&')[-1].split('\n')[0].replace(' ','') for section in input_file_sections]
    
    # Remove unused knobs
    remove_knobs = ['reinit_knobs']
    for i, knob in enumerate(input_file_knobs):
        if knob=='kt_grids_knobs':
            for line in input_file_sections[i].split('\n'):
                if 'grid_option' in line and 'range' in line:
                    remove_knobs.append('kt_grids_box_parameters')
                if 'grid_option' in line and 'box' in line:
                    remove_knobs.append('kt_grids_range_parameters')
        elif knob=='neoclassical_input':
            for line in input_file_sections[i].split('\n'):
                if 'include_neoclassical_terms' in line and 'false' in line:
                    remove_knobs.append('neoclassical_input')
                    remove_knobs.append('sfincs_input')
    input_file_sections = [input_file_sections[i] for i in range(len(input_file_sections)) if input_file_knobs[i] not in remove_knobs]
    input_file_knobs = [input_file_knobs[i] for i in range(len(input_file_knobs)) if input_file_knobs[i] not in remove_knobs]
        
    # Order the knobs
    ordered_knobs = ['physics_flags', 'species_knobs', 'species_parameters_1', 'species_parameters_2', 'species_parameters_3', \
    'parameters', 'vmec_parameters', 'kt_grids_box_parameters', 'kt_grids_range_parameters', 'zgrid_parameters', 'vpamu_grids_parameters', \
    'knobs', 'stella_diagnostics_knobs', 'init_g_knobs', 'dissipation', 'hyper', 'geo_knobs', 'kt_grids_knobs', \
    'dist_fn_knobs', 'time_advance_knobs', 'layouts_knobs', 'reinit_knobs', 'neoclassical_input', 'sfincs_input']
    input_file_sections_temp = []
    for ordered_knob in ordered_knobs:
        for i, knob in enumerate(input_file_knobs):
            if ordered_knob==knob: 
                input_file_sections_temp.append(input_file_sections[i])
                input_file_sections.pop(i)
                input_file_knobs.pop(i)
                break
    input_file_sections_temp += input_file_sections
    
    # Put the knobs together again
    input_file_sections_temp[0] = input_file_sections_temp[0][1:] if input_file_sections_temp[0][:1]=='\n' else input_file_sections_temp[0]
    input_file_text = '/'.join(input_file_sections_temp)
    input_file_text = input_file_text.replace('/&', '/\n&') 
        
    # Replace e.g. restart_file = 'restart|restart.nc' by restart_file = 'restart/restart.nc'
    input_file_text = input_file_text.replace('|','/')
    
    # Write the input file
    with open(str(path), 'w') as file:
        file.write(input_file_text) 
    return 



