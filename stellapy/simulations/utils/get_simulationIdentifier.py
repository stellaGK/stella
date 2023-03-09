
def get_simulationsIdentifier(input_file):
    ''' Returns parent_folder+input_file which is used to indentify simulations. '''
    
    # The input file is a PosixPath object: use attributes <name> and <parent>
    file_name = input_file.stem
    folder    = input_file.parent
    upfolder  = folder.parent
    upupfolder  = upfolder.parent
    
    # Get the string versions to create the simulation identifier
    input_file_str  = str(file_name)
    last_folder_str = str(folder.name) 
    
    # Make sure folder isnt $RUNS/run01, carefull of $RUNS// folders
    if "run" in last_folder_str:
        if last_folder_str.replace("run","").isdigit():
            folder = folder.parent
            upfolder = folder.parent 
            last_folder_str = str(folder.name) 
    
    # Add the parent folder to distinguish between input_files with the same name in a different folder
    simulation_id = str(upupfolder.name) +"/" + str(upfolder.name) + "/" + str(folder.name) + "/" + input_file_str
    
    # Make sure the string doesn't start with a slash
    if simulation_id.startswith('/'): simulation_id = simulation_id[1:]
    
    # Make sure there are no slashes in the middle, h5 files can't deal with slashes
    simulation_id = simulation_id.replace('/', '__')
    
    # Return the simulation identifier 
    return simulation_id
