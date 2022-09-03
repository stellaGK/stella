#================================================================
######################## UPDATE GUI #############################
#================================================================

#=================================
# Update the frames and data
#================================= 
def update_GUI(root):
    ''' Function to update all GUI elements when the selection of simulations changes

    Notes
    -----
    Calling <class_progress> will call <root.update_idletasks()> before and after each call. 
    Therefore showing the progress bar keeps the GUI up to date automatically.
    '''

    # While plotting turn the cursor into a waiting cursor
    root.config(cursor="watch")
    root.update_idletasks() 
        
    # Divide progress bar in equal bits
    numberOfStatusses = 5
    lengthDivision = int(100/numberOfStatusses)

    # Start the progress bar
    root.TabSelectedFiles.class_progress.start("Updating GUI");  status=1 

    # First time files are selected: change text on "Select folders" button
    if len(root.input_files)==0:
        root.TabSelectedFiles.class_simulations.button_selectFiles['text']   = "Select simulations" 
        root.TabSelectedFiles.class_simulations.button_selectFolders['text'] = "Select folders" 
    if len(root.input_files)>0:
        root.TabSelectedFiles.class_simulations.button_selectFiles['text']   = "Add more simulations" 
        root.TabSelectedFiles.class_simulations.button_selectFolders['text'] = "Add more folders" 
    
    # Update the scrollable canvas for the research
    root.TabSelectedFiles.class_progress.move(lengthDivision*status,"Updating research list");             status=2
    root.TabSelectedFiles.class_research.update_treeView() 
        
    # Update the scrollable canvas for the simulations
    root.TabSelectedFiles.class_progress.move(lengthDivision*status,"Updating simulation list");             status=2
    root.TabSelectedFiles.class_simulations.update_treeView()     

    # Update the input and resolution parameters shown in the frames
    root.TabSelectedFiles.class_progress.move(lengthDivision*status,"Updating input parameters");            status=4
    root.TabSelectedFiles.class_inputParameters.update_frame()
    
    # Finish Progress bar
    root.TabSelectedFiles.class_progress.finish()

    # Update the cursor
    root.config(cursor="")


