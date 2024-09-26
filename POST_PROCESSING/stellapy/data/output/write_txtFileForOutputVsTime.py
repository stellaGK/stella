
import os 
from stellapy.utils.files.get_filesInFolder import get_filesInFolder 

#===============================================================================
#                      REPLACE *.outES BY *.DT1.outES
#===============================================================================
# Make sure to run this bash command on marconi after every simulation.
#     > reduce_sizeoutes -t 1
#===============================================================================

def write_txtFileForOutputVsTime(folder, dt=1, verbose=False):
    
    # Time step
    dt = int(dt) if int(dt)==dt else dt    
    
    # Suffix of the new file
    suffix = ".dt"+str(dt)+".out" 
    
    # Get the input files
    out_files = get_filesInFolder(folder, end=".out")
    out_files = [ f for f in out_files if f.name[0:7].isdigit()] #or f.name[6:13].isdigit()  
    out_files = [ f for f in out_files if 'dt' not in f.name]  
    if out_files==[]: return 
    
    # Go through the out files
    for out_file in out_files:   
        
        # Processing status
        status = "    ("+str(out_files.index(out_file)+1)+"/"+str(len(out_files))+")  " if len(out_files)>5 else "   "
                 
        # Check whether the reduced file is older than the netcdf file
        already_written = False
        new_file = out_file.with_suffix(suffix)
        if os.path.isfile(new_file):
            time_newFile = new_file.stat().st_mtime
            time_outFile = out_file.stat().st_mtime
            if time_newFile>time_outFile:
                already_written = True

        # Only write it if we have to  
        if already_written==False or not os.path.isfile(new_file): 

            try:

                # Read the out file  
                file = open(out_file, 'r') 
                text = file.read()
                file.close()
            
                # Keep the header and bottom 
                text = text.split("------------------------------------------------------------")  
                header, text = text[0], "------------------------------------------------------------"+text[1] 
                header1 = header.split("Beginning root solves to determine theta_vmec.")[0]
                header2 = header.split("Leaving vmec_to_stella_geometry_interface.")[-1] 
                header = header1+"Beginning root solves to determine theta_vmec.\n"+" Leaving vmec_to_stella_geometry_interface."+header2
                if '############################################################' in text:
                    text = text.split('############################################################')
                    text, bottom = text[0], '############################################################'+text[1]+'############################################################'+text[2]
                else: bottom = ""
                    
                # Make the header prettier
                if "Add author names.," in header:
                    header = header.split("Add author names.,")
                    header = header[0] + "\n".join(header[1].split("\n")[5:])
                if "==> Changing code_dt to cfl_dt*cfl_cushion" in header:
                    header = header.split('==> Changing code_dt to cfl_dt*cfl_cushion')
                    header1 = header[0]+'==> Changing code_dt to cfl_dt*cfl_cushion'
                    header2 = header[1].split("\n")
                    header2[0] = header2[0]+"\n"
                    header = header1 + "\n".join(header2)
                
                # Only keep the time at each dt step 
                text = text.split("\n")
                current_time = 0
                new_text = ""
                count = 0
                for i,line in enumerate(text): 
                    if i==0: new_text += line+"\n"; continue
                    if i==1: new_text += line+"\n"; continue
                    if i==len(text)-1: new_text += line+"\n"; break 
                    if count==8:
                        count = 0
                        new_text += line+"\n" 
                        continue
                    elif count>=1 and count<=7:
                        count += 1
                        new_text += line+"\n" 
                        continue
                    elif "CHANGING TIME STEP:" in line:
                        count += 1
                        new_text += line+"\n" 
                        continue
                    data = [l for l in line.split(" ") if l!=""]
                    if data==[]: new_text += "\n"; continue 
                    try:
                        time = float(data[1])
                        if time >= current_time+dt:
                            new_text += line+"\n";
                            current_time = current_time+dt 
                    except: pass
                        
                # Full file 
                text = header + new_text + bottom
                
                # Save the nex text file
                with open(new_file, 'w') as file:
                    file.write(text) 
                    print(status+"   ---> The out(t) file is saved as", new_file.name) 
            
            except:
                print(status+"Something went wrong writing the out(t) file:", new_file.parent.name+"/"+new_file.name)      
            
        else: 
            if verbose: print(status+"The out(t) file already exists:", new_file.parent.name+"/"+new_file.name)
    
    return
