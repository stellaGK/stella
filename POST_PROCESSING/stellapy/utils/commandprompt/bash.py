
"""

#===============================================================================
#                           Bash to Python interface                           #
#===============================================================================

Interface to allow python scripts to be called from the command prompt in the 
same way that bash commands would be called. A list of commands is given through:
    >> stellapy
    
To see the available options of a command, ask for help.
    >> stellapy --help (-h)

The bash commands are defined in "stellapy/commands.sh", and an overview of commands
is given in "stellapy/utils/commandprompt/list_commands.py", which can be accessed
through the bash command ">> stellapy".

Bash options (manually input a value) and toggles (preset values) are defined through:
    bash.add_option(argument name, datatype, bash short option, bash long option, explanation) 
    bash.add_toggle(argument name, argument value, bash short option, bash long option, explanation)  
    
Examples: 
    bash.add_option('kymax', 'float', 'k', 'kymax', 'Maximum ky.')  
    bash.add_toggle('modes_id', 'stable', 's', 'stable', 'Plot the stable or unconverged modes.') 
        >> plot_gamma_vs_time --kymax 2 --stable
        >> plot_gamma_vs_time -k 2 -s 

Hanne Thienpondt 
20/01/2023

"""

import glob
import inspect
import numpy as np
import os, sys, getopt, pathlib  

#===============================================================================
#                           Bash to Python interface                           #
#===============================================================================

class Bash(): 
    
    def __init__(self, function, doc):
        
        # Save some basic information of the function
        self.function_path = function.__code__.co_filename 
        self.function_name = function.__name__ 
        self.function_doc = doc
        
        # Create bash like options
        self.short_options = ""
        self.long_options = [] 
        self.options = []
        self.toggles = []
        
        # Get the default arguments of the given function
        explanation = "Select the folders that contain simulations, default is the current working directory.\n"
        explanation += "Can be a string or list of strings e.g. '/home/path' or ['path1', 'path2']. \n"
        explanation += "Can be folder names in the current directory e.g. '[rho*, lowresolution]'."
        self.arguments = get_default_args(function)  
        self.arguments['folder'] = pathlib.Path(os.getcwd()) 
        self.add_option('folder', 'pathlib', '', '', explanation)
        pass
    
    #---------------------   
    def add_option(self, name, datatype, shortoption, longoption, explanation): 
        if longoption=="": longoption = name 
        self.options.append(Option(name, datatype, shortoption, longoption, explanation)) 
        get_argument_short = ":" if datatype not in ["True", "False", "help"] else ""
        get_argument_long = "=" if datatype not in ["True", "False", "help"] else "" 
        self.short_options += shortoption+get_argument_short if shortoption!="" else ""
        self.long_options += [longoption+get_argument_long]  
        return 
    
    #---------------------   
    def add_toggle(self, name, value, shortoption, longoption, explanation): 
        if longoption=="": longoption = value
        self.toggles.append(Toggle(name, value, shortoption, longoption, explanation)) 
        self.short_options += shortoption
        self.long_options += [longoption]
        return
    
    #---------------------   
    def add_toggleheader(self, header):  
        self.toggles.append(Toggle("HEADER", " ", " ", header, " "))  
        return
    
    #---------------------   
    def add_togglespace(self):  
        self.toggles.append(Toggle("SPACE", " ", " ", " ", " "))  
        return
    
    #---------------------   
    def add_optionheader(self, header):  
        self.options.append(Option("HEADER", " ", " ", header, " "))  
        return
    
    #---------------------   
    def add_optionspace(self):  
        self.options.append(Option("SPACE", " ", " ", " ", " "))  
        return

    #---------------------       
    def get_arguments(self):  
        
        # Avoid circular import
        from stellapy.utils.decorators.exit_program import exit_program 
        
        # Check the option lists
        short_options_temp = [i for i in self.short_options if i!=':'] 
        if (len(list(set(short_options_temp)))<len(short_options_temp)) or (len(list(set(self.long_options)))<len(self.long_options)):
            exit_reason = "One of the short of long options is used twice, please correct the script.\n"
            exit_reason += f"   Short options: {self.short_options}\n"
            exit_reason += f"   Long options: {self.long_options}\n"
            exit_program(exit_reason, Bash, sys._getframe().f_lineno)  
        
        # Add the help option last
        self.add_toggle('help', 'help', 'h', 'help', 'Print information about the command.')  
        
        # Get bash like arguments from the command prompt  
        try: opts_args = getopt.getopt(sys.argv[1:],self.short_options,self.long_options)[0] 
        except Exception as error_message: self.print_help(error=True, error_message=error_message); sys.exit()
    
        # Asign the options and arguments to variables  
        for opt, arg in opts_args:  
             
            # Initiate the arg_type
            arg_type = "undefined"
            
            # Find out whether <opt> is an option or a toggle 
            if arg_type=="undefined":
                for option in self.options:
                    if "--"+option.longoption==opt or "-"+option.shortoption==opt:
                        arg_type = "option" 
                        break
            if arg_type=="undefined":
                for toggle in self.toggles: 
                    if "--"+toggle.longoption==opt or "-"+toggle.shortoption==opt:
                        arg_type = "toggle"
                        break 
            if arg_type=="undefined" or opt=="--help":
                self.print_help()
                sys.exit() 
            if arg_type=="toggle" and opt=="-h":
                if toggle.longoption=="help":
                    self.print_help()
                    sys.exit()  
                    
            # Special treatment of the 'folder' option
            if arg_type=="option" and option.name=='folder': 
                 
                # Get the folder where the script is launched
                cwd = pathlib.Path(os.getcwd())
                found_folders = []
                 
                # <folder> can be e.g. ['nmu*', 'lowresolution'] or '/home/path'
                if '[' not in arg: given_folders = [arg]
                if '[' in arg: given_folders = [i.replace(' ','') for i in arg.split('[')[-1].split(']')[0].split(',')]
                for i, path in enumerate(given_folders): 
                    if '*' in path and '/' not in path: 
                        matching_files = [pathlib.Path(i) for i in glob.glob(str(cwd/path))]
                        found_folders += matching_files
                        if len(matching_files)==0:
                            exit_reason = f'The <folder> = {arg} could not be found.'
                            exit_program(exit_reason, Bash, sys._getframe().f_lineno)  
                    elif '*' in path and '/' in path: 
                        matching_files = [pathlib.Path(i) for i in glob.glob(path)]
                        found_folders += matching_files
                        if len(matching_files)==0:
                            exit_reason = f'The <folder> = {arg} could not be found.'
                            exit_program(exit_reason, Bash, sys._getframe().f_lineno)  
                    elif os.path.isdir('/'+path) or os.path.isfile('/'+path): 
                        found_folders.append(pathlib.Path(path))
                    elif os.path.isdir(cwd/path) or os.path.isfile(pathlib.Path(cwd/path)):
                        found_folders.append(cwd/path)
                    else: 
                        exit_reason = f'The <folder> = {arg} could not be found.'
                        exit_program(exit_reason, Bash, sys._getframe().f_lineno)  
                self.arguments[option.name] = found_folders

            # Get the arguments for an option
            elif arg_type=="option":
                if option.datatype=='str':
                    self.arguments[option.name] = str(arg)
                if option.datatype=='int':
                    self.arguments[option.name] = int(arg)
                elif option.datatype=='float':
                    self.arguments[option.name] = float(arg)
                elif option.datatype=='True':
                    self.arguments[option.name] = True
                elif option.datatype=='False':
                    self.arguments[option.name] = False
                elif option.datatype=='pathlib':
                    self.arguments[option.name] = pathlib.Path(arg) 
                elif option.datatype=='hours': 
                    if int(arg)<10:  self.arguments[option.name] = "0"+str(int(arg))+":00:00" 
                    if int(arg)>=10: self.arguments[option.name] = str(int(arg))+":00:00" 
                elif option.datatype=='minutes':
                    if int(arg)<10:  self.arguments[option.name] = "00:0"+str(int(arg))+":00" 
                    if int(arg)>=10: self.arguments[option.name] = "00:"+str(int(arg))+":00" 
                elif option.datatype=='species':
                    self.arguments[option.name] = [int(arg)]
                elif option.datatype=='range':
                    if "default" in arg:
                        self.arguments[option.name] = "default" 
                    elif "[" in arg:  
                        a = float(arg.split("[")[-1].split(",")[0])
                        b = float(arg.split("]")[0].split(",")[-1])
                        self.arguments[option.name] = [a,b]
                    elif "(" in arg:  
                        a = float(arg.split("(")[-1].split(",")[0])
                        b = float(arg.split(")")[0].split(",")[-1])
                        self.arguments[option.name] = [a,b]
                    elif "[" not in arg:
                        a = float(arg.split(" ")[0])
                        b = float(arg.split(" ")[-1])
                        self.arguments[option.name] = [a,b]
                elif option.datatype=='vector':
                    arg = arg.split("[")[-1].split("]")[0]
                    arg = [float(i) for i in arg.split(",")]
                    self.arguments[option.name] = arg
                elif option.datatype=='vector_of_tuples':
                    arg = arg.split("[")[-1].split("]")[0] 
                    arg = [i.split("(")[-1].split(")")[0] for i in arg.split(";")]
                    arg = [[float(i.split(",")[0]),float(i.split(",")[-1])] for i in arg]
                    self.arguments[option.name] = arg
                elif option.datatype=='vector_of_strings':
                    arg = arg.split("[")[-1].split("]")[0]
                    arg = [str(i).replace(' ','') for i in arg.split(",")]
                    self.arguments[option.name] = arg
                elif option.datatype=='help': 
                    self.print_help()
                    sys.exit()
                    
            # Get the arguments for a toggle
            elif arg_type=="toggle":  
                self.arguments[toggle.name] = toggle.value
     
        # Show the parameters  
        return self.arguments 
         
    #---------------------    
    def print_help(self, error=False, width=100, error_message=""):
        ''' When a bash command is called, print the function name, its descriptions and its arguments.''' 
    
        # Get the function path and name 
        path = "stellapy." + self.function_path.split("stellapy/")[-1].replace("//","/").replace("/",".").replace(".py","")
        name = self.function_name
        
        # Print the function name
        print("\n"+"".center(width,"="))
        print("\033[1m"+name.center(width," ")+"\033[0m")
        print("".center(width,"="), "\n") 
        
        # Add the function description 
        print("\033[1m DESCRIPTION \033[0m")
        print("\033[1m ----------- \033[0m", "\n")
        doc = self.function_doc; pops = []
        doc = [i if len(i)!=0 else " " for i in doc.split("\n") ]
        doc = ["    "+i for i in doc if i[0]!="#"]
        emptylines = "\n".join(doc).replace(" ","").split("\n"); number_of_lines = len(doc)
        emptylines = [i for i in range(len(emptylines)) if emptylines[i]==""] 
        for i in range(len(emptylines)):   
            if emptylines[i]==i: pops.append(i) 
            if emptylines[len(emptylines)-1-i]==(number_of_lines-1-i): pops.append(number_of_lines-1-i) 
        doc = "\n".join([doc[i] for i in range(number_of_lines) if i not in pops])
        print(doc) 
        print()
        
        # Location of the function and function arguments
        print("\033[1m FUNCTION \033[0m")
        print("\033[1m -------- \033[0m \n")
        print("    "+path+"(")
        keys = list(self.arguments.keys())
        if "folder" in keys: keys.remove("folder"); keys = ["folder"] + keys
        for name in keys: 
            value = self.arguments[name]
            if isinstance(value, str): value="'"+value+"'"
            if isinstance(value, list): value = "["+', '.join([str(i) for i in value])+"]"
            if value==None:  value="None"
            if isinstance(value, (bool)):
                if value==True:  value="True"
                if value==False: value="False"
            if name!=keys[-1]: print('       ', name+" = "+ str(value)+", ")
            if name==keys[-1]: print('       ', name+" = "+ str(value) + ")")
        print()
        
        # List the bash options and toggles
        long_option_length = 10
        for option in self.options: long_option_length = np.max([long_option_length, len(option.longoption)])
        for toggle in self.toggles: long_option_length = np.max([long_option_length, len(toggle.longoption)])
        if len(self.options)!=0:
            print("\033[1m BASH OPTIONS \033[0m")
            print("\033[1m ------------ \033[0m", "\n") 
            self.options = self.options[1:] + [self.options[0]]
            for option in self.options: self.print_option(option, long_option_length)
            print() 
        if len(self.toggles)!=0:
            print("\033[1m BASH TOGGLES \033[0m")
            print("\033[1m ------------ \033[0m", "\n") 
            for toggle in self.toggles: self.print_option(toggle, long_option_length)
            print()
            
        # Mention that we had an error
        if error:
            short_options_temp = self.short_options+"h"
            short_toggles = [short_options_temp[i] for i in range(len(short_options_temp)-1) if short_options_temp[i+1]!=":" and short_options_temp[i]!=":"]
            short_options = [short_options_temp[i] for i in range(len(short_options_temp)-1) if short_options_temp[i+1]==":"]
            long_toggles = [i for i in self.long_options if i[-1]!="="]
            long_options = [i[:-1] for i in self.long_options if i[-1]=="="]
            print("\033[1m ERROR \033[0m")
            print("\033[1m ----- \033[0m", "\n") 
            print("    One of the bash options was not defined or an argument was missing or mistakinly given. ")
            print("         Possible short toggles: ", short_toggles) 
            print("         Possible short options: ", short_options) 
            print("         Possible long toggles:  ", long_toggles)
            print("         Possible long options:  ", long_options)                    
            print() 
            if error_message!="": print("    \033[1mERROR\033[0m: "+str(error_message)+".\n")
        return 
    
    #-------------
    def print_option(self, option, long_option_length):
        
        # Sort options in sections with a certain header, and leave an empty space below
        if option.name=="HEADER": print("    \033[1m"+(option.longoption)+"\033[0m"); print(("    "+"-"*len(option.longoption))); return;
        if option.name=="HEADER": print("\033[1m"+(option.longoption).center(50)+"\033[0m"); print(("-"*len(option.longoption)).center(50)); return;
        if option.name=="SPACE": print(); return;
        
        # Turn the short and long option into a string
        long = "--"+option.longoption if option.longoption!="" else ""
        short = "-"+option.shortoption if option.shortoption!="" else ""
        long_format = "{0:<"+str(long_option_length+4)+"."+str(long_option_length+2)+"}"
        options = " "+"{0:>5}".format(short)+", "+long_format.format(long) if short!="" else " "*8+long_format.format(long)
        
        # If we have a long explanation, split it up into multiple lines
        if "\n" in option.explanation:
            parts = option.explanation.split('\n')  
            print(options, "{0:<40}".format(parts[0])) 
            for i in range(len(parts)-1):
                print(" "*(long_option_length+12), "{0:<40}".format(parts[i+1]))
                
        # If we have a short explanation, print it in one line
        else: 
            print(options, "{0:<40}".format(option.explanation))   
        return
        
#-------------
class Option():
    def __init__(self, name, datatype, shortoption, longoption, explanation):
        self.explanation = explanation
        self.shortoption = shortoption
        self.longoption = longoption
        self.datatype = datatype
        self.name = name
        return 
    
#-------------
class Toggle():
    def __init__(self, name, value, shortoption, longoption, explanation):
        self.explanation = explanation
        self.shortoption = shortoption 
        self.longoption = longoption 
        self.value = value 
        self.name = name
        return
        
#-------------------------
def get_default_args(func):
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }     
    
    
    