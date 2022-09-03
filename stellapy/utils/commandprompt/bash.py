
import inspect
import os, sys, getopt, pathlib  


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
        self.arguments = get_default_args(function)  
        self.arguments['folder'] = pathlib.Path(os.getcwd()) 
        self.add_option('folder', 'pathlib', '', 'Choose the folder where the command needs to be executed.')
        pass
    
    #---------------------   
    def add_option(self, name, datatype, shortname, explanation): 
        self.options.append(Option(name, datatype, shortname, explanation)) 
        get_argument_short = ":" if datatype not in ["True", "False", "help"] else ""
        get_argument_long = "=" if datatype not in ["True", "False", "help"] else "" 
        self.short_options += shortname+get_argument_short if shortname!="" else ""
        self.long_options += [name+get_argument_long]  
        return
    
    #---------------------   
    def add_toggle(self, name, longname, explanation): 
        self.toggles.append(Toggle(name, longname, explanation)) 
        self.long_options += [longname]
        return

    #---------------------       
    def get_arguments(self):  
        
        # Add the help option last
        self.add_option('help', 'help', 'h', 'Print information about the command.')  
        
        # Get bash like arguments from the command prompt 
        try: opts_args = getopt.getopt(sys.argv[1:],self.short_options,self.long_options)[0] 
        except getopt.GetoptError: self.print_help(error=True); sys.exit(2) 
    
        # Asign the options and arguments to variables  
        for opt, arg in opts_args:  
            
            # Initiate the arg_type
            arg_type = "undefined"
            
            # Find out whether <opt> is an option or a toggle 
            if arg_type=="undefined":
                for option in self.options:
                    if "--"+option.name==opt or "-"+option.shortname==opt:
                        arg_type = "option"
                        break
            if arg_type=="undefined":
                for toggle in self.toggles:
                    if "--"+toggle.longname==opt:
                        arg_type = "toggle"
                        break
            if arg_type=="undefined":
                self.print_help()
                sys.exit()

            # Get the arguments for an option
            if arg_type=="option":
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
                        a = arg.split("[")[-1].split(",")[0]
                        b = arg.split("]")[0].split(",")[-1]
                        self.arguments[option.name] = [a,b]
                    elif "[" not in arg:
                        a = arg.split(" ")[0]
                        b = arg.split(" ")[-1]
                        self.arguments[option.name] = [a,b]
                elif option.datatype=='help': 
                    self.print_help()
                    sys.exit()
                    
            # Get the arguments for a toggle
            if arg_type=="toggle": 
                self.arguments[toggle.name] = toggle.longname
     
        # Show the parameters
        return self.arguments 
         
    #---------------------    
    def print_help(self, error=False, width=100):
        ''' When a bash command is called, print the function name, its descriptions and its arguments.''' 
    
        # Get the function path and name 
        path = "stellapy" + self.function_path.split("stellapy/")[-1].replace("//","/").replace("/",".").replace(".py","")
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
            if value==True:  value="True"
            if value==False: value="False"
            if name!=keys[-1]: print('       ', name+" = "+ str(value))
            if name==keys[-1]: print('       ', name+" = "+ str(value) + ")")
        print()
        
        # List the bash options
        if len(self.options)!=0:
            print("\033[1m BASH OPTIONS \033[0m")
            print("\033[1m ------------ \033[0m", "\n") 
            for option in self.options:
                long = "--"+option.name if option.name!="" else ""
                short = "-"+option.shortname if option.shortname!="" else ""
                if "\n" in option.explanation:
                    parts = option.explanation.split('\n')  
                    print(" ","{0:>5}".format(short), " ", "{0:<15}".format(long), "{0:<40}".format(parts[0])) 
                    for i in range(len(parts)-1):
                        print("  ","{0:>5}".format(""), " ", "{0:<15}".format(""), "{0:<40}".format(parts[i+1]))
                else: 
                    if len(long)>15: long = long[:15]
                    print(" ","{0:>5}".format(short), " ", "{0:<15}".format(long), "{0:<40}".format(option.explanation))  
            print()
            
        if len(self.toggles)!=0:
            print("\033[1m BASH TOGGLES \033[0m")
            print("\033[1m ------------ \033[0m", "\n") 
            for toggle in self.toggles: 
                if "\n" in toggle.explanation:
                    parts = toggle.explanation.split('\n')  
                    print("  ", "{0:6}".format(""), "{0:<15}".format("--"+toggle.longname), "{0:<40}".format(parts[0])) 
                    for i in range(len(parts)-1):
                        print("  ", "{0:6}".format(""), "{0:<15}".format(""), "{0:<40}".format(parts[i+1]))
                else: 
                    print("  ", "{0:6}".format(""), "{0:<15}".format("--"+toggle.longname), "{0:<40}".format(toggle.explanation))  
            print()
            
        # Mention that we had an error
        if error:
            print("\033[1m ERROR \033[0m")
            print("\033[1m ----- \033[0m", "\n") 
            print("    One of the bash options was not defined or an argument was missing or mistakinly given. ")
            print("         Possible short arguments: ", self.short_options)
            print("         Possible long arguments:  ", self.long_options)
            print() 
        return 
        
#-------------
class Option():
    def __init__(self, name, datatype, shortname, explanation):
        self.explanation = explanation
        self.shortname = shortname
        self.datatype = datatype
        self.name = name
        return
    
#-------------
class Toggle():
    def __init__(self, name, longname, explanation):
        self.explanation = explanation
        self.longname = longname 
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
    
    
    