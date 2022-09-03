
import os, glob
import importlib 

#=====================================================================
# Make function avaliable as package.func instead of package.mod.func
#=====================================================================

def enforce_one_function_one_file():
        
    divider = '\\' if (os.name == 'nt') else '/'
    files_in_package = glob.glob(__path__[0]+'/[!_]*.py')
    packages_in_pacakge = glob.glob(__path__[0]+'/[!_]*/')
    mod_list = [file_name.split(divider)[-1].split('.')[0] for file_name in files_in_package]
    sub_pack_list = [folder_name.split(divider)[-2] for folder_name in packages_in_pacakge]
    doc_strings = []

    # Import all modules from this package
    for mod in mod_list: 
        exec('from . import ' + mod) 

    # Save the doc string from the module 
    for file_name in files_in_package:
        file_name = "stellapy." + file_name.split("/stellapy/")[-1].replace("/", ".").replace(".py","")
        mod_name, func_name = file_name.rsplit('.',1)
        module = importlib.import_module(mod_name)
        function = getattr(module, func_name)
        doc_strings.append(function.__doc__)

    # Now overwrite the modules by the functions inside
    for mod in mod_list:  
        exec('from .' + mod + ' import *') 

    # If the module has a doc string, give it to the function
    for i, file_name in enumerate(files_in_package):
        file_name = "stellapy." + file_name.split("/stellapy/")[-1].replace("/", ".").replace(".py","")
        mod_name, func_name = file_name.rsplit('.',1)
        module = importlib.import_module(mod_name)
        function = getattr(module, func_name)
        if doc_strings[i]!=None:
            function.__doc__ = doc_strings[i]
        
    # Import all subpackages
    for pack in sub_pack_list:
        exec('from . import ' + pack)
        print("pack", pack)

    # Clean up
    del glob
    try:
        del mod
    except:
        pass
    try:
        del pack
    except:
        pass
