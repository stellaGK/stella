#===============================================================================
# Package to process the output of the gyrokinetic code stella
#===============================================================================
''' 

===============================
 Diagnostic package for stella
===============================

Analysis and diagnostic tool developed for the gyrokinetic code stella.

Developed by Hanne Thienpondt.

'''

#===============================================================================
# Make every function avaliable as package.func instead of package.mod.func
#===============================================================================

def enforce_one_function_one_file(path, globals_, locals_, verbose=False):
    
    # Import packages
    import os, glob
    import importlib 
    
    # Depending on the system (windows or linux) paths are divided by / or \\  
    divider = '\\' if (os.name == 'nt') else '/'    
    
    # Get the modules (.py files) and subpackages (folders) in the <package>
    modules_in_package = glob.glob(path[0]+'/[!_]*.py')
    subpackages_in_package = glob.glob(path[0]+'/[!_]*/') 
    
    # Based on <path> find the <package> that we are loading
    if len(modules_in_package)!=0:
        function = "/"+glob.glob(path[0]+'/[!_]*.py')[0].split(divider)[-1]
        package_path = "stellapy."+glob.glob(path[0]+'/[!_]*.py')[0].split(divider+"stellapy"+divider)[-1].split(function)[0].replace(divider,".")
        if verbose: print("Function:", function, package_path)
    if len(subpackages_in_package)!=0:
        subpackage = "/"+glob.glob(path[0]+'/[!_]*/')[0].split(divider)[-2]
        package_path = "stellapy."+glob.glob(path[0]+'/[!_]*/')[0].split(divider+"stellapy"+divider)[-1].split(subpackage)[0].replace(divider,".")
        if verbose: print("Subpackage:", subpackage, package_path)

    # Only get the module and package name 
    module_names_in_package = [file_name.split(divider)[-1].split('.')[0] for file_name in modules_in_package]
    subpackage_names_in_package = [folder_name.split(divider)[-2] for folder_name in subpackages_in_package]
    
    # Keep track of the doc strings of the modules (.py files)
    doc_strings = []

    # Save the doc string from the module 
    for module_path in modules_in_package:
        module_path = "stellapy." + module_path.split("/stellapy/")[-1].replace("/", ".").replace(".py","") 
        module = importlib.import_module(module_path) 
        doc_strings.append(module.__doc__)  

    # Import all modules from this package and overwrite the modules by the functions inside
    for module_name in module_names_in_package:
            exec('from '+ package_path +' import ' + module_name, globals_, locals_) 
            if verbose: print('from '+ package_path +' import ' + module_name) 
            exec('from ' + package_path + "." + module_name + ' import *', globals_, locals_) 
            if verbose: print('from ' + package_path + "." + module_name + ' import *')

    # If the module has a doc string, give it to the function
    for i, module_path in enumerate(modules_in_package):
        module_path = "stellapy." + module_path.split("/stellapy/")[-1].replace("/", ".").replace(".py","") 
        module = importlib.import_module(module_path)
        function_name = module_path.split(".")[-1]   
        if hasattr(module, function_name):
            function = getattr(module, function_name)
            if not isinstance(function, dict):
                if doc_strings[i]!=None:
                    function.__doc__ = doc_strings[i] 

    # Import all subpackages
    for package_name in subpackage_names_in_package:
        exec('from ' + package_path + ' import ' + package_name) 
        if verbose: print('from ' + package_path + ' import ' + package_name)

    # Clean up
    del glob
    return

#===============================================================================
# Make sure the matplotlib parameters are updated whenever using stellapy
#===============================================================================
 
import stellapy.plot.utils.style.load_styleFigures 




