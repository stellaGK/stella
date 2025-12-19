#!/usr/bin/python3
import os
import shutil
import pathlib
import difflib
import platform
import subprocess
import configparser

# Package to convert input files
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'convert_input_files/convert_inputFile.py')
with open(module_path, 'r') as file: exec(file.read())

################################################################################
#                 Routines to launch a local stella simulation                 #
################################################################################
# Note that the argument of any test function is the temporary path where the 
# test is performed, i.e., test_*(tmp_path) executes stella in <tmp_path>. 
################################################################################
 
#-------------------------------------------------------------------------------
def read_nproc():
    test_directory = get_automatic_tests_directory() 
    config = configparser.ConfigParser() 
    config.read(test_directory / 'config.ini')
    nproc = config['DEFAULT']['nproc']
    return nproc
    
#-------------------------------------------------------------------------------
def get_automatic_tests_directory():
    '''Get the directory of the test modules'''
    return pathlib.Path(str(pathlib.Path(__file__)).split('AUTOMATIC_TESTS')[0]) / 'AUTOMATIC_TESTS'
 
#-------------------------------------------------------------------------------
def get_stella_expected_run_directory():
    '''Get the directory of this test file.'''
    return pathlib.Path(__file__).parent

#-------------------------------------------------------------------------------
def get_stella_path(stella_version):
    '''Returns the absolute path to the stella executable.
    Can be controlled by setting the STELLA_EXE_PATH environment variable.'''
    # Note that os.environ.get() returns <default_path> if $STELLA_EXE_PATH does not exist
    if stella_version=='master':
       default_path = get_stella_expected_run_directory() / '../../../stella'
       stella_path = pathlib.Path(os.environ.get('STELLA_EXE_PATH', default_path))
    if stella_version=='0.5':
       stella_path = get_stella_expected_run_directory() / '../../stella_releases/v0.5/stella'
    if stella_version=='0.6':
       stella_path = get_stella_expected_run_directory() / '../../stella_releases/v0.6/stella'
    if stella_version=='0.7':
       stella_path = get_stella_expected_run_directory() / '../../stella_releases/v0.7/stella'
    return stella_path.absolute()

#-------------------------------------------------------------------------------
def run_stella(stella_path, input_file, nproc=None):
    '''Run stella with a given input file.''' 
    if not nproc: nproc = read_nproc()
    print(f'Execute: mpirun --oversubscribe -np {nproc} {stella_path} {input_file}')
    subprocess.run(['mpirun','--oversubscribe', '-np', f'{nproc}', stella_path, input_file], check=True)
    return

#-------------------------------------------------------------------------------
def copy_input_file(input_file: str, destination):
    '''Copy input_file to destination directory.'''
    shutil.copyfile(get_stella_expected_run_directory() / input_file, destination / input_file)
    return destination / input_file

#-------------------------------------------------------------------------------
def copy_common_input_files(input_file: str, destination):
    '''Copy input_file to destination directory.'''
    input_files = os.listdir(get_stella_expected_run_directory()/'../common_input_files')
    if not os.path.exists(destination/'common_input_files'):
        os.makedirs(destination/'common_input_files')
    for i in input_files:
      shutil.copyfile(get_stella_expected_run_directory()/f'../common_input_files/{i}', destination/f'common_input_files/{i}')
    return

#-------------------------------------------------------------------------------
def copy_vmec_file(vmec_file: str, destination):
    '''Copy input_file to destination directory.'''
    shutil.copyfile(get_stella_expected_run_directory() / vmec_file, destination / vmec_file)
    return
    
#-------------------------------------------------------------------------------
def run_local_stella_simulation(input_file, tmp_path, stella_version, vmec_file=None, nproc=None):
    ''' Run a local stella simulation in <tmp_path>. '''
    
    # Make sure the selected stella version is implemented
    if stella_version not in ['master', '0.5', '0.6', '0.7']: 
        print(f'ABORT: Wrong stella version: {stella_version}'); sys.exit()
        
    # Copy the input file from the automatic tests folder to a temp folder
    path_input_file = copy_input_file(input_file, tmp_path)
    
    # If we want to test older stella versions, convert the input file
    # Always turn of electromagnetic effects on old stella versions
    if stella_version!='master': 
        update_inputFile(path_input_file, downgrade=True)
        path_input_file_downgraded = str(path_input_file).replace('.in', '_downgraded.in')
        input_file = pathlib.Path(path_input_file).name
        shutil.copyfile(path_input_file_downgraded, path_input_file)
        
    # Copy the VMEC files to the temp folder
    if vmec_file: copy_vmec_file(vmec_file, tmp_path)
    
    # Switch to the temp folder, and run stella from within this folder
    os.chdir(tmp_path)
    run_stella(get_stella_path(stella_version), input_file, nproc=nproc)
    
    # Return the run data
    run_data = {'input_file' : input_file, 'tmp_path' : tmp_path, 'vmec_file' : vmec_file}
    return run_data
    
################################################################################
#                         Routines to compare txt files                        #
################################################################################
    
#-------------------------------------------------------------------------------  
def compare_local_txt_with_expected_txt(local_file, expected_file, name, skiplines=0, error=False):
     
    # Check whether the files match  
    if skiplines==0:
        if os.path.getsize(local_file) != os.path.getsize(expected_file): error = True
        elif open(local_file,'r').read() != open(expected_file,'r').read(): error = True 
    if skiplines>0:
        if open(local_file,'r').readlines()[skiplines:] != open(expected_file,'r').readlines()[skiplines:]: error = True 
    
    # If the files don't match, print the differences
    if error==True:
        print(f'ERROR: {name} files do not match.')
        print_differences_in_text_files(local_file, expected_file, name=name.upper(), skiplines=skiplines)
        assert False, f'{name} files do not match.' 
     
#-------------------------------------------------------------------------------
def print_differences_in_text_files(file1, file2, name='', skiplines=0, maxlines=20): 
    print(f'\nDIFFERENCE BETWEEN {name} FILES:\n') 
    with open(file1) as f1:
        local_file = f1.readlines()[skiplines:]
    with open(file2) as f2:
        expected_file = f2.readlines()[skiplines:]
    if len(local_file)>maxlines:
        print(f'    Warning: the compared files are very long, we limit them to {maxlines} lines')
        local_file = local_file[:maxlines]; expected_file = expected_file[:maxlines]
    for line in difflib.unified_diff(local_file, expected_file, fromfile=str(file1), tofile=str(file2), lineterm=''): 
        print('    ', line) 
        
################################################################################
#                       Routines to compare netCDF files                       #
################################################################################
    
old_to_new_keys = {
        'pflx' : 'pflux_vs_s',
        'qflx' : 'qflux_vs_s',
        'vflx' : 'vflux_vs_s',
        'pflx_kxky' : 'pflux_vs_kxkyzs',
        'qflx_kxky' : 'qflux_vs_kxkyzs',
        'vflx_kxky' : 'vflux_vs_kxkyzs',
        'pflx_vs_kxky' : 'pflux_vs_kxkys',
        'qflx_vs_kxky' : 'qflux_vs_kxkys',
        'vflx_vs_kxky' : 'vflux_vs_kxkys',
        'gvmus' : 'g2_vs_vpamus',
        'gzvs' : 'g2_vs_zvpas',
        } 
new_to_old_keys = {v: k for k, v in old_to_new_keys.items()}
temp = ['pflx_vs_kxky', 'vflx_vs_kxky', 'qflx_vs_kxky', 'pflux_x', 'vflux_x', 'qflux_x']

#-------------------------------------------------------------------------------
def compare_local_netcdf_quantity_to_expected_netcdf_quantity(local_netcdf_file, expected_netcdf_file, key, error=False):

    # Check whether the potential data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
    
        # Assume the key is the same in both files
        key1 = key
        key2 = key
    
        # Check whether the key is present
        if key not in local_netcdf.keys() and key not in expected_netcdf.keys():
            print(f'\nERROR: The key "{key}" does not exist in the netCDF files.')
            assert False, f'The key "{key}" does not exist in the netCDF files.'
        elif key not in local_netcdf.keys() and key in expected_netcdf.keys():
            key1 = old_to_new_keys[key]
            if key1 not in local_netcdf.keys():
                print(f'\nERROR: The key "{key}" is present in the expected netCDF file, but not')
                print(f'         in the local netcdf file. Check whether the diagnostics has been')
                print(f'         printed and whether the key of the diagnostics has changed.')
                assert False, f'The key "{key}" does not exist in the local netCDF file.'
        elif key in local_netcdf.keys() and key not in expected_netcdf.keys():
            key2 = new_to_old_keys[key]
            if key2 not in expected_netcdf.keys():
                print(f'\nERROR: The key "{key}" is present in the local netCDF file, but not')
                print(f'         in the expected netcdf file. This mightbe a new diagnostic,')
                print(f'         hence the expected output file should be updated.')
                assert False, f'The key "{key}" does not exist in the expected netCDF file.'
        
        # Read the quantity
        local_quantity = local_netcdf[key1]
        expected_quantity = expected_netcdf[key2] 
        
        # For the fluxes the first time step gives noise
        if ('flx' in key) or ('flux' in key):
            local_quantity[0] = 0
            expected_quantity[0] = 0 
        if key in ['upar', 'spitzer2', 'upar_x', 'spitzer2_x']:
            local_quantity[0] = 0
            expected_quantity[0] = 0 
        
        # Check the operating system
        system = platform.system()
        release = platform.release()
                     
        # Check whether the quantity matches
        if not (np.allclose(local_quantity, expected_quantity, rtol=1e-8, atol=1e-100)):
        
            # For the frequency, without nonlinear interactions, we have a lot of noise on the zonal modes
            if key in ['omega']: 
                if not (np.allclose(local_quantity, expected_quantity, rtol=1e-8, atol=1e-12)):
                    print(f'\nERROR: The {key} arrays do not match in the netCDF files.'); error = True
                    print(f'Compare the {key} arrays in the local and expected netCDF files:')                     
                    compare_local_array_with_expected_array(local_quantity, expected_quantity, name=key1)   
        
            # For the moments, there are noisy elements along (t,s,tube,zed,kx,ky,ri) so check the sum
            elif key in ['density', 'upar', 'temperature', 'spitzer2']: 
                sum1 = np.sum(np.abs(local_quantity.data[:,:,:,:,:,:,:]), axis=(1,2,3,4,5,6))
                sum2 = np.sum(np.abs(expected_quantity.data[:,:,:,:,:,:,:]), axis=(1,2,3,4,5,6))
                if not (np.allclose(sum1, sum2, rtol=1e-8, atol=1e-100)):
                    print(f'\nERROR: The {key} arrays do not match in the netCDF files.'); error = True
                    print(f'Compare the {key} arrays in the local and expected netCDF files:')
                    compare_local_array_with_expected_array(local_quantity, expected_quantity, name=key1)  
        
            # The fluxes tend to have more numerical noise
            elif key in ['qflux_vs_kxkyzs']:
                if not (np.allclose(local_quantity, expected_quantity, rtol=1e-7, atol=1e-100)):
                    print(f'\nERROR: The {key} arrays do not match in the netCDF files.'); error = True
                    print(f'Compare the {key} arrays in the local and expected netCDF files:')
                    compare_local_array_with_expected_array(local_quantity, expected_quantity, name=key1)
                
            else:
                print(f'\nERROR: The {key} arrays do not match in the netCDF files.'); error = True
                print(f'Compare the {key} arrays in the local and expected netCDF files:')
                compare_local_array_with_expected_array(local_quantity, expected_quantity, name=key1)   

        else:
            print(f'The {key1} matches!' )
            compare_local_array_with_expected_array(local_quantity, expected_quantity, name=key1)   
            
        assert (not error), f'\n\nThe {key} data does not match in the netCDF files.' 
        
    return error
    
        
#-------------------------------------------------------------------------------  
def compare_local_potential_with_expected_potential(local_netcdf_file='', expected_netcdf_file='', run_data={}, error=False): 

    # Make sure we have enough data to find the files 
    if (local_netcdf_file=='' or expected_netcdf_file=='') and len(run_data.keys())==0:
        print('\nERROR: The test has been set up badly in compare_local_potential_with_expected_potential().'); error = True
        print('\n       Make sure the path to the netCDF files is defined.')
        assert (not error), f'The test has been set up badly in compare_local_potential_with_expected_potential().' 
    
    # If the input file name is given, use it to set the netCDF files
    if (local_netcdf_file=='' or expected_netcdf_file=='') and len(run_data.keys())>0:
        input_file = run_data['input_file']; tmp_path = run_data['tmp_path']
        local_netcdf_file = tmp_path / input_file.replace('.in', '.out.nc') 
        expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file.replace(".in","")}.out.nc'    
    
    # Check whether the potential data matches in the netCDF file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        
        # Read the time axis
        local_time = local_netcdf['t']
        expected_time = expected_netcdf['t']
        
        # Read the potential axis
        local_phi2 = local_netcdf['phi2']
        expected_phi2 = expected_netcdf['phi2'] 
                     
        # Check whether we have the same time and potential data
        if not (np.allclose(local_time, expected_time, equal_nan=True, rtol=1e-05, atol=1e-20)):
            print('\nERROR: The time axis does not match in the netCDF files.'); error = True
            print('\nCompare the time arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_time, expected_time)  
        if not (np.allclose(local_phi2, expected_phi2, equal_nan=True, rtol=1e-05, atol=1e-20)):
            print('\nERROR: The potential data does not match in the netCDF files.'); error = True 
            print('\nCompare the potential arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_phi2, expected_phi2) 
        assert (not error), f'The potential data does not match in the netCDF files.' 
    
    return error
        
#-------------------------------------------------------------------------------  
def compare_local_potential_with_expected_potential_em(local_netcdf_file='', expected_netcdf_file='', run_data={}, error=False, check_apar=True, check_bpar=True): 

    # Make sure we have enough data to find the files 
    if (local_netcdf_file=='' or expected_netcdf_file=='') and len(run_data.keys())==0:
        print('\nERROR: The test has been set up badly in compare_local_potential_with_expected_potential().'); error = True
        print('\n       Make sure the path to the netCDF files is defined.')
        assert (not error), f'The test has been set up badly in compare_local_potential_with_expected_potential().' 
    
    # If the input file name is given, use it to set the netCDF files
    if (local_netcdf_file=='' or expected_netcdf_file=='') and len(run_data.keys())>0:
        input_file = run_data['input_file']; tmp_path = run_data['tmp_path']
        local_netcdf_file = tmp_path / input_file.replace('.in', '.out.nc') 
        expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file.replace(".in","")}.out.nc'    
    
    # Check whether the potential data matches in the netCDF file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf:
        
        # Read the time axis
        local_time = local_netcdf['t']
        expected_time = expected_netcdf['t']
        
        # Read the potential axis
        local_phi2 = local_netcdf['phi2']
        expected_phi2 = expected_netcdf['phi2'] 

        if check_apar:
           local_apar2 = local_netcdf['apar2']
           expected_apar2 = local_netcdf['apar2']

        if check_bpar:
           local_bpar2 = local_netcdf['bpar2']
           expected_bpar2 = local_netcdf['bpar2']
        
        # Check whether we have the same time and potential data
        if not (np.allclose(local_time, expected_time, equal_nan=True, atol=1e-20)):
            print('\nERROR: The time axis does not match in the netCDF files.'); error = True
            print('\nCompare the time arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_time, expected_time)
        if not (np.allclose(local_phi2, expected_phi2, equal_nan=True, atol=1e-20)):
            print('\nERROR: The <phi potential data does not match in the netCDF files.'); error = True 
            print('\nCompare the <phi> potential arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_phi2, expected_phi2)
        if check_apar:
            if not (np.allclose(local_apar2, expected_apar2, equal_nan=True, atol=1e-20)):
                print('\nERROR: The <A_parallel> potential data does not match in the netCDF files.'); error = True 
                print('\nCompare the <A_parallel> potential arrays in the local and expected netCDF files:')
                compare_local_array_with_expected_array(local_apar2, expected_apar2)
        if check_bpar:
            if not (np.allclose(local_bpar2, expected_bpar2, equal_nan=True, atol=1e-20)):
                print('\nERROR: The <B_parallel> potential data does not match in the netCDF files.'); error = True 
                print('\nCompare the <B_parallel> potential arrays in the local and expected netCDF files:')
                compare_local_array_with_expected_array(local_bar2, expected_bpar2)
        
        assert (not error), f'The potential data does not match in the netCDF files.' 
    
    return error
    
#-------------------------------------------------------------------------------
def convert_byte_array(array):
    '''Tool to convert netcdf text arrays, to a text string that we can compare.'''
    return '\n'.join((str(line, encoding='utf-8').strip() for line in array.data))
      
#-------------------------------------------------------------------------------
def compare_local_array_with_expected_array(local_array, expected_array, name='phi2'): 
    dimensions = len(np.shape(local_array.data)) 
    print(len(name)*' '+'                                   ')
    print(len(name)*' '+f'            LOCAL       EXPECTED   ({dimensions})')
    print(len(name)*' '+'            -----       --------   ') 
    if dimensions==0:
        print(f'    {name}  =  {local_array.data:18.9e}   {expected_array.data:18.9e}')
    elif dimensions==1:
        for i in range(np.min([len(local_array.data),10])):
            print(f' {name}[{i}]  =  {local_array.data[i]:18.9e}   {expected_array.data[i]:18.9e}') 
    elif dimensions==2: # bmag 
        for i in range(np.min([len(local_array.data),10])):
            print(f' {name}[{i}]  =  {local_array.data[i,0]:18.9e}   {expected_array.data[i,0]:18.9e}')
    elif dimensions==3:  
        for i in range(np.min([len(local_array.data),10])):
            print(f' {name}[{i}]  =  {local_array.data[i,0,0]:18.9e}   {expected_array.data[i,0,0]:18.9e}')
    elif dimensions==4: # kperp, vflux_vs_kxkys
        for i in range(np.min([len(local_array.data),10])):
            print(f' {name}[{i}]  =  {np.sum(np.abs(local_array.data[i,:,:,:])):18.9e}   {np.sum(np.abs(expected_array.data[i,:,:,:])):18.9e}') 
    elif dimensions==5: 
        for i in range(np.min([len(local_array.data),10])):
            print(f' {name}[{i}]  =  {np.sum(np.abs(local_array.data[i,:,:,:,:])):18.9e}   {np.sum(np.abs(expected_array.data[i,:,:,:,:])):18.9e}') 
    elif dimensions==6: 
        for i in range(np.min([len(local_array.data),10])):
            print(f' {name}[{i}]  =  {np.sum(np.abs(local_array.data[i,:,:,:,:,:])):18.9e}   {np.sum(np.abs(expected_array.data[i,:,:,:,:,:])):18.9e}') 
    elif dimensions==7: 
        for i in range(np.min([len(local_array.data),10])):
            print(f' {name}[{i}]  =  {np.sum(np.abs(local_array.data[i,:,:,:,:,:,:])):18.9e}   {np.sum(np.abs(expected_array.data[i,:,:,:,:,:,:])):18.9e}') 
    elif dimensions==8: 
        for i in range(np.min([len(local_array.data),10])):
            print(f' {name}[{i}]  =  {np.sum(np.abs(local_array.data[i,:,:,:,:,:,:,:])):18.9e}   {np.sum(np.abs(expected_array.data[i,:,:,:,:,:,:,:])):18.9e}') 
    elif dimensions==9: 
        for i in range(np.min([len(local_array.data),10])):
            print(f' {name}[{i}]  =  {np.sum(np.abs(local_array.data[i,:,:,:,:,:,:,:,:])):18.9e}   {np.sum(np.abs(expected_array.data[i,:,:,:,:,:,:,:,:])):18.9e}') 
    else:
        print('This has not been implemented.', name, dimensions)
        print(GoImplementThis)
    print('                                       ') 

    
################################################################################
#                  Routines to compare geometry in netCDF files                #
################################################################################
    
def compare_geometry_in_netcdf_files(run_data, error=False):
    
    # File names  
    input_file = run_data['input_file']; tmp_path = run_data['tmp_path']
    input_file = input_file.replace("_v0.5","").replace("_v0.6","").replace("_v0.7","")
    local_netcdf_file = tmp_path / input_file.replace('.in', '.out.nc') 
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_file.replace(".in","")}.out.nc' 
        
    # Check the operating system
    system = platform.system()
    release = platform.release()

    # Check whether the geometry data matches in the netcdf file
    with xr.open_dataset(local_netcdf_file) as local_netcdf, xr.open_dataset(expected_netcdf_file) as expected_netcdf: 
        
        # Relevant keys for the geometry
        geometry_keys = ["bmag", "b_dot_grad_z", "gradpar", "gbdrift", "gbdrift0", "cvdrift", "cvdrift0", "kperp2", \
            "gds2", "gds21", "gds22", "grho", "jacob", "djacdrho", "q", "shat", "d2qdr2", "drhodpsi", "d2psidr2", "jtwist"]  
        for key in geometry_keys:
        
            # Compare integers and floats
            if expected_netcdf[key].shape == ():
            
                # The number of processors is allowed to differ
                if key=='nproc': continue 
                
                # The order of the calculation of <drhodpsi> was changed, so there is some numerical noise 
                if key=='drhodpsi': 
                    local_netcdf[key].data = np.round(local_netcdf[key].data, 12)
                    expected_netcdf[key].data = np.round(expected_netcdf[key].data, 12)
                    
                # Compare integers and floats
                if (local_netcdf[key] != expected_netcdf[key]):
                    print(f'\nERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                    print(f'    LOCAL:    {local_netcdf[key].data}')
                    print(f'    EXPECTED: {expected_netcdf[key].data}')
                
            # Compare texts (code_info and input_file)
            elif expected_netcdf[key].dtype.kind == 'S':
                expected_netcdf_str = convert_byte_array(expected_netcdf[key])
                local_netcdf_str = convert_byte_array(local_netcdf[key])  
                if key=='input_file': continue # The input file is allowed to differ
                if (local_netcdf_str != expected_netcdf_str):
                    print(f'\nERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                    print(f'    LOCAL:    {local_netcdf_str}')
                    print(f'    EXPECTED: {expected_netcdf_str}') 
                    
            # Compare data arrays 
            else: 
                
                # Changed definitions
                if key in ["gbdrift0", "cvdrift0", "gbdrift", "cvdrift", "gds22", "gds21", "gds2", "b_dot_grad_z", "gradpar"]:
                    if (key=="b_dot_gradz_avg"):
                        b_dot_gradz_old = expected_netcdf["gradpar"]
                        b_dot_gradz_new = local_netcdf["b_dot_gradz_avg"]
                        if not (np.allclose(b_dot_gradz_old, b_dot_gradz_new, equal_nan=True)):
                            print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                            print('\nCompare the {key} arrays in the local and expected netCDF files:')
                            compare_local_array_with_expected_array(b_dot_gradz_old, b_dot_gradz_new)
                    if (key=="b_dot_gradz"):
                        b_dot_gradz_old = expected_netcdf["b_dot_grad_z"]
                        b_dot_gradz_new = local_netcdf["b_dot_gradz"]
                        if not (np.allclose(b_dot_gradz_old, b_dot_gradz_new, equal_nan=True)):
                            print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                            print('\nCompare the {key} arrays in the local and expected netCDF files:')
                            compare_local_array_with_expected_array(b_dot_gradz_old, b_dot_gradz_new)
                    if (key=="gds2"):
                        gds2_old = expected_netcdf["gds2"]
                        gds2_new = local_netcdf["grady_dot_grady"]
                        if not (np.allclose(gds2_old, gds2_new, equal_nan=True)):
                            print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                            print('\nCompare the {key} arrays in the local and expected netCDF files:')
                            compare_local_array_with_expected_array(gds2_old, gds2_new)
                    if (key=="gds21"):
                        gds21_old =  expected_netcdf["gds21"]
                        gradx_dot_grady_new = local_netcdf["gradx_dot_grady"]
                        gds21_new = gradx_dot_grady_new * local_netcdf["shat"]
                        if not (np.allclose(gds21_old, gds21_new, equal_nan=True)):
                            print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                            print('\nCompare the {key} arrays in the local and expected netCDF files:')
                            compare_local_array_with_expected_array(gds21_old, gds21_new)
                    if (key=="gds22"):
                        gds22_old =  expected_netcdf["gds22"]
                        gradx_dot_gradx_new = local_netcdf["gradx_dot_gradx"]
                        gds22_new = gradx_dot_gradx_new * local_netcdf["shat"] * local_netcdf["shat"]
                        if not (np.allclose(gds22_old, gds22_new, equal_nan=True)):
                            print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                            print('\nCompare the {key} arrays in the local and expected netCDF files:')
                            compare_local_array_with_expected_array(gds22_old, gds22_new)
                    if (key=="gbdrift"):
                        gbdrift_old =  expected_netcdf["gbdrift"]
                        B_times_gradB_dot_grady_new = local_netcdf["B_times_gradB_dot_grady"]
                        gbdrift_new = B_times_gradB_dot_grady_new * 2
                        if not (np.allclose(gbdrift_old, gbdrift_new, equal_nan=True)):
                            print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                            print('\nCompare the {key} arrays in the local and expected netCDF files:')
                            compare_local_array_with_expected_array(gbdrift_old, gbdrift_new)
                    if (key=="cvdrift"):
                        cvdrift_old =  expected_netcdf["cvdrift"]
                        B_times_kappa_dot_grady_new = local_netcdf["B_times_kappa_dot_grady"]
                        cvdrift_new = B_times_kappa_dot_grady_new * 2
                        if not (np.allclose(cvdrift_old, cvdrift_new, equal_nan=True)):
                            print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                            print('\nCompare the {key} arrays in the local and expected netCDF files:')
                            compare_local_array_with_expected_array(cvdrift_old, cvdrift_new)
                    if (key=="gbdrift0"):
                        gbdrift0_old =  expected_netcdf["gbdrift0"]
                        B_times_gradB_dot_gradx_new = local_netcdf["B_times_gradB_dot_gradx"] 
                        gbdrift0_new = B_times_gradB_dot_gradx_new * 2 * local_netcdf["shat"]
                        if not (np.allclose(gbdrift0_old, gbdrift0_new, equal_nan=True)):
                            print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                            print('\nCompare the {key} arrays in the local and expected netCDF files:')
                            compare_local_array_with_expected_array(gbdrift0_old, gbdrift0_new)
                    elif (key=="cvdrift0"):
                        cvdrift0_old =  expected_netcdf["cvdrift0"]
                        B_times_kappa_dot_gradx_new = local_netcdf["B_times_kappa_dot_gradx"] 
                        cvdrift0_new = B_times_kappa_dot_gradx_new * 2 * local_netcdf["shat"]
                        if not (np.allclose(cvdrift0_old, cvdrift0_new, equal_nan=True)):
                            print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                            print('\nCompare the {key} arrays in the local and expected netCDF files:')
                            compare_local_array_with_expected_array(cvdrift0_old, cvdrift0_new)
            
                # Compare data arrays 
                else:
                    if not (np.allclose(local_netcdf[key], expected_netcdf[key], equal_nan=True)):
                        print(f'ERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                        print('\nCompare the {key} arrays in the local and expected netCDF files:')
                        compare_local_array_with_expected_array(local_netcdf[key], expected_netcdf[key])  
                    
        # Print "AssertionError: <error message>" if an error was encountered
        assert (not error), f'Some Miller geometry arrays in the netcdf file did not match the previous run.'  
        
    return error
    
    
################################################################################
#                  Routines to compare geometry in Miller files                #
################################################################################

def compare_geometry_files(local_geometry_file, expected_geometry_file, error=False, with_btor=True, digits=2):

    def process_error(variable):
        print(f'\nERROR: {variable} does not match in the *.geometry file.\n'); 
        return True

    # Read variables in old *.geometry file
    with open(expected_geometry_file, "r") as f:
       lines = f.readlines()
       variables = lines[1]
    variables = variables.replace('#','').replace('\n','').split(' ')
    variables = [i for i in variables if i!='']
    rhoc_old = float(variables[0])
    qinp_old = float(variables[1])
    shat_old = float(variables[2])
    rhotor_old = float(variables[3])
    aref_old = float(variables[4])
    bref_old = float(variables[5])
    dxdpsi_old = float(variables[6])
    dydalpha_old = float(variables[7])
    exb_nonlin_old = float(variables[8])
    exb_nonlin_p_old = float(variables[9])
    
    # Read arrays in old *.geometry file
    if with_btor: data = np.loadtxt(expected_geometry_file,skiprows=4,dtype='float').reshape(-1, 15)
    if not with_btor: data = np.loadtxt(expected_geometry_file,skiprows=4,dtype='float').reshape(-1, 14)
    alpha_old = data[:,0]
    zed_old = data[:,1]
    zeta_old = data[:,2]
    bmag_old = data[:,3]
    b_dot_gradz_old = data[:,4]
    gds2_old = data[:,5]
    gds21_old = data[:,6]
    gds22_old = data[:,7]
    gds23_old = data[:,8]
    gds24_old = data[:,9]
    gbdrift_old = data[:,10]
    cvdrift_old = data[:,11]
    gbdrift0_old = data[:,12]
    bmag_psi0_old = data[:,13]
    if with_btor: btor_old = data[:,14]
    
    # Read variables in new *.geometry file
    with open(local_geometry_file, "r") as f:
       lines = f.readlines()
       variables = lines[1]
    variables = variables.replace('#','').replace('\n','').split(' ')
    variables = [i for i in variables if i!='']
    rhoc_new = float(variables[0])
    qinp_new = float(variables[1])
    shat_new = float(variables[2])
    rhotor_new = float(variables[3])
    aref_new = float(variables[4])
    bref_new = float(variables[5])
    dxdpsi_new = float(variables[6])
    dydalpha_new = float(variables[7])
    exb_nonlin_new = float(variables[8])
    fluxfac_new = float(variables[9])
    
    # Read arrays in new *.geometry file
    if with_btor: data = np.loadtxt(local_geometry_file,skiprows=4,dtype='float').reshape(-1, 15)
    if not with_btor: data = np.loadtxt(local_geometry_file,skiprows=4,dtype='float').reshape(-1, 15)
    alpha_new = data[:,0]
    zed_new = data[:,1]
    zeta_new = data[:,2]
    bmag_new = data[:,3]
    b_dot_gradz_new = data[:,4]
    grady_dot_grady_new = data[:,5]
    gradx_dot_grady_new = data[:,6]
    gradx_dot_gradx_new = data[:,7]
    gds23_new = data[:,8]
    gds24_new = data[:,9]
    B_times_gradB_dot_grady_new = data[:,10]
    B_times_kappa_dot_grady_new = data[:,11]
    B_times_gradB_dot_gradx_new = data[:,12]
    bmag_psi0_new = data[:,13]
    if with_btor: btor_new = data[:,14]
    
    # New definitions
    gds2_new = grady_dot_grady_new
    gbdrift0_new = B_times_gradB_dot_gradx_new * 2 * shat_new
    gbdrift0_new = np.round(gbdrift0_new, digits)
    gbdrift0_old = np.round(gbdrift0_old, digits)
    cvdrift_new = B_times_kappa_dot_grady_new * 2
    cvdrift_new = np.round(cvdrift_new, digits)
    cvdrift_old = np.round(cvdrift_old, digits)
    gds22_new = gradx_dot_gradx_new * shat_new * shat_new
    gds22_new = np.round(gds22_new, digits)
    gds22_old = np.round(gds22_old, digits)
    if digits >= 2: digits = digits - 1
    gbdrift_new = B_times_gradB_dot_grady_new * 2
    gbdrift_new = np.round(gbdrift_new, digits)
    gbdrift_old = np.round(gbdrift_old, digits)
    if digits >= 1: digits = digits - 1
    gds21_new = gradx_dot_grady_new * shat_new
    gds21_new = np.round(gds21_new, digits)
    gds21_old = np.round(gds21_old, digits)
    
    # Compare values
    if not (np.allclose(rhoc_old, rhoc_new, equal_nan=True)): error = process_error('rhoc')
    if not (np.allclose(qinp_old, qinp_new, equal_nan=True)): error = process_error('qinp')
    if not (np.allclose(shat_old, shat_new, equal_nan=True)): error = process_error('shat')
    if not (np.allclose(aref_old, aref_new, equal_nan=True)): error = process_error('aref')
    if not (np.allclose(bref_old, bref_new, equal_nan=True)): error = process_error('bref')
    if not (np.allclose(dxdpsi_old, dxdpsi_new, equal_nan=True)): error = process_error('dxdpsi')
    if not (np.allclose(dydalpha_old, dydalpha_new, equal_nan=True)): error = process_error('dydalpha')
    if not (np.allclose(exb_nonlin_old, exb_nonlin_new, equal_nan=True)): error = process_error('exb_nonlin')
    
    # Do not compare rhotor, it was defined badly in the past
    #if not (np.allclose(rhotor_old, rhotor_new, equal_nan=True)): error = process_error('rhotor')
    
    # Compare arrays
    if not (np.allclose(alpha_old, alpha_new, equal_nan=True)): error = process_error('alpha')
    if not (np.allclose(zed_old, zed_new, equal_nan=True)): error = process_error('zed')
    if not (np.allclose(zeta_old, zeta_new, equal_nan=True)): error = process_error('zeta')
    if not (np.allclose(bmag_old, bmag_new, equal_nan=True)): error = process_error('bmag')
    if not (np.allclose(b_dot_gradz_new, b_dot_gradz_new, equal_nan=True)): error = process_error('b_dot_gradz')
    if not (np.allclose(gds2_old, gds2_new, equal_nan=True)): error = process_error('gds2')
    if not (np.allclose(gds21_old, gds21_new, equal_nan=True)): error = process_error('gds21')
    if not (np.allclose(gbdrift_old, gbdrift_new, equal_nan=True)): error = process_error('gbdrift')
    if not (np.allclose(cvdrift_old, cvdrift_new, equal_nan=True)): error = process_error('cvdrift')
    if not (np.allclose(gbdrift0_old, gbdrift0_new, equal_nan=True)): error = process_error('gbdrift0')
    if not (np.allclose(bmag_psi0_old, bmag_psi0_new, equal_nan=True)): error = process_error('bmag_psi0')
    if with_btor: 
        if not (np.allclose(btor_old, btor_new, equal_nan=True)): error = process_error('btor')
    assert (not error), f'The geometry data does not match in the *.geometry file.'
    
    # Do not compare gds23 and gds24 it was badly defined in Miller and VMEC
    #if not (np.allclose(gds23_old, gds23_new, equal_nan=True)): error = process_error('gds23')
    #if not (np.allclose(gds24_old, gds24_new, equal_nan=True)): error = process_error('gds24')
    return

#-------------------------------------------------------------------------------  
def compare_miller_input_files(local_file, expected_file, error=False):

    def process_error(variable):
        print(f'\nERROR: {variable} does not match in the *.millerlocal.input file.\n'); 
        return True
    
    # Read variables in old *.geometry file
    with open(expected_file, "r") as f:
       lines = f.readlines()
       variables1 = lines[1]
       variables2 = lines[4]
       variables3 = lines[7]
       variables4 = lines[10]
    variables1 = variables1.replace('#','').replace('\n','').split(' ')
    variables2 = variables2.replace('#','').replace('\n','').split(' ')
    variables3 = variables3.replace('#','').replace('\n','').split(' ')
    variables4 = variables4.replace('#','').replace('\n','').split(' ')
    variables1 = [i for i in variables1 if i!='']
    variables2 = [i for i in variables2 if i!='']
    variables3 = [i for i in variables3 if i!='']
    variables4 = [i for i in variables4 if i!='']
    rhoc_old = float(variables1[0])
    rmaj_old = float(variables1[1])
    rgeo_old = float(variables1[2])
    shift_old = float(variables1[3])
    qinp_old = float(variables1[4])
    shat_old = float(variables2[0])
    kappa_old = float(variables2[1])
    kapprim_old = float(variables2[2])
    tri_old = float(variables2[3])
    triprim_old = float(variables2[4])
    betaprim_old = float(variables3[0])
    dpsitordrho_old = float(variables3[1])
    rhotor_old = float(variables3[2])
    drhotordrho_old = float(variables3[3])
    d2qdr2_old = float(variables3[4])
    d2psidr2_old = float(variables4[0])
    betadbprim_old = float(variables4[1])
    psitor_lcfs_old = float(variables4[2])

    # Read variables in new *.geometry file
    with open(local_file, "r") as f:
       lines = f.readlines()
       variables1 = lines[1]
       variables2 = lines[4]
       variables3 = lines[7]
       variables4 = lines[10]
    variables1 = variables1.replace('#','').replace('\n','').split(' ')
    variables2 = variables2.replace('#','').replace('\n','').split(' ')
    variables3 = variables3.replace('#','').replace('\n','').split(' ')
    variables4 = variables4.replace('#','').replace('\n','').split(' ')
    variables1 = [i for i in variables1 if i!='']
    variables2 = [i for i in variables2 if i!='']
    variables3 = [i for i in variables3 if i!='']
    variables4 = [i for i in variables4 if i!='']
    rhoc_new = float(variables1[0])
    rmaj_new = float(variables1[1])
    rgeo_new = float(variables1[2])
    shift_new = float(variables1[3])
    qinp_new = float(variables1[4])
    shat_new = float(variables2[0])
    kappa_new = float(variables2[1])
    kapprim_new = float(variables2[2])
    tri_new = float(variables2[3])
    triprim_new = float(variables2[4])
    betaprim_new = float(variables3[0])
    dpsitordrho_new = float(variables3[1])
    rhotor_new = float(variables3[2])
    drhotordrho_new = float(variables3[3])
    d2qdr2_new = float(variables3[4])
    d2psidr2_new = float(variables4[0])
    betadbprim_new = float(variables4[1])
    psitor_lcfs_new = float(variables4[2])
    
    # Compare variables
    if not (np.allclose(rhoc_old, rhoc_new, equal_nan=True)): error = process_error('rhoc')
    if not (np.allclose(rmaj_old, rmaj_new, equal_nan=True)): error = process_error('rmaj')
    if not (np.allclose(rgeo_old, rgeo_new, equal_nan=True)): error = process_error('rgeo')
    if not (np.allclose(shift_old, shift_new, equal_nan=True)): error = process_error('shift')
    if not (np.allclose(qinp_old, qinp_new, equal_nan=True)): error = process_error('qinp')
    if not (np.allclose(shat_old, shat_new, equal_nan=True)): error = process_error('shat')
    if not (np.allclose(kappa_old, kappa_new, equal_nan=True)): error = process_error('kappa')
    if not (np.allclose(kapprim_old, kapprim_new, equal_nan=True)): error = process_error('kapprim')
    if not (np.allclose(tri_old, tri_new, equal_nan=True)): error = process_error('tri')
    if not (np.allclose(betaprim_old, betaprim_new, equal_nan=True)): error = process_error('betaprim')
    if not (np.allclose(dpsitordrho_old, dpsitordrho_new, equal_nan=True)): error = process_error('dpsitordrho')
    if not (np.allclose(drhotordrho_old, drhotordrho_new, equal_nan=True)): error = process_error('drhotordrho')
    if not (np.allclose(d2qdr2_old, d2qdr2_new, equal_nan=True)): error = process_error('d2qdr2')
    if not (np.allclose(d2psidr2_old, d2psidr2_new, equal_nan=True)): error = process_error('d2psidr2')
    if not (np.allclose(betadbprim_old, betadbprim_new, equal_nan=True)): error = process_error('betadbprim')
    if not (np.allclose(psitor_lcfs_old, psitor_lcfs_new, equal_nan=True)): error = process_error('psitor_lcfs')
    assert (not error), f'The geometry data does not match in the *.millerlocal.input file.'
    
    # Do not compare rhotor, it was defined badly in the past
    #if not (np.allclose(rhotor_old, rhotor_new, equal_nan=True)): error = process_error('rhotor')
    return shat_new
    
#-------------------------------------------------------------------------------  
def compare_miller_output_files(local_file, expected_file, shat, error=False):

    def process_error(variable):
        print(f'\nERROR: {variable} does not match in the *.millerlocal.output file.\n'); 
        return True

    # Read variables in old *.geometry file
    with open(expected_file, "r") as f:
       lines = f.readlines()
       variables = lines[0]
    variables = variables.replace('dI/dr:','').replace('d2I/dr2:','').replace('dpsi/dr:','')
    variables = variables.replace('#','').replace('\n','').split(' ')
    variables = [i for i in variables if i!='']
    dI_dr_old = float(variables[0])
    d2I_dr2_old = float(variables[1])
    dpsi_dr_old = float(variables[2])
    
    # Read arrays in old *.geometry file
    data = np.loadtxt(expected_file,skiprows=2,dtype='float').reshape(-1, 58)
    theta_old = data[:,0]
    R_old = data[:,1]
    dRdr_old = data[:,2]
    d2Rdr2_old = data[:,3]
    dR_dth_old = data[:,4]
    d2Rdrdth_old = data[:,5]
    dZ_dr_old = data[:,6]
    d2Zdr2_old = data[:,7]
    dZ_dth_old = data[:,8]
    d2Zdrdth_old = data[:,9]
    bmag_old = data[:,10]
    dBdr_old = data[:,11]
    d2Bdr2_old = data[:,12]
    dB_dth_old = data[:,13]
    d2Bdrdth_old = data[:,14]
    varthet_old = data[:,15]
    dvarthdr_old = data[:,16]
    d2varthdr2_old = data[:,17]
    jacr_old = data[:,18]
    djacrdr_old = data[:,19]
    djacdrho_old = data[:,20]
    d2jacdr2_old = data[:,21]
    grho2_old = data[:,22]
    dgr2dr_old = data[:,23]
    gthet2_old = data[:,24]
    dgt2_old = data[:,25]
    grgthet_old = data[:,26]
    dgrgt_old = data[:,27]
    galphgth_old = data[:,28]
    dgagt_old = data[:,29]
    grgalph_old = data[:,30]
    dgagr_old = data[:,31]
    galph2_old = data[:,32]
    dga2_old = data[:,33]
    cross_old = data[:,34]
    dcrossdr_old = data[:,35]
    gbdrift0_old = data[:,36]
    dgbdrift0_old = data[:,37]
    cvdrift0_old = data[:,38]
    dcvdrift0_old = data[:,39]
    gbdrift_old = data[:,40]
    dgbdrift_old = data[:,41]
    cvdrift_old = data[:,42]
    dcvdrift_old = data[:,43]
    drzdth_old = data[:,44]
    b_dot_gradz_old = data[:,45]
    dgpardr_old = data[:,46]
    b_dot_gradzB_old = data[:,47]
    dgparBdr_old = data[:,48]
    gds2_old = data[:,49]
    dgds2dr_old = data[:,50]
    gds21_old = data[:,51]
    dgds21dr_old = data[:,52]
    gds22_old = data[:,53]
    dgds22dr_old = data[:,54]
    gds23_old = data[:,55]
    gds24_old = data[:,56]
    Zr_old = data[:,57]

    # Read variables in new *.geometry file
    with open(local_file, "r") as f:
       lines = f.readlines()
       variables = lines[0]
    variables = variables.replace('dI/dr:','').replace('d2I/dr2:','').replace('dpsi/dr:','')
    variables = variables.replace('#','').replace('\n','').split(' ')
    variables = [i for i in variables if i!='']
    dI_dr_new = float(variables[0])
    d2I_dr2_new = float(variables[1])
    dpsi_dr_new = float(variables[2])
    
    # Read arrays in new *.geometry file
    data = np.loadtxt(local_file,skiprows=2,dtype='float').reshape(-1, 58)
    theta_new = data[:,0]
    R_new = data[:,1]
    dRdr_new = data[:,2]
    d2Rdr2_new = data[:,3]
    dR_dth_new = data[:,4]
    d2Rdrdth_new = data[:,5]
    dZ_dr_new = data[:,6]
    d2Zdr2_new = data[:,7]
    dZ_dth_new = data[:,8]
    d2Zdrdth_new = data[:,9]
    bmag_new = data[:,10]
    dBdr_new = data[:,11]
    d2Bdr2_new = data[:,12]
    dB_dth_new = data[:,13]
    d2Bdrdth_new = data[:,14]
    varthet_new = data[:,15]
    dvarthdr_new = data[:,16]
    d2varthdr2_new = data[:,17]
    jacr_new = data[:,18]
    djacrdr_new = data[:,19]
    djacdrho_new = data[:,20]
    d2jacdr2_new = data[:,21]
    grho2_new = data[:,22]
    dgr2dr_new = data[:,23]
    gthet2_new = data[:,24]
    dgt2_new = data[:,25]
    grgthet_new = data[:,26]
    dgrgt_new = data[:,27]
    galphgth_new = data[:,28]
    dgagt_new = data[:,29]
    grgalph_new = data[:,30]
    dgagr_new = data[:,31]
    galph2_new = data[:,32]
    dga2_new = data[:,33]
    cross_new = data[:,34]
    dcrossdr_new = data[:,35]
    B_times_gradB_dot_gradx_new = data[:,36]
    dgbdrift0_new = data[:,37]
    B_times_kappa_dot_gradx_new = data[:,38]
    dcvdrift0_new = data[:,39]
    B_times_gradB_dot_grady_new = data[:,40]
    dgbdrift_new = data[:,41]
    B_times_kappa_dot_grady_new = data[:,42]
    dcvdrift_new = data[:,43]
    drzdth_new = data[:,44]
    b_dot_gradz_new = data[:,45]
    dgpardr_new = data[:,46]
    b_dot_gradzB_new = data[:,47]
    dgparBdr_new = data[:,48]
    grady_dot_grady_new = data[:,49]
    dgds2dr_new = data[:,50]
    gradx_dot_grady_new = data[:,51]
    dgds21dr_new = data[:,52]
    gradx_dot_gradx_new = data[:,53]
    dgds22dr_new = data[:,54]
    gds23_new = data[:,55]
    gds24_new = data[:,56]
    Zr_new = data[:,57]
    
    # New definitions
    gbdrift0_new = B_times_gradB_dot_gradx_new * 2 * shat
    cvdrift0_new = B_times_kappa_dot_gradx_new * 2 * shat
    gbdrift_new = B_times_gradB_dot_grady_new * 2
    cvdrift_new = B_times_kappa_dot_grady_new * 2
    gds22_new = gradx_dot_gradx_new * shat * shat
    gds21_new = gradx_dot_grady_new * shat
    gds2_new = grady_dot_grady_new
    
    # Compare variables
    if not (np.allclose(dI_dr_old, dI_dr_new, equal_nan=True)): error = process_error('dI_dr')
    if not (np.allclose(d2I_dr2_old, d2I_dr2_new, equal_nan=True)): error = process_error('d2I_dr2')
    if not (np.allclose(dpsi_dr_old, dpsi_dr_new, equal_nan=True)): error = process_error('dpsi_dr')
           
    # Compare arrays
    if not (np.allclose(theta_old, theta_new, equal_nan=True)): error = process_error('theta')
    if not (np.allclose(R_old, R_new, equal_nan=True)): error = process_error('R')
    if not (np.allclose(dRdr_old, dRdr_new, equal_nan=True)): error = process_error('dR_dr')
    if not (np.allclose(d2Rdr2_old, d2Rdr2_new, equal_nan=True)): error = process_error('d2Rdr2')
    if not (np.allclose(dR_dth_old, dR_dth_new, equal_nan=True)): error = process_error('dR_dth')
    if not (np.allclose(d2Rdrdth_old, d2Rdrdth_new, equal_nan=True)): error = process_error('d2Rdrdth')
    if not (np.allclose(dZ_dr_old, dZ_dr_new, equal_nan=True)): error = process_error('dZ_dr')
    if not (np.allclose(d2Zdr2_old, d2Zdr2_new, equal_nan=True)): error = process_error('d2Zdr2')
    if not (np.allclose(dZ_dth_old, dZ_dth_new, equal_nan=True)): error = process_error('dZ_dth')
    if not (np.allclose(d2Zdrdth_old, d2Zdrdth_new, equal_nan=True)): error = process_error('d2Zdrdth')
    if not (np.allclose(bmag_old, bmag_new, equal_nan=True)): error = process_error('bmag')
    if not (np.allclose(dBdr_old, dBdr_new, equal_nan=True)): error = process_error('dBdr')
    if not (np.allclose(d2Bdr2_old, d2Bdr2_new, equal_nan=True)): error = process_error('d2Bdr2')
    if not (np.allclose(dB_dth_old, dB_dth_new, equal_nan=True)): error = process_error('dB_dth')
    if not (np.allclose(d2Bdrdth_old, d2Bdrdth_new, equal_nan=True)): error = process_error('d2Bdrdth')
    if not (np.allclose(varthet_old, varthet_new, equal_nan=True)): error = process_error('varthet')
    if not (np.allclose(dvarthdr_old, dvarthdr_new, equal_nan=True)): error = process_error('dvarthdr')
    if not (np.allclose(d2varthdr2_old, d2varthdr2_new, equal_nan=True)): error = process_error('d2varthdr2')  
    if not (np.allclose(jacr_old, jacr_new, equal_nan=True)): error = process_error('jacr')
    if not (np.allclose(djacrdr_old, djacrdr_new, equal_nan=True)): error = process_error('djacrdr')
    if not (np.allclose(djacdrho_old, djacdrho_new, equal_nan=True)): error = process_error('djacdrho')
    if not (np.allclose(d2jacdr2_old, d2jacdr2_new, equal_nan=True)): error = process_error('d2jacdr2')
    if not (np.allclose(grho2_old, grho2_new, equal_nan=True)): error = process_error('grho2')
    if not (np.allclose(dgr2dr_old, dgr2dr_new, equal_nan=True)): error = process_error('dgr2dr')
    if not (np.allclose(gthet2_old, gthet2_new, equal_nan=True)): error = process_error('gthet2')
    if not (np.allclose(dgt2_old, dgt2_new, equal_nan=True)): error = process_error('dgt2')
    if not (np.allclose(grgthet_old, grgthet_new, equal_nan=True)): error = process_error('grgthet')
    if not (np.allclose(dgrgt_old, dgrgt_new, equal_nan=True)): error = process_error('dgrgt')
    if not (np.allclose(galphgth_old, galphgth_new, equal_nan=True)): error = process_error('galphgth')
    if not (np.allclose(dgagt_old, dgagt_new, equal_nan=True)): error = process_error('dgagt')
    if not (np.allclose(grgalph_old, grgalph_new, equal_nan=True)): error = process_error('grgalph')
    if not (np.allclose(dgagr_old, dgagr_new, equal_nan=True)): error = process_error('dgagr')
    if not (np.allclose(galph2_old, galph2_new, equal_nan=True)): error = process_error('galph2')
    if not (np.allclose(dga2_old, dga2_new, equal_nan=True)): error = process_error('dga2')
    if not (np.allclose(cross_old, cross_new, equal_nan=True)): error = process_error('cross')
    if not (np.allclose(dcrossdr_old, dcrossdr_new, equal_nan=True)): error = process_error('dcrossdr')
    if not (np.allclose(gbdrift0_old, gbdrift0_new, equal_nan=True)): error = process_error('gbdrift0')
    if not (np.allclose(dgbdrift0_old, dgbdrift0_new, equal_nan=True)): error = process_error('dgbdrift0')
    if not (np.allclose(cvdrift0_old, cvdrift0_new, equal_nan=True)): error = process_error('cvdrift0')
    if not (np.allclose(dcvdrift0_old, dcvdrift0_new, equal_nan=True)): error = process_error('dcvdrift0')
    if not (np.allclose(gbdrift_old, gbdrift_new, equal_nan=True)): error = process_error('gbdrift')
    if not (np.allclose(dgbdrift_old, dgbdrift_new, equal_nan=True)): error = process_error('dgbdrift')
    if not (np.allclose(cvdrift_old, cvdrift_new, equal_nan=True)): error = process_error('cvdrift')
    if not (np.allclose(dcvdrift_old, dcvdrift_new, equal_nan=True)): error = process_error('dcvdrift')
    if not (np.allclose(drzdth_old, drzdth_new, equal_nan=True)): error = process_error('drzdth')
    if not (np.allclose(b_dot_gradz_old, b_dot_gradz_new, equal_nan=True)): error = process_error('b_dot_gradz')
    if not (np.allclose(dgpardr_old, dgpardr_new, equal_nan=True)): error = process_error('dgpardr')
    if not (np.allclose(b_dot_gradzB_old, b_dot_gradzB_new, equal_nan=True)): error = process_error('b_dot_gradzB')
    if not (np.allclose(dgparBdr_old, dgparBdr_new, equal_nan=True)): error = process_error('dgparBdr')
    if not (np.allclose(gds2_old, gds2_new, equal_nan=True)): error = process_error('gds2')
    if not (np.allclose(dgds2dr_old, dgds2dr_new, equal_nan=True)): error = process_error('dgds2dr')
    if not (np.allclose(gds21_old, gds21_new, equal_nan=True)): error = process_error('gds21')
    if not (np.allclose(dgds21dr_old, dgds21dr_new, equal_nan=True)): error = process_error('dgds21dr')
    if not (np.allclose(gds22_old, gds22_new, equal_nan=True)): error = process_error('gds22')
    if not (np.allclose(dgds22dr_old, dgds22dr_new, equal_nan=True)): error = process_error('dgds22dr')
    if not (np.allclose(gds23_old, gds23_new, equal_nan=True)): error = process_error('gds23')
    if not (np.allclose(gds24_old, gds24_new, equal_nan=True)): error = process_error('gds24')
    if not (np.allclose(Zr_old, Zr_new, equal_nan=True)): error = process_error('Zr')
    assert (not error), f'The geometry data does not match in the *.millerlocal.output file.'
    return
    
    
#-------------------------------------------------------------------------------
def compare_vmecgeo_files(local_geometry_file, expected_geometry_file, error=False):

    def process_error(variable):
        print(f'\nERROR: {variable} does not match in the *.vmec.geo file.\n'); 
        return True

    # Read variables in old *.geometry file
    with open(expected_geometry_file, "r") as f:
       lines = f.readlines()
       variables = lines[1]
    variables = variables.replace('#','').replace('\n','').split(' ')
    variables = [i for i in variables if i!='']
    rhotor_old = float(variables[0])
    qinp_old = float(variables[1])
    shat_old = float(variables[2])
    aref_old = float(variables[3])
    bref_old = float(variables[4])
    z_scalefac_old = float(variables[5])
    
    # Read arrays in old *.geometry file
    data = np.loadtxt(expected_geometry_file,skiprows=4,dtype='float').reshape(-1, 17)
    alpha_old = data[:,0]
    zeta_old = data[:,1]
    bmag_old = data[:,2]
    b_dot_gradz_old = data[:,3]
    bdot_grad_z_old = data[:,4]
    grad_alpha2_old = data[:,5]
    gd_alph_psi_old = data[:,6]
    grad_psi2_old = data[:,7]
    gds23_old = data[:,8]
    gds24_old = data[:,9]
    gbdrift_old = data[:,10]
    gbdrift0_old = data[:,11]
    cvdrift_old = data[:,12]
    cvdrift0_old = data[:,13]
    theta_vmec_old = data[:,14]
    B_sub_theta_old = data[:,15]
    B_sub_zeta_old = data[:,16]
    
    # Read variables in new *.geometry file
    with open(local_geometry_file, "r") as f:
       lines = f.readlines()
       variables = lines[1]
    variables = variables.replace('#','').replace('\n','').split(' ')
    variables = [i for i in variables if i!='']
    rhotor_new = float(variables[0])
    qinp_new = float(variables[1])
    shat_new = float(variables[2])
    aref_new = float(variables[3])
    bref_new = float(variables[4])
    z_scalefac_new = float(variables[5])
    
    # Read arrays in new *.geometry file
    data = np.loadtxt(local_geometry_file,skiprows=4,dtype='float').reshape(-1, 17)
    alpha_new = data[:,0]
    zeta_new = data[:,1]
    bmag_new = data[:,2]
    b_dot_gradz_new = data[:,3]
    bdot_grad_z_new = data[:,4]
    grad_alpha2_new = data[:,5]
    gd_alph_psi_new = data[:,6]
    grad_psi2_new = data[:,7]
    gds23_new = data[:,8]
    gds24_new = data[:,9]
    B_times_gradB_dot_grady_new = data[:,10]
    B_times_gradB_dot_gradx_psi_new = data[:,11]
    B_times_kappa_dot_grady_new = data[:,12]
    B_times_kappa_dot_gradx_psi_new = data[:,13]
    theta_vmec_new = data[:,14]
    B_sub_theta_new = data[:,15]
    B_sub_zeta_new = data[:,16]
    
    # New definitions (psi changed to -psi)
    gd_alph_psi_new = -gd_alph_psi_new
    B_times_gradB_dot_gradx_psi_new = - B_times_gradB_dot_gradx_psi_new
    B_times_kappa_dot_gradx_psi_new = - B_times_kappa_dot_gradx_psi_new
    
    # New definitions
    digits = 4
    gbdrift0_new = B_times_gradB_dot_gradx_psi_new * 2 * shat_new
    gbdrift0_new = np.round(gbdrift0_new, digits)
    gbdrift0_old = np.round(gbdrift0_old, digits)
    cvdrift0_new = B_times_kappa_dot_gradx_psi_new * 2 * shat_new
    cvdrift0_new = np.round(cvdrift0_new, digits)
    cvdrift0_old = np.round(cvdrift0_old, digits)
    digits = 3
    cvdrift_new = B_times_kappa_dot_grady_new * 2
    cvdrift_new = np.round(cvdrift_new, digits)
    cvdrift_old = np.round(cvdrift_old, digits)
    gbdrift_new = B_times_gradB_dot_grady_new * 2
    gbdrift_new = np.round(gbdrift_new, digits)
    gbdrift_old = np.round(gbdrift_old, digits)
    
    # Compare values
    if not (np.allclose(rhotor_old, rhotor_new, equal_nan=True)): error = process_error('rhotor')
    if not (np.allclose(qinp_old, qinp_new, equal_nan=True)): error = process_error('qinp')
    if not (np.allclose(shat_old, shat_new, equal_nan=True)): error = process_error('shat')
    if not (np.allclose(aref_old, aref_new, equal_nan=True)): error = process_error('aref')
    if not (np.allclose(bref_old, bref_new, equal_nan=True)): error = process_error('bref')
    if not (np.allclose(z_scalefac_old, z_scalefac_new, equal_nan=True)): error = process_error('z_scalefac')
    
    # Compare arrays
    if not (np.allclose(alpha_old, alpha_new, equal_nan=True)): error = process_error('alpha')
    if not (np.allclose(zeta_old, zeta_new, equal_nan=True)): error = process_error('zeta')
    if not (np.allclose(bmag_old, bmag_new, equal_nan=True)): error = process_error('bmag')
    if not (np.allclose(b_dot_gradz_old, b_dot_gradz_new, equal_nan=True)): error = process_error('b_dot_gradz')
    if not (np.allclose(bdot_grad_z_old, bdot_grad_z_new, equal_nan=True)): error = process_error('b_dot_gradz')
    if not (np.allclose(grad_alpha2_old, grad_alpha2_new, equal_nan=True)): error = process_error('grad_alpha2')
    if not (np.allclose(gd_alph_psi_old, gd_alph_psi_new, equal_nan=True)): error = process_error('gd_alph_psi')
    if not (np.allclose(grad_psi2_old, grad_psi2_new, equal_nan=True)): error = process_error('grad_psi2')
    if not (np.allclose(gbdrift_old, gbdrift_new, equal_nan=True)): error = process_error('gbdrift')
    if not (np.allclose(gbdrift0_old, gbdrift0_new, equal_nan=True)): error = process_error('gbdrift0')
    if not (np.allclose(cvdrift_old, cvdrift_new, equal_nan=True)): error = process_error('cvdrift')
    if not (np.allclose(cvdrift0_old, cvdrift0_new, equal_nan=True)): error = process_error('cvdrift0')
    if not (np.allclose(theta_vmec_old, theta_vmec_new, equal_nan=True)): error = process_error('theta_vmec')
    if not (np.allclose(B_sub_theta_old, B_sub_theta_new, equal_nan=True)): error = process_error('B_sub_theta')
    if not (np.allclose(B_sub_zeta_old, B_sub_zeta_new, equal_nan=True)): error = process_error('B_sub_zeta')
    assert (not error), f'The geometry data does not match in the *.geometry file.'
    
    # Dont compare gds23 and gds24 since it was actually broken in old stella
    #if not (np.allclose(gds23_old, gds23_new, equal_nan=True)): error = process_error('gds23')
    #if not (np.allclose(gds24_old, gds24_new, equal_nan=True)): error = process_error('gds24')
    return
