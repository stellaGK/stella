#!/usr/bin/python3   
import os
import shutil
import pathlib
import difflib
import platform
import subprocess   
import configparser

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
def get_stella_path():
    '''Returns the absolute path to the stella executable.
    Can be controlled by setting the STELLA_EXE_PATH environment variable.'''
    default_path = get_stella_expected_run_directory() / '../../../stella'
    stella_path = pathlib.Path(os.environ.get('STELLA_EXE_PATH', default_path))
    return stella_path.absolute()

#-------------------------------------------------------------------------------
def run_stella(stella_path, input_file, nproc=None):
    '''Run stella with a given input file.''' 
    if not nproc: nproc = read_nproc()
    subprocess.run(['mpirun','--oversubscribe', '-np', f'{nproc}', stella_path, input_file], check=True)
    return

#-------------------------------------------------------------------------------
def copy_input_file(input_file: str, destination):
    '''Copy input_file to destination directory.'''
    shutil.copyfile(get_stella_expected_run_directory() / input_file, destination / input_file)
    return

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
def run_local_stella_simulation(input_file, tmp_path, vmec_file=None, nproc=None):
    ''' Run a local stella simulation in <tmp_path>. '''
    copy_input_file(input_file, tmp_path)
    copy_common_input_files(input_file, tmp_path)
    if vmec_file: copy_vmec_file(vmec_file, tmp_path)
    os.chdir(tmp_path); run_stella(get_stella_path(), input_file, nproc=nproc)
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
        
            # For the moments, there are noisy elements along (t,s,tube,zed,kx,ky,ri) so check the sum
            if key in ['density', 'upar', 'temperature', 'spitzer2']: 
                sum1 = np.sum(np.abs(local_quantity.data[:,:,:,:,:,:,:]), axis=(1,2,3,4,5,6))
                sum2 = np.sum(np.abs(expected_quantity.data[:,:,:,:,:,:,:]), axis=(1,2,3,4,5,6))
                if not (np.allclose(sum1, sum2, rtol=1e-8, atol=1e-100)):
                    print(f'\nERROR: The {key} arrays do not match in the netCDF files.'); error = True
                    print(f'Compare the {key} arrays in the local and expected netCDF files:')
                    compare_local_array_with_expected_array(local_quantity, expected_quantity, name=key1)  
        
            # On macos-14, the fluxes differ slightly more 
            elif key in ['qflux_vs_kxkyzs'] and (system=='Darwin') and ('23' in release):
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
        if not (np.allclose(local_time, expected_time, equal_nan=True)):
            print('\nERROR: The time axis does not match in the netCDF files.'); error = True
            print('\nCompare the time arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_time, expected_time)  
        if not (np.allclose(local_phi2, expected_phi2, equal_nan=True)):
            print('\nERROR: The potential data does not match in the netCDF files.'); error = True 
            print('\nCompare the potential arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_phi2, expected_phi2) 
        assert (not error), f'The potential data does not match in the netCDF files.' 
    
    return error
        
#-------------------------------------------------------------------------------  
def compare_local_potential_with_expected_potential_em(local_netcdf_file='', expected_netcdf_file='', run_data={}, error=False): 

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

        local_apar2 = local_netcdf['apar2']
        expected_apar2 = local_netcdf['apar2']

        local_bpar2 = local_netcdf['bpar2']
        expected_bpar2 = local_netcdf['bpar2']
        
        # Check whether we have the same time and potential data
        if not (np.allclose(local_time, expected_time, equal_nan=True)):
            print('\nERROR: The time axis does not match in the netCDF files.'); error = True
            print('\nCompare the time arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_time, expected_time)  
        if not (np.allclose(local_phi2, expected_phi2, equal_nan=True)):
            print('\nERROR: The <phi potential data does not match in the netCDF files.'); error = True 
            print('\nCompare the <phi> potential arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_phi2, expected_phi2)
        if not (np.allclose(local_apar2, expected_apar2, equal_nan=True)):
            print('\nERROR: The <A_parallel> potential data does not match in the netCDF files.'); error = True 
            print('\nCompare the <A_parallel> potential arrays in the local and expected netCDF files:')
            compare_local_array_with_expected_array(local_apar2, expected_apar2)
        if not (np.allclose(local_bpar2, expected_bpar2, equal_nan=True)):
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
                
                # The expected output has been made on mac-os where drhodpsi has 17 digits
                # And macOS-14 is only giving 14 digits for drhodpsi, ubuntu has 16 digits
                if key=='drhodpsi': 
                    local_netcdf[key].data = np.round(local_netcdf[key].data, 15)
                    expected_netcdf[key].data = np.round(expected_netcdf[key].data, 15)
                if (key=='drhodpsi') and (system=='Darwin') and ('23' in release): 
                    local_netcdf[key].data = np.round(local_netcdf[key].data, 13)
                    expected_netcdf[key].data = np.round(expected_netcdf[key].data, 13)
                    
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
                if not (np.allclose(local_netcdf[key], expected_netcdf[key], equal_nan=True)):
                    print(f'\nERROR: The quantity <{key}> does not match in the netcdf files.'); error = True
                    print(f'Compare the {key} arrays in the local and expected netCDF files:')
                    compare_local_array_with_expected_array(local_netcdf[key], expected_netcdf[key], name=key)  
                    
        # Print "AssertionError: <error message>" if an error was encountered
        assert (not error), f'Some geometry arrays in the netcdf file did not match the previous run.'  
        
    return error
