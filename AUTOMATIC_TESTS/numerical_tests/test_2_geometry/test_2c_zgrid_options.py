
################################################################################
#               Check all the terms in the gyrokinetic equation                #
################################################################################
# Test the options for the (kx,ky) grid which are 'box' and 'range'.
################################################################################

# Python modules
import os, sys
import pathlib 
import numpy as np
import xarray as xr  

# Package to run stella 
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

#-------------------------------------------------------------------------------
#                                BOX NONLINEAR                                 #
#-------------------------------------------------------------------------------
def test_z_grid_zedequalarc(tmp_path, error=False):
 
    #---------------------------------------------------------------------------
    #                       ZED_EQUAL_ARC = FALSE, MILLER                      #
    #---------------------------------------------------------------------------

    # Run stella and compare the geometry arrays in the netcdf file
    run_data = run_local_stella_simulation('zed_equal_arc_miller_true.in', tmp_path)  
    compare_geometry_in_netcdf_files(run_data, error=False) 
    print('  -->  The zed_equal_arc=True gives the correct Miller geometry.')
    
    #---------------------------------------------------------------------------
    #                       ZED_EQUAL_ARC = TRUE, MILLER                      #
    #---------------------------------------------------------------------------

    # Run stella and compare the geometry arrays in the netcdf file
    run_data = run_local_stella_simulation('zed_equal_arc_miller_false.in', tmp_path)  
    compare_geometry_in_netcdf_files(run_data, error=False) 
    print('  -->  The zed_equal_arc=False gives the correct Miller geometry.')
    return
 
    #---------------------------------------------------------------------------
    #                        ZED_EQUAL_ARC = FALSE, VMEC                       #
    #---------------------------------------------------------------------------

    # Run stella and compare the geometry arrays in the netcdf file
    run_data = run_local_stella_simulation('zed_equal_arc_vmec_true.in', tmp_path)  
    compare_geometry_in_netcdf_files(run_data, error=False) 
    print('  -->  The zed_equal_arc=True gives the correct VMEC geometry.')
    
    #---------------------------------------------------------------------------
    #                        ZED_EQUAL_ARC = TRUE, VMEC                       #
    #---------------------------------------------------------------------------

    # Run stella and compare the geometry arrays in the netcdf file
    run_data = run_local_stella_simulation('zed_equal_arc_vmec_false.in', tmp_path)  
    compare_geometry_in_netcdf_files(run_data, error=False) 
    print('  -->  The zed_equal_arc=False gives the correct VMEC geometry.')
    return
    
