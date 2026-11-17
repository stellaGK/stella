
################################################################################
#               Check all the terms in the gyrokinetic equation                #
################################################################################
# These test are designed to test the terms in the gyrokinetic equation one 
# by one. We test the following terms:
#     - nonlinear term                 (nonlinear)
#     - parallel streaming explicit    (include_parallel_streaming, stream_implicit)
#     - parallel streaming implicit    (include_parallel_streaming, stream_implicit)
#     - mirror explicit                (include_mirror, mirror_implicit)
#     - mirror implicit                (include_mirror, mirror_implicit)
#     - magnetic drifts explicit       (xdriftknob and ydriftknob)
#     - diamagnetic drifts explicit    (wstarknob)
#     - drifts implicit                (xdriftknob, ydriftknob, wstarknob, include_parallel_streaming, drifts_implicit)
# 
# The relevant parameters in the input file are:
#       &physics_flags 
#         nonlinear = .false.
#         include_parallel_streaming = .false.
#         include_mirror = .false.
#       /
#       &time_advance_knobs
#         xdriftknob = 0
#         ydriftknob = 0
#         wstarknob = 0 
#       /
#       &dissipation
#         hyper_dissipation = .false.
#       /  
#       &knobs 
#         stream_implicit = .false.
#         mirror_implicit = .false.
#         drifts_implicit = .false.
#       /
################################################################################

# Python modules
import pytest
import os, sys
import pathlib
import numpy as np
import xarray as xr

# Package to run stella 
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                                NONLINEAR TERM                                #
#-------------------------------------------------------------------------------
def test_whether_nonlinear_term_evolves_correctly(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('nonlinear_term.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files   
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The nonlinear term is evolving correctly.')
    return
    
#-------------------------------------------------------------------------------
#                          PARALLEL STREAMING IMPLICIT                         #
#-------------------------------------------------------------------------------
def test_whether_parallel_streaming_implicit_evolves_correctly(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('parallel_streaming_implicit.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The parallel streaming term, implemented implicitly, is evolving correctly.')
    return
    
#-------------------------------------------------------------------------------
#                          PARALLEL STREAMING EXPLICIT                         #
#-------------------------------------------------------------------------------
def test_whether_parallel_streaming_explicit_evolves_correctly(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('parallel_streaming_explicit.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The parallel streaming term, implemented explicitly, is evolving correctly.')
    return
    
#-------------------------------------------------------------------------------
#                               MIRROR IMPLICIT                                #
#-------------------------------------------------------------------------------
def test_whether_mirror_implicit_evolves_correctly(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('mirror_implicit.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The mirror term, implemented implicitly, is evolving correctly.')
    return
    
#-------------------------------------------------------------------------------
#                               MIRROR EXPLICIT                                #
#-------------------------------------------------------------------------------
def test_whether_mirror_explicit_evolves_correctly(tmp_path, stella_version, error=False): 

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('mirror_explicit.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The mirror term, implemented explicitly, is evolving correctly.')
    return
    
#-------------------------------------------------------------------------------
#                           MAGNETIC DRIFTS EXPLICIT                           #
#-------------------------------------------------------------------------------
def test_whether_magnetic_drifts_explicit_evolves_correctly(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('magnetic_drifts_explicit.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The magnetic drifts, implemented explicitly, are evolving correctly.')
    return
    
#-------------------------------------------------------------------------------
#                          DIAMAGNETIC DRIFT EXPLICIT                          #
#-------------------------------------------------------------------------------
def test_whether_diamagnetic_drift_explicit_evolves_correctly(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('diamagnetic_drift_explicit.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The diamagnetic drift, implemented explicitly, is evolving correctly.')
    return
    
#-------------------------------------------------------------------------------
#                               DRIFTS IMPLICIT                                #
#-------------------------------------------------------------------------------
def test_whether_drifts_implicit_evolves_correctly(tmp_path, stella_version, error=False):

    # Run stella inside of <tmp_path>
    run_data = run_local_stella_simulation('drifts_implicit.in', tmp_path, stella_version)
     
    # Compare phi2(t) in the netCDF files
    compare_local_potential_with_expected_potential(run_data=run_data, error=False)
    print('  -->  The drifts, implemented explicitly, are evolving correctly.')
    print('  -->  WARNING: the actual implicit drifts are broken and should not be used.')
    print('  -->  We are just testing here whether the routine has changed.')
    return

