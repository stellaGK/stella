Automated numerical tests for stella
====================================

Perform some simulations and check whether the time trace of the potential matches exactly.
 
 
 
 Valid changes which required to update the tests
 ------------------------------------------------
 On 24/07/2024 ginit_option = 'default' has been changed to ensure that the reality condition
 is applied correctly. As a result, the following tests need to be updated, because they use this initialization:
    test_4_gyrokinetic_equation/test_3b_initalization_options.py::test_whether_init_default_option_is_the_same 
    test_4_gyrokinetic_equation/test_4d_parallel_streaming_implicit.py::test_whether_parallel_streaming_implicit_evolves_correctly
    test_4_gyrokinetic_equation/test_4e_parallel_streaming_explicit.py::test_whether_parallel_streaming_explicit_evolves_correctly
    test_4_gyrokinetic_equation/test_4f_mirror_implicit.py::test_whether_mirror_implicit_evolves_correctly
    test_4_gyrokinetic_equation/test_4g_mirror_explicit.py::test_whether_mirror_explicit_evolves_correctly  
    test_4_gyrokinetic_equation/test_4h_magnetic_drifts_explicit.py::test_whether_magnetic_drifts_explicit_evolves_correctly  
    test_4_gyrokinetic_equation/test_4i_diamagnetic_drift_explicit.py::test_whether_diamagnetic_drift_explicit_evolves_correctly  
    test_4_gyrokinetic_equation/test_4j_drifts_implicit.py::test_whether_drifts_implicit_evolves_correctly 
    test_5_fluxtube/test_5a_miller_linear_and_nonlinear.py::test_whether_miller_linear_evolves_correctly  

 On 14/11/2026 test 4f (test_4f_parallel_nonlinearity.py, parallel_nonlinearity.in) has been added, the first test of the parallel nonlinearity with rhostar > 0 (before, only test_all_flags2.in of test 5 switched it on, with rhostar = 0, where the term vanishes). It compares the full state and the parallel flow, relative to its own size, with the expected output, and checks that only the parallel nonlinearity drives a parallel flow from the initial noise. The two bugs of the parallel nonlinearity fixed in November 2026 change the parallel flow of this test by 1023% and 11%, and both make it fail.

 On 14/11/2026 the expected output of test 3 and of test 6a (miller_nonlinear .out.nc and .fluxes) has been recreated for the toroidal angular momentum flux of Miller geometry, which now uses R(z) instead of the constant major radius; only vflux changed.
