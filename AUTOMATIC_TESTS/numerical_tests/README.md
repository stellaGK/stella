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

 On 09/10/2026 the default of flip_flop has been changed from .false. to .true., so that the order of
 the explicit and implicit operator splitting is alternated every time step (second order in time).
 The expected output of every test which treats a term implicitly and does not set flip_flop has been
 recreated: test 3 (except flip_flop.in and stream_and_mirror_explicit.in), the implicit terms and
 kxky grid options of test 4, the Miller and VMEC flux tubes of test 5, test 6 (now using the output
 of the current stella version), the FFS runs with all terms of test 7, and the kxky grid options and
 short and full KBM/TAE runs of test 8. The final KBM frequency changed by 0.7% and the TAE
 frequency by 8% (Re[omega] = 1.077 -> 1.162), which indicates that the TAE run is not converged in delt.
