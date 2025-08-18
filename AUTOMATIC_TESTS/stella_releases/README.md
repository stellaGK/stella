

The stella releases are found at:
https://github.com/stellaGK/stella/releases

Note that for stella_v0.7 we commented line 224 in run_parameters.f90 (error = .true.), in order to allow stella to run when the variables `fapar` and `fbpar` are set in the input file. Since these variables are needed for stella_v0.6.
