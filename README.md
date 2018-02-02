# stella

INSTALLATION
1) set GK_SYSTEM='system', with system replaced by the appropriate system on which you are running.
2) make sure that fftw3 and netcdf libraries are installed.
3) set the environmental variables FFTW_LIB_DIR (directory containing libfftw3.a),
FFTW_INC_DIR (directory including fftw3.f), NETCDF_LIB_DIR (directory containing libnetcdff.a),
and NETCDF_INC_DIR (directory including netcdf.inc)
4) set the environmental variable MAKEFLAGS=-IMakefiles
