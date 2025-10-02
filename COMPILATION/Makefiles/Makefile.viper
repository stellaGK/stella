# module load gcc/14 openmpi/4.1 netcdf-mpi/4.9.2 fftw-mpi/3.3.10 mkl/2025.2 hdf5-mpi/1.14.1

# Add to bash:
# module purge
# module load gcc/14 openmpi/4.1 netcdf-mpi/4.9.2 fftw-mpi/3.3.10 mkl/2025.2 hdf5-mpi/1.14.1
# export STELLA_SYSTEM=viper
# export MAKEFLAGS=-IMakefiles

COMPILER=gcc
USE_HDF5=on
include $(COMPILATION_DIR)/Makefiles/Compilers/Makefile.$(COMPILER)

USE_LOCAL_SPFUNC=on

# FFT libraries
ifeq ($(USE_FFT),fftw3)
        FFT_INC = -I$(FFTW3_ROOT)/include
        FFT_LIB = -L$(FFTW3_ROOT)/lib -lfftw3 -lfftw3f
endif
ifeq ($(USE_FFT),fftw)
        FFT_INC = -I$(FFTW3_ROOT)/include
        FFT_LIB = -L$(FFTW3_ROOT)/lib -lrfftw -lfftw
endif

# NetCDF with HDF5
ifdef USE_NETCDF
        NETCDF_INC = -I$(NETCDF_HOME)/include -I$(HDF5_HOME)/include
        NETCDF_LIB = -L$(NETCDF_HOME)/lib -lnetcdf -lnetcdff \
                     -L$(HDF5_HOME)/lib -lhdf5_hl -lhdf5
endif

ifdef USE_LAPACK
    ifeq ($(LAPACK_LIB),)
        LAPACK_LIB = \
            -L$(MKL_HOME)/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core \
            -L$(GCC_HOME)/lib64 -lgomp \
            -lm -ldl -lm -ldl
    endif
    CPPFLAGS += -DLAPACK
endif