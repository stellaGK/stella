# makefile for NERSC Edison and Cori
# You must first load the cray-netcdf module and python module:
#   module load cray-netcdf python
# It is convenient to run
#   module unload cray-libsci
# to avoid warning messages about libsci during compiling.

ifdef GK_SYSTEM
	HOSTNAME = $(GK_SYSTEM)
else
	ifdef NERSC_HOST
		HOSTNAME = $(NERSC_HOST)
	else
		HOSTNAME="laptop"
	endif
endif

ifeq ($(HOSTNAME),edison)
	export FC = ftn
	## NERSC documentation recommends against specifying -O3
	## -mkl MUST APPEAR AT THE END!!
	export EXTRA_COMPILE_FLAGS = -openmp -mkl
	export EXTRA_LINK_FLAGS =  -openmp -mkl -Wl,-ydgemm_
	# Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.

else ifeq ($(HOSTNAME),cori)
	export FC = ftn
	## NERSC documentation recommends against specifying -O3
	## -mkl MUST APPEAR AT THE END!!
	export EXTRA_COMPILE_FLAGS = -qopenmp -mkl
	export EXTRA_LINK_FLAGS =  -qopenmp -mkl -Wl,-ydgemm_
	# Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.
else ifeq ($(HOSTNAME),archer)
	export FC = ftn
	export EXTRA_COMPILE_FLAGS = -qopenmp
	export EXTRA_LINK_FLAGS = -qopenmp -Wl,-ydgemm_
else ifeq ($(HOSTNAME), marconi)
	export FC = ifort
	export EXTRA_COMPILE_FLAGS = -r8 -double-size 128 -qopenmp -I$(NETCDFF_HOME)/include/ -I$(NETCDF_INCLUDE)
	export EXTRA_LINK_FLAGS = -qopenmp -L$(NETCDFF_HOME)/lib/ -lnetcdff -L$(NETCDF_HOME)/lib/ -lnetcdf
else ifeq ($(HOSTNAME), dirac)
	export FC = ifort
	export EXTRA_COMPILE_FLAGS = -r8 -double-size 128 -qopenmp 
	export EXTRA_LINK_FLAGS = -qopenmp -lnetcdf -lnetcdff
else ifeq ($(HOSTNAME), xula)
	export FC =mpif90
	export EXTRA_COMPILE_FLAGS = -r8 -double-size 128 -qopenmp -I$(NETCDF_HOME)/include/
	export EXTRA_LINK_FLAGS = -qopenmp -L$(NETCDF_HOME)/lib/ -lnetcdff -lnetcdf
else
	FC = mpif90
#	EXTRA_COMPILE_FLAGS = -fopenmp -I/opt/local/include -ffree-line-length-none
	EXTRA_COMPILE_FLAGS = -fopenmp -I$(NETCDF_INC_DIR) -ffree-line-length-none
	EXTRA_COMPILE_FLAGS += -fdefault-real-8 -fdefault-double-8
#	EXTRA_LINK_FLAGS =  -fopenmp -L/opt/local/lib -lnetcdff  -lnetcdf
	EXTRA_LINK_FLAGS =  -fopenmp -L$(NETCDF_LIB_DIR) -lnetcdff  -lnetcdf
endif


# The variable LIBSTELL_DIR should either be "mini_libstell", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a directory containing libstell .mod files
# elsewhere on your system.
LIBSTELL_DIR = mini_libstell
#LIBSTELL_DIR=/Users/mattland/stellopt/LIBSTELL/Release

# The variable LIBSTELL_FOR_SFINCS should either be "mini_libstell/mini_libstell.a", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a libstell.a library elsewhere on your system.
LIBSTELL_FOR_STELLA = mini_libstell/mini_libstell.a
#LIBSTELL_FOR_STELLA=/Users/mattland/stellopt/LIBSTELL/Release/libstell.a

# End of system-dependent variable assignments

# Promote real to double, as stella does:
#EXTRA_COMPILE_FLAGS += -fdefault-real-8 -fdefault-double-8

export

.PHONY: all clean

all: test_vmec_to_stella_geometry_interface

include makefile.depend

%.o: %.f90
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

%.o: %.f
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

test_vmec_to_stella_geometry_interface: $(OBJ_FILES)
	$(FC) -o test_vmec_to_stella_geometry_interface $(OBJ_FILES) $(LIBSTELL_FOR_STELLA) $(EXTRA_LINK_FLAGS)

mini_libstell/mini_libstell.a:
	$(MAKE) -C mini_libstell

clean:
	rm -f *.o *.mod *.MOD *~ test_vmec_to_stella_geometry_interface
	cd mini_libstell; rm -f *.o *.mod *.MOD *.a

test_make:
	@echo HOSTNAME is $(HOSTNAME)
	@echo FC is $(FC)
	@echo EXTRA_COMPILE_FLAGS is $(EXTRA_COMPILE_FLAGS)
	@echo EXTRA_LINK_FLAGS is $(EXTRA_LINK_FLAGS)
