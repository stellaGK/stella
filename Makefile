####################################################################
#                                                                  #
#       Makefile for the stella gyrokinetic turbulence code        # 
#                                                                  #
####################################################################
#
# Makefile written by Bill Dorland and Ryusuke Numata
#
# In the ~/.source.sh or ~/.bashrc file of your computer define, e.g., 
# export GK_SYSTEM='marconi', which matches stella/Makefiles/Makefile.$GK_SYSTEM
#
# * Available Compilers (tested on limited hosts)
#   (must be Fortran 95 Standard compliant)
#
# Intel ifort
# GNU's gfortran and g95
# IBM XL Fortran xlf90
# PathScale Compiler Suite pathf90
# The Portland Group pgf90
# NAGWare f95 (v5.1)
# Lahey/Fujitsu Fortran lf95
# 
# * Frequently Tested Hosts, Systems
#
# Standard Linux
# Standard Mac OS X with MacPorts
# Franklin at NERSC and Jaguar at NCCS (Cray XT4 with PGI)
# Bassi at NERSC (IBM Power 5 with IBM XL Fortran)
# Ranger (Sun Constellation Linux Cluster with Intel)
####################################################################


####################################################################
#                             Switches                             #
####################################################################
# In the brackets, the accepted values are shown, where "undefined" means blank.
# Note that switches with (bin) only check whether they are defined or not,
# hence the value you give it does not matter, e.g., DEBUG=off means DEBUG=on.
# The "?=" sets a variable if it does not already have a value. 
####################################################################

DEBUG ?=                   # Turns on debug mode (bin) 
TEST ?=                    # Turns on test mode (bin) 
OPT ?= on                  # Optimization (on, aggressive, undefined)
STATIC ?=                  # Prevents linking with shared libraries (bin)
DBLE ?= on                 # Promotes precisions of real and complex (bin) 
USE_NETCDF ?= on           # Uses netcdf library (bin)
USE_PARALLEL_NETCDF ?=     # Uses parallel netcdf library
USE_HDF5 ?=                # Uses hdf5 library (bin)
USE_LOCAL_RAN ?=           # Use local random number generator (mt,undefined) (see also README)
USE_LOCAL_SPFUNC ?=        # Use local special functions (bin)
USE_NAGLIB ?=              # Use nag libraray (spfunc,undefined)
USE_SFINCS ?=              # Link to sfincs library at compilation
USE_LAPACK ?= on           # Use LAPACK, needed for test particle collisions
HAS_ISO_C_BINDING ?= on    # Does the compiler support the iso_c_binding features of Fortran? (needed for local parallel LU decomposition) 

# Which FFT library to use (fftw,fftw3,mkl_fftw,undefined)
# We can't have any spaces behind 'fftw3' so we can't put the comment behind it
USE_FFT ?= fftw3

# The following variables can be set in platform-dependent Makefile's
# For example, inside stella/Makefiles/Makefile.marconi

MAKE		 = make
CPP		 = cpp
CPPFLAGS	 = -P -traditional
export FC	 = f90
export MPIFC	?= mpif90
H5FC		?= h5fc
H5FC_par	?= h5pfc
F90FLAGS	 =
F90OPTFLAGS	 =
CC		 = cc
MPICC		?= mpicc
H5CC		?= h5cc
H5CC_par	?= h5pcc
CFLAGS		 =
COPTFLAGS 	 =
LD 		 = $(FC)
LDFLAGS 	 = $(F90FLAGS)
ARCH 		 = ar
ARCHFLAGS 	 = cr
RANLIB		 = ranlib
AWK 		 = awk
PERL		 = perl
FORD       	?= ford

MPI_INC ?=
MPI_LIB ?=
FFT_INC ?=
FFT_LIB ?=
NETCDF_INC ?=
NETCDF_LIB ?=
HDF5_INC ?=
HDF5_LIB ?=
NAG_LIB ?=
NAG_PREC ?= dble
SFINCS_LIB ?=
SFINCS_INC ?=
PETSC_LIB ?=
PETSC_INC ?=
LIBSTELL_LIB ?=

# Record the top level path. Note we don't just use $(PWD) as this
# resolves to the directory from which make was invoked. The approach
# taken here ensures that GK_HEAD_DIR is the location of this
# Makefile. Note the realpath call removes the trailing slash so
# later we need to add a slash if we want to address subdirectories
GK_THIS_MAKEFILE := $(abspath $(lastword $(MAKEFILE_LIST)))
GK_HEAD_DIR := $(realpath $(dir $(GK_THIS_MAKEFILE)))

####################################################################
#                        PLATFORM DEPENDENCE                       #
####################################################################
# The compiler mode switches (DEBUG, TEST, OPT, STATIC, DBLE) must be set
# before loading Makefile.$(GK_SYSTEM) because they may affect compiler options.
# However, Makefile.local may override some options set in Makefile.$(GK_SYSTEM),
# thus it is included before and after Makefile.$(GK_SYSTEM)
####################################################################

# Here "sinclude" will include the file, without an error if it doesn't exist
sinclude Makefile.local

# Include system-dependent make variables, defined in stella/Makefiles/Makefile.$GK_SYSTEM 
# where e.g., export GK_SYSTEM='marconi' inside e.g. ~/.bashrc defined this system variable on your computer 
ifndef GK_SYSTEM
	ifdef SYSTEM
$(warning SYSTEM environment variable is obsolete)
$(warning use GK_SYSTEM instead)
	GK_SYSTEM = $(SYSTEM)
	else
$(error GK_SYSTEM is not set)
	endif
endif

# Read the stella/Makefiles/Makefile.$GK_SYSTEM file
include Makefile.$(GK_SYSTEM)

# Overwrite the file with the local file, if it exists
sinclude Makefile.local

####################################################################
#                        PLATFORM DEPENDENCE                       #
####################################################################

# Export some variables, so they become available to all child processes
export F90FLAGS
export NETCDF_INC
export NETCDF_LIB

# Export the subdirectories
export DIAG=$(GK_HEAD_DIR)/diagnostics
export COLL=$(GK_HEAD_DIR)/dissipation
export UTILS=$(GK_HEAD_DIR)/utils
export GEO=$(GK_HEAD_DIR)/geometry
export LIBSTELL=$(UTILS)/mini_libstell

# We make extra libraries <mini_libstell> and <stella_utils>
LIBSTELL_LIB=$(LIBSTELL)/mini_libstell.a

# External libraries
GIT_VERSION_DIR := $(GK_HEAD_DIR)/externals/git_version
NEASYF := $(GK_HEAD_DIR)/externals/neasyf/src

####################################################################
#                        PROCESS SOME FLAGS                        #
####################################################################

# Must invoke full functionality when we do >> make depend
ifeq ($(MAKECMDGOALS),depend)
	MAKE += USE_HDF5=on USE_FFT=fftw3 USE_NETCDF=on \
		USE_LOCAL_BESSEL=on USE_LOCAL_RAN=mt
endif

# Warning that we can't do nonlinear sims without fourier transformations
ifndef USE_FFT
$(warning USE_FFT is off)
$(warning Be sure that nonlinear run makes no sense)
endif

ifdef HAS_ISO_C_BINDING
	CPPFLAGS += -DISO_C_BINDING
endif

FC = $(MPIFC)
CC = $(MPICC)
CPPFLAGS += -DMPI

# Flags for <USE_FFT> are (fftw,fftw3,mkl_fftw,undefined)
ifeq ($(USE_FFT),fftw)
	CPPFLAGS += -DFFT=_FFTW_
	ifeq ($(FFT_LIB),)
		FFT_LIB = -lfftw -lrfftw
	endif
endif
ifeq ($(USE_FFT),fftw3)
	CPPFLAGS += -DFFT=_FFTW3_
	ifeq ($(FFT_LIB),)
		FFT_LIB = -lfftw3
	endif
endif
ifeq ($(USE_FFT),mkl_fftw)
	CPPFLAGS += -DFFT=_FFTW_
endif

ifdef USE_NETCDF
	ifeq ($(NETCDF_LIB),)
		NETCDF_LIB = -lnetcdf
        endif
	CPPFLAGS += -DNETCDF
endif

ifdef USE_LAPACK
        LAPACK_LIB ?= -llapack
	CPPFLAGS += -DLAPACK
endif

ifdef USE_HDF5
	FC = $(H5FC_par)
	CC = $(H5CC_par)
	ifdef USE_PARALLEL_NETCDF
		CPPFLAGS += -DNETCDF_PARALLEL
	endif
	CPPFLAGS += -DHDF
endif

ifeq ($(USE_LOCAL_RAN),mt)
	CPPFLAGS += -DRANDOM=_RANMT_
endif
ifdef USE_LOCAL_SPFUNC
	CPPFLAGS += -DSPFUNC=_SPLOCAL_
else
	ifeq ($(findstring spfunc,$(USE_NAGLIB)),spfunc)
		CPPFLAGS += -DSPFUNC=_SPNAG_
	endif
endif

ifdef USE_NAGLIB
	ifeq ($(NAG_PREC),dble)
		ifndef DBLE
$(warning Precision mismatch with NAG libarray)	
		endif
		CPPFLAGS += -DNAG_PREC=_NAGDBLE_
	endif
	ifeq ($(NAG_PREC),sngl)
		ifdef DBLE
$(warning Precision mismatch with NAG libarray)	
		endif
		CPPFLAGS += -DNAG_PREC=_NAGSNGL_
	endif
endif

ifdef USE_SFINCS
	CPPFLAGS += -DUSE_SFINCS
endif

# List of all the used libraries, where $LIBSTELL_LIB is our own library
LIBS	+= $(DEFAULT_LIB) $(MPI_LIB) $(FFT_LIB) $(NETCDF_LIB) $(HDF5_LIB) \
		$(NAG_LIB) $(SFINCS_LIB) $(PETSC_LIB) $(LIBSTELL_LIB) $(LAPACK_LIB)

# Include flags
INC_FLAGS= $(DEFAULT_INC) $(MPI_INC) $(FFT_INC) $(NETCDF_INC) $(HDF5_INC) \
	   $(SFINCS_INC) $(PETSC_INC) $(LAPACK_INC)

# Flags for the fortran compiler and the C compiler
F90FLAGS+= $(F90OPTFLAGS)
CFLAGS += $(COPTFLAGS)

####################################################################
#                        PROCESS DIRECTORIES                       #
####################################################################

DATE=$(shell date +%y%m%d)
TARDIR=stella_$(DATE)

# If we called make in a subdirectory, get the directory above it
TOPDIR=$(CURDIR)
ifeq ($(notdir $(CURDIR)), $(COLL))
	TOPDIR=$(subst /$(COLL),,$(CURDIR))
endif
ifeq ($(notdir $(CURDIR)), $(DIAG))
	TOPDIR=$(subst /$(DIAG),,$(CURDIR))
endif
ifeq ($(notdir $(CURDIR)), $(UTILS))
	TOPDIR=$(subst /$(UTILS),,$(CURDIR))
endif
ifeq ($(notdir $(CURDIR)), $(GEO))
	TOPDIR=$(subst /$(GEO),,$(CURDIR))
endif
ifeq ($(notdir $(CURDIR)), $(LIBSTELL))
	TOPDIR=$(subst /$(LIBSTELL),,$(CURDIR))
endif
ifneq ($(TOPDIR),$(CURDIR))
	SUBDIR=true
endif

# VPATH is a list of directories to be searched for missing source files
VPATH = $(DIAG):$(COLL):$(UTILS):$(GEO):$(LIBSTELL):$(GIT_VERSION_DIR)/src:$(NEASYF)

# Removes non-existing directories from VPATH
VPATH_tmp := $(foreach tmpvp,$(subst :, ,$(VPATH)),$(shell [ -d $(tmpvp) ] && echo $(tmpvp)))
VPATH = .:$(shell echo $(VPATH_tmp) | sed "s/ /:/g")

# If we're in a subdirectory, add the parent directory 
ifdef SUBDIR
	VPATH +=:..
endif

####################################################################
#            READ THE DEPENDENCIES OF THE SOURCE FILES             #
####################################################################

# The dependencies of one source file on other source files
# is given in the file stella/Makefile.depend
DEPEND=Makefile.depend
DEPEND_CMD=$(PERL) fortdep

# Most common include and library directories
DEFAULT_INC_LIST = . $(DIAG) $(COLL) $(UTILS) $(LIBSTELL) $(GEO) $(NEASYF)
DEFAULT_LIB_LIST =
DEFAULT_INC=$(foreach tmpinc,$(DEFAULT_INC_LIST),$(shell [ -d $(tmpinc) ] && echo -I$(tmpinc)))
DEFAULT_LIB=$(foreach tmplib,$(DEFAULT_LIB_LIST),$(shell [ -d $(tmplib) ] && echo -L$(tmplib)))

# List of intermediate f90 files generated by preprocessor
F90FROMFPP = $(patsubst %.fpp,%.f90,$(notdir $(wildcard *.fpp */*.fpp)))

# Files that we know are going to be in the submodule directories.
# This is in case the directories exist but the submodules haven't actually been initialised
GIT_VERSION_SENTINEL := $(GIT_VERSION_DIR)/Makefile
NEASYF_SENTINEL := $(NEASYF)/neasyf.f90

# If make is killed or interrupted during the execution of their recipes, the targets in ".PRECIOUS" are not deleted.
.PRECIOUS: $(GIT_VERSION_SENTINEL) $(NEASYF_SENTINEL)

# Make sure the submodules exist; they all currently share the same recipe
$(GIT_VERSION_SENTINEL) $(NEASYF_SENTINEL) submodules:
	@echo "Downloading submodules"
	git submodule update --init --recursive

# Dump the compilation flags to a file, so we can check if they change between
# invocations of `make`. The `cmp` bit checks if the file contents
# change. Adding a dependency of a file on `.compiler_flags` causes it to be
# rebuilt when the flags change. Taken from
# https://stackoverflow.com/a/3237349/2043465
COMPILER_FLAGS_CONTENTS = "FC = $(FC)\n CPPFLAGS = $(CPPFLAGS)\n F90FLAGS = $(F90FLAGS)\n INC_FLAGS = $(INC_FLAGS)\n CFLAGS = $(CFLAGS)"
COMPILER_FLAGS_CONTENTS += "\n FORTRAN_GIT_DEFS = $(FORTRAN_GIT_DEFS)"
.PHONY: force
.compiler_flags: force
	@echo -e $(COMPILER_FLAGS_CONTENTS) | cmp -s - $@ || echo -e $(COMPILER_FLAGS_CONTENTS) > $@

# Use the make rules from fortran-git-version, which requires ensuring the submodules exist
include $(GIT_VERSION_DIR)/Makefile

# We have some extra macros to add when preprocessing this file. We'll
# need the submodule to exist first
git_version_impl.f90: $(GIT_VERSION_DIR)/src/git_version_impl.F90 | $(GIT_VERSION_SENTINEL)
	$(CPP) $(CPPFLAGS) $(FORTRAN_GIT_DEFS) $< $@

# The fortdep script doesn't know about Fortran submodules, so we need
# to write the dependency ourselves
git_version_impl.o: git_version_impl.f90 git_version.o .compiler_flags
$(GIT_VERSION_DIR)/src/git_version_impl.F90: .compiler_flags


####################################################################
#                              RULES                               #
#################################################################### 

# We process the following files
.SUFFIXES: .fpp .f90 .F90 .c .o

# Use Fortran f90 compiler for ".f90.o" files, e.g., ifort = The Intel Fortran compiler
.f90.o:
	$(FC) $(F90FLAGS) $(INC_FLAGS) -c $<

# Use C preprocessor (CPP) compiler for ".fpp.f90" files
.fpp.f90:
	$(CPP) $(CPPFLAGS) $< $@

# Use Fortran f90 compiler for ".F90.o" files, but also use the C preprocessor
.F90.o:
	$(FC) $(F90FLAGS) $(CPPFLAGS) $(INC_FLAGS) -c $<

# Use a C compiler for ".c.o" files
.c.o:
	$(CC) $(CFLAGS) -c $<

####################################################################
#                              Targets                             #
####################################################################
# Targets represent executables, libraries, and utilities built by CMake
#     depend: generate dependency
#     clean: clean up
#     distclean: does "make clean" + removes platform links & executables
# "stella_all" is defined in stella/Makefile.target_stella
####################################################################

# .DEFAULT_GOAL works for GNU make 3.81 (or higher).  For 3.80 or less, see all target
# If we're in the top folder, our target is <stella_all>, otherwise it is <utils_all> or <mini_libstell_all>
.DEFAULT_GOAL := stella_all
ifeq ($(notdir $(CURDIR)),utils) 
	.DEFAULT_GOAL := utils_all
endif
ifeq ($(notdir $(CURDIR)),mini_libstell)
	.DEFAULT_GOAL := mini_libstell_all
endif

.PHONY: all stella_all

all: $(.DEFAULT_GOAL)

include $(DEPEND)

sinclude Makefile.target_stella

# Include the automated Fortran tests
# Inside 'tests/automated_fortran_tests/Makefile' we define the commands 
# >> build-pfunit-library 
# >> run-automated-fortran-tests
tests/automated_fortran_tests/Makefile:
include tests/automated_fortran_tests/Makefile

tests/automated_numerical_tests_for_stella/Makefile:
include tests/automated_numerical_tests_for_stella/Makefile
 
# Run all tests together with the 'check' command
check: run-automated-fortran-tests run-automated-numerical-tests-for-stella

############################################################### SPECIAL RULES

# comment this out to keep intermediate .f90 files
#.PRECIOUS: $(F90FROMFPP)

.INTERMEDIATE: stella_transforms.f90 stella_io.f90 stella_save.f90 \
		mp.f90 fft_work.f90 response_matrix.f90 sources.f90 \
		fields.f90 mp_lu_decomposition.f90 git_version_impl.f90

############################################################# MORE DIRECTIVES

.PHONY: depend clean distclean

depend: $(GIT_VERSION_SENTINEL) $(NEASYF)
	@$(DEPEND_CMD) -m "$(MAKE)" -1 -o -v=0 $(VPATH)


# To compile correctly it is very important to remove the previously build
# *.o and *.mod files. Here '$(MAKE) -C $(LIBSTELL) clean' will look at 
# utils/mini_libstell/makefile and it will its clean target as well.
clean:
	-rm -f *.o *.mod *.g90 *.h core */core *~ *.smod
	-rm -f $(GEO)/*.o $(GEO)/*~
	-rm -f Makefiles/*~
	-rm -f $(UTILS)/*.o $(UTILS)/*~
	-rm -f $(COLL)/*.o $(COLL)/*~
	-rm -f $(DIAG)/*.o $(DIAG)/*~
	-rm -f .compiler_flags
	$(MAKE) -C $(LIBSTELL) clean

cleanlib:
	-rm -f *.a

distclean: unlink clean cleanlib

%.o: %.mod

unlink:
	-rm -f $(F90FROMFPP)

revision:
	@LANG=C svn info | awk '{if($$1=="Revision:") printf("%20d",$$2) }' > Revision

TAGS:	*.f90 *.fpp */*.f90 */*.fpp
	etags $^

############################################################# Documentation

create_namelist_markdown:
	docs/pages/user_manual/namelist_files/combine_namelists.sh

ifneq ("$(wildcard $(shell which $(FORD) 2>/dev/null))","")
check_ford_install:
	@echo "Using ford at $(shell which $(FORD))"
else
check_ford_install:
	@echo "Ford command $(FORD) not in path -- is it installed?\\n\\tConsider installing with 'pip install --user ford' and add ${HOME}/.local/bin to PATH" ; which $(FORD)
endif

doc: docs/stella_docs.md create_namelist_markdown check_ford_install
	$(FORD) $(INC_FLAGS) -r $(GIT_VERSION) $<

cleandoc:
	@echo "FORD docs"
	-rm -rf docs/html
