#################################################################### OVERVIEW
#
#  Makefile for Trinity (GS2) / AstroGK Gyrokinetic Turbulence code 
#               GSP Gyrokinetic PIC code
#               rmhdper reduced MHD code
#
#  (requires GNU's gmake)
#
#GK_PROJECT ?= trinity
GK_PROJECT ?= gs2
#GK_PROJECT ?= agk
#GK_PROJECT ?= gsp
#GK_PROJECT ?= rmhdper
#
#  Makefile written by Bill Dorland and Ryusuke Numata
#
#  LAST UPDATE: 10/12/09
#
# * Changelogs
# 26/09/10: added targets doc_depend and doc to make doxygen documentation. EGH.
#	12/17/09: make it explicit that the codes are written in
#		  Fortran 95 Standard
#	10/12/09: share the Makefile with gsp gyrokinetic PIC code
#	04/11/09: * share the Makefile with rmhdper reduced MHD code
#		  * include progject specific target definitions
#		    Makefile.target_$(GK_PROJECT)
#	04/06/09: SYSTEM environment variable is replaced by GK_SYSTEM
#	02/25/09: gs2 replaced by trinity (MAB)
#	01/26/09: USE_C_INDEX is imported to gs2 by TT
#       12/11/08: add support for NAGWare and Lahey compilers
#       10/28/08: some non-standard macros respect environment variables
#       10/27/08: default commands for MPI and HDF5 are defined
#       09/26/08: hack for treating intermediate f90 files properly
#	 	 (this is to keep all intermediate f90 files when 
#		 .PRECIOUS is given)
#	07/01/08: new switches USE_LOCAL_SPFUNC and USE_NAGLIB
#	06/16/08: clean up unused statements related to PLATFORM_LINKS
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
#
# * Switches:
#
# Here is a list of switches with simple explanation.
# In the brackets, values accepted are shown,
# where "undefined" means blank.
# Switches with (bin) are checked if they are defined or not
# What values they have do not matter.
# Be careful that DEBUG=off means DEBUG=on.
#
# turns on debug mode (bin)
DEBUG ?=
# turns on scalasca instrumentation mode (bin)
SCAL ?=
# turns on test mode (bin)
TEST ?=
# turns on profile mode (gprof,ipm)
PROF ?=
# optimization (on,aggressive,undefined)
OPT ?= on
# prevents linking with shared libraries (bin)
STATIC ?=
# promotes precisions of real and complex (bin)
DBLE ?= on
# turns on distributed memory parallelization using MPI (bin)
USE_MPI ?= on
# turns on SHMEM parallel communications on SGI (bin)
USE_SHMEM ?=
# which FFT library to use (fftw,fftw3,mkl_fftw,undefined) 
USE_FFT ?= fftw3
# uses netcdf library (bin)
USE_NETCDF ?= on
# uses parallel netcdf library
USE_PARALLEL_NETCDF ?=
# uses hdf5 library (bin)
USE_HDF5 ?= 
# Use function pointer in layouts_indices.c (bin)
# see also README
USE_C_INDEX ?= 
# Use local random number generator (mt,undefined)
# see also README
USE_LOCAL_RAN ?=
# Use posix for command_line (bin)
USE_POSIX ?=
# Use local special functions (bin)
USE_LOCAL_SPFUNC ?= 
# Use nag libraray (spfunc,undefined)
USE_NAGLIB ?= 
# Make GS2 into a library which can be called by external programs
MAKE_LIB ?=
# Include higher-order terms in GK equation arising from low-flow physics
LOWFLOW ?=
# Use le_layout for collision operator
USE_LE_LAYOUT ?=
#
# * Targets:
#
#  depend: generate dependency
#  test_make: print values of variables
#  clean: clean up
#  distclean: does "make clean" + removes platform links & executables
#  tar: pack
#
############################################################### DEFAULT VALUE
#
# These variables can be set in platform-dependent Makefile.
#

MAKE		= make
CPP		= cpp
CPPFLAGS	= -C -P -traditional
FC		= f90
MPIFC		?= mpif90
H5FC		?= h5fc
H5FC_par	?= h5pfc
F90FLAGS	= -framework Accelerate
F90OPTFLAGS	=
CC		= cc
MPICC		?= mpicc
H5CC		?= h5cc
H5CC_par	?= h5pcc
CFLAGS		=
COPTFLAGS 	=
LD 		= $(FC)
LDFLAGS 	= $(F90FLAGS)
ARCH 		= ar
ARCHFLAGS 	= cr
RANLIB		= ranlib
AWK 		= awk
PERL		= perl

# These macros are used for the suffix problem of absoft
F90FLAGS_SFX0 =
F90FLAGS_SFX1 =
F90FLAGS_SFX2 =

MPI_INC	?=
MPI_LIB ?=
FFT_INC ?=
FFT_LIB ?=
NETCDF_INC ?=
NETCDF_LIB ?=
HDF5_INC ?=
HDF5_LIB ?=
IPM_LIB ?=
NAG_LIB ?=
NAG_PREC ?= dble
PGPLOT_LIB ?=

################################################### SET COMPILE MODE SWITCHES

ifdef TEST
$(warning TEST mode is not working yet)
	override TEST =
endif

######################################################### PLATFORM DEPENDENCE

# compile mode switches (DEBUG, TEST, PROF, OPT, STATIC, DBLE)
# must be set before loading Makefile.$(GK_SYSTEM) because they may affect
# compiler options.
# However, Makefile.local may override some options set in Makefile.$(GK_SYSTEM),
# thus it is included before and after Makefile.$(GK_SYSTEM)
sinclude Makefile.local

# include system-dependent make variables
ifndef GK_SYSTEM
	ifdef SYSTEM
$(warning SYSTEM environment variable is obsolete)
$(warning use GK_SYSTEM instead)
	GK_SYSTEM = $(SYSTEM)
	else
$(error GK_SYSTEM is not set)
	endif
endif
include Makefile.$(GK_SYSTEM)

# include Makefile.local if exists
sinclude Makefile.local

#############################################################################

UTILS=utils
GEO=geo

ifeq ($(GK_PROJECT),rmhdper)
	override USE_MPI =
	override USE_FFT = fftw
	override USE_HDF5 =
endif

ifeq ($(MAKECMDGOALS),depend)
# must invoke full functionality when make depend
	MAKE += USE_HDF5=on USE_FFT=fftw USE_NETCDF=on USE_MPI=on \
		USE_LOCAL_BESSEL=on USE_LOCAL_RAN=mt
endif

ifdef USE_SHMEM
$(warning USE_SHMEM is not working yet)
	override USE_SHMEM =	
endif

ifdef USE_HDF5
	ifndef USE_MPI
$(error Currently, USE_HDF5 works with USE_MPI)
	endif
endif

ifndef USE_FFT
$(warning USE_FFT is off)
$(warning Be sure that nonlinear run makes no sense)
endif

ifdef USE_MPI
	FC = $(MPIFC)
	CC = $(MPICC)
	CPPFLAGS += -DMPI
endif
ifdef USE_SHMEM
	CPPFLAGS += -DSHMEM
endif
ifeq ($(USE_FFT),fftw)
	CPPFLAGS += -DFFT=_FFTW_
	FFT_LIB ?= -lfftw -lrfftw
endif

ifeq ($(USE_FFT),fftw3)
	CPPFLAGS += -DFFT=_FFTW3_
	FFT_LIB ?= -lfftw -lrfftw
endif

ifeq ($(USE_FFT),mkl_fftw)
	CPPFLAGS += -DFFT=_FFTW_
endif

ifdef USE_NETCDF
	NETCDF_LIB ?= -lnetcdf
	CPPFLAGS += -DNETCDF
endif
ifdef USE_HDF5
	ifdef USE_MPI
		FC = $(H5FC_par)
		CC = $(H5CC_par)
		ifdef USE_PARALLEL_NETCDF
			CPPFLAGS += -DNETCDF_PARALLEL
		endif

	else
		FC = $(H5FC)
		CC = $(H5CC)
	endif
	CPPFLAGS += -DHDF
endif
ifdef USE_C_INDEX
	CPPFLAGS += -DUSE_C_INDEX
endif
ifeq ($(USE_LOCAL_RAN),mt)
	CPPFLAGS += -DRANDOM=_RANMT_
endif
ifdef USE_POSIX
	CPPFLAGS += -DPOSIX
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
ifndef PGPLOT_LIB
	ifeq ($(MAKECMDGOALS),agk_fields_plot)
$(error PGPLOT_LIB is not defined)
	endif
endif
ifdef MAKE_LIB
	CPPFLAGS += -DMAKE_LIB
endif
ifdef LOWFLOW
	CPPFLAGS += -DLOWFLOW
endif
ifdef USE_LE_LAYOUT
	CPPFLAGS += -DUSE_LE_LAYOUT
endif

LIBS	+= $(DEFAULT_LIB) $(MPI_LIB) $(FFT_LIB) $(NETCDF_LIB) $(HDF5_LIB) \
		$(IPM_LIB) $(NAG_LIB)
PLIBS 	+= $(LIBS) $(PGPLOT_LIB)
F90FLAGS+= $(F90OPTFLAGS) \
	   $(DEFAULT_INC) $(MPI_INC) $(FFT_INC) $(NETCDF_INC) $(HDF5_INC)
CFLAGS += $(COPTFLAGS)

DATE=$(shell date +%y%m%d)
TARDIR=$(GK_PROJECT)_$(DATE)
TOPDIR=$(CURDIR)
ifeq ($(notdir $(CURDIR)), $(UTILS))
	TOPDIR=$(subst /$(UTILS),,$(CURDIR))
endif
ifeq ($(notdir $(CURDIR)), $(GEO))
	TOPDIR=$(subst /$(GEO),,$(CURDIR))
endif
ifneq ($(TOPDIR),$(CURDIR))
	SUBDIR=true
endif

VPATH = $(UTILS):$(GEO):Aux:../$(UTILS):../$(GEO)
# this just removes non-existing directory from VPATH
VPATH_tmp := $(foreach tmpvp,$(subst :, ,$(VPATH)),$(shell [ -d $(tmpvp) ] && echo $(tmpvp)))
VPATH = .:$(shell echo $(VPATH_tmp) | sed "s/ /:/g")
#
ifdef SUBDIR
	VPATH +=:..
endif
DEPEND=Makefile.depend
DEPEND_CMD=$(PERL) fortdep

# most common include and library directories
DEFAULT_INC_LIST = . $(UTILS) $(GEO) .. ../$(UTILS) ../$(GEO)
#DEFAULT_INC_LIST = . $(UTILS) $(GEO) .. ../$(UTILS) ../$(GEO) \
#		/usr/include /usr/local/include \
#	   	/opt/local/include /sw/include
DEFAULT_LIB_LIST =
#DEFAULT_LIB_LIST = /usr/lib /usr/local/lib \
#		/opt/local/lib /sw/lib
# This default library path list can simplify the procedure of porting,
# however, I found this (actually -L/usr/lib flag) causes an error
# when linking gs2 at bassi (RS6000 with xl fortran)
DEFAULT_INC=$(foreach tmpinc,$(DEFAULT_INC_LIST),$(shell [ -d $(tmpinc) ] && echo -I$(tmpinc)))
DEFAULT_LIB=$(foreach tmplib,$(DEFAULT_LIB_LIST),$(shell [ -d $(tmplib) ] && echo -L$(tmplib)))

# list of intermediate f90 files generated by preprocessor
F90FROMFPP = $(patsubst %.fpp,%.f90,$(notdir $(wildcard *.fpp */*.fpp)))
ifdef SCAL
   FC:= scalasca -instrument $(FC)
endif
####################################################################### RULES

.SUFFIXES:
.SUFFIXES: .fpp .f90 .c .o

.f90.o: 
	$(FC) $(F90FLAGS) -c $<
.fpp.f90:
	$(CPP) $(CPPFLAGS) $< $@
.c.o:
	$(CC) $(CFLAGS) -c $<

##################################################################### TARGETS

# .DEFAULT_GOAL works for GNU make 3.81 (or higher)
# For 3.80 or less, see all target
.DEFAULT_GOAL := $(GK_PROJECT)_all
ifeq ($(notdir $(CURDIR)),utils)
	.DEFAULT_GOAL := utils_all
endif
ifeq ($(notdir $(CURDIR)),geo)
	.DEFAULT_GOAL := geo_all
endif

.PHONY: all $(GK_PROJECT)_all

all: $(.DEFAULT_GOAL)

include $(DEPEND)
#include Makefile.doc_depend

ifdef USE_C_INDEX
astrogk_mod += layouts_indices.o
gs2_mod += layouts_indices.o
endif

sinclude Makefile.target_$(GK_PROJECT)

############################################################### SPECIAL RULES

# comment this out to keep intermediate .f90 files
#.PRECIOUS: $(F90FROMFPP)

.INTERMEDIATE: $(GK_PROJECT)_transforms.f90 $(GK_PROJECT)_io.f90 $(GK_PROJECT)_save.f90 \
		mp.f90 fft_work.f90

# These are special rules for the suffix problem of absoft
# (not tested)
$(GK_PROJECT)_transforms.o: $(GK_PROJECT)_transforms.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX0) -c $<
$(GK_PROJECT)_io.o: $(GK_PROJECT)_io.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX2) -c $<
$(GK_PROJECT)_save.o: $(GK_PROJECT)_save.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX2) -c $<
mp.o: mp.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX1) -c $<
fft_work.o: fft_work.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX0) -c $<

layouts_indices.o: layouts_type.h
layouts_type.h: layouts_type.f90
	$(AWK) -f makehead.awk $^ > $@

############################################################# MORE DIRECTIVES

.PHONY: depend clean distclean tar test_make

depend:
	mv input_parameters_documentation.f90 input_parameters_documentation.f90.hidden
	@$(DEPEND_CMD) -m "$(MAKE)" -1 -o -v=0 $(VPATH)
	mv input_parameters_documentation.f90.hidden input_parameters_documentation.f90

#doc_depend:
	#echo "documentation_depend = \\" > Makefile.doc_depend
	#find | grep .fpp$$ | sed s/\.fpp/\.f90\ \\\\/ | sed s/\.\\//\\t/ >> Makefile.doc_depend
#
doc: $(F90FROMFPP)
	doxygen ../doxygen/gs2
	-rm -f $(F90FROMFPP)

sync_doc: 
	rsync -av --delete ../doxygen/doc/html/	${USER},gyrokinetics@web.sourceforge.net:htdocs/gs2_documentation/

sync_input_doc:
	coderunner cc synchronise_variables . -C gs2
	curl 'http://sourceforge.net/apps/mediawiki/gyrokinetics/index.php?title=GS2_Input_Parameters&action=edit' | sed 's/&amp;/\&/g' | sed 's/&quot;/"/g' | sed 's/&gt;/>/g' | sed 's/&lt;/</g'  | sed 's/&nbsp;/ /g' > wiki_in.tmp
	coderunner cc read_mediawiki_documentation wiki_in.tmp -C gs2
	coderunner cc write_mediawiki_documentation > wiki_out.txt -C gs2

clean:
	-rm -f *.o *.mod *.g90 *.h core */core

cleanlib:
	-rm -f *.a

distclean: unlink clean cleanlib

tar:
	@[ ! -d $(TARDIR) ] || echo "ERROR: directory $(TARDIR) exists. Stop."
	@[ -d $(TARDIR) ] || $(MAKE) tar_exec

### setting tar_exec local $(TARLIST*) variables
# expand wildcards listed $(TARLIST_wild) in ( $(TARLIST_dir) + . )
# directories and add them into TARLIST
tar_exec: TARLIST = makehead.awk fortdep AstroGK.in
tar_exec: TARLIST_dir = Makefiles utils geo Aux
tar_exec: TARLIST_wild = *.f90 *.fpp *.inc *.c Makefile Makefile.* README README.*
tar_exec: TARLIST += $(foreach dir,. $(TARLIST_dir),$(wildcard $(addprefix $(dir)/,$(TARLIST_wild))))

tar_exec:
	@mkdir $(TARDIR)
	@for dir in $(TARLIST_dir) ;\
	  do ( [ ! -d $$dir ] ||  mkdir $(TARDIR)/$$dir ; ) ;\
	done
	@for name in $(TARLIST) ;\
	  do ( [ -f $$name ] && ln $$name $(TARDIR)/$$name ; ) ;\
	done
	@tar cvf - $(TARDIR) | bzip2 -9 > $(TARDIR).tar.bz2
	@rm -rf $(TARDIR)

test_make:
	@echo GK_SYSTEM is $(GK_SYSTEM)
	@echo .DEFAULT_GOAL is $(.DEFAULT_GOAL)
	@echo VPATH is $(VPATH)
	@echo CURDIR is $(CURDIR)
	@echo TOPDIR is $(TOPDIR)
	@echo NETCDF_DIR is $(NETCDF_DIR)
	@echo FFTW_DIR is $(FFTW_DIR)
	@echo
	@echo Compile mode:
	@echo  DEBUG is $(DEBUG)
	@echo  SCAL is $(SCAL)
	@echo  TEST is $(TEST)
	@echo  PROF is $(PROF)
	@echo  OPT is $(OPT)
	@echo  STATIC is $(STATIC)
	@echo  DBLE is $(DBLE)
	@echo
	@echo Functions:
	@echo  USE_MPI is $(USE_MPI)
	@echo  USE_SHMEM is $(USE_SHMEM)
	@echo  USE_FFT is $(USE_FFT)
	@echo  USE_NETCDF is $(USE_NETCDF)
	@echo  USE_HDF5 is $(USE_HDF5)
	@echo  USE_C_INDEX is $(USE_C_INDEX)
	@echo  USE_POSIX is $(USE_POSIX)
	@echo  USE_LOCAL_RAN is $(USE_LOCAL_RAN)
	@echo  USE_LOCAL_SPFUNC is $(USE_LOCAL_SPFUNC)
	@echo  USE_NAGLIB is $(USE_NAGLIB)
	@echo  DEFAULT_LIB is $(DEFAULT_LIB)
	@echo  MPI_LIB is $(MPI_LIB)
	@echo
	@echo FC is $(FC)
	@echo F90FLAGS is $(F90FLAGS)
	@echo F90OPTFLAGS is $(F90OPTFLAGS)
	@echo CC is $(CC)
	@echo CFLAGS is $(CFLAGS)
	@echo COPTFLAGS is $(COPTFLAGS)
	@echo LD is $(LD)
	@echo LDFLAGS is $(LDFLAGS)
	@echo CPP is $(CPP)
	@echo CPPFLAGS is $(CPPFLAGS)
	@echo LIBS is $(LIBS)
	@echo PLIBS is $(PLIBS)

unlink:
	-rm -f $(F90FROMFPP) layouts_type.h

revision:
	@LANG=C svn info | awk '{if($$1=="Revision:") printf("%20d",$$2) }' > Revision

TAGS:	*.f90 *.fpp */*.f90 */*.fpp
	etags $^
