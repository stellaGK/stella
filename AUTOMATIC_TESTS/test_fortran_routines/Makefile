

# Directory
PFUNIT_TEST_DIR = $(PROJECT_DIR)/AUTOMATIC_TESTS/test_fortran_routines
PFUNIT_TEST_PATH :=.

# Find the pFUnit tests inside $(PROJECT_DIR)/AUTOMATIC_TESTS/test_fortran_routines/
PFUNIT_TESTS:= $(shell find $(PFUNIT_TEST_DIR) -name '[^.\#]*.pf')
PFUNIT_TEST_LIST := $(PFUNIT_TEST_DIR)/testSuites.inc

# Number of processors
PFUNIT_NPROC = 4

# Default major version of pFUnit
# Versions 3 and 4 require different flags, etc.
PFUNIT_VERSION_MAJOR ?= 4
PFUNIT_LOG:=pfunit.log
PFUNIT_ERROR:=pfunit.error

# Installation directory of pFUnit
# Default is the directory where `scripts/build_pfunit` puts it
# Note that the folder was called PFUNIT-4.2 when a previous version of pFUnit was used.
PFUNIT_DIR ?= $(COMPILATION_DIR)/pFUnit_build/install/PFUNIT-4.10
PFUNIT_DRIVER := $(PFUNIT_DIR)/include/driver.F90

# Preprocessor flags
CPPFLAGS += -DPFUNIT_VERSION_MAJOR=$(PFUNIT_VERSION_MAJOR)

# Flags based on pFUnit version
ifeq ($(PFUNIT_VERSION_MAJOR),4) 
  sinclude $(PFUNIT_DIR)/include/PFUNIT.mk 
  PFUNIT_MOD := pfunit
  PFUNIT_PREPROCESSOR ?= $(PFUNIT_BIN_DIR)/funitproc
  PFUNIT_INC := -I$(PFUNIT_INCLUDE_DIR)
  PFUNIT_LIB := -L$(PFUNIT_TOP_DIR)/lib -lpfunit -lfunit $(FARGPARSE_LIBRARIES)
  F90FLAGS += $(PFUNIT_EXTRA_FFLAGS)
else
  PFUNIT_MOD := pfunit_mod
  PFUNIT_PREPROCESSOR ?= $(PFUNIT_DIR)/bin/pFUnitParser.py
  PFUNIT_INC := -I$(PFUNIT_DIR)/mod
  PFUNIT_LIB := -L$(PFUNIT_DIR)/lib -lpfunit
endif

# Preprocessor flags
CPPFLAGS += -DPFUNIT_MOD=$(PFUNIT_MOD)

# We need the stella fortran objects, which should already be compiled in $(OBJ)
# Note that we need to remove stella, because it is a "main" program, if we add
# this together with $(PFUNIT_DRIVER), the compiler will not know which one to
# choose for the main program of the executable "automated-fortran-tests"
STELLA_OBJS := $(filter-out $(OBJ)/stella.o, $(stella_mod))

# To build the pFUnit executable, we need to build the 
# PFUNIT_TESTS.pf into PFUNIT_TESTS.F90 and then into PFUNIT_TESTS.o files
PFUNIT_OBJS := $(PFUNIT_TESTS:.pf=.o)
PFUNIT_INTERMEDIATES := $(PFUNIT_TESTS:.pf=.F90)
.INTERMEDIATE: $(PFUNIT_INTERMEDIATES)

###############################################################
# pFUnit library
###############################################################
# Commands to build, check and clean the <pfunit> library

.PHONY: build-pfunit-library check-pfunit-install clean-pfunit-library

build-pfunit-library:
	@echo "Building pFUnit in default location: $(COMPILATION_DIR)/pFUnit_build/install"
	@AUTOMATIC_TESTS/test_fortran_routines/build_pfunit -v ${PFUNIT_VERSION_MAJOR}

check-pfunit-install:
ifeq ("$(wildcard $(PFUNIT_PREPROCESSOR))","")
	$(info "Please install the pfunit package through:")
	$(info "  >> make build-pfunit-library")
	$(info "      ")
	$(error Could not find $(PFUNIT_PREPROCESSOR). Make sure pFUnit is installed and PFUNIT_DIR ($(PFUNIT_DIR)) is set to the current location)
endif

clean-pfunit-library:
	@echo "Removing pFUnit library"
	@AUTOMATIC_TESTS/test_fortran_routines/build_pfunit -r

distclean: clean-pfunit-library

###############################################################
# Build tests
############################################################### 

# The -I $(MOD) flag makes sure that the compiler looks for the *.mod files in the $(MOD) directory 
MOD_FLAGS_TESTS := -I $(OBJ) -I $(UTILS_OBJ)  -I $(MOD) -I $(UTILS_MOD) 
#-I $(MINILIBSTELL_MOD)

%.o: %.pf
	@echo "   (test_fortran_routines) %.pf to %.o $<"
	@$(PFUNIT_PREPROCESSOR) $< $(basename $<).F90
ifdef ONE_STEP_PP
	@$(QUIETSYM)$(FC) $(FPPFLAGS) $(F90FLAGS) -c $(basename $<).F90 $(PFUNIT_INC) -I $(PROJECT_DIR) -o $@
else
	@$(QUIETSYM)$(CPP) $(FPPFLAGS) $(basename $<).F90 -o $(basename $<).f90
	@$(QUIETSYM)$(FC) $(F90FLAGS) -c $(basename $<).f90 $(PFUNIT_INC) $(MOD_FLAGS_TESTS) -o $@
endif

$(PFUNIT_DRIVER): $(PFUNIT_TEST_LIST)
	@touch $(PFUNIT_DRIVER)

$(PFUNIT_OBJS): $(STELLA_OBJS)

###############################################################
# Test suite
###############################################################

# The test suite executable recipe to execute the Fortran tests 
GIT_OBJ :=  $(OBJ_DIR)/git_version_impl.o
automated-fortran-tests: $(GIT_OBJ) $(PFUNIT_OBJS) $(PFUNIT_DRIVER) $(STELLA_OBJS)
	@echo "                       "
	@echo "Compiling pFUnit tests executable:"
	@$(FC) -o $@ $^ -DUSE_MPI=YES $(PFUNIT_INC) $(PFUNIT_LIB) -I. $(SIMPLEDATAIO_LIB_ABS) $(LIBS) $(F90FLAGS) -I $(PFUNIT_TEST_DIR) $(MOD_FLAGS_TESTS) -D_TEST_SUITES='"$(PFUNIT_TEST_LIST)"'
	@echo "   ==> Created the executable automated-fortran-tests"

# To run the automated Fortran tests, we first need to check
# that the pfunit package has been installed correctly, and
# then we can run the automated Fortran tests
.PHONY: run-automated-fortran-tests run-automated-fortran-tests-step1
run-automated-fortran-tests:  
	@echo "                       "
	@echo "Make sure that stella is compiled:"
	@echo "  >> make clean; make -j 10"
	@echo "                       "
	@echo "In order to run the Fortran test make sure the pfunit package is build:"
	@echo "  >> make build-pfunit-library"
	@echo "                       "
	@echo "And that the python virtual environment is turned on:"
	@echo "  >> make create-test-virtualenv"
	@echo "  >> source AUTOMATIC_TESTS/venv/bin/activate"
	@echo "                       "
	@echo "Create testing scripts:"
	@$(MAKE) --no-print-directory run-automated-fortran-tests-step1 
run-automated-fortran-tests-step1: check-pfunit-install automated-fortran-tests 
	@echo "                       "
	@echo "======================================================================"
	@echo "Running: pFUnit tests"
	@mpirun -np $(PFUNIT_NPROC) $(PFUNIT_TEST_PATH)/automated-fortran-tests
	@echo "======================================================================"

###############################################################
# Coverage
###############################################################

.PHONY: coverage pfunit-coverage run-gcov

coverage:
	gcov *.f90

ifndef USE_GCOV
coverage-pfunit: coverage_gcov_fail
coverage_gcov_fail:
	$(error Cannot make coverage-pfunit as USE_GCOV undefined. Rebuild everything with USE_GCOV=on)
else
coverage-pfunit: stella automated-fortran-tests coverage
	mkdir -p coverage-pfunit
	mv *.gcov coverage-pfunit/. || echo "Couldn't move generated gcov files to coverage-pfunit. This is likely because gcov files not found in ${PROJECT_DIR}."
endif

###############################################################
# Clean up
###############################################################

.PHONY: clean-pfunit

clean-pfunit:
	@echo "   - cleaning pFUnit tests" 
	@$(RM) automated-fortran-tests $(PFUNIT_TEST_DIR)/*.o $(PFUNIT_TEST_DIR)/*.mod $(PFUNIT_TEST_DIR)/*.f90 \
      $(PFUNIT_TEST_DIR)/*~ $(PFUNIT_TEST_DIR)/*.gcov $(PFUNIT_TEST_DIR)/*.gcda \
      $(PFUNIT_TEST_DIR)/*.gcno $(PFUNIT_TEST_DIR)/*.log $(PFUNIT_TEST_DIR)/*.error \
      $(PFUNIT_OBJS) $(PFUNIT_INTERMEDIATES) ./test_*mod ./wraptest_*mod
	@echo "                     "
