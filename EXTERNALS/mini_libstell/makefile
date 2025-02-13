
# Note that we can not compile mini_libstell directly since it depends on the 
# variables $(F90FLAGS) and $(NETCDF_INC) which are set through the main stella makefile. 
# Moreover, we make use of $(MINILIBSTELL_DIR) and $(COMPILATION_DIR)

# We set the .DEFAULT_GOAL and for GNU make 3.80 or less, we use the 'all' target 
.DEFAULT_GOAL := compile-mini_libstell 
all: compile-mini_libstell

# Targets which do not create files
.PHONY: all compile-mini_libstell clean create-directories

#===================================================================
#                          COMPILER FLAGS                          #
#===================================================================

# If F90FLAGS has not been set (e.g. running makefile from geo/vmec_interface directory), 
# then set it to 'empty' so we can check for it later
F90FLAGS ?= empty

# Load up the compiler flags from the main stella makefile. 
# We need to include the NETCDF paths for mini_libstell.
EXTRA_COMPILE_FLAGS = $(F90FLAGS) $(NETCDF_INC)
EXTRA_LINK_FLAGS = $(NETCDF_LIB)

#===================================================================
#                            DIRECTORIES                           #
#===================================================================

# Location of the code 
SRC := $(MINILIBSTELL_DIR)

# Keep the compilated files in seperate folders
BUILD := $(COMPILATION_DIR)/build_make
PARENT := $(COMPILATION_DIR)/build_make/mini_libstell
MOD := $(PARENT)/modules
OBJ := $(PARENT)/objects

# Create the directories which will contain the compiled files
# These commands will be triggered if $(BUILD) is a prerequisite, but it does not exist yet.
$(PARENT) $(MOD) $(OBJ) create-directories:
	@echo "                     " 
	@echo "   Create the COMPILATION/build_make/mini_libstell/* directories"
	@echo "      - $(PARENT) "
	@echo "      - $(MOD) "
	@echo "      - $(OBJ) "
	@mkdir -p $(BUILD) 
	@mkdir -p $(PARENT) 
	@mkdir -p $(MOD)
	@mkdir -p $(OBJ)  
	@echo "                     "

#===================================================================
#                      MINI_LIBSTELL LIBRARY                       #
#===================================================================

# Gather the objects that need to be created
libstell_objects_names = ezcdf.o ezcdf_inqvar.o handle_err.o stel_constants.o ezcdf_GenGet.o \
	ezcdf_opncls.o read_wout_mod.o safe_open_mod.o stel_kinds.o vsvd0.o parse_extension.o  
libstell_objects = $(addprefix $(OBJ)/, $(libstell_objects_names)) 

# Put messages around the build command
compile-mini_libstell:  
	@$(MAKE) --no-print-directory create-directories
	@echo "   Compile mini_libstell objects:" 
	@$(MAKE) --no-print-directory -s $(BUILD)/mini_libstell.a
	@echo "          ==> Finished compiling mini_libstell objects"; echo "  "

# Create the static library, once all objects have been made 
$(BUILD)/mini_libstell.a: $(libstell_objects) | $(OBJ)
	@ar rcs $@ $(libstell_objects) 
	
# The -J $(MOD) flag makes sure that the compiler writes the *.mod files to the $(MOD) directory
# The -I $(MOD) flag makes sure that the compiler looks for the *.mod files in the $(MOD) directory
MOD_FLAGS_MINI := -I $(MOD) $(MODULE_FLAG) $(MOD)

# Pattern rules to compile the scripts
# For the preprocessed files *.F, we need to give the "NETCDF" option through the -D flag
# Without this flag, stella will compile and run but the netcdf file won't be written correctly
$(OBJ)/%.o: $(SRC)/%.f90 | $(OBJ)
	@echo "       (mini_libstell) *.f90 to *.o:  $(notdir $<)"
	@$(FC) $(EXTRA_COMPILE_FLAGS) $(MOD_FLAGS_MINI) -o $@ -c $<

$(OBJ)/%.o: $(SRC)/%.f | $(OBJ)
	@echo "       (mini_libstell) *.f to *.o:    $(notdir $<)"
	@$(FC) $(EXTRA_COMPILE_FLAGS) $(MOD_FLAGS_MINI) -o $@ -c $<

$(OBJ)/%.o: $(SRC)/%.F | $(OBJ)
	@echo "       (mini_libstell) *.F to *.o:    $(notdir $<)"
	@$(FC) $(EXTRA_COMPILE_FLAGS) $(MOD_FLAGS_MINI) -DNETCDF -o $@ -c $<

# The "ezcdf_GenGet.f90" script has a rank mistmatch warning which we will supress.
# Warning: Type mismatch between actual argument at (1) and actual argument at (2) (REAL(4)/COMPLEX(4)).
# The -w flag will suppress the warnings, while allowing for errors
$(OBJ)/ezcdf_GenGet.o: $(SRC)/ezcdf_GenGet.f90 | $(OBJ)
	@echo "       (mini_libstell) *.f90 to *.o:  $(notdir $<)"
	@$(FC) $(EXTRA_COMPILE_FLAGS) $(MOD_FLAGS_MINI) -w -o $@ -c $<
	
# Remove the compiled files 
clean:
	@echo "   - remove mini_libstell"
	@rm -f $(OBJ)/*.o $(MOD)/*.mod $(MOD)/*.MOD $(BUILD)/mini_libstell.a
	@rm -f *.o *.mod *.MOD *~ 
	@if [ -d "$(MOD)" ]; then rmdir $(MOD); fi
	@if [ -d "$(OBJ)" ]; then rmdir $(OBJ); fi
	@if [ -d "$(PARENT)" ]; then rmdir $(PARENT); fi  

# Dependencies
$(OBJ)/ezcdf.o: $(OBJ)/ezcdf_GenGet.o $(OBJ)/ezcdf_opncls.o
$(OBJ)/ezcdf_inqvar.o:
$(OBJ)/handle_err.o:
$(OBJ)/stel_constants.o: $(OBJ)/stel_kinds.o
$(OBJ)/ezcdf_GenGet.o: $(OBJ)/ezcdf_opncls.o $(OBJ)/ezcdf_inqvar.o
$(OBJ)/ezcdf_opncls.o:
$(OBJ)/read_wout_mod.o: $(OBJ)/stel_kinds.o $(OBJ)/safe_open_mod.o $(OBJ)/stel_constants.o $(OBJ)/vsvd0.o $(OBJ)/ezcdf.o $(OBJ)/parse_extension.o
$(OBJ)/safe_open_mod.o:
$(OBJ)/stel_kinds.o:
$(OBJ)/vsvd0.o: $(OBJ)/stel_kinds.o
$(OBJ)/parse_extension.o:
