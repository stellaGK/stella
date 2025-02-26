

# We set the .DEFAULT_GOAL and for GNU make 3.80 or less, we use the 'all' target 
.DEFAULT_GOAL := compile-utils
all: compile-utils

# Targets which do not create files
.PHONY: all compile-utils utils clean create-directories

#===================================================================
#                            DIRECTORIES                           #
#===================================================================

# Location of the code 
SRC := $(UTILS_DIR)

# Keep the compilated files in seperate folders
BUILD := $(COMPILATION_DIR)/build_make
PARENT := $(COMPILATION_DIR)/build_make/utils
FPP := $(UTILS_FPP)
OBJ := $(UTILS_OBJ)
MOD := $(UTILS_MOD)
OBJ_GENERAL := $(OBJ_DIR)

# Create the directories which will contain the compiled files
# These commands will be triggered if $(BUILD) is a prerequisite, but it does not exist yet.
$(PARENT) $(MOD) $(OBJ) $(FPP) create-directories:
	@echo "    " 
	@echo "   Create the COMPILATION/build_make/mini_libstell/* directories"
	@echo "      - $(PARENT) "
	@echo "      - $(MOD) "
	@echo "      - $(OBJ) "
	@echo "      - $(FPP) "
	@mkdir -p $(BUILD) 
	@mkdir -p $(PARENT) 
	@mkdir -p $(OBJ_GENERAL) 
	@mkdir -p $(MOD)
	@mkdir -p $(OBJ)  
	@mkdir -p $(FPP)  
	@echo "    "

#===================================================================
#                          UTILS SCRIPTS                           #
#===================================================================

# The 'utils' directory contains all external python scripts 
# We compile them as a seperate package, since usually we do not need to recompile it
utils_objects_names_temp = $(subst .f90,.o,$(subst .fpp,.o,$(notdir $(wildcard $(SRC)/*))))
utils_objects_names = $(subst makefile,,$(subst define.inc,,$(subst CMakeLists.txt,,$(utils_objects_names_temp))))
utils_objects = $(addprefix $(OBJ_GENERAL)/, $(utils_objects_names)) 
	
# Put messages around the build command
compile-utils:   
	@$(MAKE) --no-print-directory create-directories
	@$(MAKE) --no-print-directory ensure-that-all-f90-files-exist
	@echo "  "; echo "   Compile utils objects:" 
	@$(MAKE) --no-print-directory -s $(BUILD)/utils.a
	@echo "          ==> Finished compiling utils objects"; echo "  "

# Create the static library, once all objects have been made 
# Note that we won't actually use the library, however, to make sure we recompile
# the objects when the file has been touched, it's necessary to have a single object
# like "utils.a" which the Make.file.externals can test for, otherwise parallelisation is a mess
$(BUILD)/utils.a: $(utils_objects) | $(OBJ)
	$(ARCH) $(ARCHFLAGS) $@ $^
	$(RANLIB) $@

#===================================================================
#                           COMPILE F90                            #
#===================================================================
	
# The -J $(MOD) flag makes sure that the compiler writes the *.mod files to the $(MOD) directory
# The -I $(MOD) flag makes sure that the compiler looks for the *.mod files in the $(MOD) directory
MOD_FLAGS_UTILS := -I $(MOD) $(MODULE_FLAG) $(MOD)

# Pattern rules to compile the scripts
# Create *.o files from *.f90 files using the $(FC) compiler 
# Note that the *.f90 files can be located in $(SRC) or $(FPP)
# Where the $(FPP) directory contains the prepocessed *.f90 from *.fpp files 
$(OBJ)/%.o: $(SRC)/%.f90 | $(OBJ)
	@echo "       (utils) *.f90 to *.o: $(notdir $<)" 
	@$(FC) $(F90FLAGS) $(INC_FLAGS) $(MOD_FLAGS_UTILS) -o $@ -c $< 
$(OBJ)/%.o: $(FPP)/%.f90 | $(OBJ)
	@echo "       (utils) *.f90 to *.o: $(notdir $<)" 
	@$(FC) $(F90FLAGS) $(INC_FLAGS) $(MOD_FLAGS_UTILS) -o $@ -c $< 
	
# Make the utils objects accesible in the main objects folder
$(OBJ_GENERAL)/%.o: $(OBJ)/%.o 
	@ln -sf $< $(OBJ_GENERAL)	
	
.PRECIOUS: $(OBJ)/%.o $(OBJ_GENERAL)/%.o
	
#===================================================================
#                           COMPILE FPP                            #
#===================================================================

# Create *.f90 files from *.fpp files using the C preprocessor (CPP) compiler
# Before we can do this, the $(FPP)/%.fpp file needs to exist in the $(FPP) directory. 
# A pattern rule is used to create symbolic links from EXTERNALS/utils/*.fpp to $(FPP)/*.fpp
# Moreover the $(FPP)/define.inc needs to exist, which contains the $(CPP) settings
$(FPP)/%.f90: $(FPP)/%.fpp $(FPP)/define.inc | $(FPP)
	@echo "       (utils) *.fpp to *.f90: $(notdir $<)"
	@$(CPP) $(CPPFLAGS) $< $@

# In EXTERNALS/utils/*, we have *.fpp files, that need to be preprocessed into *.f90 files
# We will use a pattern rule to create symbolic links of these *.fpp files in the $(BUILD_F90_DIR) directory
# And we define the *.f90 file paths inside $(BUILD_F90_DIR), since we will preprocess $(BUILD_F90_DIR)/*.fpp to $(BUILD_F90_DIR)/*.f90
existing_fpp_files = $(wildcard $(SRC)/*.fpp) $(SRC)/define.inc
existing_fpp_file_names := $(notdir $(existing_fpp_files))
fpp_files_in_build_directory := $(addprefix $(FPP)/, $(existing_fpp_file_names)) 
f90_from_fpp_files = $(fpp_files_in_build_directory:.fpp=.f90)

# Make sure that all *.fpp files have been preprocessed so that the parallelisation doesn't run into issues
.PHONY: ensure-that-all-f90-files-exist ensure-that-all-f90-files-exist-step-1 ensure-that-all-f90-files-exist-step-2
ensure-that-all-f90-files-exist: | $(FPP) 
	@echo "   Preprocess *.fpp files to create *.f90 files" 
	@$(MAKE) --no-print-directory ensure-that-all-f90-files-exist-step-1
	@$(MAKE) --no-print-directory ensure-that-all-f90-files-exist-step-2
ensure-that-all-f90-files-exist-step-1: $(fpp_files_in_build_directory)  
	@echo "          ==> Created symbolic links"   
ensure-that-all-f90-files-exist-step-2: $(f90_from_fpp_files)  
	@echo "          ==> Preprocessed all *.fpp files"  
		
# Create symbolic links from EXTERNALS/utils/*.fpp to $(FPP)/*.fpp 
$(FPP)/%.fpp: $(SRC)/%.fpp | $(FPP)
	@echo "       (utils) ln -s EXTERNALS/utils/$(notdir $<) COMPILATION/utils/f90fromfpp/$(notdir $<)"	
	@ln -sf $< $(FPP)	
	
# Create a symbolic link from $(SRC)/define.inc to $(FPP)/define.inc
# Before we can do this, we need to make sure that the build directory $(FPP) exist
$(FPP)/define.inc: | $(FPP)
	@echo "       (utils) ln -s EXTERNALS/utils/define.inc COMPILATION/utils/f90fromfpp/define.inc"	
	@ln -sf $(SRC)/define.inc $(FPP)/define.inc
	
# Add the fpp folder to the VPATH for the fortdep module
export VPATH := $(VPATH):$(FPP):
	
#===================================================================
#                              CLEAN                               #
#===================================================================

# Remove the compiled files 
clean:
	@echo "   - remove utils"
	@rm -f $(OBJ)/*.o $(MOD)/*.mod $(MOD)/*.MOD $(FPP)/define.inc $(FPP)/*.f90 $(FPP)/*.fpp     
	@rm -f *.o *.mod *.MOD *~ $(BUILD)/utils.a
	@if [ -d "$(MOD)" ]; then rmdir $(MOD); fi
	@if [ -d "$(OBJ)" ]; then rmdir $(OBJ); fi
	@if [ -d "$(FPP)" ]; then rmdir $(FPP); fi  
	@if [ -d "$(PARENT)" ]; then rmdir $(PARENT); fi  


#===================================================================
#                           DEPENDENCIES                           #
#===================================================================

# Dependencies
$(OBJ)/command_line.o: 
$(OBJ)/constants.o: 
$(OBJ)/convert.o: 
$(OBJ)/fft_work.o: $(OBJ)/mp.o $(OBJ)/constants.o
$(OBJ)/file_utils.o: $(OBJ)/command_line.o 
$(OBJ)/gauss_quad.o: $(OBJ)/file_utils.o $(OBJ)/mp.o $(OBJ)/constants.o
$(OBJ)/job_manage.o: $(OBJ)/file_utils.o $(OBJ)/mp.o
$(OBJ)/linear_solve.o:  
$(OBJ)/mp.o: $(OBJ)/file_utils.o $(OBJ)/constants.o
$(OBJ)/mp_lu_decomposition.o: $(OBJ)/linear_solve.o
$(OBJ)/mt19937.o: 
$(OBJ)/netcdf_utils.o: $(OBJ)/constants.o $(OBJ)/mp.o $(OBJ)/file_utils.o
$(OBJ)/ran.o: $(OBJ)/mt19937.o
$(OBJ)/redistribute.o: $(OBJ)/mp.o $(OBJ)/job_manage.o
$(OBJ)/spfunc.o: $(OBJ)/constants.o
$(OBJ)/spl.o: $(OBJ)/constants.o
$(OBJ)/system_fortran.o: 
$(OBJ)/text_options.o: $(OBJ)/mp.o

