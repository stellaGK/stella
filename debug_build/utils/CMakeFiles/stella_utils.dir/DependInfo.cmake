
# Consider dependencies only in project.
set(CMAKE_DEPENDS_IN_PROJECT_ONLY OFF)

# The set of languages for which implicit dependencies are needed:
set(CMAKE_DEPENDS_LANGUAGES
  "Fortran"
  )
# The set of files for implicit dependencies of each language:
set(CMAKE_DEPENDS_CHECK_Fortran
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/command_line.fpp" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/command_line.fpp.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/constants.fpp" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/constants.fpp.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/convert.f90" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/convert.f90.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/fft_work.f90" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/fft_work.f90.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/file_utils.fpp" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/file_utils.fpp.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/gauss_quad.f90" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/gauss_quad.f90.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/job_manage.fpp" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/job_manage.fpp.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/linear_solve.f90" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/linear_solve.f90.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/mp.fpp" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/mp.fpp.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/mp_lu_decomposition.fpp" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/mp_lu_decomposition.fpp.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/mt19937.f90" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/mt19937.f90.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/netcdf_utils.fpp" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/netcdf_utils.fpp.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/ran.fpp" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/ran.fpp.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/redistribute.f90" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/redistribute.f90.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/smooth_step.f90" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/smooth_step.f90.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/sort.f90" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/sort.f90.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/spfunc.fpp" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/spfunc.fpp.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/spl.f90" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/spl.f90.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/system_fortran.fpp" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/system_fortran.fpp.o"
  "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/utils/text_options.f90" "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/CMakeFiles/stella_utils.dir/text_options.f90.o"
  )
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_SUBMODULE_SEP "@")
set(CMAKE_Fortran_SUBMODULE_EXT ".smod")

# Preprocessor definitions for this target.
set(CMAKE_TARGET_DEFINITIONS_Fortran
  "ANSI_CPP"
  "F200X_INTRINSICS"
  "FCOMPILER=_INTEL_"
  "FFT=_FFTW3_"
  "ISO_C_BINDING"
  "MPI"
  "NETCDF=_DEFAULT_"
  "SPFUNC=_SPF200X_"
  )

# The include file search paths:
set(CMAKE_Fortran_TARGET_INCLUDE_PATH
  "utils/mod"
  "/usr/include"
  "/cineca/prod/opt/libraries/netcdff/4.4.4/intel--pe-xe-2018--binary/include"
  "/cineca/prod/opt/libraries/netcdf/4.6.1/intel--pe-xe-2018--binary/include"
  "/marconi/prod/opt/compilers/intel/pe-xe-2018/binary/compilers_and_libraries_2018.5.274/linux/mpi/intel64/include"
  )

# The set of dependency files which are needed:
set(CMAKE_DEPENDS_DEPENDENCY_FILES
  )

# Targets to which this target links.
set(CMAKE_TARGET_LINKED_INFO_FILES
  )

# Fortran module output directory.
set(CMAKE_Fortran_TARGET_MODULE_DIR "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/utils/mod")
