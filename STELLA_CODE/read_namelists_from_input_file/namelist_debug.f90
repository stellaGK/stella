!###############################################################################
!#################### READ STELLA NAMELISTS FOR DEBUG FLAGS ####################
!###############################################################################
! 
! This module will read the namelists associated with the debug flags:
! 
!   debug_flags
!      print_extra_info_to_terminal = .true.
!      debug_all = .false.
!      stella_debug = .false.
!      fields_all_debug = .false.
!      fields_debug = .false.
!      fields_fluxtube_debug = .false.
!      fields_electromagnetic_debug = .false.
!      fields_ffs_debug = .false.
!      implicit_solve_debug = .false.
!      parallel_streaming_debug = .false.
!      mirror_terms_debug = .false.
!      neoclassical_terms_debug = .false.
!      response_matrix_debug = .false.
!      time_advance_debug = .false.
!      extended_grid_debug = .false.
!      diagnostics_all_debug = .false.
!      diagnostics_parameters = .false.
!      diagnostics_fluxes_fluxtube_debug = .false.
!      diagnostics_omega_debug = .false.
!      diagnostics_debug = .false.
!      dist_fn_debug = .false.
!      gyro_averages_debug = .false.
!      fluxes_debug = .false.
!      geometry_debug = .false.
!      const_alpha_geo = .false.
! 
! For each namelists two (or three) routines exist:
!    - set_default_parameters_<namelist>
!    - read_namelist_<namelist>
!    - check_inputs_<namelist>
! 
! First the default input parameters are set, then the default options are
! overwritten with those specified in the input file. Optionally, it is
! checked whether any input variables are clashing with each other.
! 
!###############################################################################
module namelist_debug

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_debug_flags

   private

   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                                DEBUG FLAGS                                !
   !****************************************************************************
   subroutine read_namelist_debug_flags(debug_all, stella_debug, fields_all_debug, fields_debug, &
      fields_fluxtube_debug, fields_electromagnetic_debug, fields_ffs_debug, implicit_solve_debug, & 
      parallel_streaming_debug, mirror_terms_debug, neoclassical_terms_debug, response_matrix_debug, &
      time_advance_debug, extended_grid_debug, diagnostics_all_debug, diagnostics_parameters, &
      diagnostics_fluxes_fluxtube_debug, diagnostics_omega_debug, diagnostics_debug, dist_fn_debug, &
      gyro_averages_debug, fluxes_debug, geometry_debug,  const_alpha_geo, print_extra_info_to_terminal, &
      debug_test_gyro_average, calculations_debug)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      logical, intent(out) :: debug_all, stella_debug, fields_all_debug, fields_debug, &
         fields_fluxtube_debug, fields_electromagnetic_debug, fields_ffs_debug, & 
         implicit_solve_debug, parallel_streaming_debug, mirror_terms_debug, neoclassical_terms_debug, &
         response_matrix_debug, time_advance_debug, extended_grid_debug, &
         diagnostics_all_debug, diagnostics_parameters, diagnostics_fluxes_fluxtube_debug, &
         diagnostics_omega_debug, diagnostics_debug, dist_fn_debug,&
         gyro_averages_debug, fluxes_debug, geometry_debug, const_alpha_geo, print_extra_info_to_terminal, &
         debug_test_gyro_average, calculations_debug

      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_debug_flags
      call read_input_file_debug_flags
      call check_inputs_debug_flags

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_debug_flags

         implicit none

         ! Flags that turn on all debug flags
         debug_all = .false.
         fields_all_debug = .false.
         diagnostics_all_debug = .false.

         ! Debug flags
         stella_debug = .false.
         fields_debug = .false.
         fields_fluxtube_debug = .false.
         fields_electromagnetic_debug = .false.
         fields_ffs_debug = .false.
         implicit_solve_debug = .false.
         mirror_terms_debug = .false.
         neoclassical_terms_debug = .false.
         parallel_streaming_debug = .false.
         response_matrix_debug = .false.
         time_advance_debug = .false.
         extended_grid_debug = .false.
         geometry_debug = .false.
         dist_fn_debug = .false.
         gyro_averages_debug = .false.

         ! Diagnostics debug flags
         diagnostics_debug = .false.
         diagnostics_parameters = .false.
         diagnostics_omega_debug = .false.
         diagnostics_fluxes_fluxtube_debug = .false.
         fluxes_debug = .false.

         ! Specific flags to debug certain routines
         const_alpha_geo = .false.
         debug_test_gyro_average = .false.
         calculations_debug = .false.

         ! Print stella header to the terminal
         print_extra_info_to_terminal = .true.

      end subroutine set_default_parameters_debug_flags

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_debug_flags

         use file_utils, only: input_unit_exist, error_unit

         implicit none

         ! Variables in the <debug_flags> namelist
         namelist /debug_flags/ debug_all, stella_debug, fields_all_debug, fields_debug, &
             fields_fluxtube_debug, fields_electromagnetic_debug, fields_ffs_debug, & 
             implicit_solve_debug, parallel_streaming_debug, mirror_terms_debug, neoclassical_terms_debug, &
             response_matrix_debug, time_advance_debug, extended_grid_debug, &
             diagnostics_all_debug, diagnostics_parameters, diagnostics_fluxes_fluxtube_debug, &
             diagnostics_omega_debug, diagnostics_debug, dist_fn_debug,&
             gyro_averages_debug, fluxes_debug, geometry_debug,  const_alpha_geo, print_extra_info_to_terminal, &
             debug_test_gyro_average, calculations_debug
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('debug_flags', dexist)
         if (dexist) read (unit=in_file, nml=debug_flags)

      end subroutine read_input_file_debug_flags

      !---------------------------- Check variables ----------------------------
      subroutine check_inputs_debug_flags

         use file_utils, only: error_unit

         implicit none
         
         !----------------------------------------------------------------------

         if (debug_all) then
             stella_debug = .true.
             fields_all_debug = .true.
             time_advance_debug = .true.
             implicit_solve_debug = .true.
             parallel_streaming_debug = .true.
             response_matrix_debug = .true.
             mirror_terms_debug = .true.
             neoclassical_terms_debug = .true.
             extended_grid_debug = .true.
             diagnostics_all_debug = .true.
             dist_fn_debug = .true.
             gyro_averages_debug = .true.
             geometry_debug = .true.
             print_extra_info_to_terminal = .true.
             calculations_debug = .true.
         end if

         if (diagnostics_all_debug) then
             diagnostics_debug = .true.
             diagnostics_parameters = .true.
             diagnostics_fluxes_fluxtube_debug = .true.
             diagnostics_omega_debug = .true.
             fluxes_debug = .true.
         end if

         if (fields_all_debug) then
             fields_debug = .true.
             fields_fluxtube_debug = .true.
             fields_electromagnetic_debug = .true.
             fields_ffs_debug = .true.
         end if 

      end subroutine check_inputs_debug_flags

   end subroutine read_namelist_debug_flags

end module namelist_debug
