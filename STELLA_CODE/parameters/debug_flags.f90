!###############################################################################
!############################### READ DEBUG FLAGS ##############################
!###############################################################################
! Namelist: &debug_flags
! These flags will allow you to toggle the debug flags in each file of stella.
! If you would like to debug a specific file then the corresponding flag to turn
! the debug messages in that file will be <FileNameStem_debug>.
! As an example, if you want to debug the fields.fpp debug flags then the
! corresponding debug flag in the input file will be <fields_debug = .true.>.
! If you would like all the debug flags on in every file then simply set 
! <debug_all = .true.> in the input file; this will set every debug flag to true
!###############################################################################
module debug_flags

  implicit none
  
  ! Public routines
  public :: read_debug_flags
  
  ! Debug flags
  public :: stella_debug
  public :: fields_debug
  public :: fields_fluxtube_debug
  public :: fields_electromagnetic_debug
  public :: fields_ffs_debug
  public :: time_advance_debug
  public :: implicit_solve_debug
  public :: parallel_streaming_debug
  public :: response_matrix_debug
  public :: mirror_terms_debug
  public :: neoclassical_terms_debug
  public :: extended_grid_debug
  public :: diagnostics_debug
  public :: diagnostics_parameters
  public :: diagnostics_fluxes_fluxtube_debug
  public :: diagnostics_omega_debug
  public :: fluxes_debug
  public :: geometry_debug
  public :: dist_fn_debug
  public :: gyro_averages_debug
  public :: ffs_solve_debug
  
  ! Debug flug for debugging full flux surface
  ! This will set the geometry of all field lines to be equation to
  ! the geometry on alpha = 0 (i.e. the first field line) 
  public :: const_alpha_geo
  
  private

  ! Debug flags
  logical :: debug_all 
  logical :: dist_fn_debug
  logical :: gyro_averages_debug 
  logical :: initialised = .false.
  
  ! Main stella file debug flag
  logical :: stella_debug
  
  ! Fields debug flags
  logical :: fields_all_debug
  logical :: fields_debug
  logical :: fields_fluxtube_debug
  logical :: fields_electromagnetic_debug
  logical :: fields_ffs_debug
  
  ! Gyrokinetic term debug flags
  logical :: time_advance_debug
  logical :: implicit_solve_debug
  logical :: parallel_streaming_debug
  logical :: response_matrix_debug
  logical :: mirror_terms_debug
  logical :: neoclassical_terms_debug
  logical :: extended_grid_debug
  
  ! Diagnostics debug flags
  logical :: diagnostics_all_debug
  logical :: diagnostics_debug
  logical :: diagnostics_parameters
  logical :: diagnostics_fluxes_fluxtube_debug
  logical :: diagnostics_omega_debug
  logical :: fluxes_debug
  
  ! Geometry debug flags
  logical :: geometry_debug
  
  ! For FFS
  logical :: ffs_solve_debug
  logical :: const_alpha_geo

contains

  !======================================================================
  !========================= READ DEBUG FLAGS  ==========================
  !======================================================================
  subroutine read_debug_flags

    use mp, only: proc0
    use input_file_debug, only: read_namelist_debug_flags  
    implicit none

    call read_namelist_debug_flags (debug_all, stella_debug, ffs_solve_debug, fields_all_debug, fields_debug, &
                    fields_fluxtube_debug, fields_electromagnetic_debug, fields_ffs_debug, & 
                    implicit_solve_debug, parallel_streaming_debug, mirror_terms_debug, neoclassical_terms_debug, &
                    response_matrix_debug, time_advance_debug, extended_grid_debug, &
                    diagnostics_all_debug, diagnostics_parameters, diagnostics_fluxes_fluxtube_debug, &
                    diagnostics_omega_debug, diagnostics_debug, dist_fn_debug,&
                    gyro_averages_debug, fluxes_debug, geometry_debug,  const_alpha_geo)

    call broadcast_parameters

  contains

      !**********************************************************************
      !                         BROADCAST OPTIONS                           !
      !**********************************************************************
      ! Broadcast these parameters to all the processors - necessary because
      ! the above was only done for the first processor (proc0).
      !**********************************************************************
      subroutine broadcast_parameters  

        use mp, only: broadcast

        call broadcast(stella_debug)
        call broadcast(ffs_solve_debug)
        call broadcast(fields_debug)
        call broadcast(fields_fluxtube_debug)
        call broadcast(fields_electromagnetic_debug)
        call broadcast(fields_ffs_debug)
        call broadcast(implicit_solve_debug)
        call broadcast(mirror_terms_debug)
        call broadcast(neoclassical_terms_debug)
        call broadcast(parallel_streaming_debug)
        call broadcast(response_matrix_debug)
        call broadcast(time_advance_debug)
        call broadcast(extended_grid_debug)
        call broadcast(diagnostics_debug)
        call broadcast(diagnostics_parameters)
        call broadcast(diagnostics_fluxes_fluxtube_debug)
        call broadcast(diagnostics_omega_debug) 
        call broadcast(dist_fn_debug)
        call broadcast(gyro_averages_debug)
        call broadcast(fluxes_debug) 
        call broadcast(geometry_debug)
        call broadcast(const_alpha_geo) 

      end subroutine broadcast_parameters

  end subroutine read_debug_flags
   
end module debug_flags
