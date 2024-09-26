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
  
  !> Publuc routines
  public :: read_debug_flags
  !> Public debug flags
  public :: stella_debug
  !> Fields debug flags
  public :: fields_debug
  public :: fields_fluxtube_debug
  public :: fields_electromagnetic_debug
  public :: fields_ffs_debug
  !> Gyrokinetic term debug flags
  public :: time_advance_debug
  public :: implicit_solve_debug
  public :: parallel_streaming_debug
  public :: response_matrix_debug
  public :: mirror_terms_debug
  public :: neoclassical_terms_debug
  !> Grids debug flags 
  public :: extended_grid_debug
  !> Diagnostics debug flags 
  public :: diagnostics_debug
  public :: diagnostics_parameters
  public :: diagnostics_fluxes_fluxtube_debug
  public :: diagnostics_omega_debug
  public :: fluxes_debug
  !> Geometry debug flags
  public :: geometry_debug

  public :: dist_fn_debug
  public :: gyro_averages_debug
  !> Debug flags for full flux surface
  public :: ffs_solve_debug
  !> Debug flug for debugging full flux surface
  !> This will set the geometry of all field lines to be equation to
  !> the geometry on alpha = 0 (i.e. the first field line) 
  public :: const_alpha_geo
  
  private

  logical :: debug_all
  !> Main stella file debug flag
  logical :: stella_debug
  !> Fields debug flags
  logical :: fields_all_debug
  logical :: fields_debug
  logical :: fields_fluxtube_debug
  logical :: fields_electromagnetic_debug
  logical :: fields_ffs_debug
  !> Gyrokinetic term debug flags
  logical :: time_advance_debug
  logical :: implicit_solve_debug
  logical :: parallel_streaming_debug
  logical :: response_matrix_debug
  logical :: mirror_terms_debug
  logical :: neoclassical_terms_debug
  logical :: extended_grid_debug
  !> Diagnostics debug flags
  logical :: diagnostics_all_debug
  logical :: diagnostics_debug
  logical :: diagnostics_parameters
  logical :: diagnostics_fluxes_fluxtube_debug
  logical :: diagnostics_omega_debug
  logical :: fluxes_debug
  !> Geometry debug flags
  logical :: geometry_debug

  logical :: dist_fn_debug
  logical :: gyro_averages_debug
  !> For FFS
  logical :: ffs_solve_debug
  logical :: const_alpha_geo

  logical :: initialised = .false.

contains

  !======================================================================
  !========================= READ DEBUG FLAGS  ==========================
  !======================================================================
  subroutine read_debug_flags

    use mp, only: proc0
    
    implicit none

    namelist /debug_flags/ debug_all, stella_debug, ffs_solve_debug, fields_all_debug, fields_debug, &
        fields_fluxtube_debug, fields_electromagnetic_debug, fields_ffs_debug, & 
        implicit_solve_debug, parallel_streaming_debug, mirror_terms_debug, neoclassical_terms_debug, &
        response_matrix_debug, time_advance_debug, extended_grid_debug, &
        diagnostics_all_debug, diagnostics_parameters, diagnostics_fluxes_fluxtube_debug, &
        diagnostics_omega_debug, diagnostics_debug, dist_fn_debug,&
        gyro_averages_debug, fluxes_debug, geometry_debug,  const_alpha_geo
    
    if (initialised) return
    initialised = .true.
    
    if (proc0) call set_default_parameters
    if (proc0) call read_input_file
    call broadcast_parameters


  contains

      !**********************************************************************
      !                        SET DEFAULT PARAMETERS                       !
      !**********************************************************************
      ! If not specified in the input file these are the default options that 
      ! will be set for all parameters under the namelist 
      ! &debug_flags'.
      ! The default here is that all debug flags are set to .false.
      !**********************************************************************
      subroutine set_default_parameters

        implicit none
        
        stella_debug = .false.
        ffs_solve_debug = .false.

        fields_all_debug = .false. 
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
        
        !> Diagnostircs debug
        diagnostics_all_debug = .false.
        diagnostics_debug = .false.
        diagnostics_parameters = .false.
        diagnostics_omega_debug = .false. 
        diagnostics_fluxes_fluxtube_debug = .false. 
        fluxes_debug = .false.

        geometry_debug = .false.
        
        dist_fn_debug = .false.
        gyro_averages_debug = .false.
        !###################################
        !     FOR THE PURPOSE OF DEBUGGING
        !###################################
        const_alpha_geo = .false. 

      end subroutine set_default_parameters

      !**********************************************************************
      !                         READ INPUT OPTIONS                          !
      !**********************************************************************
      ! Overwrite any default options with those specified in the input file. 
      ! Then change the other parameters consistently.
      !**********************************************************************
      subroutine read_input_file
        
        use file_utils, only: input_unit, error_unit, input_unit_exist
        
        implicit none 

        integer :: in_file
        logical :: nml_exist

        !> Overwrite the default input parameters by those specified in the input file
        !> under the heading '&numerical'  
        in_file = input_unit_exist("debug_flags", nml_exist)
        if (nml_exist) read (unit=in_file, nml=debug_flags)
        
        if(debug_all) then
          stella_debug = .true.
          ffs_solve_debug = .true.
          fields_all_debug = .true. 
          
          time_advance_debug = .true.
          implicit_solve_debug = .true.
          parallel_streaming_debug = .true.
          response_matrix_debug = .true.
          mirror_terms_debug = .true.
          neoclassical_terms_debug = .true.
          extended_grid_debug = .true. 
          !> Set all diagnostics flags to be on
          diagnostics_all_debug = .true.

          dist_fn_debug = .true.
          gyro_averages_debug = .true.
          geometry_debug = .true. 
        end if

        if(diagnostics_all_debug) then
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
      end subroutine read_input_file

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
