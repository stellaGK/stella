!###############################################################################
!######################## READ PARAMETES FOR KXKY GRIDS ########################
!###############################################################################
! Namelist: &parameters_kxky_grids
! These flags will allow you to toggle the algorithm choices in stella.
!###############################################################################

module parameters_kxky_grids

  public :: read_kxky_grid_parameters
  public :: gridopt_switch, gridopt_range, gridopt_box
  
  private

  integer :: gridopt_switch
  integer, parameter :: gridopt_range = 1, gridopt_box = 2

  integer :: kyspacing_option_switch
  integer, parameter :: kyspacing_linear = 1, kyspacing_exponential = 2

  logical :: initialised 

contains
  
  subroutine read_kxky_grid_parameters

    use mp, only: proc0, mp_abort
    use text_options, only: text_option, get_option_value
    use file_utils, only: input_unit, error_unit, input_unit_exist

    use parameters_kxky_grids_box, only : read_kxky_grids_box
    use parameters_kxky_grids_range, only : read_kxky_grids_range
    implicit none

    logical :: error = .false.

    if (initialised) return

    if (proc0) then
      call read_grid_option
      select case (gridopt_switch)
      case (gridopt_range)
         call read_kxky_grids_range
      case (gridopt_box)
         call read_kxky_grids_box
      end select
    end if

    call broadcast_parameters
    initialised = .true.

  contains

    !**********************************************************************
    !                       READ GRID OPTION FOR KXK                      !
    !**********************************************************************
    ! Read which option to select for the kxky grid layouts
    !**********************************************************************
    subroutine read_grid_option
      
      
      use file_utils, only: input_unit, error_unit, input_unit_exist
      use text_options, only: text_option, get_option_value
      
      implicit none
      
      type(text_option), dimension(5), parameter :: gridopts = &
           (/text_option('default', gridopt_range), &
           text_option('range', gridopt_range), &
           text_option('box', gridopt_box), &
           text_option('annulus', gridopt_box), &
           text_option('nonlinear', gridopt_box)/)
      
      integer :: ierr, in_file
      logical :: nml_exist
      
      character(20) :: grid_option
      
      namelist /kt_grids_knobs/ grid_option
      
      grid_option = 'default'
      
      in_file = input_unit_exist("kt_grids_knobs", nml_exist)
      if (nml_exist) read (unit=in_file, nml=kt_grids_knobs)
      
      ierr = error_unit()
      call get_option_value(grid_option, gridopts, gridopt_switch, &
           ierr, "grid_option in kt_grids_knobs")
      
    end subroutine read_grid_option
    
    subroutine broadcast_parameters
      
      use arrays_kxky
      use mp, only: broadcast
      
      implicit none

      call broadcast(gridopt_switch)
      call broadcast(centered_in_rho)
      call broadcast(periodic_variation)
      call broadcast(naky)
      call broadcast(naky_all)
      call broadcast(nakx)
      call broadcast(ny)
      call broadcast(nx)
      call broadcast(nalpha)
      call broadcast(reality)
      call broadcast(jtwist)
      call broadcast(jtwistfac)
      call broadcast(x0)
      call broadcast(y0)
      call broadcast(aky_min)
      call broadcast(aky_max)
      call broadcast(akx_min)
      call broadcast(akx_max)
      call broadcast(theta0_min)
      call broadcast(theta0_max)
      call broadcast(randomize_phase_shift)
      call broadcast(phase_shift_angle)
      call broadcast(kyspacing_option_switch)
      
    end subroutine broadcast_parameters
    
  end subroutine read_kxky_grid_parameters
  
end module parameters_kxky_grids
