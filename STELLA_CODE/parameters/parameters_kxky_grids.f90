!###############################################################################
!######################## READ PARAMETES FOR KXKY GRIDS ########################
!###############################################################################
! Namelist: &parameters_kxky_grids
! These flags will allow you to toggle the algorithm choices in stella.
!###############################################################################

module parameters_kxky_grids

  public :: read_kxky_grid_parameters

  !> Public parameters

  !> Read grid type
  public :: gridopt_switch, gridopt_range, gridopt_box

  !> For Box/Range 
  public :: naky, nakx
  public :: nx, ny
  public :: nalpha, naky_all, ikx_max
  public :: reality 
  public :: phase_shift_angle
  public :: jtwist, jtwistfac, ikx_twist_shift
  public :: centered_in_rho
  public :: periodic_variation, randomize_phase_shift
  public :: aky_min, aky_max
  public :: akx_min, akx_max
  public :: theta0_min, theta0_max
  public :: x0, y0
  
  !> ky spaction options
  public :: kyspacing_option_switch, kyspacing_linear, kyspacing_exponential
  
  private

  !> Internal 
  integer :: gridopt_switch
  integer, parameter :: gridopt_range = 1, gridopt_box = 2

  integer :: naky, nakx, nx, ny
  integer :: nalpha, naky_all, ikx_max
  logical :: reality = .false. 
  real :: phase_shift_angle
  integer :: jtwist
  real :: jtwistfac
  real :: ikx_twist_shift

  logical :: centered_in_rho, periodic_variation, randomize_phase_shift

  !> For Range
  real :: aky_min, aky_max
  real :: akx_min, akx_max
  real :: theta0_min, theta0_max
  integer :: kyspacing_option_switch
!!  character(20) :: kyspacing_option = 'default'
  integer, parameter :: kyspacing_linear = 1, kyspacing_exponential = 2

  !> For Box
  real :: x0, y0
  
  logical :: initialised 

contains
  
  subroutine read_kxky_grid_parameters

    use mp, only: proc0, mp_abort
    use text_options, only: text_option, get_option_value
    use file_utils, only: input_unit, error_unit, input_unit_exist

    use parameters_kxky_grids_box, only : read_kxky_grids_box
    use parameters_kxky_grids_range, only : read_kxky_grids_range
    
    implicit none

    if (initialised) return

    if (proc0) then
      call read_grid_option
      select case (gridopt_switch)
      case (gridopt_range)
         call read_kxky_grids_range (nalpha, naky, nakx, aky_min, aky_max, & 
              akx_min, akx_max, theta0_min, theta0_max, &
              kyspacing_option_switch, phase_shift_angle, ikx_max, naky_all)
      case (gridopt_box)
         call read_kxky_grids_box (nx, ny, ikx_max, naky_all, naky, nakx, nalpha, &
              x0, y0, jtwist, jtwistfac, phase_shift_angle, &
              centered_in_rho, randomize_phase_shift, periodic_variation, reality)
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
      
      namelist /kxky_grids_knobs/ grid_option
      
      grid_option = 'default'
      
      in_file = input_unit_exist("kxky_grids_knobs", nml_exist)
      if (nml_exist) read (unit=in_file, nml=kxky_grids_knobs)
      
      ierr = error_unit()
      call get_option_value(grid_option, gridopts, gridopt_switch, &
           ierr, "grid_option in kxky_grids_knobs")
      
    end subroutine read_grid_option
    
    subroutine broadcast_parameters
      
      use mp, only: broadcast
      
      implicit none

      call broadcast(gridopt_switch)
      call broadcast(naky)
      call broadcast(nakx)
      call broadcast(ny)
      call broadcast(nx)
      call broadcast(nalpha)
      call broadcast(naky_all)
      call broadcast(ikx_max)
      call broadcast(reality)
      call broadcast(phase_shift_angle)
      call broadcast(jtwist)
      call broadcast(jtwistfac)
      call broadcast(ikx_twist_shift)
      call broadcast(centered_in_rho)
      call broadcast(periodic_variation)
      call broadcast(randomize_phase_shift)
      call broadcast(aky_min)
      call broadcast(aky_max)
      call broadcast(akx_min)
      call broadcast(akx_max)
      call broadcast(theta0_min)
      call broadcast(theta0_max)
      call broadcast(kyspacing_option_switch)
      call broadcast(x0)
      call broadcast(y0)  

    end subroutine broadcast_parameters
    
  end subroutine read_kxky_grid_parameters
  
end module parameters_kxky_grids
