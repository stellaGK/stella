!###############################################################################
!######################## READ PARAMETES FOR KXKY GRIDS ########################
!###############################################################################
! Namelist: &parameters_kxky_grids
! These flags will allow you to toggle the algorithm choices in stella.
!###############################################################################

module parameters_kxky_grids

  public :: read_kxky_grid_parameters

  ! Public parameters

  ! Grid type
  public :: grid_option_switch, grid_option_range, grid_option_box

  ! For Box/Range 
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
  
  ! ky spaction options
  public :: kyspacing_option_switch, kyspacing_linear, kyspacing_exponential
  
  private

  integer :: grid_option_switch
  integer, parameter :: grid_option_range = 1, grid_option_box = 2

  integer :: naky, nakx, nx, ny
  integer :: nalpha, naky_all, ikx_max
  logical :: reality = .false. 
  real :: phase_shift_angle
  integer :: jtwist
  real :: jtwistfac
  real :: ikx_twist_shift

  logical :: centered_in_rho, periodic_variation, randomize_phase_shift

  ! For Range
  real :: aky_min, aky_max
  real :: akx_min, akx_max
  real :: theta0_min, theta0_max
  integer :: kyspacing_option_switch
  integer, parameter :: kyspacing_linear = 1, kyspacing_exponential = 2

  ! For Box
  real :: x0, y0
  
  ! Internal variables
  logical :: initialised 

contains
  
  subroutine read_kxky_grid_parameters

    use mp, only: proc0, mp_abort

    use input_file_kxky_grid, only: read_namelist_kxky_grid_option, &
         read_namelist_kxky_grid_box, read_namelist_kxky_grid_range
    
    implicit none

    if (initialised) return

    if (proc0) then
      call read_namelist_kxky_grid_option (grid_option_switch)
      select case (grid_option_switch)
      case (grid_option_range)
         call read_namelist_kxky_grid_range (nalpha, naky, nakx, aky_min, aky_max, & 
              akx_min, akx_max, theta0_min, theta0_max, &
              kyspacing_option_switch, phase_shift_angle, ikx_max, naky_all)
      case (grid_option_box)
         call read_namelist_kxky_grid_box (nx, ny, ikx_max, naky_all, naky, nakx, nalpha, &
              x0, y0, jtwist, jtwistfac, phase_shift_angle, &
              centered_in_rho, randomize_phase_shift, periodic_variation, reality)
      end select
    end if
    
    call broadcast_parameters
    initialised = .true.

  contains
    
    ! Broadcast parameters to all processes
    subroutine broadcast_parameters
      
      use mp, only: broadcast
      
      implicit none

      call broadcast(grid_option_switch)
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
