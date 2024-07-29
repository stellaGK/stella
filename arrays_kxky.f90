module arrays_kxky

  implicit none

  character(20) :: grid_option
  
  real, dimension(:), allocatable :: aky, akx
  real, dimension(:), allocatable :: aky_all, aky_all_ordered
  real, dimension(:, :), allocatable :: theta0, zed0
  real, dimension(:), allocatable :: x, x_d, y
  integer :: naky, nakx, nx, ny
  integer :: nalpha, naky_all, ikx_max
  logical :: reality = .false. 
  real :: dx, dy, dkx, dky, dx_d
  real :: phase_shift_angle
  logical, dimension(:), allocatable :: zonal_mode
  real, dimension(:), allocatable :: rho, rho_d, rho_clamped, rho_d_clamped
  integer :: jtwist
  real :: jtwistfac
  real :: ikx_twist_shift

  logical :: centered_in_rho, periodic_variation, randomize_phase_shift
  logical :: box
  
  integer :: boundary_size, copy_size, krook_size

  !! RANGE 
  real :: aky_min, aky_max
  real :: akx_min, akx_max
  real :: theta0_min, theta0_max
  integer :: kyspacing_option_switch
  integer, parameter :: kyspacing_linear = 1, kyspacing_exponential = 2
  
  !! BOX
  real :: x0, y0

  !! RADIAL 
  complex, dimension(:, :), allocatable :: g0x

  
end module arrays_kxky
