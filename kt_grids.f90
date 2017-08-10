module kt_grids_single
! <doc> Set up values of kx and ky for linear runs that use a single k_perp mode.
! </doc>

  implicit none

  public :: init_kt_grids_single, single_get_sizes, single_get_grids

  private

  real :: akx, aky, theta0

contains

  subroutine init_kt_grids_single
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: in_file
    logical :: exist
    namelist /kt_grids_single_parameters/ aky, theta0, akx

    aky = 0.4   ;  theta0 = 0.0   ;   akx = 0.0

    in_file = input_unit_exist ("kt_grids_single_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_single_parameters)

  end subroutine init_kt_grids_single

  subroutine single_get_sizes (naky, ntheta0, nx, ny)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny

    naky = 1  ;  ntheta0 = 1
    nx = 0    ;  ny = 0

  end subroutine single_get_sizes

  subroutine single_get_grids (aky_out, theta0_out, akx_out)
    implicit none
    real, dimension (:), intent (out) :: aky_out
    real, dimension (:,:), intent (out) :: theta0_out
    real, dimension (:), intent (out) :: akx_out

    aky_out = aky
    theta0_out = theta0
    akx_out = akx

  end subroutine single_get_grids

end module kt_grids_single

module kt_grids_range
! <doc> Set up ranges of kx and ky for linear runs.
! </doc>
  implicit none

  public :: init_kt_grids_range, range_get_sizes, range_get_grids

  private

  integer :: naky, ntheta0
  real :: aky_min, aky_max, theta0_min, theta0_max
  real :: akx_min, akx_max
  namelist /kt_grids_range_parameters/ naky, ntheta0, &
       aky_min, aky_max, theta0_min, theta0_max, akx_min, akx_max

contains

  subroutine init_kt_grids_range
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: in_file
    logical :: exist

    naky = 1          ;  ntheta0 = 1  
    aky_min = 0.0     ;  aky_max = 0.0
    theta0_min = 0.0  ;  theta0_max = 0.0
    akx_min = 0.0     ;  akx_max = 0.0

    in_file = input_unit_exist ("kt_grids_range_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_range_parameters)

  end subroutine init_kt_grids_range

  subroutine range_get_sizes (naky_x, ntheta0_x, nx, ny)
    implicit none
    integer, intent (out) :: naky_x, ntheta0_x, nx, ny

    naky_x = naky  ;  ntheta0_x = ntheta0
    nx = 0         ;  ny = 0

  end subroutine range_get_sizes

! BD: Could add some logic here to set theta0 if akx is given?  When do we need what?

  subroutine range_get_grids (aky, theta0, akx)
    use theta_grid, only: shat
    implicit none
    real, dimension (:), intent (out) :: akx, aky
    real, dimension (:,:), intent (out) :: theta0

    real :: dkx, dky, dtheta0
    integer :: i, j

    if ( size(aky) /= naky) then
       write(6,*) 'range_get_grids: size(aky) /= naky'       ;  stop
    endif

    if ( size(akx) /= ntheta0) then
       write(6,*) 'range_get_grids: size(akx) /= ntheta0'    ;  stop
    endif

    dky = 0.0
    if (naky > 1) dky = (aky_max - aky_min)/real(naky - 1)
    aky = (/ (aky_min + dky*real(i), i = 0,naky-1) /)

! set default theta0 to 0
    theta0=0.0

!
! BD: Assumption here differs from convention that abs(shat) <= 1.e-5 triggers periodic bc
!
    if (abs(shat) > epsilon(0.)) then  ! ie assumes boundary_option .eq. 'linked'
       dtheta0 = 0.0
       if (ntheta0 > 1) dtheta0 = (theta0_max - theta0_min)/real(ntheta0 - 1)

       do j = 1, naky
          theta0(j,:) &
               = (/ (theta0_min + dtheta0*real(i), i=0,ntheta0-1) /)
       end do
       akx = theta0(1,:) * shat * aky(1)

    else

!CMR, 22/9/2010:  ie here assume boundary_option .eq. 'periodic'
!new code for periodic finite kx ballooning space runs with shat=0

       dkx = 0.0
       if (ntheta0 > 1) dkx = (akx_max - akx_min)/real(ntheta0 - 1)
       akx = (/ (akx_min + dkx*real(i), i = 0,ntheta0-1) /)

    endif

  end subroutine range_get_grids

end module kt_grids_range

module kt_grids_specified
! <doc> Set up sets of (kx, ky) values for linear runs.
! </doc>
  implicit none

  public :: init_kt_grids_specified, specified_get_sizes, specified_get_grids

  private

  integer :: naky, ntheta0, nx, ny
  namelist /kt_grids_specified_parameters/ naky, ntheta0, nx, ny

contains

  subroutine init_kt_grids_specified
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: in_file
    logical :: exist

    naky = 1    ;  ntheta0 = 1
    nx = 0      ;  ny = 0

    in_file = input_unit_exist("kt_grids_specified_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_specified_parameters)

  end subroutine init_kt_grids_specified

  subroutine specified_get_sizes (naky_x, ntheta0_x, nx_x, ny_x)
    implicit none
    integer, intent (out) :: naky_x, ntheta0_x, nx_x, ny_x

    naky_x = naky  ;   ntheta0_x = ntheta0
    nx_x = nx      ;   ny_x = ny

  end subroutine specified_get_sizes

  subroutine specified_get_grids (aky, theta0, akx)
    implicit none
    real, dimension (:), intent (out) :: akx, aky
    real, dimension (:,:), intent (out) :: theta0
    real :: aky_dummy, theta0_dummy, akx_dummy
    integer :: i, naky, ntheta0

    naky = size(aky)  ;  ntheta0 = size(akx)

    do i = 1, max(naky,ntheta0)
       call read_element (i, aky_dummy, theta0_dummy, akx_dummy)
       if (i <= naky) aky(i) = aky_dummy
       if (i <= ntheta0) theta0(:,i) = theta0_dummy
       if (i <= ntheta0) akx(i) = akx_dummy
    end do

  end subroutine specified_get_grids

  subroutine read_element (i, aky_dummy, theta0_dummy, akx_dummy)
    use file_utils, only: get_indexed_namelist_unit
    implicit none
    integer, intent (in) :: i
    real, intent (out) :: aky_dummy, theta0_dummy, akx_dummy
    real :: akx, aky, theta0
    integer :: unit

    namelist /kt_grids_specified_element/ aky, theta0, akx

    aky = 0.4
    theta0 = 0.0
    akx = 0.0
    call get_indexed_namelist_unit (unit, "kt_grids_specified_element", i)
    read (unit=unit, nml=kt_grids_specified_element)
    close (unit)
    aky_dummy = aky
    theta0_dummy = theta0
    akx_dummy = akx
  end subroutine read_element

end module kt_grids_specified

module kt_grids_box
! <doc> Set the perpendicular box size and resolution for linear or nonlinear runs.
! </doc>
  implicit none

  public :: init_kt_grids_box, box_get_sizes, box_get_grids
  public :: x0, y0, jtwist !RN> Caution: these are not broadcasted!

  private

  real :: ly, y0, x0, rtwist
  integer :: naky_private, ntheta0_private, nx_private, ny_private
  integer :: nkpolar_private
  integer :: jtwist

contains

  subroutine init_kt_grids_box
    use theta_grid, only: init_theta_grid, shat
    use file_utils, only: input_unit, input_unit_exist
    use constants
    implicit none
    integer :: naky, ntheta0, nx, ny, nkpolar
    integer :: in_file
    logical :: exist
    namelist /kt_grids_box_parameters/ naky, ntheta0, ly, nx, ny, jtwist, &
         y0, rtwist, x0, nkpolar

    call init_theta_grid

    nkpolar = 0   ;   naky = 0    ;  ntheta0 = 0
    ly = 0.0      ;   y0 = 2.0    ;  x0 = 0.
    nx = 0        ;   ny = 0

    jtwist = max(int(2.0*pi*shat + 0.5),1)  ! default jtwist -- MAB
    rtwist = 0.0

    in_file = input_unit_exist("kt_grids_box_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_box_parameters)

    if (y0 < 0) y0 = -1./y0

    if (abs(ly) < epsilon(0.)) ly = 2.0*pi*y0
    if (naky == 0) naky = (ny-1)/3 + 1
    if (ntheta0 == 0) ntheta0 = 2*((nx-1)/3) + 1
    if (abs(rtwist) < epsilon(0.)) rtwist = real(jtwist)
    if (nkpolar == 0) nkpolar = int(real(naky-1.)*sqrt(2.))
    
    nkpolar_private = nkpolar
    naky_private = naky
    ntheta0_private = ntheta0
    nx_private = nx
    ny_private = ny

  end subroutine init_kt_grids_box

  subroutine box_get_sizes (naky, ntheta0, nx, ny, nkpolar)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny, nkpolar
    naky = naky_private
    ntheta0 = ntheta0_private
    nx = nx_private
    ny = ny_private
    nkpolar = nkpolar_private
  end subroutine box_get_sizes

  subroutine box_get_grids (aky, theta0, akx, ikx, iky)
    use theta_grid, only: shat
    use constants
    implicit none
    real, dimension (:), intent (out) :: akx, aky
    real, dimension (:,:), intent (out) :: theta0
    integer, dimension (:), intent (out) :: ikx, iky

    real :: dkx, dky, ratio
    integer :: i, naky, ntheta0

    naky = size(aky)    
    ntheta0 = size(akx)

    dky = 2.0*pi/ly

    if(abs(shat) <=  1.e-5) then   ! non-quantized b/c assumed to be periodic instead linked boundary conditions       

       if (abs(x0) < epsilon(0.)) then          
          
          if (rtwist > 0) then 
             ratio = rtwist
          else 
             ratio = 1. / abs(rtwist)
          end if
          
          dkx = dky / ratio
          
       else

          if (x0 > 0.) then
             dkx = 1./x0
          else
             dkx = -x0
          end if

       end if

    else
       if (jtwist /= 0) then
          dkx = dky * 2.0*pi*abs(shat)/real(jtwist)
       else
          dkx = dky
       end if
    endif

    do i = 1, naky
       iky(i) = i-1
       aky(i) = real(i-1)*dky
    end do

    do i = 1, (ntheta0+1)/2
       ikx(i) = i-1
       akx(i) = real(i-1)*dkx
    end do

    do i = (ntheta0+1)/2+1, ntheta0
       ikx(i) = i-ntheta0-1
       akx(i) = real(i-ntheta0-1)*dkx
    end do

    if (abs(shat) > epsilon(0.)) then
       do i = 1, ntheta0
          theta0(1,i) = 0.0
          theta0(2:,i) = akx(i)/(aky(2:)*shat)
       end do
    else
       do i = 1, ntheta0
          theta0(1,i) = 0.0
          theta0(2:,i) = - akx(i)/aky(2:)   ! meaningless, so be careful
       end do
    end if

  end subroutine box_get_grids

end module kt_grids_box

module kt_grids
!  <doc> Set up the perpendicular wavenumbers by calling the appropriate sub-modules. 
! </doc>
  use kt_grids_box, only: jtwist, y0
  implicit none

  public :: init_kt_grids, box, finish_kt_grids
  public :: aky, theta0, akx, xgrid
  public :: naky, ntheta0, nx, ny, reality
  public :: nkpolar 
  public :: ikx, iky, jtwist_out
  public :: gridopt_switch, grid_option
  public :: gridopt_single, gridopt_range, gridopt_specified, gridopt_box
  public :: lx, ly

  private

  real :: lx, ly
  real, dimension (:,:), allocatable :: theta0
  real, dimension (:), allocatable :: aky, akx, xgrid
  integer, dimension(:), allocatable :: ikx, iky
  integer :: naky, ntheta0, nx, ny, nkpolar, jtwist_out
  character(20) :: grid_option

  namelist /kt_grids_knobs/ grid_option

  ! internal variables
  integer :: gridopt_switch
  integer, parameter :: gridopt_single = 1, gridopt_range = 2, &
       gridopt_specified = 3, gridopt_box = 4
  logical :: reality = .false.
  logical :: box = .false.
  logical :: initialized = .false.
  logical :: nml_exist

contains

  subroutine init_kt_grids

    use geometry, only: rhoc
    use theta_grid, only: init_theta_grid, shat, qval
    use mp, only: proc0, broadcast
    use constants, only: pi

    implicit none

    integer :: ik, it

    if (initialized) return
    initialized = .true.

    call init_theta_grid

    if (proc0) then
       nkpolar = 0   ! will be set to non-zero value only in box case; only used for an MHD diagnostic
       call read_parameters
       call get_sizes
       jtwist_out = jtwist
    end if

    call broadcast (reality)
    call broadcast (box)
    call broadcast (naky)
    call broadcast (nkpolar)
    call broadcast (ntheta0)
    call broadcast (ny)
    call broadcast (nx)
    call allocate_arrays

    if (proc0) call get_grids
    call broadcast (ikx)     ! MR
    call broadcast (aky)
    call broadcast (akx)
    call broadcast (jtwist_out)
    do ik = 1, naky
       call broadcast (theta0(ik,:))
    end do

    ly = 2.*pi*y0
    lx = y0*jtwist/shat

    do it = 1, ntheta0/2
       xgrid(it) = lx*(it-1)/ntheta0
    end do
    do it = ntheta0/2+1, ntheta0
       xgrid(it) = lx*(it-1)/ntheta0 - lx
    end do
    xgrid = xgrid * rhoc/qval

  end subroutine init_kt_grids

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    implicit none
    type (text_option), dimension (6), parameter :: gridopts = &
         (/ text_option('default', gridopt_single), &
            text_option('single', gridopt_single), &
            text_option('range', gridopt_range), &
            text_option('specified', gridopt_specified), &
            text_option('box', gridopt_box), &
            text_option('nonlinear', gridopt_box) /)
    integer :: ierr, in_file

    grid_option = 'default'
    in_file = input_unit_exist ("kt_grids_knobs", nml_exist)
    if (nml_exist) read (unit=in_file, nml=kt_grids_knobs)

    ierr = error_unit()
    call get_option_value (grid_option, gridopts, gridopt_switch, &
         ierr, "grid_option in kt_grids_knobs")

  end subroutine read_parameters

  subroutine allocate_arrays
    implicit none
    allocate (akx(ntheta0))
    allocate (xgrid(ntheta0))
    allocate (aky(naky))
    allocate (theta0(naky,ntheta0))
    allocate (ikx(ntheta0))
    allocate (iky(naky))
  end subroutine allocate_arrays

  subroutine get_sizes
    use kt_grids_single, only: init_kt_grids_single, single_get_sizes
    use kt_grids_range, only: init_kt_grids_range, range_get_sizes
    use kt_grids_specified, only: init_kt_grids_specified, specified_get_sizes
    use kt_grids_box, only: init_kt_grids_box, box_get_sizes
    implicit none
    select case (gridopt_switch)
    case (gridopt_single)
       call init_kt_grids_single
       call single_get_sizes (naky, ntheta0, nx, ny)
    case (gridopt_range)
       call init_kt_grids_range
       call range_get_sizes (naky, ntheta0, nx, ny)
    case (gridopt_specified)
       call init_kt_grids_specified
       call specified_get_sizes (naky, ntheta0, nx, ny)
    case (gridopt_box)
       call init_kt_grids_box
       call box_get_sizes (naky, ntheta0, nx, ny, nkpolar)
       reality = .true.
       box = .true.
    end select
  end subroutine get_sizes

  subroutine get_grids
    use kt_grids_single, only: single_get_grids
    use kt_grids_range, only: range_get_grids
    use kt_grids_specified, only: specified_get_grids
    use kt_grids_box, only: box_get_grids
    implicit none
    select case (gridopt_switch)
    case (gridopt_single)
       call single_get_grids (aky, theta0, akx)
    case (gridopt_range)
       call range_get_grids (aky, theta0, akx)
    case (gridopt_specified)
       call specified_get_grids (aky, theta0, akx)
    case (gridopt_box)
       call box_get_grids (aky, theta0, akx, ikx, iky)
    end select
  end subroutine get_grids

  subroutine finish_kt_grids

    implicit none

    if (allocated(aky)) deallocate (akx, aky, theta0, ikx, iky, xgrid)

    reality = .false. ; box = .false.
    initialized = .false.

  end subroutine finish_kt_grids

end module kt_grids

