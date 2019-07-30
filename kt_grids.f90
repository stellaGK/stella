! Set up ranges of kx and ky for linear runs.
module kt_grids_range

  implicit none

  public :: init_kt_grids_range, range_get_sizes, range_get_grids

  private

  integer :: naky, nakx, ntheta0
  real :: aky_min, aky_max
  real :: theta0_min, theta0_max
  real :: akx_min, akx_max
  namelist /kt_grids_range_parameters/ naky, nakx, ntheta0, &
       aky_min, aky_max, theta0_min, theta0_max, akx_min, akx_max

contains

  subroutine init_kt_grids_range
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    integer :: in_file
    logical :: exist

    naky = 1 ; nakx = 1
    aky_min = 0.0 ; aky_max = 0.0
    akx_min = 0.0 ; akx_max = -1.0
    ! set these to be nonsense values
    ! so can check later if they've been set
    theta0_min = 0.0 ; theta0_max = -1.0

    in_file = input_unit_exist ("kt_grids_range_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_range_parameters)

    ! nakx is fundamental grid parameter
    ! ntheta0 is same as nakx for ballooning calculations
    ntheta0 = nakx
  end subroutine init_kt_grids_range

  subroutine range_get_sizes (naky_out, nakx_out, ntheta0_out, nx, ny)

    implicit none

    integer, intent (out) :: naky_out, nakx_out, ntheta0_out, nx, ny

    naky_out = naky
    nakx_out = nakx
    ntheta0_out = ntheta0
    nx = 0
    ny = 0

  end subroutine range_get_sizes

  subroutine range_get_grids (aky, theta0, akx)
    use mp, only: mp_abort
    use zgrid, only: shat_zero
    use stella_geometry, only: geo_surf

    implicit none

    real, dimension (:), intent (out) :: akx, aky
    real, dimension (:,:), intent (out) :: theta0

    real :: dkx, dky, dtheta0, zero
    integer :: i, j

    if (size(aky) /= naky) call mp_abort ('range_get_grids: size(aky) /= naky')
    if (size(akx) /= nakx) call mp_abort ('range_get_grids: size(akx) /= nakx')

    ! NB: we are assuming here that all ky are positive
    ! when running in range mode
    dky = 0.0
    if (naky > 1) dky = (aky_max - aky_min)/real(naky - 1)
    aky = (/ (aky_min + dky*real(i), i = 0,naky-1) /)

    ! set default akx and theta0 to 0
    akx = 0.0 ; theta0=0.0

    zero = 100.*epsilon(0.)

    ! if theta0_min and theta0_max have been specified,
    ! use them to determine akx_min and akx_max
    if (theta0_max > theta0_min-zero) then
       if (geo_surf%shat > epsilon(0.)) then
          akx_min = theta0_min * geo_surf%shat * aky(1)
          akx_max = theta0_max * geo_surf%shat * aky(1)
       else
          akx_min = theta0_max * geo_surf%shat * aky(1)
          akx_max = theta0_min * geo_surf%shat * aky(1)
       end if
    end if

    ! shat_zero is minimum shat value below which periodic BC is enforced
    if (abs(geo_surf%shat) > shat_zero) then  ! ie assumes boundary_option .eq. 'linked'
       ! if akx_min and akx_max specified in input
       ! instead of theta0_min and theta0_max,
       ! use them to get theta0_min and theta0_max
       if (theta0_min > theta0_max+zero .and. abs(aky(1)) > epsilon(0.)) then
          theta0_min = akx_min/(geo_surf%shat*aky(1))
          theta0_max = akx_max/(geo_surf%shat*aky(1))
          dtheta0 = 0.0
          if (ntheta0 > 1) dtheta0 = (theta0_max - theta0_min)/real(ntheta0 - 1)
          
          do j = 1, naky
             theta0(j,:) &
                  = (/ (theta0_min + dtheta0*real(i), i=0,ntheta0-1) /)
          end do
          akx = theta0(1,:) * geo_surf%shat * aky(1)
       else if (akx_max > akx_min-zero) then
          dkx = 0.0
          if (nakx > 1) dkx = (akx_max - akx_min)/real(nakx - 1)
          akx = (/ (akx_min + dkx*real(i), i = 0,nakx-1) /)

          dtheta0 = 0.0
          if (ntheta0 > 1) dtheta0 = (theta0_max - theta0_min)/real(ntheta0 - 1)

          if (geo_surf%shat > epsilon(0.)) then
             do j = 1, naky
                theta0(j,:) &
                     = (/ (theta0_min + dtheta0*real(i), i=0,ntheta0-1) /)
             end do
          else
             do j = 1, naky
                theta0(j,:) &
                     = (/ (theta0_min + dtheta0*real(i), i=ntheta0-1,0,-1) /)
             end do
          end if
       else
          call mp_abort ('ky=0 is inconsistent with kx_min different from kx_max. aborting.')
       end if
       
    else
       ! here assume boundary_option .eq. 'periodic'
       ! used for periodic finite kx ballooning space runs with shat=0
       dkx = 0.0
       if (nakx > 1) dkx = (akx_max - akx_min)/real(nakx - 1)
       akx = (/ (akx_min + dkx*real(i), i = 0,nakx-1) /)
    endif
    
  end subroutine range_get_grids
  
end module kt_grids_range

! Set the perpendicular box size and resolution for linear or nonlinear runs.
module kt_grids_box

  implicit none

  public :: init_kt_grids_box, box_get_sizes, box_get_grids
  public :: x0, y0, jtwist !RN> Caution: these are not broadcasted!

  private

  real :: ly, y0, x0, rtwist
  integer :: naky_box, nakx_box, ntheta0_box, nx_box, ny_box
  integer :: jtwist

contains

  subroutine init_kt_grids_box
    use zgrid, only: init_zgrid
    use file_utils, only: input_unit, input_unit_exist
    use constants
    use stella_geometry, only: twist_and_shift_geo_fac
!    use stella_geometry, only: geo_surf

    implicit none

    integer :: naky, nakx, nx, ny
    integer :: in_file
    logical :: exist
    namelist /kt_grids_box_parameters/ naky, nakx, nx, ny, &
         jtwist, rtwist, x0, y0, ly

    call init_zgrid

    nakx = 0
    naky = 0
    ly = 0.0 
    y0 = 10.0
    x0 = 0.0
    nx = 0
    ny = 0

    ! for stellarators, twist_and_shift_geo_fac = p/q (see stella_JCP)
!    jtwist = max(int(2.0*pi*geo_surf%shat*twist_and_shift_geo_fac + 0.5),1)
!    jtwist = max(1,geo_surf%shat*(gds21(1,-nzgrid)/gds22(1,nzgrid)-gds21(1,nzgrid)/gds22(1,nzgrid))+0.5)
    jtwist = max(1,int(abs(twist_and_shift_geo_fac)+0.5))
    rtwist = 0.0

    in_file = input_unit_exist("kt_grids_box_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_box_parameters)

    if (y0 < 0) y0 = -1./y0

    if (abs(ly) < epsilon(0.)) ly = 2.0*pi*y0
    if (nakx == 0) nakx = 2*((nx-1)/3) + 1
    if (naky == 0) naky = (ny-1)/3 + 1
    if (abs(rtwist) < epsilon(0.)) rtwist = real(jtwist)
    
    naky_box = naky
    nakx_box = nakx
    ntheta0_box = nakx
    nx_box = nx
    ny_box = ny

  end subroutine init_kt_grids_box

  subroutine box_get_sizes (naky, nakx, ntheta0, nx, ny)
    implicit none
    integer, intent (out) :: naky, nakx, ntheta0, nx, ny
    naky = naky_box
    nakx = nakx_box
    ntheta0 = ntheta0_box
    nx = nx_box
    ny = ny_box
  end subroutine box_get_sizes

  subroutine box_get_grids (aky, theta0, akx)

    use zgrid, only: shat_zero
    use constants
    use stella_geometry, only: geo_surf, twist_and_shift_geo_fac

    implicit none

    real, dimension (:), intent (out) :: akx, aky
    real, dimension (:,:), intent (out) :: theta0

    real :: dkx, dky, ratio
    integer :: i, naky, ntheta0, nakx
    integer :: ikx_max

    naky = naky_box
    nakx = nakx_box
    ntheta0 = ntheta0_box

    dky = 2.0*pi/ly

    ! non-quantized b/c assumed to be periodic instead 
    ! of linked boundary conditions
    if (abs(geo_surf%shat) <=  shat_zero) then   
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
          ! twist_and_shift_geo_fac = p/q (see stella_JCP)
!          dkx = dky * 2.0*pi*abs(geo_surf%shat)*twist_and_shift_geo_fac/real(jtwist)
          dkx = dky * abs(twist_and_shift_geo_fac)/real(jtwist)
       else
          dkx = dky
       end if
    endif

    ! ky goes from zero to ky_max
    do i = 1, naky
       aky(i) = real(i-1)*dky
    end do

    ikx_max = nakx/2+1

    ! kx goes from zero to kx_max ...
    do i = 1, ikx_max
       akx(i) = real(i-1)*dkx
    end do
    ! and then from -kx_max to -|kx_min|
    do i = ikx_max+1, nakx
       akx(i) = real(i-nakx-1)*dkx
    end do

    ! set theta0=0 for ky=0
    theta0(1,:) = 0.0
    if (abs(geo_surf%shat) > shat_zero) then
       do i = 1, nakx
          ! theta0 = kx/ky/shat
          theta0(2:,i) = akx(i)/(aky(2:)*geo_surf%shat)
       end do
    else
       do i = 1, nakx
          ! if shat=0, theta0 is meaningless, so be careful
          theta0(2:,i) = - akx(i)/aky(2:)
       end do
    end if

  end subroutine box_get_grids

end module kt_grids_box

! Set up the perpendicular wavenumbers by calling the appropriate sub-modules. 
module kt_grids

  use kt_grids_box, only: jtwist, y0

  implicit none

  public :: init_kt_grids, box, finish_kt_grids
  public :: aky, theta0, akx
  public :: naky, nakx, ntheta0, nx, ny, reality
  public :: jtwist_out, ikx_twist_shift
  public :: gridopt_switch, grid_option
  public :: gridopt_range, gridopt_box
  public :: lx, ly
  public :: full_flux_surface, ny_ffs
  public :: ikx_max
  public :: zonal_mode
  public :: swap_kxky, swap_kxky_back

  private

  real :: lx, ly
  real, dimension (:,:), allocatable :: theta0
  real, dimension (:), allocatable :: aky, akx
  integer :: naky, nakx, ntheta0, nx, ny
  integer :: jtwist_out, ikx_twist_shift
  integer :: ikx_max
  character(20) :: grid_option
  integer :: ny_ffs = 1
  logical :: full_flux_surface
  logical, dimension (:), allocatable :: zonal_mode

  namelist /kt_grids_knobs/ grid_option, full_flux_surface

  ! internal variables
  integer :: gridopt_switch
  integer, parameter :: gridopt_range = 1, gridopt_box = 2
  logical :: reality = .false.
  logical :: box = .false.
  logical :: initialized = .false.
  logical :: nml_exist

contains

  subroutine init_kt_grids

    use zgrid, only: init_zgrid
    use zgrid, only: shat_zero
    use mp, only: proc0, broadcast
    use constants, only: pi
    use stella_geometry, only: geo_surf, twist_and_shift_geo_fac

    implicit none

    integer :: ik

    if (initialized) return
    initialized = .true.

    call init_zgrid

    if (proc0) then
       call read_parameters
       call get_sizes
       jtwist_out = jtwist
       ! may be better to calculate ikx_twist_shift
       ! in case of stellarators
       ikx_twist_shift = -jtwist*int(sign(1.0,geo_surf%shat))
       ! get the ikx index corresponding to kx_max
       ikx_max = nakx/2+1
    end if

    call broadcast (reality)
    call broadcast (box)
    call broadcast (naky)
    call broadcast (nakx)
    call broadcast (ntheta0)
    call broadcast (ny)
    call broadcast (nx)
    call broadcast (full_flux_surface)
    call broadcast (ikx_max)
    call allocate_arrays

    if (proc0) call get_grids
    call broadcast (aky)
    call broadcast (akx)
    call broadcast (jtwist_out)
    call broadcast (ikx_twist_shift)
    do ik = 1, naky
       call broadcast (theta0(ik,:))
    end do

    ly = 2.*pi*y0
    if (abs(geo_surf%shat) > shat_zero) then
!       lx = y0*jtwist/(geo_surf%shat*twist_and_shift_geo_fac)
       lx = ly*jtwist/abs(twist_and_shift_geo_fac)
    else
       lx = ly
    end if

    ikx_max = nakx/2+1

    ! determine if iky corresponds to zonal mode
    if (.not.allocated(zonal_mode)) allocate (zonal_mode(naky))
    zonal_mode = .false.
    if (abs(aky(1)) < epsilon(0.)) zonal_mode(1) = .true.

    if (full_flux_surface) ny_ffs = ny

  end subroutine init_kt_grids

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    implicit none
    type (text_option), dimension (4), parameter :: gridopts = &
         (/ text_option('default', gridopt_range), &
            text_option('range', gridopt_range), &
            text_option('box', gridopt_box), &
            text_option('nonlinear', gridopt_box) /)
    integer :: ierr, in_file

    full_flux_surface = .false.
    grid_option = 'default'

    in_file = input_unit_exist ("kt_grids_knobs", nml_exist)
    if (nml_exist) read (unit=in_file, nml=kt_grids_knobs)

    ierr = error_unit()
    call get_option_value (grid_option, gridopts, gridopt_switch, &
         ierr, "grid_option in kt_grids_knobs")

  end subroutine read_parameters

  subroutine allocate_arrays
    implicit none
    allocate (akx(nakx))
    allocate (aky(naky))
    allocate (theta0(naky,nakx))
  end subroutine allocate_arrays

  subroutine get_sizes

    use kt_grids_range, only: init_kt_grids_range, range_get_sizes
    use kt_grids_box, only: init_kt_grids_box, box_get_sizes

    implicit none

    select case (gridopt_switch)
    case (gridopt_range)
       call init_kt_grids_range
       call range_get_sizes (naky, nakx, ntheta0, nx, ny)
    case (gridopt_box)
       call init_kt_grids_box
       call box_get_sizes (naky, nakx, ntheta0, nx, ny)
       reality = .true.
       box = .true.
    end select

  end subroutine get_sizes

  subroutine get_grids

    use kt_grids_range, only: range_get_grids
    use kt_grids_box, only: box_get_grids

    implicit none

    select case (gridopt_switch)
    case (gridopt_range)
       call range_get_grids (aky, theta0, akx)
    case (gridopt_box)
       call box_get_grids (aky, theta0, akx)
    end select

  end subroutine get_grids

  ! take an array with ky >= 0 and all kx
  ! and uses reality condition to return array
  ! with kx >= 0 and all ky
  subroutine swap_kxky (gin, gout)

    implicit none

    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    integer :: ikx, ikxneg
    integer :: iky, ikyneg

    ! first set arrays equal for ky >= 0 and kx >= 0
    gout(:naky,:) = gin(:,:ikx_max)
    ! next fill in ky < 0, kx >= 0 elements of array using reality
    ikx = 1
    ikxneg = ikx
    do iky = naky+1, 2*naky-1
       ! this is the ky index corresponding to +ky in original array
       ikyneg = 2*naky-iky+1
       gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
    end do
    do ikx = 2, ikx_max
       ikxneg = nakx-ikx+2
       do iky = naky+1, 2*naky-1
          ! this is the ky index corresponding to +ky in original array
          ikyneg = 2*naky-iky+1
          gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
       end do
    end do

  end subroutine swap_kxky

  ! take an array with ky >= 0 and all kx
  ! and uses reality condition to return array
  ! with kx >= 0 and all ky
  subroutine swap_kxky_back (gin, gout)

    implicit none

    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    integer :: ikx, ikxneg
    integer :: iky, ikyneg

    ! first set arrays equal for ky >= 0 and kx >= 0
    gout(:,:ikx_max) = gin(:naky,:)
    ! next fill in kx < 0, ky >= 0 elements of array using reality
    do ikx = ikx_max+1, nakx
       ikxneg = nakx-ikx+2
       iky = 1
       ikyneg = iky
       gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
       do iky = 2, naky
          ikyneg = 2*naky-iky+1
          gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
       end do
    end do

  end subroutine swap_kxky_back

  subroutine finish_kt_grids

    implicit none

    if (allocated(aky)) deallocate (akx, aky, theta0)

    reality = .false. ; box = .false.
    initialized = .false.

  end subroutine finish_kt_grids

end module kt_grids

