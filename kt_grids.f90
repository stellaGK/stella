! Set up the perpendicular wavenumbers by calling the appropriate sub-modules. 
module kt_grids

  implicit none

  public :: init_kt_grids, finish_kt_grids
  public :: read_kt_grids_parameters
  public :: aky, theta0, akx
  public :: naky, nakx, nx, ny, reality
  public :: dx,dy,dkx, dky
  public :: jtwist, ikx_twist_shift, x0, y0, x
  public :: nalpha
  public :: ikx_max, naky_all
  public :: zonal_mode
  public :: swap_kxky, swap_kxky_back
  public :: swap_kxky_ordered, swap_kxky_back_ordered

  private

  interface swap_kxky
     module procedure swap_kxky_real
     module procedure swap_kxky_complex
  end interface

  real, dimension (:,:), allocatable :: theta0
  real, dimension (:), allocatable :: aky, akx
  real, dimension (:), allocatable :: x
  real :: dx, dy, dkx, dky
  integer :: naky, nakx, nx, ny, nalpha
  integer :: jtwist, ikx_twist_shift
  integer :: ikx_max, naky_all
  logical :: reality = .false.
  character(20) :: grid_option
  logical, dimension (:), allocatable :: zonal_mode

  namelist /kt_grids_knobs/ grid_option

  ! internal variables
  integer :: gridopt_switch
  integer, parameter :: gridopt_range = 1, gridopt_box = 2

  real :: aky_min, aky_max
  real :: akx_min, akx_max
  real :: theta0_min, theta0_max
  real :: x0, y0
  logical :: read_kt_grids_initialized = .false.
  logical :: init_kt_grids_initialized = .false.

contains
  
  subroutine read_kt_grids_parameters

    use mp, only: proc0

    implicit none

    if (read_kt_grids_initialized) return
    read_kt_grids_initialized = .true.

    if (proc0) then
       call read_grid_option
       select case (gridopt_switch)
       case (gridopt_range)
          call read_kt_grids_range
       case (gridopt_box)
          call read_kt_grids_box
       end select
    end if

    call broadcast_input

    call allocate_arrays

  end subroutine read_kt_grids_parameters

  subroutine read_grid_option

    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value

    implicit none

    type (text_option), dimension (5), parameter :: gridopts = &
         (/ text_option('default', gridopt_range), &
         text_option('range', gridopt_range), &
         text_option('box', gridopt_box), &
         text_option('annulus', gridopt_box), &
         text_option('nonlinear', gridopt_box) /)
    
    integer :: ierr, in_file
    logical :: nml_exist
    
    grid_option = 'default'
    
    in_file = input_unit_exist ("kt_grids_knobs", nml_exist)
    if (nml_exist) read (unit=in_file, nml=kt_grids_knobs)
    
    ierr = error_unit()
    call get_option_value (grid_option, gridopts, gridopt_switch, &
         ierr, "grid_option in kt_grids_knobs")
    
  end subroutine read_grid_option

  subroutine read_kt_grids_box

    use file_utils, only: input_unit_exist
    use physics_flags, only: full_flux_surface

    implicit none

    integer :: in_file
    logical :: exist

    namelist /kt_grids_box_parameters/ nx, ny, jtwist, y0

    ! note that jtwist and y0 will possibly be modified
    ! later in init_kt_grids_box if they make it out
    ! of this subroutine with negative values
    ! it is necessary to wait until later to do this check
    ! because the values to which they may be set will
    ! depend on information from the geometry module,
    ! which itself may rely on ny from here (number of alphas)
    nx = 1
    ny = 1
    jtwist = -1
    y0 = -1.0
    nalpha = 1

    in_file = input_unit_exist("kt_grids_box_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_box_parameters)

    ! get the number of de-aliased modes in y and x
    naky = (ny-1)/3 + 1
    nakx = 2*((nx-1)/3) +  1

    reality = .true.

    if (full_flux_surface) nalpha = ny

  end subroutine read_kt_grids_box

  subroutine read_kt_grids_range

    use file_utils, only: input_unit, input_unit_exist

    implicit none

    integer :: in_file
    logical :: exist

    namelist /kt_grids_range_parameters/ naky, nakx,  &
         aky_min, aky_max, theta0_min, theta0_max, akx_min, akx_max

    nalpha = 1
    naky = 1
    nakx = 1
    aky_min = 0.0
    aky_max = 0.0
    ! set these to be nonsense values
    ! so can check later if they've been set
    akx_min = 0.0
    akx_max = -1.0
    theta0_min = 0.0
    theta0_max = -1.0    

    in_file = input_unit_exist ("kt_grids_range_parameters", exist)
    if (exist) read (in_file, nml=kt_grids_range_parameters)

  end subroutine read_kt_grids_range

  subroutine init_kt_grids (geo_surf, twist_and_shift_geo_fac)

    use common_types, only: flux_surface_type
    use zgrid, only: init_zgrid
    use zgrid, only: shat_zero
    use physics_flags, only: full_flux_surface

    implicit none

    type (flux_surface_type), intent (in) :: geo_surf
    real, intent (in) :: twist_and_shift_geo_fac

    if (init_kt_grids_initialized) return
    init_kt_grids_initialized = .true.

    call init_zgrid

    select case (gridopt_switch)
    case (gridopt_range)
       call init_kt_grids_range (geo_surf)
    case (gridopt_box)
       call init_kt_grids_box (geo_surf, twist_and_shift_geo_fac)
    end select

    ! determine if iky corresponds to zonal mode
    if (.not.allocated(zonal_mode)) allocate (zonal_mode(naky))
    zonal_mode = .false.
    if (abs(aky(1)) < epsilon(0.)) zonal_mode(1) = .true.

  end subroutine init_kt_grids

  subroutine init_kt_grids_box (geo_surf, twist_and_shift_geo_fac)

    use mp, only: mp_abort
    use common_types, only: flux_surface_type
    use constants, only: pi
    use physics_parameters, only: rhostar
    use physics_flags, only: full_flux_surface
    use zgrid, only: shat_zero

    implicit none
    
    type (flux_surface_type), intent (in) :: geo_surf
    real, intent (in) :: twist_and_shift_geo_fac

    integer :: ikx, iky

    ! set jtwist and y0 for cases where they have not been specified
    ! and for which it makes sense to set them automatically
    if (jtwist < 1) jtwist = max(1,int(abs(twist_and_shift_geo_fac)+0.5))
    ! signed version of jtwist, with sign determined by, e.g., magnetic shear
    ikx_twist_shift = -jtwist*int(sign(1.0,twist_and_shift_geo_fac))

    if (y0 < 0.) then
       if (full_flux_surface) then
          ! if simulating a flux annulus, then
          ! y0 determined by the physical
          ! extent of the device
          if (rhostar > 0.) then
             y0 = 1./(rhostar*geo_surf%rhotor)
          else
             call mp_abort ('must set rhostar if simulating a full flux surface. aborting.')
          end if
       else
          ! if simulating a flux tube
          ! makes no sense to have y0 < 0.0
          ! so abort
          call mp_abort ('y0 negative only makes sense when simulating a flux annulus.  aborting.')
       end if
    end if

    ! get the grid spacing in ky and then in kx using twist-and-shift BC
    dky = 1./y0
    ! non-quantized b/c assumed to be periodic instead 
    ! of linked boundary conditions if zero magnetic shear
    if (abs(geo_surf%shat) <= shat_zero) then
       dkx = dky / real(jtwist)
    else
       dkx = dky * abs(twist_and_shift_geo_fac) / real(jtwist)
    end if

    x0 = 1./dkx


    ! ky goes from zero to ky_max
    do iky = 1, naky
       aky(iky) = real(iky-1)*dky
    end do

    ! get the ikx index corresponding to kx_max
    ikx_max = nakx/2+1

    ! get the total number of ky values, including negative ky
    naky_all = 2*naky-1

    ! kx goes from zero to kx_max down to zero...
    do ikx = 1, ikx_max
       akx(ikx) = real(ikx-1)*dkx
    end do
    ! and then from -kx_max to -|kx_min|
    do ikx = ikx_max+1, nakx
       akx(ikx) = real(ikx-nakx-1)*dkx
    end do

    ! set theta0=0 for ky=0
    theta0(1,:) = 0.0
    if (abs(geo_surf%shat) > shat_zero) then
       do ikx = 1, nakx
          ! theta0 = kx/ky/shat
          theta0(2:,ikx) = akx(ikx)/(aky(2:)*geo_surf%shat)
       end do
    else
       do ikx = 1, nakx
          ! if shat=0, theta0 is meaningless, so be careful
          theta0(2:,ikx) = - akx(ikx)/aky(2:)
       end do
    end if

    ! for radial variation
    if(.not.allocated(x)) allocate (x(nx))

    dx = (2*pi*x0)/nx
    dy = (2*pi*y0)/ny
    do ikx = 1, nx
      x(ikx) = (ikx-0.5)*dx - pi*x0
    enddo
    
  end subroutine init_kt_grids_box

  subroutine init_kt_grids_range (geo_surf)

    use mp, only: mp_abort
    use common_types, only: flux_surface_type
    use zgrid, only: shat_zero

    implicit none

    type (flux_surface_type), intent (in) :: geo_surf

    integer :: i, j
    real :: dkx, dky, dtheta0
    real :: zero

    ! NB: we are assuming here that all ky are positive
    ! when running in range mode
    dky = 0.0
    if (naky > 1) dky = (aky_max - aky_min)/real(naky - 1)
    aky = (/ (aky_min + dky*real(i), i = 0,naky-1) /)

    ! set default akx and theta0 to 0
    akx = 0.0 ; theta0 = 0.0

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
       if (theta0_min > theta0_max+zero .and. abs(aky(1)) > zero) then
          theta0_min = akx_min/(geo_surf%shat*aky(1))
          theta0_max = akx_max/(geo_surf%shat*aky(1))
          dtheta0 = 0.0
          if (nakx > 1) dtheta0 = (theta0_max - theta0_min)/real(nakx - 1)
          
          do j = 1, naky
             theta0(j,:) &
                  = (/ (theta0_min + dtheta0*real(i), i=0,nakx-1) /)
          end do
          akx = theta0(1,:) * geo_surf%shat * aky(1)
       else if (akx_max > akx_min-zero) then
          dkx = 0.0
          if (nakx > 1) dkx = (akx_max - akx_min)/real(nakx - 1)
          akx = (/ (akx_min + dkx*real(i), i = 0,nakx-1) /)

          dtheta0 = 0.0
          if (nakx > 1) dtheta0 = (theta0_max - theta0_min)/real(nakx - 1)

          if (geo_surf%shat > epsilon(0.)) then
             do j = 1, naky
                theta0(j,:) &
                     = (/ (theta0_min + dtheta0*real(i), i=0,nakx-1) /)
             end do
          else
             do j = 1, naky
                theta0(j,:) &
                     = (/ (theta0_min + dtheta0*real(i), i=nakx-1,0,-1) /)
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

    ikx_max = nakx
    naky_all = naky

  end subroutine init_kt_grids_range

  subroutine broadcast_input

    use mp, only: broadcast

    implicit none

    call broadcast (gridopt_switch)
    call broadcast (naky)
    call broadcast (nakx)
    call broadcast (ny)
    call broadcast (nx)
    call broadcast (nalpha)
    call broadcast (reality)
    call broadcast (jtwist)
    call broadcast (y0)
    call broadcast (aky_min)
    call broadcast (aky_max)
    call broadcast (akx_min)
    call broadcast (akx_max)
    call broadcast (theta0_min)
    call broadcast (theta0_max)

  end subroutine broadcast_input

  subroutine allocate_arrays

    implicit none

    allocate (akx(nakx))
    allocate (aky(naky))
    allocate (theta0(naky,nakx))

  end subroutine allocate_arrays

  ! take an array with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
  ! and uses reality condition to return array
  ! with kx >= 0 and all ky (ordered like 0, ..., kymax, -kymax, ..., -dky)
  subroutine swap_kxky_complex (gin, gout)

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
    do iky = naky+1, naky_all
       ! this is the ky index corresponding to +ky in original array
       ikyneg = naky_all-iky+2
       gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
    end do
    do ikx = 2, ikx_max
       ikxneg = nakx-ikx+2
       do iky = naky+1, naky_all
          ! this is the ky index corresponding to +ky in original array
          ikyneg = naky_all-iky+2
          gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
       end do
    end do

  end subroutine swap_kxky_complex

  ! take an array with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
  ! and uses reality condition to return array
  ! with kx >= 0 and all ky (ordered like 0, ..., kymax, -kymax, ..., -dky)
  subroutine swap_kxky_real (gin, gout)

    implicit none

    real, dimension (:,:), intent (in) :: gin
    real, dimension (:,:), intent (out) :: gout

    integer :: ikx, ikxneg
    integer :: iky, ikyneg

    ! first set arrays equal for ky >= 0 and kx >= 0
    gout(:naky,:) = gin(:,:ikx_max)
    ! next fill in ky < 0, kx >= 0 elements of array using reality
    ikx = 1
    ikxneg = ikx
    do iky = naky+1, naky_all
       ! this is the ky index corresponding to +ky in original array
       ikyneg = naky_all-iky+2
       gout(iky,ikx) = gin(ikyneg,ikxneg)
    end do
    do ikx = 2, ikx_max
       ikxneg = nakx-ikx+2
       do iky = naky+1, naky_all
          ! this is the ky index corresponding to +ky in original array
          ikyneg = naky_all-iky+2
          gout(iky,ikx) = gin(ikyneg,ikxneg)
       end do
    end do

  end subroutine swap_kxky_real

  ! take an array with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
  ! and uses reality condition to return array
  ! with kx >= 0 and all ky (ordered like -kymax, ..., 0, ..., kymax)
  subroutine swap_kxky_ordered (gin, gout)

    implicit none

    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    integer :: ikx, ikxneg
    integer :: iky, ikyneg

    ! first set arrays equal for ky >= 0 and kx >= 0
    gout(naky:,:) = gin(:,:ikx_max)
    ! next fill in ky < 0, kx >= 0 elements of array using reality
    ikx = 1
    ikxneg = ikx
    do iky = 1, naky-1
       ! this is the ky index corresponding to +ky in original array
       ikyneg = naky-iky+1
       gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
    end do
    do ikx = 2, ikx_max
       ikxneg = nakx-ikx+2
       do iky = 1, naky-1
          ! this is the ky index corresponding to +ky in original array
          ikyneg = naky-iky+1
          gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
       end do
    end do

  end subroutine swap_kxky_ordered

  ! take an array with kx >= 0 and all ky (ordered like 0, ..., kymax, -kymax, ..., -dky)
  ! and uses reality condition to return array
  ! with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
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
          ikyneg = naky_all-iky+2
          gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
       end do
    end do

  end subroutine swap_kxky_back

  ! take an array with kx >= 0 and all ky (ordered like -kymax, ..., 0, ..., kymax)
  ! and uses reality condition to return array
  ! with ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
  subroutine swap_kxky_back_ordered (gin, gout)

    implicit none

    complex, dimension (:,:), intent (in) :: gin
    complex, dimension (:,:), intent (out) :: gout

    integer :: ikx, ikxneg
    integer :: iky, ikyneg

    ! first set arrays equal for ky >= 0 and kx >= 0
    gout(:,:ikx_max) = gin(naky:,:)
    ! next fill in kx < 0, ky >= 0 elements of array using reality
    do ikx = ikx_max+1, nakx
       ikxneg = nakx-ikx+2
       do iky = 1, naky
          ikyneg = naky-iky+1
          gout(iky,ikx) = conjg(gin(ikyneg,ikxneg))
       end do
    end do

  end subroutine swap_kxky_back_ordered

  subroutine finish_kt_grids

    implicit none

    if (allocated(aky)) deallocate (aky)
    if (allocated(akx)) deallocate (akx)
    if (allocated(theta0)) deallocate (theta0)

    if (allocated(x)) deallocate (x)

    reality = .false.
    read_kt_grids_initialized = .false.
    init_kt_grids_initialized = .false.

  end subroutine finish_kt_grids

end module kt_grids

