module vmec_geo

  implicit none

  public :: read_vmec_parameters
  public :: get_vmec_geo

  integer :: nalpha
  real :: alpha0
  real :: nzgrid_scalefac
  integer :: surface_option
  real :: nfield_periods
  real :: zeta_center, torflux
  logical :: verbose
  character (2000) :: vmec_filename
  
contains

  subroutine read_vmec_parameters (nalpha_out)

    use file_utils, only: input_unit_exist
    use mp, only: mp_abort
    use zgrid, only: zed_equal_arc

    implicit none

    integer, intent (out) :: nalpha_out

    integer :: in_file
    logical :: exist

    namelist /vmec_parameters/ nalpha, alpha0, zeta_center, nfield_periods, &
         torflux, nzgrid_scalefac, surface_option, verbose, vmec_filename

    call init_vmec_defaults

    in_file = input_unit_exist("vmec_parameters", exist)
    if (exist) read (unit=in_file, nml=vmec_parameters)

    if (nzgrid_scalefac < 1.0-epsilon(0.)) then
       write (*,*) 'nzgrid_scalefac = ', nzgrid_scalefac
       call mp_abort ('nzgrid_scalefac should always be >= 1.0.  aborting')
    else if (nzgrid_scalefac > 1.0+epsilon(0.) .and. .not.zed_equal_arc) then
       write (*,*) 'There is no reason to use nzgrid_scalefac different from 1 unless zed_equal_arc=T'
       write (*,*) 'Setting nzgrid_scalefac = 1.0'
       nzgrid_scalefac = 1.0
    end if

    nalpha_out = nalpha

  end subroutine read_vmec_parameters

  subroutine init_vmec_defaults

    use zgrid, only: zed_equal_arc

    implicit none

    vmec_filename = 'equilibria/wout_w7x_standardConfig.nc'
    nalpha = 1
    alpha0 = 0.0
    zeta_center = 0.0
    nfield_periods = -1.0
    torflux = 0.6354167d+0
    surface_option = 0
    verbose = .true.
    ! if simulating entire flux surface,
    ! must obtain vmec geo quantities
    ! on zeta grid that is longer than
    ! will ultimately be used in simulation
    ! this is related to need for gradpar to
    ! be independent of alpha
    if (zed_equal_arc) then
       nzgrid_scalefac = 2.0
    else
       nzgrid_scalefac = 1.0
    end if

  end subroutine init_vmec_defaults

  subroutine get_vmec_geo (nzgrid, surf, grho, bmag, gradpar, gds2, gds21, gds22, &
       gds23, gds24, gds25, gds26, gbdrift, gbdrift0, cvdrift, cvdrift0, sign_torflux, &
       theta_vmec, zed_scalefac, L_reference, B_reference, alpha)

    use common_types, only: flux_surface_type
    use vmec_to_gs2_geometry_interface_mod, only: vmec_to_gs2_geometry_interface

    implicit none

    integer, intent (in) :: nzgrid
    type (flux_surface_type), intent (out) :: surf
    real, dimension (:,-nzgrid:), intent (out) :: grho, bmag, gradpar, gds2, gds21, gds22, &
         gds23, gds24, gds25, gds26, gbdrift, gbdrift0, cvdrift, cvdrift0, theta_vmec
    real, dimension (:), intent (out) :: alpha
    real, intent (out) :: zed_scalefac, L_reference, B_reference
    integer, intent (out) :: sign_torflux

    integer :: i, j
    integer :: nzgrid_vmec
    real :: nfp

    real, dimension (:), allocatable :: zeta_vmec
    real, dimension (:,:), allocatable :: bmag_vmec, gradpar_vmec
    real, dimension (:,:), allocatable :: gds2_vmec, gds21_vmec, gds22_vmec
    real, dimension (:,:), allocatable :: gds23_vmec, gds24_vmec, gds25_vmec, gds26_vmec
    real, dimension (:,:), allocatable :: gbdrift_vmec, gbdrift0_vmec
    real, dimension (:,:), allocatable :: cvdrift_vmec, cvdrift0_vmec

!    real, dimension (nalpha) :: alpha
    real, dimension (-nzgrid:nzgrid) :: zeta
    real, dimension (nalpha,-nzgrid:nzgrid) :: theta

    ! nzgrid_vmec is the number of positive/negative zeta locations
    ! at which to get geometry data from vmec
    ! can be > than nzgrid for full_flux_surface case
    ! where z(zeta_max)-z(zeta_min) varies with alpha
    ! and thus a larger than usual range of zeta_max/min
    ! values are needed to avoid extrapolation
    nzgrid_vmec = nint(nzgrid*nzgrid_scalefac)

    ! allocate arrays of size 2*nzgrid_vmec+1
    allocate (zeta_vmec(-nzgrid_vmec:nzgrid_vmec))
    allocate (bmag_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))
    allocate (gradpar_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))
    allocate (gds2_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))
    allocate (gds21_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))
    allocate (gds22_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))
    allocate (gds23_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))
    allocate (gds24_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))
    allocate (gds25_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))
    allocate (gds26_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))
    allocate (gbdrift_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))
    allocate (gbdrift0_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))
    allocate (cvdrift_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))
    allocate (cvdrift0_vmec(nalpha,-nzgrid_vmec:nzgrid_vmec))

    call vmec_to_gs2_geometry_interface (vmec_filename, nalpha, alpha0, &
         nzgrid_vmec, zeta_center, nfield_periods*nzgrid_scalefac, torflux, &
         surface_option, verbose, &
         surf%rhoc, surf%qinp, surf%shat, L_reference, B_reference, nfp, &
         sign_torflux, alpha, zeta_vmec, &
         bmag_vmec, gradpar_vmec, gds2_vmec, gds21_vmec, &
         gds22_vmec, gds23_vmec, gds24_vmec, &
         gds25_vmec, gds26_vmec, gbdrift_vmec, gbdrift0_vmec, cvdrift_vmec, &
         cvdrift0_vmec, theta_vmec)

    if (nzgrid_vmec /= nzgrid) then
       ! do some stuff
    else
       zeta = zeta_vmec
       bmag = bmag_vmec
       gradpar = gradpar_vmec
       gds2 = gds2_vmec
       gds21 = gds21_vmec
       gds22 = gds22_vmec
       gds23 = gds23_vmec
       gds24 = gds24_vmec
       gds25 = gds25_vmec
       gds26 = gds26_vmec
       gbdrift = gbdrift_vmec
       gbdrift0 = gbdrift0_vmec
       cvdrift = cvdrift_vmec
       cvdrift0 = cvdrift0_vmec
    end if
    
    ! arrays over extended zeta-grid no longer needed, so deallocate
    deallocate (zeta_vmec)
    deallocate (bmag_vmec, gradpar_vmec)
    deallocate (gds2_vmec, gds21_vmec, gds22_vmec)
    deallocate (gds23_vmec, gds24_vmec, gds25_vmec, gds26_vmec)
    deallocate (gbdrift_vmec, gbdrift0_vmec)
    deallocate (cvdrift_vmec, cvdrift0_vmec)
    
    ! vmec_to_gs2_geometry_interface returns psitor/psitor_lcfs as rhoc
    ! stella uses rhoc = sqrt(psitor/psitor_lcfs) = rhotor
    surf%rhoc = sqrt(surf%rhoc)
    surf%rhotor = surf%rhoc
    ! grho = |grad rho| = |drho/dx| * |grad x|
    ! |drho/dx| = L_reference
    ! gds22 = shat^2 * |grad x|^2
    grho = sqrt(gds22/surf%shat**2)/L_reference
    surf%drhotordrho = 1.0
    surf%psitor_lcfs = 0.5*sign_torflux

    ! scale the vmec output
    ! alpha = theta_pest - iota*zeta
    ! theta_pest = theta_vmec + Lambda(psi,alpha,theta_vmec)
    ! with theta_pest a straight-field-line angle
    ! but not theta_vmec
    ! so theta is theta_pest up to constant (alpha)

    ! scale zed so that it is zeta compressed (or expanded)
    ! to the range [-pi,pi]
    ! this is 1/p from stella JCP paper
    zed_scalefac = real(nfp)/nfield_periods

!    theta = zeta/nfp/surf%qinp
    theta = spread(alpha,2,2*nzgrid+1)+spread(zeta,1,nalpha)/surf%qinp
    ! this is b . grad zed
    ! with zed = zeta scaled to run from -pi to pi
!    gradpar = gradpar/nfp/surf%qinp
    gradpar = gradpar*zed_scalefac
    gds23 = gds23*zed_scalefac
    gds24 = gds24*zed_scalefac
    ! this is the vmec theta (not straight-field-line coordinate)
    ! scaled to run between -pi and pi
    theta_vmec = theta_vmec/nfp

    open (2001,file='vmec_geo',status='unknown')
    write (2001,'(5a12)') 'rhotor', 'qinp', 'shat', 'aref', 'Bref'
    write (2001,'(5e12.4)') surf%rhoc, surf%qinp, surf%shat, L_reference, B_reference
    write (2001,*)
    write (2001,'(14a12)') '#    alpha', 'zeta', 'bmag', 'gradpar', 'gds2',&
         'gds21', 'gds22', 'gds23', 'gds24','gbdrift', 'gbdrift0', 'cvdrift',&
         'cvdrift0', 'theta_vmec'
    do j = -nzgrid, nzgrid
       do i = 1, nalpha
          write (2001,'(14e12.4)') alpha(i), zeta(j), bmag(i,j), gradpar(i,j), &
               gds2(i,j), gds21(i,j), gds22(i,j), gds23(i,j), gds24(i,j), &
               gbdrift(i,j), gbdrift0(i,j), cvdrift(i,j), cvdrift0(i,j), theta_vmec(i,j)
       end do
    end do
    close (2001)

  end subroutine get_vmec_geo

end module vmec_geo
