module vmec_geo

  implicit none

  public :: read_vmec_parameters
  public :: get_vmec_geo

  integer :: nalpha
  integer :: surface_option
  real :: nfield_periods
  real :: zeta_center, torflux
  logical :: verbose
  character (2000) :: vmec_filename
  
contains

  subroutine read_vmec_parameters (nalpha_out)

    use file_utils, only: input_unit_exist

    implicit none

    integer, intent (out) :: nalpha_out

    integer :: in_file
    logical :: exist

    namelist /vmec_parameters/ nalpha, zeta_center, nfield_periods, &
         torflux, surface_option, verbose, vmec_filename

    call init_vmec_defaults

    in_file = input_unit_exist("vmec_parameters", exist)
    if (exist) read (unit=in_file, nml=vmec_parameters)

    nalpha_out = nalpha

  end subroutine read_vmec_parameters

  subroutine init_vmec_defaults

    implicit none

    vmec_filename = 'equilibria/wout_w7x_standardConfig.nc'
    nalpha = 5
    zeta_center = 0.0
    nfield_periods = -1.0
    torflux = 0.6354167d+0
    surface_option = 0
    verbose = .true.

  end subroutine init_vmec_defaults

  subroutine get_vmec_geo (nzgrid, surf, grho, bmag, gradpar, gds2, gds21, gds22, &
       gds23, gds24, gds25, gds26, gbdrift, gbdrift0, cvdrift, cvdrift0, theta_vmec, &
       zed_scalefac)

    use common_types, only: flux_surface_type
    use vmec_to_gs2_geometry_interface_mod, only: vmec_to_gs2_geometry_interface

    implicit none

    integer, intent (in) :: nzgrid
    type (flux_surface_type), intent (out) :: surf
    real, dimension (:,-nzgrid:), intent (out) :: grho, bmag, gradpar, gds2, gds21, gds22, &
         gds23, gds24, gds25, gds26, gbdrift, gbdrift0, cvdrift, cvdrift0, theta_vmec
    real, intent (out) :: zed_scalefac

    integer :: i, j
    real :: L_reference, B_reference, nfp

    real, dimension (nalpha) :: alpha
    real, dimension (-nzgrid:nzgrid) :: zeta
    real, dimension (-nzgrid:nzgrid) :: theta

    call vmec_to_gs2_geometry_interface (vmec_filename, nalpha, nzgrid, &
         zeta_center, nfield_periods, torflux, surface_option, verbose, &
         surf%rhoc, surf%qinp, surf%shat, L_reference, B_reference, nfp, &
         alpha, zeta, bmag, gradpar, gds2, gds21, gds22, gds23, gds24, &
         gds25, gds26, gbdrift, gbdrift0, cvdrift, cvdrift0, theta_vmec)

    ! vmec_to_gs2_geometry_interface returns psitor/psitor_lcfs as rhoc
    ! stella uses rhoc = sqrt(psitor/psitor_lcfs) = rhotor
    surf%rhoc = sqrt(surf%rhoc)
    surf%rhotor = surf%rhoc
    ! grho = |grad rho| = |drho/dx| * |grad x|
    ! |drho/dx| = L_reference
    ! gds22 = shat^2 * |grad x|^2
    grho = sqrt(gds22/surf%shat**2)/L_reference
    surf%drhotordrho = 1.0
    surf%psitor_lcfs = 0.5

    ! scale the vmec output
    ! alpha = theta_pest - iota*zeta
    ! theta_pest = theta_vmec + Lambda(psi,alpha,theta_vmec)
    ! with theta_pest a straight-field-line angle
    ! but not theta_vmec
    ! so theta is theta_pest up to constant (alpha)

    ! scale zed so that it is zeta compressed (or expanded)
    ! to the range [-pi,pi]
    zed_scalefac = real(nfp)/nfield_periods

!    theta = zeta/nfp/surf%qinp
    theta = alpha+zeta/surf%qinp
    ! this is b . grad zed
    ! with zed = zeta scaled to run from -pi to pi
!    gradpar = gradpar/nfp/surf%qinp
    gradpar = gradpar*zed_scalefac
    gds23 = gds23*zed_scalefac
    gds24 = gds24*zed_scalefac
    ! this is the vmec theta (not straight-field-line coordinate)
    ! scaled to run between -pi and pi
    theta_vmec = theta_vmec/nfp

    open (2001,file='vmec.geo',status='unknown')
    write (2001,'(5a12)') 'rhotor', 'qinp', 'shat', 'aref', 'Bref'
    write (2001,'(5e12.4)') surf%rhoc, surf%qinp, surf%shat, L_reference, B_reference
    write (2001,*)
    write (2001,'(9a12)') 'alpha', 'zeta', 'bmag', 'gradpar', 'gds2', 'gds21', 'gds22', 'gds23', 'gds24'
    do j = -nzgrid, nzgrid
       do i = 1, nalpha
          write (2001,'(9e12.4)') alpha(i), zeta(j), bmag(i,j), gradpar(i,j), &
               gds2(i,j), gds21(i,j), gds22(i,j), gds23(i,j), gds24(i,j)
       end do
    end do
    write (2001,*)
    write (2001,'(7a12)') 'alpha', 'zeta', 'gbdrift', 'gbdrift0', 'cvdrift', 'cvdrift0', 'theta_vmec'
    do j = -nzgrid, nzgrid
       do i = 1, nalpha
          write (2001,'(7e12.4)') alpha(i), zeta(j), gbdrift(i,j), gbdrift0(i,j), cvdrift(i,j), cvdrift0(i,j), theta_vmec(i,j)
       end do
    end do
    close (2001)

  end subroutine get_vmec_geo

end module vmec_geo
