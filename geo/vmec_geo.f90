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
    nfield_periods = 1.0
    torflux = 0.6354167d+0
    surface_option = 0
    verbose = .true.

  end subroutine init_vmec_defaults

  subroutine get_vmec_geo (nzgrid, surf, bmag, gradpar, gds2, gds21, gds22, &
       gbdrift, gbdrift0, cvdrift, cvdrift0)

    use common_types, only: flux_surface_type
    use vmec_to_gs2_geometry_interface_mod, only: vmec_to_gs2_geometry_interface

    implicit none

    integer, intent (in) :: nzgrid
    type (flux_surface_type), intent (out) :: surf
    real, dimension (:,-nzgrid:), intent (out) :: bmag, gradpar, gds2, gds21, gds22, &
         gbdrift, gbdrift0, cvdrift, cvdrift0

    integer :: i, j
    real :: L_reference, B_reference

    real, dimension (nalpha) :: alpha
    real, dimension (-nzgrid:nzgrid) :: zeta

    call vmec_to_gs2_geometry_interface (vmec_filename, nalpha, nzgrid, &
         zeta_center, nfield_periods, torflux, surface_option, verbose, &
         surf%rhoc, surf%qinp, surf%shat, L_reference, B_reference, &
         alpha, zeta, bmag, gradpar, gds2, gds21, gds22, &
         gbdrift, gbdrift0, cvdrift, cvdrift0)

    open (2001,file='vmec.geo',status='unknown')
    write (2001,'(5a12)') 'torflux', 'qinp', 'shat', 'aref', 'Bref'
    write (2001,'(5e12.4)') surf%rhoc, surf%qinp, surf%shat, L_reference, B_reference
    write (2001,*)
    write (2001,'(7a12)') 'alpha', 'zeta', 'bmag', 'gradpar', 'gds2', 'gds21', 'gds22'
    do j = -nzgrid, nzgrid
       do i = 1, nalpha
          write (2001,'(7e12.4)') alpha(i), zeta(j), bmag(i,j), gradpar(i,j), gds2(i,j), gds21(i,j), gds22(i,j)
       end do
    end do
    write (2001,*)
    write (2001,'(6a12)') 'alpha', 'zeta', 'gbdrift', 'gbdrift0', 'cvdrift', 'cvdrift0'
    do j = -nzgrid, nzgrid
       do i = 1, nalpha
          write (2001,'(6e12.4)') alpha(i), zeta(j), gbdrift(i,j), gbdrift0(i,j), cvdrift(i,j), cvdrift0(i,j)
       end do
    end do
    close (2001)

  end subroutine get_vmec_geo

end module vmec_geo
