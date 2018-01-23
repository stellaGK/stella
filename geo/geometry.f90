module geometry

  use common_types, only: flux_surface_type

  implicit none

  public :: init_geometry, finish_geometry
  public :: grho
  public :: bmag, dbdzed, btor
  public :: gradpar
  public :: cvdrift, cvdrift0
  public :: gbdrift, gbdrift0
  public :: gds2, gds21, gds22, gds23, gds24
  public :: jacob
  public :: drhodpsi
  public :: dl_over_b
  public :: dBdrho, d2Bdrdth, dgradpardrho, dIdrho
  public :: geo_surf
  public :: Rmajor
  public :: nalpha

  private

  type (flux_surface_type) :: geo_surf

  real :: dIdrho
  real :: drhodpsi, rhotor, drhotordrho, shat, qinp, rgeo
  real, dimension (:), allocatable :: grho
  real, dimension (:,:), allocatable :: gradpar
  real, dimension (:,:), allocatable :: bmag, dbdzed
  real, dimension (:,:), allocatable :: cvdrift, cvdrift0
  real, dimension (:,:), allocatable :: gbdrift, gbdrift0
  real, dimension (:,:), allocatable :: gds2, gds21, gds22, gds23, gds24
  real, dimension (:), allocatable :: jacob
  real, dimension (:), allocatable :: dl_over_b
  real, dimension (:), allocatable :: dBdrho, d2Bdrdth, dgradpardrho
  real, dimension (:), allocatable :: btor, Rmajor

  integer :: geo_option_switch
  integer, parameter :: geo_option_local = 1
  integer, parameter :: geo_option_inputprof = 2
  integer, parameter :: geo_option_vmec = 3

  ! number of field line labels to include
  ! default is one (only > 1 for alpha_global = .true.)
  integer :: nalpha = 1

  logical :: geoinit = .false.

contains

  subroutine init_geometry

    use mp, only: proc0
    use millerlocal, only: read_local_parameters, get_local_geo
    use vmec_geo, only: read_vmec_parameters, get_vmec_geo
    use inputprofiles_interface, only: read_inputprof_geo
    use zgrid, only: nzed, nzgrid
    use zgrid, only: zed, delzed
    use zgrid, only: shat_zero
    use zgrid, only: boundary_option_switch, boundary_option_self_periodic

    implicit none

    real :: dpsidrho
    integer :: iy

    if (geoinit) return
    geoinit = .true.

    if (proc0) then
       call read_parameters
       select case (geo_option_switch)
       case (geo_option_local)
          ! read in Miller local parameters
          call read_local_parameters (geo_surf)
          ! allocate geometry arrays
          call allocate_arrays (nalpha, nzgrid)
          ! use Miller local parameters to get 
          ! geometric coefficients needed by stella
          call get_local_geo (nzed, nzgrid, zed, &
               dpsidrho, dIdrho, grho, bmag(1,:), &
               gds2(1,:), gds21(1,:), gds22(1,:), &
               gds23(1,:), gds24(1,:), gradpar(1,:), &
               gbdrift0(1,:), gbdrift(1,:), cvdrift0(1,:), cvdrift(1,:), &
               dBdrho, d2Bdrdth, dgradpardrho, btor, &
               rmajor)
          drhodpsi = 1./dpsidrho
       case (geo_option_inputprof)
          ! first read in some local parameters
          ! only thing needed really is rhoc
          call read_local_parameters (geo_surf)
          ! allocate geometry arrays
          call allocate_arrays (nalpha, nzgrid)
          ! now overwrite local parameters
          ! with those from input.profiles file
          ! use rhoc from input as surface
          call read_inputprof_geo (geo_surf)
          call get_local_geo (nzed, nzgrid, zed, &
               dpsidrho, dIdrho, grho, bmag(1,:), &
               gds2(1,:), gds21(1,:), gds22(1,:), &
               gds23(1,:), gds24(1,:), gradpar(1,:), &
               gbdrift0(1,:), gbdrift(1,:), cvdrift0(1,:), cvdrift(1,:), &
               dBdrho, d2Bdrdth, dgradpardrho, btor, &
               rmajor)
          drhodpsi = 1./dpsidrho
       case (geo_option_vmec)
          ! read in input parameters for vmec
          ! nalpha may be specified via input file
          call read_vmec_parameters (nalpha)
          ! allocate geometry arrays
          call allocate_arrays (nalpha, nzgrid)
          ! get geometry coefficients from vmec
          call get_vmec_geo (nzgrid, geo_surf, bmag, gradpar, gds2, gds21, gds22, &
               gbdrift, gbdrift0, cvdrift, cvdrift0)
          ! FLAG -- NOT SURE IF THIS IS CORRECT
          drhodpsi = 1.0
          ! FLAG -- NEED TO SEE IF MATT CAN PROVIDE THESE
          gds23 = 0. ; gds24 = 0.
       end select
    end if

    call broadcast_arrays

    ! FLAG -- THIS SHOULD BE GENERALIZED TO ACCOUNT FOR ALPHA VARIATION
    jacob = 1.0/(drhodpsi*gradpar(1,:)*bmag(1,:))
    
    dl_over_b = delzed*jacob
    dl_over_b = dl_over_b / sum(dl_over_b)

    do iy = 1, nalpha
       call get_dzed (nzgrid, delzed, bmag(iy,:), dbdzed(iy,:))
    end do

    ! if magnetic shear almost zero, override parallel
    ! boundary condition so that it is periodic
    if(abs(geo_surf%shat) <=  shat_zero) &
         boundary_option_switch = boundary_option_self_periodic

  end subroutine init_geometry

  subroutine allocate_arrays (nalpha, nzgrid)

    implicit none

    integer, intent (in) :: nalpha, nzgrid

    if (.not.allocated(bmag)) allocate (bmag(nalpha,-nzgrid:nzgrid))
    if (.not.allocated(gds2)) allocate (gds2(nalpha,-nzgrid:nzgrid))
    if (.not.allocated(gds21)) allocate (gds21(nalpha,-nzgrid:nzgrid))
    if (.not.allocated(gds22)) allocate (gds22(nalpha,-nzgrid:nzgrid))
    if (.not.allocated(gds23)) allocate (gds23(nalpha,-nzgrid:nzgrid))
    if (.not.allocated(gds24)) allocate (gds24(nalpha,-nzgrid:nzgrid))
    if (.not.allocated(gbdrift)) allocate (gbdrift(nalpha,-nzgrid:nzgrid))
    if (.not.allocated(gbdrift0)) allocate (gbdrift0(nalpha,-nzgrid:nzgrid))
    if (.not.allocated(cvdrift)) allocate (cvdrift(nalpha,-nzgrid:nzgrid))
    if (.not.allocated(cvdrift0)) allocate (cvdrift0(nalpha,-nzgrid:nzgrid))
    if (.not.allocated(dbdzed)) allocate (dbdzed(nalpha,-nzgrid:nzgrid))

    ! FLAG - NEED TO SORT OUT 1D VS 2D FOR GRADPAR
    if (.not.allocated(gradpar)) allocate (gradpar(nalpha,-nzgrid:nzgrid))

    if (.not.allocated(grho)) allocate (grho(-nzgrid:nzgrid))
    if (.not.allocated(btor)) allocate (btor(-nzgrid:nzgrid))
    if (.not.allocated(rmajor)) allocate (rmajor(-nzgrid:nzgrid))
    if (.not.allocated(jacob)) allocate (jacob(-nzgrid:nzgrid))
    if (.not.allocated(dl_over_b)) allocate (dl_over_b(-nzgrid:nzgrid))
    if (.not.allocated(dBdrho)) allocate (dBdrho(-nzgrid:nzgrid))
    if (.not.allocated(d2Bdrdth)) allocate (d2Bdrdth(-nzgrid:nzgrid))
    if (.not.allocated(dgradpardrho)) allocate (dgradpardrho(-nzgrid:nzgrid))

  end subroutine allocate_arrays

  subroutine read_parameters

    use text_options
    use file_utils, only: error_unit, input_unit_exist

    implicit none

    integer :: in_file, ierr
    logical :: exist

    character (20) :: geo_option
    type (text_option), dimension (5), parameter :: geoopts = (/ &
         text_option('default', geo_option_local), &
         text_option('miller', geo_option_local), &
         text_option('local', geo_option_local), &
         text_option('input.profiles', geo_option_inputprof), &
         text_option('vmec', geo_option_vmec) /)

    namelist /geo_knobs/ geo_option

    geo_option = 'local'

    in_file = input_unit_exist("geo_knobs", exist)
    if (exist) read (unit=in_file, nml=geo_knobs)

    ierr = error_unit()
    call get_option_value &
         (geo_option, geoopts, geo_option_switch, &
         ierr, "geo_option in geo_knobs")
    
  end subroutine read_parameters

  subroutine broadcast_arrays

    use mp, only: broadcast

    implicit none

    call broadcast (qinp)
    call broadcast (shat)
    call broadcast (drhodpsi)
    call broadcast (dIdrho)
    call broadcast (grho)
    call broadcast (bmag)
    call broadcast (btor)
    call broadcast (rmajor)
    call broadcast (gradpar)
    call broadcast (gds2)
    call broadcast (gds21)
    call broadcast (gds22)
    call broadcast (gds23)
    call broadcast (gds24)
    call broadcast (gbdrift0)
    call broadcast (gbdrift)
    call broadcast (cvdrift0)
    call broadcast (cvdrift)
    call broadcast (dBdrho)
    call broadcast (d2Bdrdth)
    call broadcast (dgradpardrho)

    call broadcast (geo_surf%rmaj)
    call broadcast (geo_surf%rgeo)
    call broadcast (geo_surf%kappa)
    call broadcast (geo_surf%kapprim)
    call broadcast (geo_surf%tri)
    call broadcast (geo_surf%triprim)
    call broadcast (geo_surf%rhoc)
    call broadcast (geo_surf%dr)
    call broadcast (geo_surf%shift)
    call broadcast (geo_surf%qinp)
    call broadcast (geo_surf%shat)
    call broadcast (geo_surf%betaprim)
    call broadcast (geo_surf%betadbprim)
    call broadcast (geo_surf%d2qdr2)
    call broadcast (geo_surf%d2psidr2)
    call broadcast (geo_surf%dpsitordrho)
    call broadcast (geo_surf%rhotor)
    call broadcast (geo_surf%drhotordrho)

  end subroutine broadcast_arrays

  ! given function f(z:-pi->pi), calculate z derivative
  ! second order accurate, with equal grid spacing assumed
  ! assumes periodic in z -- may need to change this in future
  subroutine get_dzed (nz, dz, f, df)

    implicit none

    integer, intent (in) :: nz
    real, dimension (-nz:), intent (in) :: dz, f
    real, dimension (-nz:), intent (out) :: df

    df(-nz+1:nz-1) = (f(-nz+2:)-f(:nz-2))/(dz(:nz-2)+dz(-nz+1:nz-1))

    ! FLAG -- THIS MAY NEED TO BE CHANGED
    ! use periodicity at boundary
    df(-nz) = (f(-nz+1)-f(nz-1))/(dz(-nz)+dz(nz-1))
    df(nz) = df(-nz)

  end subroutine get_dzed

  subroutine finish_geometry

    implicit none

    if (allocated(grho)) deallocate (grho)
    if (allocated(bmag)) deallocate (bmag)
    if (allocated(btor)) deallocate (btor)
    if (allocated(rmajor)) deallocate (rmajor)
    if (allocated(dbdzed)) deallocate (dbdzed)
    if (allocated(jacob)) deallocate (jacob)
    if (allocated(gradpar)) deallocate (gradpar)
    if (allocated(dl_over_b)) deallocate (dl_over_b)
    if (allocated(gds2)) deallocate (gds2)
    if (allocated(gds21)) deallocate (gds21)
    if (allocated(gds22)) deallocate (gds22)
    if (allocated(gds23)) deallocate (gds23)
    if (allocated(gds24)) deallocate (gds24)
    if (allocated(gbdrift)) deallocate (gbdrift)
    if (allocated(gbdrift0)) deallocate (gbdrift0)
    if (allocated(cvdrift)) deallocate (cvdrift)
    if (allocated(cvdrift0)) deallocate (cvdrift0)
    if (allocated(dBdrho)) deallocate (dBdrho)
    if (allocated(d2Bdrdth)) deallocate (d2Bdrdth)
    if (allocated(dgradpardrho)) deallocate (dgradpardrho)

    geoinit = .false.

  end subroutine finish_geometry

end module geometry
