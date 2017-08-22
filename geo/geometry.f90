module geometry

  implicit none

  public :: init_geometry, finish_geometry
  public :: grho
  public :: bmag, dbdthet
  public :: gradpar
  public :: cvdrift, cvdrift0
  public :: gbdrift, gbdrift0
  public :: gds2, gds21, gds22
  public :: jacob
  public :: drhodpsi
  public :: dl_over_b
  public :: shat, qinp

  private

  real :: drhodpsi, shat, qinp
  real, dimension (:), allocatable :: grho, bmag, dbdthet
  real, dimension (:), allocatable :: cvdrift, cvdrift0
  real, dimension (:), allocatable :: gbdrift, gbdrift0
  real, dimension (:), allocatable :: gds2, gds21, gds22
  real, dimension (:), allocatable :: jacob, gradpar
  real, dimension (:), allocatable :: dl_over_b

  integer :: geo_option_switch
  integer, parameter :: geo_option_local = 1
  integer, parameter :: geo_option_inputprof = 2

  logical :: geoinit = .false.

contains

  subroutine init_geometry (nzed, nzgrid, zed, dz)

    use mp, only: proc0
    use millerlocal, only: read_local_parameters,  get_local_geo
!    use inputprofiles_interface, only: read_inputprof

    implicit none

    integer, intent (in) :: nzed, nzgrid
    real, dimension (-nzgrid:), intent (in) :: zed, dz

    real :: dpsidrho

    if (geoinit) return
    geoinit = .true.

    call allocate_arrays (nzgrid)

    if (proc0) then
       call read_parameters
       select case (geo_option_switch)
       case (geo_option_local)
          call read_local_parameters (qinp, shat)
          call get_local_geo (nzed, nzgrid, zed, &
               dpsidrho, grho, bmag, &
               gds2, gds21, gds22, gradpar, &
               gbdrift0, gbdrift, cvdrift0, cvdrift)
          drhodpsi = 1./dpsidrho
       case (geo_option_inputprof)
!          call read_inputprof
!          call get_local_geo (nzed, nzgrid, zed, &
!               dpsidrho, grho, bmag, &
!               gds2, gds21, gds22, gradpar, &
!               gbdrift0, gbdrift, cvdrift0, cvdrift)
          drhodpsi = 1./dpsidrho
       end select
    end if

    call broadcast_arrays

    jacob = 1.0/(drhodpsi*gradpar*bmag)
    
    dl_over_b = dz*jacob
    dl_over_b = dl_over_b / sum(dl_over_b)

    call get_dthet (nzgrid, dz, bmag, dbdthet)

  end subroutine init_geometry

  subroutine allocate_arrays (nzgrid)

    implicit none

    integer, intent (in) :: nzgrid

    if (.not.allocated(grho)) allocate (grho(-nzgrid:nzgrid))
    if (.not.allocated(bmag)) allocate (bmag(-nzgrid:nzgrid))
    if (.not.allocated(dbdthet)) allocate (dbdthet(-nzgrid:nzgrid))
    if (.not.allocated(jacob)) allocate (jacob(-nzgrid:nzgrid))
    if (.not.allocated(gradpar)) allocate (gradpar(-nzgrid:nzgrid))
    if (.not.allocated(dl_over_b)) allocate (dl_over_b(-nzgrid:nzgrid))
    if (.not.allocated(gds2)) allocate (gds2(-nzgrid:nzgrid))
    if (.not.allocated(gds21)) allocate (gds21(-nzgrid:nzgrid))
    if (.not.allocated(gds22)) allocate (gds22(-nzgrid:nzgrid))
    if (.not.allocated(gbdrift)) allocate (gbdrift(-nzgrid:nzgrid))
    if (.not.allocated(gbdrift0)) allocate (gbdrift0(-nzgrid:nzgrid))
    if (.not.allocated(cvdrift)) allocate (cvdrift(-nzgrid:nzgrid))
    if (.not.allocated(cvdrift0)) allocate (cvdrift0(-nzgrid:nzgrid))

  end subroutine allocate_arrays

  subroutine read_parameters

    use text_options
    use file_utils, only: error_unit, input_unit_exist

    implicit none

    integer :: in_file, ierr
    logical :: exist

    character (20) :: geo_option
    type (text_option), dimension (4), parameter :: geoopts = (/ &
         text_option('default', geo_option_local), &
         text_option('miller', geo_option_local), &
         text_option('local', geo_option_local), &
         text_option('input.profiles', geo_option_inputprof) /)

    namelist /geo_knobs/ geo_option

    geo_option = 'local'

    in_file = input_unit_exist("theta_grid_knobs", exist)
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
    call broadcast (grho)
    call broadcast (bmag)
    call broadcast (gradpar)
    call broadcast (gds2)
    call broadcast (gds21)
    call broadcast (gds22)
    call broadcast (gbdrift0)
    call broadcast (gbdrift)
    call broadcast (cvdrift0)
    call broadcast (cvdrift)

  end subroutine broadcast_arrays

  ! given function f(theta:-pi->pi), calculate theta derivative
  ! second order accurate, with equal grid spacing assumed
  ! assumes periodic in theta -- may need to change this in future
  subroutine get_dthet (nz, dz, f, df)

    implicit none

    integer, intent (in) :: nz
    real, dimension (-nz:), intent (in) :: dz, f
    real, dimension (-nz:), intent (out) :: df

    df(-nz+1:nz-1) = (f(-nz+2:)-f(:nz-2))/(dz(:nz-2)+dz(-nz+1:nz-1))

    ! use periodicity at boundary
    df(-nz) = (f(-nz+1)-f(nz-1))/(dz(-nz)+dz(nz-1))
    df(nz) = df(-nz)

  end subroutine get_dthet

  subroutine finish_geometry

    implicit none

    if (allocated(grho)) deallocate (grho)
    if (allocated(bmag)) deallocate (bmag)
    if (allocated(dbdthet)) deallocate (dbdthet)
    if (allocated(jacob)) deallocate (jacob)
    if (allocated(gradpar)) deallocate (gradpar)
    if (allocated(dl_over_b)) deallocate (dl_over_b)
    if (allocated(gds2)) deallocate (gds2)
    if (allocated(gds21)) deallocate (gds21)
    if (allocated(gds22)) deallocate (gds22)
    if (allocated(gbdrift)) deallocate (gbdrift)
    if (allocated(gbdrift0)) deallocate (gbdrift0)
    if (allocated(cvdrift)) deallocate (cvdrift)
    if (allocated(cvdrift0)) deallocate (cvdrift0)

    geoinit = .false.

  end subroutine finish_geometry

end module geometry
