module zgrid

  implicit none

  public :: init_zgrid, finish_zgrid
  public :: nzed, nzgrid, nperiod
  public :: nztot, nz2pi
  public :: zed
  public :: delzed
  public :: shat_zero
  public :: boundary_option_switch
  public :: boundary_option_zero
  public :: boundary_option_self_periodic
  public :: boundary_option_linked

  private

  integer :: nzed, nzgrid, nperiod, nztot, nz2pi
  real :: shat_zero
  real, dimension (:), allocatable :: zed, delzed

  integer :: boundary_option_switch
  integer, parameter :: boundary_option_zero = 1, &
       boundary_option_self_periodic = 2, &
       boundary_option_linked = 3

  logical :: zgridinit = .false.

contains

  subroutine init_zgrid

    use mp, only: proc0
    use constants, only: pi

    implicit none

    integer :: i

    if (zgridinit) return
    zgridinit = .true.

    if (proc0) then
       call read_parameters
    end if
    call broadcast_parameters

    if (.not.allocated(zed)) allocate (zed(-nzgrid:nzgrid))
    if (.not.allocated(delzed)) allocate (delzed(-nzgrid:nzgrid))

    zed = (/ (i*pi/real(nzed/2), i=-nzgrid, nzgrid ) /)
    delzed(:nzgrid-1) = zed(-nzgrid+1:) - zed(:nzgrid-1)
    delzed(nzgrid) = delzed(-nzgrid)

    nztot = 2*nzgrid+1
    nz2pi = 2*(nzed/2)+1

  end subroutine init_zgrid

  subroutine read_parameters

    use file_utils, only: input_unit_exist, error_unit
    use text_options, only: text_option, get_option_value

    implicit none

    integer :: in_file, ierr
    logical :: exist

    type (text_option), dimension (6), parameter :: boundaryopts = &
         (/ text_option('default', boundary_option_zero), &
            text_option('zero', boundary_option_zero), &
            text_option('unconnected', boundary_option_zero), &
            text_option('self-periodic', boundary_option_self_periodic), &
            text_option('periodic', boundary_option_self_periodic), &
            text_option('linked', boundary_option_linked) /)
    character(20) :: boundary_option

    namelist /zgrid_parameters/ nzed, nperiod, shat_zero, boundary_option

    nzed = 24
    nperiod = 1
    boundary_option = 'default'
    ! set minimum shat value below which we assume
    ! periodic BC
    shat_zero = 1.e-5

    in_file = input_unit_exist("zgrid_parameters", exist)
    if (exist) read (unit=in_file, nml=zgrid_parameters)

    ierr = error_unit()
    call get_option_value &
         (boundary_option, boundaryopts, boundary_option_switch, &
         ierr, "boundary_option in dist_fn_knobs")

    ! note that boundary_option may be changed to self-periodic later
    ! if magnetic shear is smaller than shat_zero
    ! cannot do this here as magnetic shear has yet to be input

    nzgrid = nzed/2 + (nperiod-1)*nzed

  end subroutine read_parameters

  subroutine broadcast_parameters

    use mp, only: broadcast

    implicit none

    call broadcast (nzed)
    call broadcast (nzgrid)
    call broadcast (nperiod)
    call broadcast (shat_zero)
    call broadcast (boundary_option_switch)

  end subroutine broadcast_parameters

  subroutine finish_zgrid

    implicit none

    if (allocated(zed)) deallocate (zed)
    if (allocated(delzed)) deallocate (delzed)

    zgridinit = .false.

  end subroutine finish_zgrid

end module zgrid
