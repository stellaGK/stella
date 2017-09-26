module zgrid

  implicit none

  public :: init_zgrid, finish_zgrid
  public :: nzed, nzgrid, nperiod
  public :: nztot
  public :: zed
  public :: delzed
  public :: shat_zero

  private

  integer :: nzed, nzgrid, nperiod, nztot
  real :: shat_zero
  real, dimension (:), allocatable :: zed, delzed

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

  end subroutine init_zgrid

  subroutine read_parameters

    use file_utils, only: input_unit_exist

    implicit none

    integer :: in_file
    logical :: exist

    namelist /zgrid_parameters/ nzed, nperiod, shat_zero

    nzed = 32
    nperiod = 1
    ! set minimum shat value below which we assume
    ! periodic BC
    shat_zero = 1.e-5

    in_file = input_unit_exist("zgrid_parameters", exist)
    if (exist) read (unit=in_file, nml=zgrid_parameters)

    nzgrid = nzed/2 + (nperiod-1)*nzed

  end subroutine read_parameters

  subroutine broadcast_parameters

    use mp, only: broadcast

    implicit none

    call broadcast (nzed)
    call broadcast (nzgrid)
    call broadcast (nperiod)
    call broadcast (shat_zero)

  end subroutine broadcast_parameters

  subroutine finish_zgrid

    implicit none

    if (allocated(zed)) deallocate (zed)
    if (allocated(delzed)) deallocate (delzed)

    zgridinit = .false.

  end subroutine finish_zgrid

end module zgrid
