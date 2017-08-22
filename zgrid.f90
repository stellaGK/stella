module zgrid

  implicit none

  public :: init_theta_grid, finish_theta_grid
  public :: ntheta, ntgrid, nperiod
  public :: theta
  public :: delthet
  public :: shat_zero

  private

  integer :: ntheta, ntgrid, nperiod
  real :: shat_zero
  real, dimension (:), allocatable :: theta, delthet

  logical :: zgridinit = .false.

contains

  subroutine init_theta_grid

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

    if (.not.allocated(theta)) allocate (theta(-ntgrid:ntgrid))
    if (.not.allocated(delthet)) allocate (delthet(-ntgrid:ntgrid))

    theta = (/ (i*pi/real(ntheta/2), i=-ntgrid, ntgrid ) /)
    delthet(:ntgrid-1) = theta(-ntgrid+1:) - theta(:ntgrid-1)
    delthet(ntgrid) = delthet(-ntgrid)

  end subroutine init_theta_grid

  subroutine read_parameters

    use file_utils, only: input_unit_exist

    implicit none

    integer :: in_file
    logical :: exist

    namelist /zgrid_parameters/ ntheta, nperiod, shat_zero

    ntheta = 32
    nperiod = 1
    ! set minimum shat value below which we assume
    ! periodic BC
    shat_zero = 1.e-5

    in_file = input_unit_exist("zgrid_parameters", exist)
    if (exist) read (unit=in_file, nml=zgrid_parameters)

    ntgrid = ntheta/2 + (nperiod-1)*ntheta

  end subroutine read_parameters

  subroutine broadcast_parameters

    use mp, only: broadcast

    implicit none

    call broadcast (ntheta)
    call broadcast (ntgrid)
    call broadcast (nperiod)
    call broadcast (shat_zero)

  end subroutine broadcast_parameters

  subroutine finish_theta_grid

    implicit none

    if (allocated(theta)) deallocate (theta)
    if (allocated(delthet)) deallocate (delthet)

    zgridinit = .false.

  end subroutine finish_theta_grid

end module zgrid
