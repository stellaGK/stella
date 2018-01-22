module physics_parameters

  implicit none

  public :: init_physics_parameters
  public :: finish_physics_parameters
  public :: beta, zeff, tite, nine, rhostar, vnew_ref

  private

  real :: beta, zeff, tite, nine, rhostar, vnew_ref
  logical :: initialized = .false.

contains

  subroutine init_physics_parameters

    implicit none

    if (initialized) return
    initialized = .true.

    call read_parameters

  end subroutine init_physics_parameters

  subroutine read_parameters

    use file_utils, only: input_unit_exist
    use mp, only: proc0, broadcast

    implicit none

    integer :: in_file
    logical :: rpexist
    
    namelist /parameters/ beta, zeff, tite, nine, rhostar, vnew_ref

    if (proc0) then
       beta = 0.0
       vnew_ref = 0.01
       rhostar = 3.0e-3 ! = m_ref * vt_ref / (e * B_ref * a_ref), with refs in SI
       zeff = 1.0
       tite = 1.0
       nine = 1.0

       in_file = input_unit_exist("parameters", rpexist)
       if (rpexist) read (unit=in_file,nml=parameters)
    end if

    call broadcast (beta)
    call broadcast (vnew_ref)
    call broadcast (zeff)
    call broadcast (rhostar)
    call broadcast (tite)
    call broadcast (nine)

  end subroutine read_parameters

  subroutine finish_physics_parameters

    implicit none

    initialized = .false.

  end subroutine finish_physics_parameters

end module physics_parameters
