module dissipation

  implicit none

  public :: hyper_dissipation
  public :: init_dissipation
  public :: advance_hyper_dissipation

  private

  logical :: hyper_dissipation
  real :: D_hyper

contains

  subroutine init_dissipation

    implicit none

    call read_parameters

  end subroutine init_dissipation

  subroutine read_parameters

    use file_utils, only: input_unit_exist
    use mp, only: proc0, broadcast

    implicit none

    namelist /dissipation/ hyper_dissipation, D_hyper

    integer :: in_file
    logical :: dexist

    if (proc0) then
       hyper_dissipation = .false.
       D_hyper = 0.05

       in_file = input_unit_exist("dissipation", dexist)
       if (dexist) read (unit=in_file, nml=dissipation)
    end if

    call broadcast (hyper_dissipation)
    call broadcast (D_hyper)

  end subroutine read_parameters

  subroutine advance_hyper_dissipation (g)

    use stella_time, only: code_dt
    use zgrid, only: nzgrid
    use stella_layouts, only: vmu_lo
    use dist_fn_arrays, only: kperp2

    implicit none

    complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: g

    integer :: ivmu
    real :: k2max

    k2max = maxval(kperp2)

    ! add in hyper-dissipation of form dg/dt = -D*(k/kmax)^4*g
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       g(:,:,:,ivmu) = g(:,:,:,ivmu)/(1.+code_dt*(kperp2/k2max)**2*D_hyper)
    end do

  end subroutine advance_hyper_dissipation

end module dissipation
