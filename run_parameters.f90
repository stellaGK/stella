!> This module is basically a store for the input parameters that are specified in the namelists \a knobs and \a parameters. In general, the names of the public variables in this module are the same as the name of the input parameter they correspond to.
 

module run_parameters

  implicit none

  public :: init_run_parameters, finish_run_parameters
  public :: fphi, fapar, fbpar
  public :: nonlinear
  public :: code_delt_max
  public :: nstep
  public :: cfl_cushion
  public :: avail_cpu_time
  
  private

  real :: cfl_cushion
  real :: fphi, fapar, fbpar
  real :: delt, code_delt_max
  logical :: nonlinear
  real :: avail_cpu_time
  integer :: nstep
  integer, public :: delt_option_switch
  integer, public, parameter :: delt_option_hand = 1, delt_option_auto = 2
  logical :: initialized = .false.
  logical :: rpexist, knexist

contains

  subroutine init_run_parameters

    use stella_time, only: init_delt
    
    implicit none

    if (initialized) return
    initialized = .true.

    call read_parameters

    call init_delt (delt)

  end subroutine init_run_parameters

  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use mp, only: proc0, broadcast
    use stella_save, only: init_dt
    use text_options, only: text_option, get_option_value
    implicit none
    type (text_option), dimension (3), parameter :: deltopts = &
         (/ text_option('default', delt_option_hand), &
            text_option('set_by_hand', delt_option_hand), &
            text_option('check_restart', delt_option_auto) /)
    character(20) :: delt_option

    integer :: ierr, istatus, in_file
    real :: delt_saved

    namelist /knobs/ fphi, fapar, fbpar, delt, nstep, &
         delt_option, &
         avail_cpu_time, nonlinear, cfl_cushion

    if (proc0) then
       fphi = 1.0
       fapar = 1.0
       fbpar = -1.0
       nonlinear = .false.
       delt_option = 'default'
       avail_cpu_time = 1.e10
       cfl_cushion = 0.5

       in_file = input_unit_exist("knobs", knexist)
       if (knexist) read (unit=in_file, nml=knobs)

       ierr = error_unit()
       call get_option_value &
            (delt_option, deltopts, delt_option_switch, ierr, &
            "delt_option in knobs")

    end if

    call broadcast (delt_option_switch)
    call broadcast (delt)
    call broadcast (cfl_cushion)
    call broadcast (fphi)
    call broadcast (fapar)
    call broadcast (fbpar)
    call broadcast (nonlinear)
    call broadcast (nstep)
    call broadcast (avail_cpu_time)
    
    code_delt_max = delt

    delt_saved = delt
    if (delt_option_switch == delt_option_auto) then
       call init_dt (delt_saved, istatus)
       if (istatus == 0) delt  = delt_saved
    endif

  end subroutine read_parameters

  subroutine finish_run_parameters

    implicit none

    initialized = .false.

  end subroutine finish_run_parameters

end module run_parameters
