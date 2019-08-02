!> This module is basically a store for the input parameters that are specified in the namelists \a knobs and \a parameters. In general, the names of the public variables in this module are the same as the name of the input parameter they correspond to.
 
module run_parameters

  implicit none

  public :: init_run_parameters, finish_run_parameters
  public :: fphi, fapar, fbpar
  public :: include_parallel_streaming
  public :: include_mirror
  public :: include_collisions
  public :: nonlinear
  public :: code_delt_max
  public :: nstep
  public :: cfl_cushion, delt_adjust
  public :: avail_cpu_time
  public :: stream_implicit, mirror_implicit
  public :: driftkinetic_implicit
  public :: fully_explicit
  public :: stream_cell, stream_matrix_inversion
  public :: mirror_semi_lagrange, mirror_linear_interp
  public :: zed_upwind, vpa_upwind, time_upwind
  public :: fields_kxkyz
  
  private

  real :: cfl_cushion, delt_adjust
  real :: fphi, fapar, fbpar
  real :: delt, code_delt_max
  real :: zed_upwind, vpa_upwind, time_upwind
  logical :: include_parallel_streaming
  logical :: include_mirror, include_collisions
  logical :: nonlinear
  logical :: stream_implicit, mirror_implicit
  logical :: driftkinetic_implicit
  logical :: fully_explicit
  logical :: stream_cell, stream_matrix_inversion
  logical :: mirror_semi_lagrange, mirror_linear_interp
  logical :: fields_kxkyz
  real :: avail_cpu_time
  integer :: nstep
  integer, public :: delt_option_switch
  integer, public, parameter :: delt_option_hand = 1, delt_option_auto = 2
  logical :: initialized = .false.
  logical :: knexist

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
         delt_option, nonlinear, &
         avail_cpu_time, cfl_cushion, delt_adjust, &
         stream_implicit, mirror_implicit, driftkinetic_implicit, &
         stream_cell, stream_matrix_inversion, &
         mirror_semi_lagrange, mirror_linear_interp, &
         include_parallel_streaming, include_mirror, &
         include_collisions, &
         zed_upwind, vpa_upwind, time_upwind, &
         fields_kxkyz

    if (proc0) then
       fphi = 1.0
       fapar = 1.0
       fbpar = -1.0
       fields_kxkyz = .false.
       stream_implicit = .true.
       mirror_implicit = .true.
       driftkinetic_implicit = .false.
       nonlinear = .false.
       include_parallel_streaming = .true.
       include_mirror = .true.
       include_collisions = .false.
       mirror_semi_lagrange = .true.
       mirror_linear_interp = .false.
       stream_cell = .true.
       stream_matrix_inversion = .false.
       delt_option = 'default'
       zed_upwind = 0.02
       vpa_upwind = 0.02
       time_upwind = 0.02
       avail_cpu_time = 1.e10
       cfl_cushion = 0.5
       delt_adjust = 2.0

       in_file = input_unit_exist("knobs", knexist)
       if (knexist) read (unit=in_file, nml=knobs)

       ierr = error_unit()
       call get_option_value &
            (delt_option, deltopts, delt_option_switch, ierr, &
            "delt_option in knobs")

    end if

    call broadcast (fields_kxkyz)
    call broadcast (delt_option_switch)
    call broadcast (delt)
    call broadcast (cfl_cushion)
    call broadcast (delt_adjust)
    call broadcast (fphi)
    call broadcast (fapar)
    call broadcast (fbpar)
    call broadcast (stream_implicit)
    call broadcast (mirror_implicit)
    call broadcast (driftkinetic_implicit)
    call broadcast (nonlinear)
    call broadcast (include_parallel_streaming)
    call broadcast (include_mirror)
    call broadcast (include_collisions)
    call broadcast (mirror_semi_lagrange)
    call broadcast (mirror_linear_interp)
    call broadcast (stream_cell)
    call broadcast (stream_matrix_inversion)
    call broadcast (zed_upwind)
    call broadcast (vpa_upwind)
    call broadcast (time_upwind)
    call broadcast (nstep)
    call broadcast (avail_cpu_time)
    
    if (.not.include_mirror) mirror_implicit = .false.

    code_delt_max = delt

    delt_saved = delt
    if (delt_option_switch == delt_option_auto) then
       call init_dt (delt_saved, istatus)
       if (istatus == 0) delt  = delt_saved
    endif

    if (driftkinetic_implicit) stream_implicit = .false.

    if (mirror_implicit .or. stream_implicit .or. driftkinetic_implicit) then
       fully_explicit = .false.
    else
       fully_explicit = .true.
    end if

  end subroutine read_parameters

  subroutine finish_run_parameters

    implicit none

    initialized = .false.

  end subroutine finish_run_parameters

end module run_parameters
