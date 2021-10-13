module physics_flags

  implicit none

  public :: init_physics_flags
  public :: finish_physics_flags
  public :: full_flux_surface
  public :: radial_variation
  public :: include_parallel_nonlinearity
  public :: include_parallel_streaming
  public :: include_drifts
  public :: include_mirror
  public :: prp_shear_enabled
  public :: hammett_flow_shear
  public :: include_pressure_variation
  public :: include_geometric_variation
  public :: nonlinear
  public :: override_vexb
  public :: vexb_x
  public :: vexb_y

  private

  logical :: full_flux_surface
  logical :: radial_variation
  logical :: include_parallel_nonlinearity
  logical :: include_parallel_streaming
  logical :: include_drifts
  logical :: include_pressure_variation
  logical :: include_geometric_variation
  logical :: include_mirror
  logical :: nonlinear
  logical :: prp_shear_enabled
  logical :: hammett_flow_shear
  logical :: override_vexb
  real :: vexb_x
  real :: vexb_y
  logical :: initialized = .false.

contains

  subroutine init_physics_flags

    implicit none

    if (initialized) return
    initialized = .true.

    call read_parameters

  end subroutine init_physics_flags

  subroutine read_parameters

    use file_utils, only: input_unit_exist
    use mp, only: proc0, broadcast

    implicit none

    integer :: in_file
    logical :: rpexist

    namelist /physics_flags/ full_flux_surface, radial_variation, &
         include_parallel_nonlinearity, include_parallel_streaming, &
         include_mirror, nonlinear, &
         include_pressure_variation, include_geometric_variation, include_drifts, &
         override_vexb, vexb_x, vexb_y

    if (proc0) then
       full_flux_surface = .false.
       radial_variation = .false.
       include_pressure_variation = .true.
       include_geometric_variation = .true.
       include_parallel_nonlinearity = .false.
       include_parallel_streaming = .true.
       include_mirror = .true.
       include_drifts = .true.
       nonlinear = .false.
       override_vexb = .false.
       vexb_x = 0
       vexb_y = 0
       in_file = input_unit_exist("physics_flags", rpexist)
       if (rpexist) read (unit=in_file,nml=physics_flags)
    end if

    prp_shear_enabled = .false.
    hammett_flow_shear = .true.

    call broadcast (full_flux_surface)
    call broadcast (radial_variation)
    call broadcast (include_pressure_variation)
    call broadcast (include_geometric_variation)
    call broadcast (include_parallel_nonlinearity)
    call broadcast (include_parallel_streaming)
    call broadcast (include_mirror)
    call broadcast (include_drifts)
    call broadcast (nonlinear)
    call broadcast (override_vexb)
    call broadcast (vexb_x)
    call broadcast (vexb_y)

  end subroutine read_parameters

  subroutine finish_physics_flags

    implicit none

    initialized = .false.

  end subroutine finish_physics_flags

end module physics_flags
