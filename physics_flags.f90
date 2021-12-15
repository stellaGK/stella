module physics_flags

   implicit none

   public :: init_physics_flags
   public :: finish_physics_flags
   public :: full_flux_surface
   public :: radial_variation
   public :: include_parallel_nonlinearity
   public :: include_parallel_streaming
   public :: include_mirror
   public :: prp_shear_enabled
   public :: hammett_flow_shear
   public :: include_pressure_variation
   public :: include_geometric_variation
   public :: nonlinear
   public :: adiabatic_option_switch
   public :: adiabatic_option_fieldlineavg
   public :: const_alpha_geo

   private

   logical :: full_flux_surface
   logical :: radial_variation
   logical :: include_parallel_nonlinearity
   logical :: include_parallel_streaming
   logical :: include_pressure_variation
   logical :: include_geometric_variation
   logical :: include_mirror
   logical :: nonlinear
   logical :: prp_shear_enabled
   logical :: hammett_flow_shear
   logical :: const_alpha_geo

   integer :: adiabatic_option_switch
   integer, parameter :: adiabatic_option_default = 1, &
                         adiabatic_option_zero = 2, &
                         adiabatic_option_fieldlineavg = 3

   logical :: initialized = .false.

contains

   subroutine init_physics_flags

      implicit none

      if (initialized) return
      initialized = .true.

      call read_parameters

   end subroutine init_physics_flags

   subroutine read_parameters

      use file_utils, only: input_unit_exist, error_unit
      use mp, only: proc0, broadcast
      use text_options, only: text_option, get_option_value

      implicit none

      integer :: in_file, ierr
      logical :: rpexist

      type(text_option), dimension(6), parameter :: adiabaticopts = &
                                                    (/text_option('default', adiabatic_option_default), &
                                                      text_option('no-field-line-average-term', adiabatic_option_default), &
                                                      text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
                                                      text_option('iphi00=0', adiabatic_option_default), &
                                                      text_option('iphi00=1', adiabatic_option_default), &
                                                      text_option('iphi00=2', adiabatic_option_fieldlineavg)/)
      character(30) :: adiabatic_option

      namelist /physics_flags/ full_flux_surface, radial_variation, &
         include_parallel_nonlinearity, include_parallel_streaming, &
         include_mirror, nonlinear, &
         include_pressure_variation, include_geometric_variation, &
         adiabatic_option, const_alpha_geo

      if (proc0) then
         full_flux_surface = .false.
         radial_variation = .false.
         include_pressure_variation = .true.
         include_geometric_variation = .true.
         include_parallel_nonlinearity = .false.
         include_parallel_streaming = .true.
         include_mirror = .true.
         nonlinear = .false.
         adiabatic_option = 'default'
         const_alpha_geo = .false.

         in_file = input_unit_exist("physics_flags", rpexist)
         if (rpexist) read (unit=in_file, nml=physics_flags)

         ierr = error_unit()
         call get_option_value &
            (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
             ierr, "adiabatic_option in physics_flags")
      end if

      prp_shear_enabled = .false.
      hammett_flow_shear = .true.

      call broadcast(full_flux_surface)
      call broadcast(radial_variation)
      call broadcast(include_pressure_variation)
      call broadcast(include_geometric_variation)
      call broadcast(include_parallel_nonlinearity)
      call broadcast(include_parallel_streaming)
      call broadcast(include_mirror)
      call broadcast(nonlinear)
      call broadcast(adiabatic_option_switch)
      call broadcast(const_alpha_geo)

   end subroutine read_parameters

   subroutine finish_physics_flags

      implicit none

      initialized = .false.

   end subroutine finish_physics_flags

end module physics_flags
