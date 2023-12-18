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
   public :: quasilinear_ExB
   public :: remove_phiZ_NL
   public :: remove_gZ_NL
   public :: remove_gZ_even_NL
   public :: remove_gZ_odd_NL
   public :: scale_NL_term
   public :: maxwellianize_gZ_NL
   public :: remove_uparZ_NL
   public :: remove_tempZ_NL
   public :: scale_parallel_streaming_zonal
   public :: adiabatic_option_switch
   public :: adiabatic_option_fieldlineavg
   public :: const_alpha_geo
   public :: suppress_zonal_interaction
   public :: suppress_zonal_kmin
   public :: suppress_zonal_kmax

   private

   logical :: full_flux_surface
   logical :: radial_variation
   logical :: include_parallel_nonlinearity
   logical :: include_parallel_streaming
   logical :: include_pressure_variation
   logical :: include_geometric_variation
   logical :: include_mirror
   logical :: nonlinear
   logical :: quasilinear_ExB
   logical :: remove_phiZ_NL
   logical :: remove_gZ_NL
   logical :: remove_gZ_even_NL
   logical :: remove_gZ_odd_NL
   logical :: maxwellianize_gZ_NL
   logical :: remove_uparZ_NL
   logical :: remove_tempZ_NL
   logical :: prp_shear_enabled
   logical :: hammett_flow_shear
   logical :: const_alpha_geo
   logical :: suppress_zonal_interaction

   integer :: adiabatic_option_switch
   integer :: scale_NL_term
   real :: scale_parallel_streaming_zonal
   real :: suppress_zonal_kmin
   real :: suppress_zonal_kmax
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
         include_mirror, nonlinear, quasilinear_ExB, scale_parallel_streaming_zonal, &
         remove_gZ_NL, remove_phiZ_NL, remove_gZ_even_NL, remove_gZ_odd_NL, scale_NL_term,&
         maxwellianize_gZ_NL, remove_uparZ_NL, remove_tempZ_NL,&
         include_pressure_variation, include_geometric_variation, &
         adiabatic_option, const_alpha_geo, &
         suppress_zonal_interaction, suppress_zonal_kmin, suppress_zonal_kmax

      if (proc0) then
         full_flux_surface = .false.
         radial_variation = .false.
         include_pressure_variation = .true.
         include_geometric_variation = .true.
         include_parallel_nonlinearity = .false.
         include_parallel_streaming = .true.
         include_mirror = .true.
         suppress_zonal_interaction= .false.
         suppress_zonal_kmin = -1
         suppress_zonal_kmax = 1e10
         nonlinear = .false.
         quasilinear_ExB = .false.
         remove_gZ_NL    = .false.
         remove_phiZ_NL    = .false.
         remove_gZ_even_NL = .false.
         remove_gZ_odd_NL  = .false.
         scale_NL_term = 1.0
         maxwellianize_gZ_NL  = .false.
         remove_uparZ_NL = .false.
         remove_tempZ_NL = .false.
         scale_parallel_streaming_zonal = 1.0
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
      call broadcast(quasilinear_ExB)
      call broadcast(remove_gZ_NL)
      call broadcast(remove_phiZ_NL)
      call broadcast(remove_gZ_even_NL)
      call broadcast(remove_gZ_odd_NL)
      call broadcast(scale_NL_term)
      call broadcast(maxwellianize_gZ_NL)
      call broadcast(remove_uparZ_NL)
      call broadcast(remove_tempZ_NL)
      call broadcast(scale_parallel_streaming_zonal)
      call broadcast(adiabatic_option_switch)
      call broadcast(const_alpha_geo)
      call broadcast(suppress_zonal_interaction)
      call broadcast(suppress_zonal_kmin)
      call broadcast(suppress_zonal_kmax)

   end subroutine read_parameters

   subroutine finish_physics_flags

      implicit none

      initialized = .false.

   end subroutine finish_physics_flags

end module physics_flags
