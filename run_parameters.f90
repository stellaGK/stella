!> This module is basically a store for the input parameters that are specified in the namelists \a knobs and \a parameters. In general, the names of the public variables in this module are the same as the name of the input parameter they correspond to.

module run_parameters

   implicit none

   public :: init_run_parameters, finish_run_parameters
   public :: fphi, fapar, fbpar
   public :: nstep, tend, delt
   public :: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
   public :: delt_max, delt_min
   public :: avail_cpu_time
   public :: stream_implicit, mirror_implicit
   public :: drifts_implicit
   public :: driftkinetic_implicit
   public :: fully_explicit, fully_implicit
   public :: ky_solve_radial, ky_solve_real
   public :: maxwellian_inside_zed_derivative
   public :: stream_matrix_inversion
   public :: mirror_semi_lagrange, mirror_linear_interp
   public :: zed_upwind, vpa_upwind, time_upwind
   public :: fields_kxkyz, mat_gen, mat_read
   public :: rng_seed
   public :: use_deltaphi_for_response_matrix
   public :: maxwellian_normalization
   public :: time_upwind_plus, time_upwind_minus
   public :: zed_upwind_plus, zed_upwind_minus

   private

   real :: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
   real :: fphi, fapar, fbpar
   real :: delt, tend, delt_max, delt_min
   real :: vpa_upwind
   real :: time_upwind, time_upwind_plus, time_upwind_minus
   real :: zed_upwind, zed_upwind_plus, zed_upwind_minus
   logical :: stream_implicit, mirror_implicit, drifts_implicit
   logical :: driftkinetic_implicit
   logical :: fully_explicit, fully_implicit
   logical :: maxwellian_inside_zed_derivative
   logical :: stream_matrix_inversion
   logical :: mirror_semi_lagrange, mirror_linear_interp
   logical :: fields_kxkyz, mat_gen, mat_read
   logical :: ky_solve_real
   logical :: use_deltaphi_for_response_matrix
   logical :: maxwellian_normalization
   real :: avail_cpu_time
   integer :: nstep, ky_solve_radial
   integer :: rng_seed
   integer, public :: delt_option_switch, lu_option_switch
   integer, public, parameter :: delt_option_hand = 1, delt_option_auto = 2
   integer, public, parameter :: lu_option_none = 1, &
                                 lu_option_local = 2, &
                                 lu_option_global = 3
   logical :: initialized = .false.
   logical :: knexist

contains

   subroutine init_run_parameters

      implicit none

      if (initialized) return
      initialized = .true.

      call read_parameters

   end subroutine init_run_parameters

   subroutine read_parameters

      use file_utils, only: input_unit, error_unit, input_unit_exist
      use mp, only: mp_abort, proc0, broadcast
      use text_options, only: text_option, get_option_value
      use physics_flags, only: include_mirror, full_flux_surface, radial_variation, include_parallel_streaming
      use physics_flags, only: nonlinear

      implicit none

      type(text_option), dimension(3), parameter :: deltopts = &
                                                    (/text_option('default', delt_option_auto), &
                                                      text_option('set_by_hand', delt_option_hand), &
                                                      text_option('check_restart', delt_option_auto)/)
      type(text_option), dimension(4), parameter :: lu_opts = &
                                                    (/text_option('default', lu_option_none), &
                                                      text_option('none', lu_option_none), &
                                                      text_option('local', lu_option_local), &
                                                      text_option('global', lu_option_global)/)

      character(20) :: delt_option, lu_option
      logical :: error = .false.
      integer :: ierr, in_file

      namelist /knobs/ fphi, fapar, fbpar, delt, nstep, tend, &
         delt_option, lu_option, &
         avail_cpu_time, delt_max, delt_min, &
         cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower, &
         stream_implicit, mirror_implicit, driftkinetic_implicit, &
         drifts_implicit, use_deltaphi_for_response_matrix, &
         maxwellian_normalization, &
         stream_matrix_inversion, maxwellian_inside_zed_derivative, &
         mirror_semi_lagrange, mirror_linear_interp, &
         zed_upwind, vpa_upwind, time_upwind, &
         fields_kxkyz, mat_gen, mat_read, rng_seed, &
         ky_solve_radial, ky_solve_real

      if (proc0) then

         ! Default parameters in namelist <knobs>
         fphi = 1.0
         fapar = 1.0
         fbpar = -1.0
         fields_kxkyz = .false.
         stream_implicit = .true.
         mirror_implicit = .true.
         drifts_implicit = .true.
         driftkinetic_implicit = .false.
         maxwellian_inside_zed_derivative = .false.
         mirror_semi_lagrange = .true.
         mirror_linear_interp = .false.
         stream_matrix_inversion = .false.
         use_deltaphi_for_response_matrix = .false.
         maxwellian_normalization = .false.
         delt_option = 'default'
         zed_upwind = 0.02
         vpa_upwind = 0.02
         time_upwind = 0.02
         rng_seed = -1 !negative values use current time as seed
         ky_solve_radial = 0
         ky_solve_real = .false.
         mat_gen = .false.
         mat_read = .false.

         ! Stella runs until t*v_{th,i}/a=tend or until istep=nstep
         tend = -1.0
         nstep = -1

         ! Set the available wall time in seconds, 5 minutes before the wall, stella will make a clean exit
         avail_cpu_time = 1.e10

         ! code_dt needs to stay within [cfl_dt*cfl_cushion_upper, cfl_dt*cfl_cushion_lower]
         ! code_dt can be increased if cfl_dt increases, however, never increase above delt_max (=delt by default)
         ! Exit stella if code_dt < delt_min (e.g. when the code blows up)
         cfl_cushion_upper = 0.5       ! Stay a factor of 2 under the CFL condition, otherwise it might run out of control
         cfl_cushion_middle = 0.25     ! If code_dt>cfl_dt/2 or code_dt<cfl_dt/100000, set code_dt to cfl_dt/4
         cfl_cushion_lower = 0.00001   ! Default is very low to not trigger it.
         delt_max = -1
         delt_min = 1.e-10

         ! The response matrix is solved with a none, local or global scheme, local seems to be the most efficient
         lu_option = 'default'

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist("knobs", knexist)
         if (knexist) read (unit=in_file, nml=knobs)

         ierr = error_unit()
         call get_option_value &
            (delt_option, deltopts, delt_option_switch, ierr, &
             "delt_option in knobs")
         call get_option_value &
            (lu_option, lu_opts, lu_option_switch, ierr, &
             "lu_option in knobs")

         ! Abort if neither tend nor nstep are set
         if (tend < 0 .and. nstep < 0) then
            ierr = error_unit()
            write (ierr, *) ''
            write (ierr, *) 'Please specify either <nstep> or <tend> in the <knobs> namelist.'
            write (ierr, *) 'Aborting.'
            write (*, *) ''
            write (*, *) 'Please specify either <nstep> or <tend> in the <knobs> namelist.'
            write (*, *) 'Aborting.'
            error = .true.
         end if

         ! Abort if cfl_cushion_lower>cfl_cushion_upper or if cfl_cushion_lower==cfl_cushion_upper
         if ((cfl_cushion_lower > cfl_cushion_upper - 0.001) &
             .or. (cfl_cushion_middle > cfl_cushion_upper - 0.001) &
             .or. (cfl_cushion_middle < cfl_cushion_lower + 0.001)) then
            ierr = error_unit()
            write (ierr, *) ''
            write (ierr, *) 'Please make sure that <cfl_cushion_upper> is bigger than <cfl_cushion_lower>,'
            write (ierr, *) 'and that <cfl_cushion_middle> lies in between <cfl_cushion_upper> and <cfl_cushion_lower>.'
            write (ierr, *) 'Aborting.'
            write (*, *) ''
            write (*, *) 'Please make sure that <cfl_cushion_upper> is bigger than <cfl_cushion_lower>,'
            write (*, *) 'and that <cfl_cushion_middle> lies in between <cfl_cushion_upper> and <cfl_cushion_lower>.'
            write (*, *) 'Aborting.'
            error = .true.
         end if

         if (radial_variation .and. maxwellian_normalization) then
            write (*, *) '!!!WARNING!!!'
            write (*, *) 'maxwellian_normalization is not currently supported for use with radial_variation.'
            write (*, *) 'forcing maxwellian_normalization = F.'
            write (*, *) '!!!WARNING!!!'
            maxwellian_normalization = .false.
         end if

         if (maxwellian_normalization .and. mirror_semi_lagrange) then
            write (*, *) '!!!WARNING!!!'
            write (*, *) 'maxwellian_normalization is not consistent with mirror_semi_lagrange = T.'
            write (*, *) 'forcing mirror_semi_lagrange = F.'
            write (*, *) '!!!WARNING!!!'
            mirror_semi_lagrange = .false.
         end if

         if (drifts_implicit) then
            if (.not. stream_implicit) then
               write (*, *) '!!!WARNING!!!'
               write (*, *) 'drifts_implicit = T requires stream_implicit = T.'
               write (*, *) 'forcing drifts_implicit = F.'
               drifts_implicit = .false.
            else if (.not. include_parallel_streaming) then
               write (*, *) '!!!WARNING!!!'
               write (*, *) 'drifts_implicit = T requires include_parallel_streaming = T.'
               write (*, *) 'forcing drifts_implicit = F.'
               drifts_implicit = .false.
            end if
         end if
         
      end if

      ! Exit stella if we ran into an error
      call broadcast(error)
      if (error) call mp_abort('Aborting in run_parameters.f90')

      call broadcast(fields_kxkyz)
      call broadcast(delt_option_switch)
      call broadcast(delt)
      call broadcast(lu_option_switch)
      call broadcast(cfl_cushion_upper)
      call broadcast(cfl_cushion_middle)
      call broadcast(cfl_cushion_lower)
      call broadcast(delt_max)
      call broadcast(delt_min)
      call broadcast(fphi)
      call broadcast(fapar)
      call broadcast(fbpar)
      call broadcast(stream_implicit)
      call broadcast(mirror_implicit)
      call broadcast(drifts_implicit)
      call broadcast(driftkinetic_implicit)
      call broadcast(maxwellian_inside_zed_derivative)
      call broadcast(mirror_semi_lagrange)
      call broadcast(mirror_linear_interp)
      call broadcast(stream_matrix_inversion)
      call broadcast(use_deltaphi_for_response_matrix)
      call broadcast(maxwellian_normalization)
      call broadcast(zed_upwind)
      call broadcast(vpa_upwind)
      call broadcast(time_upwind)
      call broadcast(nstep)
      call broadcast(tend)
      call broadcast(avail_cpu_time)
      call broadcast(rng_seed)
      call broadcast(ky_solve_radial)
      call broadcast(ky_solve_real)
      call broadcast(mat_gen)
      call broadcast(mat_read)

      ! calculate some useful derived quantities that are used repeatedly across modules
      time_upwind_plus = 0.5 * (1.0 + time_upwind)
      time_upwind_minus = 0.5 * (1.0 - time_upwind)
      zed_upwind_plus = 0.5 * (1.0 + zed_upwind)
      zed_upwind_minus = 0.5 * (1.0 - zed_upwind)

      if (.not. include_mirror) mirror_implicit = .false.

      if (driftkinetic_implicit) then
         stream_implicit = .false.
      else if (stream_implicit .and. full_flux_surface) then
         stream_implicit = .false.
         write (*, *)
         write (*, *) "!!!WARNING!!!"
         write (*, *) "The option stream_implicit=T is not supported for full_flux_surface=T."
         write (*, *) "Setting driftkinetic_implicit=T instead."
         write (*, *) "!!!WARNING!!!"
         write (*, *)
         driftkinetic_implicit = .true.
      end if

      if (mirror_implicit .or. stream_implicit .or. driftkinetic_implicit .or. drifts_implicit) then
         fully_explicit = .false.
      else
         fully_explicit = .true.
      end if

      if (mirror_implicit .and. stream_implicit .and. drifts_implicit .and. .not. nonlinear) then
         fully_implicit = .true.
      else
         fully_implicit = .false.
      end if

      !> print warning messages and override inconsistent or unsupported options for full_flux_surface = T
      if (full_flux_surface) then
         if (fields_kxkyz) then
            write (*, *)
            write (*, *) '!!!WARNING!!!'
            write (*, *) 'The option fields_kxkyz=T is not currently supported for full_flux_surface=T.'
            write (*, *) 'Forcing fields_kxkyz=F.'
            write (*, *) '!!!WARNING!!!'
            write (*, *)
            fields_kxkyz = .false.
         end if
         if (mirror_semi_lagrange) then
            write (*, *)
            write (*, *) '!!!WARNING!!!'
            write (*, *) 'The option mirror_semi_lagrange=T is not consistent with full_flux_surface=T.'
            write (*, *) 'Forcing mirror_semi_lagrange=F.'
            write (*, *) '!!!WARNING!!!'
            mirror_semi_lagrange = .false.
         end if
         ! the full flux surface implementation relies on the use of gnorm = g / F_Maxwellian
         ! as the evolved pdf
         if (full_flux_surface) then
            maxwellian_normalization = .true.
         end if
      end if

   end subroutine read_parameters

   subroutine finish_run_parameters

      implicit none

      initialized = .false.

   end subroutine finish_run_parameters

end module run_parameters
