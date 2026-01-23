!###############################################################################
!########################## CALCULATIONS - TIME STEP ###########################
!###############################################################################
! 
! This module is used for CFL calculations, and resetting the time step, dt.
! 
!###############################################################################
module calculations_timestep

   ! Load debug flags
   use debug_flags, only: debug => time_advance_debug

   implicit none

   ! Make routines accesible to other modules
   public :: init_cfl
   public :: reset_dt

   private
    
contains
    
!###############################################################################
!########################### INITIALISE CFL CONDITION ##########################
!###############################################################################
   subroutine init_cfl

      ! Parallelisation
      use mp, only: proc0, nproc, max_allreduce, min_allreduce
      use mp, only: scope, allprocs, subprocs
      use debug_flags, only: print_extra_info_to_terminal
      
      ! Grids
      use grids_kxky, only: nx
      use grids_z, only: delzed
      use grids_velocity, only: dvpa
      use grids_kxky, only: akx, aky, rho
      
      ! CFL parameters
      use grids_time, only: code_dt, write_dt, cfl_dt_linear
      use parameters_numerical, only: cfl_cushion_upper
      use parameters_numerical, only: cfl_cushion_middle
      use parameters_numerical, only: cfl_cushion_lower
            
      ! Gyrokinetic equation 
      use parameters_numerical, only: stream_implicit
      use parameters_numerical, only: mirror_implicit
      use parameters_numerical, only: drifts_implicit
      use parameters_physics, only: include_drive
      use parameters_physics, only: include_xdrift, include_ydrift
      use gk_parallel_streaming, only: stream
      use gk_parallel_streaming, only: stream_rad_var1
      use gk_parallel_streaming, only: stream_rad_var2
      use gk_mirror, only: mirror
      
      ! Arrays.
      use arrays, only: wdriftx_g, wdrifty_g
      use arrays, only: wstar
      
      ! NEO's neoclassical corrections.
      use arrays, only: neo_chi_coeff, neo_apar_coeff, neo_dchidz_coeff
      use arrays, only: wstar1, wpol
      ! use arrays, only: neomagx, neomagy
      use arrays, only: neocurvx, neocurvy
      
      ! Collisions
      use dissipation_and_collisions, only: include_collisions, collisions_implicit
      use dissipation_and_collisions, only: cfl_dt_vpadiff, cfl_dt_mudiff
      
      ! Radial variation - multibox - flow shear
      use parameters_physics, only: radial_variation
      use gk_flow_shear, only: prl_shear, shift_times
      use gk_flow_shear, only: prp_shear_enabled
      use file_utils, only: runtype_option_switch, runtype_multibox

      ! For NEO's neoclassical corrections.
      use neoclassical_terms_neo, only: neoclassical_is_enabled

      implicit none

      ! Local variables
      real :: cfl_dt_mirror, cfl_dt_stream, cfl_dt_shear
      real :: cfl_dt_wdriftx, cfl_dt_wdrifty, cfl_dt_wstar
      real :: wdriftx_max, wdrifty_max, wstar_max, zero
      
      ! For NEO's neoclassical corrections. 

      real :: cfl_dt_neo_chi_coeff, cfl_dt_neo_apar_coeff, cfl_dt_neo_dchidz_coeff
      real :: cfl_dt_wstar1, cfl_dt_wpol
      ! real :: cfl_dt_neomagx, cfl_dt_neomagy
      real :: cfl_dt_neocurvx, cfl_dt_neocurvy
      real :: neo_chi_coeff_max, neo_apar_coeff_max, neo_dchidz_coeff_max
      real :: wstar1_max, wpol_max
      ! real :: neomagx_max, neomagy_max
      real :: neocurvx_max, neocurvy_max      


      !-------------------------------------------------------------------------

      ! Avoid divide by zero in cfl_dt terms below
      zero = 100.*epsilon(0.)

      ! FLAG -- assuming equal spacing in zed!
      
      if (cfl_dt_linear < 0) cfl_dt_linear = code_dt / cfl_cushion_upper

      if (.not. drifts_implicit .and. include_xdrift) then
      
         ! Get the local max value of wdriftx on each processor
         wdriftx_max = maxval(abs(wdriftx_g))
         ! Compare these max values across processors to get global max
         if (nproc > 1) then
            call max_allreduce(wdriftx_max)
         end if
         
         ! NB: wdriftx_g has code_dt built-in, which accounts for code_dt factor here
         cfl_dt_wdriftx = abs(code_dt) / max(maxval(abs(akx)) * wdriftx_max, zero)
         cfl_dt_linear = cfl_dt_wdriftx
         
      end if

      cfl_dt_shear = abs(code_dt) / max(maxval(abs(aky)) * maxval(abs(prl_shear)), zero)
      cfl_dt_linear = min(cfl_dt_linear, cfl_dt_shear)

      if (prp_shear_enabled) then
         cfl_dt_shear = minval(shift_times)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_shear)
      end if

      ! Stream has code_dt built-in, which accounts for code_dt factor here
      if (.not. stream_implicit) then
         cfl_dt_stream = abs(code_dt) * delzed(0) / max(maxval(abs(stream)), zero)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_stream)
      end if

      ! TODO:GA- add correct CFL condition
      ! if (driftkinetic_implicit) then
      !    cfl_dt_stream = abs(code_dt) * delzed(0) / max(maxval(abs(stream_correction)), zero)
      !    cfl_dt_linear = min(cfl_dt_linear, cfl_dt_stream)
      ! end if

      if (.not. mirror_implicit) then
         ! NB: mirror has code_dt built-in, which accounts for code_dt factor here
         cfl_dt_mirror = abs(code_dt) * dvpa / max(maxval(abs(mirror)), zero)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_mirror)
      end if

      ! While other quantities should go here, parallel streaming with electrons is what will limit us
      if (radial_variation) then
         cfl_dt_stream = abs(code_dt) * delzed(0) / max(maxval(abs(stream_rad_var1)), zero)
         cfl_dt_stream = cfl_dt_stream / abs(rho(nx) + zero)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_stream)

         cfl_dt_stream = abs(code_dt) * delzed(0) / max(maxval(abs(stream_rad_var2)), zero)
         cfl_dt_stream = cfl_dt_stream / abs(rho(nx) + zero)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_stream)
      end if

      if (include_collisions .and. .not. collisions_implicit) then
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_vpadiff)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_mudiff)
      end if

      if (.not. drifts_implicit .and. include_ydrift) then
      
         ! Get the local max value of wdrifty on each processor
         wdrifty_max = maxval(abs(wdrifty_g))
         
         ! Compare these max values across processors to get global max
         if (nproc > 1) then
            call max_allreduce(wdrifty_max)
         end if
         
         ! wdrifty_g has code_dt built-in, which accounts for code_dt factor here
         cfl_dt_wdrifty = abs(code_dt) / max(maxval(abs(aky)) * wdrifty_max, zero)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_wdrifty)
         
      end if
   
      
     if (.not. drifts_implicit .and. include_drive) then

         wstar_max = maxval(abs(wstar))
         if (nproc > 1) then
             call max_allreduce(wstar_max)
         end if
         
         cfl_dt_wstar = abs(code_dt) / max(maxval(abs(aky)) * wstar_max, zero)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_wstar)

      end if

      ! Check that the introduction of the neoclassical chi coeffecient doesn't break the CFL condition. 
      ! if (neoclassical_is_enabled()) then
          ! neo_chi_coeff_max = maxval(abs(neo_chi_coeff))
          ! if (nproc > 1) then
              ! call max_allreduce(neo_chi_coeff_max)
          ! end if

          ! cfl_dt_neo_chi_coeff = abs(code_dt) / max(neo_chi_coeff_max, zero)
          ! cfl_dt_linear = min(cfl_dt_linear, cfl_dt_neo_chi_coeff)
      ! end if

      ! Check that the introduction of the neoclassical apar coeffecient doesn't break the CFL condition.
      if (neoclassical_is_enabled()) then
          neo_apar_coeff_max = maxval(abs(neo_apar_coeff))
          if (nproc > 1) then
              call max_allreduce(neo_apar_coeff_max)
          end if

          cfl_dt_neo_apar_coeff = abs(code_dt) / max(neo_apar_coeff_max, zero)
          cfl_dt_linear = min(cfl_dt_linear, cfl_dt_neo_apar_coeff)
      end if

      ! Check that the introduction of the neoclassical dchi/dz coeffecient doesn't break the CFL condition.
      if (neoclassical_is_enabled()) then
          neo_dchidz_coeff_max = maxval(abs(neo_dchidz_coeff))
          if (nproc > 1) then
              call max_allreduce(neo_dchidz_coeff_max)
          end if

          cfl_dt_neo_dchidz_coeff = abs(code_dt) / max(neo_dchidz_coeff_max, zero)
          cfl_dt_linear = min(cfl_dt_linear, cfl_dt_neo_dchidz_coeff)
      end if

      ! Check that the introduction of the neoclassical wstar1 drive doesn't break the CFL condition.
      if (neoclassical_is_enabled()) then
          wstar1_max = maxval(abs(wstar1))
          if (nproc > 1) then
              call max_allreduce(wstar1_max)
          end if

          cfl_dt_wstar1 = abs(code_dt) / max(maxval(abs(aky)) * wstar1_max, zero)
          cfl_dt_linear = min(cfl_dt_linear, cfl_dt_wstar1)
      end if

      ! Check that the introduction of the neoclassical wpol drive doesn't break the CFL condition.
      if (neoclassical_is_enabled()) then
          ! Only calculate the CFL constaint if there are non-zero akx present. 
          if (maxval(abs(akx)) > epsilon(0.0)) then
              wpol_max = maxval(abs(wpol))
              if (nproc > 1) then
                  call max_allreduce(wpol_max)
              end if

              cfl_dt_wpol = abs(code_dt) / max(maxval(abs(akx)) * wpol_max, zero) 
              cfl_dt_linear = min(cfl_dt_linear, cfl_dt_wpol)
          end if
      end if

      ! Check that the introduction of the neoclassical neomagy doesn't break the CFL condition.
      ! if (neoclassical_is_enabled()) then
          ! neomagy_max = maxval(abs(neomagy))
          ! if (nproc > 1) then
              ! call max_allreduce(neomagy_max)
          ! end if

          ! cfl_dt_neomagy = abs(code_dt) / max(maxval(abs(aky)) * neomagy_max, zero)
          ! cfl_dt_linear = min(cfl_dt_linear, cfl_dt_neomagy)
      ! end if

      ! Check that the introduction of the neoclassical neomagx drive doesn't break the CFL condition.
      ! if (neoclassical_is_enabled()) then
          ! Only calculate the CFL constaint if there are non-zero akx present. 
          ! if (maxval(abs(akx)) > epsilon(0.0)) then
              ! neomagx_max = maxval(abs(neomagx))
              ! if (nproc > 1) then
                  ! call max_allreduce(neomagx_max)
              ! end if

              ! cfl_dt_neomagx = abs(code_dt) / max(maxval(abs(akx)) * neomagx_max, zero)
              ! cfl_dt_linear = min(cfl_dt_linear, cfl_dt_neomagx)
          ! end if
      ! end if

      ! Check that the introduction of the neoclassical neocurvy doesn't break the CFL condition.
      if (neoclassical_is_enabled()) then
          neocurvy_max = maxval(abs(neocurvy))
          if (nproc > 1) then
              call max_allreduce(neocurvy_max)
          end if

          cfl_dt_neocurvy = abs(code_dt) / max(maxval(abs(aky)) * neocurvy_max, zero)
          cfl_dt_linear = min(cfl_dt_linear, cfl_dt_neocurvy)
      end if

      ! Check that the introduction of the neoclassical neocurvx drive doesn't break the CFL condition.
      if (neoclassical_is_enabled()) then
          ! Only calculate the CFL constaint if there are non-zero akx present. 
          if (maxval(abs(akx)) > epsilon(0.0)) then
              neocurvx_max = maxval(abs(neocurvx))
              if (nproc > 1) then
                  call max_allreduce(neocurvx_max)
              end if

              cfl_dt_neocurvx = abs(code_dt) / max(maxval(abs(akx)) * neocurvx_max, zero)
              cfl_dt_linear = min(cfl_dt_linear, cfl_dt_neocurvx)
          end if
      end if

      ! Get the minimum value of <cfl_dt_linear> on all processors
      if (runtype_option_switch == runtype_multibox) call scope(allprocs)
      call min_allreduce(cfl_dt_linear)
      if (runtype_option_switch == runtype_multibox) call scope(subprocs)

      ! Print the CFL condition to the command prompt
      if (proc0 .and. print_extra_info_to_terminal) then
         write (*, '(A)') "############################################################"
         write (*, '(A)') "                        CFL CONDITION"
         write (*, '(A)') "############################################################"
         write (*, '(A16)') 'LINEAR CFL_DT: '
         if (.not. drifts_implicit .and. include_xdrift) write (*, '(A12,ES12.4)') '   wdriftx: ', cfl_dt_wdriftx
         if (.not. drifts_implicit .and. include_ydrift) write (*, '(A12,ES12.4)') '   wdrifty: ', cfl_dt_wdrifty
         if (.not. drifts_implicit .and. include_drive) write (*, '(A12,ES12.4)') '   wstar: ', cfl_dt_wstar
         if (.not. stream_implicit) write (*, '(A12,ES12.4)') '   stream: ', cfl_dt_stream
         if (.not. mirror_implicit) write (*, '(A12,ES12.4)') '   mirror: ', cfl_dt_mirror
         ! if (neoclassical_is_enabled()) write (*, '(A12,ES12.4)') '   neo_chi: ', cfl_dt_neo_chi_coeff
         if (neoclassical_is_enabled()) write (*, '(A12,ES12.4)') '  neo_apar: ', cfl_dt_neo_apar_coeff
         if (neoclassical_is_enabled()) write (*, '(A12,ES12.4)') 'neo_dchidz: ', cfl_dt_neo_dchidz_coeff
         if (neoclassical_is_enabled()) write (*, '(A12,ES12.4)') 'wstar1: ', cfl_dt_wstar1
         if (neoclassical_is_enabled() .and. maxval(abs(akx)) > epsilon(0.0)) write (*, '(A12,ES12.4)') 'wpol: ', cfl_dt_wpol
         ! if (neoclassical_is_enabled()) write (*, '(A12,ES12.4)') 'neomagy: ', cfl_dt_neomagy
         ! if (neoclassical_is_enabled() .and. maxval(abs(akx)) > epsilon(0.0)) write (*, '(A12,ES12.4)') 'neomagx: ', cfl_dt_neomagx
         if (neoclassical_is_enabled()) write (*, '(A12,ES12.4)') 'neocurvy: ', cfl_dt_neocurvy
         if (neoclassical_is_enabled() .and. maxval(abs(akx)) > epsilon(0.0)) write (*, '(A12,ES12.4)') 'neocurvx: ', cfl_dt_neocurvx 
         write (*, '(A12,ES12.4)') '   total: ', cfl_dt_linear
         write (*, *)
      end if

      ! Reduce the time step if it's much larger than the CFL condition
      if (abs(code_dt) > cfl_dt_linear * cfl_cushion_upper) then
         if (proc0) then
            write (*, *) 'CHANGING TIME STEP:'
            write (*, '(A22, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
            write (*, '(A22, ES10.2E2)') "   cfl_dt_linear:"//REPEAT(' ', 50), cfl_dt_linear
            write (*, '(A22, ES10.2E2)') "   cfl_cushion_upper:"//REPEAT(' ', 50), cfl_cushion_upper
            write (*, '(A22, ES10.2E2)') "   cfl_cushion_middle:"//REPEAT(' ', 50), cfl_cushion_middle
            write (*, '(A22, ES10.2E2)') "   cfl_cushion_lower:"//REPEAT(' ', 50), cfl_cushion_lower
            write (*, '(A70)') '     ==> User-specified delt is larger than cfl_dt*cfl_cushion_upper.'//REPEAT(' ', 50)
            write (*, '(A55,ES12.4)') '     ==> Changing code_dt to cfl_dt*cfl_cushion_upper ='//REPEAT(' ', 50), cfl_dt_linear * cfl_cushion_upper
            write (*, *)
         end if
         code_dt = sign(1.0, code_dt) * cfl_dt_linear * cfl_cushion_upper
         call reset_dt
      else if (proc0 .and. print_extra_info_to_terminal) then
         call write_dt
         write (*, *)
      end if

    end subroutine init_cfl

   !****************************************************************************
   !                            Reset the time step                             
   !****************************************************************************
   subroutine reset_dt

      use parameters_numerical, only: stream_implicit, driftkinetic_implicit
      use parameters_physics, only: radial_variation
      
      use arrays, only: initialised_radial_variation, initialised_implicit_drifts
      use arrays, only: initialised_wdrift, initialised_wstar
      
      use gk_mirror, only: init_mirror
      use gk_flow_shear, only: init_flow_shear
      use gk_drive, only: init_wstar
      use gk_magnetic_drift, only: init_wdrift
      use gk_radial_variation, only: init_radial_variation
      use gk_sources, only: init_source_timeaverage
      use gk_sources, only: init_quasineutrality_source
      use gk_parallel_streaming, only: init_parallel_streaming
      use response_matrix, only: init_response_matrix
      use dissipation_and_collisions, only: init_collisions
      use dissipation_and_collisions, only: include_collisions
      
      ! Flags -> Logicals for initalisation (true or false)
      use gk_parallel_streaming, only: initialised_parallel_streaming
      use dissipation_and_collisions, only: initialised_collisions
      use response_matrix, only: initialised_response_matrix
      use gk_flow_shear, only: initialised_flow_shear
      use gk_sources, only: initialised_qn_source
      use gk_mirror, only: initialised_mirror

      implicit none

      !----------------------------------------------------------------------
         
      ! Need to recompute mirror and streaming terms to account for updated code_dt
      initialised_wdrift = .false.
      initialised_wstar = .false.
      initialised_radial_variation = .false.
      initialised_implicit_drifts = .false.
      initialised_flow_shear = .false.
      initialised_mirror = .false.
      initialised_parallel_streaming = .false.
      initialised_qn_source = .false.

      ! Re-initialise all routines that depend on the time step
      if (debug) write (6, *) 'time_advance::reset_dt::init_wstar'
      call init_wstar 
      if (debug) write (6, *) 'time_advance::reset_dt::init_wdrift'
      call init_wdrift 
      if (debug) write (6, *) 'time_advance::reset_dt::init_mirror'
      call init_mirror
      if (debug) write (6, *) 'time_advance::reset_dt::init_parallel_streaming'
      call init_parallel_streaming
      if (debug) write (6, *) 'time_advance::reset_dt::init_flow_shear'
      call init_flow_shear
      if (debug) write (6, *) 'time_advance::reset_dt::init_source_timeaverage'
      call init_source_timeaverage
      if (debug) write (6, *) 'time_advance::reset_dt::init_quasineutrality_source'
      call init_quasineutrality_source
      if (radial_variation) then
         if (debug) write (6, *) 'time_advance::reset_dt::init_radial_variation'
         call init_radial_variation
      end if
      if (include_collisions) then
         if (debug) write (6, *) 'time_advance::reset_dt::init_collisions'
         initialised_collisions = .false.
         call init_collisions
      end if
      
      ! Do not try to re-init response matrix before it has been initialised the first time
      if ((stream_implicit .or. driftkinetic_implicit) .and. initialised_response_matrix) then
         initialised_response_matrix = .false.
         if (debug) write (6, *) 'time_advance::reset_dt::init_response_matrix'
         call init_response_matrix
      end if

   end subroutine reset_dt

end module calculations_timestep
