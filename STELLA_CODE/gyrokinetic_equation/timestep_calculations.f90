!###############################################################################
!###############################################################################
!###############################################################################
! 
! This module ...
! 
!###############################################################################
module timestep_calculations

   use debug_flags, only: debug => time_advance_debug

   implicit none

   public :: init_cfl
   public :: reset_dt

   private
    
contains
    

   subroutine init_cfl

      use mp, only: proc0, nproc, max_allreduce, min_allreduce
      use mp, only: scope, allprocs, subprocs
      use file_utils, only: runtype_option_switch, runtype_multibox
      use arrays_store_useful, only: wdriftx_g, wdrifty_g
      use stella_time, only: code_dt, write_dt, cfl_dt_linear
      use parameters_numerical, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
      use parameters_numerical, only: stream_implicit, mirror_implicit, drifts_implicit
      use parameters_physics, only: radial_variation, prp_shear_enabled
      use parameters_kxky_grid, only: nx
      use grids_z, only: delzed
      use grids_velocity, only: dvpa
      use grids_kxky, only: akx, aky, rho
      use gk_parallel_streaming, only: stream
      use gk_parallel_streaming, only: stream_rad_var1, stream_rad_var2
      use gk_mirror, only: mirror
      use gk_flow_shear, only: prl_shear, shift_times
      use dissipation, only: include_collisions, collisions_implicit
      use dissipation, only: cfl_dt_vpadiff, cfl_dt_mudiff
      use debug_flags, only: print_extra_info_to_terminal

      implicit none

      real :: cfl_dt_mirror, cfl_dt_stream, cfl_dt_shear
      real :: cfl_dt_wdriftx, cfl_dt_wdrifty
      real :: zero
      real :: wdriftx_max, wdrifty_max

      ! Avoid divide by zero in cfl_dt terms below
      zero = 100.*epsilon(0.)

      ! FLAG -- assuming equal spacing in zed!

      if (cfl_dt_linear < 0) cfl_dt_linear = code_dt / cfl_cushion_upper

      if (.not. drifts_implicit) then
         ! Get the local max value of wdriftx on each processor
         wdriftx_max = maxval(abs(wdriftx_g))
         ! Xompare these max values across processors to get global max
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

      if (.not. stream_implicit) then
         ! NB: stream has code_dt built-in, which accounts for code_dt factor here
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

      if (radial_variation) then
         ! While other quantities should go here, parallel streaming with electrons
         ! is what will limit us
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

      if (.not. drifts_implicit) then
         ! Get the local max value of wdrifty on each processor
         wdrifty_max = maxval(abs(wdrifty_g))
         ! Compare these max values across processors to get global max
         if (nproc > 1) then
            call max_allreduce(wdrifty_max)
         end if
         ! NB: wdrifty_g has code_dt built-in, which accounts for code_dt factor here
         cfl_dt_wdrifty = abs(code_dt) / max(maxval(abs(aky)) * wdrifty_max, zero)
         cfl_dt_linear = min(cfl_dt_linear, cfl_dt_wdrifty)
      end if

      if (runtype_option_switch == runtype_multibox) call scope(allprocs)
      call min_allreduce(cfl_dt_linear)
      if (runtype_option_switch == runtype_multibox) call scope(subprocs)

      if (proc0 .and. print_extra_info_to_terminal) then
         write (*, '(A)') "############################################################"
         write (*, '(A)') "                        CFL CONDITION"
         write (*, '(A)') "############################################################"
         write (*, '(A16)') 'LINEAR CFL_DT: '
         if (.not. drifts_implicit) write (*, '(A12,ES12.4)') '   wdriftx: ', cfl_dt_wdriftx
         if (.not. drifts_implicit) write (*, '(A12,ES12.4)') '   wdrifty: ', cfl_dt_wdrifty
         if (.not. stream_implicit) write (*, '(A12,ES12.4)') '   stream: ', cfl_dt_stream
         if (.not. mirror_implicit) write (*, '(A12,ES12.4)') '   mirror: ', cfl_dt_mirror
         write (*, '(A12,ES12.4)') '   total: ', cfl_dt_linear
         write (*, *)
      end if

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

    subroutine reset_dt

        use dissipation, only: init_collisions, collisions_initialized, include_collisions
        use parameters_numerical, only: stream_implicit, driftkinetic_implicit
        use parameters_physics, only: radial_variation
        use response_matrix, only: response_matrix_initialized
        use response_matrix, only: init_response_matrix
        use gk_parallel_streaming, only: parallel_streaming_initialized
        use gk_parallel_streaming, only: init_parallel_streaming
        use gk_mirror, only: mirror_initialized
        use gk_mirror, only: init_mirror
        use gk_flow_shear, only: flow_shear_initialized
        use gk_flow_shear, only: init_flow_shear
        use gk_drive, only: init_wstar
        use gk_magnetic_drift, only: init_wdrift
        use gk_radial_variation, only: init_radial_variation
        use gk_sources, only: init_source_timeaverage
        use gk_sources, only: init_quasineutrality_source, qn_source_initialized

        use arrays_store_useful, only: wdriftinit, wstarinit, &
                        radialinit, driftimpinit
            

        ! need to recompute mirror and streaming terms
        ! to account for updated code_dt
        wdriftinit = .false.
        wstarinit = .false.
        radialinit = .false.
        driftimpinit = .false.
        flow_shear_initialized = .false.
        mirror_initialized = .false.
        parallel_streaming_initialized = .false.
        qn_source_initialized = .false.

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
            collisions_initialized = .false.
            call init_collisions
        end if
        ! do not try to re-init response matrix
        ! before it has been initialized the first time
        if ((stream_implicit .or. driftkinetic_implicit) .and. response_matrix_initialized) then
            response_matrix_initialized = .false.
            if (debug) write (6, *) 'time_advance::reset_dt::init_response_matrix'
            call init_response_matrix
        end if

   end subroutine reset_dt

end module timestep_calculations
