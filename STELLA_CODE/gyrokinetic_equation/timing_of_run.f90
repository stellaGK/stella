module timing_of_run
    use debug_flags, only: debug => time_advance_debug

    implicit none

    public :: reset_dt
    private
    
contains
    
    subroutine reset_dt

        use parallel_streaming, only: parallel_streaming_initialized
        use parallel_streaming, only: init_parallel_streaming
        use dissipation, only: init_collisions, collisions_initialized, include_collisions
        use parameters_numerical, only: stream_implicit, driftkinetic_implicit
        use response_matrix, only: response_matrix_initialized
        use response_matrix, only: init_response_matrix
        use mirror_terms, only: mirror_initialized
        use mirror_terms, only: init_mirror
        use flow_shear, only: flow_shear_initialized
        use flow_shear, only: init_flow_shear
        use parameters_physics, only: radial_variation
        use sources, only: init_source_timeaverage
        use sources, only: init_quasineutrality_source, qn_source_initialized
        use arrays_drifts, only: init_wdrift, init_wstar
        use radial_variation_time_advance, only: init_radial_variation

        use store_arrays_useful, only: wdriftinit, wstarinit, parnlinit, &
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

end module timing_of_run