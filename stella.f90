program stella

  use redistribute, only: scatter
  use job_manage, only: time_message, checkstop
  use run_parameters, only: nstep, fphi, fapar
  use stella_time, only: update_time, code_time, code_dt
  use dist_redistribute, only: kxkyz2vmu
  use time_advance, only: advance_stella
  use stella_diagnostics, only: diagnose_stella, nsave
  use stella_save, only: stella_save_for_restart
  use dist_fn_arrays, only: gnew, gvmu

  implicit none

  logical :: debug = .false.
  logical :: stop_stella = .false.

  integer :: istep
  integer :: istatus
  real, dimension (2) :: time_init = 0.
  real, dimension (2) :: time_diagnostics = 0.
  real, dimension (2) :: time_total = 0.

  call init_stella

  if (debug) write (*,*) 'stella::diagnose_stella'
  call diagnose_stella (0)

  if (debug) write (*,*) 'stella::advance_stella'
  do istep = 1, nstep
     if (debug) write (*,*) 'istep = ', istep
     call advance_stella (istep)
!     call advance_stella
     if (nsave > 0 .and. mod(istep,nsave)==0) then
        call scatter (kxkyz2vmu, gnew, gvmu)
        call stella_save_for_restart (gvmu, code_time, code_dt, istatus, fphi, fapar)
     end if
     call update_time
     call time_message(.false.,time_diagnostics,' diagnostics')
     call diagnose_stella (istep)
     call time_message(.false.,time_diagnostics,' diagnostics')
     if (mod(istep,10)==0) call checkstop (stop_stella)
     if (stop_stella) exit
  end do

  if (debug) write (*,*) 'stella::finish_stella'
  call finish_stella

contains

  subroutine init_stella

    use mp, only: init_mp, broadcast
    use mp, only: proc0
    use file_utils, only: init_file_utils
    use file_utils, only: run_name
    use job_manage, only: checktime, time_message
    use physics_parameters, only: init_physics_parameters
    use physics_flags, only: init_physics_flags
    use run_parameters, only: init_run_parameters
    use run_parameters, only: avail_cpu_time, nstep
    use run_parameters, only: stream_implicit, driftkinetic_implicit
    use species, only: init_species
    use zgrid, only: init_zgrid
    use stella_geometry, only: init_geometry
    use stella_geometry, only: geo_surf, twist_and_shift_geo_fac
    use stella_layouts, only: init_stella_layouts
    use response_matrix, only: init_response_matrix
    use init_g, only: ginit, init_init_g
    use fields, only: init_fields, advance_fields
    use stella_time, only: init_tstart
    use init_g, only: tstart
    use stella_diagnostics, only: init_stella_diagnostics
    use fields_arrays, only: phi, apar
    use dist_fn_arrays, only: gnew
    use dist_fn, only: init_gxyz, init_dist_fn
    use time_advance, only: init_time_advance
    use extended_zgrid, only: init_extended_zgrid
    use kt_grids, only: init_kt_grids, read_kt_grids_parameters
    use vpamu_grids, only: init_vpamu_grids

    implicit none

    logical :: exit, list, restarted
    character (500), target :: cbuff

    ! initialize mpi message passing
    call init_mp
    debug = debug .and. proc0
    if (debug) write (*,*) 'stella::init_stella::check_time'
    ! initialize timer
    call checktime(avail_cpu_time,exit)
    if (proc0) then
       if (debug) write (*,*) 'stella::init_stella::write_start_message'
       ! write message to screen with useful info
       ! regarding start of simulation
       call write_start_message
       if (debug) write (*,*) 'stella::init_stella::init_file_utils'
       ! initialize file i/o
       call init_file_utils (list)
       call time_message(.false.,time_total,' Total')
       call time_message(.false.,time_init,' Initialization')
       cbuff = trim(run_name)
    end if

    call broadcast (cbuff)
    if (.not. proc0) run_name => cbuff

    if (debug) write(6,*) "stella::init_stella::init_zgrid"
    call init_zgrid
    if (debug) write (6,*) "stella::init_stella::read_kt_grids_parameters"
    call read_kt_grids_parameters
    if (debug) write(6,*) "stella::init_stella::init_geometry"
    call init_geometry
    if (debug) write(6,*) "stella::init_stella::init_physics_parameters"
    call init_physics_parameters
    if (debug) write(6,*) "stella::init_stella::init_physics_flags"
    call init_physics_flags
    if (debug) write (6,*) 'stella::init_stella::init_species'
    call init_species
    if (debug) write(6,*) "stella::init_stella::init_init_g"
    call init_init_g
    if (debug) write(6,*) "stella::init_stella::init_run_parameters"
    call init_run_parameters
    if (debug) write (6,*) 'stella::init_stella::init_stella_layouts'
    call init_stella_layouts
    if (debug) write (6,*) 'stella::init_stella::init_kt_grids'
    call init_kt_grids (geo_surf, twist_and_shift_geo_fac)
    if (debug) write (6,*) 'stella::init_stella::init_vpamu_grids'
    call init_vpamu_grids
    if (debug) write(6,*) "stella::init_stella::init_dist_fn"
    call init_dist_fn
    if (debug) write (6,*) 'stella::init_stella::init_extended_zgrid'
    call init_extended_zgrid
    if (debug) write (6,*) 'stella::init_stella::init_fields'
    call init_fields
    if (debug) write (6,*) 'stella::init_stella::init_time_advance'
    call init_time_advance
    if (debug) write(6,*) "stella::init_stella::ginit"
    call ginit (restarted)
    if (debug) write(6,*) "stella::init_stella::init_gxyz"
    call init_gxyz
    if (debug) write(6,*) "stella::init_stella::init_response_matrix"
    if (stream_implicit .or. driftkinetic_implicit) call init_response_matrix

    if (.not.restarted) then
       if (debug) write (6,*) 'stella::init_stella::get_fields'
       ! get initial field from initial distribution function
       call advance_fields (gnew, phi, apar, dist='gbar')
    end if

    if (debug) write (6,*) 'stella::init_stella::init_stella_diagnostics'
    call init_stella_diagnostics (nstep)
    if (debug) write (6,*) 'stella::init_stella::init_tstart'
    call init_tstart (tstart)

    if (proc0) call time_message(.false.,time_init,' Initialization')

  end subroutine init_stella

  subroutine write_start_message

    use mp, only: proc0, nproc

    implicit none

    if (proc0) then
       if (nproc==1) then
          write (*,*) "Running on ", nproc, " processor"
       else
          write (*,*) "Running on ", nproc, " processors"
       end if
       write (*,*)
    end if

  end subroutine write_start_message

  subroutine finish_stella

    use mp, only: finish_mp
    use mp, only: proc0
    use file_utils, only: finish_file_utils
    use job_manage, only: time_message
    use physics_parameters, only: finish_physics_parameters
    use physics_flags, only: finish_physics_flags
    use run_parameters, only: finish_run_parameters
    use zgrid, only: finish_zgrid
    use species, only: finish_species
    use time_advance, only: time_gke, time_parallel_nl
    use time_advance, only: finish_time_advance
    use parallel_streaming, only: time_parallel_streaming
    use mirror_terms, only: time_mirror
    use dissipation, only: time_collisions
    use dist_fn, only: finish_dist_fn
    use fields, only: finish_fields
    use fields, only: time_field_solve
    use stella_diagnostics, only: finish_stella_diagnostics
    use response_matrix, only: finish_response_matrix
    use stella_geometry, only: finish_geometry
    use extended_zgrid, only: finish_extended_zgrid
    use vpamu_grids, only: finish_vpamu_grids
    use kt_grids, only: finish_kt_grids

    implicit none

    if (debug) write (*,*) 'stella::finish_stella::finish_stella_diagnostics'
    call finish_stella_diagnostics
    if (debug) write (*,*) 'stella::finish_stella::finish_response_matrix'
    call finish_response_matrix
    if (debug) write (*,*) 'stella::finish_stella::finish_fields'
    call finish_fields
    if (debug) write (*,*) 'stella::finish_stella::finish_time_advance'
    call finish_time_advance
    if (debug) write (*,*) 'stella::finish_stella::finish_extended_zgrid'
    call finish_extended_zgrid
    if (debug) write (*,*) 'stella::finish_stella::finish_dist_fn'
    call finish_dist_fn
    if (debug) write (*,*) 'stella::finish_stella::finish_vpamu_grids'
    call finish_vpamu_grids
    if (debug) write (*,*) 'stella::finish_stella::finish_kt_grids'
    call finish_kt_grids
    if (debug) write (*,*) 'stella::finish_stella::finish_run_parameters'
    call finish_run_parameters
    if (debug) write (*,*) 'stella::finish_stella::finish_species'
    call finish_species
    if (debug) write (*,*) 'stella::finish_stella::finish_physics_flags'
    call finish_physics_flags
    if (debug) write (*,*) 'stella::finish_stella::finish_physics_parameters'
    call finish_physics_parameters
    if (debug) write (*,*) 'stella::finish_stella::finish_geometry'
    call finish_geometry
    if (debug) write (*,*) 'stella::finish_stella::finish_zgrid'
    call finish_zgrid
    if (debug) write (*,*) 'stella::finish_stella::finish_file_utils'
    if (proc0) then
       call finish_file_utils
       call time_message(.false.,time_total,' Total')
       write (*,*)
       write (*,fmt=101) 'initialization:', time_init(1)/60., 'min'
       write (*,fmt=101) 'diagnostics:', time_diagnostics(1)/60., 'min'
       write (*,fmt=101) 'fields:', time_field_solve(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_field_solve(1,2)/60., 'min'
       write (*,fmt=101) 'mirror:', time_mirror(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_mirror(1,2)/60., 'min'
       write (*,fmt=101) 'stream:', time_parallel_streaming(1)/60., 'min'
       write (*,fmt=101) 'dgdx:', time_gke(1,5)/60., 'min'
       write (*,fmt=101) 'dgdy:', time_gke(1,4)/60., 'min'
       write (*,fmt=101) 'wstar:', time_gke(1,6)/60., 'min'
       write (*,fmt=101) 'collisions:', time_collisions(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_collisions(1,2)/60., 'min'
       write (*,fmt=101) 'ExB nonlin:', time_gke(1,7)/60., 'min'
       write (*,fmt=101) 'parallel nonlin:', time_parallel_nl(1,1)/60., 'min'
       write (*,fmt=101) '(redistribute):', time_parallel_nl(1,2)/60., 'min'
       write (*,fmt=101) 'total implicit: ', time_gke(1,9)/60., 'min'
       write (*,fmt=101) 'total explicit: ', time_gke(1,8)/60., 'min'
       write (*,fmt=101) 'total:', time_total(1)/60., 'min'
       write (*,*)
    end if
101 format (a17,0pf8.2,a4)

    if (debug) write (*,*) 'stella::finish_stella::finish_mp'
    ! finish (clean up) mpi message passing
    call finish_mp

  end subroutine finish_stella

end program stella
