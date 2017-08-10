program stella

  use job_manage, only: time_message
  use run_parameters, only: nstep
  use stella_time, only: update_time
  use dist_fn, only: advance_stella
  use stella_diagnostics, only: diagnose_stella

  implicit none

  logical :: debug = .false.

  integer :: istep
  real, dimension (2) :: time_init = 0.
  real, dimension (2) :: time_diagnostics = 0.
  real, dimension (2) :: time_total = 0.

  call init_stella

  if (debug) write (*,*) 'stella::diagnose_stella'
  call diagnose_stella (0)

  if (debug) write (*,*) 'stella::advance_stella'
  do istep = 1, nstep
     if (debug) write (*,*) 'istep = ', istep
     call advance_stella
     call update_time
     call time_message(.false.,time_diagnostics,' diagnostics')
     call diagnose_stella (istep)
     call time_message(.false.,time_diagnostics,' diagnostics')
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
    use run_parameters, only: avail_cpu_time, nstep
    use fields, only: init_fields
    use stella_time, only: init_tstart
    use init_g, only: tstart
    use stella_diagnostics, only: init_stella_diagnostics

    implicit none

    logical :: exit, list
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

    if (debug) write (*,*) 'stella::init_stella::init_fields'
    call init_fields
    if (debug) write (*,*) 'stella::init_stella::init_stella_diagnostics'
    call init_stella_diagnostics (nstep)
    if (debug) write (*,*) 'stella::init_stella::init_tstart'
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
    use fields, only: finish_fields
    use dist_fn, only: time_gke
    use stella_diagnostics, only: finish_stella_diagnostics

    implicit none

    if (debug) write (*,*) 'stella::finish_stella::finish_stella_diagnostics'
    call finish_stella_diagnostics
    if (debug) write (*,*) 'stella::finish_stella::finish_fields'
    call finish_fields
    if (debug) write (*,*) 'stella::finish_stella::finish_file_utils'
    if (proc0) then
       call finish_file_utils
       call time_message(.false.,time_total,' Total')
       write (*,*)
       write (*,fmt=101) 'initialization:', time_init(1)/60., 'min'
       write (*,fmt=101) 'diagnostics:', time_diagnostics(1)/60., 'min'
       write (*,fmt=101) 'fields:', time_gke(1,2)/60., 'min'
       write (*,fmt=101) 'mirror:', time_gke(1,2)/60., 'min'
       write (*,fmt=101) 'stream:', time_gke(1,3)/60., 'min'
       write (*,fmt=101) 'dgdx:', time_gke(1,5)/60., 'min'
       write (*,fmt=101) 'dgdy:', time_gke(1,4)/60., 'min'
       write (*,fmt=101) 'wstar:', time_gke(1,6)/60., 'min'
       write (*,fmt=101) 'total:', time_total(1)/60., 'min'
       write (*,*)
    end if
101 format (a16,0pf8.2,a4)

    if (debug) write (*,*) 'stella::finish_stella::finish_mp'
    ! finish (clean up) mpi message passing
    call finish_mp

  end subroutine finish_stella

end program stella
