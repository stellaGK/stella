module time_advance

  public :: init_time_advance, finish_time_advance
  public :: advance_stella
  public :: time_gke, time_parallel_nl
  public :: checksum

  private
  
  interface get_dgdy
     module procedure get_dgdy_2d
     module procedure get_dgdy_3d
     module procedure get_dgdy_4d
  end interface

  interface get_dgdx
     module procedure get_dgdx_2d
     module procedure get_dgdx_3d
     module procedure get_dgdx_4d
  end interface

  interface get_dchidy
     module procedure get_dchidy_4d
     module procedure get_dchidy_2d
  end interface

  interface checksum
     module procedure checksum_field
     module procedure checksum_dist
  end interface

  logical :: time_advance_initialized = .false.

  logical :: wdriftinit = .false.
  logical :: wstarinit = .false.
  logical :: parnlinit = .false.
  logical :: readinit = .false.

  ! if .true., dist fn is represented on alpha grid
  ! if .false., dist fn is given on k-alpha grid
  ! default is .false.; will only ever be set to
  ! .true. during full_flux_surface simulations
!  logical :: alpha_space = .false.

  integer :: explicit_option_switch
  integer, parameter :: explicit_option_rk3 = 1, &
       explicit_option_rk2 = 2, &
       explicit_option_rk4 = 3

!  logical :: wdrifty_explicit, wdrifty_implicit
!  logical :: wstar_explicit, wstar_implicit
  real :: xdriftknob, ydriftknob, wstarknob
  logical :: flip_flop

!  complex, dimension (:,:,:), allocatable :: gamtot_wstar, apar_denom_wstar
!  complex, dimension (:,:), allocatable :: gamtot3_wstar

  ! factor multiplying parallel nonlinearity
  real, dimension (:,:), allocatable :: par_nl_fac
  ! factor multiplying higher order linear term in parallel acceleration
  real, dimension (:,:), allocatable :: par_nl_curv
  real, dimension (:,:), allocatable :: par_nl_driftx, par_nl_drifty

  ! needed for timing various pieces of gke solve
  real, dimension (2,9) :: time_gke = 0.
  real, dimension (2,2) :: time_parallel_nl = 0.

  logical :: debug = .false.

contains

  subroutine init_time_advance

    use mp, only: proc0
    use stella_transforms, only: init_transforms
    use physics_flags, only: full_flux_surface
    use physics_flags, only: include_parallel_nonlinearity
    use physics_flags, only: nonlinear
    use neoclassical_terms, only: init_neoclassical_terms
    use dissipation, only: init_dissipation
    use parallel_streaming, only: init_parallel_streaming
    use dist_redistribute, only: init_redistribute
    use mirror_terms, only: init_mirror

    implicit none

    if (time_advance_initialized) return
    time_advance_initialized = .true.

    debug = debug .and. proc0
    
    if (debug) write (6,*) 'time_advance::init_time_advance::read_parameters'
    call read_parameters
    if (debug) write (6,*) 'time_advance::init_time_advance::allocate_arrays'
    call allocate_arrays
    if (debug) write (6,*) 'time_advance::init_time_advance::init_neoclassical_terms'
    call init_neoclassical_terms
    if (debug) write (6,*) 'time_advance::init_time_advance::init_mirror'
    call init_mirror
    if (debug) write (6,*) 'time_advance::init_time_advance::init_parstream'
    call init_parallel_streaming
    if (debug) write (6,*) 'time_advance::init_time_advance::init_wdrift'
    call init_wdrift
    if (debug) write (6,*) 'time_advance::init_time_advance::init_wstar'
    call init_wstar
    if (debug) write (6,*) 'time_advance::init_time_advance::init_parallel_nonlinearity'
    if (include_parallel_nonlinearity) call init_parallel_nonlinearity
    if (debug) write (6,*) 'time_advance::init_time_advance::init_dissipation'
    call init_dissipation
    if (debug) write (6,*) 'time_advance::init_time_advance::init_redistribute'
    call init_redistribute
    if (debug) write (6,*) 'time_advance::init_time_advance::init_cfl'
    call init_cfl
!    if (nonlinear .or. full_flux_surface .or. include_parallel_nonlinearity) then
!       if (debug) write (6,*) 'time_advance::init_time_advance::init_transforms'
!       call init_transforms
!    end if

  end subroutine init_time_advance

  subroutine read_parameters

    use file_utils, only: error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast
    use run_parameters, only: fully_explicit

    implicit none

    logical :: taexist

    type (text_option), dimension (4), parameter :: explicitopts = &
         (/ text_option('default', explicit_option_rk3), &
         text_option('rk3', explicit_option_rk3), &
         text_option('rk2', explicit_option_rk2), &
         text_option('rk4', explicit_option_rk4) /)
    character(10) :: explicit_option

    namelist /time_advance_knobs/ xdriftknob, ydriftknob, wstarknob, explicit_option, flip_flop

    integer :: ierr, in_file

    if (readinit) return
    readinit = .true.

    if (proc0) then
       explicit_option = 'default'
       xdriftknob = 1.0
       ydriftknob = 1.0
       wstarknob = 1.0
       flip_flop = .false.
!       wdrifty_explicit = .true.
!       wstar_explicit = .true.

       in_file = input_unit_exist("time_advance_knobs", taexist)
       if (taexist) read (unit=in_file, nml=time_advance_knobs)

       ierr = error_unit()
       call get_option_value &
            (explicit_option, explicitopts, explicit_option_switch, &
            ierr, "explicit_option in time_advance_knobs")
    end if

    call broadcast (explicit_option_switch)
    call broadcast (xdriftknob)
    call broadcast (ydriftknob)
    call broadcast (wstarknob)
    call broadcast (flip_flop)
!    call broadcast (wdrifty_explicit)
!    call broadcast (wstar_explicit)

    if (fully_explicit) flip_flop = .false.

  end subroutine read_parameters 

  subroutine init_wdrift

    use dist_fn_arrays, only: wdriftx_g, wdrifty_g
    use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use species, only: spec
    use zgrid, only: nzgrid
    use kt_grids, only: nalpha
    use stella_geometry, only: cvdrift, gbdrift
    use stella_geometry, only: cvdrift0, gbdrift0
    use stella_geometry, only: gds23, gds24
    use stella_geometry, only: geo_surf
    use stella_geometry, only: dxdpsi, drhodpsi, dydalpha
    use vpamu_grids, only: vpa, vperp2, ztmax, maxwell_mu
    use neoclassical_terms, only: include_neoclassical_terms
    use neoclassical_terms, only: dphineo_dzed, dphineo_drho, dphineo_dalpha
    use neoclassical_terms, only: dfneo_dvpa, dfneo_dzed, dfneo_dalpha

    implicit none

    integer :: ivmu, iv, imu, is
    real :: fac
    real, dimension (:,:), allocatable :: wcvdrifty, wgbdrifty
    real, dimension (:,:), allocatable :: wcvdriftx, wgbdriftx

    if (wdriftinit) return
    wdriftinit = .true.

    if (.not.allocated(wdriftx_phi)) &
         allocate (wdriftx_phi(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    if (.not.allocated(wdrifty_phi)) &
         allocate (wdrifty_phi(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    if (.not.allocated(wdriftx_g)) &
         allocate (wdriftx_g(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    if (.not.allocated(wdrifty_g)) &
         allocate (wdrifty_g(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    allocate (wcvdrifty(nalpha,-nzgrid:nzgrid))
    allocate (wgbdrifty(nalpha,-nzgrid:nzgrid))
    allocate (wcvdriftx(nalpha,-nzgrid:nzgrid))
    allocate (wgbdriftx(nalpha,-nzgrid:nzgrid))

    ! FLAG -- need to deal with shat=0 case.  ideally move away from q as x-coordinate
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)

       fac = -ydriftknob*0.5*code_dt*spec(is)%tz
       ! this is the curvature drift piece of wdrifty with missing factor of vpa
       ! vpa factor is missing to avoid singularity when including 
       ! non-Maxwellian corrections to equilibrium
       wcvdrifty = fac*cvdrift*vpa(iv)
       ! this is the grad-B drift piece of wdrifty
       wgbdrifty = fac*gbdrift*0.5*vperp2(:,:,imu)
       wdrifty_g(:,:,ivmu) = wcvdrifty*vpa(iv) + wgbdrifty
       ! if including neoclassical correction to equilibrium Maxwellian,
       ! then add in v_E^{nc} . grad y dg/dy coefficient here
       if (include_neoclassical_terms) then
          wdrifty_g(:,:,ivmu) = wdrifty_g(:,:,ivmu)+code_dt*0.5*(gds23*dphineo_dzed &
               + drhodpsi*dydalpha*dphineo_drho)
       end if

       wdrifty_phi(:,:,ivmu) = ztmax(iv,is)*maxwell_mu(:,:,imu) &
            * (wgbdrifty + wcvdrifty*vpa(iv))
       ! if including neoclassical corrections to equilibrium,
       ! add in -(Ze/m) * v_curv/vpa . grad y d<phi>/dy * dF^{nc}/dvpa term
       ! and v_E . grad z dF^{nc}/dz (here get the dphi/dy part of v_E)
       if (include_neoclassical_terms) then
          wdrifty_phi(:,:,ivmu) = wdrifty_phi(:,:,ivmu) &
               - 0.5*spec(is)%zt*dfneo_dvpa(:,:,ivmu)*wcvdrifty &
               - code_dt*0.5*dfneo_dzed(:,:,ivmu)*gds23
       end if

       fac = -xdriftknob*0.5*code_dt*spec(is)%tz/geo_surf%shat
       ! this is the curvature drift piece of wdriftx with missing factor of vpa
       ! vpa factor is missing to avoid singularity when including 
       ! non-Maxwellian corrections to equilibrium
       wcvdriftx = fac*cvdrift0*vpa(iv)
       ! this is the grad-B drift piece of wdriftx
       wgbdriftx = fac*gbdrift0*0.5*vperp2(:,:,imu)
       wdriftx_g(:,:,ivmu) = wcvdriftx*vpa(iv) + wgbdriftx
       ! if including neoclassical correction to equilibrium Maxwellian,
       ! then add in v_E^{nc} . grad x dg/dx coefficient here
       if (include_neoclassical_terms) then
          wdriftx_g(:,:,ivmu) = wdriftx_g(:,:,ivmu)+code_dt*0.5*(gds24*dphineo_dzed &
               - dxdpsi*dphineo_dalpha)
       end if
       wdriftx_phi(:,:,ivmu) = ztmax(iv,is)*maxwell_mu(:,:,imu) &
            * (wgbdriftx + wcvdriftx*vpa(iv))
       ! if including neoclassical corrections to equilibrium,
       ! add in (Ze/m) * v_curv/vpa . grad x d<phi>/dx * dF^{nc}/dvpa term
       ! and v_E . grad z dF^{nc}/dz (here get the dphi/dx part of v_E)
       ! and v_E . grad alpha dF^{nc}/dalpha (dphi/dx part of v_E)
       if (include_neoclassical_terms) then
          wdriftx_phi(:,:,ivmu) = wdriftx_phi(:,:,ivmu) &
               - 0.5*spec(is)%zt*dfneo_dvpa(:,:,ivmu)*wcvdriftx &
               + code_dt*0.5*(dfneo_dalpha(:,:,ivmu)*dxdpsi-dfneo_dzed(:,:,ivmu)*gds24)
       end if
    end do

    deallocate (wcvdriftx, wgbdriftx, wcvdrifty, wgbdrifty)

  end subroutine init_wdrift

  subroutine init_wstar

    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use species, only: spec
    use zgrid, only: nzgrid
    use kt_grids, only: nalpha
    use stella_geometry, only: dydalpha, drhodpsi
    use vpamu_grids, only: maxwell_vpa, maxwell_mu
    use vpamu_grids, only: vperp2, vpa
    use dist_fn_arrays, only: wstar
    use neoclassical_terms, only: include_neoclassical_terms
    use neoclassical_terms, only: dfneo_drho

    implicit none

    integer :: is, imu, iv, ivmu
    real, dimension (:,:), allocatable :: energy

    if (wstarinit) return
    wstarinit = .true.

    if (.not.allocated(wstar)) &
         allocate (wstar(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc)) ; wstar = 0.0

    allocate (energy(nalpha,-nzgrid:nzgrid))

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       iv = iv_idx(vmu_lo,ivmu)
       energy = vpa(iv)**2 + vperp2(:,:,imu)
       if (include_neoclassical_terms) then
          wstar(:,:,ivmu) = dydalpha*drhodpsi*wstarknob*0.5*code_dt &
               * (maxwell_vpa(iv)*maxwell_mu(:,:,imu) &
               * (spec(is)%fprim+spec(is)%tprim*(energy-1.5)) &
                - dfneo_drho(:,:,ivmu))
       else
          wstar(:,:,ivmu) = dydalpha*drhodpsi*wstarknob*0.5*code_dt &
               * maxwell_vpa(iv)*maxwell_mu(:,:,imu) &
               * (spec(is)%fprim+spec(is)%tprim*(energy-1.5))
       end if
    end do

    deallocate (energy)

  end subroutine init_wstar

  subroutine init_parallel_nonlinearity

    use physics_parameters, only: rhostar
    use species, only: spec, nspec
    use zgrid, only: nztot, nzgrid
    use kt_grids, only: akx, aky
    use kt_grids, only: nakx, naky
    use stella_geometry, only: geo_surf, drhodpsi
    use stella_geometry, only: gradpar, dbdzed, bmag
    use stella_geometry, only: cvdrift, cvdrift0

    implicit none

    if (.not. allocated(par_nl_fac)) allocate (par_nl_fac(-nzgrid:nzgrid,nspec))
    ! this is the factor multiplying -dphi/dz * dg/dvpa in the parallel nonlinearity
    par_nl_fac = 0.5*rhostar*spread(spec%stm*spec%zt,1,nztot)*spread(gradpar,2,nspec)

    if (.not. allocated(par_nl_curv)) allocate (par_nl_curv(-nzgrid:nzgrid,nspec))
    ! ydriftknob is here because this term comes from bhat x curvature . grad B
    par_nl_curv = -ydriftknob*rhostar*geo_surf%rgeo*geo_surf%betaprim*drhodpsi &
         *spread(dbdzed(1,:)*gradpar/bmag(1,:),2,nspec)/spread(spec%zt,1,nztot)

    if (.not. allocated(par_nl_drifty)) allocate (par_nl_drifty(naky,-nzgrid:nzgrid))
    par_nl_drifty = 0.25*rhostar*spread(aky,2,nztot)*spread(cvdrift(1,:),1,naky)
    if (.not. allocated(par_nl_drifty)) allocate (par_nl_driftx(nakx,-nzgrid:nzgrid))
    par_nl_driftx = 0.25*rhostar*spread(akx/geo_surf%shat,2,nztot)*spread(cvdrift0(1,:),1,nakx)

    parnlinit = .true.

  end subroutine init_parallel_nonlinearity

  subroutine allocate_arrays

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx
    use dist_fn_arrays, only: g1, g2, g3

    implicit none

    if (.not.allocated(g1)) &
         allocate (g1(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    g1 = 0.
    if (.not.allocated(g2) .and. explicit_option_switch/=explicit_option_rk2) then
       allocate (g2(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       g2 = 0.
    else
       allocate (g2(1,1,1,1,1))
    end if
    if (.not.allocated(g3) .and. explicit_option_switch==explicit_option_rk4) then
       allocate (g3(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       g3 = 0.
    else
       allocate (g3(1,1,1,1,1))
    end if

  end subroutine allocate_arrays

  subroutine init_cfl
    
    use mp, only: proc0, nproc, max_allreduce
    use dist_fn_arrays, only: wdriftx_g, wdrifty_g
    use stella_time, only: cfl_dt, code_dt, write_dt
    use run_parameters, only: cfl_cushion
    use zgrid, only: delzed
    use vpamu_grids, only: dvpa
    use kt_grids, only: akx, aky
    use run_parameters, only: stream_implicit, mirror_implicit
    use parallel_streaming, only: stream
    use mirror_terms, only: mirror

    implicit none
    
    real :: cfl_dt_mirror, cfl_dt_stream
    real :: cfl_dt_wdriftx, cfl_dt_wdrifty
    real :: zero
    real :: wdriftx_max, wdrifty_max

    ! avoid divide by zero in cfl_dt terms below
    zero = 100.*epsilon(0.)

    ! FLAG -- assuming equal spacing in zed!

    ! get the local max value of wdriftx on each processor
    wdriftx_max = maxval(abs(wdriftx_g))
    ! compare these max values across processors to get global max
    if (nproc > 1) call max_allreduce (wdriftx_max)
    ! NB: wdriftx_g has code_dt built-in, which accounts for code_dt factor here
    cfl_dt_wdriftx = code_dt/max(maxval(abs(akx))*wdriftx_max,zero)
    cfl_dt = cfl_dt_wdriftx

    if (.not.stream_implicit) then
       ! NB: stream has code_dt built-in, which accounts for code_dt factor here
       cfl_dt_stream = code_dt*delzed(0)/max(maxval(abs(stream)),zero)
       cfl_dt = min(cfl_dt,cfl_dt_stream)
    end if

    if (.not.mirror_implicit) then
       ! NB: mirror has code_dt built-in, which accounts for code_dt factor here
       cfl_dt_mirror = code_dt*dvpa/max(maxval(abs(mirror)),zero)
       cfl_dt = min(cfl_dt,cfl_dt_mirror)
    end if

    ! get the local max value of wdrifty on each processor
    wdrifty_max = maxval(abs(wdrifty_g))
    ! compare these max values across processors to get global max
    if (nproc > 1) call max_allreduce (wdrifty_max)
    ! NB: wdrifty_g has code_dt built-in, which accounts for code_dt factor here
    cfl_dt_wdrifty = code_dt/max(maxval(abs(aky))*wdrifty_max,zero)
    cfl_dt = min(cfl_dt,cfl_dt_wdrifty)
    
    if (proc0) then
       write (*,'(a16)') 'LINEAR CFL_DT: '
       write (*,'(a9,e12.4)') 'wdriftx: ', cfl_dt_wdriftx
       write (*,'(a9,e12.4)') 'wdrifty: ', cfl_dt_wdrifty
!       if (wdrifty_explicit) write (*,'(a9,e12.4)') 'wdrifty: ', cfl_dt_wdrifty
       if (.not.stream_implicit) write (*,'(a9,e12.4)') 'stream: ', cfl_dt_stream
       if (.not.mirror_implicit) write (*,'(a9,e12.4)') 'mirror: ', cfl_dt_mirror
       write (*,*)
    end if

    if (code_dt > cfl_dt*cfl_cushion) then
       if (proc0) then
          write (*,'(a21,e12.4,a35)') 'user-specified delt =', code_dt, 'is larger than cfl_dt*cfl_cushion.'
          write (*,'(a41,e12.4)') 'changing code_dt to cfl_dt*cfl_cushion =', cfl_dt*cfl_cushion
       end if
       code_dt = cfl_dt*cfl_cushion
       call reset_dt
    else if (proc0) then
       call write_dt
       write (*,*)
    end if
    
  end subroutine init_cfl

  subroutine reset_dt

    use parallel_streaming, only: parallel_streaming_initialized
    use parallel_streaming, only: init_parallel_streaming
    use run_parameters, only: stream_implicit, driftkinetic_implicit
    use response_matrix, only: response_matrix_initialized
    use response_matrix, only: init_response_matrix
    use mirror_terms, only: mirror_initialized
    use mirror_terms, only: init_mirror

    implicit none

    ! need to recompute mirror and streaming terms
    ! to account for updated code_dt
    wdriftinit = .false.
    wstarinit = .false.
    mirror_initialized = .false.
    parallel_streaming_initialized = .false.
    call init_wstar
    call init_wdrift
    call init_mirror
    call init_parallel_streaming
    ! do not try to re-init response matrix
    ! before it has been initialized the first time
    if ((stream_implicit.or.driftkinetic_implicit) .and. response_matrix_initialized) then
       response_matrix_initialized = .false.
       call init_response_matrix
    end if
!     if (wstar_implicit) then
!        get_fields_wstar_initialized = .false.
!        ! call to init_get_fields only needed for initial step
!        call init_get_fields
!        call init_get_fields_wstar
!     end if

  end subroutine reset_dt

  subroutine advance_stella (istep)

    use dist_fn_arrays, only: gold, gnew
    use fields_arrays, only: phi, apar
    use fields_arrays, only: phi_old
    use run_parameters, only: fully_explicit

    implicit none

    integer, intent (in) :: istep

    ! save value of phi
    ! for use in diagnostics (to obtain frequency)
    phi_old = phi

    ! reverse the order of operations every time step
    ! as part of alternating direction operator splitting
    ! this is needed to ensure 2nd order accuracy in time
    if (mod(istep,2)==1 .or. .not.flip_flop) then
       ! advance the explicit parts of the GKE
!    call advance_explicit (phi, apar, gnew)
       call advance_explicit (gnew)

       ! enforce periodicity for zonal mode
!    if (zonal_mode(1)) gnew(1,:,-nzgrid,:) = gnew(1,:,nzgrid,:)

       ! use operator splitting to separately evolve
       ! all terms treated implicitly
       if (.not.fully_explicit) call advance_implicit (istep, phi, apar, gnew)
    else
       if (.not.fully_explicit) call advance_implicit (istep, phi, apar, gnew)
       call advance_explicit (gnew)
    end if

!    else
!       ! use operator splitting to separately evolve
!       ! all terms treated implicitly
!       if (.not.fully_explicit) call advance_implicit (istep, phi, apar, gnew)
!
!       ! advance the explicit parts of the GKE
!       if (.not.fully_implicit) call advance_explicit (phi, apar, gnew)
!    end if

    ! next line is likely unnecessary
    gold = gnew

  end subroutine advance_stella

!  subroutine advance_explicit (phi, apar, g)
  subroutine advance_explicit (g)

    use mp, only: proc0
    use job_manage, only: time_message
    use zgrid, only: nzgrid
    use kt_grids, only: zonal_mode
    use stella_layouts, only: vmu_lo, iv_idx
    use parallel_streaming, only: stream_sign

    implicit none

!    complex, dimension (:,:,-nzgrid:), intent (in out) :: phi, apar
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g

    integer :: ivmu, iv, sgn

    ! start the timer for the explicit part of the solve
    if (proc0) call time_message(.false.,time_gke(:,8),' explicit')

    select case (explicit_option_switch)
    case (explicit_option_rk2)
       ! SSP RK2
       call advance_explicit_rk2 (g)
    case (explicit_option_rk3)
       ! default is SSP RK3
       call advance_explicit_rk3 (g)
    case (explicit_option_rk4)
       ! RK4
       call advance_explicit_rk4 (g)
    end select
    
    ! enforce periodicity for zonal modes
    if (zonal_mode(1)) then
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          iv = iv_idx(vmu_lo,ivmu)
          ! stream_sign > 0 corresponds to dz/dt < 0
          sgn = stream_sign(iv)
          g(1,:,sgn*nzgrid,:,ivmu) = g(1,:,-sgn*nzgrid,:,ivmu)
       end do
    end if

    ! stop the timer for the explicit part of the solve
    if (proc0) call time_message(.false.,time_gke(:,8),' explicit')

  end subroutine advance_explicit

  ! strong stability-preserving RK2
  subroutine advance_explicit_rk2 (g)

    use dist_fn_arrays, only: g1, gold
    use zgrid, only: nzgrid
    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g

    integer :: icnt
    logical :: restart_time_step

    ! if CFL condition is violated by nonlinear term
    ! then must modify time step size and restart time step
    ! assume false and test
    restart_time_step = .false.

    gold = g

    icnt = 1
    ! SSP rk3 algorithm to advance explicit part of code
    ! if GK equation written as dg/dt = rhs - vpar . grad h,
    ! solve_gke returns rhs*dt
    do while (icnt <= 2)

       select case (icnt)
       case (1)
          call solve_gke (gold, g1, restart_time_step)
       case (2)
          g1 = gold + g1
          call solve_gke (g1, g, restart_time_step)
       end select
       if (restart_time_step) then
          icnt = 1
       else
          icnt = icnt + 1
       end if
    end do

    ! this is gbar at intermediate time level
    g = 0.5*gold + 0.5*(g1 + g)

  end subroutine advance_explicit_rk2

  ! strong stability-preserving RK3
  subroutine advance_explicit_rk3 (g)

    use dist_fn_arrays, only: g1, g2, gold
    use zgrid, only: nzgrid
    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g

    integer :: icnt
    logical :: restart_time_step

    ! if CFL condition is violated by nonlinear term
    ! then must modify time step size and restart time step
    ! assume false and test
    restart_time_step = .false.

    gold = g

    icnt = 1
    ! SSP rk3 algorithm to advance explicit part of code
    ! if GK equation written as dg/dt = rhs - vpar . grad h,
    ! solve_gke returns rhs*dt
    do while (icnt <= 3)
       select case (icnt)
       case (1)
          call solve_gke (gold, g1, restart_time_step)
       case (2)
          g1 = gold + g1
          call solve_gke (g1, g2, restart_time_step)
       case (3)
          g2 = g1 + g2
          call solve_gke (g2, g, restart_time_step)
       end select
       if (restart_time_step) then
          icnt = 1
       else
          icnt = icnt + 1
       end if
    end do

    ! this is gbar at intermediate time level
    g = gold/3. + 0.5*g1 + (g2 + g)/6.

  end subroutine advance_explicit_rk3

  subroutine advance_explicit_rk4 (g)

    use dist_fn_arrays, only: g1, g2, g3, gold
    use zgrid, only: nzgrid
    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g

    integer :: icnt
    logical :: restart_time_step

    ! if CFL condition is violated by nonlinear term
    ! then must modify time step size and restart time step
    ! assume false and test
    restart_time_step = .false.

    gold = g

    icnt = 1
    ! SSP rk3 algorithm to advance explicit part of code
    ! if GK equation written as dg/dt = rhs - vpar . grad h,
    ! solve_gke returns rhs*dt
    do while (icnt <= 4)
       select case (icnt)
       case (1)
          call solve_gke (gold, g1, restart_time_step)
       case (2)
          ! g1 is h*k1
          g3 = gold + 0.5*g1
          call solve_gke (g3, g2, restart_time_step)
          g1 = g1 + 2.*g2
       case (3)
          ! g2 is h*k2
          g2 = gold+0.5*g2
          call solve_gke (g2, g3, restart_time_step)
          g1 = g1 + 2.*g3
       case (4)
          ! g3 is h*k3
          g3 = gold+g3
          call solve_gke (g3, g, restart_time_step)
          g1 = g1 + g
       end select
       if (restart_time_step) then
          icnt = 1
       else
          icnt = icnt + 1
       end if
    end do

    ! this is gbar at intermediate time level
    g = gold + g1/6.

  end subroutine advance_explicit_rk4

  subroutine solve_gke (gin, rhs_ky, restart_time_step)

    use job_manage, only: time_message
    use fields_arrays, only: phi, apar
    use stella_layouts, only: vmu_lo
    use stella_transforms, only: transform_y2ky
    use redistribute, only: gather, scatter
    use physics_flags, only: include_parallel_nonlinearity
    use physics_flags, only: include_parallel_streaming
    use physics_flags, only: include_mirror
    use physics_flags, only: nonlinear
    use run_parameters, only: fphi, fapar
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, ny
    use physics_flags, only: full_flux_surface
    use run_parameters, only: stream_implicit, mirror_implicit
    use dissipation, only: include_collisions, advance_collisions_explicit, collisions_implicit
    use parallel_streaming, only: advance_parallel_streaming_explicit
    use fields, only: advance_fields, fields_updated
    use mirror_terms, only: advance_mirror_explicit

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: gin
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (out), target :: rhs_ky
    logical, intent (out) :: restart_time_step

    complex, dimension (:,:,:,:,:), allocatable, target :: rhs_y
    complex, dimension (:,:,:,:,:), pointer :: rhs

    rhs_ky = 0.

    if (full_flux_surface) then
       allocate (rhs_y(ny,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_proc))
       rhs_y = 0.
       rhs => rhs_y
    else
       rhs => rhs_ky
    end if

    ! start with gbar in k-space and (ky,kx,z) local
    ! obtain fields corresponding to gbar
    call advance_fields (gin, phi, apar, dist='gbar')

    ! default is to continue with same time step size
    ! if estimated CFL condition for nonlinear terms is violated
    ! then restart_time_step will be set to .true.
    restart_time_step = .false.
    ! calculate and add ExB nonlinearity to RHS of GK eqn
    ! do this first, as the CFL condition may require a change in time step
    ! and thus recomputation of mirror, wdrift, wstar, and parstream
    if (nonlinear) call advance_ExB_nonlinearity (gin, rhs, restart_time_step)

    if (include_parallel_nonlinearity .and. .not.restart_time_step) &
         call advance_parallel_nonlinearity (gin, rhs, restart_time_step)

    if (.not.restart_time_step) then

       ! calculate and add mirror term to RHS of GK eqn
       if (include_mirror.and..not.mirror_implicit) then
          call advance_mirror_explicit (gin, rhs)
       end if

       ! calculate and add alpha-component of magnetic drift term to RHS of GK eqn
       call advance_wdrifty_explicit (gin, phi, rhs)

       ! calculate and add psi-component of magnetic drift term to RHS of GK eqn
       call advance_wdriftx (gin, phi, rhs)

       ! calculate and add omega_* term to RHS of GK eqn
       call advance_wstar_explicit (rhs)

       if (include_collisions.and..not.collisions_implicit) call advance_collisions_explicit (gin, phi, rhs)
       
       if (full_flux_surface) then
          call transform_y2ky (rhs_y, rhs_ky)
          deallocate (rhs_y)
       end if
       
       ! calculate and add parallel streaming term to RHS of GK eqn
       if (include_parallel_streaming.and.(.not.stream_implicit)) &
            call advance_parallel_streaming_explicit (gin, rhs_ky)

       ! calculate and add omega_* term to RHS of GK eqn
!       if (wstar_explicit) call advance_wstar_explicit (rhs_ky)
!       call advance_wstar_explicit (rhs_ky)
       ! calculate and add collision term to RHS of GK eqn
       !    call advance_collisions
    end if

    fields_updated = .false.

    nullify (rhs)

  end subroutine solve_gke

  subroutine advance_wstar_explicit (gout)

    use mp, only: proc0, mp_abort
    use job_manage, only: time_message
    use fields_arrays, only: phi, apar
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx
    use physics_flags, only: full_flux_surface

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout

    complex, dimension (:,:,:,:,:), allocatable :: g0

    if (proc0) call time_message(.false.,time_gke(:,6),' wstar advance')

    allocate (g0(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    ! omega_* stays in ky,kx,z space with ky,kx,z local
    ! get d<chi>/dy
    if (debug) write (*,*) 'time_advance::solve_gke::get_dchidy'
    if (full_flux_surface) then
       ! FLAG -- ADD SOMETHING HERE
       call mp_abort ('wstar term not yet setup for full_flux_surface = .true. aborting.')
    else
       call get_dchidy (phi, apar, g0)
       ! multiply with omega_* coefficient and add to source (RHS of GK eqn)
       if (debug) write (*,*) 'time_advance::solve_gke::add_wstar_term'
       call add_wstar_term (g0, gout)
    end if
    deallocate (g0)

    if (proc0) call time_message(.false.,time_gke(:,6),' wstar advance')

  end subroutine advance_wstar_explicit

  subroutine advance_wdrifty_explicit (g, phi, gout)

    use mp, only: proc0
    use stella_layouts, only: vmu_lo
    use job_manage, only: time_message
    use stella_transforms, only: transform_ky2y
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky, ny
    use physics_flags, only: full_flux_surface
    use gyro_averages, only: gyro_average
    use dist_fn_arrays, only: wdrifty_g, wdrifty_phi

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phi
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout

    integer :: ivmu
    complex, dimension (:,:,:,:), allocatable :: dphidy
    complex, dimension (:,:,:,:,:), allocatable :: g0k, g0y

    ! alpha-component of magnetic drift (requires ky -> y)
    if (proc0) call time_message(.false.,time_gke(:,4),' dgdy advance')

    allocate (dphidy(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (g0k(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    if (debug) write (*,*) 'time_advance::solve_gke::get_dgdy'
    call get_dgdy (g, g0k)
    call get_dgdy (phi, dphidy)

    if (full_flux_surface) then
       allocate (g0y(ny,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       ! transform dg/dy from k-space to y-space
       call transform_ky2y (g0k, g0y)
       ! add vM . grad y dg/dy term to equation
       call add_dg_term_global (g0y, wdrifty_g, gout)
       ! get <dphi/dy> in k-space
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          call gyro_average (dphidy, ivmu, g0k(:,:,:,:,ivmu))
       end do
       ! transform d<phi>/dy from k-space to y-space
       call transform_ky2y (g0k, g0y)
       ! add vM . grad y d<phi>/dy term to equation
       call add_dphi_term_global (g0y, wdrifty_phi, gout)
       deallocate (g0y)
    else
       if (debug) write (*,*) 'time_advance::solve_gke::add_dgdy_term'
       ! add vM . grad y dg/dy term to equation
       call add_dg_term (g0k, wdrifty_g(1,:,:), gout)
       ! get <dphi/dy> in k-space
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          call gyro_average (dphidy, ivmu, g0k(:,:,:,:,ivmu))
       end do
       ! add vM . grad y d<phi>/dy term to equation
       call add_dphi_term (g0k, wdrifty_phi(1,:,:), gout)
    end if
    deallocate (g0k, dphidy)

    if (proc0) call time_message(.false.,time_gke(:,4),' dgdy advance')

  end subroutine advance_wdrifty_explicit

  subroutine advance_wdriftx (g, phi, gout)

    use mp, only: proc0
    use stella_layouts, only: vmu_lo
    use job_manage, only: time_message
    use stella_transforms, only: transform_ky2y
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky, ny, akx
    use physics_flags, only: full_flux_surface
    use gyro_averages, only: gyro_average
    use dist_fn_arrays, only: wdriftx_g, wdriftx_phi

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phi
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout

    integer :: ivmu
    complex, dimension (:,:,:,:), allocatable :: dphidx
    complex, dimension (:,:,:,:,:), allocatable :: g0k, g0y

    ! psi-component of magnetic drift (requires ky -> y)
    if (proc0) call time_message(.false.,time_gke(:,5),' dgdx advance')

    ! do not calculate if wdriftx terms are all zero
    if (maxval(abs(akx))<epsilon(0.)) then
       if (proc0) call time_message(.false.,time_gke(:,5),' dgdx advance')
       return
    end if

    allocate (dphidx(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (g0k(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    if (debug) write (*,*) 'time_advance::solve_gke::get_dgdx'
    call get_dgdx (g, g0k)
    call get_dgdx (phi, dphidx)

    if (full_flux_surface) then
       allocate (g0y(ny,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       ! transform dg/dx from k-space to y-space
       call transform_ky2y (g0k, g0y)
       ! add vM . grad x dg/dx term to equation
       call add_dg_term_global (g0y, wdriftx_g, gout)
       ! get <dphi/dx> in k-space
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          call gyro_average (dphidx, ivmu, g0k(:,:,:,:,ivmu))
       end do
       ! transform d<phi>/dx from k-space to y-space
       call transform_ky2y (g0k, g0y)
       ! add vM . grad x d<phi>/dx term to equation
       call add_dphi_term_global (g0y, wdriftx_phi, gout)
       deallocate (g0y)
    else
       if (debug) write (*,*) 'time_advance::solve_gke::add_dgdx_term'
       ! add vM . grad x dg/dx term to equation
       call add_dg_term (g0k, wdriftx_g(1,:,:), gout)
       ! get <dphi/dx> in k-space
       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          call gyro_average (dphidx, ivmu, g0k(:,:,:,:,ivmu))
       end do
       ! add vM . grad x d<phi>/dx term to equation
       call add_dphi_term (g0k, wdriftx_phi(1,:,:), gout)
    end if
    deallocate (g0k, dphidx)

    if (proc0) call time_message(.false.,time_gke(:,5),' dgdx advance')

  end subroutine advance_wdriftx

  subroutine advance_ExB_nonlinearity (g, gout, restart_time_step)

    use mp, only: proc0, min_allreduce
    use stella_layouts, only: vmu_lo
    use job_manage, only: time_message
    use fields_arrays, only: phi, apar
    use stella_transforms, only: transform_ky2y, transform_y2ky
    use stella_transforms, only: transform_kx2x, transform_x2kx
    use stella_time, only: cfl_dt, code_dt, code_dt_max
    use run_parameters, only: cfl_cushion, delt_adjust
    use zgrid, only: nzgrid, ntubes
    use stella_geometry, only: exb_nonlin_fac
    use kt_grids, only: nakx, naky, nx, ny, ikx_max
    use kt_grids, only: akx, aky
    use physics_flags, only: full_flux_surface
    use kt_grids, only: swap_kxky, swap_kxky_back
    use constants, only: pi

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout
    logical, intent (out) :: restart_time_step
    
    integer :: ivmu, iz, it

    complex, dimension (:,:), allocatable :: g0k, g0k_swap
    complex, dimension (:,:), allocatable :: g0kxy
    real, dimension (:,:), allocatable :: g0xy, g1xy, bracket

    ! alpha-component of magnetic drift (requires ky -> y)
    if (proc0) call time_message(.false.,time_gke(:,7),' ExB nonlinear advance')

    restart_time_step = .false.

    allocate (g0k(naky,nakx))
    allocate (g0k_swap(2*naky-1,ikx_max))
    allocate (g0kxy(ny,ikx_max))
    allocate (g0xy(ny,nx))
    allocate (g1xy(ny,nx))
    allocate (bracket(ny,nx))

    if (debug) write (*,*) 'time_advance::solve_gke::advance_ExB_nonlinearity::get_dgdy'
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do it = 1, ntubes
          do iz = -nzgrid, nzgrid
             call get_dgdy (g(:,:,iz,it,ivmu), g0k)
             call swap_kxky (g0k, g0k_swap)

             ! we have i*ky*g(kx,ky) for ky >= 0 and all kx
             ! want to do 1D complex to complex transform in y
             ! which requires i*ky*g(kx,ky) for all ky and kx >= 0
             ! use g(kx,-ky) = conjg(g(-kx,ky))
             ! so i*(-ky)*g(kx,-ky) = -i*ky*conjg(g(-kx,ky)) = conjg(i*ky*g(-kx,ky))
             ! and i*kx*g(kx,-ky) = i*kx*conjg(g(-kx,ky)) = conjg(i*(-kx)*g(-kx,ky))
             ! and i*(-ky)*J0(kx,-ky)*phi(kx,-ky) = conjg(i*ky*J0(-kx,ky)*phi(-kx,ky))
             ! and i*kx*J0(kx,-ky)*phi(kx,-ky) = conjg(i*(-kx)*J0(-kx,ky)*phi(-kx,ky))
             ! i.e., can calculate dg/dx, dg/dy, d<phi>/dx and d<phi>/dy 
             ! on stella (kx,ky) grid, then conjugate and flip sign of (kx,ky)
             ! NB: J0(kx,ky) = J0(-kx,-ky)
             call transform_ky2y (g0k_swap, g0kxy)
             call transform_kx2x (g0kxy, g0xy)

             call get_dchidx (iz, ivmu, phi(:,:,iz,it), apar(:,:,iz,it), g0k)
             call swap_kxky (g0k, g0k_swap)
             call transform_ky2y (g0k_swap, g0kxy)
             call transform_kx2x (g0kxy, g1xy)
             g1xy = g1xy*exb_nonlin_fac
             bracket = g0xy*g1xy
             cfl_dt = min(cfl_dt,2.*pi/(maxval(abs(g1xy))*aky(naky)))

             call get_dgdx (g(:,:,iz,it,ivmu), g0k)
             call swap_kxky (g0k, g0k_swap)
             call transform_ky2y (g0k_swap, g0kxy)
             call transform_kx2x (g0kxy, g0xy)
             call get_dchidy (iz, ivmu, phi(:,:,iz,it), apar(:,:,iz,it), g0k)
             call swap_kxky (g0k, g0k_swap)
             call transform_ky2y (g0k_swap, g0kxy)
             call transform_kx2x (g0kxy, g1xy)
             g1xy = g1xy*exb_nonlin_fac
             bracket = bracket - g0xy*g1xy
             cfl_dt = min(cfl_dt,2.*pi/(maxval(abs(g1xy))*akx(ikx_max)))

             call transform_x2kx (bracket, g0kxy)

             if (full_flux_surface) then
                gout(:,:,iz,it,ivmu) = g0kxy
             else
                call transform_y2ky (g0kxy, g0k_swap)
                call swap_kxky_back (g0k_swap, gout(:,:,iz,it,ivmu))
             end if
          end do
       end do
       ! enforce periodicity for zonal mode
       ! FLAG -- THIS IS PROBABLY NOT NECESSARY (DONE AT THE END OF EXPLICIT ADVANCE)
       ! AND MAY INDEED BE THE WRONG THING TO DO
       gout(1,:,-nzgrid,:,ivmu) = 0.5*(gout(1,:,nzgrid,:,ivmu)+gout(1,:,-nzgrid,:,ivmu))
       gout(1,:,nzgrid,:,ivmu) = gout(1,:,-nzgrid,:,ivmu)
    end do
    deallocate (g0k, g0k_swap, g0kxy, g0xy, g1xy, bracket)

    call min_allreduce (cfl_dt)

    if (code_dt > cfl_dt*cfl_cushion) then
       if (proc0) then
          write (*,*) 'code_dt= ', code_dt, 'larger than cfl_dt*cfl_cushion= ', cfl_dt*cfl_cushion
          write (*,*) 'setting code_dt=cfl_dt*cfl_cushion/delt_adjust and restarting time step'
       end if
       code_dt = cfl_dt*cfl_cushion/delt_adjust
       call reset_dt
       restart_time_step = .true.
    else if (code_dt < min(cfl_dt*cfl_cushion/delt_adjust,code_dt_max)) then
       if (proc0) then
          write (*,*) 'code_dt= ', code_dt, 'smaller than cfl_dt*cfl_cushion/delt_adjust= ', cfl_dt*cfl_cushion/delt_adjust
          write (*,*) 'setting code_dt=min(cfl_dt*cfl_cushion/delt_adjust,delt) and restarting time step'
       end if
       code_dt = min(cfl_dt*cfl_cushion/delt_adjust,code_dt_max)
       call reset_dt
       ! FLAG -- NOT SURE THIS IS CORRECT
       gout = code_dt*gout
    else
       gout = code_dt*gout
    end if

    if (proc0) call time_message(.false.,time_gke(:,7),' ExB nonlinear advance')

  end subroutine advance_ExB_nonlinearity

  subroutine advance_parallel_nonlinearity (g, gout, restart_time_step)

    use constants, only: zi
    use mp, only: proc0, min_allreduce
    use stella_layouts, only: vmu_lo, xyz_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use job_manage, only: time_message
    use finite_differences, only: second_order_centered_zed
    use finite_differences, only: third_order_upwind
    use redistribute, only: gather, scatter
    use fields_arrays, only: phi
    use stella_transforms, only: transform_ky2y, transform_y2ky
    use stella_transforms, only: transform_kx2x, transform_x2kx
    use stella_time, only: cfl_dt, code_dt, code_dt_max
    use run_parameters, only: cfl_cushion, delt_adjust
    use zgrid, only: nzgrid, delzed, ntubes
    use extended_zgrid, only: neigen, nsegments, ikxmod
    use extended_zgrid, only: iz_low, iz_up
    use physics_flags, only: full_flux_surface
    use kt_grids, only: nakx, naky, nx, ny, ikx_max
    use kt_grids, only: zonal_mode
    use kt_grids, only: swap_kxky, swap_kxky_back
    use vpamu_grids, only: nvpa, nmu
    use vpamu_grids, only: dvpa, vpa, mu
    use gyro_averages, only: gyro_average
    use parallel_streaming, only: stream_sign
    use dist_redistribute, only: xyz2vmu
    use extended_zgrid, only: fill_zed_ghost_zones

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout
    logical, intent (out) :: restart_time_step
    
    integer :: ivmu, ixyz
    integer :: iz, it, iv, imu, is
    integer :: iky, ie, iseg
    integer :: advect_sign
    real, dimension (:), allocatable :: dgdv
    real, dimension (:,:,:,:,:), allocatable :: g0xy
    real, dimension (:,:,:), allocatable :: gxy_vmulocal
    real, dimension (:,:), allocatable :: advect_speed
    complex, dimension (2) :: gleft, gright
    complex, dimension (:,:,:,:), allocatable :: phi_gyro, dphi
    complex, dimension (:,:), allocatable :: g0k, g0kxy, dphi_swap
    complex, dimension (:,:), allocatable :: tmp

    ! alpha-component of magnetic drift (requires ky -> y)
    if (proc0) call time_message(.false.,time_parallel_nl(:,1),' parallel nonlinearity advance')

    restart_time_step = .false.

    ! overview:
    ! need g and d<phi>/dz in (x,y) space in
    ! order to upwind dg/dvpa
    ! 1) transform d<phi>/dz from (kx,ky) to (x,y). layout: vmu_lo
    ! 2) need sign of parnl advection in xyz_lo (since dg/dvpa 
    !    requires vpa local), so d<phi>/dz(vmu_lo) --> d<phi>/dz(xyz_lo)
    ! 3) transform g from (kx,ky) to (x,y). layout: vmu_lo
    ! 4) dg/dvpa requires vpa local, so g(vmu_lo) --> g(xyz_lo)
    ! 5) calculate dg/dvpa
    ! 6) multiply dg/dvpa with d<phi>/dz
    ! 7) product(xyz_lo) --> product(vmu_lo)
    ! 8) inverse transform product(vmu_lo)

    allocate (g0xy(ny,nx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
    allocate (g0kxy(ny,ikx_max))
    allocate (phi_gyro(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (dphi(naky,nakx,-nzgrid:nzgrid,ntubes))
    allocate (dphi_swap(2*naky-1,ikx_max))
    allocate (tmp(size(gout,1),size(gout,2)))

    ! get d<phi>/dz in vmu_lo
    ! we will need to transform it to real-space
    ! as its sign is needed for upwinding of dg/dvpa
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       iv = iv_idx(vmu_lo,ivmu)
       imu = imu_idx(vmu_lo,ivmu)
       is = is_idx(vmu_lo,ivmu)
       ! construct <phi>
       call gyro_average (phi, ivmu, phi_gyro)
       do iky = 1, naky
          do it = 1, ntubes
             do ie = 1, neigen(iky)
                do iseg = 1, nsegments(ie,iky)
                   ! first fill in ghost zones at boundaries in g(z)
                   call fill_zed_ghost_zones (it, iseg, ie, iky, phi_gyro, gleft, gright)
                   ! now get d<phi>/dz
                   call second_order_centered_zed (iz_low(iseg), iseg, nsegments(ie,iky), &
                        phi_gyro(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it), &
                        delzed(0), stream_sign(iv), gleft, gright, zonal_mode(iky), &
                        dphi(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it))
                end do
             end do
          end do
       end do

       dphi = dphi*spread(spread(spread(par_nl_fac(:,is),1,naky),2,nakx),4,ntubes) &
            + spread(spread(spread(par_nl_curv(:,is),1,naky),2,nakx),4,ntubes)*vpa(iv)*mu(imu) &
            + zi*vpa(iv)*phi_gyro*spread(spread(par_nl_driftx,1,naky)+spread(par_nl_drifty,2,nakx),4,ntubes)

       do it = 1, ntubes
          do iz = -nzgrid, nzgrid
             ! use reality to swap from ky >= 0, all kx to kx >= 0 , all ky
             call swap_kxky (dphi(:,:,iz,it), dphi_swap)
             ! transform in y
             call transform_ky2y (dphi_swap, g0kxy)
             ! transform in x
             call transform_kx2x (g0kxy, g0xy(:,:,iz,it,ivmu))
!          ! get dvpa/dt
!          g0xy(:,:,iz,ivmu) = g0xy(:,:,iz,ivmu)*par_nl_fac(iz,is) + par_nl_curv(iz,is)*vpa(iv)*mu(imu)
          end do
       end do
    end do

    ! do not need phi_gyro or dphidz  again so deallocate
    deallocate (phi_gyro)
    deallocate (dphi, dphi_swap)

    allocate (gxy_vmulocal(nvpa,nmu,xyz_lo%llim_proc:xyz_lo%ulim_alloc))
    allocate (advect_speed(nmu,xyz_lo%llim_proc:xyz_lo%ulim_alloc))

    ! we now have the advection velocity in vpa in (x,y) space
    ! next redistribute it so that (vpa,mu) are local
    if (proc0) call time_message(.false.,time_parallel_nl(:,2),' parallel nonlinearity redist')
    call scatter (xyz2vmu, g0xy, gxy_vmulocal)
    if (proc0) call time_message(.false.,time_parallel_nl(:,2),' parallel nonlinearity redist')
    ! advect_speed does not depend on vpa
    advect_speed = gxy_vmulocal(1,:,:)

    allocate (g0k(2*naky-1,ikx_max))

    ! transform g from (kx,ky) to (x,y)
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do it = 1, ntubes
          do iz = -nzgrid, nzgrid
             call swap_kxky (g(:,:,iz,it,ivmu), g0k)
             ! transform in y
             call transform_ky2y (g0k, g0kxy)
             ! transform in x
             call transform_kx2x (g0kxy, g0xy(:,:,iz,it,ivmu))
          end do
       end do
    end do

    ! finished with g0k
    deallocate (g0k)

    ! redistribute so that (vpa,mu) local
    if (proc0) call time_message(.false.,time_parallel_nl(:,2),' parallel nonlinearity redist')
    call scatter (xyz2vmu, g0xy, gxy_vmulocal)
    if (proc0) call time_message(.false.,time_parallel_nl(:,2),' parallel nonlinearity redist')
    
    allocate (dgdv(nvpa))

    ! we now need to form dg/dvpa and obtain product of dg/dvpa with advection speed
    do ixyz = xyz_lo%llim_proc, xyz_lo%ulim_proc
       do imu = 1, nmu
          ! advect_sign set to +/- 1 depending on sign of the parallel nonlinearity 
          ! advection velocity
          ! NB: advect_sign = -1 corresponds to positive advection velocity
          advect_sign = int(sign(1.0,advect_speed(imu,ixyz)))
          call third_order_upwind (1,gxy_vmulocal(:,imu,ixyz),dvpa,advect_sign,dgdv)
          gxy_vmulocal(:,imu,ixyz) = dgdv*advect_speed(imu,ixyz)
          cfl_dt = min(cfl_dt,dvpa/abs(advect_speed(imu,ixyz)))
       end do
    end do

    ! finished with dgdv and advect_speed
    deallocate (dgdv, advect_speed)
    
    ! now that we have the full parallel nonlinearity in (x,y)-space
    ! need to redistribute so that (x,y) local for transforms
    if (proc0) call time_message(.false.,time_parallel_nl(:,2),' parallel nonlinearity redist')
    call gather (xyz2vmu, gxy_vmulocal, g0xy)
    if (proc0) call time_message(.false.,time_parallel_nl(:,2),' parallel nonlinearity redist')

    ! finished with gxy_vmulocal
    deallocate (gxy_vmulocal)

    allocate (g0k(2*naky-1,ikx_max))

    ! g0xy is parallel nonlinearity term with (x,y) on processor
    ! need to inverse Fourier transform
    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       do it = 1, ntubes
          do iz = -nzgrid, nzgrid
             call transform_x2kx (g0xy(:,:,iz,it,ivmu), g0kxy)
             if (full_flux_surface) then
                gout(:,:,iz,it,ivmu) = gout(:,:,iz,it,ivmu) + code_dt*g0kxy
             else
                call transform_y2ky (g0kxy, g0k)
                call swap_kxky_back (g0k, tmp)
                gout(:,:,iz,it,ivmu) = gout(:,:,iz,it,ivmu) + code_dt*tmp
             end if
          end do
       end do
    end do
    deallocate (g0k, g0kxy, g0xy)

    call min_allreduce (cfl_dt)

    if (code_dt > cfl_dt*cfl_cushion) then
       if (proc0) then
          write (*,*) 'code_dt= ', code_dt, 'larger than cfl_dt*cfl_cushion= ', cfl_dt
          write (*,*) 'setting code_dt=cfl_dt*cfl_cushion/delt_adjust and restarting time step'
       end if
       code_dt = cfl_dt*cfl_cushion/delt_adjust
       call reset_dt
       restart_time_step = .true.
    else if (code_dt < min(cfl_dt*cfl_cushion/delt_adjust,code_dt_max)) then
       if (proc0) then
          write (*,*) 'code_dt= ', code_dt, 'smaller than cfl_dt*cfl_cushion/delt_adjust= ', cfl_dt*cfl_cushion/delt_adjust
          write (*,*) 'setting code_dt=min(cfl_dt*cfl_cushion/delt_adjust,delt) and restarting time step'
       end if
       code_dt = min(cfl_dt*cfl_cushion/delt_adjust,code_dt_max)
       call reset_dt
!    else
!       gout = code_dt*gout
    end if

    if (proc0) call time_message(.false.,time_parallel_nl(:,1),' parallel nonlinearity advance')

  end subroutine advance_parallel_nonlinearity

  subroutine get_dgdy_2d (g, dgdy)

    use constants, only: zi
    use kt_grids, only: nakx, aky

    implicit none

    complex, dimension (:,:), intent (in) :: g
    complex, dimension (:,:), intent (out) :: dgdy

    dgdy = zi*spread(aky,2,nakx)*g

  end subroutine get_dgdy_2d

  subroutine get_dgdy_3d (g, dgdy)

    use constants, only: zi
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, aky

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: dgdy

    dgdy = zi*spread(spread(spread(aky,2,nakx),3,2*nzgrid+1),4,ntubes)*g
    
  end subroutine get_dgdy_3d

  subroutine get_dgdy_4d (g, dgdy)

    use constants, only: zi
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, aky

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (out) :: dgdy

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       dgdy(:,:,:,:,ivmu) = zi*spread(spread(spread(aky,2,nakx),3,2*nzgrid+1),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine get_dgdy_4d

  subroutine add_dg_term (g, wdrift_in, src)

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (-nzgrid:,vmu_lo%llim_proc:), intent (in) :: wdrift_in
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,:,ivmu) = src(:,:,:,:,ivmu) &
            + spread(spread(spread(wdrift_in(:,ivmu),1,naky),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_dg_term

  subroutine add_dg_term_global (g, wdrift_in, src)

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: wdrift_in
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,:,ivmu) = src(:,:,:,:,ivmu) - spread(spread(wdrift_in(:,:,ivmu),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_dg_term_global

  subroutine get_dgdx_2d (g, dgdx)

    use constants, only: zi
    use kt_grids, only: naky, akx

    implicit none

    complex, dimension (:,:), intent (in) :: g
    complex, dimension (:,:), intent (out) :: dgdx

    dgdx = zi*spread(akx,1,naky)*g

  end subroutine get_dgdx_2d

  subroutine get_dgdx_3d (g, dgdx)

    use constants, only: zi
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, akx
    
    implicit none
    
    complex, dimension (:,:,-nzgrid:,:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:), intent (out) :: dgdx

    dgdx = zi*spread(spread(spread(akx,1,naky),3,2*nzgrid+1),4,ntubes)*g

  end subroutine get_dgdx_3d

  subroutine get_dgdx_4d (g, dgdx)

    use constants, only: zi
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, akx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (out) :: dgdx

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       dgdx(:,:,:,:,ivmu) = zi*spread(spread(spread(akx,1,naky),3,2*nzgrid+1),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine get_dgdx_4d

  subroutine add_dphi_term (g, wdrift, src)

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (-nzgrid:,vmu_lo%llim_proc:), intent (in) :: wdrift
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,:,ivmu) = src(:,:,:,:,ivmu) &
            + spread(spread(spread(wdrift(:,ivmu),1,naky),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_dphi_term

  subroutine add_dphi_term_global (g, wdrift, src)

    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    real, dimension (:,-nzgrid:,vmu_lo%llim_proc:), intent (in) :: wdrift
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,:,ivmu) = src(:,:,:,:,ivmu) &
            + spread(spread(wdrift(:,:,ivmu),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_dphi_term_global

  subroutine get_dchidy_4d (phi, apar, dchidy)

    use constants, only: zi
    use gyro_averages, only: gyro_average
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: is_idx, iv_idx
    use run_parameters, only: fphi, fapar
    use species, only: spec
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: vpa
    use kt_grids, only: nakx, aky, naky

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: phi, apar
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (out) :: dchidy

    integer :: ivmu, iv, is
    complex, dimension (:,:,:,:), allocatable :: field

    allocate (field(naky,nakx,-nzgrid:nzgrid,ntubes))

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       is = is_idx(vmu_lo,ivmu)
       iv = iv_idx(vmu_lo,ivmu)
       field = zi*spread(spread(spread(aky,2,nakx),3,2*nzgrid+1),4,ntubes) &
            * ( fphi*phi - fapar*vpa(iv)*spec(is)%stm*apar )
       call gyro_average (field, ivmu, dchidy(:,:,:,:,ivmu))
    end do

    deallocate (field)

  end subroutine get_dchidy_4d

  subroutine get_dchidy_2d (iz, ivmu, phi, apar, dchidy)

    use constants, only: zi
    use gyro_averages, only: gyro_average
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: is_idx, iv_idx
    use run_parameters, only: fphi, fapar
    use species, only: spec
    use vpamu_grids, only: vpa
    use kt_grids, only: nakx, aky, naky

    implicit none

    integer, intent (in) :: ivmu, iz
    complex, dimension (:,:), intent (in) :: phi, apar
    complex, dimension (:,:), intent (out) :: dchidy

    integer :: iv, is
    complex, dimension (:,:), allocatable :: field

    allocate (field(naky,nakx))

    is = is_idx(vmu_lo,ivmu)
    iv = iv_idx(vmu_lo,ivmu)
    field = zi*spread(aky,2,nakx) &
         * ( fphi*phi - fapar*vpa(iv)*spec(is)%stm*apar )
    call gyro_average (field, iz, ivmu, dchidy)
    
    deallocate (field)

  end subroutine get_dchidy_2d

  subroutine get_dchidx (iz, ivmu, phi, apar, dchidx)

    use constants, only: zi
    use gyro_averages, only: gyro_average
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: is_idx, iv_idx
    use run_parameters, only: fphi, fapar
    use species, only: spec
    use vpamu_grids, only: vpa
    use kt_grids, only: akx, naky, nakx

    implicit none

    integer, intent (in) :: ivmu, iz
    complex, dimension (:,:), intent (in) :: phi, apar
    complex, dimension (:,:), intent (out) :: dchidx

    integer :: iv, is
    complex, dimension (:,:), allocatable :: field

    allocate (field(naky,nakx))

    is = is_idx(vmu_lo,ivmu)
    iv = iv_idx(vmu_lo,ivmu)
    field = zi*spread(akx,1,naky) &
         * ( fphi*phi - fapar*vpa(iv)*spec(is)%stm*apar )
    call gyro_average (field, iz, ivmu, dchidx)
    
    deallocate (field)

  end subroutine get_dchidx

  subroutine add_wstar_term (g, src)

    use dist_fn_arrays, only: wstar
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: src

    integer :: ivmu

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       src(:,:,:,:,ivmu) = src(:,:,:,:,ivmu) &
            + spread(spread(spread(wstar(1,:,ivmu),1,naky),2,nakx),4,ntubes)*g(:,:,:,:,ivmu)
    end do

  end subroutine add_wstar_term

  subroutine advance_implicit (istep, phi, apar, g)
!  subroutine advance_implicit (phi, apar, g)

    use mp, only: proc0
    use job_manage, only: time_message
    use run_parameters, only: fphi, fapar
    use stella_layouts, only: vmu_lo
    use zgrid, only: nzgrid
    use dissipation, only: hyper_dissipation, advance_hyper_dissipation
    use physics_flags, only: include_parallel_streaming
    use physics_flags, only: include_mirror
    use run_parameters, only: stream_implicit, mirror_implicit
    use parallel_streaming, only: advance_parallel_streaming_implicit
    use fields, only: advance_fields, fields_updated
    use mirror_terms, only: advance_mirror_implicit
    use dissipation, only: collisions_implicit, include_collisions
    use dissipation, only: advance_collisions_implicit
    use run_parameters, only: driftkinetic_implicit

    implicit none

    integer, intent (in) :: istep
    complex, dimension (:,:,-nzgrid:,:), intent (in out) :: phi, apar
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
!    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out), target :: g

!    complex, dimension (:,:,:,:,:), pointer :: gk, gy
!    complex, dimension (:,:,:,:,:), allocatable, target :: g_dual

!    ! the 'g' that enters this subroutine may be in alpha-space or kalpha-space
!    ! figure out which it is
!    if (size(g,1) == naky) then
!       alpha_space = .false.
!       gk => g
!       if (full_flux_surface) then
!          allocate (g_dual(nalpha,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
!          gy => g_dual
!       end if
!    else
!       alpha_space = .true.
!       allocate (g_dual(naky,nakx,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
!       gy => g
!       gk => g_dual
!    end if

    ! start the timer for the implicit part of the solve
    if (proc0) call time_message(.false.,time_gke(:,9),' implicit')

    ! reverse the order of operations every time step
    ! as part of alternating direction operator splitting
    ! this is needed to ensure 2nd order accuracy in time
!    if (mod(istep,2)==0) then
       ! g^{*} (coming from explicit solve) is input
       ! get g^{**}, with g^{**}-g^{*} due to mirror term

    if (mod(istep,2)==1 .or. .not.flip_flop) then

       if (hyper_dissipation) then
!          ! for hyper-dissipation, need to be in k-alpha space
!          if (alpha_space) call transform_y2ky (gy, gk)
          call advance_hyper_dissipation (g)
          fields_updated = .false.
       end if

       if (collisions_implicit .and. include_collisions) then
          call advance_fields (g, phi, apar, dist='gbar')
          call advance_collisions_implicit (mirror_implicit, phi, apar, g)
          fields_updated = .false.
       end if

       if (mirror_implicit .and. include_mirror) then
!          if (full_flux_surface) then
!             allocate (gy(ny,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
!             if (.not.alpha_space) call transform_ky2y (g, gy)
!          else
!             g_mirror => g
!          end if
          call advance_mirror_implicit (collisions_implicit, g)
          fields_updated = .false.
       end if
       
       ! get updated fields corresponding to advanced g
       ! note that hyper-dissipation and mirror advances
       ! depended only on g and so did not need field update
       call advance_fields (g, phi, apar, dist='gbar')

       ! g^{**} is input
       ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
       if ((stream_implicit.or.driftkinetic_implicit) .and. include_parallel_streaming) &
            call advance_parallel_streaming_implicit (g, phi, apar)

    else

       ! get updated fields corresponding to advanced g
       ! note that hyper-dissipation and mirror advances
       ! depended only on g and so did not need field update
       call advance_fields (g, phi, apar, dist='gbar')

       ! g^{**} is input
       ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
       if ((stream_implicit.or.driftkinetic_implicit) .and. include_parallel_streaming) &
            call advance_parallel_streaming_implicit (g, phi, apar)
       
       if (mirror_implicit .and. include_mirror) then
          call advance_mirror_implicit (collisions_implicit, g)
          fields_updated = .false.
       end if

       if (collisions_implicit .and. include_collisions) then
          call advance_fields (g, phi, apar, dist='gbar')
          call advance_collisions_implicit (mirror_implicit, phi, apar, g)
          fields_updated = .false.
       end if

       if (hyper_dissipation) then
          call advance_hyper_dissipation (g)
          fields_updated = .false.
       end if

    end if

!       call checksum (phi, phitot)
!       call checksum (g, gtot)
!       if (proc0) write (*,*) 'stream', phitot, gtot
       
!!       if (wdrifty_implicit) then
!!          call advance_wdrifty_implicit (g)
!!          call advance_fields (g, phi, apar, dist='gbar')
!!       end if
!!    if (wstar_implicit) call advance_wstar_implicit (g, phi, apar)
!    else
!!       if (wdrifty_implicit) then
!!          call advance_wdrifty_implicit (g)
!!          call advance_fields (g, phi, apar, dist='gbar')
!!       end if
!!    if (wstar_implicit) call advance_wstar_implicit (g, phi, apar)

!        ! g^{**} is input
!        ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
!        if (stream_implicit) call advance_parallel_streaming_implicit (g, phi, apar)

!        ! g^{*} (coming from explicit solve) is input
!        ! get g^{**}, with g^{**}-g^{*} due to mirror term
!        if (mirror_implicit) then
!           call advance_mirror_implicit (g)
!           ! get updated fields corresponding to mirror-advanced g
!           call advance_fields (g, phi, apar, dist='gbar')
!        end if
!     end if

    ! stop the timer for the implict part of the solve
    if (proc0) call time_message(.false.,time_gke(:,9),' implicit')

  end subroutine advance_implicit

!   subroutine advance_wdrifty_implicit (g)

!     use constants, only: zi
!     use stella_layouts, only: vmu_lo
!     use zgrid, only: nzgrid
!     use kt_grids, only: aky, naky
!     use dist_fn_arrays, only: wdrifty_g

!     implicit none

!     complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: g

!     integer :: ivmu, iky, iz
!     complex :: tmp

!     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!        do iz = -nzgrid, nzgrid
!           do iky = 1, naky
!              tmp = 0.5*zi*aky(iky)*wdrifty_g(1,iz,ivmu)
!              g(iky,:,iz,ivmu) = g(iky,:,iz,ivmu) * (1.0 + tmp) / (1.0 - tmp)
! !             g(iky,:,iz,ivmu) = g(iky,:,iz,ivmu) * (1.0 + 2.0*tmp)
!           end do
!        end do
!     end do

!   end subroutine advance_wdrifty_implicit

!   subroutine advance_wstar_implicit (g, phi, apar)

!     use constants, only: zi
!     use stella_layouts, only: vmu_lo
!     use stella_layouts, only: iv_idx, is_idx
!     use run_parameters, only: fphi, fapar
!     use zgrid, only: nzgrid
!     use dist_fn_arrays, only: gstar_to_g

!     implicit none

!     complex, dimension (:,:,-nzgrid:,vmu_lo%llim_proc:), intent (in out) :: g
!     complex, dimension (:,:,-nzgrid:), intent (in out) :: phi, apar

!     ! given g^{*}, obtain phi^{*} and apar^{*}
! !    call advance_fields (g, phi, apar, dist='gbar')

!     ! solve g^{**} = g^{*}+i*wstar*ky*J0*(chi^{**}+chi^{*})/2
!     ! define gstar^{**} = g^{**} - i*wstar*ky*J0*chi^{**}/2
!     ! so that gstar^{**} = g^{*} + i*wstar*ky*J0*chi^{*}/2

!     ! BACKWARDS DIFFERENCE FLAG
! !     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
! !        iv = iv_idx(vmu_lo,ivmu)
! !        is = is_idx(vmu_lo,ivmu)
! !        do iz = -nzgrid, nzgrid
! !           do ikx = 1, nakx
! !              do iky = 1, naky
! !                 g(iky,ikx,iz,ivmu) = g(iky,ikx,iz,ivmu) &
! !                      + zi*aky(iky)*wstar(iz,ivmu)*aj0x(iky,ikx,iz,ivmu) &
! !                      * (fphi*phi(iky,ikx,iz) - fapar*vpa(iv)*spec(is)%stm*apar(iky,ikx,iz))
! !              end do
! !           end do
! !        end do
! !     end do

!     ! now that we have gstar^{**}=g^{*}, obtain corresponding fields, phi^{**} and apar^{**}
!     call advance_fields (g, phi, apar, dist='gstar')
!     ! convert from gstar^{**} to g^{**}
!     call gstar_to_g (g, phi, apar, fphi, fapar)

!   end subroutine advance_wstar_implicit

  subroutine checksum_field (field, total)

    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky
    use extended_zgrid, only: neigen, nsegments, ikxmod
    use extended_zgrid, only: iz_low, iz_up

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: field
    real, intent (out) :: total

    integer :: it, iky, ie, iseg
    integer :: ikx

    total = 0.

    do iky = 1, naky
       do it = 1, ntubes
          do ie = 1, neigen(iky)
             iseg = 1
             ikx = ikxmod(iseg,ie,iky)
             total = total + sum(cabs(field(iky,ikx,iz_low(iseg):iz_up(iseg),it)))
             if (nsegments(ie,iky) > 1) then
                do iseg = 2, nsegments(ie,iky)
                   ikx = ikxmod(iseg,ie,iky)
                   total = total + sum(cabs(field(iky,ikx,iz_low(iseg)+1:iz_up(iseg),it)))
                end do
             end if
          end do
       end do
    end do

  end subroutine checksum_field

  subroutine checksum_dist (dist, total)

    use mp, only: sum_allreduce
    use zgrid, only: nzgrid
    use stella_layouts, only: vmu_lo

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: dist
    real, intent (out) :: total

    integer :: ivmu
    real :: subtotal

    total = 0.

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
       call checksum (dist(:,:,:,:,ivmu), subtotal)
       total = total + subtotal
    end do

    call sum_allreduce (total)

  end subroutine checksum_dist

!   subroutine init_get_fields_wstar

!     use constants, only: zi
!     use mp, only: sum_allreduce, mp_abort
!     use stella_layouts, only: kxkyz_lo
!     use stella_layouts, onlY: iz_idx, ikx_idx, iky_idx, is_idx
!     use dist_fn_arrays, only: aj0v
!     use run_parameters, only: fphi, fapar
!     use run_parameters, only: tite, nine, beta
!     use stella_time, only: code_dt
!     use species, only: spec, has_electron_species
!     use stella_geometry, only: dl_over_b
!     use zgrid, only: nzgrid
!     use vpamu_grids, only: nvpa, nvgrid, nmu
!     use vpamu_grids, only: energy
!     use vpamu_grids, only: maxwell_vpa, integrate_vmu
!     use species, only: spec
!     use kt_grids, only: naky, nakx, aky, akx

!     implicit none

!     integer :: ikxkyz, ig, ikx, iky, is
!     complex :: tmp
!     complex, dimension (:,:), allocatable :: g0

!     if (get_fields_wstar_initialized) return
!     get_fields_wstar_initialized = .true.

!     if (.not.allocated(gamtot_wstar)) &
!          allocate (gamtot_wstar(naky,nakx,-nzgrid:nzgrid))
!     gamtot_wstar = 0.
!     if (.not.allocated(gamtot3_wstar)) &
!          allocate (gamtot3_wstar(nakx,-nzgrid:nzgrid))
!     gamtot3_wstar = 0.
!     if (.not.allocated(apar_denom_wstar)) &
!          allocate (apar_denom_wstar(naky,nakx,-nzgrid:nzgrid))
!     apar_denom_wstar = 0.

!     if (fphi > epsilon(0.0)) then
!        allocate (g0(-nvgrid:nvgrid,nmu))
!        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
!           iky = iky_idx(kxkyz_lo,ikxkyz)
!           ikx = ikx_idx(kxkyz_lo,ikxkyz)
!           ig = iz_idx(kxkyz_lo,ikxkyz)
!           is = is_idx(kxkyz_lo,ikxkyz)
!           g0 = -zi*aky(iky)*spec(is)%z*spec(is)%dens*spread(aj0v(:,ikxkyz)**2,1,nvpa) &
!                ! BACKWARDS DIFFERENCE FLAG
! !               *anon(ig,:,:)*0.25*code_dt &
!                *spread(maxwell_vpa,2,nmu)*0.5*code_dt &
!                *(spec(is)%fprim+spec(is)%tprim*(energy(ig,:,:)-1.5))
!           call integrate_vmu (g0, ig, tmp)
!           gamtot_wstar(iky,ikx,ig) = gamtot_wstar(iky,ikx,ig) + tmp
!        end do
!        call sum_allreduce (gamtot_wstar)

!        gamtot_wstar = gamtot_wstar + gamtot

!        deallocate (g0)

!        if (.not.has_electron_species(spec)) then
!           ! no need to do anything extra for ky /= 0 because
!           ! already accounted for in gamtot_h
!           if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
!              if (abs(aky(1)) < epsilon(0.)) then
!                 do ikx = 1, nakx
!                    ! do not need kx=ky=0 mode
!                    if (abs(akx(ikx)) < epsilon(0.)) cycle
!                    tmp = nine/tite-sum(dl_over_b/gamtot_wstar(1,ikx,:))
!                    gamtot3_wstar(ikx,:) = 1./(gamtot_wstar(1,ikx,:)*tmp)
!                 end do
!              end if
!           end if
!        end if
!     end if

!     ! FLAG -- NEED TO SORT OUT FINITE FAPAR FOR GSTAR
!      if (fapar > epsilon(0.)) then
!         write (*,*) 'APAR NOT SETUP FOR GSTAR YET. aborting'
!         call mp_abort ('APAR NOT SETUP FOR GSTAR YET. aborting')
!      end if
        
! !        allocate (g0(-nvgrid:nvgrid,nmu))
! !        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
! !           iky = iky_idx(kxkyz_lo,ikxkyz)
! !           ikx = ikx_idx(kxkyz_lo,ikxkyz)
! !           ig = iz_idx(kxkyz_lo,ikxkyz)
! !           is = is_idx(kxkyz_lo,ikxkyz)
! !           g0 = spread(vpa**2,2,nmu)*spread(aj0v(:,ikxkyz)**2,1,nvpa)*anon(ig,:,:)
! !           wgt = 2.0*beta*spec(is)%z*spec(is)%z*spec(is)%dens/spec(is)%mass
! !           call integrate_vmu (g0, ig, tmp)
! !           apar_denom(iky,ikx,ig) = apar_denom(iky,ikx,ig) + tmp*wgt
! !        end do
! !        call sum_allreduce (apar_denom)
! !        apar_denom = apar_denom + kperp2

! !        deallocate (g0)
! !     end if
    
!   end subroutine init_get_fields_wstar

  subroutine finish_time_advance

    use stella_transforms, only: finish_transforms
    use physics_flags, only: full_flux_surface
    use extended_zgrid, only: finish_extended_zgrid
    use parallel_streaming, only: finish_parallel_streaming
    use mirror_terms, only: finish_mirror
    use dist_redistribute, only: finish_redistribute
    use neoclassical_terms, only: finish_neoclassical_terms
    use dissipation, only: finish_dissipation

    implicit none

    if (full_flux_surface) call finish_transforms
    call finish_redistribute
    call finish_dissipation
    call finish_parallel_nonlinearity
    call finish_wstar
    call finish_wdrift
    call finish_parallel_streaming
    call finish_mirror
    call finish_neoclassical_terms
    call deallocate_arrays

    time_advance_initialized = .false.
    readinit = .false.

  end subroutine finish_time_advance

  subroutine finish_parallel_nonlinearity

    implicit none

    if (allocated(par_nl_fac)) deallocate (par_nl_fac)
    if (allocated(par_nl_curv)) deallocate (par_nl_curv)
    if (allocated(par_nl_driftx)) deallocate (par_nl_driftx)
    if (allocated(par_nl_drifty)) deallocate (par_nl_drifty)

    parnlinit = .false.

  end subroutine finish_parallel_nonlinearity

  subroutine finish_wdrift

    use dist_fn_arrays, only: wdriftx_g, wdrifty_g
    use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi

    implicit none

    if (allocated(wdriftx_g)) deallocate (wdriftx_g)
    if (allocated(wdrifty_g)) deallocate (wdrifty_g)
    if (allocated(wdriftx_phi)) deallocate (wdriftx_phi)
    if (allocated(wdrifty_phi)) deallocate (wdrifty_phi)

    wdriftinit = .false.

  end subroutine finish_wdrift

  subroutine finish_wstar

    use dist_fn_arrays, only: wstar

    implicit none

    if (allocated(wstar)) deallocate (wstar)

    wstarinit = .false.

  end subroutine finish_wstar

  subroutine deallocate_arrays

    use dist_fn_arrays, only: g1, g2, g3

    implicit none

    if (allocated(g1)) deallocate (g1)
    if (allocated(g2)) deallocate (g2)
    if (allocated(g3)) deallocate (g3)

  end subroutine deallocate_arrays

end module time_advance
