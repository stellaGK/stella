
module time_advance

   use debug_flags, only: debug => time_advance_debug

   implicit none 
   
   public :: init_time_advance, finish_time_advance
   public :: advance_stella

   private

   logical :: time_advance_initialized = .false.
   ! logical :: wdriftinit = .false.
   ! logical :: wstarinit = .false.
   ! logical :: parnlinit = .false.
   logical :: readinit = .false.
   ! logical :: radialinit = .false.
   ! logical :: driftimpinit = .false.

contains

   !****************************************************************************
   !                          INITIALISE TIME ADVANCE                          !
   !****************************************************************************

   subroutine init_time_advance

      use mp, only: proc0
      use parameters_physics, only: radial_variation
      use parameters_physics, only: include_parallel_nonlinearity
      use neoclassical_terms, only: init_neoclassical_terms
      use dissipation, only: init_collisions, include_collisions
      use parallel_streaming, only: init_parallel_streaming
      use mirror_terms, only: init_mirror
      use flow_shear, only: init_flow_shear
      use sources, only: init_quasineutrality_source, init_source_timeaverage

      use arrays_drifts, only: init_wdrift, init_wstar
      use parallel_nonlinearity, only: init_parallel_nonlinearity
      use store_arrays_useful, only: wdriftinit, wstarinit, parnlinit, &
            radialinit, driftimpinit

      implicit none

      if (time_advance_initialized) return
      time_advance_initialized = .true.

      debug = debug .and. proc0

      !> read time_advance_knobs namelist from the input file;
      !> sets the explicit time advance option, as well as allows for scaling of
      !> the x and y components of the magnetic drifts and of the drive term
      !> allocate distribution function sized arrays needed, e.g., for Runge-Kutta time advance
      if (debug) write (6, *) 'time_advance::init_time_advance::allocate_arrays'
      call allocate_arrays
      !> set up neoclassical corrections to the equilibrium Maxwellian;
      !> only calculated/needed when simulating higher order terms in rhostar for intrinsic rotation
      if (debug) write (6, *) 'time_advance::init_time_advance::init_neoclassical_terms'
      call init_neoclassical_terms
      !> calculate the term multiplying dg/dvpa in the mirror term
      !> and set up either the semi-Lagrange machinery or the tridiagonal matrix to be inverted
      !> if solving implicitly
      if (debug) write (6, *) 'time_advance::init_time_advance::init_mirror'
      call init_mirror
      !> calculate the term multiplying dg/dz in the parallel streaming term
      !> and set up the tridiagonal matrix to be inverted if solving implicitly
      if (debug) write (6, *) 'time_advance::init_time_advance::init_parstream'
      call init_parallel_streaming
      !> allocate and calculate the factors multiplying dg/dx, dg/dy, dphi/dx and dphi/dy
      !> in the magnetic drift terms
      if (debug) write (6, *) 'time_advance::init_time_advance::init_wdrift'
      call init_wdrift 
      !> allocate and calculate the factor multiplying dphi/dy in the gradient drive term
      if (debug) write (6, *) 'time_advance::init_time_advance::init_wstar'
      call init_wstar 
      if (debug) write (6, *) 'time_advance::init_time_advance::init_flow_shear'
      call init_flow_shear
      if (debug) write (6, *) 'time_advance::init_time_advance::init_parallel_nonlinearity'
      if (include_parallel_nonlinearity) call init_parallel_nonlinearity 
      if (debug) write (6, *) 'time_advance::init_time_advance::init_radial_variation'
      if (radial_variation) call init_radial_variation
      if (include_collisions) then
         if (debug) write (6, *) 'time_advance::init_time_advance::init_collisions'
         call init_collisions
      end if
      if (debug) write (6, *) 'time_advance::init_time_advance::init_cfl'
      call init_cfl

      if (debug) write (6, *) 'time_advance::init_time_advance::init_source_timeaverage'
      call init_source_timeaverage
      if (debug) write (6, *) 'time_advance::init_time_advance::init_quasineutrality_source'
      call init_quasineutrality_source

   end subroutine init_time_advance


   subroutine init_radial_variation
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use stella_time, only: code_dt
      use species, only: spec, pfac
      use z_grid, only: nzgrid
      use parameters_kxky_grid, only: nalpha
      use geometry, only: drhodpsi, dydalpha, gfac
      use geometry, only: dBdrho, geo_surf, q_as_x
      use geometry, only: dcvdriftdrho, dcvdrift0drho
      use geometry, only: dgbdriftdrho, dgbdrift0drho
      use velocity_grids, only: vperp2, vpa, mu
      use velocity_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use store_arrays_distribution_fn, only: wstarp
      use store_arrays_distribution_fn, only: wdriftx_phi, wdrifty_phi
      use store_arrays_distribution_fn, only: wdriftpx_g, wdriftpy_g
      use store_arrays_distribution_fn, only: wdriftpx_phi, wdriftpy_phi!, adiabatic_phi

      use parameters_physics, only: xdriftknob, ydriftknob, wstarknob

      use store_arrays_useful, only: radialinit
      
!   use neoclassical_terms, only: include_neoclassical_terms

      implicit none

      integer :: is, imu, iv, ivmu
      real :: fac
      real, dimension(:, :), allocatable :: energy

      real, dimension(:, :), allocatable :: wcvdrifty, wgbdrifty
      real, dimension(:, :), allocatable :: wcvdriftx, wgbdriftx

!wstar

      if (radialinit) return
      radialinit = .true.

      allocate (wcvdrifty(nalpha, -nzgrid:nzgrid))
      allocate (wgbdrifty(nalpha, -nzgrid:nzgrid))
      allocate (wcvdriftx(nalpha, -nzgrid:nzgrid))
      allocate (wgbdriftx(nalpha, -nzgrid:nzgrid))

      if (.not. allocated(wstarp)) &
         allocate (wstarp(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstarp = 0.0
      if (.not. allocated(wdriftpx_phi)) &
         allocate (wdriftpx_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
!   if (.not.allocated(adiabatic_phi)) &
!      allocate (adiabatic_phi(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not. allocated(wdriftpy_phi)) &
         allocate (wdriftpy_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not. allocated(wdriftpx_g)) &
         allocate (wdriftpx_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not. allocated(wdriftpy_g)) &
         allocate (wdriftpy_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      allocate (energy(nalpha, -nzgrid:nzgrid))

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         is = is_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         iv = iv_idx(vmu_lo, ivmu)
         energy = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
         !FLAG DSO - THIS NEEDS TO BE ADDED SOMEDAY!
         !if (include_neoclassical_terms) then
         !   wstarp(:,:,ivmu) = dydalpha*drhodpsi*wstarknob*0.5*code_dt &
         !        * (maxwell_vpa(iv)*maxwell_mu(:,:,imu) &
         !        * (spec(is)%fprim+spec(is)%tprim*(energy-1.5)) &
         !        - dfneo_drho(:,:,ivmu))
         !else
         !recall that fprim = -dn/dr and trpim = -dt/dr

         wstarp(:, :, ivmu) = -wstarknob * 0.5 * code_dt &
                              * dydalpha * drhodpsi * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is) &
                              * (pfac * (spec(is)%d2ndr2 - (spec(is)%fprim)**2 - (spec(is)%tprim)**2 * energy) &
                                 + pfac * (spec(is)%d2Tdr2 - (spec(is)%tprim)**2) * (energy - 1.5) &
                                 - gfac * 2 * spec(is)%tprim * mu(imu) * spread(dBdrho, 1, nalpha) &
                                 + (spec(is)%fprim + spec(is)%tprim * (energy - 1.5)) &
                                 * (pfac * (spec(is)%fprim + spec(is)%tprim * (energy - 1.5)) &
                                    + gfac * 2 * mu(imu) * spread(dBdrho, 1, nalpha) &
                                    + gfac * drhodpsi * geo_surf%d2psidr2))

         !end if

         !wdrift
         fac = -ydriftknob * 0.5 * code_dt * spec(is)%tz_psi0
         ! this is the curvature drift piece of wdrifty with missing factor of vpa
         ! vpa factor is missing to avoid singularity when including
         ! non-Maxwellian corrections to equilibrium
         wcvdrifty = gfac * fac * dcvdriftdrho * vpa(iv)
         ! this is the grad-B drift piece of wdrifty
         wgbdrifty = gfac * fac * dgbdriftdrho * 0.5 * vperp2(:, :, imu)
         wdriftpy_g(:, :, ivmu) = wcvdrifty * vpa(iv) + wgbdrifty

         wdriftpy_phi(:, :, ivmu) = spec(is)%zt * (wgbdrifty + wcvdrifty * vpa(iv)) &
                                    * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is) &
                                    - wdrifty_phi(:, :, ivmu) * (pfac * (spec(is)%fprim + spec(is)%tprim * (energy - 2.5)) &
                                                                 + gfac * 2.*mu(imu) * spread(dBdrho, 1, nalpha))

         if (q_as_x) then
            fac = -xdriftknob * 0.5 * code_dt * spec(is)%tz_psi0
         else
            fac = -xdriftknob * 0.5 * code_dt * spec(is)%tz_psi0 / geo_surf%shat
         end if
         ! this is the curvature drift piece of wdriftx with missing factor of vpa
         ! vpa factor is missing to avoid singularity when including
         ! non-Maxwellian corrections to equilibrium
         wcvdriftx = gfac * fac * dcvdrift0drho * vpa(iv)
         ! this is the grad-B drift piece of wdriftx
         wgbdriftx = gfac * fac * dgbdrift0drho * 0.5 * vperp2(:, :, imu)
         wdriftpx_g(:, :, ivmu) = wgbdriftx + wcvdriftx * vpa(iv)

         wdriftpx_phi(:, :, ivmu) = spec(is)%zt * (wgbdriftx + wcvdriftx * vpa(iv)) &
                                    * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is) &
                                    - wdriftx_phi(:, :, ivmu) * (pfac * (spec(is)%fprim + spec(is)%tprim * (energy - 2.5)) &
                                                                 + gfac * 2.*mu(imu) * spread(dBdrho, 1, nalpha))

!      !the next piece is everything under the x derivative, as this needs to be
!      !transformed separately
!      wdriftpx_phi(:,:,ivmu) = spec(is)%zt*(wgbdriftx + wcvdriftx*vpa(iv))  &
!           * maxwell_vpa(iv,is)*maxwell_mu(:,:,imu,is)*maxwell_fac(is)

!      !this is variation in the Maxwellian part of the adiabatic response of phi,
!      !which needs to be transformed separately before differentiation wrt x
!      !the gyroaveraging and quasineutrality is already done in fields
!      adiabatic_phi(:,:,ivmu) = -(pfac*(spec(is)%fprim+spec(is)%tprim*(energy-2.5)) &
!                                 +gfac*2.*mu(imu)*spread(dBdrho,1,nalpha))

      end do

      deallocate (energy, wcvdriftx, wgbdriftx, wcvdrifty, wgbdrifty)

   end subroutine init_radial_variation

   subroutine allocate_arrays

      use stella_layouts, only: vmu_lo
      use z_grid, only: nzgrid, ntubes
      use parameters_kxky_grid, only: naky, nakx
      use store_arrays_distribution_fn, only: g0, g1, g2, g3
      use parameters_numerical, only: explicit_algorithm_switch, explicit_algorithm_rk3, &
           explicit_algorithm_rk2, explicit_algorithm_rk4, explicit_algorithm_euler

      implicit none

      if (.not. allocated(g0)) &
         allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g0 = 0.
      if (.not. allocated(g1)) &
         allocate (g1(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g1 = 0.
      if (.not. allocated(g2)) &
         allocate (g2(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g2 = 0.
      if (.not. allocated(g3) .and. explicit_algorithm_switch == explicit_algorithm_rk4) then
         allocate (g3(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         g3 = 0.
      else
         allocate (g3(1, 1, 1, 1, 1))
      end if

   end subroutine allocate_arrays

   subroutine init_cfl

      use mp, only: proc0, nproc, max_allreduce, min_allreduce
      use mp, only: scope, allprocs, subprocs
      use store_arrays_distribution_fn, only: wdriftx_g, wdrifty_g
      use stella_time, only: code_dt, write_dt, cfl_dt_linear
      use parameters_numerical, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
      use parameters_physics, only: radial_variation, prp_shear_enabled
      use z_grid, only: delzed
      use velocity_grids, only: dvpa
      use grids_kxky, only: akx, aky, rho
      use parameters_kxky_grid, only: nx
      use parameters_numerical, only: stream_implicit, mirror_implicit, drifts_implicit
      use parallel_streaming, only: stream
      use parallel_streaming, only: stream_rad_var1, stream_rad_var2
      use mirror_terms, only: mirror
      use flow_shear, only: prl_shear, shift_times
      use file_utils, only: runtype_option_switch, runtype_multibox
      use dissipation, only: include_collisions, collisions_implicit
      use dissipation, only: cfl_dt_vpadiff, cfl_dt_mudiff
      use debug_flags, only: print_extra_info_to_terminal

      use timing_of_run, only: reset_dt

      implicit none

      real :: cfl_dt_mirror, cfl_dt_stream, cfl_dt_shear
      real :: cfl_dt_wdriftx, cfl_dt_wdrifty
      real :: zero
      real :: wdriftx_max, wdrifty_max

      ! avoid divide by zero in cfl_dt terms below
      zero = 100.*epsilon(0.)

      ! FLAG -- assuming equal spacing in zed!

      if (cfl_dt_linear < 0) cfl_dt_linear = code_dt / cfl_cushion_upper

      if (.not. drifts_implicit) then
         ! get the local max value of wdriftx on each processor
         wdriftx_max = maxval(abs(wdriftx_g))
         ! compare these max values across processors to get global max
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

      !> TODO:GA- add correct CFL condition 
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
         !while other quantities should go here, parallel streaming with electrons
         !is what will limit us
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
         ! get the local max value of wdrifty on each processor
         wdrifty_max = maxval(abs(wdrifty_g))
         ! compare these max values across processors to get global max
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

   ! subroutine reset_dt

   !    use parallel_streaming, only: parallel_streaming_initialized
   !    use parallel_streaming, only: init_parallel_streaming
   !    use dissipation, only: init_collisions, collisions_initialized, include_collisions
   !    use parameters_numerical, only: stream_implicit, driftkinetic_implicit
   !    use response_matrix, only: response_matrix_initialized
   !    use response_matrix, only: init_response_matrix
   !    use mirror_terms, only: mirror_initialized
   !    use mirror_terms, only: init_mirror
   !    use flow_shear, only: flow_shear_initialized
   !    use flow_shear, only: init_flow_shear
   !    use parameters_physics, only: radial_variation
   !    use sources, only: init_source_timeaverage
   !    use sources, only: init_quasineutrality_source, qn_source_initialized
   !    use arrays_drifts, only: init_wdrift, init_wstar
   !    use store_arrays_useful, only: wdriftinit, wstarinit, radialinit, driftimpinit, &
   !                      flow_shear_initialized, mirror_initialized, &
   !                      parallel_streaming_initialized, qn_source_initialized
   !    implicit none

   !    ! need to recompute mirror and streaming terms
   !    ! to account for updated code_dt
   !    wdriftinit = .false.
   !    wstarinit = .false.
   !    radialinit = .false.
   !    driftimpinit = .false.
   !    flow_shear_initialized = .false.
   !    mirror_initialized = .false.
   !    parallel_streaming_initialized = .false.
   !    qn_source_initialized = .false.

   !    if (debug) write (6, *) 'time_advance::reset_dt::init_wstar'
   !    call init_wstar 
   !    if (debug) write (6, *) 'time_advance::reset_dt::init_wdrift'
   !    call init_wdrift 
   !    if (debug) write (6, *) 'time_advance::reset_dt::init_mirror'
   !    call init_mirror
   !    if (debug) write (6, *) 'time_advance::reset_dt::init_parallel_streaming'
   !    call init_parallel_streaming
   !    if (debug) write (6, *) 'time_advance::reset_dt::init_flow_shear'
   !    call init_flow_shear
   !    if (debug) write (6, *) 'time_advance::reset_dt::init_source_timeaverage'
   !    call init_source_timeaverage
   !    if (debug) write (6, *) 'time_advance::reset_dt::init_quasineutrality_source'
   !    call init_quasineutrality_source
   !    if (radial_variation) then
   !       if (debug) write (6, *) 'time_advance::reset_dt::init_radial_variation'
   !       call init_radial_variation
   !    end if
   !    if (include_collisions) then
   !       if (debug) write (6, *) 'time_advance::reset_dt::init_collisions'
   !       collisions_initialized = .false.
   !       call init_collisions
   !    end if
   !    ! do not try to re-init response matrix
   !    ! before it has been initialized the first time
   !    if ((stream_implicit .or. driftkinetic_implicit) .and. response_matrix_initialized) then
   !       response_matrix_initialized = .false.
   !       if (debug) write (6, *) 'time_advance::reset_dt::init_response_matrix'
   !       call init_response_matrix
   !    end if

   ! end subroutine reset_dt

   !****************************************************************************
   !****************************************************************************
   !                        MAIN TIME ADVANCE OF STELLA                        !
   !****************************************************************************
   !****************************************************************************

   subroutine advance_stella(istep, stop_stella)

      use store_arrays_distribution_fn, only: gold, gnew
      use store_arrays_fields, only: phi, apar, bpar
      use store_arrays_fields, only: phi_old, apar_old
      use fields, only: advance_fields, fields_updated
      use parameters_numerical, only: fully_explicit, fully_implicit
      use parameters_multibox, only: rk_step
      use sources, only: include_qn_source, update_quasineutrality_source
      use sources, only: source_option_switch, source_option_projection
      use sources, only: source_option_krook
      use sources, only: update_tcorr_krook, project_out_zero
      use parameters_physics, only: include_apar
      use mp, only: proc0, broadcast

      use parameters_numerical, only: flip_flop
      implicit none

      integer, intent(in) :: istep
      logical, intent(in out) :: stop_stella

      logical :: restart_time_step, time_advance_successful
      integer :: count_restarts

      !> unless running in multibox mode, no need to worry about
      !> mb_communicate calls as the subroutine is immediately exited
      !> if not in multibox mode.
      if (.not. rk_step) then
         if (debug) write (*, *) 'time_advance::multibox'
         call mb_communicate(gnew)
      end if

      !> save value of phi & apar
      !> for use in diagnostics (to obtain frequency)
      phi_old = phi
      if (include_apar) apar_old = apar

      ! Flag which is set to true once we've taken a step without needing to
      ! reset dt (which can be done by the nonlinear term(s))
      time_advance_successful = .false.

      ! If cfl_cushion_lower is chosen too close to cfl_cushion_upper, then
      ! we might get stuck restarting the time step over and over, so exit stella
      count_restarts = 1

      ! Attempt the Lie or flip-flop time advance until we've done it without the
      ! timestep changing.
      do while (.not. time_advance_successful)

         ! If we've already attempted a time advance then we've updated gnew, so reset it.
         gnew = gold

         ! Ensure fields are consistent with gnew.
         call advance_fields(gnew, phi, apar, bpar, dist='g')

         ! Keep track whether any routine wants to modify the time step
         restart_time_step = .false.

         !> reverse the order of operations every time step
         !> as part of alternating direction operator splitting
         !> this is needed to ensure 2nd order accuracy in time
         if (mod(istep, 2) == 1 .or. .not. flip_flop) then

            !> Advance the explicit parts of the GKE
            if (debug) write (*, *) 'time_advance::advance_explicit'
            if (.not. fully_implicit) call advance_explicit(gnew, restart_time_step, istep)
            if (debug) write (*, *) 'time_advance::advance_implicit'
            !> Use operator splitting to separately evolve all terms treated implicitly
            if (.not. restart_time_step .and. .not. fully_explicit) call advance_implicit(istep, phi, apar, bpar, gnew)
         else
            if (debug) write (*, *) 'time_advance::advance_implicit'
            if (.not. fully_explicit) call advance_implicit(istep, phi, apar, bpar, gnew)
            if (debug) write (*, *) 'time_advance::advance_explicit'
            if (.not. fully_implicit) call advance_explicit(gnew, restart_time_step, istep)
         end if

         ! If the time step has not been restarted, the time advance was succesfull
         ! Otherwise, discard changes to gnew and start the time step again, fields
         ! will have to be recalculated
         if (.not. restart_time_step) then
            time_advance_successful = .true.
         else
            count_restarts = count_restarts + 1
            fields_updated = .false.
         end if

         ! At some point, give up on restarting the time step
         if (count_restarts > 5) then
            stop_stella = .true.
            call broadcast(stop_stella)
            gnew = gold
            fields_updated = .false.
            if (proc0) then
               write (*, *)
               write (*, *) 'EXITING STELLA BECAUSE WE ALREADY RESTARTED THE TIME STEP 5 TIMES.'
               write (*, *) 'CHANGE CFL_CUSHION_UPPER AND CFL_CUSHION_LOWER AND RESTART THE SIMULATION.'
               write (*, *)
            end if
            exit
         end if

      end do

      ! presumably this is to do with the radially global version of the code?
      ! perhaps it could be packaged together with thee update_delay_krook code
      ! below and made into a single call where all of this happens so that
      ! users of the flux tube version of the code need not worry about it.
      if (source_option_switch == source_option_projection) then
         call project_out_zero(gold, gnew)
         fields_updated = .false.
      end if

      gold = gnew

      !> Ensure fields are updated so that omega calculation is correct.
      call advance_fields(gnew, phi, apar, bpar, dist='g')

      !update the delay parameters for the Krook operator
      if (source_option_switch == source_option_krook) call update_tcorr_krook(gnew)
      if (include_qn_source) call update_quasineutrality_source

   end subroutine advance_stella

   !-------------------------------------------------------------------------------
   !                       EXPLICIT TIME ADVANCE SUBROUTINES
   !-------------------------------------------------------------------------------
   !> advance_explicit takes as input the guiding centre distribution function
   !> in k-space and updates it to account for all of the terms in the GKE that
   !> are advanced explicitly in time
   subroutine advance_explicit(g, restart_time_step, istep)

      use mp, only: proc0
      use job_manage, only: time_message
      use z_grid, only: nzgrid
      use extended_zgrid, only: periodic, phase_shift
      use parameters_kxky_grid, only: naky
      use stella_layouts, only: vmu_lo, iv_idx
      use parameters_physics, only: include_apar
      use parallel_streaming, only: stream_sign
      use store_arrays_fields, only: phi, apar, bpar
      use fields, only: advance_fields
      use g_tofrom_h, only: gbar_to_g

      use parameters_numerical, only: explicit_algorithm_switch, explicit_algorithm_rk3, &
           explicit_algorithm_rk2, explicit_algorithm_rk4, explicit_algorithm_euler

      use store_arrays_useful, only: time_gke
      
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      integer :: ivmu, iv, sgn, iky

      !> start the timer for the explicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 8), ' explicit')

      ! incoming pdf is g = h - (Z F0/T) (J0 phi + 4 mu (T/Z) (J1/gamma) bpar)
      ! if include_apar = T, convert from g to gbar = g + Z F0/T (2J0 vpa vth apar),
      ! as gbar appears in time derivative
      if (include_apar) then
         ! if the fields are not already updated, then update them
         call advance_fields(g, phi, apar, bpar, dist='g')
         call gbar_to_g(g, apar, -1.0)
      end if

      select case (explicit_algorithm_switch)
      case (explicit_algorithm_euler)
         !> forward Euler
         call advance_explicit_euler(g, restart_time_step, istep)
      case (explicit_algorithm_rk2)
         !> SSP RK2
         call advance_explicit_rk2(g, restart_time_step, istep)
      case (explicit_algorithm_rk3)
         !> default is SSP RK3
         call advance_explicit_rk3(g, restart_time_step, istep)
      case (explicit_algorithm_rk4)
         !> RK4
         call advance_explicit_rk4(g, restart_time_step, istep)
      end select

      if (include_apar) then
         ! if the fields are not already updated, then update them
         call advance_fields(g, phi, apar, bpar, dist='gbar')
         ! implicit solve will use g rather than gbar for advance,
         ! so convert from gbar to g
         call gbar_to_g(g, apar, 1.0)
      end if

      !> enforce periodicity for periodic (including zonal) modes
      do iky = 1, naky
         if (periodic(iky)) then
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               !> stream_sign > 0 corresponds to dz/dt < 0
               sgn = stream_sign(iv)
               g(iky, :, sgn * nzgrid, :, ivmu) = &
                  g(iky, :, -sgn * nzgrid, :, ivmu) * phase_shift(iky)**(-sgn)
            end do
         end if
      end do

      !> stop the timer for the explicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 8), ' explicit')

   end subroutine advance_explicit

   !-------------------------------------------------------------------------------
   !                        EXPLICIT EULER TIME ADVANCE SUBROUTINE
   !-------------------------------------------------------------------------------
   !> advance_explicit_euler uses forward Euler to advance one time step
   subroutine advance_explicit_euler(g, restart_time_step, istep)

      use store_arrays_distribution_fn, only: g0
      use z_grid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      !> rk_step only true if running in multibox mode
      if (rk_step) call mb_communicate(g)

      g0 = g

      call solve_gke(g0, g, restart_time_step, istep)

      g = g0 + g

   end subroutine advance_explicit_euler

   !-------------------------------------------------------------------------------
   !                         EXPLICIT RK2 TIME ADVANCE SUBROUTINE
   !-------------------------------------------------------------------------------
   !> advance_expliciit_rk2 uses strong stability-preserving rk2 to advance one time step
   subroutine advance_explicit_rk2(g, restart_time_step, istep)

      use store_arrays_distribution_fn, only: g0, g1
      use z_grid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      integer :: icnt

      !> rk_step only true if running in multibox mode
      if (rk_step) call mb_communicate(g)

      g0 = g
      icnt = 1

      !> SSP rk2 algorithm to advance explicit part of code
      !> if GK equation written as dg/dt = rhs - vpar . grad h,
      !> solve_gke returns rhs*dt
      do while (icnt <= 2)
         select case (icnt)
         case (1)
            call solve_gke(g0, g1, restart_time_step, istep)
         case (2)
            g1 = g0 + g1
            if (rk_step) call mb_communicate(g1)
            call solve_gke(g1, g, restart_time_step, istep)
         end select
         if (restart_time_step) then
            ! If the code_dt is reset, we need to quit this loop and restart the timestep again
            icnt = 10
         else
            icnt = icnt + 1
         end if
      end do

      !> this is g at intermediate time level
      g = 0.5 * g0 + 0.5 * (g1 + g)

   end subroutine advance_explicit_rk2

   !-------------------------------------------------------------------------------
   !                         EXPLICIT RK3 TIME ADVANCE SUBROUTINE
   !-------------------------------------------------------------------------------
   !> strong stability-preserving RK3
   subroutine advance_explicit_rk3(g, restart_time_step, istep)

      use store_arrays_distribution_fn, only: g0, g1, g2
      use z_grid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      integer :: icnt

      !> rk_STEP = false unless in multibox mode
      if (rk_step) call mb_communicate(g)

      g0 = g
      icnt = 1

      !> SSP rk3 algorithm to advance explicit part of code
      !> if GK equation written as dg/dt = rhs - vpar . grad h,
      !> solve_gke returns rhs*dt
      do while (icnt <= 3)
         select case (icnt)
         case (1)
            call solve_gke(g0, g1, restart_time_step, istep)
         case (2)
            g1 = g0 + g1
            if (rk_step) call mb_communicate(g1)
            call solve_gke(g1, g2, restart_time_step, istep)
         case (3)
            g2 = g1 + g2
            if (rk_step) call mb_communicate(g2)
            call solve_gke(g2, g, restart_time_step, istep)
         end select
         if (restart_time_step) then
            ! If the code_dt is reset, we need to quit this loop and restart the timestep again
            icnt = 10
         else
            icnt = icnt + 1
         end if
      end do

      !> this is g at intermediate time level
      g = g0 / 3.+0.5 * g1 + (g2 + g) / 6.

   end subroutine advance_explicit_rk3

   !-------------------------------------------------------------------------------
   !                         EXPLICIT RK4 TIME ADVANCE SUBROUTINE
   !-------------------------------------------------------------------------------
   !> standard RK4
   subroutine advance_explicit_rk4(g, restart_time_step, istep)

      use store_arrays_distribution_fn, only: g0, g1, g2, g3
      use z_grid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      integer :: icnt

      !> rk_step is false unless in multibox mode
      if (rk_step) call mb_communicate(g)

      g0 = g
      icnt = 1

      !> RK4 algorithm to advance explicit part of code
      !> if GK equation written as dg/dt = rhs - vpar . grad h,
      !> solve_gke returns rhs*dt
      do while (icnt <= 4)
         select case (icnt)
         case (1)
            call solve_gke(g0, g1, restart_time_step, istep)
         case (2)
            ! g1 is h*k1
            g3 = g0 + 0.5 * g1
            if (rk_step) call mb_communicate(g3)
            call solve_gke(g3, g2, restart_time_step, istep)
            g1 = g1 + 2.*g2
         case (3)
            ! g2 is h*k2
            g2 = g0 + 0.5 * g2
            if (rk_step) call mb_communicate(g2)
            call solve_gke(g2, g3, restart_time_step, istep)
            g1 = g1 + 2.*g3
         case (4)
            ! g3 is h*k3
            g3 = g0 + g3
            if (rk_step) call mb_communicate(g3)
            call solve_gke(g3, g, restart_time_step, istep)
            g1 = g1 + g
         end select
         if (restart_time_step) then
            ! If the code_dt is reset, we need to quit this loop and restart the timestep again
            icnt = 10
         else
            icnt = icnt + 1
         end if
      end do

      !> this is g at intermediate time level
      g = g0 + g1 / 6.

   end subroutine advance_explicit_rk4

   !-------------------------------------------------------------------------------
   !                  NEEDED FOR ALL EXPLICIT TIME ADVANCE SUBROUTINES
   !-------------------------------------------------------------------------------
   !> solve_gke accepts as argument pdf, the guiding centre distribution function in k-space,
   !> and returns rhs_ky, the right-hand side of the gyrokinetic equation in k-space;
   !> i.e., if dg/dt = r, then rhs_ky = r*dt;
   !> note that if include_apar = T, then the input pdf is actually gbar = g + (Ze/T)*(vpa/c)*<Apar>*F0
   subroutine solve_gke(pdf, rhs_ky, restart_time_step, istep)

      use job_manage, only: time_message
      use store_arrays_fields, only: phi, apar, bpar
      use stella_layouts, only: vmu_lo
      use stella_transforms, only: transform_y2ky
      use parameters_physics, only: include_parallel_nonlinearity
      use parameters_physics, only: include_parallel_streaming
      use parameters_physics, only: include_mirror, include_apar
      use parameters_physics, only: include_nonlinear, include_bpar
      use parameters_physics, only: full_flux_surface, radial_variation
      use parameters_physics, only: g_exb
      use z_grid, only: nzgrid, ntubes
      use parameters_kxky_grid, only: ikx_max, ny, naky_all
      use calculations_kxky, only: swap_kxky_back
      use grids_kxky, only: zonal_mode, akx
      use parameters_numerical, only: stream_implicit, mirror_implicit, drifts_implicit
      use dissipation, only: include_collisions, advance_collisions_explicit, collisions_implicit
      use sources, only: source_option_switch, source_option_krook
      use sources, only: add_krook_operator
      use parallel_streaming, only: advance_parallel_streaming_explicit
      use fields, only: fields_updated, advance_fields
      use fields_radial_variation, only: get_radial_correction
      use mirror_terms, only: advance_mirror_explicit
      use flow_shear, only: advance_parallel_flow_shear
      use parameters_multibox, only: include_multibox_krook
      use multibox, only: add_multibox_krook
      use store_arrays_distribution_fn, only: g_scratch
      use gyro_averages, only: gyro_average
      use arrays_gyro_averages, only: j0_ffs
      use g_tofrom_h, only: gbar_to_g 
      use dissipation, only: hyper_dissipation
      ! TMP FOR TESTING -- MAB
      use fields, only: fields_updated
      use parameters_physics, only: xdriftknob, ydriftknob

      use advance_explicit_drifts, only: advance_wstar_explicit
      use advance_explicit_drifts, only: advance_wdriftx_explicit, advance_wdrifty_explicit
      use parallel_nonlinearity, only: advance_parallel_nonlinearity

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: pdf
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out), target :: rhs_ky
      logical, intent(out) :: restart_time_step
      integer, intent(in) :: istep

      complex, dimension(:, :, :, :, :), allocatable, target :: rhs_y
      complex, dimension(:, :, :, :, :), pointer :: rhs
      complex, dimension(:, :), allocatable :: rhs_ky_swap

      integer :: iz, it, ivmu

      rhs_ky = 0.

      !> if full_flux_surface = .true., then initially obtain the RHS of the GKE in alpha-space;
      !> will later inverse Fourier transform to get RHS in k_alpha-space
      if (full_flux_surface) then
         !> rhs_ky will always be needed as the array returned by the subroutine,
         !> but intermediate array rhs_y (RHS of gke in alpha-space) only needed for full_flux_surface = .true.
         allocate (rhs_y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         rhs_y = 0.
         !> rhs is array referred to for both flux tube and full-flux-surface simulations;
         !> for full-flux-surface it should point to rhs_y
         rhs => rhs_y
      else
         !> rhs is array referred to for both flux tube and full-flux-surface simulations;
         !> for flux tube it should point to rhs_ky
         rhs => rhs_ky
      end if

      !> start with g in k-space and (ky,kx,z) local
      !> obtain fields corresponding to g
      if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_fields'

      ! if advancing apar, then gbar is evolved in time rather than g
      if (include_apar) then
         call advance_fields(pdf, phi, apar, bpar, dist='gbar')

         ! convert from gbar to g = h - (Z F0/T)( J0 phi + 4 mu (T/Z) (J1/gamma) bpar),
         ! as all terms on RHS of GKE use g rather than gbar
         call gbar_to_g(pdf, apar, 1.0)
      else
         call advance_fields(pdf, phi, apar, bpar, dist='g')
      end if

      if (radial_variation) call get_radial_correction(pdf, phi, dist='gbar')

      !> obtain the gyro-average of the electrostatic potential phi and store in g_scratch;
      !> this can be a particularly costly operation when simulating a full flux surface
      !> due to the coupling of different k-alphas inherent in the gyro-average;
      !> calculate once here to avoid repeated calculation later
      !> TODO-GA : can this be spec up??
      if (full_flux_surface) call gyro_average(phi, g_scratch, j0_ffs)

      !! INSERT TEST HERE TO SEE IF dg/dy, dg/dx, d<phi>/dy, d<phi>/dx WILL BE NEEDED
      !! IF SO, PRE-COMPUTE ONCE HERE

      !> default is to continue with same time step size.
      !> if estimated CFL condition for nonlinear terms is violated
      !> then restart_time_step will be set to .true.
      restart_time_step = .false.
      !> calculate and add ExB nonlinearity to RHS of GK eqn
      !> do this first, as the CFL condition may require a change in time step
      !> and thus recomputation of mirror, wdrift, wstar, and parstream
      if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_ExB_nonlinearity'
      if (include_nonlinear) call advance_ExB_nonlinearity(pdf, rhs, restart_time_step, istep)

      !> include contribution from the parallel nonlinearity (aka turbulent acceleration)
      if (include_parallel_nonlinearity .and. .not. restart_time_step) &
         call advance_parallel_nonlinearity(pdf, rhs, restart_time_step)

      if (.not. restart_time_step) then

         !> include contribution from perp flow shear in the parallel component of the toroidal flow
         if ((g_exb**2) > epsilon(0.0)) call advance_parallel_flow_shear(rhs)

         !> calculate and add mirror term to RHS of GK eqn
         if (include_mirror .and. .not. mirror_implicit) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_mirror_explicit'
            call advance_mirror_explicit(pdf, rhs)
         end if

         if (.not. drifts_implicit) then
            !> calculate and add alpha-component of magnetic drift term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wdrifty_explicit'
            if (abs(ydriftknob) > epsilon(0.0)) then
               call advance_wdrifty_explicit(pdf, phi, bpar, rhs)
            end if

            !> calculate and add psi-component of magnetic drift term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wdriftx_explicit'
            if (abs(xdriftknob) > epsilon(0.0)) then
               call advance_wdriftx_explicit(pdf, phi, bpar, rhs)
            end if
            
            !> calculate and add omega_* term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wstar_explicit'
            call advance_wstar_explicit(phi, rhs)
         end if

         !> calculate and add contribution from collisions to RHS of GK eqn
         if (include_collisions .and. .not. collisions_implicit) call advance_collisions_explicit(pdf, phi, bpar, rhs)

         !> calculate and add parallel streaming term to RHS of GK eqn
         if (include_parallel_streaming .and. (.not. stream_implicit)) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_parallel_streaming_explicit'
            call advance_parallel_streaming_explicit(pdf, phi, bpar, rhs)
         end if
         
         if (hyper_dissipation) then
            call advance_hyper_explicit(pdf, rhs)
         end if

         !> if simulating a full flux surface (flux annulus), all terms to this point have been calculated
         !> in real-space in alpha (y); transform to kalpha (ky) space before adding to RHS of GKE.
         !> NB: it may be that for fully explicit calculation, this transform can be eliminated with additional code changes
         if (full_flux_surface) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::transform_y2ky'
            allocate (rhs_ky_swap(naky_all, ikx_max))
            it = 1
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               do iz = -nzgrid, nzgrid
                  call transform_y2ky(rhs_y(:, :, iz, it, ivmu), rhs_ky_swap)
                  call swap_kxky_back(rhs_ky_swap, rhs_ky(:, :, iz, it, ivmu))
               end do
               ! ensure that the kx=ky=0 mode is zeroed out
               if (zonal_mode(1) .and. akx(1) < epsilon(0.)) then
                  rhs_ky(1, 1, :, it, ivmu) = 0.0
               end if
            end do
            deallocate (rhs_ky_swap)
         end if

         if (radial_variation) call advance_radial_variation(pdf, rhs)

         if (source_option_switch == source_option_krook) call add_krook_operator(pdf, rhs)

         if (include_multibox_krook) call add_multibox_krook(pdf, rhs)

      end if

      ! if advancing apar, need to convert input pdf back from g to gbar
      if (include_apar) call gbar_to_g(pdf, apar, -1.0)

      fields_updated = .false.

      if (allocated(rhs_y)) deallocate (rhs_y)
      nullify (rhs)

   end subroutine solve_gke

   subroutine advance_hyper_explicit(gin, gout)

      use stella_layouts, only: vmu_lo
      use z_grid, only: nzgrid, ntubes
      use parameters_kxky_grid, only: naky, nakx
      use hyper, only: advance_hyper_vpa, advance_hyper_zed
      use hyper, only: hyp_zed, hyp_vpa

      implicit none
      
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      complex, dimension(:, :, :, :, :), allocatable :: dg

      allocate (dg(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); dg = 0.0

      if (hyp_zed) then
         call advance_hyper_zed(gin, dg)
         gout = gout + dg
      end if
      if (hyp_vpa) then
         call advance_hyper_vpa(gin, dg)
      end if
      deallocate (dg)
      
    end subroutine advance_hyper_explicit

   subroutine advance_ExB_nonlinearity(g, gout, restart_time_step, istep)

      use mp, only: proc0, min_allreduce
      use mp, only: scope, allprocs, subprocs
      use stella_layouts, only: vmu_lo, imu_idx, is_idx
      use job_manage, only: time_message
      use gyro_averages, only: gyro_average
      use calculations_kxky_derivatives, only: get_dchidx, get_dchidy
      use store_arrays_fields, only: phi, apar, bpar, shift_state
      use store_arrays_fields, only: phi_corr_QN, phi_corr_GA

      use stella_transforms, only: transform_y2ky, transform_x2kx
      use stella_transforms, only: transform_y2ky_xfirst, transform_x2kx_xfirst
      use stella_time, only: cfl_dt_ExB, cfl_dt_linear, code_dt, code_dt_max
      use parameters_numerical, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
      use parameters_physics, only: g_exb, g_exbfac, fphi
      use z_grid, only: nzgrid, ntubes
      use geometry, only: exb_nonlin_fac, exb_nonlin_fac_p, gfac
      use parameters_kxky_grid, only: nakx, ikx_max, naky, naky_all, nx, ny
      use grids_kxky, only: akx, aky, rho_clamped
      use parameters_physics, only: full_flux_surface, radial_variation
      use parameters_physics, only: prp_shear_enabled, hammett_flow_shear
      use parameters_physics, only: include_apar, include_bpar
      use grids_kxky, only: x
      use calculations_kxky, only: swap_kxky, swap_kxky_back
      use constants, only: pi, zi
      use file_utils, only: runtype_option_switch, runtype_multibox
      use parameters_physics, only: suppress_zonal_interaction
      use store_arrays_distribution_fn, only: g_scratch
      use g_tofrom_h, only: g_to_h

      use calculations_kxky_derivatives, only: get_dgdy, get_dgdx

      use store_arrays_useful, only: time_gke
      use timing_of_run, only: reset_dt

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
      logical, intent(out) :: restart_time_step
      integer, intent(in) :: istep

      complex, dimension(:, :), allocatable :: g0k, g0a, g0k_swap
      complex, dimension(:, :), allocatable :: g0kxy, g0xky, prefac
      real, dimension(:, :), allocatable :: g0xy, g1xy, bracket

      real :: zero, cfl_dt
      integer :: ivmu, iz, it, imu, is
      logical :: yfirst

      ! alpha-component of magnetic drift (requires ky -> y)
      if (proc0) call time_message(.false., time_gke(:, 7), ' ExB nonlinear advance')

      if (debug) write (*, *) 'time_advance::solve_gke::advance_ExB_nonlinearity::get_dgdy'

      ! avoid divide by zero in cfl_dt terms below
      zero = 100.*epsilon(0.)

      ! Initialize cfl_dt_ExB
      cfl_dt_ExB = 10000000.

      restart_time_step = .false.
      ! this statement seems to imply that flow shear is not compatible with FFS
      ! need to check
      yfirst = .not. prp_shear_enabled

      allocate (g0k(naky, nakx))
      allocate (g0a(naky, nakx))
      allocate (g0xy(ny, nx))
      allocate (g1xy(ny, nx))
      allocate (bracket(ny, nx))
      allocate (prefac(naky, nx))

      if (yfirst) then
         allocate (g0k_swap(naky_all, ikx_max))
         allocate (g0kxy(ny, ikx_max))
      else
         allocate (g0xky(naky, nx))
      end if

      !> compute phase factor needed when running with equilibrium flow shear
      prefac = 1.0
      if (prp_shear_enabled .and. hammett_flow_shear) then
         prefac = exp(-zi * g_exb * g_exbfac * spread(x, 1, naky) * spread(aky * shift_state, 2, nx))
      end if

      ! incoming pdf is g = <f>
      ! for EM simulations, the pdf entering the ExB nonlinearity needs to be
      ! the non-Boltzmann part of f (h = f + (Ze/T)*phi*F0)
      if (include_apar .or. include_bpar) call g_to_h(g, phi, bpar, fphi)

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               !> compute i*ky*g
               call get_dgdy(g(:, :, iz, it, ivmu), g0k)
               !> FFT to get dg/dy in (y,x) space
               call forward_transform(g0k, g0xy)
               !> compute i*kx*<chi>
               if (full_flux_surface) then
                  call get_dgdx(g_scratch(:, :, iz, it, ivmu), g0k)
               else
                  call get_dchidx(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k)
               end if
               !> zero out the zonal contribution to d<chi>/dx if requested
               if (suppress_zonal_interaction) then
                  g0k(1,:) = 0.0
               end if
               !> if running with equilibrium flow shear, make adjustment to
               !> the term multiplying dg/dy
               if (prp_shear_enabled .and. hammett_flow_shear) then
                  call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0a)
                  g0k = g0k - g_exb * g_exbfac * spread(shift_state, 2, nakx) * g0a
               end if
               !> FFT to get d<chi>/dx in (y,x) space
               call forward_transform(g0k, g1xy)
               !> multiply by the geometric factor appearing in the Poisson bracket;
               !> i.e., (dx/dpsi*dy/dalpha)*0.5
               g1xy = g1xy * exb_nonlin_fac
               !> compute the contribution to the Poisson bracket from dg/dy*d<chi>/dx
               bracket = g0xy * g1xy

               !> estimate the CFL dt due to the above contribution
               cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * aky(naky), zero))

               if (radial_variation) then
                  bracket = bracket + gfac * g0xy * g1xy * exb_nonlin_fac_p * spread(rho_clamped, 1, ny)
                  call gyro_average(phi_corr_QN(:, :, iz, it), iz, ivmu, g0a)
                  g0a = fphi * (g0a + phi_corr_GA(:, :, iz, it, ivmu))
                  call get_dgdx(g0a, g0k)
                  call forward_transform(g0k, g1xy)
                  g1xy = g1xy * exb_nonlin_fac
                  bracket = bracket + g0xy * g1xy
                  !> estimate the CFL dt due to the above contribution
                  cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * aky(naky), zero))
               end if

               !> compute dg/dx in k-space (= i*kx*g)
               call get_dgdx(g(:, :, iz, it, ivmu), g0k)
               !> zero out the zonal contribution to dg/dx if requested
               if (suppress_zonal_interaction) then
                  g0k(1,:) = 0.0
               end if
               !> if running with equilibrium flow shear, correct dg/dx term
               if (prp_shear_enabled .and. hammett_flow_shear) then
                  call get_dgdy(g(:, :, iz, it, ivmu), g0a)
                  g0k = g0k - g_exb * g_exbfac * spread(shift_state, 2, nakx) * g0a
               end if
               !> FFT to get dg/dx in (y,x) space
               call forward_transform(g0k, g0xy)
               !> compute d<chi>/dy in k-space
               if (full_flux_surface) then
                  call get_dgdy(g_scratch(:, :, iz, it, ivmu), g0k)
               else
                  call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k)
               end if
               !> FFT to get d<chi>/dy in (y,x) space
               call forward_transform(g0k, g1xy)
               !> multiply by the geometric factor appearing in the Poisson bracket;
               !> i.e., (dx/dpsi*dy/dalpha)*0.5
               g1xy = g1xy * exb_nonlin_fac
               !> compute the contribution to the Poisson bracket from dg/dy*d<chi>/dx
               bracket = bracket - g0xy * g1xy

               !> estimate the CFL dt due to the above contribution
               cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * akx(ikx_max), zero))

               if (radial_variation) then
                  bracket = bracket - gfac * g0xy * g1xy * exb_nonlin_fac_p * spread(rho_clamped, 1, ny)
                  call gyro_average(phi_corr_QN(:, :, iz, it), iz, ivmu, g0a)
                  g0a = fphi * (g0a + phi_corr_GA(:, :, iz, it, ivmu))
                  call get_dgdy(g0a, g0k)
                  call forward_transform(g0k, g1xy)
                  g1xy = g1xy * exb_nonlin_fac
                  bracket = bracket - g0xy * g1xy
                  !> estimate the CFL dt due to the above contribution
                  cfl_dt_ExB = min(cfl_dt_ExB, 2.*pi / max(maxval(abs(g1xy)) * akx(ikx_max), zero))
               end if

               if (yfirst) then
                  call transform_x2kx(bracket, g0kxy)
                  if (full_flux_surface) then
                     gout(:, :, iz, it, ivmu) = g0kxy
                  else
                     call transform_y2ky(g0kxy, g0k_swap)
                     call swap_kxky_back(g0k_swap, gout(:, :, iz, it, ivmu))
                  end if
               else
                  call transform_y2ky_xfirst(bracket, g0xky)
                  g0xky = g0xky / prefac
                  call transform_x2kx_xfirst(g0xky, gout(:, :, iz, it, ivmu))
               end if
            end do
         end do
      end do

      ! convert back from h to g = <f> (only needed for EM sims)
      if (include_apar .or. include_bpar) call g_to_h(g, phi, bpar, -fphi)

      deallocate (g0k, g0a, g0xy, g1xy, bracket)
      if (allocated(g0k_swap)) deallocate (g0k_swap)
      if (allocated(g0xky)) deallocate (g0xky)
      if (allocated(g0kxy)) deallocate (g0kxy)

      if (runtype_option_switch == runtype_multibox) call scope(allprocs)

      call min_allreduce(cfl_dt_ExB)

      if (runtype_option_switch == runtype_multibox) call scope(subprocs)

      !> check estimated cfl_dt to see if the time step size needs to be changed
      cfl_dt = min(cfl_dt_ExB, cfl_dt_linear)
      if (code_dt > cfl_dt * cfl_cushion_upper) then
         if (proc0) then
            write (*, *) ' '
            write (*, '(A30,I0,A1)') 'CHANGING TIME STEP: (istep = ', istep, ')'
            write (*, '(A22, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
            write (*, '(A22, ES10.2E2)') "   cfl_dt_ExB:"//REPEAT(' ', 50), cfl_dt_ExB
            write (*, '(A22, ES10.2E2)') "   cfl_dt_linear:"//REPEAT(' ', 50), cfl_dt_linear
            write (*, '(A22, ES10.2E2)') "   cfl_cushion_upper:"//REPEAT(' ', 50), cfl_cushion_upper
            write (*, '(A22, ES10.2E2)') "   cfl_cushion_middle:"//REPEAT(' ', 50), cfl_cushion_middle
            write (*, '(A22, ES10.2E2)') "   cfl_cushion_lower:"//REPEAT(' ', 50), cfl_cushion_lower
            write (*, '(A62)') '     ==> The code_dt is larger than cfl_dt*cfl_cushion_upper.'
            write (*, '(A59,ES11.4)') '      ==> Decreasing code_dt to cfl_dt*cfl_cushion_middle =', cfl_dt * cfl_cushion_middle
            write (*, *) ' '
         end if
         code_dt = cfl_dt * cfl_cushion_middle
         call reset_dt
         restart_time_step = .true.
      else if (code_dt < min(cfl_dt * cfl_cushion_lower, code_dt_max)) then
         if (proc0) then
            write (*, *) ' '
            write (*, '(A30,I0,A1)') 'CHANGING TIME STEP: (istep = ', istep, ')'
            write (*, '(A22, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
            write (*, '(A22, ES10.2E2)') "   cfl_dt_ExB:"//REPEAT(' ', 50), cfl_dt_ExB
            write (*, '(A22, ES10.2E2)') "   cfl_dt_linear:"//REPEAT(' ', 50), cfl_dt_linear
            write (*, '(A22, ES10.2E2)') "   cfl_cushion_upper:"//REPEAT(' ', 50), cfl_cushion_upper
            write (*, '(A22, ES10.2E2)') "   cfl_cushion_middle:"//REPEAT(' ', 50), cfl_cushion_middle
            write (*, '(A22, ES10.2E2)') "   cfl_cushion_lower:"//REPEAT(' ', 50), cfl_cushion_lower
            write (*, '(A63)') '     ==> The code_dt is smaller than cfl_dt*cfl_cushion_lower.'
            write (*, '(A59,ES11.4)') '      ==> Increasing code_dt to cfl_dt*cfl_cushion_middle =', cfl_dt * cfl_cushion_middle
            write (*, *) ' '
         end if
         code_dt = min(cfl_dt * cfl_cushion_middle, code_dt_max)
         call reset_dt
         restart_time_step = .true.
      else
         gout = code_dt * gout
      end if

      if (proc0) call time_message(.false., time_gke(:, 7), ' ExB nonlinear advance')

   contains

      subroutine forward_transform(gk, gx)

         use stella_transforms, only: transform_ky2y, transform_kx2x
         use stella_transforms, only: transform_ky2y_xfirst, transform_kx2x_xfirst

         implicit none

         complex, dimension(:, :), intent(in) :: gk
         real, dimension(:, :), intent(out) :: gx

         if (yfirst) then
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
            ! TODO DSO: coordinate change for shearing
            call swap_kxky(gk, g0k_swap)
            call transform_ky2y(g0k_swap, g0kxy)
            call transform_kx2x(g0kxy, gx)
         else
            call transform_kx2x_xfirst(gk, g0xky)
            g0xky = g0xky * prefac
            call transform_ky2y_xfirst(g0xky, gx)
         end if

      end subroutine forward_transform

   end subroutine advance_ExB_nonlinearity

!    subroutine advance_parallel_nonlinearity(g, gout, restart_time_step)

!       use constants, only: zi
!       use mp, only: proc0, min_allreduce, mp_abort
!       use mp, only: scope, allprocs, subprocs
!       use stella_layouts, only: vmu_lo, xyz_lo
!       use stella_layouts, only: iv_idx, imu_idx, is_idx
!       use job_manage, only: time_message
!       use finite_differences, only: second_order_centered_zed
!       use finite_differences, only: third_order_upwind
!       use redistribute, only: gather, scatter
!       use store_arrays_fields, only: phi, phi_corr_QN, phi_corr_GA
!       use stella_transforms, only: transform_ky2y, transform_y2ky
!       use stella_transforms, only: transform_kx2x, transform_x2kx
!       use stella_time, only: cfl_dt_parallel, cfl_dt_linear, code_dt, code_dt_max
!       use parameters_numerical, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
!       use z_grid, only: nzgrid, delzed, ntubes
!       use extended_zgrid, only: neigen, nsegments, ikxmod
!       use extended_zgrid, only: iz_low, iz_up
!       use extended_zgrid, only: periodic
!       use parameters_physics, only: full_flux_surface, radial_variation
!       use grids_kxky, only: akx, aky, rho_clamped
!       use parameters_kxky_grid, only: nakx, naky, nx, ny, ikx_max
!       use calculations_kxky, only: swap_kxky, swap_kxky_back
!       use velocity_grids, only: nvpa, nmu
!       use velocity_grids, only: dvpa, vpa, mu
!       use gyro_averages, only: gyro_average
!       use parallel_streaming, only: stream_sign
!       use dist_redistribute, only: xyz2vmu
!       use file_utils, only: runtype_option_switch, runtype_multibox
!       use extended_zgrid, only: fill_zed_ghost_zones

!       implicit none

!       complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
!       complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
!       logical, intent(out) :: restart_time_step

!       integer :: ivmu, ixyz
!       integer :: iz, it, iv, imu, is
!       integer :: iky, ie, iseg
!       integer :: advect_sign
!       real :: cfl_dt
!       real, dimension(:), allocatable :: dgdv
!       real, dimension(:, :, :, :, :), allocatable :: g0xy
!       real, dimension(:, :, :), allocatable :: gxy_vmulocal
!       real, dimension(:, :), allocatable :: g1xy, advect_speed
!       complex, dimension(2) :: gleft, gright
!       complex, dimension(:, :, :, :), allocatable :: phi_gyro, dphidz
!       complex, dimension(:, :), allocatable :: g0k, g0kxy, g0k_swap
!       complex, dimension(:, :), allocatable :: tmp
      
!       ! WARNING this routine will probably break if neigen_max = 0
!       ! which happens when be set grid_option = 'range'

!       ! alpha-component of magnetic drift (requires ky -> y)
!       if (proc0) call time_message(.false., time_parallel_nl(:, 1), ' parallel nonlinearity advance')

!       ! Initialize cfl_dt_parallel
!       cfl_dt_parallel = 10000000.

!       restart_time_step = .false.

!       ! overview:
!       ! need g and d<phi>/dz in (x,y) space in
!       ! order to upwind dg/dvpa
!       ! 1) transform d<phi>/dz from (kx,ky) to (x,y). layout: vmu_lo
!       ! 2) need sign of parnl advection in xyz_lo (since dg/dvpa
!       !    requires vpa local), so d<phi>/dz(vmu_lo) --> d<phi>/dz(xyz_lo)
!       ! 3) transform g from (kx,ky) to (x,y). layout: vmu_lo
!       ! 4) dg/dvpa requires vpa local, so g(vmu_lo) --> g(xyz_lo)
!       ! 5) calculate dg/dvpa
!       ! 6) multiply dg/dvpa with d<phi>/dz
!       ! 7) product(xyz_lo) --> product(vmu_lo)
!       ! 8) inverse transform product(vmu_lo)

!       allocate (g0k(naky, nakx))
!       allocate (g0xy(ny, nx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
!       allocate (g0kxy(ny, ikx_max))
!       if (radial_variation) allocate (g1xy(ny, nx))
!       allocate (phi_gyro(naky, nakx, -nzgrid:nzgrid, ntubes))
!       allocate (dphidz(naky, nakx, -nzgrid:nzgrid, ntubes))
!       allocate (g0k_swap(2 * naky - 1, ikx_max))
!       allocate (tmp(size(gout, 1), size(gout, 2)))

!       ! get d<phi>/dz in vmu_lo
!       ! we will need to transform it to real-space
!       ! as its sign is needed for upwinding of dg/dvpa
!       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!          iv = iv_idx(vmu_lo, ivmu)
!          imu = imu_idx(vmu_lo, ivmu)
!          is = is_idx(vmu_lo, ivmu)

!          ! construct <phi>
!          dphidz = phi
!          if (radial_variation) dphidz = dphidz + phi_corr_QN
!          call gyro_average(dphidz, ivmu, phi_gyro)
!          if (radial_variation) phi_gyro = phi_gyro + phi_corr_GA(:, :, :, :, ivmu)

!          do iky = 1, naky
!             do it = 1, ntubes
!                do ie = 1, neigen(iky)
!                   do iseg = 1, nsegments(ie, iky)
!                      ! first fill in ghost zones at boundaries in g(z)
!                      call fill_zed_ghost_zones(it, iseg, ie, iky, phi_gyro, gleft, gright)
!                      ! now get d<phi>/dz
!                      call second_order_centered_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
!                                                     phi_gyro(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
!                                                     delzed(0), stream_sign(iv), gleft, gright, periodic(iky), &
!                                                     dphidz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
!                   end do
!                end do
!             end do
!          end do

!          if (radial_variation) then
!             do it = 1, ntubes
!                do iz = -nzgrid, nzgrid
!                   ! use reality to swap from ky >= 0, all kx to kx >= 0 , all ky
!                   call swap_kxky(dphidz(:, :, iz, it), g0k_swap)
!                   ! transform in y
!                   call transform_ky2y(g0k_swap, g0kxy)
!                   ! transform in x
!                   call transform_kx2x(g0kxy, g1xy)
!                   g0xy(:, :, iz, it, ivmu) = g1xy * (par_nl_fac(iz, is) + d_par_nl_fac_dr(iz, is) * spread(rho_clamped, 1, ny))

!                   g0k = zi * spread(aky, 2, nakx) * phi_gyro(:, :, iz, it)
!                   call swap_kxky(g0k, g0k_swap)
!                   call transform_ky2y(g0k_swap, g0kxy)
!                   call transform_kx2x(g0kxy, g1xy)
!                   g0xy(:, :, iz, it, ivmu) = g0xy(:, :, iz, it, ivmu) &
!                                              + vpa(iv) * g1xy * (par_nl_drifty(iz) + d_par_nl_drifty_dr(iz) * spread(rho_clamped, 1, ny))

!                   g0k = zi * spread(akx, 1, naky) * phi_gyro(:, :, iz, it)
!                   call swap_kxky(g0k, g0k_swap)
!                   call transform_ky2y(g0k_swap, g0kxy)
!                   call transform_kx2x(g0kxy, g1xy)
!                   g0xy(:, :, iz, it, ivmu) = g0xy(:, :, iz, it, ivmu) &
!                                              + vpa(iv) * g1xy * (par_nl_driftx(iz) + d_par_nl_driftx_dr(iz) * spread(rho_clamped, 1, ny))

!                   g0xy(:, :, iz, it, ivmu) = g0xy(:, :, iz, it, ivmu) &
!                                              + vpa(iv) * mu(imu) * (par_nl_curv(iz, is) + d_par_nl_curv_dr(iz, is) * spread(rho_clamped, 1, ny))

!                end do
!             end do
!          else
!             do it = 1, ntubes
!                do iz = -nzgrid, nzgrid
!                   g0k = dphidz(:, :, iz, it) * par_nl_fac(iz, is) + vpa(iv) * mu(imu) * par_nl_curv(iz, is) &
!                         + zi * vpa(iv) * phi_gyro(:, :, iz, it) * (spread(akx, 1, naky) * par_nl_driftx(iz) &
!                                                                    + spread(aky, 2, nakx) * par_nl_drifty(iz))
!                   ! use reality to swap from ky >= 0, all kx to kx >= 0 , all ky
!                   call swap_kxky(g0k, g0k_swap)
!                   ! transform in y
!                   call transform_ky2y(g0k_swap, g0kxy)
!                   ! transform in x
!                   call transform_kx2x(g0kxy, g0xy(:, :, iz, it, ivmu))
!                end do
!             end do
!          end if
!       end do

!       ! do not need phi_gyro or dphidz  again so deallocate
!       deallocate (phi_gyro, dphidz)
!       deallocate (g0k)
!       if (allocated(g1xy)) deallocate (g1xy)

!       allocate (gxy_vmulocal(nvpa, nmu, xyz_lo%llim_proc:xyz_lo%ulim_alloc))
!       allocate (advect_speed(nmu, xyz_lo%llim_proc:xyz_lo%ulim_alloc))

!       ! we now have the advection velocity in vpa in (x,y) space
!       ! next redistribute it so that (vpa,mu) are local
!       if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
!       call scatter(xyz2vmu, g0xy, gxy_vmulocal)
!       if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
!       ! advect_speed does not depend on vpa
!       advect_speed = gxy_vmulocal(1, :, :)

!       ! transform g from (kx,ky) to (x,y)
!       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!          do it = 1, ntubes
!             do iz = -nzgrid, nzgrid
!                call swap_kxky(g(:, :, iz, it, ivmu), g0k_swap)
!                ! transform in y
!                call transform_ky2y(g0k_swap, g0kxy)
!                ! transform in x
!                call transform_kx2x(g0kxy, g0xy(:, :, iz, it, ivmu))
!             end do
!          end do
!       end do

!       ! redistribute so that (vpa,mu) local
!       if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
!       call scatter(xyz2vmu, g0xy, gxy_vmulocal)
!       if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')

!       allocate (dgdv(nvpa))

!       ! we now need to form dg/dvpa and obtain product of dg/dvpa with advection speed
!       do ixyz = xyz_lo%llim_proc, xyz_lo%ulim_proc
!          do imu = 1, nmu
!             ! advect_sign set to +/- 1 depending on sign of the parallel nonlinearity
!             ! advection velocity
!             ! NB: advect_sign = -1 corresponds to positive advection velocity
!             advect_sign = int(sign(1.0, advect_speed(imu, ixyz)))
!             call third_order_upwind(1, gxy_vmulocal(:, imu, ixyz), dvpa, advect_sign, dgdv)
!             gxy_vmulocal(:, imu, ixyz) = dgdv * advect_speed(imu, ixyz)
!             cfl_dt_parallel = min(cfl_dt_parallel, dvpa / abs(advect_speed(imu, ixyz)))
!          end do
!       end do

!       ! finished with dgdv and advect_speed
!       deallocate (dgdv, advect_speed)

!       ! now that we have the full parallel nonlinearity in (x,y)-space
!       ! need to redistribute so that (x,y) local for transforms
!       if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
!       call gather(xyz2vmu, gxy_vmulocal, g0xy)
!       if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')

!       ! finished with gxy_vmulocal
!       deallocate (gxy_vmulocal)

!       ! g0xy is parallel nonlinearity term with (x,y) on processor
!       ! need to inverse Fourier transform
!       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!          do it = 1, ntubes
!             do iz = -nzgrid, nzgrid
!                call transform_x2kx(g0xy(:, :, iz, it, ivmu), g0kxy)
!                if (full_flux_surface) then
!                   gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + code_dt * g0kxy
!                else
!                   call transform_y2ky(g0kxy, g0k_swap)
!                   call swap_kxky_back(g0k_swap, tmp)
!                   gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + code_dt * tmp
!                end if
!             end do
!          end do
!       end do
!       deallocate (g0k_swap, g0kxy, g0xy)

!       if (runtype_option_switch == runtype_multibox) call scope(allprocs)

!       call min_allreduce(cfl_dt_parallel)

!       if (runtype_option_switch == runtype_multibox) call scope(subprocs)

!       !> check estimated cfl_dt to see if the time step size needs to be changed
!       cfl_dt = min(cfl_dt_parallel, cfl_dt_linear)
!       if (code_dt > cfl_dt * cfl_cushion_upper) then
!          if (proc0) then
!             write (*, *) ' '
!             write (*, *) 'CHANGING TIME STEP:'
!             write (*, '(A22, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
!             write (*, '(A22, ES10.2E2)') "   cfl_dt_parallel:"//REPEAT(' ', 50), cfl_dt_parallel
!             write (*, '(A22, ES10.2E2)') "   cfl_dt_linear:"//REPEAT(' ', 50), cfl_dt_linear
!             write (*, '(A22, ES10.2E2)') "   cfl_cushion_upper:"//REPEAT(' ', 50), cfl_cushion_upper
!             write (*, '(A22, ES10.2E2)') "   cfl_cushion_middle:"//REPEAT(' ', 50), cfl_cushion_middle
!             write (*, '(A22, ES10.2E2)') "   cfl_cushion_lower:"//REPEAT(' ', 50), cfl_cushion_lower
!             write (*, '(A62)') '     ==> The code_dt is larger than cfl_dt*cfl_cushion_upper.'
!             write (*, '(A59,ES11.4)') '      ==> Decreasing code_dt to cfl_dt*cfl_cushion_middle =', cfl_dt * cfl_cushion_middle
!             write (*, *) ' '
!          end if
!          code_dt = cfl_dt * cfl_cushion_middle
!          call reset_dt
!          restart_time_step = .true.
!       else if (code_dt < min(cfl_dt * cfl_cushion_lower, code_dt_max)) then
!          if (proc0) then
!             write (*, *) ' '
!             write (*, *) 'CHANGING TIME STEP:'
!             write (*, '(A22, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
!             write (*, '(A22, ES10.2E2)') "   cfl_dt_parallel:"//REPEAT(' ', 50), cfl_dt_parallel
!             write (*, '(A22, ES10.2E2)') "   cfl_dt_linear:"//REPEAT(' ', 50), cfl_dt_linear
!             write (*, '(A22, ES10.2E2)') "   cfl_cushion_upper:"//REPEAT(' ', 50), cfl_cushion_upper
!             write (*, '(A22, ES10.2E2)') "   cfl_cushion_middle:"//REPEAT(' ', 50), cfl_cushion_middle
!             write (*, '(A22, ES10.2E2)') "   cfl_cushion_lower:"//REPEAT(' ', 50), cfl_cushion_lower
!             write (*, '(A63)') '     ==> The code_dt is smaller than cfl_dt*cfl_cushion_lower.'
!             write (*, '(A59,ES11.4)') '      ==> Increasing code_dt to cfl_dt*cfl_cushion_middle =', cfl_dt * cfl_cushion_middle
!             write (*, *) ' '
!          end if
!          code_dt = min(cfl_dt * cfl_cushion_middle, code_dt_max)
!          call reset_dt
!          restart_time_step = .true.
! !    else
! !       gout = code_dt*gout
!       end if

!       if (proc0) call time_message(.false., time_parallel_nl(:, 1), ' parallel nonlinearity advance')

!    end subroutine advance_parallel_nonlinearity

   subroutine advance_radial_variation(g, gout)

      use mp, only: mp_abort, proc0
      use job_manage, only: time_message
      use calculations_kxky_derivatives, only: get_dchidy
      use store_arrays_fields, only: phi, apar, bpar
      use store_arrays_fields, only: phi_corr_QN, phi_corr_GA
!   use store_arrays_fields, only: apar_corr_QN, apar_corr_GA
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
      use z_grid, only: nzgrid, ntubes
      use parameters_kxky_grid, only: nakx, naky
      use calculations_kxky, only: multiply_by_rho
      use gyro_averages, only: gyro_average, gyro_average_j1
      use parameters_physics, only: full_flux_surface, fphi
      use parameters_physics, only: include_parallel_streaming, include_mirror
      use store_arrays_distribution_fn, only: wdriftx_phi, wdrifty_phi
      use store_arrays_distribution_fn, only: wdriftpx_g, wdriftpy_g
      use store_arrays_distribution_fn, only: wdriftpx_phi, wdriftpy_phi !, adiabatic_phi
      use store_arrays_distribution_fn, only: wstar, wstarp
      use mirror_terms, only: add_mirror_radial_variation
      use flow_shear, only: prl_shear, prl_shear_p
      use parallel_streaming, only: add_parallel_streaming_radial_variation

      use calculations_kxky_derivatives, only: get_dgdy, get_dgdx

      use store_arrays_useful, only: time_gke

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      integer :: ia, ivmu, iv, imu, is, iz, it

      complex, dimension(:, :), allocatable :: g0k, g1k, g0a
      complex, dimension(:, :, :, :, :), allocatable :: g_corr

      allocate (g0k(naky, nakx))

      allocate (g1k(naky, nakx))
      allocate (g0a(naky, nakx))

      if (debug) write (*, *) 'time_advance::solve_gke::advance_radial_variation'

      if (proc0) call time_message(.false., time_gke(:, 10), ' radial variation advance')

      if (include_mirror .or. include_parallel_streaming) then
         allocate (g_corr(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         g_corr = 0.
      end if

      !grab the mirror and parallel streaming corrections here to save on FFTs
      if (include_mirror) then
         call add_mirror_radial_variation(g, g_corr)
      end if
      if (include_parallel_streaming) then
         call add_parallel_streaming_radial_variation(g, g_corr, gout)
      end if

      if (full_flux_surface) then
         ! FLAG -- ADD SOMETHING HERE
         call mp_abort('wstarp term not yet setup for full_flux_surface = .true. aborting.')
      end if

      ia = 1
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0k = 0.

               !wstar variation
               call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0a)
               g0k = g0k + g0a * wstarp(ia, iz, ivmu)

               !radial variation in ExB nonlinearity is handled in advance_ExB_nonlinearity

               !wdrift(x/y) - g

               call get_dgdx(g(:, :, iz, it, ivmu), g0a)
               g0k = g0k + g0a * wdriftpx_g(ia, iz, ivmu)

               call get_dgdy(g(:, :, iz, it, ivmu), g0a)
               g0k = g0k + g0a * wdriftpy_g(ia, iz, ivmu)

               !wdrift - phi
               call get_dgdx(phi(:, :, iz, it), g1k)
               !wdriftx variation
               call gyro_average(g1k, iz, ivmu, g0a)
               g0k = g0k + g0a * wdriftpx_phi(ia, iz, ivmu)

               call get_dgdy(phi(:, :, iz, it), g1k)
               !wdrifty variation
               call gyro_average(g1k, iz, ivmu, g0a)
               g0k = g0k + g0a * wdriftpy_phi(ia, iz, ivmu)

               !prl_shear variation
               g0k = g0k + g0a * prl_shear_p(ia, iz, ivmu)

               !mirror term and/or parallel streaming
               if (include_mirror .or. include_parallel_streaming) then
                  g0k = g0k + g_corr(:, :, iz, it, ivmu)
               end if

               !inverse and forward transforms
               call multiply_by_rho(g0k)

               !
               !quasineutrality/gyroaveraging
               !
               call gyro_average(phi_corr_QN(:, :, iz, it), iz, ivmu, g0a)
               g0a = fphi * (g0a + phi_corr_GA(:, :, iz, it, ivmu))

               !wstar - gyroaverage/quasineutrality variation
               call get_dgdy(g0a, g1k)
               g0k = g0k + g1k * wstar(ia, iz, ivmu)

               !wdrifty gyroaverage/quasineutrality variation
               g0k = g0k + g1k * wdrifty_phi(ia, iz, ivmu)

               !prl_shear gyroaverage/quasineutrality variation
               g0k = g0k + g1k * prl_shear(ia, iz, ivmu)

               !wdriftx gyroaverage/quasineutrality variation
               call get_dgdx(g0a, g1k)
               g0k = g0k + g1k * wdriftx_phi(ia, iz, ivmu)

!              !wdriftx F_M/T_s variation
!              call gyro_average (phi(:, :, iz, it), iz, ivmu, g0a)
!              g0a = adiabatic_phi(ia, iz, ivmu) * g0a
!              call multiply_by_rho(g0a)
!              call get_dgdx(g0a, g1k)
!              g0k = g0k + g1k * wdriftx_phi(ia, iz, ivmu)

               gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + g0k
            end do
         end do
      end do

      deallocate (g0k, g1k, g0a)
      if (allocated(g_corr)) deallocate (g_corr)

      if (proc0) call time_message(.false., time_gke(:, 10), ' radial variation advance')

   end subroutine advance_radial_variation

   !-------------------------------------------------------------------------------
   !                           IMPLICIT TIME ADVANCE SUBROUTINE
   !-------------------------------------------------------------------------------

   subroutine advance_implicit(istep, phi, apar, bpar, g)

      use mp, only: proc0
      use job_manage, only: time_message
      use stella_layouts, only: vmu_lo
      use z_grid, only: nzgrid
      use dissipation, only: hyper_dissipation
      use hyper, only: advance_hyper_dissipation
      use parameters_physics, only: include_parallel_streaming
      use parameters_physics, only: radial_variation, full_flux_surface
      use parameters_physics, only: include_mirror, prp_shear_enabled
      use parameters_numerical, only: stream_implicit, mirror_implicit, drifts_implicit
      use implicit_solve, only: advance_implicit_terms
      use fields, only: advance_fields, fields_updated
      use mirror_terms, only: advance_mirror_implicit
      use dissipation, only: collisions_implicit, include_collisions
      use dissipation, only: advance_collisions_implicit
      use flow_shear, only: advance_perp_flow_shear
      use parameters_multibox, only: rk_step
      use parameters_numerical, only: flip_flop
      use store_arrays_useful, only: time_gke
      implicit none

      integer, intent(in) :: istep
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g

      ! start the timer for the implicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 9), ' implicit')

      ! reverse the order of operations every time step
      ! as part of alternating direction operator splitting
      ! this is needed to ensure 2nd order accuracy in time
      ! if (mod(istep,2)==0) then
      ! g^{*} (coming from explicit solve) is input
      ! get g^{**}, with g^{**}-g^{*} due to mirror term

      if (rk_step) call mb_communicate(g)

      if (mod(istep, 2) == 1 .or. .not. flip_flop) then

         if (prp_shear_enabled) then
            call advance_perp_flow_shear(g)
            fields_updated = .false.
         end if

         if (hyper_dissipation) then
            call advance_hyper_dissipation(g)
            fields_updated = .false.
         end if

         if (collisions_implicit .and. include_collisions) then
            call advance_fields(g, phi, apar, bpar, dist='g')
            call advance_collisions_implicit(mirror_implicit, phi, apar, bpar, g)
            fields_updated = .false.
         end if

         if (mirror_implicit .and. include_mirror) then
            call advance_mirror_implicit(collisions_implicit, g, apar)
            fields_updated = .false.
         end if

         call advance_fields(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 
         ! g^{**} is input
         ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
         if (stream_implicit .and. include_parallel_streaming) then
            call advance_implicit_terms(g, phi, apar, bpar)
            if (radial_variation .or. full_flux_surface) fields_updated = .false.
         end if

         ! update the fields if not already updated
         call advance_fields(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 
      else

         ! get updated fields corresponding to advanced g
         ! note that hyper-dissipation and mirror advances
         ! depended only on g and so did not need field update
         call advance_fields(g, phi, apar, bpar, dist='g')
         fields_updated = .true. 

         ! g^{**} is input
         ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
         if (stream_implicit .and. include_parallel_streaming) then
            call advance_implicit_terms(g, phi, apar, bpar)
            if (radial_variation .or. full_flux_surface) fields_updated = .false.
         end if

         if (mirror_implicit .and. include_mirror) then
            call advance_mirror_implicit(collisions_implicit, g, apar)
            fields_updated = .false.
         end if

         if (collisions_implicit .and. include_collisions) then
            call advance_fields(g, phi, apar, bpar, dist='g')
            call advance_collisions_implicit(mirror_implicit, phi, apar, bpar, g)
            fields_updated = .false.
         end if

         if (hyper_dissipation) then
            call advance_hyper_dissipation(g)
            fields_updated = .false.
         end if

         if (prp_shear_enabled) then
            call advance_perp_flow_shear(g)
            fields_updated = .false.
         end if

      end if

      ! stop the timer for the implict part of the solve
      if (proc0) call time_message(.false., time_gke(:, 9), ' implicit')

   end subroutine advance_implicit

   !******************************************************************************
   !                           MULTIBOX COMMUNICATION SUBROUTINE
   !******************************************************************************

   subroutine mb_communicate(g_in)

      use mp, only: job
      use stella_layouts, only: vmu_lo
      use z_grid, only: nzgrid
      use multibox, only: multibox_communicate, apply_radial_boundary_conditions
      use parameters_multibox, only: use_dirichlet_bc
      use fields, only: fields_updated, advance_fields
      use store_arrays_fields, only: phi, apar, bpar
      use file_utils, only: runtype_option_switch, runtype_multibox

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g_in

      if (runtype_option_switch == runtype_multibox) then
         if (job /= 1) then
            call advance_fields(g_in, phi, apar, bpar, dist='g')
         end if

         call multibox_communicate(g_in)

         if (job == 1) then
            fields_updated = .false.
            call advance_fields(g_in, phi, apar, bpar, dist='g')
         end if
      else if (use_dirichlet_BC) then
         call apply_radial_boundary_conditions(g_in)
         fields_updated = .false.
         call advance_fields(g_in, phi, apar, bpar, dist='g')
      end if

   end subroutine mb_communicate

   !******************************************************************************
   !                           FINISH TIME ADVANCE SUBROUTINE
   !******************************************************************************
   subroutine finish_time_advance

      use stella_transforms, only: finish_transforms
      use parameters_physics, only: full_flux_surface
      use extended_zgrid, only: finish_extended_zgrid
      use parallel_streaming, only: finish_parallel_streaming
      use mirror_terms, only: finish_mirror
      use flow_shear, only: finish_flow_shear
      use neoclassical_terms, only: finish_neoclassical_terms
      use dissipation, only: finish_dissipation
      use arrays_drifts, only: finish_wstar, finish_wdrift
      use parallel_nonlinearity, only: finish_parallel_nonlinearity
      implicit none

      if (full_flux_surface) call finish_transforms
      call finish_dissipation
      call finish_parallel_nonlinearity 
      call finish_wstar 
      call finish_wdrift 
      call finish_parallel_streaming
      call finish_flow_shear
      call finish_mirror
      call finish_neoclassical_terms
      call deallocate_arrays

      time_advance_initialized = .false.
      readinit = .false.

   end subroutine finish_time_advance


   subroutine deallocate_arrays

      use store_arrays_distribution_fn, only: g0, g1, g2, g3

      implicit none

      if (allocated(g0)) deallocate (g0)
      if (allocated(g1)) deallocate (g1)
      if (allocated(g2)) deallocate (g2)
      if (allocated(g3)) deallocate (g3)

   end subroutine deallocate_arrays

end module time_advance
