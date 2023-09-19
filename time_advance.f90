
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

   interface checksum
      module procedure checksum_field
      module procedure checksum_dist
   end interface

   logical :: time_advance_initialized = .false.

   logical :: wdriftinit = .false.
   logical :: wstarinit = .false.
   logical :: parnlinit = .false.
   logical :: readinit = .false.
   logical :: radialinit = .false.
   logical :: driftimpinit = .false.

   ! if .true., dist fn is represented on alpha grid
   ! if .false., dist fn is given on k-alpha grid
   ! default is .false.; will only ever be set to
   ! .true. during full_flux_surface simulations
!  logical :: alpha_space = .false.

   integer :: explicit_option_switch
   integer, parameter :: explicit_option_rk3 = 1, &
                         explicit_option_rk2 = 2, &
                         explicit_option_rk4 = 3

   real :: xdriftknob, ydriftknob, wstarknob
   logical :: flip_flop

   ! factor multiplying parallel nonlinearity
   real, dimension(:, :), allocatable :: par_nl_fac, d_par_nl_fac_dr
   ! factor multiplying higher order linear term in parallel acceleration
   real, dimension(:, :), allocatable :: par_nl_curv, d_par_nl_curv_dr
   real, dimension(:), allocatable :: par_nl_driftx, par_nl_drifty
   real, dimension(:), allocatable :: d_par_nl_driftx_dr, d_par_nl_drifty_dr

   ! needed for timing various pieces of gke solve
   real, dimension(2, 10) :: time_gke = 0.
   real, dimension(2, 2) :: time_parallel_nl = 0.

   logical :: debug = .false.

contains

   subroutine init_time_advance

      use mp, only: proc0
      use physics_flags, only: radial_variation
      use physics_flags, only: include_parallel_nonlinearity
      use neoclassical_terms, only: init_neoclassical_terms
      use dissipation, only: init_collisions, include_collisions
      use parallel_streaming, only: init_parallel_streaming
      use mirror_terms, only: init_mirror
      use flow_shear, only: init_flow_shear
      use sources, only: init_quasineutrality_source, init_source_timeaverage

      implicit none

      if (time_advance_initialized) return
      time_advance_initialized = .true.

      debug = debug .and. proc0

      !> read time_advance_knobs namelist from the input file;
      !> sets the explicit time advance option, as well as allows for scaling of
      !> the x and y components of the magnetic drifts and of the drive term
      if (debug) write (6, *) 'time_advance::init_time_advance::read_parameters'
      call read_parameters
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

   subroutine read_parameters

      use file_utils, only: error_unit, input_unit_exist
      use text_options, only: text_option, get_option_value
      use mp, only: proc0, broadcast
      use run_parameters, only: fully_explicit

      implicit none

      logical :: taexist

      type(text_option), dimension(4), parameter :: explicitopts = &
                                                    (/text_option('default', explicit_option_rk3), &
                                                      text_option('rk3', explicit_option_rk3), &
                                                      text_option('rk2', explicit_option_rk2), &
                                                      text_option('rk4', explicit_option_rk4)/)
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

         in_file = input_unit_exist("time_advance_knobs", taexist)
         if (taexist) read (unit=in_file, nml=time_advance_knobs)

         ierr = error_unit()
         call get_option_value &
            (explicit_option, explicitopts, explicit_option_switch, &
             ierr, "explicit_option in time_advance_knobs")
      end if

      call broadcast(explicit_option_switch)
      call broadcast(xdriftknob)
      call broadcast(ydriftknob)
      call broadcast(wstarknob)
      call broadcast(flip_flop)

      if (fully_explicit) flip_flop = .false.

   end subroutine read_parameters

   subroutine init_wdrift

      use mp, only: mp_abort
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
      use stella_geometry, only: geo_surf, q_as_x
      use stella_geometry, only: dxdXcoord, drhodpsi, dydalpha
      use vpamu_grids, only: vpa, vperp2
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use neoclassical_terms, only: include_neoclassical_terms
      use neoclassical_terms, only: dphineo_dzed, dphineo_drho, dphineo_dalpha
      use neoclassical_terms, only: dfneo_dvpa, dfneo_dzed, dfneo_dalpha

      implicit none

      integer :: ivmu, iv, imu, is
      real :: fac
      real, dimension(:, :), allocatable :: wcvdrifty, wgbdrifty
      real, dimension(:, :), allocatable :: wcvdriftx, wgbdriftx

      if (wdriftinit) return
      wdriftinit = .true.

      !> allocate wdriftx_phi, the factor multiplying dphi/dx in the magnetic drift term
      if (.not. allocated(wdriftx_phi)) then
         allocate (wdriftx_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         wdriftx_phi = 0.0
      end if
      !> allocate wdrifty_phi, the factor multiplying dphi/dy in the magnetic drift term
      if (.not. allocated(wdrifty_phi)) then
         allocate (wdrifty_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         wdrifty_phi = 0.0
      end if
      !> allocate wdriftx_g, the factor multiplying dg/dx in the magnetic drift term
      if (.not. allocated(wdriftx_g)) then
         allocate (wdriftx_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         wdriftx_g = 0.0
      end if
      !> allocate wdrifty_g, the factor multiplying dg/dy in the magnetic drift term
      if (.not. allocated(wdrifty_g)) then
         allocate (wdrifty_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         wdrifty_g = 0.0
      end if

      allocate (wcvdrifty(nalpha, -nzgrid:nzgrid))
      allocate (wgbdrifty(nalpha, -nzgrid:nzgrid))
      allocate (wcvdriftx(nalpha, -nzgrid:nzgrid))
      allocate (wgbdriftx(nalpha, -nzgrid:nzgrid))

      ! FLAG -- need to deal with shat=0 case.  ideally move away from q as x-coordinate
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         fac = -ydriftknob * 0.5 * code_dt * spec(is)%tz_psi0
         !> this is the curvature drift piece of wdrifty with missing factor of vpa
         !> vpa factor is missing to avoid singularity when including
         !> non-Maxwellian corrections to equilibrium
         wcvdrifty = fac * cvdrift * vpa(iv)
         !> this is the grad-B drift piece of wdrifty
         wgbdrifty = fac * gbdrift * 0.5 * vperp2(:, :, imu)
         wdrifty_g(:, :, ivmu) = wcvdrifty * vpa(iv) + wgbdrifty
         !> if including neoclassical correction to equilibrium Maxwellian,
         !> then add in v_E^{nc} . grad y dg/dy coefficient here
         if (include_neoclassical_terms) then
            wdrifty_g(:, :, ivmu) = wdrifty_g(:, :, ivmu) + code_dt * 0.5 * (gds23 * dphineo_dzed &
                                                                             + drhodpsi * dydalpha * dphineo_drho)
         end if

         wdrifty_phi(:, :, ivmu) = spec(is)%zt * (wgbdrifty + wcvdrifty * vpa(iv))

         !> if maxwwellian_normalization = .true., evolved distribution function is normalised by a Maxwellian
         !> otherwise, it is not; a Maxwellian weighting factor must thus be included
         wdrifty_phi(:, :, ivmu) = wdrifty_phi(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
         !> if including neoclassical corrections to equilibrium,
         !> add in -(Ze/m) * v_curv/vpa . grad y d<phi>/dy * dF^{nc}/dvpa term
         !> and v_E . grad z dF^{nc}/dz (here get the dphi/dy part of v_E)
         if (include_neoclassical_terms) then
            !> NB: the below neoclassical correction needs to be divided by an equilibrium Maxwellian
            wdrifty_phi(:, :, ivmu) = wdrifty_phi(:, :, ivmu) &
                                      - 0.5 * spec(is)%zt * dfneo_dvpa(:, :, ivmu) * wcvdrifty &
                                      - code_dt * 0.5 * dfneo_dzed(:, :, ivmu) * gds23
         end if

         if (q_as_x) then
            fac = -xdriftknob * 0.5 * code_dt * spec(is)%tz_psi0
         else
            fac = -xdriftknob * 0.5 * code_dt * spec(is)%tz_psi0 / geo_surf%shat
         end if
         !> this is the curvature drift piece of wdriftx with missing factor of vpa
         !> vpa factor is missing to avoid singularity when including
         !> non-Maxwellian corrections to equilibrium
         wcvdriftx = fac * cvdrift0 * vpa(iv)
         !> this is the grad-B drift piece of wdriftx
         wgbdriftx = fac * gbdrift0 * 0.5 * vperp2(:, :, imu)
         wdriftx_g(:, :, ivmu) = wcvdriftx * vpa(iv) + wgbdriftx
         !> if including neoclassical correction to equilibrium Maxwellian,
         !> then add in v_E^{nc} . grad x dg/dx coefficient here
         if (include_neoclassical_terms) then
            wdriftx_g(:, :, ivmu) = wdriftx_g(:, :, ivmu) + code_dt * 0.5 * (gds24 * dphineo_dzed &
                                                                             - dxdXcoord * dphineo_dalpha)
         end if
         wdriftx_phi(:, :, ivmu) = spec(is)%zt * (wgbdriftx + wcvdriftx * vpa(iv))
         !> otherwise, it is not; a Maxwellian weighting factor must thus be included
         wdriftx_phi(:, :, ivmu) = wdriftx_phi(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
         !> if including neoclassical corrections to equilibrium,
         !> add in (Ze/m) * v_curv/vpa . grad x d<phi>/dx * dF^{nc}/dvpa term
         !> and v_E . grad z dF^{nc}/dz (here get the dphi/dx part of v_E)
         !> and v_E . grad alpha dF^{nc}/dalpha (dphi/dx part of v_E)
         if (include_neoclassical_terms) then
            !> NB: the below neoclassical correction needs to be divided by an equilibrium Maxwellian
            wdriftx_phi(:, :, ivmu) = wdriftx_phi(:, :, ivmu) &
                                      - 0.5 * spec(is)%zt * dfneo_dvpa(:, :, ivmu) * wcvdriftx &
                                      + code_dt * 0.5 * (dfneo_dalpha(:, :, ivmu) * dxdXcoord - dfneo_dzed(:, :, ivmu) * gds24)
         end if

      end do

      deallocate (wcvdriftx, wgbdriftx, wcvdrifty, wgbdrifty)

   end subroutine init_wdrift

   subroutine init_wstar

      use mp, only: mp_abort
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use stella_time, only: code_dt
      use species, only: spec
      use zgrid, only: nzgrid
      use kt_grids, only: nalpha
      use stella_geometry, only: dydalpha, drhodpsi
      use vpamu_grids, only: vperp2, vpa
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use dist_fn_arrays, only: wstar
      use neoclassical_terms, only: include_neoclassical_terms
      use neoclassical_terms, only: dfneo_drho

      implicit none

      integer :: is, imu, iv, ivmu
      real, dimension(:, :), allocatable :: energy

      if (wstarinit) return
      wstarinit = .true.

      if (.not. allocated(wstar)) &
         allocate (wstar(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstar = 0.0

      allocate (energy(nalpha, -nzgrid:nzgrid))

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         is = is_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         iv = iv_idx(vmu_lo, ivmu)
         energy = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
         if (include_neoclassical_terms) then
            wstar(:, :, ivmu) = dydalpha * drhodpsi * wstarknob * 0.5 * code_dt &
                                * (maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is) &
                                   * (spec(is)%fprim + spec(is)%tprim * (energy - 1.5)) &
                                   - dfneo_drho(:, :, ivmu))
         else
            wstar(:, :, ivmu) = dydalpha * drhodpsi * wstarknob * 0.5 * code_dt &
                                * (spec(is)%fprim + spec(is)%tprim * (energy - 1.5)) &
                                * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
         end if
      end do

      deallocate (energy)

   end subroutine init_wstar

   subroutine init_parallel_nonlinearity

      use physics_parameters, only: rhostar
      use species, only: spec, nspec
      use zgrid, only: nztot, nzgrid
      use stella_geometry, only: geo_surf, drhodpsi, q_as_x
      use stella_geometry, only: gradpar, dbdzed, bmag
      use stella_geometry, only: cvdrift, cvdrift0
      use stella_geometry, only: dIdrho, dgradpardrho, dBdrho, d2Bdrdth
      use stella_geometry, only: dcvdriftdrho, dcvdrift0drho
      use physics_flags, only: radial_variation

      implicit none

      if (.not. allocated(par_nl_fac)) allocate (par_nl_fac(-nzgrid:nzgrid, nspec))
      ! this is the factor multiplying -dphi/dz * dg/dvpa in the parallel nonlinearity
      par_nl_fac = 0.5 * rhostar * spread(spec%stm_psi0 * spec%zt_psi0, 1, nztot) * spread(gradpar, 2, nspec)

      if (.not. allocated(par_nl_curv)) allocate (par_nl_curv(-nzgrid:nzgrid, nspec))
      ! ydriftknob is here because this term comes from bhat x curvature . grad B
      par_nl_curv = -ydriftknob * rhostar * geo_surf%rgeo * geo_surf%betaprim * drhodpsi &
                    * spread(dbdzed(1, :) * gradpar / bmag(1, :), 2, nspec) / spread(spec%zt_psi0, 1, nztot)

      if (.not. allocated(par_nl_drifty)) allocate (par_nl_drifty(-nzgrid:nzgrid))
      par_nl_drifty = 0.25 * rhostar * cvdrift(1, :)
      if (.not. allocated(par_nl_driftx)) allocate (par_nl_driftx(-nzgrid:nzgrid))
      if (q_as_x) then
         par_nl_driftx = 0.25 * rhostar * cvdrift0(1, :)
      else
         par_nl_driftx = 0.25 * rhostar * cvdrift0(1, :) / geo_surf%shat
      end if

      if (radial_variation) then
         if (.not. allocated(d_par_nl_fac_dr)) allocate (d_par_nl_fac_dr(-nzgrid:nzgrid, nspec))
         ! this is the factor multiplying -dphi/dz * dg/dvpa in the parallel nonlinearity
         d_par_nl_fac_dr = 0.5 * rhostar * spread(spec%stm_psi0 * spec%zt_psi0, 1, nztot) * spread(dgradpardrho, 2, nspec)

         if (.not. allocated(d_par_nl_curv_dr)) allocate (d_par_nl_curv_dr(-nzgrid:nzgrid, nspec))
         ! ydriftknob is here because this term comes from bhat x curvature . grad B
         ! handle terms with no zeroes
         d_par_nl_curv_dr = par_nl_curv * (dIdrho / geo_surf%rgeo - drhodpsi * geo_surf%d2psidr2 &
                                           - spread(dBdrho / bmag(1, :) + dgradpardrho / gradpar, 2, nspec))
         ! handle terms with possible zeroes
         d_par_nl_curv_dr = d_par_nl_curv_dr &
                            - ((ydriftknob * rhostar * geo_surf%rgeo * drhodpsi * spread(gradpar / bmag(1, :), 2, nspec)) &
                               / spread(spec%zt_psi0, 1, nztot)) &
                            * (geo_surf%betadbprim * spread(dbdzed(1, :), 2, nspec) &
                               + geo_surf%betaprim * spread(d2Bdrdth, 2, nspec))

         if (.not. allocated(d_par_nl_drifty_dr)) allocate (d_par_nl_drifty_dr(-nzgrid:nzgrid))
         d_par_nl_drifty_dr = 0.25 * rhostar * dcvdriftdrho(1, :)
         if (.not. allocated(d_par_nl_drifty_dr)) allocate (d_par_nl_driftx_dr(-nzgrid:nzgrid))
         if (q_as_x) then
            d_par_nl_driftx_dr = 0.25 * rhostar * dcvdrift0drho(1, :)
         else
            d_par_nl_driftx_dr = 0.25 * rhostar * dcvdrift0drho(1, :) / geo_surf%shat
         end if
      end if

      parnlinit = .true.

   end subroutine init_parallel_nonlinearity

   subroutine init_radial_variation
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use stella_time, only: code_dt
      use species, only: spec, pfac
      use zgrid, only: nzgrid
      use kt_grids, only: nalpha
      use stella_geometry, only: drhodpsi, dydalpha, gfac
      use stella_geometry, only: dBdrho, geo_surf, q_as_x
      use stella_geometry, only: dcvdriftdrho, dcvdrift0drho
      use stella_geometry, only: dgbdriftdrho, dgbdrift0drho
      use vpamu_grids, only: vperp2, vpa, mu
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use dist_fn_arrays, only: wstarp
      use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
      use dist_fn_arrays, only: wdriftpx_g, wdriftpy_g
      use dist_fn_arrays, only: wdriftpx_phi, wdriftpy_phi!, adiabatic_phi
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
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use dist_fn_arrays, only: g0, g1, g2, g3

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
      if (.not. allocated(g3) .and. explicit_option_switch == explicit_option_rk4) then
         allocate (g3(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         g3 = 0.
      else
         allocate (g3(1, 1, 1, 1, 1))
      end if

   end subroutine allocate_arrays

   subroutine init_cfl

      use mp, only: proc0, nproc, max_allreduce, min_allreduce
      use mp, only: scope, allprocs, subprocs
      use dist_fn_arrays, only: wdriftx_g, wdrifty_g
      use stella_time, only: code_dt, write_dt, cfl_dt_linear
      use run_parameters, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
      use physics_flags, only: radial_variation, prp_shear_enabled
      use zgrid, only: delzed
      use vpamu_grids, only: dvpa
      use kt_grids, only: akx, aky, nx, rho
      use run_parameters, only: stream_implicit, mirror_implicit, drifts_implicit
      use parallel_streaming, only: stream
      use parallel_streaming, only: stream_rad_var1, stream_rad_var2
      use mirror_terms, only: mirror
      use flow_shear, only: prl_shear, shift_times
      use file_utils, only: runtype_option_switch, runtype_multibox
      use dissipation, only: include_collisions, collisions_implicit
      use dissipation, only: cfl_dt_vpadiff, cfl_dt_mudiff

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

      if (proc0) then
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
      else if (proc0) then
         call write_dt
         write (*, *)
      end if

   end subroutine init_cfl

   subroutine reset_dt

      use parallel_streaming, only: parallel_streaming_initialized
      use parallel_streaming, only: init_parallel_streaming
      use dissipation, only: init_collisions, collisions_initialized, include_collisions
      use run_parameters, only: stream_implicit, driftkinetic_implicit
      use response_matrix, only: response_matrix_initialized
      use response_matrix, only: init_response_matrix
      use mirror_terms, only: mirror_initialized
      use mirror_terms, only: init_mirror
      use flow_shear, only: flow_shear_initialized
      use flow_shear, only: init_flow_shear
      use physics_flags, only: radial_variation
      use sources, only: init_source_timeaverage
      use sources, only: init_quasineutrality_source, qn_source_initialized

      implicit none

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

   subroutine advance_stella(istep, stop_stella)

      use dist_fn_arrays, only: gold, gnew
      use fields_arrays, only: phi, apar
      use fields_arrays, only: phi_old
      use fields, only: advance_fields, fields_updated
      use run_parameters, only: fully_explicit, fully_implicit
      use multibox, only: RK_step
      use sources, only: include_qn_source, update_quasineutrality_source
      use sources, only: source_option_switch, source_option_projection
      use sources, only: source_option_krook
      use sources, only: update_tcorr_krook, project_out_zero
      use mp, only: proc0, broadcast

      implicit none

      integer, intent(in) :: istep
      logical, intent(in out) :: stop_stella

      logical :: restart_time_step, time_advance_successful
      integer :: count_restarts

      !> unless running in multibox mode, no need to worry about
      !> mb_communicate calls as the subroutine is immediately exited
      !> if not in multibox mode.
      if (.not. RK_step) then
         if (debug) write (*, *) 'time_advance::multibox'
         call mb_communicate(gnew)
      end if

      !> save value of phi
      !> for use in diagnostics (to obtain frequency)
      phi_old = phi

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
         call advance_fields(gnew, phi, apar, dist='gbar')

         ! Keep track whether any routine wants to modify the time step
         restart_time_step = .false.

         !> reverse the order of operations every time step
         !> as part of alternating direction operator splitting
         !> this is needed to ensure 2nd order accuracy in time
         if (mod(istep, 2) == 1 .or. .not. flip_flop) then

            !> Advance the explicit parts of the GKE
            if (debug) write (*, *) 'time_advance::advance_explicit'
            if (.not. fully_implicit) call advance_explicit(gnew, restart_time_step, istep)

            !> Use operator splitting to separately evolve all terms treated implicitly
            if (.not. restart_time_step .and. .not. fully_explicit) call advance_implicit(istep, phi, apar, gnew)
         else
            if (.not. fully_explicit) call advance_implicit(istep, phi, apar, gnew)
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
      call advance_fields(gnew, phi, apar, dist='gbar')

      !update the delay parameters for the Krook operator
      if (source_option_switch == source_option_krook) call update_tcorr_krook(gnew)
      if (include_qn_source) call update_quasineutrality_source

   end subroutine advance_stella

   !> advance_explicit takes as input the guiding centre distribution function
   !> in k-space and updates it to account for all of the terms in the GKE that
   !> are advanced explicitly in time
   subroutine advance_explicit(g, restart_time_step, istep)

      use mp, only: proc0
      use job_manage, only: time_message
      use zgrid, only: nzgrid
      use extended_zgrid, only: periodic, phase_shift
      use kt_grids, only: naky
      use stella_layouts, only: vmu_lo, iv_idx
      use parallel_streaming, only: stream_sign

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      integer :: ivmu, iv, sgn, iky

      !> start the timer for the explicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 8), ' explicit')

      select case (explicit_option_switch)
      case (explicit_option_rk2)
         !> SSP RK2
         call advance_explicit_rk2(g, restart_time_step, istep)
      case (explicit_option_rk3)
         !> default is SSP RK3
         call advance_explicit_rk3(g, restart_time_step, istep)
      case (explicit_option_rk4)
         !> RK4
         call advance_explicit_rk4(g, restart_time_step, istep)
      end select

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

   !> advance_expliciit_rk2 uses strong stability-preserving RK2 to advance one time step
   subroutine advance_explicit_rk2(g, restart_time_step, istep)

      use dist_fn_arrays, only: g0, g1
      use zgrid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use multibox, only: RK_step

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      integer :: icnt

      !> RK_step only true if running in multibox mode
      if (RK_step) call mb_communicate(g)

      g0 = g
      icnt = 1

      !> SSP rk3 algorithm to advance explicit part of code
      !> if GK equation written as dg/dt = rhs - vpar . grad h,
      !> solve_gke returns rhs*dt
      do while (icnt <= 2)
         select case (icnt)
         case (1)
            call solve_gke(g0, g1, restart_time_step, istep)
         case (2)
            g1 = g0 + g1
            if (RK_step) call mb_communicate(g1)
            call solve_gke(g1, g, restart_time_step, istep)
         end select
         if (restart_time_step) then
            ! If the code_dt is reset, we need to quit this loop and restart the timestep again
            icnt = 10
         else
            icnt = icnt + 1
         end if
      end do

      !> this is gbar at intermediate time level
      g = 0.5 * g0 + 0.5 * (g1 + g)

   end subroutine advance_explicit_rk2

   !> strong stability-preserving RK3
   subroutine advance_explicit_rk3(g, restart_time_step, istep)

      use dist_fn_arrays, only: g0, g1, g2
      use zgrid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use multibox, only: RK_step

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      integer :: icnt

      !> RK_STEP = false unless in multibox mode
      if (RK_step) call mb_communicate(g)

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
            if (RK_step) call mb_communicate(g1)
            call solve_gke(g1, g2, restart_time_step, istep)
         case (3)
            g2 = g1 + g2
            if (RK_step) call mb_communicate(g2)
            call solve_gke(g2, g, restart_time_step, istep)
         end select
         if (restart_time_step) then
            ! If the code_dt is reset, we need to quit this loop and restart the timestep again
            icnt = 10
         else
            icnt = icnt + 1
         end if
      end do

      !> this is gbar at intermediate time level
      g = g0 / 3.+0.5 * g1 + (g2 + g) / 6.

   end subroutine advance_explicit_rk3

   !> standard RK4
   subroutine advance_explicit_rk4(g, restart_time_step, istep)

      use dist_fn_arrays, only: g0, g1, g2, g3
      use zgrid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use multibox, only: RK_step

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      integer :: icnt

      !> RK_step is false unless in multibox mode
      if (RK_step) call mb_communicate(g)

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
            if (RK_step) call mb_communicate(g3)
            call solve_gke(g3, g2, restart_time_step, istep)
            g1 = g1 + 2.*g2
         case (3)
            ! g2 is h*k2
            g2 = g0 + 0.5 * g2
            if (RK_step) call mb_communicate(g2)
            call solve_gke(g2, g3, restart_time_step, istep)
            g1 = g1 + 2.*g3
         case (4)
            ! g3 is h*k3
            g3 = g0 + g3
            if (RK_step) call mb_communicate(g3)
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

      !> this is gbar at intermediate time level
      g = g0 + g1 / 6.

   end subroutine advance_explicit_rk4

   !> solve_gke accepts as argument gin, the guiding centre distribution function in k-space,
   !> and returns rhs_ky, the right-hand side of the gyrokinetic equation in k-space;
   !> i.e., if dg/dt = r, then rhs_ky = r*dt
   subroutine solve_gke(gin, rhs_ky, restart_time_step, istep)

      use job_manage, only: time_message
      use fields_arrays, only: phi, apar
      use stella_layouts, only: vmu_lo
      use stella_transforms, only: transform_y2ky
      use redistribute, only: gather, scatter
      use physics_flags, only: include_parallel_nonlinearity
      use physics_flags, only: include_parallel_streaming
      use physics_flags, only: include_mirror
      use physics_flags, only: nonlinear
      use physics_flags, only: full_flux_surface, radial_variation
      use physics_parameters, only: g_exb
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: ikx_max, ny, naky_all
      use kt_grids, only: swap_kxky_back
      use run_parameters, only: stream_implicit, mirror_implicit, drifts_implicit
      use dissipation, only: include_collisions, advance_collisions_explicit, collisions_implicit
      use sources, only: source_option_switch, source_option_krook
      use sources, only: add_krook_operator
      use parallel_streaming, only: advance_parallel_streaming_explicit
      use fields, only: advance_fields, fields_updated, get_radial_correction
      use mirror_terms, only: advance_mirror_explicit
      use flow_shear, only: advance_parallel_flow_shear
      use multibox, only: include_multibox_krook, add_multibox_krook

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
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

      !> start with gbar in k-space and (ky,kx,z) local
      !> obtain fields corresponding to gbar
      if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_fields'
      call advance_fields(gin, phi, apar, dist='gbar')

      if (radial_variation) call get_radial_correction(gin, phi, dist='gbar')

      !> default is to continue with same time step size.
      !> if estimated CFL condition for nonlinear terms is violated
      !> then restart_time_step will be set to .true.
      restart_time_step = .false.
      !> calculate and add ExB nonlinearity to RHS of GK eqn
      !> do this first, as the CFL condition may require a change in time step
      !> and thus recomputation of mirror, wdrift, wstar, and parstream
      if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_ExB_nonlinearity'
      if (nonlinear) call advance_ExB_nonlinearity(gin, rhs, restart_time_step, istep)

      !> include contribution from the parallel nonlinearity (aka turbulent acceleration)
      if (include_parallel_nonlinearity .and. .not. restart_time_step) &
         call advance_parallel_nonlinearity(gin, rhs, restart_time_step)

      if (.not. restart_time_step) then

         !> include contribution from perp flow shear in the parallel component of the toroidal flow
         if ((g_exb**2) > epsilon(0.0)) call advance_parallel_flow_shear(rhs)

         !> calculate and add mirror term to RHS of GK eqn
         if (include_mirror .and. .not. mirror_implicit) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_mirror_explicit'
            call advance_mirror_explicit(gin, rhs)
         end if

         if (.not. drifts_implicit) then
            !> calculate and add alpha-component of magnetic drift term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wdrifty_explicit'
            call advance_wdrifty_explicit(gin, phi, rhs)

            !> calculate and add psi-component of magnetic drift term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wdriftx_explicit'
            call advance_wdriftx_explicit(gin, phi, rhs)

            !> calculate and add omega_* term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wstar_explicit'
            call advance_wstar_explicit(phi, rhs)
         end if

         !> calculate and add contribution from collisions to RHS of GK eqn
         if (include_collisions .and. .not. collisions_implicit) call advance_collisions_explicit(gin, phi, rhs)

         !> calculate and add parallel streaming term to RHS of GK eqn
         if (include_parallel_streaming .and. (.not. stream_implicit)) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_parallel_streaming_explicit'
            call advance_parallel_streaming_explicit(gin, phi, rhs)
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
            end do
            deallocate (rhs_ky_swap)
         end if

         if (radial_variation) call advance_radial_variation(gin, rhs)

         if (source_option_switch == source_option_krook) call add_krook_operator(gin, rhs)

         if (include_multibox_krook) call add_multibox_krook(gin, rhs)

      end if

      fields_updated = .false.

      if (allocated(rhs_y)) deallocate (rhs_y)
      nullify (rhs)

   end subroutine solve_gke

   subroutine advance_wstar_explicit(phi, gout)

      use mp, only: proc0, mp_abort
      use job_manage, only: time_message
      use fields, only: get_dchidy
      use fields_arrays, only: apar
      use stella_layouts, only: vmu_lo
      use stella_transforms, only: transform_ky2y
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, naky_all, nakx, ikx_max, ny
      use kt_grids, only: swap_kxky
      use physics_flags, only: full_flux_surface
      use dist_fn_arrays, only: wstar

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      complex, dimension(:, :, :, :, :), allocatable :: g0, g0y
      complex, dimension(:, :), allocatable :: g0_swap

      integer :: iz, it, ivmu

      !> start timing the time advance due to the driving gradients
      if (proc0) call time_message(.false., time_gke(:, 6), ' wstar advance')

      allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      if (debug) write (*, *) 'time_advance::solve_gke::get_dchidy'
      !> get d<chi>/dy in k-space
      call get_dchidy(phi, apar, g0)

      if (full_flux_surface) then
         !> assume only a single flux surface simulated
         it = 1
         allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (g0_swap(naky_all, ikx_max))
         !> transform d<chi>/dy from k-space to y-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do iz = -nzgrid, nzgrid
               call swap_kxky(g0(:, :, iz, it, ivmu), g0_swap)
               call transform_ky2y(g0_swap, g0y(:, :, iz, it, ivmu))
            end do
         end do
         !> multiply d<chi>/dy with omega_* coefficient and add to source (RHS of GK eqn)
         !       call add_wstar_term_ffs (g0y, gout)
         call add_explicit_term_ffs(g0y, wstar, gout)
         deallocate (g0y, g0_swap)
      else
         !> omega_* stays in ky,kx,z space with ky,kx,z local
         !> multiply d<chi>/dy with omega_* coefficient and add to source (RHS of GK eqn)
         if (debug) write (*, *) 'time_advance::solve_gke::add_wstar_term'
         !       call add_wstar_term (g0, gout)
         call add_explicit_term(g0, wstar(1, :, :), gout)
      end if
      deallocate (g0)

      !> stop timing the time advance due to the driving gradients
      if (proc0) call time_message(.false., time_gke(:, 6), ' wstar advance')

   end subroutine advance_wstar_explicit

   !> advance_wdrifty_explicit subroutine calculates and adds the y-component of the
   !> magnetic drift term to the RHS of the GK equation
   subroutine advance_wdrifty_explicit(g, phi, gout)

      use mp, only: proc0
      use stella_layouts, only: vmu_lo
      use job_manage, only: time_message
      use stella_transforms, only: transform_ky2y
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, ikx_max, naky, naky_all, ny
      use kt_grids, only: swap_kxky
      use physics_flags, only: full_flux_surface
      use gyro_averages, only: gyro_average
      use dist_fn_arrays, only: wdrifty_g, wdrifty_phi

      use gyro_averages, only: j0_ffs
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      integer :: ivmu, iz, it
      complex, dimension(:, :, :, :), allocatable :: dphidy
      complex, dimension(:, :, :, :, :), allocatable :: g0k, g0y
      complex, dimension(:, :), allocatable :: g0k_swap

      !> start the timing of the y component of the magnetic drift advance
      if (proc0) call time_message(.false., time_gke(:, 4), ' dgdy advance')

      allocate (dphidy(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g0k(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wdrifty_explicit::get_dgdy'
      !> calculate dg/dy in (ky,kx) space
      call get_dgdy(g, g0k)
      !> calculate dphi/dy in (ky,kx) space
      call get_dgdy(phi, dphidy)

      if (full_flux_surface) then
         !> assume only a single flux surface simulated
         it = 1
         allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (g0k_swap(naky_all, ikx_max))
         !> transform dg/dy from k-space to y-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do iz = -nzgrid, nzgrid
               call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
               call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
            end do
         end do
         !> add vM . grad y dg/dy term to equation
         call add_explicit_term_ffs(g0y, wdrifty_g, gout)

         !> get <dphi/dy> in k-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            call gyro_average(dphidy, g0k(:, :, :, :, ivmu), j0_ffs(:, :, :, ivmu))
         end do

         !> transform d<phi>/dy from k-space to y-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do iz = -nzgrid, nzgrid
               call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
               call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
            end do
         end do

         !> add vM . grad y d<phi>/dy term to equation
         call add_explicit_term_ffs(g0y, wdrifty_phi, gout)

         deallocate (g0y, g0k_swap)
      else
         if (debug) write (*, *) 'time_advance::solve_gke::add_dgdy_term'
         ! add vM . grad y dg/dy term to equation
         call add_explicit_term(g0k, wdrifty_g(1, :, :), gout)

         ! get <dphi/dy> in k-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            call gyro_average(dphidy, ivmu, g0k(:, :, :, :, ivmu))
         end do

         ! add vM . grad y d<phi>/dy term to equation
         call add_explicit_term(g0k, wdrifty_phi(1, :, :), gout)
      end if
      deallocate (g0k, dphidy)

      !> stop the timing of the y component of the magnetic drift advance
      if (proc0) call time_message(.false., time_gke(:, 4), ' dgdy advance')

   end subroutine advance_wdrifty_explicit

   !> advance_wdriftx_explicit subroutine calculates and adds the x-component of the
   !> magnetic drift term to the RHS of the GK equation
   subroutine advance_wdriftx_explicit(g, phi, gout)

      use mp, only: proc0
      use stella_layouts, only: vmu_lo
      use job_manage, only: time_message
      use stella_transforms, only: transform_ky2y
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, ikx_max, naky, naky_all, ny, akx
      use kt_grids, only: swap_kxky
      use physics_flags, only: full_flux_surface
      use gyro_averages, only: gyro_average
      use dist_fn_arrays, only: wdriftx_g, wdriftx_phi

      use gyro_averages, only: j0_ffs
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      integer :: ivmu, iz, it
      complex, dimension(:, :, :, :), allocatable :: dphidx
      complex, dimension(:, :, :, :, :), allocatable :: g0k, g0y
      complex, dimension(:, :), allocatable :: g0k_swap

      !> start the timing of the x component of the magnetic drift advance
      if (proc0) call time_message(.false., time_gke(:, 5), ' dgdx advance')

      !> do not calculate if wdriftx terms are all zero
      if (maxval(abs(akx)) < epsilon(0.)) then
         if (proc0) call time_message(.false., time_gke(:, 5), ' dgdx advance')
         return
      end if

      allocate (dphidx(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g0k(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      if (debug) write (*, *) 'time_advance::solve_gke::get_dgdx'
      !> calculate dg/dx in (ky,kx) space
      call get_dgdx(g, g0k)
      !> calculate dphi/dx in (ky,kx) space
      call get_dgdx(phi, dphidx)

      if (full_flux_surface) then
         !> assume a single flux surface is simulated
         it = 1
         allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (g0k_swap(naky_all, ikx_max))
         !> transform dg/dx from k-space to y-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do iz = -nzgrid, nzgrid
               call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
               call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
            end do
         end do
         !> add vM . grad x dg/dx term to equation
         call add_explicit_term_ffs(g0y, wdriftx_g, gout)
         !> get <dphi/dx> in k-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            call gyro_average(dphidx, g0k(:, :, :, :, ivmu), j0_ffs(:, :, :, ivmu))
         end do
         !> transform d<phi>/dx from k-space to y-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do iz = -nzgrid, nzgrid
               call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
               call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
            end do
         end do
         !> add vM . grad x d<phi>/dx term to equation
         call add_explicit_term_ffs(g0y, wdriftx_phi, gout)
         deallocate (g0y, g0k_swap)
      else
         if (debug) write (*, *) 'time_advance::solve_gke::add_dgdx_term'
         !> add vM . grad x dg/dx term to equation
         call add_explicit_term(g0k, wdriftx_g(1, :, :), gout)
         !> get <dphi/dx> in k-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            call gyro_average(dphidx, ivmu, g0k(:, :, :, :, ivmu))
         end do
         !> add vM . grad x d<phi>/dx term to equation
         call add_explicit_term(g0k, wdriftx_phi(1, :, :), gout)
      end if
      deallocate (g0k, dphidx)

      !> stop the timing of the x component of the magnetic drift advance
      if (proc0) call time_message(.false., time_gke(:, 5), ' dgdx advance')

   end subroutine advance_wdriftx_explicit

   subroutine advance_ExB_nonlinearity(g, gout, restart_time_step, istep)

      use mp, only: proc0, min_allreduce
      use mp, only: scope, allprocs, subprocs
      use stella_layouts, only: vmu_lo, imu_idx, is_idx
      use job_manage, only: time_message
      use gyro_averages, only: gyro_average
      use fields, only: get_dchidx, get_dchidy
      use fields_arrays, only: phi, apar, shift_state
      use fields_arrays, only: phi_corr_QN, phi_corr_GA
!   use fields_arrays, only: apar_corr_QN, apar_corr_GA
      use stella_transforms, only: transform_y2ky, transform_x2kx
      use stella_transforms, only: transform_y2ky_xfirst, transform_x2kx_xfirst
      use stella_time, only: cfl_dt_ExB, cfl_dt_linear, code_dt, code_dt_max
      use run_parameters, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower, fphi
      use physics_parameters, only: g_exb, g_exbfac
      use zgrid, only: nzgrid, ntubes
      use stella_geometry, only: exb_nonlin_fac, exb_nonlin_fac_p, gfac
      use kt_grids, only: nakx, ikx_max, naky, naky_all, nx, ny
      use kt_grids, only: akx, aky, rho_clamped
      use physics_flags, only: full_flux_surface, radial_variation
      use physics_flags, only: prp_shear_enabled, hammett_flow_shear
      use kt_grids, only: x, swap_kxky, swap_kxky_back
      use constants, only: pi, zi
      use file_utils, only: runtype_option_switch, runtype_multibox

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
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
               call get_dchidx(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), g0k)
               !> if running with equilibrium flow shear, make adjustment to
               !> the term multiplying dg/dy
               if (prp_shear_enabled .and. hammett_flow_shear) then
                  call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), g0a)
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
               !> if running with equilibrium flow shear, correct dg/dx term
               if (prp_shear_enabled .and. hammett_flow_shear) then
                  call get_dgdy(g(:, :, iz, it, ivmu), g0a)
                  g0k = g0k - g_exb * g_exbfac * spread(shift_state, 2, nakx) * g0a
               end if
               !> FFT to get dg/dx in (y,x) space
               call forward_transform(g0k, g0xy)
               !> compute d<chi>/dy in k-space
               call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), g0k)
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

   subroutine advance_parallel_nonlinearity(g, gout, restart_time_step)

      use constants, only: zi
      use mp, only: proc0, min_allreduce, mp_abort
      use mp, only: scope, allprocs, subprocs
      use stella_layouts, only: vmu_lo, xyz_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use job_manage, only: time_message
      use finite_differences, only: second_order_centered_zed
      use finite_differences, only: third_order_upwind
      use redistribute, only: gather, scatter
      use fields_arrays, only: phi, phi_corr_QN, phi_corr_GA
      use stella_transforms, only: transform_ky2y, transform_y2ky
      use stella_transforms, only: transform_kx2x, transform_x2kx
      use stella_time, only: cfl_dt_parallel, cfl_dt_linear, code_dt, code_dt_max
      use run_parameters, only: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower
      use zgrid, only: nzgrid, delzed, ntubes
      use extended_zgrid, only: neigen, nsegments, ikxmod
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: periodic
      use physics_flags, only: full_flux_surface, radial_variation
      use kt_grids, only: akx, aky, nakx, naky, nx, ny, ikx_max
      use kt_grids, only: swap_kxky, swap_kxky_back, rho_clamped
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: dvpa, vpa, mu
      use gyro_averages, only: gyro_average
      use parallel_streaming, only: stream_sign
      use dist_redistribute, only: xyz2vmu
      use file_utils, only: runtype_option_switch, runtype_multibox
      use extended_zgrid, only: fill_zed_ghost_zones

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
      logical, intent(out) :: restart_time_step

      integer :: ivmu, ixyz
      integer :: iz, it, iv, imu, is
      integer :: iky, ie, iseg
      integer :: advect_sign
      real :: cfl_dt
      real, dimension(:), allocatable :: dgdv
      real, dimension(:, :, :, :, :), allocatable :: g0xy
      real, dimension(:, :, :), allocatable :: gxy_vmulocal
      real, dimension(:, :), allocatable :: g1xy, advect_speed
      complex, dimension(2) :: gleft, gright
      complex, dimension(:, :, :, :), allocatable :: phi_gyro, dphidz
      complex, dimension(:, :), allocatable :: g0k, g0kxy, g0k_swap
      complex, dimension(:, :), allocatable :: tmp

      ! alpha-component of magnetic drift (requires ky -> y)
      if (proc0) call time_message(.false., time_parallel_nl(:, 1), ' parallel nonlinearity advance')

      ! Initialize cfl_dt_parallel
      cfl_dt_parallel = 10000000.

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

      allocate (g0k(naky, nakx))
      allocate (g0xy(ny, nx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      allocate (g0kxy(ny, ikx_max))
      if (radial_variation) allocate (g1xy(ny, nx))
      allocate (phi_gyro(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dphidz(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g0k_swap(2 * naky - 1, ikx_max))
      allocate (tmp(size(gout, 1), size(gout, 2)))

      ! get d<phi>/dz in vmu_lo
      ! we will need to transform it to real-space
      ! as its sign is needed for upwinding of dg/dvpa
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         ! construct <phi>
         dphidz = phi
         if (radial_variation) dphidz = dphidz + phi_corr_QN
         call gyro_average(dphidz, ivmu, phi_gyro)
         if (radial_variation) phi_gyro = phi_gyro + phi_corr_GA(:, :, :, :, ivmu)

         do iky = 1, naky
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  do iseg = 1, nsegments(ie, iky)
                     ! first fill in ghost zones at boundaries in g(z)
                     call fill_zed_ghost_zones(it, iseg, ie, iky, phi_gyro, gleft, gright)
                     ! now get d<phi>/dz
                     call second_order_centered_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
                                                    phi_gyro(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                                                    delzed(0), stream_sign(iv), gleft, gright, periodic(iky), &
                                                    dphidz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
                  end do
               end do
            end do
         end do

         if (radial_variation) then
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  ! use reality to swap from ky >= 0, all kx to kx >= 0 , all ky
                  call swap_kxky(dphidz(:, :, iz, it), g0k_swap)
                  ! transform in y
                  call transform_ky2y(g0k_swap, g0kxy)
                  ! transform in x
                  call transform_kx2x(g0kxy, g1xy)
                  g0xy(:, :, iz, it, ivmu) = g1xy * (par_nl_fac(iz, is) + d_par_nl_fac_dr(iz, is) * spread(rho_clamped, 1, ny))

                  g0k = zi * spread(aky, 2, nakx) * phi_gyro(:, :, iz, it)
                  call swap_kxky(g0k, g0k_swap)
                  call transform_ky2y(g0k_swap, g0kxy)
                  call transform_kx2x(g0kxy, g1xy)
                  g0xy(:, :, iz, it, ivmu) = g0xy(:, :, iz, it, ivmu) &
                                             + vpa(iv) * g1xy * (par_nl_drifty(iz) + d_par_nl_drifty_dr(iz) * spread(rho_clamped, 1, ny))

                  g0k = zi * spread(akx, 1, naky) * phi_gyro(:, :, iz, it)
                  call swap_kxky(g0k, g0k_swap)
                  call transform_ky2y(g0k_swap, g0kxy)
                  call transform_kx2x(g0kxy, g1xy)
                  g0xy(:, :, iz, it, ivmu) = g0xy(:, :, iz, it, ivmu) &
                                             + vpa(iv) * g1xy * (par_nl_driftx(iz) + d_par_nl_driftx_dr(iz) * spread(rho_clamped, 1, ny))

                  g0xy(:, :, iz, it, ivmu) = g0xy(:, :, iz, it, ivmu) &
                                             + vpa(iv) * mu(imu) * (par_nl_curv(iz, is) + d_par_nl_curv_dr(iz, is) * spread(rho_clamped, 1, ny))

               end do
            end do
         else
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g0k = dphidz(:, :, iz, it) * par_nl_fac(iz, is) + vpa(iv) * mu(imu) * par_nl_curv(iz, is) &
                        + zi * vpa(iv) * phi_gyro(:, :, iz, it) * (spread(akx, 1, naky) * par_nl_driftx(iz) &
                                                                   + spread(aky, 2, nakx) * par_nl_drifty(iz))
                  ! use reality to swap from ky >= 0, all kx to kx >= 0 , all ky
                  call swap_kxky(g0k, g0k_swap)
                  ! transform in y
                  call transform_ky2y(g0k_swap, g0kxy)
                  ! transform in x
                  call transform_kx2x(g0kxy, g0xy(:, :, iz, it, ivmu))
               end do
            end do
         end if
      end do

      ! do not need phi_gyro or dphidz  again so deallocate
      deallocate (phi_gyro, dphidz)
      deallocate (g0k)
      if (allocated(g1xy)) deallocate (g1xy)

      allocate (gxy_vmulocal(nvpa, nmu, xyz_lo%llim_proc:xyz_lo%ulim_alloc))
      allocate (advect_speed(nmu, xyz_lo%llim_proc:xyz_lo%ulim_alloc))

      ! we now have the advection velocity in vpa in (x,y) space
      ! next redistribute it so that (vpa,mu) are local
      if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
      call scatter(xyz2vmu, g0xy, gxy_vmulocal)
      if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
      ! advect_speed does not depend on vpa
      advect_speed = gxy_vmulocal(1, :, :)

      ! transform g from (kx,ky) to (x,y)
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               call swap_kxky(g(:, :, iz, it, ivmu), g0k_swap)
               ! transform in y
               call transform_ky2y(g0k_swap, g0kxy)
               ! transform in x
               call transform_kx2x(g0kxy, g0xy(:, :, iz, it, ivmu))
            end do
         end do
      end do

      ! redistribute so that (vpa,mu) local
      if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
      call scatter(xyz2vmu, g0xy, gxy_vmulocal)
      if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')

      allocate (dgdv(nvpa))

      ! we now need to form dg/dvpa and obtain product of dg/dvpa with advection speed
      do ixyz = xyz_lo%llim_proc, xyz_lo%ulim_proc
         do imu = 1, nmu
            ! advect_sign set to +/- 1 depending on sign of the parallel nonlinearity
            ! advection velocity
            ! NB: advect_sign = -1 corresponds to positive advection velocity
            advect_sign = int(sign(1.0, advect_speed(imu, ixyz)))
            call third_order_upwind(1, gxy_vmulocal(:, imu, ixyz), dvpa, advect_sign, dgdv)
            gxy_vmulocal(:, imu, ixyz) = dgdv * advect_speed(imu, ixyz)
            cfl_dt_parallel = min(cfl_dt_parallel, dvpa / abs(advect_speed(imu, ixyz)))
         end do
      end do

      ! finished with dgdv and advect_speed
      deallocate (dgdv, advect_speed)

      ! now that we have the full parallel nonlinearity in (x,y)-space
      ! need to redistribute so that (x,y) local for transforms
      if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')
      call gather(xyz2vmu, gxy_vmulocal, g0xy)
      if (proc0) call time_message(.false., time_parallel_nl(:, 2), ' parallel nonlinearity redist')

      ! finished with gxy_vmulocal
      deallocate (gxy_vmulocal)

      ! g0xy is parallel nonlinearity term with (x,y) on processor
      ! need to inverse Fourier transform
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               call transform_x2kx(g0xy(:, :, iz, it, ivmu), g0kxy)
               if (full_flux_surface) then
                  gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + code_dt * g0kxy
               else
                  call transform_y2ky(g0kxy, g0k_swap)
                  call swap_kxky_back(g0k_swap, tmp)
                  gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + code_dt * tmp
               end if
            end do
         end do
      end do
      deallocate (g0k_swap, g0kxy, g0xy)

      if (runtype_option_switch == runtype_multibox) call scope(allprocs)

      call min_allreduce(cfl_dt_parallel)

      if (runtype_option_switch == runtype_multibox) call scope(subprocs)

      !> check estimated cfl_dt to see if the time step size needs to be changed
      cfl_dt = min(cfl_dt_parallel, cfl_dt_linear)
      if (code_dt > cfl_dt * cfl_cushion_upper) then
         if (proc0) then
            write (*, *) ' '
            write (*, *) 'CHANGING TIME STEP:'
            write (*, '(A22, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
            write (*, '(A22, ES10.2E2)') "   cfl_dt_parallel:"//REPEAT(' ', 50), cfl_dt_parallel
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
            write (*, *) 'CHANGING TIME STEP:'
            write (*, '(A22, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
            write (*, '(A22, ES10.2E2)') "   cfl_dt_parallel:"//REPEAT(' ', 50), cfl_dt_parallel
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
!    else
!       gout = code_dt*gout
      end if

      if (proc0) call time_message(.false., time_parallel_nl(:, 1), ' parallel nonlinearity advance')

   end subroutine advance_parallel_nonlinearity

   subroutine advance_radial_variation(g, gout)

      use mp, only: mp_abort, proc0
      use job_manage, only: time_message
      use fields, only: get_dchidy
      use fields_arrays, only: phi, apar
      use fields_arrays, only: phi_corr_QN, phi_corr_GA
!   use fields_arrays, only: apar_corr_QN, apar_corr_GA
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, naky, multiply_by_rho
      use gyro_averages, only: gyro_average, gyro_average_j1
      use run_parameters, only: fphi
      use physics_flags, only: full_flux_surface
      use physics_flags, only: include_parallel_streaming, include_mirror
      use dist_fn_arrays, only: wdriftx_phi, wdrifty_phi
      use dist_fn_arrays, only: wdriftpx_g, wdriftpy_g
      use dist_fn_arrays, only: wdriftpx_phi, wdriftpy_phi !, adiabatic_phi
      use dist_fn_arrays, only: wstar, wstarp
      use mirror_terms, only: add_mirror_radial_variation
      use flow_shear, only: prl_shear, prl_shear_p
      use parallel_streaming, only: add_parallel_streaming_radial_variation

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
               call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), g0a)
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

   !> compute dg/dy in k-space
   !> accepts g(ky,kx)
   subroutine get_dgdy_2d(g, dgdy)

      use constants, only: zi
      use kt_grids, only: nakx, aky

      implicit none

      complex, dimension(:, :), intent(in) :: g
      complex, dimension(:, :), intent(out) :: dgdy

      dgdy = zi * spread(aky, 2, nakx) * g

   end subroutine get_dgdy_2d

   !> compute dg/dy in k-space
   !> accepts g(ky,kx,z,tube)
   subroutine get_dgdy_3d(g, dgdy)

      use constants, only: zi
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, aky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdy

      integer :: it, iz, ikx

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               dgdy(:, ikx, iz, it) = zi * aky(:) * g(:, ikx, iz, it)
            end do
         end do
      end do

   end subroutine get_dgdy_3d

   !> compute dg/dy in k-space
   !> accepts g(ky,kx,z,tube,(vpa,mu,spec))
   subroutine get_dgdy_4d(g, dgdy)

      use constants, only: zi
      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, aky

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dgdy

      integer :: ivmu, ikx, iz, it

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  dgdy(:, ikx, iz, it, ivmu) = zi * aky(:) * g(:, ikx, iz, it, ivmu)
               end do
            end do
         end do
      end do

   end subroutine get_dgdy_4d

   !> compute dg/dx in k-space
   !> accepts g(ky,kx)
   subroutine get_dgdx_2d(g, dgdx)

      use constants, only: zi
      use kt_grids, only: naky, akx

      implicit none

      complex, dimension(:, :), intent(in) :: g
      complex, dimension(:, :), intent(out) :: dgdx

      dgdx = zi * spread(akx, 1, naky) * g

   end subroutine get_dgdx_2d

   !> compute dg/dx in k-space
   !> accepts g(ky,kx,z,tube)
   subroutine get_dgdx_3d(g, dgdx)

      use constants, only: zi
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: akx, nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdx

      integer :: ikx, iz, it

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               dgdx(:, ikx, iz, it) = zi * akx(ikx) * g(:, ikx, iz, it)
            end do
         end do
      end do

   end subroutine get_dgdx_3d

   !> compute dg/dx in k-space
   !> accepts g(ky,kx,z,tube,(vpa,mu,spec))
   subroutine get_dgdx_4d(g, dgdx)

      use constants, only: zi
      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: akx, nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dgdx

      integer :: ivmu, ikx, iz, it

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  dgdx(:, ikx, iz, it, ivmu) = zi * akx(ikx) * g(:, ikx, iz, it, ivmu)
               end do
            end do
         end do
      end do

   end subroutine get_dgdx_4d

   subroutine add_explicit_term(g, pre_factor, src)

      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(-nzgrid:, vmu_lo%llim_proc:), intent(in) :: pre_factor
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: src

      integer :: ivmu
      integer :: iky, ikx, iz, it

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  do iky = 1, naky
                     src(iky, ikx, iz, it, ivmu) = src(iky, ikx, iz, it, ivmu) + pre_factor(iz, ivmu) * g(iky, ikx, iz, it, ivmu)
                  end do
               end do
            end do
         end do
      end do

   end subroutine add_explicit_term

   !> add vM . grad y d<phi>/dy or vM . grad x d<phi>/dx (or equivalents with g) or omega_* * d<phi>/dy term to RHS of GK equation
   subroutine add_explicit_term_ffs(g, pre_factor, src)

      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: ikx_max, nalpha

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(in) :: pre_factor
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: src

      integer :: ivmu
      integer :: ia, ikx, iz, it

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, ikx_max
                  do ia = 1, nalpha
                     src(ia, ikx, iz, it, ivmu) = src(ia, ikx, iz, it, ivmu) + pre_factor(ia, iz, ivmu) * g(ia, ikx, iz, it, ivmu)
                  end do
               end do
            end do
         end do
      end do

   end subroutine add_explicit_term_ffs

   subroutine advance_implicit(istep, phi, apar, g)

      use mp, only: proc0
      use job_manage, only: time_message
      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid
      use dissipation, only: hyper_dissipation
      use hyper, only: advance_hyper_dissipation
      use physics_flags, only: include_parallel_streaming
      use physics_flags, only: radial_variation, full_flux_surface
      use physics_flags, only: include_mirror, prp_shear_enabled
      use run_parameters, only: stream_implicit, mirror_implicit, drifts_implicit
      use implicit_solve, only: advance_implicit_terms
      use fields, only: advance_fields, fields_updated
      use mirror_terms, only: advance_mirror_implicit
      use dissipation, only: collisions_implicit, include_collisions
      use dissipation, only: advance_collisions_implicit
      use run_parameters, only: driftkinetic_implicit
      use flow_shear, only: advance_perp_flow_shear
      use multibox, only: RK_step

      implicit none

      integer, intent(in) :: istep
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
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
      if (proc0) call time_message(.false., time_gke(:, 9), ' implicit')

      ! reverse the order of operations every time step
      ! as part of alternating direction operator splitting
      ! this is needed to ensure 2nd order accuracy in time
!    if (mod(istep,2)==0) then
      ! g^{*} (coming from explicit solve) is input
      ! get g^{**}, with g^{**}-g^{*} due to mirror term

      if (RK_step) call mb_communicate(g)

      if (mod(istep, 2) == 1 .or. .not. flip_flop) then

         if (prp_shear_enabled) then
            call advance_perp_flow_shear(g)
            fields_updated = .false.
         end if

         if (hyper_dissipation) then
!          ! for hyper-dissipation, need to be in k-alpha space
!          if (alpha_space) call transform_y2ky (gy, gk)
            call advance_hyper_dissipation(g)
            fields_updated = .false.
         end if

         if (collisions_implicit .and. include_collisions) then
            call advance_fields(g, phi, apar, dist='gbar')
            call advance_collisions_implicit(mirror_implicit, phi, apar, g)
            fields_updated = .false.
         end if

         if (mirror_implicit .and. include_mirror) then
!          if (full_flux_surface) then
!             allocate (gy(ny,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
!             if (.not.alpha_space) call transform_ky2y (g, gy)
!          else
!             g_mirror => g
!          end if
            call advance_mirror_implicit(collisions_implicit, g)
            fields_updated = .false.
         end if

         ! get updated fields corresponding to advanced g
         ! note that hyper-dissipation and mirror advances
         ! depended only on g and so did not need field update
         call advance_fields(g, phi, apar, dist='gbar')

         ! g^{**} is input
         ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
         if ((stream_implicit .or. driftkinetic_implicit) .and. include_parallel_streaming) then
            call advance_implicit_terms(g, phi, apar)
            if (radial_variation .or. full_flux_surface) fields_updated = .false.
         end if

         call advance_fields(g, phi, apar, dist='gbar')

      else

         ! get updated fields corresponding to advanced g
         ! note that hyper-dissipation and mirror advances
         ! depended only on g and so did not need field update
         call advance_fields(g, phi, apar, dist='gbar')

         ! g^{**} is input
         ! get g^{***}, with g^{***}-g^{**} due to parallel streaming term
         if ((stream_implicit .or. driftkinetic_implicit) .and. include_parallel_streaming) then
            call advance_implicit_terms(g, phi, apar)
            if (radial_variation .or. full_flux_surface) fields_updated = .false.
         end if

         if (mirror_implicit .and. include_mirror) then
            call advance_mirror_implicit(collisions_implicit, g)
            fields_updated = .false.
         end if

         if (collisions_implicit .and. include_collisions) then
            call advance_fields(g, phi, apar, dist='gbar')
            call advance_collisions_implicit(mirror_implicit, phi, apar, g)
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

   subroutine mb_communicate(g_in)

      use mp, only: job
      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid
      use multibox, only: multibox_communicate, use_dirichlet_bc, apply_radial_boundary_conditions
      use fields, only: fields_updated, advance_fields
      use fields_arrays, only: phi, apar
      use file_utils, only: runtype_option_switch, runtype_multibox

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g_in

      if (runtype_option_switch == runtype_multibox) then
         if (job /= 1) then
            call advance_fields(g_in, phi, apar, dist='gbar')
         end if

         call multibox_communicate(g_in)

         if (job == 1) then
            fields_updated = .false.
            call advance_fields(g_in, phi, apar, dist='gbar')
         end if
      else if (use_dirichlet_BC) then
         call apply_radial_boundary_conditions(g_in)
         fields_updated = .false.
         call advance_fields(g_in, phi, apar, dist='gbar')
      end if

   end subroutine mb_communicate

   subroutine checksum_field(field, total)

      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky
      use extended_zgrid, only: neigen, nsegments, ikxmod
      use extended_zgrid, only: iz_low, iz_up

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: field
      real, intent(out) :: total

      integer :: it, iky, ie, iseg
      integer :: ikx

      total = 0.

      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               iseg = 1
               ikx = ikxmod(iseg, ie, iky)
               total = total + sum(cabs(field(iky, ikx, iz_low(iseg):iz_up(iseg), it)))
               if (nsegments(ie, iky) > 1) then
                  do iseg = 2, nsegments(ie, iky)
                     ikx = ikxmod(iseg, ie, iky)
                     total = total + sum(cabs(field(iky, ikx, iz_low(iseg) + 1:iz_up(iseg), it)))
                  end do
               end if
            end do
         end do
      end do

   end subroutine checksum_field

   subroutine checksum_dist(dist, total, norm)

      use mp, only: sum_allreduce
      use zgrid, only: nzgrid, ntubes
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use kt_grids, only: naky, nakx
      use vpamu_grids, only: maxwell_vpa, maxwell_mu

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: dist
      real, intent(out) :: total
      logical, intent(in), optional :: norm

      integer :: ivmu, iv, imu, is
      integer :: iky, ikx, it
      real :: subtotal

      complex, dimension(:, :, :, :), allocatable :: dist_single

      total = 0.

      allocate (dist_single(naky, nakx, -nzgrid:nzgrid, ntubes))
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         dist_single = dist(:, :, :, :, ivmu)
         if (present(norm)) then
            if (norm) then
               iv = iv_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               is = is_idx(vmu_lo, ivmu)
               do it = 1, ntubes
                  do ikx = 1, nakx
                     do iky = 1, naky
                        dist_single(iky, ikx, :, it) = dist_single(iky, ikx, :, it) * maxwell_vpa(iv, is) * maxwell_mu(1, :, imu, is)
                     end do
                  end do
               end do
            else
            end if
         end if
         call checksum(dist_single, subtotal)
         total = total + subtotal
      end do
      deallocate (dist_single)

      call sum_allreduce(total)

   end subroutine checksum_dist

   subroutine finish_time_advance

      use stella_transforms, only: finish_transforms
      use physics_flags, only: full_flux_surface
      use extended_zgrid, only: finish_extended_zgrid
      use parallel_streaming, only: finish_parallel_streaming
      use mirror_terms, only: finish_mirror
      use flow_shear, only: finish_flow_shear
      use neoclassical_terms, only: finish_neoclassical_terms
      use dissipation, only: finish_dissipation

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
      use dist_fn_arrays, only: wdriftpx_g, wdriftpy_g
      use dist_fn_arrays, only: wdriftpx_phi, wdriftpy_phi
!   use dist_fn_arrays, only: adiabatic_phi

      implicit none

      if (allocated(wdriftx_g)) deallocate (wdriftx_g)
      if (allocated(wdrifty_g)) deallocate (wdrifty_g)
      if (allocated(wdriftx_phi)) deallocate (wdriftx_phi)
      if (allocated(wdrifty_phi)) deallocate (wdrifty_phi)
      if (allocated(wdriftpx_g)) deallocate (wdriftpx_g)
      if (allocated(wdriftpy_g)) deallocate (wdriftpy_g)
      if (allocated(wdriftpx_phi)) deallocate (wdriftpx_phi)
      if (allocated(wdriftpy_phi)) deallocate (wdriftpy_phi)
!   if (allocated(adiabatic_phi)) deallocate (adiabatic_phi)

      wdriftinit = .false.

   end subroutine finish_wdrift

   subroutine finish_wstar

      use dist_fn_arrays, only: wstar, wstarp

      implicit none

      if (allocated(wstar)) deallocate (wstar)
      if (allocated(wstarp)) deallocate (wstarp)

      wstarinit = .false.

   end subroutine finish_wstar

   subroutine deallocate_arrays

      use dist_fn_arrays, only: g0, g1, g2, g3

      implicit none

      if (allocated(g0)) deallocate (g0)
      if (allocated(g1)) deallocate (g1)
      if (allocated(g2)) deallocate (g2)
      if (allocated(g3)) deallocate (g3)

   end subroutine deallocate_arrays

end module time_advance
