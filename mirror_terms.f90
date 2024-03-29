module mirror_terms

   implicit none

   public :: mirror_initialized
   public :: init_mirror, finish_mirror
   public :: mirror
   public :: advance_mirror_explicit, advance_mirror_implicit
   public :: add_mirror_radial_variation
   public :: time_mirror

   private

!  interface checksum
!     module procedure checksum_field
!     module procedure checksum_dist
!  end interface

   logical :: mirror_initialized = .false.
   real, dimension(2, 2) :: time_mirror = 0.

   integer, dimension(:, :), allocatable :: mirror_sign
   real, dimension(:, :, :, :), allocatable :: mirror
   real, dimension(:, :, :, :), allocatable :: mirror_rad_var
   real, dimension(:, :, :), allocatable :: mirror_tri_a, mirror_tri_b, mirror_tri_c
   real, dimension(:, :, :), allocatable :: mirror_int_fac
   real, dimension(:, :, :, :), allocatable :: mirror_interp_loc
   integer, dimension(:, :, :, :), allocatable :: mirror_interp_idx_shift

contains

   subroutine init_mirror

      use stella_time, only: code_dt
      use species, only: spec, nspec
      use vpamu_grids, only: nmu
      use vpamu_grids, only: mu
      use zgrid, only: nzgrid, nztot
      use kt_grids, only: nalpha
      use stella_geometry, only: dbdzed, b_dot_grad_z, gfac
      use stella_geometry, only: d2Bdrdth, dgradpardrho
      use neoclassical_terms, only: include_neoclassical_terms
      use neoclassical_terms, only: dphineo_dzed
      use run_parameters, only: mirror_implicit, mirror_semi_lagrange
      use physics_flags, only: include_mirror, radial_variation

      implicit none

      integer :: iz, ia, imu
      real, dimension(:, :), allocatable :: neoclassical_term

      if (mirror_initialized) return
      mirror_initialized = .true.

      if (.not. allocated(mirror)) allocate (mirror(nalpha, -nzgrid:nzgrid, nmu, nspec)); mirror = 0.
      if (.not. allocated(mirror_sign)) allocate (mirror_sign(nalpha, -nzgrid:nzgrid)); mirror_sign = 0

      allocate (neoclassical_term(-nzgrid:nzgrid, nspec))
      if (include_neoclassical_terms) then
         neoclassical_term = spread(dphineo_dzed(1, :), 2, nspec) * spread(spec%zt_psi0, 1, nztot) * 0.5
      else
         neoclassical_term = 0.
      end if

      !> mirror has sign consistent with being on RHS of GKE;
      !> it is the factor multiplying dg/dvpa in the mirror term
      if (include_mirror) then
         do imu = 1, nmu
            do ia = 1, nalpha
               do iz = -nzgrid, nzgrid
                  mirror(ia, iz, imu, :) = code_dt * spec%stm_psi0 * b_dot_grad_z(ia, iz) &
                                           * (mu(imu) * dbdzed(ia, iz) + neoclassical_term(iz, :))
               end do
            end do
         end do
      else
         mirror = 0.
      end if

      deallocate (neoclassical_term)

      if (radial_variation) then
         if (.not. allocated(mirror_rad_var)) then
            allocate (mirror_rad_var(nalpha, -nzgrid:nzgrid, nmu, nspec)); 
            mirror_rad_var = 0.
         end if
         !FLAG should include neoclassical corrections here?
         do imu = 1, nmu
            do ia = 1, nalpha
               do iz = -nzgrid, nzgrid
                  mirror_rad_var(ia, iz, imu, :) = code_dt * spec%stm_psi0 * mu(imu) * gfac &
                                                   * (dgradpardrho(iz) * dbdzed(ia, iz) &
                                                      + b_dot_grad_z(ia, iz) * d2Bdrdth(iz))
               end do
            end do
         end do
      end if

      do ia = 1, nalpha
         !> mirror_sign set to +/- 1 depending on the sign of the mirror term.
         !> NB: mirror_sign = -1 corresponds to positive advection velocity
         do iz = -nzgrid, nzgrid
            mirror_sign(ia, iz) = int(sign(1.0, mirror(ia, iz, 1, 1)))
         end do
      end do

      if (mirror_implicit) then
         if (mirror_semi_lagrange) then
            call init_mirror_semi_lagrange
         else
            !> set up the tridiagonal matrix that must be inverted
            !> for the implicit treatment of the mirror operator
            call init_invert_mirror_operator
         end if
      end if

   end subroutine init_mirror

   subroutine init_mirror_semi_lagrange

      use zgrid, only: nzgrid
      use vpamu_grids, only: nmu, dvpa
      use species, only: nspec
      use kt_grids, only: nalpha

      implicit none

      if (.not. allocated(mirror_interp_idx_shift)) &
         allocate (mirror_interp_idx_shift(nalpha, -nzgrid:nzgrid, nmu, nspec))
      if (.not. allocated(mirror_interp_loc)) &
         allocate (mirror_interp_loc(nalpha, -nzgrid:nzgrid, nmu, nspec))

      mirror_interp_idx_shift = int(mirror / dvpa)
      mirror_interp_loc = abs(mod(mirror, dvpa)) / dvpa

      ! f at shifted vpa
      ! is f(iv+idx_shift)*(1-mirror_interp_loc)
      ! + f(iv+idx_shift + mirror_sign)*mirror_interp_loc

   end subroutine init_mirror_semi_lagrange

   subroutine init_invert_mirror_operator

      use mp, only: mp_abort
      use stella_layouts, only: kxkyz_lo, kxyz_lo, vmu_lo
      use stella_layouts, only: iz_idx, is_idx, imu_idx, iv_idx, iy_idx
      use zgrid, only: nzgrid
      use vpamu_grids, only: dvpa, vpa, mu
      use vpamu_grids, only: nvpa, nmu
      use physics_flags, only: full_flux_surface
      use species, only: spec
      use kt_grids, only: nalpha
      use stella_geometry, only: dbdzed
      use neoclassical_terms, only: include_neoclassical_terms
      use neoclassical_terms, only: dphineo_dzed
      use run_parameters, only: vpa_upwind, time_upwind
      use run_parameters, only: maxwellian_normalization

      implicit none

      integer :: sgn
      integer :: iy, iz, is, imu, iv
      integer :: ivmu, ikxkyz, ikxyz
      integer :: llim, ulim
      real :: tupwndfac, zero
      real, dimension(:, :), allocatable :: a, b, c

      zero = 100.*epsilon(0.)

      !> mirror_int_fac = exp(vpa^2 * (mu*dB/dz)/(mu*dB/dz + Z*e*dpihnc/dz))
      !> is the integrating factor needed to turn the dg/dvpa part of the GKE advance
      !> into an advection equation
      if (.not. allocated(mirror_int_fac)) then
         if (include_neoclassical_terms) then
            allocate (mirror_int_fac(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               is = is_idx(vmu_lo, ivmu)
               do iy = 1, nalpha
                  where (abs(mu(imu) * dbdzed(iy, :)) > zero)
                     mirror_int_fac(iy, :, ivmu) = exp(vpa(iv)**2 * mu(imu) * dbdzed(iy, :) &
                                                       / (mu(imu) * dbdzed(iy, :) + spec(is)%zt_psi0 * dphineo_dzed(1, :) * 0.5))
                  elsewhere
                     mirror_int_fac(iy, :, ivmu) = 1.0
                  end where
               end do
            end do
         else
            ! mirror_int_fac should never be used in this case
            allocate (mirror_int_fac(1, 1, 1)); mirror_int_fac = 0.
         end if
      end if

      !> a, b and c contain the sub-, main- and super-diagonal terms, respectively
      allocate (a(nvpa, -1:1)); a = 0.
      allocate (b(nvpa, -1:1)); b = 0.
      allocate (c(nvpa, -1:1)); c = 0.

      if (.not. allocated(mirror_tri_a)) then
         !> if running in full-flux-surface mode, solve mirror advance
         !> in y-space rather than ky-space due to alpha-dependence of coefficients
         if (full_flux_surface) then
            llim = kxyz_lo%llim_proc
            ulim = kxyz_lo%ulim_alloc
         else
            llim = kxkyz_lo%llim_proc
            ulim = kxkyz_lo%ulim_alloc
         end if

         allocate (mirror_tri_a(nvpa, nmu, llim:ulim)); mirror_tri_a = 0.
         allocate (mirror_tri_b(nvpa, nmu, llim:ulim)); mirror_tri_b = 0.
         allocate (mirror_tri_c(nvpa, nmu, llim:ulim)); mirror_tri_c = 0.
      end if

      !> corresponds to sign of mirror term positive on RHS of equation
      a(2:, 1) = -0.5 * (1.0 - 2.0 * vpa_upwind) / dvpa
      b(2:, 1) = -2.0 * vpa_upwind / dvpa
      c(2:nvpa - 1, 1) = 0.5 * (1.0 + 2.0 * vpa_upwind) / dvpa
      !> must treat boundary carefully
      !> assumes fully upwinded at outgoing boundary
      b(1, 1) = -1.0 / dvpa
      c(1, 1) = 1.0 / dvpa

      !> corresponds to sign of mirror term negative on RHS of equation
      a(2:nvpa - 1, -1) = -0.5 * (1.0 + 2.0 * vpa_upwind) / dvpa
      b(:nvpa - 1, -1) = 2.0 * vpa_upwind / dvpa
      c(:nvpa - 1, -1) = 0.5 * (1.0 - 2.0 * vpa_upwind) / dvpa
      !> must treat boundary carefully
      !> assumes fully upwinded at outgoing boundary
      a(nvpa, -1) = -1.0 / dvpa
      b(nvpa, -1) = 1.0 / dvpa

      !> time_upwind = 0.0 corresponds to centered in time
      !> time_upwind = 1.0 corresponds to fully implicit (upwinded)
      tupwndfac = 0.5 * (1.0 + time_upwind)
      a = a * tupwndfac
      c = c * tupwndfac
      if (maxwellian_normalization) then
         !> account for fact that we have expanded d(gnorm)/dvpa, where gnorm = g/exp(-v^s);
         !> this gives rise to d(gnorm*exp(-vpa^2))/dvpa + 2*vpa*gnorm*exp(-vpa^2) term
         !> we solve for gnorm*exp(-vpa^2) and later multiply by exp(vpa^2) to get gnorm
         b = b + spread(2.0 * vpa, 2, 3)
      end if

      if (full_flux_surface) then
         do ikxyz = kxyz_lo%llim_proc, kxyz_lo%ulim_proc
            iy = iy_idx(kxyz_lo, ikxyz)
            iz = iz_idx(kxyz_lo, ikxyz)
            is = is_idx(kxyz_lo, ikxyz)
            sgn = mirror_sign(iy, iz)
            do imu = 1, nmu
               do iv = 1, nvpa
                  mirror_tri_a(iv, imu, ikxyz) = -a(iv, sgn) * mirror(iy, iz, imu, is)
                  mirror_tri_b(iv, imu, ikxyz) = 1.0 - b(iv, sgn) * mirror(iy, iz, imu, is) * tupwndfac
                  mirror_tri_c(iv, imu, ikxyz) = -c(iv, sgn) * mirror(iy, iz, imu, is)
               end do
            end do
         end do
      else
         !> multiply by mirror coefficient
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iy = 1
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            sgn = mirror_sign(iy, iz)
            do imu = 1, nmu
               do iv = 1, nvpa
                  mirror_tri_a(iv, imu, ikxkyz) = -a(iv, sgn) * mirror(iy, iz, imu, is)
                  mirror_tri_b(iv, imu, ikxkyz) = 1.0 - b(iv, sgn) * mirror(iy, iz, imu, is) * tupwndfac
                  mirror_tri_c(iv, imu, ikxkyz) = -c(iv, sgn) * mirror(iy, iz, imu, is)
               end do
            end do
         end do
      end if

      deallocate (a, b, c)

   end subroutine init_invert_mirror_operator

   !> advance_mirror_explicit calculates the contribution to the RHS of the gyrokinetic equation
   !> due to the mirror force term; it treats all terms explicitly in time
   subroutine advance_mirror_explicit(g, gout)

      use mp, only: proc0
      use redistribute, only: gather, scatter
      use dist_fn_arrays, only: gvmu
      use job_manage, only: time_message
      use stella_layouts, only: kxyz_lo, kxkyz_lo, vmu_lo
      use stella_layouts, only: iv_idx, is_idx
      use stella_transforms, only: transform_ky2y
      use zgrid, only: nzgrid, ntubes
      use physics_flags, only: full_flux_surface
      use kt_grids, only: nakx, naky, naky_all, ny, ikx_max
      use kt_grids, only: swap_kxky
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: vpa, maxwell_vpa
      use run_parameters, only: fields_kxkyz, maxwellian_normalization
      use dist_redistribute, only: kxkyz2vmu, kxyz2vmu

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      complex, dimension(:, :, :), allocatable :: g0v
      complex, dimension(:, :, :, :, :), allocatable :: g0x
      complex, dimension(:, :), allocatable :: dgdv, g_swap

      integer :: ikxyz, iz, it
      integer :: ivmu, iv, imu, is

      !> start the timer for this subroutine
      if (proc0) call time_message(.false., time_mirror(:, 1), ' Mirror advance')

      if (full_flux_surface) then
         !> assume we are simulating a single flux surface
         it = 1

         allocate (g0v(nvpa, nmu, kxyz_lo%llim_proc:kxyz_lo%ulim_alloc))
         allocate (g0x(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (dgdv(nvpa, nmu))
         allocate (g_swap(naky_all, ikx_max))

         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do iz = -nzgrid, nzgrid
               !> swap from ky >= 0 and all kx to kx >= 0 and all ky
               !> needed for ky2y transform below
               call swap_kxky(g(:, :, iz, it, ivmu), g_swap)
               !> for upwinding of vpa, need to evaluate dg/dvpa in y-space
               !> this is necessary because the advection speed contains dB/dz, which depends on y
               !> first must take g(ky,kx) and transform to g(y,kx)
               call transform_ky2y(g_swap, g0x(:, :, iz, it, ivmu))
            end do
         end do
         !> remap g so velocities are local
         call scatter(kxyz2vmu, g0x, g0v)
         !> next, calculate dg/dvpa;
         !> we enforce a boundary condition on <f>, but with full_flux_surface = T,
         !> g = <f> / F, so we use the chain rule to get two terms:
         !> one with exp(vpa^2)*d<f>/dvpa and another that is proportional to exp(vpa^2) * <f>/F * d ln F /dvpa
         do ikxyz = kxyz_lo%llim_proc, kxyz_lo%ulim_proc
            is = is_idx(kxyz_lo, ikxyz)
            !> remove exp(-vpa^2) normalisation from g before differentiating
            do imu = 1, nmu
               dgdv(:, imu) = g0v(:, imu, ikxyz) * maxwell_vpa(:, is)
            end do
            !> get d <f> / dvpa
            call get_dgdvpa_ffs(dgdv, ikxyz)
            do iv = 1, nvpa
               g0v(iv, :, ikxyz) = dgdv(iv, :) / maxwell_vpa(iv, is) + 2.0 * vpa(iv) * g0v(iv, :, ikxyz)
            end do
         end do
         !> then take the results and remap again so y,kx,z local.
         call gather(kxyz2vmu, g0v, g0x)

         !> finally add the mirror term to the RHS of the GK eqn
         call add_mirror_term_ffs(g0x, gout)
         deallocate (dgdv, g_swap)
      else
         allocate (g0v(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
         allocate (g0x(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

         if (.not. fields_kxkyz) then
            if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
            call scatter(kxkyz2vmu, g, gvmu)
            if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         end if
         ! get dg/dvpa and store in g0v
         g0v = gvmu
         ! remove exp(-vpa^2) normalization from pdf before differentiating
         if (maxwellian_normalization) then
            do iv = 1, nvpa
               g0v(iv, :, :) = g0v(iv, :, :) * maxwell_vpa(iv, 1)
            end do
         end if
         call get_dgdvpa_explicit(g0v)
         if (maxwellian_normalization) then
            do iv = 1, nvpa
               g0v(iv, :, :) = g0v(iv, :, :) / maxwell_vpa(iv, 1) + 2.0 * vpa(iv) * gvmu(iv, :, :)
            end do
         end if
         ! swap layouts so that (z,kx,ky) are local
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         call gather(kxkyz2vmu, g0v, g0x)
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         ! get mirror term and add to source
         call add_mirror_term(g0x, gout)
      end if
      deallocate (g0x, g0v)

      if (proc0) call time_message(.false., time_mirror(:, 1), ' Mirror advance')

   end subroutine advance_mirror_explicit

   subroutine add_mirror_radial_variation(g, gout)

      use mp, only: proc0
      use redistribute, only: gather, scatter
      use dist_fn_arrays, only: gvmu
      use job_manage, only: time_message
      use stella_layouts, only: kxkyz_lo, vmu_lo
      use stella_layouts, only: is_idx, imu_idx
      use zgrid, only: nzgrid, ntubes
      use physics_flags, only: full_flux_surface
      use vpamu_grids, only: nvpa, nmu
      use run_parameters, only: fields_kxkyz
      use dist_redistribute, only: kxkyz2vmu

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      complex, dimension(:, :, :), allocatable :: g0v

      integer :: iz, it, imu, is, ivmu, ia

      allocate (g0v(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

      !if (proc0) call time_message(.false.,time_mirror(:,1),' Mirror global variation advance')

      ia = 1

      ! the mirror term is most complicated of all when doing full flux surface
      if (full_flux_surface) then
         ! FLAG DSO - Someday one should be able to do full global
      else
         if (.not. fields_kxkyz) then
            if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
            call scatter(kxkyz2vmu, g, gvmu)
            if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         end if
         ! get dg/dvpa and store in g0v
         g0v = gvmu
         call get_dgdvpa_explicit(g0v)
         ! swap layouts so that (z,kx,ky) are local
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         call gather(kxkyz2vmu, g0v, gout)
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')

         ! get mirror term and add to source
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) &
                                             * mirror_rad_var(ia, iz, imu, is)
               end do
            end do
         end do
      end if

      deallocate (g0v)

      !if (proc0) call time_message(.false.,time_mirror(:,1),' Mirror global variation advance')

   end subroutine add_mirror_radial_variation

   subroutine get_dgdvpa_ffs(g, ikxyz)

      use finite_differences, only: third_order_upwind
      use stella_layouts, only: kxyz_lo, iz_idx, iy_idx, is_idx
      use vpamu_grids, only: nvpa, nmu, dvpa

      implicit none

      complex, dimension(:, :), intent(in out) :: g
      integer, intent(in) :: ikxyz

      integer :: imu, iz, iy, is
      complex, dimension(:), allocatable :: tmp

      allocate (tmp(nvpa))
      iz = iz_idx(kxyz_lo, ikxyz)
      iy = iy_idx(kxyz_lo, ikxyz)
      is = is_idx(kxyz_lo, ikxyz)
      do imu = 1, nmu
         ! tmp is dg/dvpa
         call third_order_upwind(1, g(:, imu), dvpa, mirror_sign(iy, iz), tmp)
         g(:, imu) = tmp
      end do
      deallocate (tmp)

   end subroutine get_dgdvpa_ffs

   subroutine get_dgdvpa_explicit(g)

      use finite_differences, only: third_order_upwind
      use stella_layouts, only: kxkyz_lo, iz_idx, is_idx
      use vpamu_grids, only: nvpa, nmu, dvpa

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in out) :: g

      integer :: ikxkyz, imu, iz, is
      complex, dimension(:), allocatable :: tmp

      allocate (tmp(nvpa))
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            call third_order_upwind(1, g(:, imu, ikxkyz), dvpa, mirror_sign(1, iz), tmp)
            g(:, imu, ikxkyz) = tmp
         end do
      end do

      deallocate (tmp)

   end subroutine get_dgdvpa_explicit

   subroutine add_mirror_term(g, src)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: imu_idx, is_idx
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: src

      integer :: imu, is, ivmu
      integer :: it, iz, ikx

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  src(:, ikx, iz, it, ivmu) = src(:, ikx, iz, it, ivmu) + mirror(1, iz, imu, is) * g(:, ikx, iz, it, ivmu)
               end do
            end do
         end do
      end do

   end subroutine add_mirror_term

   subroutine add_mirror_term_ffs(g, src)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: imu_idx, is_idx
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: ikx_max

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: src

      integer :: imu, is, ivmu
      integer :: it, iz, ikx

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, ikx_max
                  src(:, ikx, iz, it, ivmu) = src(:, ikx, iz, it, ivmu) + mirror(:, iz, imu, is) * g(:, ikx, iz, it, ivmu)
               end do
            end do
         end do
      end do

   end subroutine add_mirror_term_ffs

   ! advance mirror implicit solve dg/dt = mu/m * bhat . grad B (dg/dvpa + m*vpa/T * g)
   subroutine advance_mirror_implicit(collisions_implicit, g)

      use constants, only: zi
      use mp, only: proc0
      use job_manage, only: time_message
      use redistribute, only: gather, scatter
      use finite_differences, only: fd_variable_upwinding_vpa
      use stella_layouts, only: vmu_lo, kxyz_lo, kxkyz_lo
      use stella_layouts, only: iz_idx, is_idx, iv_idx, imu_idx
      use stella_transforms, only: transform_ky2y, transform_y2ky
      use zgrid, only: nzgrid, ntubes
      use dist_fn_arrays, only: gvmu
      use physics_flags, only: full_flux_surface
      use kt_grids, only: ny, nakx
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: dvpa, maxwell_vpa, vpa
      use neoclassical_terms, only: include_neoclassical_terms
      use run_parameters, only: vpa_upwind, time_upwind
      use run_parameters, only: mirror_semi_lagrange, maxwellian_normalization
      use dist_redistribute, only: kxkyz2vmu, kxyz2vmu

      implicit none

      logical, intent(in) :: collisions_implicit
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g

      integer :: ikxyz, ikxkyz, ivmu
      integer :: iv, imu, iz, is, ikx, it
      real :: tupwnd
      complex, dimension(:, :, :), allocatable :: g0v
      complex, dimension(:, :, :, :, :), allocatable :: g0x

      if (proc0) call time_message(.false., time_mirror(:, 1), ' Mirror advance')

      tupwnd = (1.0 - time_upwind) * 0.5

      ! FLAG -- STILL NEED TO IMPLEMENT VARIABLE TIME UPWINDING
      ! FOR FULL_FLUX_SURFACE

      ! now that we have g^{*}, need to solve
      ! g^{n+1} = g^{*} - dt*mu*bhat . grad B d((h^{n+1}+h^{*})/2)/dvpa
      ! define A_0^{-1} = dt*mu*bhat.gradB/2
      ! so that (A_0 + I)h^{n+1} = (A_0-I)h^{*}
      ! will need (I-A_0^{-1})h^{*} in Sherman-Morrison approach
      ! to invert and obtain h^{n+1}
      if (full_flux_surface) then
         ! if implicit treatment of collisions, then already have updated gvmu in kxkyz_lo
         if (.not. collisions_implicit) then
            ! get g^{*} with v-space on processor
            if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
            call scatter(kxkyz2vmu, g, gvmu)
            if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         end if

         allocate (g0v(nvpa, nmu, kxyz_lo%llim_proc:kxyz_lo%ulim_alloc))
         allocate (g0x(ny, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         ! for upwinding, need to evaluate dg^{*}/dvpa in y-space
         ! first must take g^{*}(ky) and transform to g^{*}(y)
         call transform_ky2y(g, g0x)

         write (*, *) 'WARNING: full_flux_surface not working in implicit_mirror advance!'

         ! convert g to g*(integrating factor), as this is what is being advected
         ! integrating factor = exp(m*vpa^2/2T * (mu*dB/dz) / (mu*dB/dz + Z*e*dphinc/dz))
         ! if dphinc/dz=0, simplifies to exp(m*vpa^2/2T)
         if (include_neoclassical_terms) then
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
                     do ikx = 1, nakx
                        g0x(:, ikx, iz, it, ivmu) = g0x(:, ikx, iz, it, ivmu) * mirror_int_fac(:, iz, ivmu)
                     end do
                  end do
               end do
            end do
         else
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               is = is_idx(vmu_lo, ivmu)
               g0x(:, :, :, :, ivmu) = g0x(:, :, :, :, ivmu) / maxwell_vpa(iv, is)
            end do
         end if

         ! second, remap g so velocities are local
         call scatter(kxyz2vmu, g0x, g0v)

         do ikxyz = kxyz_lo%llim_proc, kxyz_lo%ulim_proc
            do imu = 1, nmu
               call invert_mirror_operator(imu, ikxyz, g0v(:, imu, ikxyz))
            end do
         end do

         ! then take the results and remap again so y,kx,z local.
         call gather(kxyz2vmu, g0v, g0x)

         ! convert back from g*(integrating factor) to g
         ! integrating factor = exp(m*vpa^2/2T * (mu*dB/dz) / (mu*dB/dz + Z*e*dphinc/dz))
         ! if dphinc/dz=0, simplifies to exp(m*vpa^2/2T)
         if (include_neoclassical_terms) then
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
                     do ikx = 1, nakx
                        g0x(:, ikx, iz, it, ivmu) = g0x(:, ikx, iz, it, ivmu) / mirror_int_fac(:, iz, ivmu)
                     end do
                  end do
               end do
            end do
         else
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               is = is_idx(vmu_lo, ivmu)
               g0x(:, :, :, :, ivmu) = g0x(:, :, :, :, ivmu) * maxwell_vpa(iv, is)
            end do
         end if

         ! finally transform back from y to ky space
         call transform_y2ky(g0x, g)
      else
         ! if implicit treatment of collisions, then already have updated gvmu in kxkyz_lo
         if (.not. collisions_implicit) then
            ! get g^{*} with v-space on processor
            if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
            call scatter(kxkyz2vmu, g, gvmu)
            if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         end if

         allocate (g0v(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
         allocate (g0x(1, 1, 1, 1, 1))

         if (mirror_semi_lagrange) then
            call vpa_interpolation(gvmu, g0v)
         else
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               iz = iz_idx(kxkyz_lo, ikxkyz)
               is = is_idx(kxkyz_lo, ikxkyz)
               do imu = 1, nmu
                  ! remove exp(-vpa^2) normalization from pdf before differentiating
                  if (maxwellian_normalization) then
                     gvmu(:, imu, ikxkyz) = gvmu(:, imu, ikxkyz) * maxwell_vpa(:, is)
                  end if

                  ! calculate dg/dvpa
                  call fd_variable_upwinding_vpa(1, gvmu(:, imu, ikxkyz), dvpa, &
                                                 mirror_sign(1, iz), vpa_upwind, g0v(:, imu, ikxkyz))
                  ! construct RHS of GK equation for mirror advance;
                  ! i.e., (1-(1+alph)/2*dt*mu/m*b.gradB*(d/dv+m*vpa/T))*g^{n+1}
                  ! = RHS = (1+(1-alph)/2*dt*mu/m*b.gradB*(d/dv+m*vpa/T))*g^{n}
                  if (maxwellian_normalization) then
                     g0v(:, imu, ikxkyz) = gvmu(:, imu, ikxkyz) + tupwnd * mirror(1, iz, imu, is) * (g0v(:, imu, ikxkyz) &
                                                                                                     + 2.0 * vpa * gvmu(:, imu, ikxkyz))
                  else
                     g0v(:, imu, ikxkyz) = gvmu(:, imu, ikxkyz) + tupwnd * mirror(1, iz, imu, is) * g0v(:, imu, ikxkyz)
                  end if

                  ! invert_mirror_operator takes rhs of equation and
                  ! returns g^{n+1}
                  call invert_mirror_operator(imu, ikxkyz, g0v(:, imu, ikxkyz))

                  if (maxwellian_normalization) g0v(:, imu, ikxkyz) = g0v(:, imu, ikxkyz) / maxwell_vpa(:, is)
               end do
            end do
         end if

         ! then take the results and remap again so ky,kx,z local.
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         call gather(kxkyz2vmu, g0v, g)
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
      end if

      deallocate (g0x, g0v)

      if (proc0) call time_message(.false., time_mirror, ' Mirror advance')

   end subroutine advance_mirror_implicit

   subroutine vpa_interpolation(grid, interp)

      use vpamu_grids, only: nvpa, nmu
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, is_idx
      use run_parameters, only: mirror_linear_interp

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: grid
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(out) :: interp

      integer :: ikxkyz, iz, is, iv, imu
      integer :: shift, sgn, llim, ulim
      real :: fac0, fac1, fac2, fac3
      real :: tmp0, tmp1, tmp2, tmp3

      if (mirror_linear_interp) then
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            do imu = 1, nmu
               shift = mirror_interp_idx_shift(1, iz, imu, is)
               sgn = mirror_sign(1, iz)
               if (sgn > 0) then
                  llim = 1
                  ulim = nvpa - 1 - shift
               else
                  llim = nvpa
                  ulim = 2 - shift
               end if
               fac1 = 1.0 - mirror_interp_loc(1, iz, imu, is)
               fac2 = mirror_interp_loc(1, iz, imu, is)
               do iv = llim, ulim, sgn
                  interp(iv, imu, ikxkyz) = grid(iv + shift, imu, ikxkyz) * fac1 &
                                            + grid(iv + shift + sgn, imu, ikxkyz) * fac2
               end do
               ! either assume BC for g is zero beyond grid extrema
               ! or dg/dvpa is zero beyond grid extrema
               ! at boundary cell, use zero incoming BC for point just beyond boundary
               interp(ulim + sgn, imu, ikxkyz) = grid(ulim + shift + sgn, imu, ikxkyz) * fac1
               ! use zero incoming BC for cells beyond +/- nvgrid
               if (shift > 0) then
!                   interp(nvpa-shift+1:,imu,ikxkyz) = 0.

                  do iv = nvpa - shift, nvpa
                     interp(iv, imu, ikxkyz) = grid(nvpa, imu, ikxkyz) * real(nvpa - iv) / real(shift + 1)
                  end do

               else if (shift < 0) then
!                   interp(:-shift,imu,ikxkyz) = 0.

                  do iv = 1, -shift
                     interp(iv, imu, ikxkyz) = grid(1, imu, ikxkyz) * real(iv - 1) / real(-shift)
                  end do

               end if
            end do
         end do
      else
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            do imu = 1, nmu
               tmp0 = mirror_interp_loc(1, iz, imu, is)
               tmp1 = tmp0 - 2.0
               tmp2 = tmp0 - 1.0
               tmp3 = tmp0 + 1.0

               shift = mirror_interp_idx_shift(1, iz, imu, is)
               sgn = mirror_sign(1, iz)

               if (shift == 0) then
                  ! deal with boundary point near outgoing BC
                  ! using 2-point (linear) interpolation
                  ! could do 3-point to improve accuracy
                  fac1 = 1.0 - tmp0
                  fac2 = tmp0
                  if (sgn > 0) then
                     llim = 2
                     ulim = nvpa - 2 - shift
                  else
                     llim = nvpa - 1
                     ulim = 3 - shift
                  end if
                  iv = llim - sgn
                  interp(iv, imu, ikxkyz) = grid(iv + shift, imu, ikxkyz) * fac1 &
                                            + grid(iv + shift + sgn, imu, ikxkyz) * fac2
               else
                  if (sgn > 0) then
                     llim = 1
                     ulim = nvpa - 2 - shift
                  else
                     llim = nvpa
                     ulim = 3 - shift
                  end if
               end if

               ! if llim > ulim (for sgn > 0) or llim < ulim (for sgn < 0)
               ! then there are no elements to be interpolated
               if (sgn * ulim < sgn * llim) then
                  interp(:, imu, ikxkyz) = 0.
               else
                  ! coefficient multiplying g_{iv-1}
                  fac0 = -tmp0 * tmp1 * tmp2
                  ! coefficient multiplying g_{iv}
                  fac1 = 3.*tmp1 * tmp2 * tmp3
                  ! coefficient multiplying g_{iv+1}
                  fac2 = -3.*tmp0 * tmp1 * tmp3
                  ! coefficient multiplying g_{iv+2}
                  fac3 = tmp0 * tmp2 * tmp3
                  do iv = llim, ulim, sgn
                     interp(iv, imu, ikxkyz) = (grid(iv + shift - sgn, imu, ikxkyz) * fac0 &
                                                + grid(iv + shift, imu, ikxkyz) * fac1 &
                                                + grid(iv + shift + sgn, imu, ikxkyz) * fac2 &
                                                + grid(iv + shift + 2 * sgn, imu, ikxkyz) * fac3) / 6.
                  end do

                  ! at boundary cell, use zero incoming BC for point just beyond boundary
                  interp(ulim + sgn, imu, ikxkyz) = (grid(ulim + shift, imu, ikxkyz) * fac0 &
                                                     + grid(ulim + shift + sgn, imu, ikxkyz) * fac1 &
                                                     + grid(ulim + shift + 2 * sgn, imu, ikxkyz) * fac2) / 6.
                  interp(ulim + 2 * sgn, imu, ikxkyz) = (grid(ulim + shift + sgn, imu, ikxkyz) * fac0 &
                                                         + grid(ulim + shift + 2 * sgn, imu, ikxkyz) * fac1) / 6.
                  ! use zero incoming BC for cells beyond +/- nvgrid
                  if (shift > 0) then
                     interp(nvpa - shift + 1:, imu, ikxkyz) = 0.
                  else if (shift < 0) then
                     interp(:-shift, imu, ikxkyz) = 0.
                  end if
               end if
            end do
         end do
      end if

   end subroutine vpa_interpolation

   subroutine invert_mirror_operator(imu, ilo, g)

      use finite_differences, only: tridag

      implicit none

      integer, intent(in) :: imu, ilo
      complex, dimension(:), intent(in out) :: g

      call tridag(1, mirror_tri_a(:, imu, ilo), mirror_tri_b(:, imu, ilo), mirror_tri_c(:, imu, ilo), g)

   end subroutine invert_mirror_operator

   subroutine finish_mirror

      use run_parameters, only: mirror_implicit, mirror_semi_lagrange

      implicit none

      if (allocated(mirror)) deallocate (mirror)
      if (allocated(mirror_sign)) deallocate (mirror_sign)
      if (allocated(mirror_rad_var)) deallocate (mirror_rad_var)

      if (mirror_implicit) then
         if (mirror_semi_lagrange) then
            call finish_mirror_semi_lagrange
         else
            call finish_invert_mirror_operator
         end if
      end if

      mirror_initialized = .false.

   end subroutine finish_mirror

   subroutine finish_mirror_semi_lagrange

      implicit none

      if (allocated(mirror_interp_loc)) deallocate (mirror_interp_loc)
      if (allocated(mirror_interp_idx_shift)) deallocate (mirror_interp_idx_shift)

   end subroutine finish_mirror_semi_lagrange

   subroutine finish_invert_mirror_operator

      implicit none

      if (allocated(mirror_tri_a)) then
         deallocate (mirror_tri_a)
         deallocate (mirror_tri_b)
         deallocate (mirror_tri_c)
      end if

      if (allocated(mirror_int_fac)) deallocate (mirror_int_fac)

   end subroutine finish_invert_mirror_operator

   ! subroutine checksum_field (field, total)

   !   use zgrid, only: nzgrid, ntubes
   !   use kt_grids, only: naky
   !   use extended_zgrid, only: neigen, nsegments, ikxmod
   !   use extended_zgrid, only: iz_low, iz_up

   !   implicit none

   !   complex, dimension (:,:,-nzgrid:,:), intent (in) :: field
   !   real, intent (out) :: total

   !   integer :: it, iky, ie, iseg
   !   integer :: ikx

   !   total = 0.

   !   do iky = 1, naky
   !      do it = 1, ntubes
   !         do ie = 1, neigen(iky)
   !            iseg = 1
   !            ikx = ikxmod(iseg,ie,iky)
   !            total = total + sum(cabs(field(iky,ikx,iz_low(iseg):iz_up(iseg),it)))
   !            if (nsegments(ie,iky) > 1) then
   !               do iseg = 2, nsegments(ie,iky)
   !                  ikx = ikxmod(iseg,ie,iky)
   !                  total = total + sum(cabs(field(iky,ikx,iz_low(iseg)+1:iz_up(iseg),it)))
   !               end do
   !            end if
   !         end do
   !      end do
   !   end do

   ! end subroutine checksum_field

   ! subroutine checksum_dist (dist, total, norm)

   !   use mp, only: sum_allreduce
   !   use zgrid, only: nzgrid, ntubes
   !   use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
   !   use kt_grids, only: naky, nakx
   !   use vpamu_grids, only: maxwell_vpa, maxwell_mu

   !   implicit none

   !   complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: dist
   !   real, intent (out) :: total
   !   logical, intent (in), optional :: norm

   !   integer :: ivmu, iv, imu, is
   !   integer :: iky, ikx, it
   !   real :: subtotal

   !   complex, dimension (:,:,:,:), allocatable :: dist_single

   !   total = 0.

   !   allocate (dist_single(naky,nakx,-nzgrid:nzgrid,ntubes))
   !   do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
   !      dist_single = dist(:,:,:,:,ivmu)
   !      if (present(norm)) then
   !         if (norm) then
   !            iv = iv_idx(vmu_lo,ivmu)
   !            imu = imu_idx(vmu_lo,ivmu)
   !            is = is_idx(vmu_lo,ivmu)
   !            do it = 1, ntubes
   !               do ikx = 1, nakx
   !                  do iky = 1, naky
   !                     dist_single(iky,ikx,:,it) = dist_single(iky,ikx,:,it) * maxwell_vpa(iv,is) * maxwell_mu(1,:,imu,is)
   !                  end do
   !               end do
   !            end do
   !         else
   !         end if
   !      end if
   !      call checksum (dist_single, subtotal)
   !      total = total + subtotal
   !   end do
   !   deallocate (dist_single)

   !   call sum_allreduce (total)

   ! end subroutine checksum_dist

end module mirror_terms
