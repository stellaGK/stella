!###############################################################################
!################### MIRROR TERM OF THE GYROKINETIC EQUATION ###################
!###############################################################################
! 
! This module evolves the magnetic mirror term:
!     - v_{th,s} mu_s b . ∇B ∂g_{k,s} / ∂v_parallel
! 
! Define the mirror coefficient: 
!        <mirror> = v_{th,s} * mu_s * b . ∇z * dB/dz * code_dt/dvpar
! If neoclassical effects are included this is also added to <mirror> here. 
!###############################################################################
module gk_mirror

   ! Load debug flags
   use debug_flags, only: debug => mirror_terms_debug
   
   implicit none

   ! Make routines available to other modules
   public :: initialised_mirror
   public :: init_mirror, finish_mirror
   public :: mirror
   public :: advance_mirror_explicit, advance_mirror_implicit
   public :: add_mirror_radial_variation

   private

   integer, dimension(:, :), allocatable :: mirror_sign
   real, dimension(:, :, :, :), allocatable :: mirror
   real, dimension(:, :, :, :), allocatable :: mirror_rad_var
   real, dimension(:, :, :), allocatable :: mirror_tri_a, mirror_tri_b, mirror_tri_c
   real, dimension(:, :, :), allocatable :: mirror_int_fac
   real, dimension(:, :, :, :), allocatable :: mirror_interp_loc
   integer, dimension(:, :, :, :), allocatable :: mirror_interp_idx_shift
   complex, dimension(:, :, :, :), allocatable :: response_apar_denom
   
   ! Only initialise once
   logical :: initialised_mirror = .false.

contains

   !****************************************************************************
   !                              INITIALISE MIRROR
   !****************************************************************************
  subroutine init_mirror
    
      use grids_time, only: code_dt
      use grids_species, only: spec, nspec
      use grids_velocity, only: nmu
      use grids_velocity, only: mu
      use grids_z, only: nzgrid, nztot
      use grids_kxky, only: nalpha
      use geometry, only: dbdzed, b_dot_gradz, gfac
      use geometry, only: d2Bdrdth, d_bdotgradz_drho
      use neoclassical_terms, only: include_neoclassical_terms
      use neoclassical_terms, only: dphineo_dzed
      use parameters_numerical, only: mirror_implicit, mirror_semi_lagrange
      use parameters_physics, only: include_apar
      use parameters_physics, only: include_mirror, radial_variation
      use parameters_physics, only: mirrorknob

      implicit none

      integer :: iz, ia, imu
      real, dimension(:, :), allocatable :: neoclassical_term

      !-------------------------------------------------------------------------

      if (initialised_mirror) return
      initialised_mirror = .true.
 
      if (debug) write(*,*) 'no debug messages for mirror_terms.f90 yet'
      
      if (.not. allocated(mirror)) allocate (mirror(nalpha, -nzgrid:nzgrid, nmu, nspec)); mirror = 0.
      if (.not. allocated(mirror_sign)) allocate (mirror_sign(nalpha, -nzgrid:nzgrid)); mirror_sign = 0

      allocate (neoclassical_term(-nzgrid:nzgrid, nspec))
      if (include_neoclassical_terms) then
         neoclassical_term = spread(dphineo_dzed(1, :), 2, nspec) * spread(spec%zt_psi0, 1, nztot) * 0.5
      else
         neoclassical_term = 0.
      end if

      ! Mirror has sign consistent with being on RHS of GKE;
      ! it is the factor multiplying dg/dvpa in the mirror term
      if (include_mirror) then
         do imu = 1, nmu
            do ia = 1, nalpha
               do iz = -nzgrid, nzgrid
                  mirror(ia, iz, imu, :) = mirrorknob * code_dt * spec%stm_psi0 * b_dot_gradz(ia, iz) &
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
         !FLAG - should include neoclassical corrections here?
         do imu = 1, nmu
            do ia = 1, nalpha
               do iz = -nzgrid, nzgrid
                  mirror_rad_var(ia, iz, imu, :) = code_dt * spec%stm_psi0 * mu(imu) * gfac &
                                                   * (d_bdotgradz_drho(iz) * dbdzed(ia, iz) &
                                                      + b_dot_gradz(ia, iz) * d2Bdrdth(iz))
               end do
            end do
         end do
      end if

      do ia = 1, nalpha
         ! Mirror_sign set to +/- 1 depending on the sign of the mirror term.
         ! NB: mirror_sign = -1 corresponds to positive advection velocity
         do iz = -nzgrid, nzgrid
            mirror_sign(ia, iz) = int(sign(1.0, mirror(ia, iz, 1, 1)))
         end do
      end do

      if (mirror_implicit) then
         if (mirror_semi_lagrange) then
            call init_mirror_semi_lagrange
         else
            ! Set up the tridiagonal matrix that must be inverted
            ! for the implicit treatment of the mirror operator
            call init_invert_mirror_operator

            ! If advancing apar, need to get response of pdf to unit impulse in apar
            if (include_apar) call init_mirror_response
         end if

      end if

   end subroutine init_mirror

   !****************************************************************************
   !                       INITIALISE THE SEMI LAGRANGIAN ALGORITHM 
   !****************************************************************************
   ! If semi-lagrange is chosen as the mirror advance algorithm then this
   ! is the initialisaiton for the mirror terms
   !****************************************************************************
   subroutine init_mirror_semi_lagrange

      use grids_z, only: nzgrid
      use grids_velocity, only: nmu, dvpa
      use grids_species, only: nspec
      use grids_kxky, only: nalpha

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

   !****************************************************************************
   !                     INITIALISE IMPLICIT MIRROR MATRIX 
   !****************************************************************************
   ! If mirror advance option is implicit then we need to create the matrix that
   ! is to be inverted.
   ! This is a tri-diagonal matrix in vpa that is inverted. The tri-diagonal 
   ! arises from the vpar derivative. 
   ! Note that this is an option for FFS as well. 
   !****************************************************************************
   subroutine init_invert_mirror_operator

      ! Parallelisation
      use mp, only: mp_abort
      use parallelisation_layouts, only: kxkyz_lo, kxyz_lo, vmu_lo
      use parallelisation_layouts, only: iz_idx, is_idx, imu_idx, iv_idx, iy_idx
      
      ! Parameters
      use parameters_physics, only: full_flux_surface
      use parameters_numerical, only: vpa_upwind, time_upwind

      ! Grids
      use grids_z, only: nzgrid
      use grids_velocity, only: dvpa, vpa, mu
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: spec
      use grids_kxky, only: nalpha

      ! Geometry + Neoclassical
      use geometry, only: dbdzed
      use neoclassical_terms, only: include_neoclassical_terms
      use neoclassical_terms, only: dphineo_dzed

      implicit none

      integer :: sgn
      integer :: iy, iz, is, imu, iv
      integer :: ivmu, ikxkyz, ikxyz
      integer :: llim, ulim
      real :: tupwndfac, zero
      real, dimension(:, :), allocatable :: a, b, c

      !-------------------------------------------------------------------------

      zero = 100.*epsilon(0.)

      ! mirror_int_fac = exp(vpa^2 * (mu*dB/dz)/(mu*dB/dz + Z*e*dpihnc/dz))
      ! is the integrating factor needed to turn the dg/dvpa part of the GKE advance
      ! into an advection equation
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

      ! a, b and c contain the sub-, main- and super-diagonal terms, respectively
      allocate (a(nvpa, -1:1)); a = 0.
      allocate (b(nvpa, -1:1)); b = 0.
      allocate (c(nvpa, -1:1)); c = 0.

      if (.not. allocated(mirror_tri_a)) then
         ! If running in full-flux-surface mode, solve mirror advance
         ! in y-space rather than ky-space due to alpha-dependence of coefficients
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

      ! Corresponds to sign of mirror term positive on RHS of equation
      a(2:, 1) = -0.5 * (1.0 - 2.0 * vpa_upwind) / dvpa
      b(2:, 1) = -2.0 * vpa_upwind / dvpa
      c(2:nvpa - 1, 1) = 0.5 * (1.0 + 2.0 * vpa_upwind) / dvpa
      ! We must treat boundary carefully
      ! Assume fully upwinded at outgoing boundary
      b(1, 1) = -1.0 / dvpa
      c(1, 1) = 1.0 / dvpa

      ! This corresponds to sign of mirror term negative on RHS of equation
      a(2:nvpa - 1, -1) = -0.5 * (1.0 + 2.0 * vpa_upwind) / dvpa
      b(:nvpa - 1, -1) = 2.0 * vpa_upwind / dvpa
      c(:nvpa - 1, -1) = 0.5 * (1.0 - 2.0 * vpa_upwind) / dvpa
      ! We must treat boundary carefully
      ! Assume fully upwinded at outgoing boundary
      a(nvpa, -1) = -1.0 / dvpa
      b(nvpa, -1) = 1.0 / dvpa

      ! time_upwind = 0.0 corresponds to centered in time
      ! time_upwind = 1.0 corresponds to fully implicit (upwinded)
      tupwndfac = 0.5 * (1.0 + time_upwind)
      a = a * tupwndfac
      c = c * tupwndfac

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
         ! Multiply by mirror coefficient
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

   !****************************************************************************
   !                          ELECTROMAGNETIC MIRROR RESPONSE
   !****************************************************************************
   ! If electromagnetic effects with apar are included then we need to get 
   ! response of the distrubution function to a unit impulse in apar.
   ! This is needed when the mirror terms are being treated implicitly. 
   !****************************************************************************
   subroutine init_mirror_response

      use parallelisation_layouts, only: kxkyz_lo
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: nmu, nvpa
      use grids_kxky, only: naky, nakx
      use field_equations_electromagnetic, only: advance_apar

      implicit none

      complex :: apar
      integer :: ikxkyz, imu
      complex, dimension(:), allocatable :: rhs
      complex, dimension(:, :, :), allocatable :: mirror_response_g
      character(5) :: dist

      !-------------------------------------------------------------------------

      allocate (rhs(nvpa))
      allocate (mirror_response_g(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      if (.not. allocated(response_apar_denom)) then
         allocate (response_apar_denom(naky, nakx, -nzgrid:nzgrid, ntubes))
      end if

      ! Give unit impulse to apar and solve for g (stored in mirror_response_g)
      apar = 1.0
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         do imu = 1, nmu
            ! Calculate the rhs of the 'homogeneous' mirror advance equation
            call get_mirror_rhs_apar_contribution(rhs, apar, imu, ikxkyz)
            ! Invert the mirror operator and replace rhs with the solution
            call invert_mirror_operator(imu, ikxkyz, rhs)
            ! Store the solution in mirror_response_g
            mirror_response_g(:, imu, ikxkyz) = rhs
         end do
      end do

      ! Calculate the contribution to apar from mirror_response_g
      dist = 'g'
      call advance_apar(mirror_response_g, dist, response_apar_denom)

      ! response_apar_denom is the thing that the inhomogeneous contribution to apar
      ! must be divided by to obtain the total apar; i.e., apar = apar_inh + apar_h,
      ! and apar_h is proportional to apar, so (1 - apar_h / apar) * apar = apar_inh
      response_apar_denom = 1.0 - response_apar_denom

      deallocate (rhs, mirror_response_g)

   end subroutine init_mirror_response

   !****************************************************************************
   !                           EXPLICIT MIRROR ADVANCE
   !****************************************************************************
   ! advance_mirror_explicit calculates the contribution to the RHS of the 
   ! gyrokinetic equation due to the mirror force term.
   ! This treats all terms explicitly in time. This is fine for adiabatic/modified
   ! boltzmann electrons. 
   !****************************************************************************
   subroutine advance_mirror_explicit(g, gout)

      ! Parallelisation
      use mp, only: proc0
      use job_manage, only: time_message
      use timers, only: time_mirror
      use initialise_redistribute, only: kxkyz2vmu, kxyz2vmu
      use redistribute, only: gather, scatter
      use parallelisation_layouts, only: fields_kxkyz
      use parallelisation_layouts, only: iv_idx, is_idx
      use parallelisation_layouts, only: kxyz_lo, kxkyz_lo, vmu_lo
      
      ! Parameters
      use parameters_physics, only: full_flux_surface

      ! Grids + Arrays
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: vpa, maxwell_vpa
      use grids_kxky, only: nakx, naky, naky_all, ny, ikx_max
      use arrays_distribution_function, only: gvmu

      ! Calculations
      use calculations_transforms, only: transform_ky2y
      use calculations_kxky, only: swap_kxky      

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      complex, dimension(:, :, :), allocatable :: g0v
      complex, dimension(:, :, :, :, :), allocatable :: g0x
      complex, dimension(:, :), allocatable :: dgdv, g_swap
      integer :: ikxyz, iz, it
      integer :: ivmu, iv, is

      !-------------------------------------------------------------------------

      ! Start the timer for this subroutine
      if (proc0) call time_message(.false., time_mirror(:, 1), ' Mirror advance')

      if (full_flux_surface) then
         ! Assume we are simulating a single flux train
         it = 1

         allocate (g0v(nvpa, nmu, kxyz_lo%llim_proc:kxyz_lo%ulim_alloc))
         allocate (g0x(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (dgdv(nvpa, nmu))
         allocate (g_swap(naky_all, ikx_max))

         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do iz = -nzgrid, nzgrid
               ! Swap from ky >= 0 and all kx to kx >= 0 and all ky
               ! This is needed for the ky2y transform below
               call swap_kxky(g(:, :, iz, it, ivmu), g_swap)
               ! For upwinding of vpa, need to evaluate dg/dvpa in y-space.
               ! This is necessary because the advection speed contains dB/dz, which depends on y.
               ! We first must take g(ky,kx) and transform to g(y,kx)
               call transform_ky2y(g_swap, g0x(:, :, iz, it, ivmu))
            end do
         end do
         ! Remap g so velocities are local
         call scatter(kxyz2vmu, g0x, g0v)
         ! Next, calculate dg/dvpa;
         ! We enforce a boundary condition on <f>, but with full_flux_surface = T,
         ! g = <f> / F, so we use the chain rule to get two terms:
         ! One with exp(vpa^2)*d<f>/dvpa and another that is proportional to exp(vpa^2) * <f>/F * d ln F /dvpa
         do ikxyz = kxyz_lo%llim_proc, kxyz_lo%ulim_proc
            is = is_idx(kxyz_lo, ikxyz)
            ! Remove exp(-vpa^2) normalisation from g before differentiating
            dgdv(:, :) = g0v(:, :, ikxyz)
            ! Get d <f> / dvpa
            call get_dgdvpa_ffs(dgdv, ikxyz)
            g0v(:, :, ikxyz) = dgdv(:, :)
         end do
         ! Then take the results and remap again so y,kx,z local.
         call gather(kxyz2vmu, g0v, g0x)

         ! Finally add the mirror term to the RHS of the GK eqn
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

         ! Incoming gvmu is g = <f>
         ! Get dg/dvpa and store in g0v
         g0v = gvmu

         call get_dgdvpa_explicit(g0v)

         ! Swap layouts so that (z,kx,ky) are local
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         call gather(kxkyz2vmu, g0v, g0x)
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         ! Get mirror term and add to source
         call add_mirror_term(g0x, gout)
      end if
      deallocate (g0x, g0v)

      if (proc0) call time_message(.false., time_mirror(:, 1), ' Mirror advance')

   end subroutine advance_mirror_explicit

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine get_dgdvpa_ffs(g, ikxyz)

      use calculations_finite_differences, only: third_order_upwind
      use parallelisation_layouts, only: kxyz_lo, iz_idx, iy_idx, is_idx
      use grids_velocity, only: nvpa, nmu, dvpa

      implicit none

      complex, dimension(:, :), intent(in out) :: g
      integer, intent(in) :: ikxyz

      integer :: imu, iz, iy, is
      complex, dimension(:), allocatable :: tmp

      !-------------------------------------------------------------------------

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

   !****************************************************************************
   !                        ADVANCE MIRROR EXPLICIT
   !****************************************************************************
   ! Solve mirror term using an explicit algorithm: 
   !                 dg/dt = mu/m * hat{b} . grad B (dg/dvpa)
   ! The g term on the RHS is treated at the previous time step {n} and we update
   ! for g at the next intermediate time step:
   !                       dg/dt = (g^{n+1} - g{n})/delta t  
   ! RHS : mu/m * hat{b} . grad B * dg^{n}/dvpa )
   ! Note that we allow for upwinding in vpar seperately. This is an independent toggle
   subroutine get_dgdvpa_explicit(g)

      use calculations_finite_differences, only: third_order_upwind
      use parallelisation_layouts, only: kxkyz_lo, iz_idx, is_idx
      use grids_velocity, only: nvpa, nmu, dvpa

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in out) :: g

      integer :: ikxkyz, imu, iz, is
      complex, dimension(:), allocatable :: tmp

      !-------------------------------------------------------------------------

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

   !****************************************************************************
   !                            Add Mirror Term
   !****************************************************************************
   ! add_mirror_term is used to add mirror term to the RHS in explicit advance
   subroutine add_mirror_term(g, src)

      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: imu_idx, is_idx
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: src

      integer :: imu, is, ivmu
      integer :: it, iz, ikx

      !-------------------------------------------------------------------------

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

   !****************************************************************************
   !                       Add Mirror Term - Full Flux Surface
   !****************************************************************************
   subroutine add_mirror_term_ffs(g, src)

      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: imu_idx, is_idx
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: ikx_max

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

   !****************************************************************************
   !                     ADVANE MIRROR IMPLICIT 
   !****************************************************************************
   ! Advance mirror implicit solve:
   !              dg/dt = mu/m * bhat . grad B (dg/dvpa + m*vpa/T * g)
   !****************************************************************************
   subroutine advance_mirror_implicit(collisions_implicit, g, apar)

      ! Parallelisation
      use constants, only: zi
      use mp, only: proc0
      use job_manage, only: time_message
      use timers, only: time_mirror
      use initialise_redistribute, only: kxkyz2vmu, kxyz2vmu
      use redistribute, only: gather, scatter
      use parallelisation_layouts, only: iy_idx
      use parallelisation_layouts, only: vmu_lo, kxyz_lo, kxkyz_lo
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use parallelisation_layouts, only: iv_idx, imu_idx

      ! Grids + Arrays
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: ny, nakx
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: maxwell_vpa
      use grids_kxky, only: naky_all, ikx_max
      use grids_velocity, only: dvpa
      use arrays_distribution_function, only: gvmu

      ! Parameters
      use parameters_physics, only: full_flux_surface
      use parameters_physics, only: include_apar
      use parameters_numerical, only: time_upwind
      use parameters_numerical, only: vpa_upwind, time_upwind
      use parameters_numerical, only: mirror_semi_lagrange

      ! Calculations
      use calculations_tofrom_ghf, only: gbar_to_g
      use calculations_kxky, only: swap_kxky, swap_kxky_back
      use calculations_finite_differences, only: fd_variable_upwinding_vpa
      use calculations_transforms, only: transform_ky2y, transform_y2ky

      ! Fields + Neoclassical
      use field_equations, only: fields_updated
      use field_equations_electromagnetic, only: advance_apar
      use neoclassical_terms, only: include_neoclassical_terms
      
      implicit none

      logical, intent(in) :: collisions_implicit
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: apar

      integer :: ikxyz, ikxkyz, ivmu
      integer :: iky, ikx, iz, it, is
      integer :: imu
      character(5) :: dist
      real :: tupwnd
      complex, dimension(:), allocatable :: rhs
      complex, dimension(:, :, :), allocatable :: g0v
      complex, dimension(:, :, :, :, :), allocatable :: g0x

      ! Arays for FFS
      complex, dimension(:, :, :), allocatable :: dgdvpa
      integer :: iy
      complex, dimension(:, :), allocatable :: g_swap

      !-------------------------------------------------------------------------

      if (proc0) call time_message(.false., time_mirror(:, 1), ' Mirror advance')

      tupwnd = (1.0 - time_upwind) * 0.5
      
      ! Incoming pdf is g = <f>
      dist = 'g'

      ! Now that we have g^{*}, we need to solve
      !    g^{n+1} = g^{*} - dt*mu*bhat . grad B d((h^{n+1}+h^{*})/2)/dvpa
      ! Define A_0^{-1} = dt*mu*bhat.gradB/2 so that (A_0 + I)h^{n+1} = (A_0-I)h^{*}
      ! will need (I-A_0^{-1})h^{*} in Sherman-Morrison approach to invert and obtain h^{n+1}
      
      ! ------------------------------------------------------------------------
      !                                Flux Tube
      ! ------------------------------------------------------------------------
      
      ! Flux tube simulations
      if (.not. full_flux_surface) then
      
         ! For the mirror term, we need the velocity data to be local, therefore, scatter the
         ! g(naky, nakx, -nzgrid:nzgrid, ntubes, vmu-layout) data to gvmu(nvpa, nmu, kxkyz-layout) 
         if (.not. collisions_implicit) then 
            if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
            call scatter(kxkyz2vmu, g, gvmu)
            if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         end if

         allocate (g0v(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
         allocate (g0x(1, 1, 1, 1, 1))

         if (mirror_semi_lagrange) then
            call vpa_interpolation(gvmu, g0v)
         else
            allocate (rhs(nvpa))

            ! If fields are not updated, then update apar before converting from g to gbar
            ! in get_mirror_rhs_g_contribution below
            if (include_apar .and. .not. fields_updated) call advance_apar(gvmu, dist, apar)

            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               iky = iky_idx(kxkyz_lo, ikxkyz)
               ikx = ikx_idx(kxkyz_lo, ikxkyz)
               iz = iz_idx(kxkyz_lo, ikxkyz)
               it = it_idx(kxkyz_lo, ikxkyz)

               do imu = 1, nmu
                  ! calculate the contribution to the mirror advance due to source terms
                  ! involving the particle distribution function
                  call get_mirror_rhs_g_contribution(gvmu(:, imu, ikxkyz), apar(iky, ikx, iz, it), imu, ikxkyz, rhs)

                  ! invert_mirror_operator takes rhs of equation and
                  ! returns g^{n+1}
                  call invert_mirror_operator(imu, ikxkyz, rhs)
                  g0v(:, imu, ikxkyz) = rhs
               end do
            end do

            ! if not advancing apar, g0v now contains the updated pdf, g.
            ! if advanceing apar, g0v contains the 'inhomogeneous' pdf, g_{inh}
            if (include_apar) then
               ! if advancing apar, need to calculate the contribution to apar from g0v, which is the solution to the 'inhomogenous' equation
               call advance_apar(g0v, dist, apar)

               ! the total apar is the above contribution from the 'inhomogeneous' apar / response_apar_denom,
               ! with the denominator pre-calculated using a response matrix approach; in this case,
               ! the response matrix is diagonal, so it is just a scalar divide
               apar = apar / response_apar_denom

               do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                  iky = iky_idx(kxkyz_lo, ikxkyz)
                  ikx = ikx_idx(kxkyz_lo, ikxkyz)
                  iz = iz_idx(kxkyz_lo, ikxkyz)
                  it = it_idx(kxkyz_lo, ikxkyz)
                  do imu = 1, nmu
                     ! now that we have apar at the future time level, use it to solve for g;
                     ! first get the rhs of the 'homogeneous' equation, which only has
                     ! contributions from terms proportional to apar
                     call get_mirror_rhs_apar_contribution(rhs, apar(iky, ikx, iz, it), imu, ikxkyz)
                     ! invert the mirror operator to find the 'homogeneous' solution
                     call invert_mirror_operator(imu, ikxkyz, rhs)
                     ! add the 'homogeneous' solution to the 'inhomogeneous' one found above
                     ! to get g = <f> at the future time step
                     g0v(:, imu, ikxkyz) = g0v(:, imu, ikxkyz) + rhs
                  end do
               end do
            end if

            deallocate (rhs)
         end if

         ! then take the results and remap again so ky,kx,z local.
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         call gather(kxkyz2vmu, g0v, g)
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         
      ! ------------------------------------------------------------------------
      !                            Full Flux Surface                            
      ! ------------------------------------------------------------------------

      ! Full flux surface simulations
      elseif (full_flux_surface) then

         allocate (g0v(nvpa, nmu, kxyz_lo%llim_proc:kxyz_lo%ulim_alloc))
         allocate (g0x(ny, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

         ! For upwinding, need to evaluate dg^{*}/dvpa in y-space
         ! first must take g^{*}(ky) and transform to g^{*}(y)
         allocate (g_swap(naky_all, ikx_max))
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call swap_kxky(g(:, :, iz, it, ivmu), g_swap)
                  call transform_ky2y(g_swap, g0x(:, :, iz, it, ivmu))
               end do
            end do
         end do

         ! Convert g to g*(integrating factor), as this is what is being advected
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
         end if

         ! Second, remap g so velocities are local
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         call scatter(kxyz2vmu, g0x, g0v)
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')

         allocate (dgdvpa(nvpa, nmu, kxyz_lo%llim_proc:kxyz_lo%ulim_alloc))
         do ikxyz = kxyz_lo%llim_proc, kxyz_lo%ulim_proc
            iz = iz_idx(kxyz_lo, ikxyz)
            is = is_idx(kxyz_lo, ikxyz)
            iy = iy_idx(kxyz_lo, ikxyz)
            do imu = 1, nmu
               call fd_variable_upwinding_vpa(1, g0v(:, imu, ikxyz), dvpa, &
                                              mirror_sign(iy, iz), vpa_upwind, dgdvpa(:, imu, ikxyz))
               dgdvpa(:, imu, ikxyz) = g0v(:, imu, ikxyz) + tupwnd * mirror(iy, iz, imu, is) * &
                                       dgdvpa(:, imu, ikxyz)
               call invert_mirror_operator(imu, ikxyz, dgdvpa(:, imu, ikxyz))
            end do
         end do
         g0v = dgdvpa
         deallocate (dgdvpa)

         ! Then take the results and remap again so y,kx,z local.
         call gather(kxyz2vmu, g0v, g0x)

         ! Convert back from g*(integrating factor) to g
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
         end if

         ! Finally transform back from y to ky space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call transform_y2ky(g0x(:, :, iz, it, ivmu), g_swap)
                  call swap_kxky_back(g_swap, g(:, :, iz, it, ivmu))
               end do
            end do
         end do
         deallocate (g_swap)
      end if

      deallocate (g0x, g0v)

      if (proc0) call time_message(.false., time_mirror, ' Mirror advance')

   end subroutine advance_mirror_implicit

   !****************************************************************************
   !                           ADD MIRROR RHS G CONTRIBUTION
   !****************************************************************************
   subroutine get_mirror_rhs_g_contribution(g_in, apar, imu, ikxkyz, rhs)

      use parameters_physics, only: include_apar
      use parameters_numerical, only: vpa_upwind, time_upwind_minus
      use calculations_tofrom_ghf, only: gbar_to_g
      use parallelisation_layouts, only: kxkyz_lo, iz_idx, is_idx
      use calculations_finite_differences, only: fd_variable_upwinding_vpa
      use grids_velocity, only: dvpa, vpa, nvpa

      implicit none

      ! Arguments
      complex, dimension(:), intent(in) :: g_in
      complex, intent(in) :: apar
      integer, intent(in) :: imu, ikxkyz
      complex, dimension(:), intent(out) :: rhs

      ! Local variables
      integer :: iz, is
      complex, dimension(:), allocatable :: dgdv

      !-------------------------------------------------------------------------

      iz = iz_idx(kxkyz_lo, ikxkyz)
      is = is_idx(kxkyz_lo, ikxkyz)

      ! The incoming pdf is g = <f>
      rhs = g_in
      ! When advancing apar, need to compute g^{n+1} = gbar^{n} + dt * dg/dvpa * (...)
      ! the vpa derivative appearing on the RHS of the mirror equation
      ! should be operating on g, so need to have both gbar and g.
      ! NB: changes may need to be made to this call to gbar_to_g 
      if (include_apar) then
         ! RHS is converted from g to gbar
         call gbar_to_g(rhs, apar, imu, ikxkyz, -1.0)
      end if

      ! Calculate dg/dvpa
      allocate (dgdv(nvpa))
      call fd_variable_upwinding_vpa(1, g_in, dvpa, &
                                     mirror_sign(1, iz), vpa_upwind, dgdv)

      ! Construct RHS of GK equation for mirror advance;
      ! i.e., (1-(1+alph)/2*dt*mu/m*b.gradB*(d/dv+m*vpa/T))*g^{n+1}
      ! = RHS = (1+(1-alph)/2*dt*mu/m*b.gradB*(d/dv+m*vpa/T))*g^{n}

      rhs = rhs + time_upwind_minus * mirror(1, iz, imu, is) * dgdv

      deallocate (dgdv)

   end subroutine get_mirror_rhs_g_contribution

   !****************************************************************************
   !                        ADD MIRROR CONTRBUTION FROM APAR
   !****************************************************************************
   subroutine get_mirror_rhs_apar_contribution(rhs, apar, imu, ikxkyz)

      use grids_species, only: spec
      use grids_velocity, only: nvpa
      use grids_velocity, only: maxwell_vpa, maxwell_mu, vpa
      use parallelisation_layouts, only: kxkyz_lo, is_idx, iz_idx
      use calculations_gyro_averages, only: gyro_average

      implicit none

      ! Arguments
      complex, dimension(:), intent(out) :: rhs
      complex, intent(in) :: apar
      integer, intent(in) :: imu, ikxkyz

      ! Local variables
      integer :: ia, iz, is
      real :: pre_factor
      complex, dimension(:), allocatable :: vpa_scratch

      !-------------------------------------------------------------------------

      allocate (vpa_scratch(nvpa))

      is = is_idx(kxkyz_lo, ikxkyz)
      iz = iz_idx(kxkyz_lo, ikxkyz)

      ia = 1
      pre_factor = -2.0 * spec(is)%zt * spec(is)%stm_psi0
      call gyro_average(pre_factor * vpa * apar, imu, ikxkyz, vpa_scratch)

      vpa_scratch = vpa_scratch * maxwell_vpa(:, is) * maxwell_mu(ia, iz, imu, is)

      rhs = vpa_scratch

      deallocate (vpa_scratch)

   end subroutine get_mirror_rhs_apar_contribution

   !****************************************************************************
   !                       ADD MIRROR TERM FOR RADIAL VARIATION
   !****************************************************************************
   subroutine add_mirror_radial_variation(g, gout)

      ! Parallelisation
      use mp, only: proc0
      use initialise_redistribute, only: kxkyz2vmu
      use job_manage, only: time_message
      use timers, only: time_mirror
      use redistribute, only: gather, scatter
      use parallelisation_layouts, only: kxkyz_lo, vmu_lo
      use parallelisation_layouts, only: is_idx, imu_idx
      use parallelisation_layouts, only: fields_kxkyz

      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: nvpa, nmu
      use arrays_distribution_function, only: gvmu
      use parameters_physics, only: full_flux_surface
      
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      complex, dimension(:, :, :), allocatable :: g0v
      integer :: iz, it, imu, is, ivmu, ia

      !-------------------------------------------------------------------------

      allocate (g0v(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

      ! if (proc0) call time_message(.false.,time_mirror(:,1),' Mirror global variation advance')

      ia = 1

      if (full_flux_surface) then
         ! FLAG DSO - Someday one should be able to do full global
      else
         if (.not. fields_kxkyz) then
            if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
            call scatter(kxkyz2vmu, g, gvmu)
            if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         end if
         ! Get dg/dvpa and store in g0v
         g0v = gvmu
         call get_dgdvpa_explicit(g0v)
         ! Swap layouts so that (z,kx,ky) are local
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')
         call gather(kxkyz2vmu, g0v, gout)
         if (proc0) call time_message(.false., time_mirror(:, 2), ' mirror_redist')

         ! Get mirror term and add to source
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

      ! if (proc0) call time_message(.false.,time_mirror(:,1),' Mirror global variation advance')

   end subroutine add_mirror_radial_variation

   !****************************************************************************
   !                               VPA INTERPOLAITON
   !****************************************************************************
   ! This is needed if advancing mirror using semi-lagrange
   !****************************************************************************
   subroutine vpa_interpolation(grid, interp)

      use grids_velocity, only: nvpa, nmu
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: iz_idx, is_idx
      use parameters_numerical, only: mirror_linear_interp

      implicit none

      ! Arguments
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: grid
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(out) :: interp

      ! Local variables
      integer :: ikxkyz, iz, is, iv, imu
      integer :: shift, sgn, llim, ulim
      real :: fac0, fac1, fac2, fac3
      real :: tmp0, tmp1, tmp2, tmp3

      !-------------------------------------------------------------------------

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
               ! Either assume BC for g is zero beyond grid extrema
               ! or dg/dvpa is zero beyond grid extrema
               ! at boundary cell, use zero incoming BC for point just beyond boundary
               interp(ulim + sgn, imu, ikxkyz) = grid(ulim + shift + sgn, imu, ikxkyz) * fac1
               ! Use zero incoming BC for cells beyond +/- nvgrid
               if (shift > 0) then
               ! interp(nvpa-shift+1:,imu,ikxkyz) = 0.

                  do iv = nvpa - shift, nvpa
                     interp(iv, imu, ikxkyz) = grid(nvpa, imu, ikxkyz) * real(nvpa - iv) / real(shift + 1)
                  end do

               else if (shift < 0) then
               ! interp(:-shift,imu,ikxkyz) = 0.

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
                  ! Deal with boundary point near outgoing BC using 2-point
                  ! (linear) interpolation -> could do 3-point to improve accuracy
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

               ! If llim > ulim (for sgn > 0) or llim < ulim (for sgn < 0)
               ! then there are no elements to be interpolated
               if (sgn * ulim < sgn * llim) then
                  interp(:, imu, ikxkyz) = 0.
               else
                  ! Coefficient multiplying g_{iv-1}
                  fac0 = -tmp0 * tmp1 * tmp2
                  ! Coefficient multiplying g_{iv}
                  fac1 = 3.*tmp1 * tmp2 * tmp3
                  ! Coefficient multiplying g_{iv+1}
                  fac2 = -3.*tmp0 * tmp1 * tmp3
                  ! Coefficient multiplying g_{iv+2}
                  fac3 = tmp0 * tmp2 * tmp3
                  do iv = llim, ulim, sgn
                     interp(iv, imu, ikxkyz) = (grid(iv + shift - sgn, imu, ikxkyz) * fac0 &
                        + grid(iv + shift, imu, ikxkyz) * fac1 &
                        + grid(iv + shift + sgn, imu, ikxkyz) * fac2 &
                        + grid(iv + shift + 2 * sgn, imu, ikxkyz) * fac3) / 6.
                  end do

                  ! At boundary cell, use zero incoming BC for point just beyond boundary
                  interp(ulim + sgn, imu, ikxkyz) = (grid(ulim + shift, imu, ikxkyz) * fac0 &
                     + grid(ulim + shift + sgn, imu, ikxkyz) * fac1 &
                     + grid(ulim + shift + 2 * sgn, imu, ikxkyz) * fac2) / 6.
                  interp(ulim + 2 * sgn, imu, ikxkyz) = (grid(ulim + shift + sgn, imu, ikxkyz) * fac0 &
                     + grid(ulim + shift + 2 * sgn, imu, ikxkyz) * fac1) / 6.
                  ! Use zero incoming BC for cells beyond +/- nvgrid
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

   !****************************************************************************
   !                         INVERT MIRROR OPERATOR 
   !****************************************************************************
   ! Use the Tridiagonal matrix algorithm to invert the mirror operator for 
   ! implicit mirror advance
   !****************************************************************************
   subroutine invert_mirror_operator(imu, ilo, g)

      use calculations_finite_differences, only: tridag

      implicit none

      integer, intent(in) :: imu, ilo
      complex, dimension(:), intent(in out) :: g

      !-------------------------------------------------------------------------

      call tridag(1, mirror_tri_a(:, imu, ilo), mirror_tri_b(:, imu, ilo), mirror_tri_c(:, imu, ilo), g)

   end subroutine invert_mirror_operator
   
!###############################################################################
!################################# FINISH MIRROR ###############################
!###############################################################################

   !****************************************************************************
   !                          FINISH MIRROR - DEALLOCATE
   !****************************************************************************
   subroutine finish_mirror

      use parameters_numerical, only: mirror_implicit, mirror_semi_lagrange

      implicit none

      !-------------------------------------------------------------------------

      if (allocated(mirror)) deallocate (mirror)
      if (allocated(mirror_sign)) deallocate (mirror_sign)
      if (allocated(mirror_rad_var)) deallocate (mirror_rad_var)

      if (mirror_implicit) then
         if (mirror_semi_lagrange) then
            call finish_mirror_semi_lagrange
         else
            call finish_invert_mirror_operator
            call finish_mirror_response
         end if
      end if

      initialised_mirror = .false.

   end subroutine finish_mirror

   !****************************************************************************
   !                        FINISH MIRROR: SEMI LAGRANGE
   !****************************************************************************
   subroutine finish_mirror_semi_lagrange

      implicit none

      if (allocated(mirror_interp_loc)) deallocate (mirror_interp_loc)
      if (allocated(mirror_interp_idx_shift)) deallocate (mirror_interp_idx_shift)

   end subroutine finish_mirror_semi_lagrange

   !****************************************************************************
   !                       FINISH MIRROR: IMPLICIT
   !****************************************************************************
   subroutine finish_invert_mirror_operator

      implicit none

      if (allocated(mirror_tri_a)) then
         deallocate (mirror_tri_a)
         deallocate (mirror_tri_b)
         deallocate (mirror_tri_c)
      end if

      if (allocated(mirror_int_fac)) deallocate (mirror_int_fac)

   end subroutine finish_invert_mirror_operator

   !****************************************************************************
   !                    FINISH MIRROR: IMPLICIT ELECTROMAGNETIC
   !****************************************************************************
   subroutine finish_mirror_response

      implicit none

      if (allocated(response_apar_denom)) deallocate (response_apar_denom)

   end subroutine finish_mirror_response

end module gk_mirror
