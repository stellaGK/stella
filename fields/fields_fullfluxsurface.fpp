module fields_ffs

   use common_types, only: coupled_alpha_type, gam0_ffs_type

   implicit none

   public :: init_fields_ffs
   public :: get_fields_ffs
   public :: finish_fields_ffs

   private

   !> arrays allocated/used if simulating a full flux surface
   type(coupled_alpha_type), dimension(:, :, :), allocatable :: gam0_ffs
   type(gam0_ffs_type), dimension(:, :), allocatable :: lu_gam0_ffs
   complex, dimension(:), allocatable :: adiabatic_response_factor

   logical :: debug = .false.

contains
   

!###############################################################################
!###################### ADVANCE FULL FLUX SURFACE FIELDS #######################
!###############################################################################


!> get_fields_ffs accepts as input the guiding centre distribution function g
!> and calculates/returns the electronstatic potential phi for full_flux_surface simulations
   subroutine get_fields_ffs(g, phi, apar)

      use mp, only: mp_abort
      use physics_parameters, only: nine, tite
      use stella_layouts, only: vmu_lo
      use run_parameters, only: fphi, fapar
      use species, only: modified_adiabatic_electrons, adiabatic_electrons
      use zgrid, only: nzgrid
      use kt_grids, only: nakx, ikx_max, naky, naky_all
      use kt_grids, only: swap_kxky_ordered
      use volume_averages, only: flux_surface_average_ffs

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar

      integer :: iz, ikx
      complex, dimension(:), allocatable :: phi_fsa
      complex, dimension(:, :, :), allocatable :: phi_swap, source

      if (fphi > epsilon(0.0)) then
         allocate (source(naky, nakx, -nzgrid:nzgrid))
         !> calculate the contribution to quasineutrality coming from the velocity space
         !> integration of the guiding centre distribution function g;
         !> the sign is consistent with phi appearing on the RHS of the eqn and int g appearing on the LHS.
         !> this is returned in source
         if (debug) write (*, *) 'fields::advance_fields::get_fields_ffs::get_g_integral_contribution'
         call get_g_integral_contribution(g, source)
         !> use sum_s int d3v <g> and QN to solve for phi
         !> NB: assuming here that ntubes = 1 for FFS sim
         if (debug) write (*, *) 'fields::advance_fields::get_phi_ffs'
         call get_phi_ffs(source, phi(:, :, :, 1))
         !> if using a modified Boltzmann response for the electrons, then phi
         !> at this stage is the 'inhomogeneous' part of phi.
         if (modified_adiabatic_electrons) then
            !> first must get phi on grid that includes positive and negative ky (but only positive kx)
            allocate (phi_swap(naky_all, ikx_max, -nzgrid:nzgrid))
            if (debug) write (*, *) 'fields::advance_fields::get_fields_ffs::swap_kxky_ordered'
            do iz = -nzgrid, nzgrid
               call swap_kxky_ordered(phi(:, :, iz, 1), phi_swap(:, :, iz))
            end do
            !> calculate the flux surface average of this phi_inhomogeneous
            allocate (phi_fsa(nakx))
            if (debug) write (*, *) 'fields::advance_fields::get_fields_ffs::flux_surface_average_ffs'
            do ikx = 1, nakx
               call flux_surface_average_ffs(phi_swap(:, ikx, :), phi_fsa(ikx))
            end do
            !> use the flux surface average of phi_inhomogeneous, together with the
            !> adiabatic_response_factor, to obtain the flux-surface-averaged phi
            phi_fsa = phi_fsa * adiabatic_response_factor
            !> use the computed flux surface average of phi as an additional sosurce in quasineutrality
            !> to obtain the electrostatic potential; only affects the ky=0 component of QN
            do ikx = 1, nakx
               source(1, ikx, :) = source(1, ikx, :) + phi_fsa(ikx) * tite / nine
            end do
            if (debug) write (*, *) 'fields::advance_fields::get_fields_ffs::get_phi_ffs2s'
            call get_phi_ffs(source, phi(:, :, :, 1))
            deallocate (phi_swap, phi_fsa)
         end if
         deallocate (source)
      else if (.not. adiabatic_electrons) then
         !> if adiabatic electrons are not employed, then
         !> no explicit equation for the ky=kx=0 component of phi;
         !> hack for now is to set it equal to zero.
         phi(1, 1, :, :) = 0.
      end if

      apar = 0.
      if (fapar > epsilon(0.0)) then
         call mp_abort('apar not yet supported for full_flux_surface = T. aborting.')
      end if

   contains

      subroutine get_g_integral_contribution(g, source)

         use mp, only: sum_allreduce
         use stella_layouts, only: vmu_lo
         use species, only: spec
         use zgrid, only: nzgrid
         use kt_grids, only: naky, nakx
         use vpamu_grids, only: integrate_species_ffs
         use gyro_averages, only: gyro_average, j0_B_maxwell_ffs

         implicit none

         complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
         complex, dimension(:, :, -nzgrid:), intent(in out) :: source

         integer :: it, iz, ivmu
         complex, dimension(:, :, :), allocatable :: gyro_g

         !> assume there is only a single flux surface being simulated
         it = 1
         allocate (gyro_g(naky, nakx, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         !> loop over zed location within flux tube
         do iz = -nzgrid, nzgrid
            !> loop over super-index ivmu, which include vpa, mu and spec
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               !> gyroaverage the distribution function g at each phase space location
               call gyro_average(g(:, :, iz, it, ivmu), gyro_g(:, :, ivmu), j0_B_maxwell_ffs(:, :, iz, ivmu))
            end do
            !> integrate <g> over velocity space and sum over species within each processor
            !> as v-space and species possibly spread over processors, wlil need to
            !> gather sums from each proceessor and sum them all together below
            call integrate_species_ffs(gyro_g, spec%z * spec%dens_psi0, source(:, :, iz), reduce_in=.false.)
         end do
         !> gather sub-sums from each processor and add them together
         !> store result in phi, which will be further modified below to account for polarization term
         call sum_allreduce(source)
         !> no longer need <g>, so deallocate
         deallocate (gyro_g)

      end subroutine get_g_integral_contribution

   end subroutine get_fields_ffs

   subroutine get_phi_ffs(rhs, phi)

      use zgrid, only: nzgrid
      use kt_grids, only: swap_kxky_ordered, swap_kxky_back_ordered
      use kt_grids, only: naky_all, ikx_max
      use gyro_averages, only: band_lu_solve_ffs

      implicit none

      complex, dimension(:, :, -nzgrid:), intent(in) :: rhs
      complex, dimension(:, :, -nzgrid:), intent(out) :: phi

      integer :: iz
      complex, dimension(:, :, :), allocatable :: rhs_swap

      allocate (rhs_swap(naky_all, ikx_max, -nzgrid:nzgrid))

      !> change from rhs defined on grid with ky >=0 and kx from 0,...,kxmax,-kxmax,...,-dkx
      !> to rhs_swap defined on grid with ky = -kymax,...,kymax and kx >= 0
      do iz = -nzgrid, nzgrid
         call swap_kxky_ordered(rhs(:, :, iz), rhs_swap(:, :, iz))
      end do

      !> solve sum_s Z_s int d^3v <g> = gam0*phi
      !> where sum_s Z_s int d^3v <g> is initially passed in as rhs_swap
      !> and then rhs_swap is over-written with the solution to the linear system
      call band_lu_solve_ffs(lu_gam0_ffs, rhs_swap)

      !> swap back from the ordered grid in ky to the original (kx,ky) grid
      do iz = -nzgrid, nzgrid
         call swap_kxky_back_ordered(rhs_swap(:, :, iz), phi(:, :, iz))
      end do

      deallocate (rhs_swap)

   end subroutine get_phi_ffs

!###############################################################################
!############################ INITALISE AND FINALISE ###########################
!###############################################################################

   !============================================================================
   !============================= INITALISE ARRAYS =============================
   !============================================================================ 

   subroutine init_fields_ffs

      use species, only: modified_adiabatic_electrons

      implicit none

      !> calculate and LU factorise the matrix multiplying the electrostatic potential in quasineutrality
      !> this involves the factor 1-Gamma_0(kperp(alpha))
      call init_gamma0_factor_ffs

      !> if using a modified Boltzmann response for the electrons
      if (modified_adiabatic_electrons) then
         !> obtain the response of phi_homogeneous to a unit perturbation in flux-surface-averaged phi
         call init_adiabatic_response_factor
      end if

   end subroutine init_fields_ffs

   !> calculate and LU factorise the matrix multiplying the electrostatic potential in quasineutrality
   !> this involves the factor 1-Gamma_0(kperp(alpha))
   subroutine init_gamma0_factor_ffs

      use mp, only: proc0
      use spfunc, only: j0
      use dist_fn_arrays, only: kperp2
      use stella_transforms, only: transform_alpha2kalpha
      use physics_parameters, only: nine, tite
      use species, only: spec, nspec
      use species, only: adiabatic_electrons
      use zgrid, only: nzgrid
      use stella_geometry, only: bmag
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use kt_grids, only: nalpha, ikx_max, naky_all, naky
      use kt_grids, only: swap_kxky_ordered
      use vpamu_grids, only: vperp2, maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: integrate_species
      use gyro_averages, only: band_lu_factorisation_ffs

      implicit none

      integer :: iky, ikx, iz, ia
      integer :: ivmu, iv, imu, is
      real :: arg

      real, dimension(:, :, :), allocatable :: kperp2_swap
      real, dimension(:), allocatable :: aj0_alpha, gam0_alpha
      real, dimension(:), allocatable :: wgts
      complex, dimension(:), allocatable :: gam0_kalpha

      debug = debug .and. proc0

      if (debug) write (*, *) 'fields::init_fields::init_gamm0_factor_ffs'

      allocate (kperp2_swap(naky_all, ikx_max, nalpha))
      allocate (aj0_alpha(vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      allocate (gam0_alpha(nalpha))
      allocate (gam0_kalpha(naky))
      !> wgts are species-dependent factors appearing in Gamma0 factor
      allocate (wgts(nspec))
      wgts = spec%dens * spec%z**2 / spec%temp
      !> allocate gam0_ffs array, which will contain the Fourier coefficients in y
      !> of the Gamma0 factor that appears in quasineutrality
      if (.not. allocated(gam0_ffs)) then
         allocate (gam0_ffs(naky_all, ikx_max, -nzgrid:nzgrid))
      end if

      do iz = -nzgrid, nzgrid
         !> in calculating the Fourier coefficients for Gamma_0, change loop orders
         !> so that inner loop is over ivmu super-index;
         !> this is done because we must integrate over v-space and sum over species,
         !> and we want to minimise memory usage where possible (so, e.g., aj0_alpha need
         !> only be a function of ivmu and can be over-written for each (ia,iky,ikx)).
         do ia = 1, nalpha
            call swap_kxky_ordered(kperp2(:, :, ia, iz), kperp2_swap(:, :, ia))
         end do
         do ikx = 1, ikx_max
            do iky = 1, naky_all
               do ia = 1, nalpha
                  !> get J0 for all vpar, mu, spec values
                  do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                     is = is_idx(vmu_lo, ivmu)
                     imu = imu_idx(vmu_lo, ivmu)
                     iv = iv_idx(vmu_lo, ivmu)
                     !> calculate the argument of the Bessel function J0
                     arg = spec(is)%bess_fac * spec(is)%smz_psi0 * sqrt(vperp2(ia, iz, imu) * kperp2_swap(iky, ikx, ia)) / bmag(ia, iz)
                     !> compute J0 corresponding to the given argument arg
                     aj0_alpha(ivmu) = j0(arg)
                     !> form coefficient needed to calculate 1-Gamma_0
                     aj0_alpha(ivmu) = (1.0 - aj0_alpha(ivmu)**2)
                  end do

                  !> calculate gamma0(kalpha,alpha,...) = sum_s Zs^2 * ns / Ts int d3v (1-J0^2)*F_{Maxwellian}
                  !> note that v-space Jacobian contains alpha-dependent factor, B(z,alpha),
                  !> but this is not a problem as we have yet to transform from alpha to k_alpha
                  call integrate_species(aj0_alpha, iz, wgts, gam0_alpha(ia), ia)
                  !> if Boltzmann response used, account for non-flux-surface-averaged component of electron density
                  if (adiabatic_electrons) then
                     gam0_alpha(ia) = gam0_alpha(ia) + tite / nine
                  else if (ikx == 1 .and. iky == naky) then
                     !> if kx = ky = 0, 1-Gam0 factor is zero;
                     !> this leads to eqn of form 0 * phi_00 = int d3v g.
                     !> hack for now is to set phi_00 = 0, as above inversion is singular.
                     !> to avoid singular inversion, set gam0_alpha = 1.0
                     gam0_alpha(ia) = 1.0
                  end if
               end do
               !> fourier transform Gamma_0(alpha) from alpha to k_alpha space
               call transform_alpha2kalpha(gam0_alpha, gam0_kalpha)
               gam0_ffs(iky, ikx, iz)%max_idx = naky
               !> allocate array to hold the Fourier coefficients
               if (.not. associated(gam0_ffs(iky, ikx, iz)%fourier)) &
                  allocate (gam0_ffs(iky, ikx, iz)%fourier(gam0_ffs(iky, ikx, iz)%max_idx))
               !> fill the array with the requisite coefficients
               gam0_ffs(iky, ikx, iz)%fourier = gam0_kalpha(:gam0_ffs(iky, ikx, iz)%max_idx)
!                call test_ffs_bessel_coefs (gam0_ffs(iky,ikx,iz)%fourier, gam0_alpha, iky, ikx, iz, gam0_ffs_unit)
            end do
         end do
      end do

      !> LU factorise array of gam0, using the LAPACK zgbtrf routine for banded matrices
      if (.not. allocated(lu_gam0_ffs)) then
         allocate (lu_gam0_ffs(ikx_max, -nzgrid:nzgrid))
!          call test_band_lu_factorisation (gam0_ffs, lu_gam0_ffs)
         call band_lu_factorisation_ffs(gam0_ffs, lu_gam0_ffs)
      end if

      deallocate (wgts)
      deallocate (kperp2_swap)
      deallocate (aj0_alpha, gam0_alpha)
      deallocate (gam0_kalpha)

   end subroutine init_gamma0_factor_ffs

   !> solves Delta * phi_hom = -delta_{ky,0} * ne/Te for phi_hom
   !> this is the vector describing the response of phi_hom to a unit impulse in phi_fsa
   !> it is the sum over ky and integral over kx of this that is needed, and this
   !> is stored in adiabatic_response_factor
   subroutine init_adiabatic_response_factor

      use physics_parameters, only: nine, tite
      use zgrid, only: nzgrid
      use stella_transforms, only: transform_alpha2kalpha
      use kt_grids, only: naky, naky_all, ikx_max
      use gyro_averages, only: band_lu_solve_ffs
      use volume_averages, only: flux_surface_average_ffs

      implicit none

      integer :: ikx
      complex, dimension(:, :, :), allocatable :: adiabatic_response_vector

      allocate (adiabatic_response_vector(naky_all, ikx_max, -nzgrid:nzgrid))
      if (.not. allocated(adiabatic_response_factor)) allocate (adiabatic_response_factor(ikx_max))

      !> adiabatic_response_vector is initialised to be the rhs of the equation for the
      !> 'homogeneous' part of phi, with a unit impulse assumed for the flux-surface-averaged phi
      !> only the ky=0 component contributes to the flux-surface-averaged potential
      adiabatic_response_vector = 0.0
      adiabatic_response_vector(naky, :, :) = tite / nine
      !> pass in the rhs and overwrite with the solution for phi_homogeneous
      call band_lu_solve_ffs(lu_gam0_ffs, adiabatic_response_vector)

      !> obtain the flux surface average of the response vector
      do ikx = 1, ikx_max
         call flux_surface_average_ffs(adiabatic_response_vector(:, ikx, :), adiabatic_response_factor(ikx))
      end do
      adiabatic_response_factor = 1.0 / (1.0 - adiabatic_response_factor)

      deallocate (adiabatic_response_vector)

   end subroutine init_adiabatic_response_factor

   !============================================================================
   !============================= ALLOCATE ARRAYS ==============================
   !============================================================================

   !> TODO-GA: add allocate fields subroutine 

   !============================================================================
   !========================== FINISH THE FFS FIELDS ===========================
   !============================================================================
   subroutine finish_fields_ffs

      implicit none

      !> arrays only allocated/used if simulating a full flux surface
      if (allocated(gam0_ffs)) deallocate (gam0_ffs)
      if (allocated(lu_gam0_ffs)) deallocate (lu_gam0_ffs)
      if (allocated(adiabatic_response_factor)) deallocate (adiabatic_response_factor)

   end subroutine finish_fields_ffs

end module fields_ffs
