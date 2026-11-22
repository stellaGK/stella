!###############################################################################
!                   Linearised Fokker-Planck Collision Model                    
!###############################################################################
! 
! This module adds collisions to the gyrokinetic equation, based on gyrokinetic, 
! linearised Fokker–Planck collision model
! 
! See [2024 - Von Boetticher - Linearised Fokker–Planck collision model for gyrokinetic simulations]
! 
! The basic Spitzer collision frequency nu^ab is given by equation (12) in [2024 - Von Boetticher]:
!     nu^ab = n_b * e_a^2 * e_b^2 * ln Lambda / (4 * pi * varepsilon_0^2 * m_a^2 * v_{th,a}^3)
! 
! The operator C^ab[h_a, h_b] = C_test^ab[h_a] + C_field^ab[h_b] of species a colliding with species b has
!   - the test-particle part C_test: pitch-angle scattering (nuD), energy diffusion (nupa), the mixed
!     vpa-mu derivatives (nux = nupa - nuD) and the classical (finite-Larmor-radius) diffusion;
!   - the field-particle part C_field: an expansion in Legendre polynomials (l <= lmax) and associated
!     Laguerre polynomials (j <= jmax) of the field-particle operator of Hirshman and Sigmar, with the
!     functions Delta_j^l(v) (deltaj) and the moments psi^{ab,jlm} of h_b, which restore the momentum and
!     energy exchange; with lmax = jmax = 1 it holds the density, momentum and energy terms.
! The options are read in namelist_dissipation.f90 (<collisions_fokker_planck> for users and
! <collisions_fokker_planck_developers> for tests of the operator); the variables of this module have the
! names of these options, except that the scale_* factors also hold the include_* switches (0 = off).
!
! Implicit algorithm (advance_collisions_fp_implicit, the default), backward Euler in time:
!   init_collisions_fp builds, on every (ky, kx, z, tube) point of kxkyzt_lo, the block-tridiagonal matrix
!   1 - Deltat C_test (blocks aa_blcs, bb_blcs, cc_blcs of nmu x nmu per vpa row, in band storage
!   cdiffmat_band, LU factorised; init_fp_diffmatrix), the functions Delta_j^l (init_deltaj_vmu) and the
!   response matrix fp_response of phi and the psi moments (init_fp_conserve). Every time step,
!   advance_implicit_fp
!     1. converts g to h and solves [1 - Deltat C_test] h_inh = h (solve_test_particle_operator);
!     2. takes phi_inh and the psi_inh moments of h_inh (get_psi);
!     3. solves the response matrix for phi^{n+1} and psi^{n+1};
!     4. solves [1 - Deltat C_test] h = h + Ze/T J0 phi^{n+1} F0 + sum_jlm psi^{n+1} Delta_j^l and converts
!        h back to g.
!   Without the field-particle part, the response matrix only holds phi.
! Explicit algorithm (advance_collisions_fp_explicit, a developer option): the test-particle part of
!   like-species collisions only, applied to h on kxkyz_lo, for comparisons with the implicit algorithm.
!
! exact_conservation conserves momentum and energy to machine precision; see init_collisions_fp for its
! two algorithms (block products on the equally spaced mu grid, the test-particle matrix on the
! Gauss-Laguerre grid).
!
! The electron-ion options (scale_energy_diffusion_ei, scale_pitch_angle_scattering_ei,
! mass_ratio_approximation_ei) and the Spitzer problem take the ions as species 1 and the electrons as
! species 2; stella stops otherwise (init_collisions_fp).
! 
! TODO - stella changed a lot since collisions have been implemented
!        compare the collisions with an older stella version!
!
! POSSIBLE IMPROVEMENTS (from the commented-out code removed in Nov 2026, see the git history)
!   - The functions Delta_j^l should satisfy the self-adjointness of the field-particle operator,
!     int v^l L_j Delta_j^{l,ab} d^3v = (m_a/m_b)^{7/2} int v^l L_j Delta_j^{l,ba} d^3v (up to the
!     normalisation), which was checked by hand with printed integrals in init_deltaj_vmu. This could
!     become an automatic test of the field-particle part on the grid.
!   - The mu derivatives at the ends of the mu grid (imu = 1 and nmu) are first order. Second-order
!     one-sided stencils were written but commented out in the explicit operator (and in an old version
!     of the implicit matrix); they would have to be used in both algorithms together, and checked for
!     the density conservation of conservative_test_particle_density.
!   - Several parts are near copies that could share code: the matrix assembly for iv = 1, nvpa and
!     the interior, the half-point collision frequencies (written inline about 20 times, and in
!     init_nusDpa), the mm = 0, mm > 0 and mm < 0 branches of the psi moments, and the isb = 0 and
!     isb /= 0 branches of get_psi.
!
! OPEN QUESTION - The vpa integration weights are not the same in the whole implicit step (Nov 2026)
!        The vpa integrals use the weights that set_vpa_weights set last: the conservative weights
!        (dvpa/sqrt(pi)) or the weights of the velocity grid (Simpson's rule by default). The collision
!        step of dissipation_and_collisions sets the conservative weights before it calls
!        advance_collisions_fp_implicit, and the Fokker-Planck routines switch them several times. As a
!        result, with the field-particle part, the inhomogeneous potential phi_inh of the response-matrix
!        solve and the first density correction use the conservative weights, while the response matrix
!        (fp_response, built in the initialisation) and the rest of the step use the weights of the grid.
!        Using one set of weights in the whole step (the conservative weights if conservative_test_particle_density,
!        otherwise those of the grid, as in the initialisation) changes phi^2 by about 1e-7 to 1e-9 in the
!        tests with the field-particle part, and by about 5e-4 in the implicit test without it. It has
!        not been shown which version is more accurate, so the code is unchanged for now. A test that
!        settles it: compare the growth rate and the conservation of density, momentum and energy of
!        both versions against a converged reference (e.g. GS2) at increasing nvpa.
!
!###############################################################################
! PARALLELISATION NOTES (kxkyzt_lo refactor)
!
! The implicit collision step in this module uses a banded LU factorisation of
! a (nvpa*nmu) x (nvpa*nmu) operator that couples all mu values together
! (and all vpa values together) within each (kx, ky, z, tube, species)
! point. The matrix arrays aa_blcs / bb_blcs / cc_blcs (and the band-storage
! cdiffmat_band) therefore require both vpa and mu to be local on the
! processor that performs the solve.
!
! Until May 2026 the matrix data was kept on kxkyz_lo (which parallelises
! kx, ky, z, tubes AND species).
! 
! The kxkyzt_lo layout (parallelises kx, ky, z, tubes only, and keeps vpa, mu, 
! and species all local on each processor) is the natural layout for this
! solve. Each rank that owns an (iky,ikx,iz,it) point owns its full per-
! species band matrix; cross-species accesses become local operations and
! the global replication of cdiffmat_band / ipiv disappears. The implicit
! solve still keeps the entire (vpa,mu) axis on a single rank, so the
! existing LAPACK zgbtrf / zgbtrs band routines apply unchanged.
!
! Redistribute kxkyzt2vmu (declared in initialise_redistribute) moves data
! between
!   g_kxkyzt(nvpa, nmu, nspec, kxkyzt-layout)
! and
!   g(naky, nakx, -nzgrid:nzgrid, ntubes, vmu-layout)
! and replaces the previous kxkyz2vmu hop for the implicit collision step.
!
!###############################################################################
module collisions_fokkerplanck

   implicit none

   ! Make routines available to other modules
   public :: read_parameters_fp
   public :: init_collisions_fp, finish_collisions_fp
   public :: reset_dt_collisions_fp
   public :: advance_collisions_fp_explicit, advance_collisions_fp_implicit

   private

   logical :: include_vpa_derivatives, include_mu_derivatives
   logical :: conservative_test_particle_density, exact_field_particle_density
   logical :: field_particle_density_from_test_particle
   ! exact_conservation is the user's switch; init_collisions_fp turns it into one of the two algorithms
   logical :: exact_conservation, exact_conservation_block_products, exact_conservation_test_particle_matrix
   logical :: spitzer_problem, no_j1l1, no_j1l2, no_j0l2
   logical :: include_field_particle_part, include_test_particle_part
   logical :: include_unlike_species_collisions, include_like_species_collisions
   ! The algorithm of exact_conservation: 'automatic' (from the mu grid), 'block_products' or 'test_particle_matrix'
   character(len=30) :: exact_conservation_algorithm

   integer :: nresponse = 1
   real :: scale_flr_terms
   real :: scale_mixed_derivative_terms
   real :: scale_nu_ii, scale_nu_ie, scale_nu_ee, scale_nu_ei
   ! scale_pitch_angle_scattering is 1 or 0 (include_pitch_angle_scattering)
   real :: scale_energy_diffusion_ei, scale_pitch_angle_scattering_ei, scale_pitch_angle_scattering
   logical :: mass_ratio_approximation_ei
   integer :: jmax = 1
   integer :: lmax = 1

   complex, dimension(:, :, :), allocatable :: fp_response
   integer, dimension(:, :), allocatable :: diff_idx

   ! aa_blcs/bb_blcs/cc_blcs are distributed on kxkyzt_lo with explicit
   ! (target_is, source_isb) axes, shape:
   !   (nvpa, nmu, nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec, nspec)
   ! Was previously (nvpa, nmu, nmu, kxkyz_lo, nspec) with target species
   ! implicit in kxkyz_lo.
   complex, dimension(:, :, :, :, :, :), allocatable :: aa_blcs, cc_blcs
   complex, dimension(:, :, :, :, :, :), allocatable :: bb_blcs
   ! cdiffmat_band has shape
   !   (3*(nmu+1)+1, nmu*nvpa, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec)
   complex, dimension(:, :, :, :), allocatable :: cdiffmat_band
   ! blockmatrix has shape
   !   (nvpa*nmu, nvpa*nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec, nspec)
   complex, dimension(:, :, :, :, :), allocatable :: blockmatrix
   ! blockmatrix_sum has shape
   !   (nvpa*nmu, nvpa*nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec)
   complex, dimension(:, :, :, :), allocatable :: blockmatrix_sum
   ! ipiv has shape
   !   (nvpa*nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec)
   integer, dimension(:, :, :), allocatable :: ipiv
   real, dimension(:, :, :, :, :), allocatable :: nuD, nupa, nux
   real, dimension(:, :, :, :), allocatable :: mw, modmw
   real, dimension(:, :, :), allocatable :: velvpamu
   integer :: info

   real, dimension(:, :, :, :, :, :, :, :), allocatable :: deltaj, deltaj_tp

   ! The z points at which this processor stores deltaj, deltaj_tp, legendre_vpamu and jm, see init_local_z:
   ! izl(iz) = 1, ..., nz_local is the slot of iz in these arrays (0 if it is not stored here), iz_of_izl(izl) the z
   ! index of a slot, local_z(iz) whether iz is stored here (also the z points on which calc_deltaj_vmu, calc_psi_vmu
   ! and vLj_vmu compute), and iz_share_lo to iz_share_hi the share of the z points on which this processor computes psijnorm
   integer :: nz_local, iz_share_lo, iz_share_hi
   integer, dimension(:), allocatable :: izl, iz_of_izl
   logical, dimension(:), allocatable :: local_z
   ! The test-particle operator -code_dt C^{ab} applied to vpa F0a (1) and v^2 F0a (2), at the point (naky, 1, iz, 1)
   ! of every z point of this processor that it owns on kxkyzt_lo (exact_conservation, see apply_testpart_matrix)
   real, dimension(:, :, :, :, :), allocatable :: tp_f0moments
   real, dimension(:, :, :, :, :), allocatable :: psijnorm
   real, dimension(:, :, :, :, :), allocatable :: legendre_vpamu
   real, dimension(:, :, :, :, :, :), allocatable :: jm
   real, dimension(:), allocatable :: mwnorm
   real, dimension(:), allocatable :: modmwnorm
   ! With kinetic electrons, restore the density of the test-particle solve on the velocity grid (see solve_test_particle_operator)
   logical :: density_correction_solve
   real :: spitzer_scale_epar, spitzer_scale_grad_t

   ! Only initialise once
   logical :: initialised_fokker_planck = .false.

contains

!###############################################################################
!################################ READ PARAMETERS ##############################
!###############################################################################

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine read_parameters_fp(collisions_implicit)

      use mp, only: broadcast
      use namelist_dissipation, only: read_namelist_collisions_fokker_planck

      implicit none

      ! In: whether <dissipation_and_collisions_options> asks for implicit collisions. Out: .false. if the
      ! explicit algorithm is used, which is chosen with use_explicit_algorithm in the developer namelist
      logical, intent(in out) :: collisions_implicit

      logical :: explicit_algorithm

      !----------------------------------------------------------------------

      ! Read the namelists <collisions_fokker_planck> and <collisions_fokker_planck_developers>; the options are
      ! returned under the internal names of this module, see read_namelist_collisions_fokker_planck
      call read_namelist_collisions_fokker_planck(collisions_implicit, include_test_particle_part, &
               include_field_particle_part, lmax, jmax, include_unlike_species_collisions, &
               include_like_species_collisions, scale_nu_ii, scale_nu_ie, scale_nu_ee, scale_nu_ei, &
               scale_energy_diffusion_ei, scale_pitch_angle_scattering_ei, scale_pitch_angle_scattering, &
               mass_ratio_approximation_ei, spitzer_problem, conservative_test_particle_density, &
               exact_field_particle_density, field_particle_density_from_test_particle, &
               exact_conservation, exact_conservation_algorithm, include_vpa_derivatives, include_mu_derivatives, &
               scale_flr_terms, scale_mixed_derivative_terms, spitzer_scale_epar, spitzer_scale_grad_t, &
               no_j1l1, no_j1l2, no_j0l2, explicit_algorithm)

      ! Broadcast the parameters to all processors
      call broadcast(include_field_particle_part)
      call broadcast(include_test_particle_part)
      call broadcast(include_unlike_species_collisions)
      call broadcast(include_like_species_collisions)
      call broadcast(scale_nu_ii)
      call broadcast(scale_nu_ie)
      call broadcast(scale_nu_ee)
      call broadcast(scale_nu_ei)
      call broadcast(scale_energy_diffusion_ei)
      call broadcast(scale_pitch_angle_scattering)
      call broadcast(mass_ratio_approximation_ei)
      call broadcast(scale_pitch_angle_scattering_ei)
      call broadcast(conservative_test_particle_density)
      call broadcast(exact_field_particle_density)
      call broadcast(field_particle_density_from_test_particle)
      call broadcast(exact_conservation)
      call broadcast(exact_conservation_algorithm)
      call broadcast(spitzer_problem)
      call broadcast(include_vpa_derivatives)
      call broadcast(include_mu_derivatives)
      call broadcast(scale_flr_terms)
      call broadcast(scale_mixed_derivative_terms)
      call broadcast(jmax)
      call broadcast(lmax)
      call broadcast(spitzer_scale_epar)
      call broadcast(spitzer_scale_grad_t)
      call broadcast(no_j1l1)
      call broadcast(no_j1l2)
      call broadcast(no_j0l2)
      call broadcast(explicit_algorithm)

      if (explicit_algorithm) collisions_implicit = .false.

      ! The user's exact_conservation is turned into one of the two algorithms in init_collisions_fp, once the
      ! mu grid is known
      exact_conservation_block_products = .false.
      exact_conservation_test_particle_matrix = .false.

   end subroutine read_parameters_fp

!###############################################################################
!########################## INITIALISE COLLISION MODEL #########################
!###############################################################################

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_collisions_fp(collisions_implicit, cfl_dt_collisions)

      use grids_species, only: spec, nspec, has_electron_species
      use grids_velocity, only: dvpa, dmu, mu, nmu, equally_spaced_mu_grid, conservative_wgts_vpa
      use geometry, only: bmag
      use parallelisation_layouts
      use parameters_numerical, only: fully_explicit
      use common_types, only: spec_type
      use mp, only: mp_abort, nproc, proc0

      implicit none

      logical, intent(in) :: collisions_implicit
      real, intent(out) :: cfl_dt_collisions

      integer :: is, is2
      integer, parameter :: ion_species = 1
      integer, parameter :: electron_species = 2
      real :: vnew_max

      !-------------------------------------------------------------------------

      if (initialised_fokker_planck) return
      initialised_fokker_planck = .true.

      ! The field-particle part is only implemented in the implicit algorithm: the explicit algorithm
      ! (advance_collisions_fp_explicit) applies the test-particle operator only, so stop instead of
      ! silently running without the field-particle part
      if (include_field_particle_part .and. .not. collisions_implicit) then
         call mp_abort('The field-particle part of the Fokker-Planck operator (include_field_particle_part = .true.) ' &
            //'is only implemented in the implicit algorithm. Aborting.')
      end if

      ! The explicit algorithm only has like-species collisions (it uses vnew(is) of the species itself), so stop
      ! instead of silently dropping the collisions between different species
      if (.not. collisions_implicit .and. include_unlike_species_collisions .and. nspec > 1) then
         call mp_abort('The explicit Fokker-Planck algorithm only has collisions between particles of the same ' &
            //'species: with several species, set include_unlike_species_collisions = .false. in ' &
            //'&collisions_fokker_planck_developers. Aborting.')
      end if

      ! exact_conservation (conserve momentum and energy to machine precision) has two algorithms, chosen from
      ! the mu grid unless the developer option exact_conservation_algorithm says otherwise:
      !  - on the equally spaced mu grid, 'block_products' (exact_conservation_block_products): the test-particle
      !    operator applied to vpa F0 and v^2 F0 is taken from the blocks of the matrix (cheap, MPI parallel);
      !  - on the Gauss-Laguerre mu grid, 'test_particle_matrix' (exact_conservation_test_particle_matrix): the
      !    field-particle terms are computed with the dense test-particle matrix (expensive, and with several
      !    species not MPI parallel).
      ! Both need the discretisation of the test-particle part that conserves density
      ! (conservative_test_particle_density), set here: without it the energy is only conserved to about 1e-4
      ! (equally spaced) or 1e-2 (Gauss-Laguerre). On the Gauss-Laguerre grid that discretisation alone does not
      ! conserve the density, so the density term of field_particle_density_from_test_particle is added too.
      ! With ions and electrons (test 4g), the momentum and energy are then conserved to about 1e-12 on both
      ! grids, the density to 1e-11 (equally spaced) and 1e-8 (Gauss-Laguerre).
      if (exact_conservation) then
         if (.not. include_field_particle_part) call mp_abort('exact_conservation needs the field-particle ' &
            //'part of the Fokker-Planck operator (include_field_particle_part = .true.). Aborting.')
         if (.not. include_test_particle_part) call mp_abort('exact_conservation needs the test-particle ' &
            //'part of the Fokker-Planck operator (include_test_particle_part = .true.). Aborting.')
         if (lmax /= 1 .or. jmax /= 1) call mp_abort('exact_conservation only works with ' &
            //'field_particle_lmax = field_particle_jmax = 1. Aborting.')
         ! The exact conservation terms assume that the test-particle operator has no mixed (nux) derivative
         ! terms, whose discretisation does not conserve; with them, exact_conservation conserves the energy
         ! worse than without it. Dropping them changes the operator, so it is left to the user.
         if (scale_mixed_derivative_terms /= 0.) call mp_abort('exact_conservation only works without the ' &
            //'mixed-derivative terms of the test-particle part: set include_mixed_derivative_terms = .false. in ' &
            //'&collisions_fokker_planck. Aborting.')
         ! Without the conservative vpa weights the energy is only conserved to 1e-5 to 1e-2
         if (.not. conservative_wgts_vpa) call mp_abort('exact_conservation needs the conservative vpa ' &
            //'integration weights: set conservative_wgts_vpa = .true. in &velocity_grids. Aborting.')
         if (exact_conservation_algorithm == 'test_particle_matrix' .or. &
             (exact_conservation_algorithm == 'automatic' .and. .not. equally_spaced_mu_grid)) then
            exact_conservation_test_particle_matrix = .true.
            if (.not. equally_spaced_mu_grid) field_particle_density_from_test_particle = .true.
            if (proc0) write (*, '(A)') 'Warning: exact_conservation on the Gauss-Laguerre mu grid uses the ' &
               //'dense test-particle matrix, (nvpa*nmu)^2 complex numbers per point and species pair, and a ' &
               //'matrix product every time step; the equally spaced mu grid (equally_spaced_mu_grid = .true.) ' &
               //'is much cheaper.'
         else
            if (.not. equally_spaced_mu_grid) call mp_abort('exact_conservation_algorithm = "block_products" ' &
               //'needs the equally spaced mu grid (equally_spaced_mu_grid = .true.). Aborting.')
            exact_conservation_block_products = .true.
         end if
         conservative_test_particle_density = .true.
         if (proc0) then
            if (exact_conservation_block_products) then
               write (*, '(A)') 'exact_conservation: using the block-products algorithm, with ' &
                  //'conservative_test_particle_density.'
            else
               if (field_particle_density_from_test_particle) then
                  write (*, '(A)') 'exact_conservation: using the test-particle-matrix algorithm, with ' &
                     //'conservative_test_particle_density and field_particle_density_from_test_particle.'
               else
                  write (*, '(A)') 'exact_conservation: using the test-particle-matrix algorithm, with ' &
                     //'conservative_test_particle_density.'
               end if
            end if
         end if
      end if

      ! With several species, these options use the test-particle matrices (blockmatrix) of another
      ! species at the same (ky,kx,z) point (get_psi, get_testpart_density), which are only available
      ! when that block is on the same processor; otherwise the term is silently skipped
      if (include_field_particle_part .and. (exact_conservation_test_particle_matrix &
          .or. field_particle_density_from_test_particle) .and. nspec > 1 .and. nproc > 1) then
         call mp_abort('exact_conservation on the Gauss-Laguerre mu grid and field_particle_density_from_test_particle ' &
            //'are not MPI parallel for several species: run on one processor, or use exact_conservation on the ' &
            //'equally spaced mu grid (equally_spaced_mu_grid = .true.). Aborting.')
      end if

      ! include_vpa_derivatives is only used by the explicit algorithm: the implicit matrix always has the vpa terms
      if (.not. include_vpa_derivatives .and. collisions_implicit) then
         call mp_abort('include_vpa_derivatives = .false. is only implemented in the explicit Fokker-Planck ' &
            //'algorithm (use_explicit_algorithm = .true.). Aborting.')
      end if

      ! The electron-ion developer options and the Spitzer problem take the ions as species 1 and the electrons as
      ! species 2 (is == 2, isb == 1 in init_nusDpa, init_fp_diffmatrix and advance_implicit_fp)
      if ((scale_energy_diffusion_ei /= 1. .or. scale_pitch_angle_scattering_ei /= 1. .or. mass_ratio_approximation_ei &
           .or. spitzer_problem) .and. nspec > 1) then
         if (spec(1)%type /= ion_species .or. spec(2)%type /= electron_species) then
            call mp_abort('The electron-ion options of the Fokker-Planck operator (scale_energy_diffusion_ei, ' &
               //'scale_pitch_angle_scattering_ei, mass_ratio_approximation_ei) and spitzer_problem need the ions ' &
               //'as species 1 and the electrons as species 2. Aborting.')
         end if
      end if

      ! disable inter-species collisions if include_unlike_species_collisions==false
      if (.not. include_unlike_species_collisions) then
         do is = 1, nspec
            do is2 = 1, nspec
               if (is /= is2) then
                  spec(is)%vnew(is2) = 0.
               end if
            end do
         end do
      end if

      ! disable intra-species collisions if include_like_species_collisions==false
      if (.not. include_like_species_collisions) then
         do is = 1, nspec
            do is2 = 1, nspec
               if (is == is2) then
                  spec(is)%vnew(is2) = 0.
               end if
            end do
         end do
      end if

      ! control inter-species collisions
      do is = 1, nspec
         do is2 = 1, nspec
            if ((spec(is)%type == ion_species) .and. (spec(is2)%type == ion_species)) then
               spec(is)%vnew(is2) = spec(is)%vnew(is2) * scale_nu_ii
            else if ((spec(is)%type == ion_species) .and. (spec(is2)%type == electron_species)) then
               spec(is)%vnew(is2) = spec(is)%vnew(is2) * scale_nu_ie
            else if ((spec(is)%type == electron_species) .and. (spec(is2)%type == electron_species)) then
               spec(is)%vnew(is2) = spec(is)%vnew(is2) * scale_nu_ee
            else if ((spec(is)%type == electron_species) .and. (spec(is2)%type == ion_species)) then
               spec(is)%vnew(is2) = spec(is)%vnew(is2) * scale_nu_ei
            end if
            ! AVB: to do - add impurity collision control
         end do
      end do

      ! initialise speed dependent collision frequencies
      call init_nusDpa

      density_correction_solve = has_electron_species(spec) .and. (.not. conservative_test_particle_density)

      if (collisions_implicit) then
         fully_explicit = .false.

         call init_local_z
         call init_legendre
         call init_bessel_fn
         call init_fp_diffmatrix
         call init_deltaj_vmu
         call init_fp_conserve
      else
         vnew_max = 0.0
         do is = 1, nspec
            vnew_max = max(vnew_max, maxval(spec(is)%vnew))
         end do
         cfl_dt_collisions = min(2.0 * dvpa**2 / vnew_max, minval(bmag) / (vnew_max * maxval(mu(2:) / dmu(:nmu - 1)**2)))
      end if

   end subroutine init_collisions_fp

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_nusDpa

      ! AVB: compute the collision frequencies nuD and nupa

      use grids_z, only: nzgrid
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_velocity, only: nmu, mu, vpa, nvpa
      use geometry, only: bmag
      use grids_species, only: spec, nspec
      use spfunc, only: erf => erf_ext
      use calculations_finite_differences, only: fd3pt
      use grids_velocity, only: maxwell_mu, maxwell_vpa
      use constants, only: pi

      implicit none

      real, dimension(-nzgrid:nzgrid) :: v2mwint, v4mwint
      integer :: ia, imu, iv, iz, is, isb
      real :: x, Gf, massr

      !-------------------------------------------------------------------------

      if (.not. allocated(nuD)) allocate (nuD(nvpa, nmu, -nzgrid:nzgrid, nspec, nspec))
      if (.not. allocated(nupa)) allocate (nupa(nvpa, nmu, -nzgrid:nzgrid, nspec, nspec))
      if (.not. allocated(nux)) allocate (nux(nvpa, nmu, -nzgrid:nzgrid, nspec, nspec))
      if (.not. allocated(mw)) allocate (mw(nvpa, nmu, -nzgrid:nzgrid, nspec))
      if (.not. allocated(modmw)) allocate (modmw(nvpa, nmu, -nzgrid:nzgrid, nspec))
      if (.not. allocated(velvpamu)) allocate (velvpamu(nvpa, nmu, -nzgrid:nzgrid))

      ia = 1

      do is = 1, nspec
         do isb = 1, nspec

            massr = spec(is)%mass / spec(isb)%mass

            do iz = -nzgrid, nzgrid
               do iv = 1, nvpa
                  do imu = 1, nmu
                     x = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mu(imu))
                     Gf = (erf(x / sqrt(massr)) - x / sqrt(massr) * (2 / sqrt(pi)) * exp(-x**2 / massr)) / (2 * x**2 / massr)
                     nuD(iv, imu, iz, is, isb) = scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(x / sqrt(massr)) - Gf) / x**3
                     nupa(iv, imu, iz, is, isb) = spec(is)%vnew(isb) * 2 * Gf / x**3
                     velvpamu(iv, imu, iz) = x
                     mw(iv, imu, iz, is) = maxwell_vpa(iv, is) * maxwell_mu(1, iz, imu, is)
                  end do
               end do
            end do

            ! electron-ion collisions
            ! approximation of Lorentz operator using mass-ratio expansion
            if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
               do iz = -nzgrid, nzgrid
                  do iv = 1, nvpa
                     do imu = 1, nmu
                        x = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mu(imu))
                        nuD(iv, imu, iz, is, isb) = scale_pitch_angle_scattering * spec(is)%vnew(isb) * 1./x**3
                     end do
                  end do
               end do
            end if

         end do
      end do

      ! get a function with vanishing energy moment, modmw
      do is = 1, nspec
         do iz = -nzgrid, nzgrid
            call integrate_vmu(velvpamu(:, :, iz)**2 * mw(:, :, iz, is), iz, v2mwint(iz))
            call integrate_vmu(velvpamu(:, :, iz)**4 * mw(:, :, iz, is), iz, v4mwint(iz))
            modmw(:, :, iz, is) = mw(:, :, iz, is) - velvpamu(:, :, iz)**2 * v2mwint(iz) / v4mwint(iz) * mw(:, :, iz, is)
         end do
      end do

      nux = scale_mixed_derivative_terms * (nupa - scale_pitch_angle_scattering * nuD)

      if (nspec > 1) then
         ! scale_energy_diffusion_ei controls e-i energy diffusion, note that it is also used in blockmatrix
         nux(:, :, :, 2, 1) = scale_mixed_derivative_terms * (scale_energy_diffusion_ei * nupa(:, :, :, 2, 1) - scale_pitch_angle_scattering_ei * scale_pitch_angle_scattering * nuD(:, :, :, 2, 1))
         nuD(:, :, :, 2, 1) = scale_pitch_angle_scattering_ei * nuD(:, :, :, 2, 1)
      end if

   end subroutine init_nusDpa

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_fp_diffmatrix

      use grids_time, only: code_dt
      use grids_species, only: nspec, spec
      use grids_velocity, only: dvpa, vpa, nvpa, mu, nmu, maxwell_mu, maxwell_vpa, dmu
      use parallelisation_layouts, only: kxkyzt_lo
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use parallelisation_layouts, only: idx
      use geometry, only: bmag
      use arrays, only: kperp2
      use constants, only: pi
      use common_types, only: spec_type
      use spfunc, only: erf => erf_ext
      use mp, only: sum_allreduce
      use banded_lu, only: band_ldab, band_row, band_lu_factorise, band_lu_factorise_threaded

      implicit none

      integer :: ikxkyzt, iky, ikx, iz, is, isb
      integer :: imu, ia, iv, ivv, imu2
      integer :: nc, nb, bm_colind, bm_rowind
      integer :: ikxkyzt_first, ikxkyzt_last, npoints_chunk
      logical :: factorise_threaded
      real :: vpap, vpam, mum, mup
      real :: xpv, xmv, nupapv, nupamv, nuDpv, nuDmv, mwpv, mwmv, gam_mu, gam_mum, gam_mup
      real :: mwm, mwp, nuDm, nuDp, nupam, nupap, xm, xp
      real :: massr, eiediff, eidefl

      integer, parameter :: ion_species = 1
      integer, parameter :: electron_species = 2

      !-------------------------------------------------------------------------

      if (.not. allocated(cdiffmat_band)) &
         allocate (cdiffmat_band(band_ldab(nmu + 1, nmu + 1), nmu * nvpa, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec))
      if (.not. allocated(ipiv)) &
         allocate (ipiv(nvpa * nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec))

      ! The blocks aa_blcs, bb_blcs and cc_blcs are (nvpa,nmu,nmu) per point and species pair, i.e. larger than
      ! cdiffmat_band, so they are not stored for all points at once: they are calculated for a chunk of at most
      ! <npoints_chunk> points (about 128 MB of blocks), added to cdiffmat_band (and blockmatrix), and freed.
      ! The blocks keep the index ikxkyzt of the point, so the calculation below is the same for every chunk.
      npoints_chunk = max(1, int(128.e6 / (3. * 16. * nvpa * nmu**2 * nspec**2)))
      cdiffmat_band = 0.
      if ((exact_conservation_test_particle_matrix) .or. (field_particle_density_from_test_particle)) then
         if (.not. allocated(blockmatrix)) &
            allocate (blockmatrix(nvpa * nmu, nvpa * nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec, nspec))
         if (.not. allocated(blockmatrix_sum)) &
            allocate (blockmatrix_sum(nvpa * nmu, nvpa * nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec))
         blockmatrix = 0.
      end if
      ! exact_conservation only needs the test-particle operator applied to vpa F0 and v^2 F0 (see apply_testpart_matrix),
      ! which are taken from the blocks of every chunk below, without building blockmatrix
      if (exact_conservation_block_products) then
         if (.not. allocated(tp_f0moments)) allocate (tp_f0moments(nvpa * nmu, 2, nz_local, nspec, nspec))
         tp_f0moments = 0.
      end if

      do ikxkyzt_first = kxkyzt_lo%llim_proc, kxkyzt_lo%ulim_proc, npoints_chunk
      ikxkyzt_last = min(ikxkyzt_first + npoints_chunk - 1, kxkyzt_lo%ulim_proc)
      allocate (aa_blcs(nvpa, nmu, nmu, ikxkyzt_first:ikxkyzt_last, nspec, nspec))
      allocate (bb_blcs(nvpa, nmu, nmu, ikxkyzt_first:ikxkyzt_last, nspec, nspec))
      allocate (cc_blcs(nvpa, nmu, nmu, ikxkyzt_first:ikxkyzt_last, nspec, nspec))

      ! AVB: calculate the discretisation matrix -\Delta t C_test^{ab}
      ! because of mixed vpa-mu derivatives in the test particle operator
      ! this matrix is block tri-diagonal, with dimension nmu*nvpa x nmu*nvpa
      ! store and operate with the matrix in band format

      ! aa_blcs stores subdiagonal blocks, bb_blcs diagonal blocks and cc_blcs superdiagonal blocks
      ! aa_blcs(1,:,:) and cc_blcs(nvpa,:,:) are never used
      ! mu-derivatives are contained within blocks, thus blocks have dimension nmu x nmu

      ia = 1
      aa_blcs = 0.
      bb_blcs = 0.
      cc_blcs = 0.

      !$omp parallel do default(none) &
      !$omp   firstprivate(kxkyzt_lo, ikxkyzt_first, ikxkyzt_last, ia, nspec, nvpa, nmu) &
      !$omp   private(imu, iv, ikxkyzt, iky, ikx, iz, is, isb, &
      !$omp           vpap, vpam, mum, mup, xpv, xmv, &
      !$omp           nupapv, nupamv, nuDpv, nuDmv, mwpv, mwmv, &
      !$omp           gam_mu, gam_mum, gam_mup, mwm, mwp, &
      !$omp           nuDm, nuDp, nupam, nupap, xm, xp, &
      !$omp           massr, eiediff, eidefl) &
      !$omp   shared(aa_blcs, bb_blcs, cc_blcs, nux, mw, maxwell_mu, maxwell_vpa, &
      !$omp          vpa, mu, bmag, dmu, dvpa, spec, code_dt, &
      !$omp          spitzer_problem, conservative_test_particle_density, include_mu_derivatives,nupa, &
      !$omp          mass_ratio_approximation_ei, scale_energy_diffusion_ei, scale_pitch_angle_scattering_ei, scale_pitch_angle_scattering, nuD)
      do imu = 1, nmu
         do iv = 1, nvpa
            do ikxkyzt = ikxkyzt_first, ikxkyzt_last
               iky = iky_idx(kxkyzt_lo, ikxkyzt)
               ikx = ikx_idx(kxkyzt_lo, ikxkyzt)
               iz = iz_idx(kxkyzt_lo, ikxkyzt)

               do is = 1, nspec

               if (spitzer_problem) then
                  if (.not. (spec(is)%type == electron_species)) cycle ! add eon-eon and eon-ion collisions only for Spitzer problem
               end if

               do isb = 1, nspec
                  ! the factors on the e-i energy diffusion and pitch-angle scattering (e.g. scale_energy_diffusion_ei = 0
                  ! for the Spitzer problem); electrons are species 2 and ions species 1, see init_collisions_fp
                  if ((is == 2) .and. (isb == 1)) then
                     eiediff = scale_energy_diffusion_ei
                     eidefl = scale_pitch_angle_scattering_ei
                  else
                     eiediff = 1.
                     eidefl = 1.
                  end if

                  massr = spec(is)%mass / spec(isb)%mass

                  if (iv == 1) then
                     vpap = 0.5 * (vpa(iv) + vpa(iv + 1))
                     mwpv = exp(-vpap**2) * maxwell_mu(1, iz, imu, is)
                     xpv = sqrt(vpap**2 + 2 * bmag(ia, iz) * mu(imu))
                     nupapv = spec(is)%vnew(isb) * 2 * (erf(xpv / sqrt(massr)) &
                                          - xpv / sqrt(massr) * (2 / sqrt(pi)) * exp(-(xpv / sqrt(massr))**2)) / (2 * (xpv / sqrt(massr))**2) / xpv**3

                     if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
                        nuDpv = spec(is)%vnew(isb) / xpv**3
                     else
                        nuDpv = spec(is)%vnew(isb) * (erf(xpv / sqrt(massr)) - (erf(xpv / sqrt(massr)) &
                                       - (xpv / sqrt(massr)) * (2 / sqrt(pi)) * exp(-(xpv / sqrt(massr))**2)) / (2 * (xpv / sqrt(massr))**2)) / xpv**3
                     end if

                     if (imu == 1) then
                        ! using ghost cell at mu=0:
                        cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
       - 0.5*code_dt*vpa(iv+1)*0.5*(mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)+mu(imu+1)*nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))&
                                                                 / mw(iv + 1, imu + 1, iz, is) / (dvpa) / dmu(imu)
                        cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv / dvpa**2 / mw(iv + 1, imu, iz, is) &
       + 0.5*code_dt*vpa(iv+1)*0.5*(mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)+mu(imu+1)*nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))&
                                                             / mw(iv + 1, imu, iz, is) / (dvpa) / dmu(imu)
                        bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                 - 0.5*code_dt*vpa(iv)*0.5*(mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)+mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                                                 / mw(iv, imu + 1, iz, is) / (dvpa) / dmu(imu)
                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                 + 0.5*code_dt*vpa(iv)*0.5*(mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)+mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                                             / mw(iv, imu, iz, is) / (dvpa) / dmu(imu)
                        !
                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
      + code_dt * (0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv) / dvpa**2 / mw(iv, imu, iz, is)

                        ! mu operator
                        if (include_mu_derivatives) then
                           ! use ghost cell at mu_{0} = 0
                           mup = 0.5 * (mu(imu) + mu(imu + 1))
                           mwp = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mup)
                           xp = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mup)
                           if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
                              nuDp = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xp**3
                           else
                              nuDp = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xp / sqrt(massr)) &
                                    - (erf(xp / sqrt(massr)) - xp / sqrt(massr) * (2 / sqrt(pi)) * exp(-xp**2 / massr)) / (2 * xp**2 / massr)) / xp**3
                           end if
      nupap = spec(is)%vnew(isb) * 2 * (erf(xp / sqrt(massr)) - xp / sqrt(massr) * (2 / sqrt(pi)) * exp(-xp**2 / massr)) / (2 * xp**2 / massr) / xp**3
      gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*scale_pitch_angle_scattering*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                           gam_mup = 2 * (eiediff * nupap * mup**2 + eidefl * scale_pitch_angle_scattering * nuDp * vpa(iv)**2 / (2 * bmag(ia, iz)) * mup) * mwp
                   bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) - code_dt * (gam_mu*(-1) / dmu(imu) * dmu(imu) / 2./mu(imu) &
                         + (gam_mup*(-1) / dmu(imu) - gam_mu*(-1) / dmu(imu)) * mu(imu) / (dmu(imu) / 2.)) / mw(iv, imu, iz, is) / (dmu(imu) / 2.+mu(imu))
          bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) - code_dt * (gam_mu * 1 / dmu(imu) * dmu(imu) / 2./mu(imu) &
                   + (gam_mup * 1 / dmu(imu) - gam_mu * 1 / dmu(imu)) * mu(imu) / (dmu(imu) / 2.)) / mw(iv, imu + 1, iz, is) / (dmu(imu) / 2.+mu(imu))
                           ! mixed derivative:
                           if (conservative_test_particle_density) then
                              ! to ensure density conservation, change discretisation of mixed derivative term at imu=1
                              ! from explicit routine: Dmuh1 = (vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz,is)*Dvpah &
                              cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                         - 0.5 * 1.0 * code_dt * (0.5 * (vpa(iv) * mu(imu) * nux(iv, imu, iz, is, isb) * mw(iv, imu, iz, is) + vpa(iv + 1) * mu(imu) &
                                     * nux(iv + 1, imu, iz, is, isb) * mw(iv + 1, imu, iz, is)) * 1 / mw(iv + 1, imu, iz, is)) / (2 * dmu(imu)) / dvpa
                              cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
         - 0.5 * 1.0 * code_dt * (0.5 * (vpa(iv) * mu(imu + 1) * nux(iv, imu + 1, iz, is, isb) * mw(iv, imu + 1, iz, is) + vpa(iv + 1) * mu(imu + 1) &
                         * nux(iv + 1, imu + 1, iz, is, isb) * mw(iv + 1, imu + 1, iz, is)) * 1 / mw(iv + 1, imu + 1, iz, is)) / (2 * dmu(imu)) / dvpa
                              bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                         + 0.5 * 1.0 * code_dt * (0.5 * (vpa(iv) * mu(imu) * nux(iv, imu, iz, is, isb) * mw(iv, imu, iz, is) + vpa(iv + 1) * mu(imu) &
                                         * nux(iv + 1, imu, iz, is, isb) * mw(iv + 1, imu, iz, is)) * 1 / mw(iv, imu, iz, is)) / (2 * dmu(imu)) / dvpa
                              bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
         + 0.5 * 1.0 * code_dt * (0.5 * (vpa(iv) * mu(imu + 1) * nux(iv, imu + 1, iz, is, isb) * mw(iv, imu + 1, iz, is) + vpa(iv + 1) * mu(imu + 1) &
                             * nux(iv + 1, imu + 1, iz, is, isb) * mw(iv + 1, imu + 1, iz, is)) * 1 / mw(iv, imu + 1, iz, is)) / (2 * dmu(imu)) / dvpa
                           else
                              cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                      - 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv+1,imu,iz,is)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu)))&
                                                                   / (mu(imu) + dmu(imu)) / dvpa
                              cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                - 1.0*code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv+1,imu+1,iz,is)* mu(imu)/dmu(imu))&
                                                                       / (mu(imu) + dmu(imu)) / dvpa
                              bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                        + 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv,imu,iz,is)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu)))&
                                                                   / (mu(imu) + dmu(imu)) / dvpa
                              bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                  + 1.0*code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv,imu+1,iz,is)* mu(imu)/dmu(imu))&
                                                                       / (mu(imu) + dmu(imu)) / dvpa
                           end if
                        end if
                     else if (imu == nmu) then
                        ! AVB: one-sided difference for mu-derivative at imu=nmu:
                        if (conservative_test_particle_density) then
                           cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
       + 0.5*code_dt*vpa(iv+1)*0.5*(mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is)+mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is))&
                                                                    / mw(iv + 1, imu - 1, iz, is) / (dvpa) / dmu(nmu - 1)
                           cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv / dvpa**2 / mw(iv + 1, imu, iz, is) &
       - 0.5*code_dt*vpa(iv+1)*0.5*(mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is)+mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is))&
                                                                / mw(iv + 1, imu, iz, is) / (dvpa) / dmu(nmu - 1)
                           bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                 + 0.5*code_dt*vpa(iv)*0.5*(mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is))&
                                                                    / mw(iv, imu - 1, iz, is) / (dvpa) / dmu(nmu - 1)
                           bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                 - 0.5*code_dt*vpa(iv)*0.5*(mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is))&
                                                                / mw(iv, imu, iz, is) / (dvpa) / dmu(nmu - 1)
                        else
                           cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                          + 1.0*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                           cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv / dvpa**2 / mw(iv + 1, imu, iz, is) &
                            - 1.0*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu,iz,is) / (2*dvpa) / dmu(nmu-1)
                        end if
                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
      + code_dt * (0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv) / dvpa**2 / mw(iv, imu, iz, is)
                        ! mu operator
                        if (include_mu_derivatives) then
                           mum = 0.5 * (mu(imu) + mu(imu - 1))
                           mwm = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mum)
                           xm = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mum)
                           if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
                              nuDm = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xm**3
                           else
                              nuDm = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xm / sqrt(massr)) &
                                    - (erf(xm / sqrt(massr)) - xm / sqrt(massr) * (2 / sqrt(pi)) * exp(-xm**2 / massr)) / (2 * xm**2 / massr)) / xm**3
                           end if
      nupam = spec(is)%vnew(isb) * 2 * (erf(xm / sqrt(massr)) - xm / sqrt(massr) * (2 / sqrt(pi)) * exp(-xm**2 / massr)) / (2 * xm**2 / massr) / xm**3
      gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*scale_pitch_angle_scattering*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                           gam_mum = 2 * (eiediff * nupam * mum**2 + eidefl * scale_pitch_angle_scattering * nuDm * vpa(iv)**2 / (2 * bmag(ia, iz)) * mum) * mwm
                           if (conservative_test_particle_density) then
                              ! to ensure density conservation we assume that the argument of the outer derivative vanishes at nmu+1/2, ie
                              ! d/dmu[...]_{nmu} = ([...]_{nmu+1/2} - [...]_{nmu-1/2})/(dmu_{nmu-1}/2+dmu(nmu-1)/2)
                              ! where [...]_{nmu+1/2} is a ghost cell at mu = mu_{nmu} + dmu(nmu-1)/2, with [...]_{nmu+1/2} = 0.
        bb_blcs(iv,imu,imu,ikxkyzt,is,isb)  = bb_blcs(iv,imu,imu,ikxkyzt,is,isb) - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)) &
                               + (-gam_mu / dmu(imu - 1)) * dmu(imu - 1) / dmu(imu - 1)) / mw(iv, imu, iz, is) / (dmu(imu - 1) / 2.+dmu(imu - 1) / 2.)
bb_blcs(iv,imu,imu-1,ikxkyzt,is,isb)= bb_blcs(iv,imu,imu-1,ikxkyzt,is,isb) - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * (-1)/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)) &
                           - gam_mu*(-1) / dmu(imu - 1) * dmu(imu - 1) / dmu(imu - 1)) / mw(iv, imu - 1, iz, is) / (dmu(imu - 1) / 2.+dmu(imu - 1) / 2.)
                           else
     bb_blcs(iv,imu,imu,ikxkyzt,is,isb)  = bb_blcs(iv,imu,imu,ikxkyzt,is,isb) - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                 + (-gam_mu / dmu(imu - 1)) * dmu(imu - 1) / 2./dmu(imu - 1)) / mw(iv, imu, iz, is) / (dmu(imu - 1) / 2.+dmu(imu - 1))
                                     bb_blcs(iv,imu,imu-1,ikxkyzt,is,isb)= bb_blcs(iv,imu,imu-1,ikxkyzt,is,isb) - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * (-1)/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                             - gam_mu*(-1) / dmu(imu - 1) * dmu(imu - 1) / 2./dmu(imu - 1)) / mw(iv, imu - 1, iz, is) / (dmu(imu - 1) / 2.+dmu(imu - 1))
                           end if
                           ! mixed derivative:
                           if (conservative_test_particle_density) then
                              cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)+vpa(iv+1)*mu(imu  )*nux(iv+1,imu  ,iz,is,isb)*mw(iv+1,imu  ,iz,is))&
                                                                                      * 1 / mw(iv + 1, imu, iz, is)) / (2 * dmu(imu - 1)) / dvpa
                              cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+vpa(iv+1)*mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is))&
                                                                                        * 1 / mw(iv + 1, imu - 1, iz, is)) / (2 * dmu(imu - 1)) / dvpa
                              bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)+vpa(iv+1)*mu(imu  )*nux(iv+1,imu  ,iz,is,isb)*mw(iv+1,imu  ,iz,is))&
                                                                                      * 1 / mw(iv, imu, iz, is)) / (2 * dmu(imu - 1)) / dvpa
                              bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+vpa(iv+1)*mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is))&
                                                                                          * 1 / mw(iv, imu - 1, iz, is)) / (2 * dmu(imu - 1)) / dvpa
                           else
                              cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
            - 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv+1,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1)))&
                                                                   / (dmu(imu - 1) + dmu(imu - 1)) / (dvpa)
                              cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                           + 1.0*code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv+1,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1))&
                                                                       / (dmu(imu - 1) + dmu(imu - 1)) / (dvpa)
                              bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
              + 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1)))&
                                                                   / (dmu(imu - 1) + dmu(imu - 1)) / (dvpa)
                              bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                             - 1.0*code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1))&
                                                                       / (dmu(imu - 1) + dmu(imu - 1)) / (dvpa)
                           end if

                        end if
                     else
                        ! interior mu points for iv = 1
                        ! use ghost cell for derivative
                        ! d/dvpa [vpa*mu*nux*F0*dhdmu] = [ ]_nvpa+0.5-[ ]_nvpa-0.5 / dvpa = -[ ]_nvpa-0.5 / dvpa
                        ! = -0.5*([ ]_nvpa-1 + [ ]_nvpa) / dvpa
                        ! could be cleared up, by moving non-nux part of aa_blcs down to non-nux part of bb_bcls
                        cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                + 0.5*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu-1,iz,is)*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                        cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                - 0.5*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu+1,iz,is)*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                        cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv / dvpa**2 / mw(iv + 1, imu, iz, is) &
             - 0.5*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                        bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
     + 0.5*code_dt*vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)/mw(iv,imu-1,iz,is)*dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                        bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
     - 0.5*code_dt*vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)/mw(iv,imu+1,iz,is)*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                 - 0.5*code_dt*vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)

                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
      + code_dt * (0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv) / dvpa**2 / mw(iv, imu, iz, is)

                        ! mu operator
                        if (include_mu_derivatives) then
                           ! quantities at mu_{i+1/2} and mu_{i-1/2}:
                           mup = 0.5 * (mu(imu) + mu(imu + 1))
                           mum = 0.5 * (mu(imu) + mu(imu - 1))
                           mwp = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mup)
                           mwm = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mum)
                           xp = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mup)
                           xm = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mum)
                           if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
                              nuDp = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xp**3
                              nuDm = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xm**3
                           else
                              nuDp = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xp / sqrt(massr)) &
                                    - (erf(xp / sqrt(massr)) - xp / sqrt(massr) * (2 / sqrt(pi)) * exp(-xp**2 / massr)) / (2 * xp**2 / massr)) / xp**3
                              nuDm = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xm / sqrt(massr)) &
                                    - (erf(xm / sqrt(massr)) - xm / sqrt(massr) * (2 / sqrt(pi)) * exp(-xm**2 / massr)) / (2 * xm**2 / massr)) / xm**3
                           end if
      nupap = spec(is)%vnew(isb) * 2 * (erf(xp / sqrt(massr)) - xp / sqrt(massr) * (2 / sqrt(pi)) * exp(-xp**2 / massr)) / (2 * xp**2 / massr) / xp**3
      nupam = spec(is)%vnew(isb) * 2 * (erf(xm / sqrt(massr)) - xm / sqrt(massr) * (2 / sqrt(pi)) * exp(-xm**2 / massr)) / (2 * xm**2 / massr) / xm**3
      gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*scale_pitch_angle_scattering*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                           gam_mum = 2 * (eiediff * nupam * mum**2 + eidefl * scale_pitch_angle_scattering * nuDm * vpa(iv)**2 / (2 * bmag(ia, iz)) * mum) * mwm
                           gam_mup = 2 * (eiediff * nupap * mup**2 + eidefl * scale_pitch_angle_scattering * nuDp * vpa(iv)**2 / (2 * bmag(ia, iz)) * mup) * mwp

                           bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                  - code_dt*((gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mum / dmu(imu-1))*dmu(imu)/dmu(imu-1) &
       + (-gam_mu * (dmu(imu) / dmu(imu - 1) - dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu)) - gam_mup / dmu(imu)) * dmu(imu - 1) / dmu(imu)) &
                                                                / mw(iv, imu, iz, is) * 2./(dmu(imu - 1) + dmu(imu))
                           bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                - code_dt * ((gam_mu*(-1) * dmu(imu) / dmu(imu - 1) / (dmu(imu - 1) + dmu(imu)) - gam_mum*(-1) / dmu(imu - 1)) * dmu(imu) / dmu(imu - 1) &
 - gam_mu*(-1) * dmu(imu) / dmu(imu - 1) / (dmu(imu - 1) + dmu(imu)) * dmu(imu - 1) / dmu(imu)) / mw(iv, imu - 1, iz, is) * 2./(dmu(imu - 1) + dmu(imu))
                           bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                                 - code_dt * (gam_mu * dmu(imu - 1) / dmu(imu) / (dmu(imu - 1) + dmu(imu)) * dmu(imu) / dmu(imu - 1) &
                                    + (gam_mup / dmu(imu) - gam_mu * dmu(imu - 1) / dmu(imu) / (dmu(imu - 1) + dmu(imu))) * dmu(imu - 1) / dmu(imu)) &
                                                                    / mw(iv, imu + 1, iz, is) * 2 / (dmu(imu - 1) + dmu(imu))
                           ! mixed derivative:
                           cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                    - 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)+vpa(iv+1)*mu(imu  )*nux(iv+1,imu  ,iz,is,isb)*mw(iv+1,imu  ,iz,is))&
                             * 1 / mw(iv + 1, imu, iz, is) * (dmu(imu) / dmu(imu - 1) - dmu(imu - 1) / dmu(imu))) / (dmu(imu - 1) + dmu(imu)) / (dvpa)
                           cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                    + 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+vpa(iv+1)*mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is))&
                                                     * 1 / mw(iv + 1, imu - 1, iz, is) * dmu(imu) / dmu(imu - 1)) / (dmu(imu - 1) + dmu(imu)) / (dvpa)
                           cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                    - 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)+vpa(iv+1)*mu(imu+1)*nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))&
                                                     * 1 / mw(iv + 1, imu + 1, iz, is) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu)) / (dvpa)
                           bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                    + 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is)+vpa(iv+1)*mu(imu  )*nux(iv+1,imu  ,iz,is,isb)*mw(iv+1,imu  ,iz,is))&
                                 * 1 / mw(iv, imu, iz, is) * (dmu(imu) / dmu(imu - 1) - dmu(imu - 1) / dmu(imu))) / (dmu(imu - 1) + dmu(imu)) / (dvpa)
                           bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                    - 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+vpa(iv+1)*mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is))&
                                                         * 1 / mw(iv, imu - 1, iz, is) * dmu(imu) / dmu(imu - 1)) / (dmu(imu - 1) + dmu(imu)) / (dvpa)
                           bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                    + 0.5*1.0*code_dt*(0.5*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)+vpa(iv+1)*mu(imu+1)*nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))&
                                                         * 1 / mw(iv, imu + 1, iz, is) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu)) / (dvpa)
                        end if
                     end if

                  else if (iv == nvpa) then
                     vpam = 0.5 * (vpa(iv) + vpa(iv - 1))
                     mwmv = exp(-vpam**2) * maxwell_mu(1, iz, imu, is)
                     xmv = sqrt(vpam**2 + 2 * bmag(ia, iz) * mu(imu))
                   nupamv = spec(is)%vnew(isb)*2*(erf(xmv/sqrt(massr))-xmv/sqrt(massr)*(2/sqrt(pi))*exp(-xmv**2/massr)) / (2*xmv**2/massr)/xmv**3
                     if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
                        nuDmv = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xmv**3
                     else
                        nuDmv = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xmv / sqrt(massr)) &
                               - (erf(xmv / sqrt(massr)) - xmv / sqrt(massr) * (2 / sqrt(pi)) * exp(-xmv**2 / massr)) / (2 * xmv**2 / massr)) / xmv**3
                     end if
                     if (imu == 1) then
                        ! one-sided difference for mu-derivative at imu=1:
                        ! ie d/dvpa [vpa*mu*nux*F0*dhdmu] = [ ]_nvpa+0.5-[ ]_nvpa-0.5 / dvpa = -[ ]_nvpa-0.5 / dvpa
                        ! = -0.5*([ ]_nvpa-1 + [ ]_nvpa) / dvpa
                        ! could be cleared up, by moving non-nux part of aa_blcs down to non-nux part of bb_bcls
                        aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
       + 0.5*code_dt*vpa(iv-1)*0.5*(mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)+mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is))&
                                                                 / mw(iv - 1, imu + 1, iz, is) / (dvpa) / dmu(imu)
                        aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) - code_dt * 0.5 * (eiediff * nupamv * vpam**2 &
                                                + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv / dvpa**2 / mw(iv - 1, imu, iz, is) &
       - 0.5*code_dt*vpa(iv-1)*0.5*(mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)+mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is))&
                                                             / mw(iv - 1, imu, iz, is) / (dvpa) / dmu(imu)
                        bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                 + 0.5*code_dt*vpa(iv)*0.5*(mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)+mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                                                 / mw(iv, imu + 1, iz, is) / (dvpa) / dmu(imu)
                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                 - 0.5*code_dt*vpa(iv)*0.5*(mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)+mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                                             / mw(iv, imu, iz, is) / (dvpa) / dmu(imu)

                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
      + code_dt * (0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv) / dvpa**2 / mw(iv, imu, iz, is)

                        ! mu operator
                        if (include_mu_derivatives) then
                           ! use ghost cell at mu_{0} = 0
                           mup = 0.5 * (mu(imu) + mu(imu + 1))
                           mwp = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mup)
                           xp = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mup)
                           if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
                              nuDp = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xp**3
                           else
                              nuDp = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xp / sqrt(massr)) &
                                    - (erf(xp / sqrt(massr)) - xp / sqrt(massr) * (2 / sqrt(pi)) * exp(-xp**2 / massr)) / (2 * xp**2 / massr)) / xp**3
                           end if
      nupap = spec(is)%vnew(isb) * 2 * (erf(xp / sqrt(massr)) - xp / sqrt(massr) * (2 / sqrt(pi)) * exp(-xp**2 / massr)) / (2 * xp**2 / massr) / xp**3
      gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*scale_pitch_angle_scattering*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                           gam_mup = 2 * (eiediff * nupap * mup**2 + eidefl * scale_pitch_angle_scattering * nuDp * vpa(iv)**2 / (2 * bmag(ia, iz)) * mup) * mwp

                   bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) - code_dt * (gam_mu*(-1) / dmu(imu) * dmu(imu) / 2./mu(imu) &
                         + (gam_mup*(-1) / dmu(imu) - gam_mu*(-1) / dmu(imu)) * mu(imu) / (dmu(imu) / 2.)) / mw(iv, imu, iz, is) / (dmu(imu) / 2.+mu(imu))
          bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) - code_dt * (gam_mu * 1 / dmu(imu) * dmu(imu) / 2./mu(imu) &
                   + (gam_mup * 1 / dmu(imu) - gam_mu * 1 / dmu(imu)) * mu(imu) / (dmu(imu) / 2.)) / mw(iv, imu + 1, iz, is) / (dmu(imu) / 2.+mu(imu))
                           ! mixed derivative:
                           if (conservative_test_particle_density) then
                              ! to ensure density conservation, change discretisation of mixed derivative term at imu=1
                              ! from explicit routine: Dmuh1 = (vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz,is)*Dvpah + vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz)*mw(iv,imu+1,iz,is)*Dvpah_p) / (2.*dmu(imu))
                              aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu  )*nux(iv-1,imu  ,iz,is,isb)*mw(iv-1,imu  ,iz,is)+vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is))&
                                                                                      * 1 / mw(iv - 1, imu, iz, is)) / (2 * dmu(imu)) / (dvpa)
                              aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is)+vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                                                                          * 1 / mw(iv - 1, imu + 1, iz, is)) / (2 * dmu(imu)) / (dvpa)
                              bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu  )*nux(iv-1,imu  ,iz,is,isb)*mw(iv-1,imu  ,iz,is)+vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is))&
                                                                                      * 1 / mw(iv, imu, iz, is)) / (2 * dmu(imu)) / (dvpa)
                              bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is)+vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                                                                          * 1 / mw(iv, imu + 1, iz, is)) / (2 * dmu(imu)) / (dvpa)
                           else
                              aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                      + 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv-1,imu,iz,is)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu)))&
                                                                   / (mu(imu) + dmu(imu)) / (dvpa)
                              aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                + 1.0*code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv-1,imu+1,iz,is)* mu(imu)/dmu(imu))&
                                                                       / (mu(imu) + dmu(imu)) / (dvpa)
                              bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                        - 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv,imu,iz,is)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu)))&
                                                                   / (mu(imu) + dmu(imu)) / (dvpa)
                              bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                  - 1.0*code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv,imu+1,iz,is)* mu(imu)/dmu(imu))&
                                                                       / (mu(imu) + dmu(imu)) / (dvpa)
                           end if
                        end if
                     else if (imu == nmu) then
                        ! one-sided difference for mu-derivative at imu=nmu:
                        if (conservative_test_particle_density) then
                           ! take vpa derivative using ghost cell at nvpa+0.5
                           ! ie d/dvpa [vpa*mu*nux*F0*dhdmu] = [ ]_nvpa+0.5-[ ]_nvpa-0.5 / dvpa = -[ ]_nvpa-0.5 / dvpa
                           ! = -0.5*([ ]_nvpa-1 + [ ]_nvpa) / dvpa
                           ! note this could be cleared up, by moving the non-nux part of aa_blcs down to the non-nux part of bb_bcls
                           aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
       - 0.5*code_dt*vpa(iv-1)*0.5*(mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is))&
                                                                    / mw(iv - 1, imu - 1, iz, is) / (dvpa) / dmu(nmu - 1)
                           aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv / dvpa**2 / mw(iv - 1, imu, iz, is) &
       + 0.5*code_dt*vpa(iv-1)*0.5*(mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is))&
                                                                / mw(iv - 1, imu, iz, is) / (dvpa) / dmu(nmu - 1)
                           bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                 - 0.5*code_dt*vpa(iv)*0.5*(mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is))&
                                                                    / mw(iv, imu - 1, iz, is) / (dvpa) / dmu(nmu - 1)
                           bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                 + 0.5*code_dt*vpa(iv)*0.5*(mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)+mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is))&
                                                                / mw(iv, imu, iz, is) / (dvpa) / dmu(nmu - 1)
                        else
                           aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                          - 1.0*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                           aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv / dvpa**2 / mw(iv - 1, imu, iz, is) &
                            + 1.0*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu,iz,is) / (2*dvpa) / dmu(nmu-1)
                        end if
                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
      + code_dt * (0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv) / dvpa**2 / mw(iv, imu, iz, is)

                        ! mu operator
                        if (include_mu_derivatives) then
                           mum = 0.5 * (mu(imu) + mu(imu - 1))
                           mwm = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mum)
                           xm = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mum)
                           if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
                              nuDm = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xm**3
                           else
                              nuDm = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xm / sqrt(massr)) &
                                    - (erf(xm / sqrt(massr)) - xm / sqrt(massr) * (2 / sqrt(pi)) * exp(-xm**2 / massr)) / (2 * xm**2 / massr)) / xm**3
                           end if
      nupam = spec(is)%vnew(isb) * 2 * (erf(xm / sqrt(massr)) - xm / sqrt(massr) * (2 / sqrt(pi)) * exp(-xm**2 / massr)) / (2 * xm**2 / massr) / xm**3
     gam_mu  = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*scale_pitch_angle_scattering*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                           gam_mum = 2 * (eiediff * nupam * mum**2 + eidefl * scale_pitch_angle_scattering * nuDm * vpa(iv)**2 / (2 * bmag(ia, iz)) * mum) * mwm
                           if (conservative_test_particle_density) then
                              ! to ensure density conservation we assume that the argument of the outer derivative vanishes at nmu+1/2, ie
                              ! d/dmu[...]_{nmu} = ([...]_{nmu+1/2} - [...]_{nmu-1/2})/(dmu_{nmu-1}/2+dmu(nmu-1)/2), where [...]_{nmu+1/2} is a ghost cell at mu = mu_{nmu} + dmu(nmu-1)/2, with [...]_{nmu+1/2} = 0.
         bb_blcs(iv,imu,imu,ikxkyzt,is,isb) = bb_blcs(iv,imu,imu,ikxkyzt,is,isb) - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)) &
                               + (-gam_mu / dmu(imu - 1)) * dmu(imu - 1) / dmu(imu - 1)) / mw(iv, imu, iz, is) / (dmu(imu - 1) / 2.+dmu(imu - 1) / 2.)
                                     bb_blcs(iv,imu,imu-1,ikxkyzt,is,isb) = bb_blcs(iv,imu,imu-1,ikxkyzt,is,isb) - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * (-1)/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)) &
                           - gam_mu*(-1) / dmu(imu - 1) * dmu(imu - 1) / dmu(imu - 1)) / mw(iv, imu - 1, iz, is) / (dmu(imu - 1) / 2.+dmu(imu - 1) / 2.)
                           else
      bb_blcs(iv,imu,imu,ikxkyzt,is,isb) = bb_blcs(iv,imu,imu,ikxkyzt,is,isb) - code_dt*((gam_mu/dmu(imu-1) - gam_mum/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                                 + (-gam_mu / dmu(imu - 1)) * dmu(imu - 1) / 2./dmu(imu - 1)) / mw(iv, imu, iz, is) / (dmu(imu - 1) / 2.+dmu(imu - 1))
                                     bb_blcs(iv,imu,imu-1,ikxkyzt,is,isb) = bb_blcs(iv,imu,imu-1,ikxkyzt,is,isb) - code_dt*((-gam_mu/dmu(imu-1) - gam_mum * (-1)/dmu(imu-1))*dmu(imu-1)/(dmu(imu-1)/2.) &
                             - gam_mu*(-1) / dmu(imu - 1) * dmu(imu - 1) / 2./dmu(imu - 1)) / mw(iv, imu - 1, iz, is) / (dmu(imu - 1) / 2.+dmu(imu - 1))
                           end if
                           ! mixed derivative:
                           if (conservative_test_particle_density) then
                              ! 04.03. removed an extra bracket at end of lines in following block
                              aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is))&
                                                                                        * 1 / mw(iv - 1, imu - 1, iz, is)) / (2 * dmu(imu - 1)) / dvpa
                              aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                        - 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu  )*nux(iv-1,imu  ,iz,is,isb)*mw(iv-1,imu  ,iz,is)+vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is))&
                                                                                      * 1 / mw(iv - 1, imu, iz, is)) / (2 * dmu(imu - 1)) / dvpa
                              bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is))&
                                                                                          * 1 / mw(iv, imu - 1, iz, is)) / (2 * dmu(imu - 1)) / dvpa
                              bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                        + 0.5*1.0*code_dt*(0.5*(vpa(iv-1)*mu(imu  )*nux(iv-1,imu  ,iz,is,isb)*mw(iv-1,imu  ,iz,is)+vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is))&
                                                                                      * 1 / mw(iv, imu, iz, is)) / (2 * dmu(imu - 1)) / dvpa
                           else
                              aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
            + 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv-1,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1)))&
                                                                   / (dmu(imu - 1) + dmu(imu - 1)) / (dvpa)
                              aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                           - 1.0*code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv-1,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1))&
                                                                       / (dmu(imu - 1) + dmu(imu - 1)) / (dvpa)
                              bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
              - 1.0*code_dt*(vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1)))&
                                                                   / (dmu(imu - 1) + dmu(imu - 1)) / (dvpa)
                              bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                             + 1.0*code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1))&
                                                                       / (dmu(imu - 1) + dmu(imu - 1)) / (dvpa)
                           end if
                        end if
                     else
                        ! interior mu points for iv=nvpa
                        ! take vpa derivative using ghost cell at nvpa+0.5
                        ! ie d/dvpa [vpa*mu*nux*F0*dhdmu] = [ ]_nvpa+0.5-[ ]_nvpa-0.5 / dvpa = -[ ]_nvpa-0.5 / dvpa
                        ! = -0.5*([ ]_nvpa-1 + [ ]_nvpa) / dvpa

                        ! note this could be cleared up, by moving the non-nux part of aa_blcs down to the non-nux part of bb_bcls
                        aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                - 0.5*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu-1,iz,is) * dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                        aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                + 0.5*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu+1,iz,is) * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                        aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv / dvpa**2 / mw(iv - 1, imu, iz, is) &
             + 0.5*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                        bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
   - 0.5*code_dt*vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)/mw(iv,imu-1,iz,is) * dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                        bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
   + 0.5*code_dt*vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)/mw(iv,imu+1,iz,is) * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (dvpa)
                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                 + 0.5*code_dt*vpa(iv)*mu(imu)*nux(iv,imu,iz,is,isb)*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (dvpa)

                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
      + code_dt * (0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv) / dvpa**2 / mw(iv, imu, iz, is)

                        ! mu operator
                        if (include_mu_derivatives) then
                           ! quantities at mu_{i+1/2} and mu_{i-1/2}:
                           mup = 0.5 * (mu(imu) + mu(imu + 1))
                           mum = 0.5 * (mu(imu) + mu(imu - 1))
                           mwp = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mup)
                           mwm = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mum)
                           xp = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mup)
                           xm = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mum)
                           if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
                              nuDp = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xp**3
                              nuDm = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xm**3
                           else
                                     nuDp = eidefl*scale_pitch_angle_scattering*spec(is)%vnew(isb)*(erf(xp/sqrt(massr))-(erf(xp/sqrt(massr))-xp/sqrt(massr)*(2/sqrt(pi))*exp(-xp**2/massr)) / (2*xp**2/massr)) / xp**3
                                     nuDm = eidefl*scale_pitch_angle_scattering*spec(is)%vnew(isb)*(erf(xm/sqrt(massr))-(erf(xm/sqrt(massr))-xm/sqrt(massr)*(2/sqrt(pi))*exp(-xm**2/massr)) / (2*xm**2/massr)) / xm**3
                           end if
      nupap = spec(is)%vnew(isb) * 2 * (erf(xp / sqrt(massr)) - xp / sqrt(massr) * (2 / sqrt(pi)) * exp(-xp**2 / massr)) / (2 * xp**2 / massr) / xp**3
      nupam = spec(is)%vnew(isb) * 2 * (erf(xm / sqrt(massr)) - xm / sqrt(massr) * (2 / sqrt(pi)) * exp(-xm**2 / massr)) / (2 * xm**2 / massr) / xm**3
      gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*scale_pitch_angle_scattering*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                           gam_mum = 2 * (eiediff * nupam * mum**2 + eidefl * scale_pitch_angle_scattering * nuDm * vpa(iv)**2 / (2 * bmag(ia, iz)) * mum) * mwm
                           gam_mup = 2 * (eiediff * nupap * mup**2 + eidefl * scale_pitch_angle_scattering * nuDp * vpa(iv)**2 / (2 * bmag(ia, iz)) * mup) * mwp
                           bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                  - code_dt*((gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mum / dmu(imu-1))*dmu(imu)/dmu(imu-1) &
       + (-gam_mu * (dmu(imu) / dmu(imu - 1) - dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu)) - gam_mup / dmu(imu)) * dmu(imu - 1) / dmu(imu)) &
                                                                / mw(iv, imu, iz, is) * 2./(dmu(imu - 1) + dmu(imu))
                           bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                - code_dt * ((gam_mu*(-1) * dmu(imu) / dmu(imu - 1) / (dmu(imu - 1) + dmu(imu)) - gam_mum*(-1) / dmu(imu - 1)) * dmu(imu) / dmu(imu - 1) &
 - gam_mu*(-1) * dmu(imu) / dmu(imu - 1) / (dmu(imu - 1) + dmu(imu)) * dmu(imu - 1) / dmu(imu)) / mw(iv, imu - 1, iz, is) * 2./(dmu(imu - 1) + dmu(imu))
                           bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                                 - code_dt * (gam_mu * dmu(imu - 1) / dmu(imu) / (dmu(imu - 1) + dmu(imu)) * dmu(imu) / dmu(imu - 1) &
                                    + (gam_mup / dmu(imu) - gam_mu * dmu(imu - 1) / dmu(imu) / (dmu(imu - 1) + dmu(imu))) * dmu(imu - 1) / dmu(imu)) &
                                                                    / mw(iv, imu + 1, iz, is) * 2 / (dmu(imu - 1) + dmu(imu))
                           ! mixed derivative, one-sided difference in vpa at iv = nvpa:
                           ! use second order accurate treatment for vpa derivative at nvpa
                           aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
 + 0.5*code_dt*(0.5*(vpa(iv-1)*mu(imu  )*nux(iv-1,imu  ,iz,is,isb)*mw(iv-1,imu  ,iz,is)+vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is))&
                             * 1 / mw(iv - 1, imu, iz, is) * (dmu(imu) / dmu(imu - 1) - dmu(imu - 1) / dmu(imu))) / (dmu(imu - 1) + dmu(imu)) / (dvpa)
                           aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
 - 0.5*code_dt*(0.5*(vpa(iv-1)*mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is))&
                                                     * 1 / mw(iv - 1, imu - 1, iz, is) * dmu(imu) / dmu(imu - 1)) / (dmu(imu - 1) + dmu(imu)) / (dvpa)
                           aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
 + 0.5*code_dt*(0.5*(vpa(iv-1)*mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is)+vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                                     * 1 / mw(iv - 1, imu + 1, iz, is) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu)) / (dvpa)
                           bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
 - 0.5*code_dt*(0.5*(vpa(iv-1)*mu(imu  )*nux(iv-1,imu  ,iz,is,isb)*mw(iv-1,imu  ,iz,is)+vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu  ,iz,is))&
                                 * 1 / mw(iv, imu, iz, is) * (dmu(imu) / dmu(imu - 1) - dmu(imu - 1) / dmu(imu))) / (dmu(imu - 1) + dmu(imu)) / (dvpa)
                           bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
 + 0.5*code_dt*(0.5*(vpa(iv-1)*mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is))&
                                                         * 1 / mw(iv, imu - 1, iz, is) * dmu(imu) / dmu(imu - 1)) / (dmu(imu - 1) + dmu(imu)) / (dvpa)
                           bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
 - 0.5*code_dt*(0.5*(vpa(iv-1)*mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is)+vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is))&
                                                         * 1 / mw(iv, imu + 1, iz, is) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu)) / (dvpa)
                        end if
                     end if

                  else ! interior vpa points
                     vpam = 0.5 * (vpa(iv) + vpa(iv - 1))
                     vpap = 0.5 * (vpa(iv) + vpa(iv + 1))
                     mwmv = exp(-vpam**2) * maxwell_mu(1, iz, imu, is)
                     mwpv = exp(-vpap**2) * maxwell_mu(1, iz, imu, is)
                     xpv = sqrt(vpap**2 + 2 * bmag(ia, iz) * mu(imu))
                     xmv = sqrt(vpam**2 + 2 * bmag(ia, iz) * mu(imu))
                   nupamv = spec(is)%vnew(isb)*2*(erf(xmv/sqrt(massr))-xmv/sqrt(massr)*(2/sqrt(pi))*exp(-xmv**2/massr)) / (2*xmv**2/massr)/xmv**3
                   nupapv = spec(is)%vnew(isb)*2*(erf(xpv/sqrt(massr))-xpv/sqrt(massr)*(2/sqrt(pi))*exp(-xpv**2/massr)) / (2*xpv**2/massr)/xpv**3
                     if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
                        nuDmv = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xmv**3
                        nuDpv = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xpv**3
                     else
                        nuDmv = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xmv / sqrt(massr)) &
                               - (erf(xmv / sqrt(massr)) - xmv / sqrt(massr) * (2 / sqrt(pi)) * exp(-xmv**2 / massr)) / (2 * xmv**2 / massr)) / xmv**3
                        nuDpv = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xpv / sqrt(massr)) &
                               - (erf(xpv / sqrt(massr)) - xpv / sqrt(massr) * (2 / sqrt(pi)) * exp(-xpv**2 / massr)) / (2 * xpv**2 / massr)) / xpv**3
                     end if

                     if (imu == 1) then
                        ! one-sided difference for mu-derivative at imu=1:
                        if (.not. conservative_test_particle_density) then
                           aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                + code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu+1,iz,is) / (2*dvpa) / dmu(imu)
                           aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv / dvpa**2 / mw(iv - 1, imu, iz, is) &
                                                      - code_dt * vpa(iv - 1) * mu(imu) * nux(iv - 1, imu, iz, is, isb) / (2 * dvpa) / dmu(imu)
                           cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                - code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu+1,iz,is) / (2*dvpa) / dmu(imu)
                           cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv / dvpa**2 / mw(iv + 1, imu, iz, is) &
                                                      + code_dt * vpa(iv + 1) * mu(imu) * nux(iv + 1, imu, iz, is, isb) / (2 * dvpa) / dmu(imu)
                        else if (conservative_test_particle_density) then
                           ! to ensure density conservation, assume that nux*F0 vanishes at iv=1 and iv=nvpa
                           aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
           + 0.5*code_dt*vpa(iv-1)*(mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)+mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is))&
                                                                    / mw(iv - 1, imu + 1, iz, is) / (2 * dvpa) / dmu(imu)
                           aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv / dvpa**2 / mw(iv - 1, imu, iz, is) &
           - 0.5*code_dt*vpa(iv-1)*(mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)+mu(imu+1)*nux(iv-1,imu+1,iz,is,isb)*mw(iv-1,imu+1,iz,is))&
                                                                / mw(iv - 1, imu, iz, is) / (2 * dvpa) / dmu(imu)
                           cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
           - 0.5*code_dt*vpa(iv+1)*(mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)+mu(imu+1)*nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))&
                                                                    / mw(iv + 1, imu + 1, iz, is) / (2 * dvpa) / dmu(imu)
                           cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv / dvpa**2 / mw(iv + 1, imu, iz, is) &
           + 0.5*code_dt*vpa(iv+1)*(mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)+mu(imu+1)*nux(iv+1,imu+1,iz,is,isb)*mw(iv+1,imu+1,iz,is))&
                                                                / mw(iv + 1, imu, iz, is) / (2 * dvpa) / dmu(imu)
                        end if
                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                     + code_dt * (0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv &
                 + 0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv) / dvpa**2 / mw(iv, imu, iz, is)

                        ! mu operator
                        if (include_mu_derivatives) then
                           ! use ghost cell at mu_{0} = 0, where term vanishes, so dmu(0) = mu(1).
                           mup = 0.5 * (mu(imu) + mu(imu + 1))
                           mwp = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mup)
                           xp = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mup)
                           if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
                              nuDp = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xp**3
                           else
                              nuDp = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xp / sqrt(massr)) &
                                    - (erf(xp / sqrt(massr)) - xp / sqrt(massr) * (2 / sqrt(pi)) * exp(-xp**2 / massr)) / (2 * xp**2 / massr)) / xp**3
                           end if
      nupap = spec(is)%vnew(isb) * 2 * (erf(xp / sqrt(massr)) - xp / sqrt(massr) * (2 / sqrt(pi)) * exp(-xp**2 / massr)) / (2 * xp**2 / massr) / xp**3
      gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*scale_pitch_angle_scattering*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                           gam_mup = 2 * (eiediff * nupap * mup**2 + eidefl * scale_pitch_angle_scattering * nuDp * vpa(iv)**2 / (2 * bmag(ia, iz)) * mup) * mwp
                   bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) - code_dt * (gam_mu*(-1) / dmu(imu) * dmu(imu) / 2./mu(imu) &
                         + (gam_mup*(-1) / dmu(imu) - gam_mu*(-1) / dmu(imu)) * mu(imu) / (dmu(imu) / 2.)) / mw(iv, imu, iz, is) / (dmu(imu) / 2.+mu(imu))
          bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) - code_dt * (gam_mu * 1 / dmu(imu) * dmu(imu) / 2./mu(imu) &
                   + (gam_mup * 1 / dmu(imu) - gam_mu * 1 / dmu(imu)) * mu(imu) / (dmu(imu) / 2.)) / mw(iv, imu + 1, iz, is) / (dmu(imu) / 2.+mu(imu))
                           ! mixed derivative:
                           if (conservative_test_particle_density) then
                              ! to ensure density conservation, we change discretisation of mixed derivative term at imu=1
                              ! from explicit routine: Dmuh1 = (vpa(iv)*mu(imu)*nux(iv,imu,iz)*mw(iv,imu,iz,is)*Dvpah &
                              aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
         + code_dt * (vpa(iv) * mu(imu) * nux(iv, imu, iz, is, isb) * mw(iv, imu, iz, is) * 1 / mw(iv - 1, imu, iz, is)) / (2 * dmu(imu)) / (2 * dvpa)
                              aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                             + code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv-1,imu+1,iz,is)) / (2*dmu(imu)) / (2*dvpa)
                              cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
         - code_dt * (vpa(iv) * mu(imu) * nux(iv, imu, iz, is, isb) * mw(iv, imu, iz, is) * 1 / mw(iv + 1, imu, iz, is)) / (2 * dmu(imu)) / (2 * dvpa)
                              cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                             - code_dt*(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,isb)*mw(iv,imu+1,iz,is)*1/mw(iv+1,imu+1,iz,is)) / (2*dmu(imu)) / (2*dvpa)
                           else
                              aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                        + code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv-1,imu,iz,is)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu)))&
                                                                   / (mu(imu) + dmu(imu)) / (2 * dvpa)
                              aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
+ code_dt * (vpa(iv) * mu(imu + 1) * nux(iv, imu + 1, iz, is, isb) * mw(iv, imu + 1, iz, is) * 1 / mw(iv - 1, imu + 1, iz, is) * mu(imu) / dmu(imu)) &
                                                                       / (mu(imu) + dmu(imu)) / (2 * dvpa)
                              cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                       - code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv+1,imu,iz,is)*(dmu(imu)/mu(imu)-mu(imu)/dmu(imu))) &
                                                                   / (mu(imu) + dmu(imu)) / (2 * dvpa)
                              cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
- code_dt * (vpa(iv) * mu(imu + 1) * nux(iv, imu + 1, iz, is, isb) * mw(iv, imu + 1, iz, is) * 1 / mw(iv + 1, imu + 1, iz, is) * mu(imu) / dmu(imu)) &
                                                                       / (mu(imu) + dmu(imu)) / (2 * dvpa)
                           end if
                        end if

                     else if (imu == nmu) then
                        ! one-sided difference for mu-derivative at imu=nmu:
                        ! to be consistent with treatment of mixed mu operator; assume that nux(imu)=0.
                        if (conservative_test_particle_density) then
                           aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
       - 1.0*code_dt*vpa(iv-1)*0.5*(mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is))&
                                                                    / mw(iv - 1, imu - 1, iz, is) / (2 * dvpa) / dmu(nmu - 1)
                           aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv / dvpa**2 / mw(iv - 1, imu, iz, is) &
       + 1.0*code_dt*vpa(iv-1)*0.5*(mu(imu-1)*nux(iv-1,imu-1,iz,is,isb)*mw(iv-1,imu-1,iz,is)+mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is))&
                                                                / mw(iv - 1, imu, iz, is) / (2 * dvpa) / dmu(nmu - 1)
                           cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
       + 1.0*code_dt*vpa(iv+1)*0.5*(mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is)+mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is))&
                                                                    / mw(iv + 1, imu - 1, iz, is) / (2 * dvpa) / dmu(nmu - 1)
                           cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv / dvpa**2 / mw(iv + 1, imu, iz, is) &
       - 1.0*code_dt*vpa(iv+1)*0.5*(mu(imu-1)*nux(iv+1,imu-1,iz,is,isb)*mw(iv+1,imu-1,iz,is)+mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is))&
                                                                / mw(iv + 1, imu, iz, is) / (2 * dvpa) / dmu(nmu - 1)
                        else
                           aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                          - 1.0*code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                           aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv / dvpa**2 / mw(iv - 1, imu, iz, is) &
                                            + 1.0 * code_dt * vpa(iv - 1) * mu(imu) * nux(iv - 1, imu, iz, is, isb) / (2 * dvpa) / dmu(nmu - 1)
                           cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                          + 1.0*code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu-1,iz,is) / (2*dvpa) / dmu(nmu-1)
                           cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv / dvpa**2 / mw(iv + 1, imu, iz, is) &
                                            - 1.0 * code_dt * vpa(iv + 1) * mu(imu) * nux(iv + 1, imu, iz, is, isb) / (2 * dvpa) / dmu(nmu - 1)
                        end if

                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                     + code_dt * (0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv &
                 + 0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv) / dvpa**2 / mw(iv, imu, iz, is)

                        ! mu operator
                        if (include_mu_derivatives) then
                           mum = 0.5 * (mu(imu) + mu(imu - 1))
                           mwm = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mum)
                           xm = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mum)
                           if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
                              nuDm = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xm**3
                           else
                              nuDm = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xm / sqrt(massr)) &
                                    - (erf(xm / sqrt(massr)) - xm / sqrt(massr) * (2 / sqrt(pi)) * exp(-xm**2 / massr)) / (2 * xm**2 / massr)) / xm**3
                           end if
      nupam = spec(is)%vnew(isb) * 2 * (erf(xm / sqrt(massr)) - xm / sqrt(massr) * (2 / sqrt(pi)) * exp(-xm**2 / massr)) / (2 * xm**2 / massr) / xm**3
      gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*scale_pitch_angle_scattering*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                           gam_mum = 2 * (eiediff * nupam * mum**2 + eidefl * scale_pitch_angle_scattering * nuDm * vpa(iv)**2 / (2 * bmag(ia, iz)) * mum) * mwm

                           if (conservative_test_particle_density) then
                              ! to ensure density conservation assume that argument of outer derivative vanishes at nmu+1/2, ie
                              ! d/dmu[...]_{nmu} = ([...]_{nmu+1/2} - [...]_{nmu-1/2})/(dmu_{nmu-1}/2+dmu(nmu-1)/2)
                              ! where [...]_{nmu+1/2} is a ghost cell at mu = mu_{nmu} + dmu(nmu-1)/2, with [...]_{nmu+1/2} = 0.
                              bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                                       - code_dt * ((gam_mu / dmu(imu - 1) - gam_mum / dmu(imu - 1)) * dmu(imu - 1) / (dmu(imu - 1)) &
                               + (-gam_mu / dmu(imu - 1)) * dmu(imu - 1) / dmu(imu - 1)) / mw(iv, imu, iz, is) / (dmu(imu - 1) / 2.+dmu(imu - 1) / 2.)
                              bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                                   - code_dt * ((-gam_mu / dmu(imu - 1) - gam_mum*(-1) / dmu(imu - 1)) * dmu(imu - 1) / (dmu(imu - 1)) &
                           - gam_mu*(-1) / dmu(imu - 1) * dmu(imu - 1) / dmu(imu - 1)) / mw(iv, imu - 1, iz, is) / (dmu(imu - 1) / 2.+dmu(imu - 1) / 2.)
                           else
                              bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                                  - code_dt * ((gam_mu / dmu(imu - 1) - gam_mum / dmu(imu - 1)) * dmu(imu - 1) / (dmu(imu - 1) / 2.) &
                                 + (-gam_mu / dmu(imu - 1)) * dmu(imu - 1) / 2./dmu(imu - 1)) / mw(iv, imu, iz, is) / (dmu(imu - 1) / 2.+dmu(imu - 1))
                              bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                              - code_dt * ((-gam_mu / dmu(imu - 1) - gam_mum*(-1) / dmu(imu - 1)) * dmu(imu - 1) / (dmu(imu - 1) / 2.) &
                             - gam_mu*(-1) / dmu(imu - 1) * dmu(imu - 1) / 2./dmu(imu - 1)) / mw(iv, imu - 1, iz, is) / (dmu(imu - 1) / 2.+dmu(imu - 1))
                           end if
                           ! no distinction here between density_conserving and default
                           if (conservative_test_particle_density) then
                              aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                           - code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv-1,imu-1,iz,is)) / (2*dmu(imu-1)) / (2*dvpa)
                              aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
     - code_dt * (vpa(iv) * mu(imu) * nux(iv, imu, iz, is, isb) * mw(iv, imu, iz, is) * 1 / mw(iv - 1, imu, iz, is)) / (2 * dmu(imu - 1)) / (2 * dvpa)
                              cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                           + code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv+1,imu-1,iz,is)) / (2*dmu(imu-1)) / (2*dvpa)
                              cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
     + code_dt * (vpa(iv) * mu(imu) * nux(iv, imu, iz, is, isb) * mw(iv, imu, iz, is) * 1 / mw(iv + 1, imu, iz, is)) / (2 * dmu(imu - 1)) / (2 * dvpa)
                           else
                              aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
           + code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv-1,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) &
                                                                   / (dmu(imu - 1) + dmu(imu - 1)) / (2 * dvpa)
                              aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                              - code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv-1,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1)) &
                                                                       / (dmu(imu - 1) + dmu(imu - 1)) / (2 * dvpa)
                              cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
           - code_dt*(vpa(iv)*mu(imu  )*nux(iv,imu  ,iz,is,isb)*mw(iv,imu,iz,is)*1/mw(iv+1,imu,iz,is)*(dmu(imu-1)/dmu(imu-1)-dmu(imu-1)/dmu(imu-1))) &
                                                                   / (dmu(imu - 1) + dmu(imu - 1)) / (2 * dvpa)
                              cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                              + code_dt*(vpa(iv)*mu(imu-1)*nux(iv,imu-1,iz,is,isb)*mw(iv,imu-1,iz,is)*1/mw(iv+1,imu-1,iz,is)* dmu(imu-1)/dmu(imu-1)) &
                                                                       / (dmu(imu - 1) + dmu(imu - 1)) / (2 * dvpa)
                           end if

                        end if
                     else ! interior mu points
                        if (conservative_test_particle_density) then
                           aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv / dvpa**2 / mw(iv - 1, imu, iz, is) &
              +code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb) * (dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                           cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv / dvpa**2 / mw(iv + 1, imu, iz, is) &
              -code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb) * (dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                           ! vpa operator, mixed (interior treatment):
                           aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                    - code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu-1,iz,is) * dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                           aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                    + code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu+1,iz,is) * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                           cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                    + code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu-1,iz,is) * dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                           cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                    - code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu+1,iz,is) * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                        else
                           aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupamv * vpam**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv / dvpa**2 / mw(iv - 1, imu, iz, is) &
              +code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb) * (dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                           cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
  - code_dt * 0.5 * (eiediff * nupapv * vpap**2 + eidefl * scale_pitch_angle_scattering * 2 * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv / dvpa**2 / mw(iv + 1, imu, iz, is) &
              -code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb) * (dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                           ! vpa operator, mixed (interior treatment):
                           aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                    - code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu-1,iz,is) * dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                           aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                    + code_dt*vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,isb)*mw(iv-1,imu,iz,is)/mw(iv-1,imu+1,iz,is) * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                           cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                                    + code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu-1,iz,is) * dmu(imu)/dmu(imu-1) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                           cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                    - code_dt*vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,isb)*mw(iv+1,imu,iz,is)/mw(iv+1,imu+1,iz,is) * dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu)) / (2*dvpa)
                        end if

                        bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                     + code_dt * (0.5 * (eiediff * nupapv * vpap**2 + 2 * eidefl * scale_pitch_angle_scattering * nuDpv * bmag(ia, iz) * mu(imu)) * mwpv &
                 + 0.5 * (eiediff * nupamv * vpam**2 + 2 * eidefl * scale_pitch_angle_scattering * nuDmv * bmag(ia, iz) * mu(imu)) * mwmv) / dvpa**2 / mw(iv, imu, iz, is)

                        ! mu operator
                        if (include_mu_derivatives) then
                           ! quantities at mu_{i+1/2} and mu_{i-1/2}:
                           mup = 0.5 * (mu(imu) + mu(imu + 1))
                           mum = 0.5 * (mu(imu) + mu(imu - 1))
                           mwp = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mup)
                           mwm = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mum)
                           xp = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mup)
                           xm = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mum)
                           if ((is == 2) .and. (isb == 1) .and. (mass_ratio_approximation_ei)) then
                              nuDp = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xp**3
                              nuDm = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) / xm**3
                           else
                              nuDp = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xp / sqrt(massr)) &
                                    - (erf(xp / sqrt(massr)) - xp / sqrt(massr) * (2 / sqrt(pi)) * exp(-xp**2 / massr)) / (2 * xp**2 / massr)) / xp**3
                              nuDm = eidefl * scale_pitch_angle_scattering * spec(is)%vnew(isb) * (erf(xm / sqrt(massr)) &
                                    - (erf(xm / sqrt(massr)) - xm / sqrt(massr) * (2 / sqrt(pi)) * exp(-xm**2 / massr)) / (2 * xm**2 / massr)) / xm**3
                           end if
      nupap = spec(is)%vnew(isb) * 2 * (erf(xp / sqrt(massr)) - xp / sqrt(massr) * (2 / sqrt(pi)) * exp(-xp**2 / massr)) / (2 * xp**2 / massr) / xp**3
      nupam = spec(is)%vnew(isb) * 2 * (erf(xm / sqrt(massr)) - xm / sqrt(massr) * (2 / sqrt(pi)) * exp(-xm**2 / massr)) / (2 * xm**2 / massr) / xm**3
      gam_mu = 2*(eiediff*nupa(iv,imu,iz,is,isb)*mu(imu)**2+eidefl*scale_pitch_angle_scattering*nuD(iv,imu,iz,is,isb)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)
                           gam_mum = 2 * (eiediff * nupam * mum**2 + eidefl * scale_pitch_angle_scattering * nuDm * vpa(iv)**2 / (2 * bmag(ia, iz)) * mum) * mwm
                           gam_mup = 2 * (eiediff * nupap * mup**2 + eidefl * scale_pitch_angle_scattering * nuDp * vpa(iv)**2 / (2 * bmag(ia, iz)) * mup) * mwp
                           ! include_mu_derivatives (interior treatment):
                           bb_blcs(iv, imu, imu, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                  - code_dt*((gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mum / dmu(imu-1))*dmu(imu)/dmu(imu-1) &
                                    +(-gam_mu*(dmu(imu)/dmu(imu-1) - dmu(imu-1)/dmu(imu)) / (dmu(imu-1)+dmu(imu)) - gam_mup / dmu(imu))* dmu(imu-1)/dmu(imu)) / mw(iv,imu,iz,is) * 2./(dmu(imu-1)+dmu(imu))
                           bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                - code_dt * ((gam_mu*(-1) * dmu(imu) / dmu(imu - 1) / (dmu(imu - 1) + dmu(imu)) - gam_mum*(-1) / dmu(imu - 1)) * dmu(imu) / dmu(imu - 1) &
 - gam_mu*(-1) * dmu(imu) / dmu(imu - 1) / (dmu(imu - 1) + dmu(imu)) * dmu(imu - 1) / dmu(imu)) / mw(iv, imu - 1, iz, is) * 2./(dmu(imu - 1) + dmu(imu))
                           bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = bb_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                                                 - code_dt * (gam_mu * dmu(imu - 1) / dmu(imu) / (dmu(imu - 1) + dmu(imu)) * dmu(imu) / dmu(imu - 1) &
       + (gam_mup/dmu(imu) - gam_mu*dmu(imu-1)/dmu(imu) / (dmu(imu-1)+dmu(imu))) * dmu(imu-1)/dmu(imu)) / mw(iv,imu+1,iz,is) * 2/(dmu(imu-1)+dmu(imu))
                           ! mu operator, mixed (interior treatment):
                           aa_blcs(iv, imu, imu, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                      + code_dt * (vpa(iv) * mu(imu) * nux(iv, imu, iz, is, isb) * mw(iv, imu, iz, is) * 1 / mw(iv - 1, imu, iz, is) &
                                                       * (dmu(imu) / dmu(imu - 1) - dmu(imu - 1) / dmu(imu))) / (dmu(imu - 1) + dmu(imu)) / (2 * dvpa)
                           aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                      - code_dt * (vpa(iv) * mu(imu - 1) * nux(iv, imu - 1, iz, is, isb) * mw(iv, imu - 1, iz, is) * 1 / mw(iv - 1, imu - 1, iz, is) &
                                                                                 * dmu(imu) / dmu(imu - 1)) / (dmu(imu - 1) + dmu(imu)) / (2 * dvpa)
                           aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = aa_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                      + code_dt * (vpa(iv) * mu(imu + 1) * nux(iv, imu + 1, iz, is, isb) * mw(iv, imu + 1, iz, is) * 1 / mw(iv - 1, imu + 1, iz, is) &
                                                                                 * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu)) / (2 * dvpa)
                           cc_blcs(iv, imu, imu, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu, ikxkyzt, is, isb) &
                                      - code_dt * (vpa(iv) * mu(imu) * nux(iv, imu, iz, is, isb) * mw(iv, imu, iz, is) * 1 / mw(iv + 1, imu, iz, is) &
                                                       * (dmu(imu) / dmu(imu - 1) - dmu(imu - 1) / dmu(imu))) / (dmu(imu - 1) + dmu(imu)) / (2 * dvpa)
                           cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu - 1, ikxkyzt, is, isb) &
                      + code_dt * (vpa(iv) * mu(imu - 1) * nux(iv, imu - 1, iz, is, isb) * mw(iv, imu - 1, iz, is) * 1 / mw(iv + 1, imu - 1, iz, is) &
                                                                                 * dmu(imu) / dmu(imu - 1)) / (dmu(imu - 1) + dmu(imu)) / (2 * dvpa)
                           cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) = cc_blcs(iv, imu, imu + 1, ikxkyzt, is, isb) &
                      - code_dt * (vpa(iv) * mu(imu + 1) * nux(iv, imu + 1, iz, is, isb) * mw(iv, imu + 1, iz, is) * 1 / mw(iv + 1, imu + 1, iz, is) &
                                                                                 * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu)) / (2 * dvpa)
                        end if
                     end if
                  end if            
               end do
            end do
         end do
      end do
      end do
      !$omp end parallel do

      ! Without the test-particle part (include_test_particle_part = .false., to test the field-particle part alone) all blocks
      ! are zero: the identity of 1 - C^{ab} is added once below, after the blocks of all species are summed.
      ! (Before Nov 2026 the diagonal blocks were set to 1 for every species pair, which gave (nspec + 1) times
      ! the identity.) The classical (FLR) diffusion below is part of the test-particle part and left out too.
      if (.not. include_test_particle_part) then
         aa_blcs = 0.
         cc_blcs = 0.
         bb_blcs = 0.
      end if

      ! construct full block matrix
      if ((exact_conservation_test_particle_matrix) .or. (field_particle_density_from_test_particle)) then
         ! this is memory intensive, operating with blockmatrix is slow
         ! currently used for exact conservation scheme on non-uniform grids
         ! AVB: to do - replace this with band matrix operations

         do is = 1, nspec
            do isb = 1, nspec
               do ikxkyzt = ikxkyzt_first, ikxkyzt_last
                  do iv = 1, nvpa
                     ! diagonal blocks:
                     blockmatrix(nmu * (iv - 1) + 1:nmu * iv, nmu * (iv - 1) + 1:nmu * iv, ikxkyzt, is, isb) = bb_blcs(iv, :, :, ikxkyzt, is, isb)
                     if (iv < nvpa) then
                        ! subdiagonal blocks:
                        blockmatrix(nmu * iv + 1:nmu * (iv + 1), nmu * (iv - 1) + 1:nmu * iv, ikxkyzt, is, isb) = aa_blcs(iv + 1, :, :, ikxkyzt, is, isb)
                        ! superdiagonal blocks:
                        blockmatrix(nmu * (iv - 1) + 1:nmu * iv, nmu * iv + 1:nmu * (iv + 1), ikxkyzt, is, isb) = cc_blcs(iv, :, :, ikxkyzt, is, isb)
                     end if
                  end do
               end do
            end do
         end do

      end if

      ! switch to band-storage for the band LU routines of utils/banded_lu (see band_row):
      ! and sum the interspecies and intraspecies operators for each species
      ! a_ij is stored in aband(kl+ku+1+i-j,j) for $\max(1,j-ku) \leq i \leq \min(m,j+kl)$
      !
      ! All inputs (aa/bb/cc_blcs) and the output (cdiffmat_band) live on
      ! kxkyzt_lo with all species local on each rank, so the assembly is
      ! fully local: no global temporary, no sum_allreduce, no copy step.

      !$omp parallel do default(none) &
      !$omp   firstprivate(kxkyzt_lo, ikxkyzt_first, ikxkyzt_last, nmu, nvpa, nspec) &
      !$omp   private(ikxkyzt, iky, ikx, iz, is, iv, imu, imu2, isb, bm_colind, bm_rowind) &
      !$omp   shared(cdiffmat_band, aa_blcs, bb_blcs, cc_blcs)
      do ikxkyzt = ikxkyzt_first, ikxkyzt_last
         iky = iky_idx(kxkyzt_lo, ikxkyzt)
         ikx = ikx_idx(kxkyzt_lo, ikxkyzt)
         iz = iz_idx(kxkyzt_lo, ikxkyzt)

         do is = 1, nspec  ! target species

            ! loop through aa_blcs, bb_blcs, cc_blcs
            ! find corresponding index of blockmatrix at each location within blcs
            ! calculate index of cdiffmat_band

            do iv = 1, nvpa

               ! bb_blcs (intra-species: target=source=is)
               do imu = 1, nmu
                  bm_rowind = (iv - 1) * nmu + imu
                  do imu2 = 1, nmu
                     bm_colind = (iv - 1) * nmu + imu2
                     if ((max(1, bm_colind - (nmu + 1)) <= bm_rowind) .and. (bm_rowind <= min(nvpa * nmu, bm_colind + (nmu + 1)))) then
                        cdiffmat_band(band_row(bm_rowind, bm_colind, nmu + 1, nmu + 1), bm_colind, ikxkyzt, is) = bb_blcs(iv, imu, imu2, ikxkyzt, is, is)
                     end if
                  end do
               end do
               ! cc_blcs (intra-species)
               if (iv < nvpa) then ! aa_blcs and cc_blcs contain only (nvpa-1) blocks, since they are off-diagonal
                  do imu = 1, nmu
                     bm_rowind = (iv - 1) * nmu + imu
                     do imu2 = 1, nmu
                        bm_colind = nmu + (iv - 1) * nmu + imu2 ! nvpa*nmu - nmu + imu2
                        if ((max(1, bm_colind - (nmu + 1)) <= bm_rowind) .and. (bm_rowind <= min(nvpa * nmu, bm_colind + (nmu + 1)))) then
                           cdiffmat_band(band_row(bm_rowind, bm_colind, nmu + 1, nmu + 1), bm_colind, ikxkyzt, is) = cc_blcs(iv, imu, imu2, ikxkyzt, is, is)
                        end if
                     end do
                  end do
                  ! aa_blcs (intra-species)
                  do imu = 1, nmu
                     bm_rowind = nmu + (iv - 1) * nmu + imu
                     do imu2 = 1, nmu
                        bm_colind = (iv - 1) * nmu + imu2
                        if ((max(1, bm_colind - (nmu + 1)) <= bm_rowind) .and. (bm_rowind <= min(nvpa * nmu, bm_colind + (nmu + 1)))) then
                           cdiffmat_band(band_row(bm_rowind, bm_colind, nmu + 1, nmu + 1), bm_colind, ikxkyzt, is) = aa_blcs(1 + iv, imu, imu2, ikxkyzt, is, is)
                        end if
                     end do
                  end do
               end if

               ! inter-species test particle contributions:
               do isb = 1, nspec
                  if (isb == is) cycle
                  ! bb_blcs
                  do imu = 1, nmu
                     bm_rowind = (iv - 1) * nmu + imu
                     do imu2 = 1, nmu
                        bm_colind = (iv - 1) * nmu + imu2
                        if ((max(1, bm_colind - (nmu + 1)) <= bm_rowind) .and. (bm_rowind <= min(nvpa * nmu, bm_colind + (nmu + 1)))) then
                           cdiffmat_band(band_row(bm_rowind, bm_colind, nmu + 1, nmu + 1), bm_colind, ikxkyzt, is) = &
                          cdiffmat_band(band_row(bm_rowind, bm_colind, nmu + 1, nmu + 1), bm_colind, ikxkyzt, is) + bb_blcs(iv, imu, imu2, ikxkyzt, is, isb)
                        end if
                     end do
                  end do
                  ! cc_blcs
                  if (iv < nvpa) then
                     do imu = 1, nmu
                        bm_rowind = (iv - 1) * nmu + imu
                        do imu2 = 1, nmu
                           bm_colind = (iv - 1) * nmu + nmu + imu2
                           if ((max(1, bm_colind - (nmu + 1)) <= bm_rowind) .and. (bm_rowind <= min(nvpa * nmu, bm_colind + (nmu + 1)))) then
                              cdiffmat_band(band_row(bm_rowind, bm_colind, nmu + 1, nmu + 1), bm_colind, ikxkyzt, is) = &
                             cdiffmat_band(band_row(bm_rowind, bm_colind, nmu + 1, nmu + 1), bm_colind, ikxkyzt, is) + cc_blcs(iv, imu, imu2, ikxkyzt, is, isb)
                           end if
                        end do
                     end do
                     ! aa_blcs
                     do imu = 1, nmu
                        bm_rowind = (iv - 1) * nmu + nmu + imu
                        do imu2 = 1, nmu
                           bm_colind = (iv - 1) * nmu + imu2
                           if ((max(1, bm_colind - (nmu + 1)) <= bm_rowind) .and. (bm_rowind <= min(nvpa * nmu, bm_colind + (nmu + 1)))) then
                              cdiffmat_band(band_row(bm_rowind, bm_colind, nmu + 1, nmu + 1), bm_colind, ikxkyzt, is) = &
                             cdiffmat_band(band_row(bm_rowind, bm_colind, nmu + 1, nmu + 1), bm_colind, ikxkyzt, is) + aa_blcs(1 + iv, imu, imu2, ikxkyzt, is, isb)
                           end if
                        end do
                     end do
                  end if
               end do

            end do  ! iv

         end do  ! is (target species)

      end do  ! ikxkyzt
      !$omp end parallel do

      if (exact_conservation_block_products) call get_tp_f0moments(ikxkyzt_first, ikxkyzt_last)

      deallocate (aa_blcs, bb_blcs, cc_blcs)
      end do  ! chunks of ikxkyzt

      ! add the gyro-diffusive term to cdiffmat (diagonal)
      !$omp parallel do default(none) &
      !$omp   firstprivate(kxkyzt_lo, nmu, nvpa, nspec, ia) &
      !$omp   private(ikxkyzt, iky, ikx, iz, is, iv, imu, ivv) &
      !$omp   shared(cdiffmat_band, kperp2, code_dt, scale_flr_terms, include_test_particle_part)
      do ikxkyzt = kxkyzt_lo%llim_proc, kxkyzt_lo%ulim_proc
         iky = iky_idx(kxkyzt_lo, ikxkyzt)
         ikx = ikx_idx(kxkyzt_lo, ikxkyzt)
         iz = iz_idx(kxkyzt_lo, ikxkyzt)
         do is = 1, nspec
            do iv = 1, nvpa
               do imu = 1, nmu
                  ivv = nmu * (iv - 1) + imu
                  if (include_test_particle_part) cdiffmat_band(band_row(ivv, ivv, nmu + 1, nmu + 1), ivv, ikxkyzt, is) = cdiffmat_band(band_row(ivv, ivv, nmu + 1, nmu + 1), ivv, ikxkyzt, is) &
                     + code_dt * scale_flr_terms * kperp2(iky, ikx, ia, iz) * gyro_diffusion_coefficient(iv, imu, iz, is)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! add 1 to the diagonal, since the matrix operator is 1 - C^{ab}[h_a].
      ! All ranks operate on their own kxkyzt slices independently -- no
      ! sum_allreduce or temporary global is needed because there is no
      ! overlapping ownership across ranks.
      do is = 1, nspec
         do iv = 1, nmu * nvpa
            do imu = 1, nmu * nvpa
               if ((max(1, iv - (nmu + 1)) <= imu) .and. (imu <= min(nvpa * nmu, iv + (nmu + 1)))) then
                  if (iv == imu) then
                     cdiffmat_band(band_row(imu, iv, nmu + 1, nmu + 1), iv, :, is) = cdiffmat_band(band_row(imu, iv, nmu + 1, nmu + 1), iv, :, is) + 1.
                  end if
               end if
            end do
         end do
      end do

      ! AVB: LU factorise cdiffmat, using the band LU routines of utils/banded_lu (LAPACK's zgbtrf, or PLASMA,
      ! which parallelises each factorisation itself, so then the loop is serial; inside the parallel loop
      ! every factorisation uses zgbtrf). Each rank factorises only the (ikxkyzt, is) slices it owns.
      nc = nvpa * nmu
      nb = nmu + 1
      factorise_threaded = band_lu_factorise_threaded()
      !$omp parallel do default(none) collapse(2) if(.not. factorise_threaded) &
      !$omp   firstprivate(kxkyzt_lo, nc, nb, nspec) &
      !$omp   private(ikxkyzt, is, info) &
      !$omp   shared(cdiffmat_band, ipiv)
      do is = 1, nspec
         do ikxkyzt = kxkyzt_lo%llim_proc, kxkyzt_lo%ulim_proc
            call band_lu_factorise(nc, nb, nb, cdiffmat_band(:, :, ikxkyzt, is), ipiv(:, ikxkyzt, is), info)
         end do
      end do
      !$omp end parallel do

   end subroutine init_fp_diffmatrix

   !****************************************************************************
   !                  The z points of the coefficients on this processor
   !****************************************************************************
   ! The coefficients of the field-particle operator (deltaj, deltaj_tp), the Legendre
   ! functions (legendre_vpamu) and the Bessel functions (jm) are only used at the
   ! z points of the kxkyz_lo and kxkyzt_lo points of this processor, so every processor
   ! stores them only there, which needs no communication. These z points are numbered
   ! izl(iz) = 1, ..., nz_local (0 if iz is not stored here), with iz = iz_of_izl(izl).
   ! They also include the share iz_share_lo to iz_share_hi of the z points on which this
   ! processor computes psijnorm, which is needed on every z point (see init_deltaj_vmu).
   !****************************************************************************
   subroutine init_local_z

      use mp, only: nproc, iproc
      use grids_z, only: nzgrid
      use parallelisation_layouts, only: kxkyz_lo, kxkyzt_lo, iz_idx

      implicit none

      integer :: i, iz, nz_total

      !-------------------------------------------------------------------------

      ! The map does not change when the time step changes (reset_dt_collisions_fp)
      if (allocated(izl)) return

      allocate (local_z(-nzgrid:nzgrid))
      local_z = .false.
      do i = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         local_z(iz_idx(kxkyz_lo, i)) = .true.
      end do
      do i = kxkyzt_lo%llim_proc, kxkyzt_lo%ulim_proc
         local_z(iz_idx(kxkyzt_lo, i)) = .true.
      end do
      nz_total = 2 * nzgrid + 1
      iz_share_lo = -nzgrid + (iproc * nz_total) / nproc
      iz_share_hi = -nzgrid + ((iproc + 1) * nz_total) / nproc - 1
      local_z(iz_share_lo:iz_share_hi) = .true.

      allocate (izl(-nzgrid:nzgrid), iz_of_izl(max(1, count(local_z))))
      izl = 0
      nz_local = 0
      do iz = -nzgrid, nzgrid
         if (local_z(iz)) then
            nz_local = nz_local + 1
            izl(iz) = nz_local
            iz_of_izl(nz_local) = iz
         end if
      end do

   end subroutine init_local_z

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_legendre

      use grids_velocity, only: mu, nmu, vpa, nvpa
      use geometry, only: bmag
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx

      implicit none

      integer :: iv, imu, iz, izs, ia, mm, ll
      double precision :: xi

      !-------------------------------------------------------------------------

      ! Only on the z points of this processor, see init_local_z
      allocate (legendre_vpamu(0:lmax, -lmax:lmax, nvpa, nmu, nz_local))

      ! note lmin = 0, lmax = nsph-1
      ! mmin = -lmax, mmax = lmax

      legendre_vpamu = 0.

      ia = 1

      !$omp parallel do default(none) &
      !$omp   firstprivate(nvpa, nmu, nz_local, lmax, ia) &
      !$omp   private(iv, imu, iz, izs, ll, mm, xi) &
      !$omp   shared(legendre_vpamu, vpa, mu, bmag, iz_of_izl)
      do iv = 1, nvpa
         do imu = 1, nmu
            do izs = 1, nz_local
               iz = iz_of_izl(izs)
               do ll = 0, lmax
                  do mm = -lmax, lmax
                     xi = dble(vpa(iv) / sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mu(imu)))
                     legendre_vpamu(ll, mm, iv, imu, izs) = associated_legendre(ll, mm, xi)
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine init_legendre

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_bessel_fn
      use grids_velocity, only: nmu, vperp2
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use spfunc, only: j0
      use grids_species, only: spec, nspec
      use geometry, only: bmag
      use grids_kxky, only: naky, nakx
      use arrays, only: kperp2
      use mp, only: sum_allreduce

      implicit none

      integer :: iky, ikx, iz, izs, is, ia, mm, imu
      real :: arg

      !-------------------------------------------------------------------------

      ! Only on the z points of this processor, see init_local_z
      allocate (jm(nmu, 0:lmax, naky, nakx, nz_local, nspec))

      jm = 0
      ia = 1

      ! The Bessel functions are used at the kxkyz_lo and kxkyzt_lo points of every species (fp_response,
      ! get_psi_response, get_psi and the conservation terms), so every processor computes them on all
      ! (ky, kx, species) of its z points, J0 with the same argument as aj0v in arrays_gyro_averages
      do izs = 1, nz_local
         iz = iz_of_izl(izs)
         do is = 1, nspec
            do ikx = 1, nakx
               do iky = 1, naky
                  do imu = 1, nmu
                     arg = spec(is)%bess_fac * spec(is)%smz_psi0 * sqrt(vperp2(ia, iz, imu) * kperp2(iky, ikx, ia, iz)) / bmag(ia, iz)
                     jm(imu, 0, iky, ikx, izs, is) = j0(arg)
                  end do
                  ! J_m of m >= 1 have the same argument as J_0, with bess_fac (0 for a drift-kinetic species)
                  do mm = 1, lmax
                     do imu = 1, nmu
                        arg = spec(is)%bess_fac * spec(is)%smz_psi0 * sqrt(vperp2(ia, iz, imu) * kperp2(iky, ikx, ia, iz)) / bmag(ia, iz)
                        jm(imu, mm, iky, ikx, izs, is) = bessel_jn(mm, arg)
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! The m /= 0 terms of the field-particle part (J1 for l = 1: the perpendicular momentum) restore what the
      ! classical (gyro-)diffusion of the test-particle part removes, so scale them by the same <scale_flr_terms>. With scale_flr_terms = 0
      ! (no classical diffusion) they would otherwise inject perpendicular momentum, and the growth rate runs away with
      ! the collision frequency (benchmark C1, Oct 2026: gamma = 0.28 instead of 0.08 at nu_ii = 1). J0 is kept on the
      ! m = 0 terms, so that the momentum and energy at the particle position are conserved, as in GS2, GENE and CGYRO.
      ! This replaces the input <cfac2>, which removed J0 and J1 of the field-particle part together.
      if (lmax > 0) then
         do mm = 1, lmax
            jm(:, mm, :, :, :, :) = scale_flr_terms * jm(:, mm, :, :, :, :)
         end do
      end if

   end subroutine init_bessel_fn

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_deltaj_vmu

      use grids_species, only: nspec, spec
      use grids_z, only: nzgrid
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_velocity, only: nmu, vpa, nvpa, set_vpa_weights
      use constants, only: pi
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, it_idx
      use mp, only: sum_allreduce

      implicit none

      real, dimension(0:lmax, 0:jmax, nvpa, nmu, 1, -nzgrid:nzgrid) :: vlaguerre_vmu
      real, dimension(:, :, :), allocatable :: deltaj_share
      integer :: ll, iv, jj, ia, iz, izs, isa, isb, ix
      real, dimension(-nzgrid:nzgrid) :: deltajint, deltajint_tp
      real, dimension(nvpa*nmu, -nzgrid:nzgrid) :: vpaF0vec, v2F0vec

      logical :: conservative_wgts

      !-------------------------------------------------------------------------

      ia = 1

      ! Only on the z points of this processor, see init_local_z
      allocate (deltaj(0:lmax, 0:jmax, nspec, nspec, nvpa, nmu, ia, nz_local))
      allocate (deltaj_tp(0:lmax, 0:jmax, nspec, nspec, nvpa, nmu, ia, nz_local))
      allocate (psijnorm(0:lmax, 0:jmax, nspec, nspec, -nzgrid:nzgrid))

      allocate (mwnorm(-nzgrid:nzgrid))
      allocate (modmwnorm(-nzgrid:nzgrid))

      if (conservative_test_particle_density) then
         conservative_wgts = .true.
         call set_vpa_weights(conservative_wgts)
      else if (exact_conservation_test_particle_matrix) then
         conservative_wgts = .false.
         call set_vpa_weights(conservative_wgts)
      else
         conservative_wgts = .false.
         call set_vpa_weights(conservative_wgts)
      end if

      ! AVB: to do - option for cons_wgts when exact_conservation is true

      ! get Delta_j^{l,ab}[x_b^l L_j^{l+0.5}(x_b^2)F_{0b}](x_a)
      ! and Delta_j^{l,ba}[x_a^l L_j^{l+0.5}(x_a^2)F_{0a}](x_b)

      deltaj = 0

      ! construct an identity matrix required below

      ! Delta_j is local in z, and its calculation (with the recursion of calc_deltaj_vmu and the incomplete gamma
      ! functions of calc_delta0) is expensive, so every processor computes it only on its own z points (local_z, see
      ! init_local_z), with the (vpa, mu, z) loop of calc_delta0 shared among the threads
      allocate (deltaj_share(nvpa, nmu, -nzgrid:nzgrid))
      do isa = 1, nspec
         do isb = 1, nspec
            do ll = 0, lmax
               do jj = 0, jmax
                  call calc_deltaj_vmu(jj, jj, ll, isa, isb, deltaj_share)
                  do izs = 1, nz_local
                     deltaj(ll, jj, isa, isb, :, :, ia, izs) = deltaj_share(:, :, iz_of_izl(izs))
                  end do
               end do
            end do
         end do
      end do
      deallocate (deltaj_share)

      do isa = 1, nspec
         do isb = 1, nspec
            do ll = 0, lmax
               do jj = 0, jmax

                  call vLj_vmu(jj, ll, vlaguerre_vmu(ll, jj, :, :, ia, :))

                  if (spitzer_problem) then
                     if ((exact_conservation_block_products) .and. (ll == 1) .and. (jj == 0) &
                         .and. .not. ((isa == 1) .and. (isb == 1)) &
                         .and. .not. ((isa == 1) .and. (isb == 2)) &
                         ) then

                        ! to ensure conservation of momentum to machine precision

                        ! calculate delta_{j,l}^{ab} using the differential test particle operator C^{ab} acting on m_a v_\parallel F0a, on every z point
                        call apply_testpart_matrix(vpaF0vec, 1, isa, isb, (spec(isa)%mass / spec(isb)%mass)**2)

                        do izs = 1, nz_local
                           do ix = 1, nvpa
                              deltaj_tp(ll, jj, isa, isb, ix, :, ia, izs) = vpaF0vec(nmu * (ix - 1) + 1:nmu * (ix - 1) + nmu, iz_of_izl(izs))
                           end do
                        end do

                        do izs = 1, nz_local
                           deltaj_tp(ll, jj, isa, isb, :, :, ia, izs) = deltaj_tp(ll, jj, isa, isb, :, :, ia, izs) * velvpamu(:, :, iz_of_izl(izs)) / spread(vpa, 2, nmu)
                        end do

                     end if

                  else
                     if ((exact_conservation_block_products) .and. (ll == 1) .and. (jj == 0)) then

                        if (spec(isa)%vnew(isb) == 0) cycle

                        ! calculate delta_{j,l}^{ab} using the differential test particle operator C^{ab} acting on m_a v_\parallel F0a, on every z point
                        call apply_testpart_matrix(vpaF0vec, 1, isa, isb, (spec(isa)%mass / spec(isb)%mass)**2)

                        do izs = 1, nz_local
                           do ix = 1, nvpa
                              deltaj_tp(ll, jj, isa, isb, ix, :, ia, izs) = vpaF0vec(nmu * (ix - 1) + 1:nmu * (ix - 1) + nmu, iz_of_izl(izs))
                           end do
                        end do

                        do izs = 1, nz_local
                           deltaj_tp(ll, jj, isa, isb, :, :, ia, izs) = deltaj_tp(ll, jj, isa, isb, :, :, ia, izs) * velvpamu(:, :, iz_of_izl(izs)) / spread(vpa, 2, nmu)
                        end do

                     end if
                  end if

                  if (spitzer_problem) then
                     if ((exact_conservation_block_products) .and. (ll == 0) .and. (jj == 1) &
                         .and. .not. ((isa == 1) .and. (isb == 1)) &
                         .and. .not. ((isa == 1) .and. (isb == 2)) &
                         ) then

                        ! energy conservation to machine precision

                        ! calculate delta_{j,l}^{ab} using the differential test particle operator C^{ab} acting on m_a v^2 F0a, on every z point
                        call apply_testpart_matrix(v2F0vec, 2, isa, isb, -(spec(isa)%mass / spec(isb)%mass)**1.5)

                        do izs = 1, nz_local
                           do ix = 1, nvpa
                              deltaj_tp(ll, jj, isa, isb, ix, :, ia, izs) = v2F0vec(nmu * (ix - 1) + 1:nmu * (ix - 1) + nmu, iz_of_izl(izs))
                           end do
                        end do

                     end if

                  else
                     if ((exact_conservation_block_products) .and. (ll == 0) .and. (jj == 1)) then

                        if (spec(isa)%vnew(isb) == 0) cycle

                        ! calculate delta_{j,l}^{ab} using the differential test particle operator C^{ab} acting on m_a v^2 F0a, on every z point
                        call apply_testpart_matrix(v2F0vec, 2, isa, isb, -(spec(isa)%mass / spec(isb)%mass)**1.5)

                        do izs = 1, nz_local
                           do ix = 1, nvpa
                              deltaj_tp(ll, jj, isa, isb, ix, :, ia, izs) = v2F0vec(nmu * (ix - 1) + 1:nmu * (ix - 1) + nmu, iz_of_izl(izs))
                           end do
                        end do

                     end if
                  end if

                  ! required integrals to ensure density conservation in field particle operator:
                  do iz = -nzgrid, nzgrid
                     call integrate_vmu(mw(:, :, iz, isa), iz, mwnorm(iz))
                  end do

                  do iz = -nzgrid, nzgrid
                     call integrate_vmu(modmw(:, :, iz, isa), iz, modmwnorm(iz))
                  end do

                  ! The l = 0 terms of the field-particle part (energy exchange, j >= 1) must not change the density.
                  ! For electrons colliding with ions, Delta_j^{ab} depends on the ion speed, which is sqrt(m_i/m_e) = 43
                  ! times the electron speed, so part of it lies inside the first cells of the electron grid. On the grid
                  ! its density moment was then as large as its energy moment (benchmark C2, Nov 2026), so the energy
                  ! exchange created electron density and the growth rate ran away at large collisionality. Remove the
                  ! density moment on the grid with modmw, which has no energy moment, so the energy exchange is unchanged.
                  if ((ll == 0) .and. (jj >= 1)) then
                     do izs = 1, nz_local
                        iz = iz_of_izl(izs)
                        call integrate_vmu(deltaj(ll, jj, isa, isb, :, :, ia, izs), iz, deltajint(iz))
                        deltaj(ll, jj, isa, isb, :, :, ia, izs) = deltaj(ll, jj, isa, isb, :, :, ia, izs) &
                                                - modmw(:, :, iz, isa) / modmwnorm(iz) * deltajint(iz)
                     end do
                  end if

                  ! Likewise, with exact conservation the l = 1 terms with j >= 1 must not change the parallel momentum, which
                  ! only the l = 1, j = 0 term restores. Analytically their Laguerre polynomials make them orthogonal to the
                  ! momentum, but not on the velocity grid: with 8 equally spaced mu points the l = 1, j = 1 term restored 4% more
                  ! momentum than the test-particle part removed (Nov 2026). Remove the momentum moment of Delta_j on the grid
                  ! with v F0, which times P_1 is vpa F0, with no density and (odd in vpa) no energy.
                  if ((exact_conservation_block_products .or. exact_conservation_test_particle_matrix) .and. (ll == 1) .and. (jj >= 1)) then
                     do izs = 1, nz_local
                        iz = iz_of_izl(izs)
                        call integrate_vmu(spread(vpa, 2, nmu) * legendre_vpamu(1, 0, :, :, izs) * deltaj(ll, jj, isa, isb, :, :, ia, izs), &
                                           iz, deltajint(iz))
                        call integrate_vmu(spread(vpa, 2, nmu)**2 * mw(:, :, iz, isa), iz, deltajint_tp(iz))
                        deltaj(ll, jj, isa, isb, :, :, ia, izs) = deltaj(ll, jj, isa, isb, :, :, ia, izs) &
                                                - velvpamu(:, :, iz) * mw(:, :, iz, isa) * deltajint(iz) / deltajint_tp(iz)
                     end do
                  end if

                  ! to ensure number conservation to machine precision, \Delta_{j=1} -> \Delta_{j=1} - F0/\int F0 dv * \int \Delta_{j=1} dv

                  ! AVB: to do - need to generalise this to higher-order terms in the field particle operator

                  if ((exact_field_particle_density) .and. (ll == 0) .and. (jj == 1)) then

                     do izs = 1, nz_local
                        iz = iz_of_izl(izs)
                        call integrate_vmu(deltaj(ll, jj, isa, isb, :, :, ia, izs), iz, deltajint(iz))
                        call integrate_vmu(deltaj_tp(ll, jj, isa, isb, :, :, ia, izs), iz, deltajint_tp(iz))
                     end do

                     ! remove the density moment of deltaj and deltaj_tp

                     do izs = 1, nz_local
                        iz = iz_of_izl(izs)
                        deltaj(ll, jj, isa, isb, :, :, ia, izs) = deltaj(ll, jj, isa, isb, :, :, ia, izs) &
                                                - mw(:, :, iz, isa) / mwnorm(iz) * deltajint(iz)
                        deltaj_tp(ll, jj, isa, isb, :, :, ia, izs) = deltaj_tp(ll, jj, isa, isb, :, :, ia, izs) &
                                             - mw(:, :, iz, isa) / mwnorm(iz) * deltajint_tp(iz)
                     end do

                  end if

               end do
            end do
         end do
      end do

      ! get normalisations for psijls; psijnorm is needed on every z point (get_psi), so every processor computes it
      ! on its share of the z points (iz_share_lo to iz_share_hi, see init_local_z) and the sum over the processors
      ! gives it everywhere
      psijnorm = 0.
      !$omp parallel do default(none) collapse(4) &
      !$omp   firstprivate(nspec, lmax, jmax, ia, iz_share_lo, iz_share_hi) &
      !$omp   private(isa, isb, ll, iz, jj) &
      !$omp   shared(psijnorm, legendre_vpamu, vlaguerre_vmu, deltaj, velvpamu, spec, izl, &
      !$omp          exact_conservation_block_products, exact_conservation_test_particle_matrix)
      do isa = 1, nspec
         do isb = 1, nspec
            do ll = 0, lmax
               do iz = iz_share_lo, iz_share_hi
                  do jj = 0, jmax
                     if ((ll == 0) .and. (jj == 0)) then
                        psijnorm(ll, jj, isa, isb, iz) = 1 ! never used
                     else

                        ! note self-adjointness of deltaj after normalisation and without coll freqs (\Delta_j') is
                        ! \int x_a^l L_j^{l+1/2}(x_a^2) \Delta_j'^{ab}[\tilde{f}_b/F0b * F0b] x_a^2 dx_a
                        ! = ma^0.5/mb^0.5 * ma^3/mb^3 \int f_b/F0b \Delta_j^{ba}'[x_a^l L_j^{l+1/2}(x_a^2) \tilde{F0a}] x_b^2 dx_b

                        if ((exact_conservation_block_products) .and. (ll == 1) .and. (jj == 0)) then
                           ! momentum conservation term for uniform mu grid
            call integrate_vmu(3 * legendre_vpamu(ll, 0, :, :, izl(iz))**2 * vlaguerre_vmu(ll, jj, :, :, ia, iz) * deltaj(ll, jj, isa, isb, :, :, ia, izl(iz)) &
                            * (spec(isa)%mass / spec(isb)%mass)**(-3) * (spec(isa)%mass / spec(isb)%mass)**(-0.5), iz, psijnorm(ll, jj, isa, isb, iz))
                        else if ((exact_conservation_block_products) .and. (ll == 0) .and. (jj == 1)) then
                           ! energy conservation term for uniform mu grid
                call integrate_vmu(legendre_vpamu(ll, 0, :, :, izl(iz))**2 * vlaguerre_vmu(ll, jj, :, :, ia, iz) * deltaj(ll, jj, isa, isb, :, :, ia, izl(iz)) &
                            * (spec(isa)%mass / spec(isb)%mass)**(-3) * (spec(isa)%mass / spec(isb)%mass)**(-0.5), iz, psijnorm(ll, jj, isa, isb, iz))
                        else if ((exact_conservation_test_particle_matrix) .and. (ll == 0) .and. (jj == 1)) then
                           ! energy conservation term for non-uniform mu grid
                call integrate_vmu(legendre_vpamu(ll, 0, :, :, izl(iz))**2 * vlaguerre_vmu(ll, jj, :, :, ia, iz) * deltaj(ll, jj, isa, isb, :, :, ia, izl(iz)) &
                            * (spec(isa)%mass / spec(isb)%mass)**(-3) * (spec(isa)%mass / spec(isb)%mass)**(-0.5), iz, psijnorm(ll, jj, isa, isb, iz))
                        else if ((exact_conservation_test_particle_matrix) .and. (ll == 1) .and. (jj == 0)) then
                           ! momentum conservation term for non-uniform mu grid
            call integrate_vmu(3 * legendre_vpamu(ll, 0, :, :, izl(iz))**2 * vlaguerre_vmu(ll, jj, :, :, ia, iz) * deltaj(ll, jj, isa, isb, :, :, ia, izl(iz)) &
                            * (spec(isa)%mass / spec(isb)%mass)**(-3) * (spec(isa)%mass / spec(isb)%mass)**(-0.5), iz, psijnorm(ll, jj, isa, isb, iz))
                        else if (.not. ((exact_conservation_block_products) .or. (exact_conservation_test_particle_matrix)) .and. (ll == 0) .and. (jj == 1)) then
                           ! non-exact conservation of energy
                           if (spec(isa)%mass / spec(isb)%mass < 1.) then
                              ! use self-adjointness to avoid resolution problems
                             call integrate_vmu((-velvpamu(:, :, iz)**2) * deltaj(ll, jj, isb, isa, :, :, ia, izl(iz)), iz, psijnorm(ll, jj, isa, isb, iz))
                           else
                      call integrate_vmu((-velvpamu(:, :, iz)**2) * deltaj(ll, jj, isa, isb, :, :, ia, izl(iz)) * (spec(isa)%mass / spec(isb)%mass)**(-3) &
                                                 * (spec(isa)%mass / spec(isb)%mass)**(-0.5), iz, psijnorm(ll, jj, isa, isb, iz))
                           end if
                        else
                           ! non-exact conservation of momentum, and higher-order terms
                           if (spec(isa)%mass / spec(isb)%mass < 1.) then
                              ! use self-adjointness to avoid resolution problems
                  call integrate_vmu(vlaguerre_vmu(ll, jj, :, :, ia, iz) * deltaj(ll, jj, isb, isa, :, :, ia, izl(iz)), iz, psijnorm(ll, jj, isa, isb, iz))
                           else
           call integrate_vmu(vlaguerre_vmu(ll, jj, :, :, ia, iz) * deltaj(ll, jj, isa, isb, :, :, ia, izl(iz)) * (spec(isa)%mass / spec(isb)%mass)**(-3) &
                                                 * (spec(isa)%mass / spec(isb)%mass)**(-0.5), iz, psijnorm(ll, jj, isa, isb, iz))
                           end if
                        end if

                     end if
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do
      call sum_allreduce(psijnorm)

      psijnorm = psijnorm / (4 * pi) ! account for theta and phi integrals included in integrate_vmu

      conservative_wgts = .false.
      call set_vpa_weights(conservative_wgts)

   end subroutine init_deltaj_vmu

   !****************************************************************************
   !        The test-particle operator applied to a moment of F0 on every z
   !****************************************************************************
   ! vec(:, iz) = factor * C^{isa, isb}[F0 moment] = factor * tp_f0moments / code_dt / vnew
   ! on every z point, with moment 1 for vpa F0a and 2 for v^2 F0a, for the exact
   ! conservation terms of init_deltaj_vmu, which do not depend on (ky, kx). The operator
   ! of one point per z point, (iky = naky, ikx = 1, tube 1), is used (get_tp_f0moments);
   ! it is on one processor (kxkyzt_lo), and the sum over the processors gives the result
   ! on every processor.
   !****************************************************************************
   subroutine apply_testpart_matrix(vec, moment, isa, isb, factor)

      use mp, only: sum_allreduce
      use grids_z, only: nzgrid
      use grids_kxky, only: naky
      use grids_species, only: spec
      use grids_time, only: code_dt
      use parallelisation_layouts, only: kxkyzt_lo, idx

      implicit none

      real, dimension(:, -nzgrid:), intent(out) :: vec
      integer, intent(in) :: moment, isa, isb
      real, intent(in) :: factor

      integer :: iz, ikxkyzt

      !-------------------------------------------------------------------------

      do iz = -nzgrid, nzgrid
         ikxkyzt = idx(kxkyzt_lo, naky, 1, iz, 1)
         ! without a-b collisions (vnew = 0) the operator and the term are zero
         if (ikxkyzt >= kxkyzt_lo%llim_proc .and. ikxkyzt <= kxkyzt_lo%ulim_proc .and. spec(isa)%vnew(isb) /= 0.) then
            vec(:, iz) = factor * tp_f0moments(:, moment, izl(iz), isa, isb) / code_dt / spec(isa)%vnew(isb)
         else
            vec(:, iz) = 0.
         end if
      end do
      call sum_allreduce(vec)

   end subroutine apply_testpart_matrix

   !****************************************************************************
   !         The test-particle operator applied to vpa F0 and v^2 F0
   !****************************************************************************
   ! At the point (iky = naky, ikx = 1, tube 1) of every z point among the points
   ! ikxkyzt_first to ikxkyzt_last of kxkyzt_lo, apply the block-tridiagonal matrix
   ! -code_dt C^{ab} of the chunk of blocks aa_blcs, bb_blcs and cc_blcs (see
   ! init_fp_diffmatrix) to vpa F0a and v^2 F0a, and store the real part in tp_f0moments.
   ! The block row iv has the blocks aa_blcs(iv) (column iv - 1), bb_blcs(iv) (column iv)
   ! and cc_blcs(iv) (column iv + 1), as blockmatrix in init_fp_diffmatrix.
   !****************************************************************************
   subroutine get_tp_f0moments(ikxkyzt_first, ikxkyzt_last)

      use grids_species, only: nspec
      use grids_velocity, only: vpa, nvpa, mu, nmu
      use grids_kxky, only: naky
      use geometry, only: bmag
      use parallelisation_layouts, only: kxkyzt_lo, iz_idx, idx

      implicit none

      integer, intent(in) :: ikxkyzt_first, ikxkyzt_last

      integer :: ikxkyzt, iz, is, isb, iv, imu, im, ia
      real, dimension(nvpa, nmu, 2) :: x
      complex :: prod

      !-------------------------------------------------------------------------

      ia = 1
      do ikxkyzt = ikxkyzt_first, ikxkyzt_last
         iz = iz_idx(kxkyzt_lo, ikxkyzt)
         if (ikxkyzt /= idx(kxkyzt_lo, naky, 1, iz, 1)) cycle
         do is = 1, nspec
            do iv = 1, nvpa
               do imu = 1, nmu
                  x(iv, imu, 1) = vpa(iv) * mw(iv, imu, iz, is)
                  x(iv, imu, 2) = (vpa(iv)**2 + 2 * bmag(ia, iz) * mu(imu)) * mw(iv, imu, iz, is)
               end do
            end do
            do isb = 1, nspec
               do im = 1, 2
                  do iv = 1, nvpa
                     do imu = 1, nmu
                        prod = 0.
                        if (iv > 1) prod = prod + sum(aa_blcs(iv, imu, :, ikxkyzt, is, isb) * x(iv - 1, :, im))
                        prod = prod + sum(bb_blcs(iv, imu, :, ikxkyzt, is, isb) * x(iv, :, im))
                        if (iv < nvpa) prod = prod + sum(cc_blcs(iv, imu, :, ikxkyzt, is, isb) * x(iv + 1, :, im))
                        tp_f0moments(nmu * (iv - 1) + imu, im, izl(iz), is, isb) = real(prod)
                     end do
                  end do
               end do
            end do
         end do
      end do

   end subroutine get_tp_f0moments

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_fp_conserve

      use linear_solve, only: lu_decomposition
      use grids_time, only: code_dt
      use grids_species, only: nspec
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: ztmax, maxwell_mu, nmu, nvpa, set_vpa_weights
      use grids_kxky, only: naky, nakx
      use parallelisation_layouts, only: kxkyz_lo, iky_idx, ikx_idx, iz_idx, is_idx, it_idx
      use parallelisation_layouts, only: kxkyzt_lo
      use initialise_redistribute, only: kxkyz2kxkyzt
      use redistribute, only: gather, scatter
      use arrays_distribution_function, only: gvmu
      use field_equations_fluxtube, only: advance_fields_fluxtube_using_field_equations
      use field_equations_collisions, only: get_fields_by_spec_idx
      use constants, only: pi

      implicit none

      integer :: ikxkyz, ikxkyzt, iky, ikx, iz, is, it, iv, imu, ix, ia, idx1, idx2, il, im, ij, mm, ll, jj
      integer :: il1, im1, ij1, mm1, ll1, jj1
      integer :: il2, im2, ij2, mm2, ll2, jj2, isa, isb
      logical :: conservative_wgts
      real :: dum2

      complex, dimension(:, :, :, :), allocatable :: dum1, dum3
      complex, dimension(:, :, :, :, :), allocatable :: field
      complex, dimension(:), allocatable :: ghrs
      ! Temporary kxkyzt-shaped buffer for the implicit solves; species axis is
      ! local on kxkyzt_lo so we add an explicit nspec dimension.
      complex, dimension(:, :, :, :), allocatable :: gvmu_kxkyzt

      !-------------------------------------------------------------------------

      if (.not. allocated(fp_response)) then
         if (include_field_particle_part) then
            nresponse = 1 + (jmax + 1) * (lmax + 1)**2 * nspec**2
         else
            nresponse = 1
         end if
         allocate (fp_response(nresponse, nresponse, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc)); fp_response = 0.
         allocate (diff_idx(nresponse, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      end if

      allocate (dum1(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dum3(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      allocate (ghrs(nmu * nvpa)); ghrs = 0.

      ia = 1

      ! set wgts to uniform to ensure exact conservation properties
      if (conservative_test_particle_density) then
         conservative_wgts = .true.
         call set_vpa_weights(conservative_wgts)
      end if
      if (exact_conservation_test_particle_matrix) then
         conservative_wgts = .false.
         call set_vpa_weights(conservative_wgts)
      end if

      ! phi response. Solve runs on kxkyzt_lo (where cdiffmat_band lives), then
      ! result is scattered back to kxkyz_lo for the rest of the routine.
      allocate (gvmu_kxkyzt(nvpa, nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec))
      gvmu_kxkyzt = 0.
      do ikxkyzt = kxkyzt_lo%llim_proc, kxkyzt_lo%ulim_proc
         iky = iky_idx(kxkyzt_lo, ikxkyzt)
         ikx = ikx_idx(kxkyzt_lo, ikxkyzt)
         iz = iz_idx(kxkyzt_lo, ikxkyzt)
         do is = 1, nspec
   ghrs = reshape(transpose(spread(ztmax(:,is),2,nmu)*spread(maxwell_mu(1,iz,:,is),1,nvpa)*spread(jm(:,0,iky,ikx,izl(iz),is),1,nvpa)), shape=(/ nmu*nvpa /))
            call solve_test_particle_operator(ghrs, ikxkyzt, is)
            gvmu_kxkyzt(:, :, ikxkyzt, is) = transpose(reshape(ghrs, shape=(/nmu, nvpa/)))
         end do
      end do
      call scatter(kxkyz2kxkyzt, gvmu_kxkyzt, gvmu)
      deallocate (gvmu_kxkyzt)

      ! gvmu contains dhs/dphi
      ! for phi equation, need 1-P[dhs/dphi]
      ! note that get_fields sums over species, as required in response matrix
      call advance_fields_fluxtube_using_field_equations(gvmu, field(:, :, :, :, 1), dum1, dum3, dist='h') 

      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         fp_response(1, 1, ikxkyz) = 1.0 - field(iky, ikx, iz, it, 1)
      end do

      ! field particle operator

      ! collect all tp terms for species a in blockmatrix_sum(isa)
      ! required for tp density conservation to machine precision
      ! to do - avoid operations with blockmatrix or blockmatrix_sum, use band storage
      ! Now distributed on kxkyzt_lo with explicit target species:
      !   blockmatrix_sum(:,:,ikxkyzt,is) = sum_{isb} blockmatrix(:,:,ikxkyzt,is,isb)
      if (field_particle_density_from_test_particle) then
         do ikxkyzt = kxkyzt_lo%llim_proc, kxkyzt_lo%ulim_proc
            do is = 1, nspec
               blockmatrix_sum(:, :, ikxkyzt, is) = blockmatrix(:, :, ikxkyzt, is, is)
               do isb = 1, nspec
                  if (isb == is) cycle
                  blockmatrix_sum(:, :, ikxkyzt, is) = blockmatrix_sum(:, :, ikxkyzt, is) + blockmatrix(:, :, ikxkyzt, is, isb)
               end do
            end do
         end do
      end if

      ! field particle terms
      if (include_field_particle_part) then
         idx1 = 1 ! first row (phi equation)
         do idx2 = 1, (jmax + 1) * (lmax + 1)**2 * nspec ! column indices

            isa = 1 + int((idx2 - 1) / float((jmax + 1) * (lmax + 1)**2))
            ij = 1 + mod(1 + mod(idx2 - 1, (jmax + 1) * (lmax + 1)**2) - 1, jmax + 1)
            il = 1 + int(sqrt(1.0 * (1 + mod(idx2 - 1, (jmax + 1) * (lmax + 1)**2) - ij) / (jmax + 1)))
            im = 1 + (1 + mod(idx2 - 1, (jmax + 1) * (lmax + 1)**2) - ij) / (jmax + 1) - (il - 1)**2
            ll = il - 1
            mm = -ll + im - 1
            jj = ij - 1

            ! get psi response, dh/dpsi_{ikn}, need responses:
            ! dh/dpsi_aa, dh/dpsi_ab, ...
            ! dh/dpsi_bb, dh/dpsi_ba, ...

            ! jj=0, ll=0 term is zero because Delta_j^l = 0 for jj=0, ll=0
            ! optionally replace this term with density conserving term
            ! that ensures conservation of density of the test particle operator to machine precision when using non-uniform mu-grid

            if ((jj == 0) .and. (ll == 0) .and. (field_particle_density_from_test_particle)) then

               ! get testpart_den response for kperp = 0, only need dh_isa / d testpart_den_isa
               ! since this includes all interspecies test-particle operators.
               ! Solve runs on kxkyzt_lo for is = isa only; other species are
               ! left at zero (matching the gvmu = 0. + cycle of the original).
               allocate (gvmu_kxkyzt(nvpa, nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec))
               gvmu_kxkyzt = 0.
               is = isa
               do ikxkyzt = kxkyzt_lo%llim_proc, kxkyzt_lo%ulim_proc
                  iky = iky_idx(kxkyzt_lo, ikxkyzt)
                  ikx = ikx_idx(kxkyzt_lo, ikxkyzt)
                  iz = iz_idx(kxkyzt_lo, ikxkyzt)

                  do iv = 1, nvpa
                     do imu = 1, nmu
                        ghrs(nmu * (iv - 1) + imu) = -code_dt * modmw(iv, imu, iz, is) / modmwnorm(iz)
                     end do
                  end do
            call solve_test_particle_operator(ghrs, ikxkyzt, is)
                  do ix = 1, nvpa
                     gvmu_kxkyzt(ix, :, ikxkyzt, is) = ghrs(nmu * (ix - 1) + 1:nmu * (ix - 1) + nmu)
                  end do

               end do
               gvmu = 0.
               call scatter(kxkyz2kxkyzt, gvmu_kxkyzt, gvmu)
               deallocate (gvmu_kxkyzt)

               ! get phi_isa(dh_is2a/dh_is1), phi_isa(dh_is2a/dh_is2) ...
               call get_fields_by_spec_idx(isa, gvmu, field) ! AVB: check - using by_spec_idx instead of by_spec_mod

               do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                  iky = iky_idx(kxkyz_lo, ikxkyz)
                  ikx = ikx_idx(kxkyz_lo, ikxkyz)
                  iz = iz_idx(kxkyz_lo, ikxkyz)
                  it = it_idx(kxkyz_lo, ikxkyz)
                  fp_response(idx1, 2 + (idx2 - 1) * nspec:1 + idx2 * nspec, ikxkyz) = -field(iky, ikx, iz, it, :)
               end do
            else
               ! get the responses dh_is2a/dh_is1, dh_is2a/dh_is2, dh_is2a/dh_is3 ...
               call get_psi_response(ll, mm, jj, isa, gvmu)

               ! get phi_isa(dh_is2a/dh_is1), phi_isa(dh_is2a/dh_is2) ...
               call get_fields_by_spec_idx(isa, gvmu, field) ! AVB: check - using by_spec_idx instead of by_spec_mod

               do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                  iky = iky_idx(kxkyz_lo, ikxkyz)
                  ikx = ikx_idx(kxkyz_lo, ikxkyz)
                  iz = iz_idx(kxkyz_lo, ikxkyz)
                  it = it_idx(kxkyz_lo, ikxkyz)
                  fp_response(idx1, 2 + (idx2 - 1) * nspec:1 + idx2 * nspec, ikxkyz) = -field(iky, ikx, iz, it, :)
               end do

            end if

         end do

         idx2 = 1 ! first column
         do idx1 = 1, (jmax + 1) * (lmax + 1)**2 * nspec ! row indices

            isa = 1 + int((idx1 - 1) / float((jmax + 1) * (lmax + 1)**2))
            ij = 1 + mod(1 + mod(idx1 - 1, (jmax + 1) * (lmax + 1)**2) - 1, jmax + 1)
            il = 1 + int(sqrt(1.0 * (1 + mod(idx1 - 1, (jmax + 1) * (lmax + 1)**2) - ij) / (jmax + 1)))
            im = 1 + (1 + mod(idx1 - 1, (jmax + 1) * (lmax + 1)**2) - ij) / (jmax + 1) - (il - 1)**2
            ll = il - 1
            mm = -ll + im - 1
            jj = ij - 1

            ! get phi responses, dh_{isa}/dphi, dh_{isb}/dphi, dh_{isc}/dphi ... store in gvmu
            ! Solve runs on kxkyzt_lo (where cdiffmat_band lives), then result
            ! is scattered back to kxkyz_lo for the rest of the routine.
            allocate (gvmu_kxkyzt(nvpa, nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec))
            gvmu_kxkyzt = 0.
            do ikxkyzt = kxkyzt_lo%llim_proc, kxkyzt_lo%ulim_proc
               iky = iky_idx(kxkyzt_lo, ikxkyzt)
               ikx = ikx_idx(kxkyzt_lo, ikxkyzt)
               iz = iz_idx(kxkyzt_lo, ikxkyzt)
               do is = 1, nspec
                  do iv = 1, nvpa
                     do imu = 1, nmu
                        ghrs(nmu * (iv - 1) + imu) = ztmax(iv, is) * maxwell_mu(1, iz, imu, is) * jm(imu, 0, iky, ikx, izl(iz), is)
                     end do
                  end do
                  call solve_test_particle_operator(ghrs, ikxkyzt, is)
                  do ix = 1, nvpa
                     gvmu_kxkyzt(ix, :, ikxkyzt, is) = ghrs(nmu * (ix - 1) + 1:nmu * (ix - 1) + nmu)
                  end do
               end do
            end do
            call scatter(kxkyz2kxkyzt, gvmu_kxkyzt, gvmu)
            deallocate (gvmu_kxkyzt)

            if ((jj == 0) .and. (ll == 0) .and. (field_particle_density_from_test_particle)) then

               ! get C_testpart[isa,isa+isb+isc...][dh_a/dphi], store in field(:,:,:,:,isa)
               field = 0.
               call get_testpart_density(isa, 0, gvmu, field)

               do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                  iky = iky_idx(kxkyz_lo, ikxkyz)
                  ikx = ikx_idx(kxkyz_lo, ikxkyz)
                  iz = iz_idx(kxkyz_lo, ikxkyz)
                  it = it_idx(kxkyz_lo, ikxkyz)
                  fp_response(2 + (idx1 - 1) * nspec:1 + idx1 * nspec, idx2, ikxkyz) = -field(iky, ikx, iz, it, :)
               end do
            else
               ! get psi_{isa,isa}[dh_{isa}/dphi], psi_{isa,isb}[dh_{isb}/dphi], psi_{isa,isc}[dh_{isc}/dphi] ...
               call get_psi(gvmu, field, isa, 0, ll, mm, jj)
               do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                  iky = iky_idx(kxkyz_lo, ikxkyz)
                  ikx = ikx_idx(kxkyz_lo, ikxkyz)
                  iz = iz_idx(kxkyz_lo, ikxkyz)
                  it = it_idx(kxkyz_lo, ikxkyz)
                  fp_response(2 + (idx1 - 1) * nspec:1 + idx1 * nspec, idx2, ikxkyz) = -field(iky, ikx, iz, it, :)
               end do
            end if

         end do

         ! interior entries
         do idx1 = 1, (jmax + 1) * (lmax + 1)**2 * nspec ! row indices
            do idx2 = 1, (jmax + 1) * (lmax + 1)**2 * nspec ! column indices

               isa = 1 + int((idx1 - 1) / float((jmax + 1) * (lmax + 1)**2))
               ij1 = 1 + mod(1 + mod(idx1 - 1, (jmax + 1) * (lmax + 1)**2) - 1, jmax + 1)
               il1 = 1 + int(sqrt(1.0 * (1 + mod(idx1 - 1, (jmax + 1) * (lmax + 1)**2) - ij1) / (jmax + 1)))
               im1 = 1 + (1 + mod(idx1 - 1, (jmax + 1) * (lmax + 1)**2) - ij1) / (jmax + 1) - (il1 - 1)**2
               ll1 = il1 - 1
               mm1 = -ll1 + im1 - 1
               jj1 = ij1 - 1

               isb = 1 + int((idx2 - 1) / float((jmax + 1) * (lmax + 1)**2))
               ij2 = 1 + mod(1 + mod(idx2 - 1, (jmax + 1) * (lmax + 1)**2) - 1, jmax + 1)
               il2 = 1 + int(sqrt(1.0 * (1 + mod(idx2 - 1, (jmax + 1) * (lmax + 1)**2) - ij2) / (jmax + 1)))
               im2 = 1 + (1 + mod(idx2 - 1, (jmax + 1) * (lmax + 1)**2) - ij2) / (jmax + 1) - (il2 - 1)**2
               ll2 = il2 - 1
               mm2 = -ll2 + im2 - 1
               jj2 = ij2 - 1

               if ((jj2 == 0) .and. (ll2 == 0) .and. (field_particle_density_from_test_particle)) then
                  ! Solve runs on kxkyzt_lo for is = isb only; other species are
                  ! left at zero (matching the gvmu = 0. + cycle of the original).
                  allocate (gvmu_kxkyzt(nvpa, nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec))
                  gvmu_kxkyzt = 0.
                  is = isb
                  do ikxkyzt = kxkyzt_lo%llim_proc, kxkyzt_lo%ulim_proc
                     iky = iky_idx(kxkyzt_lo, ikxkyzt)
                     ikx = ikx_idx(kxkyzt_lo, ikxkyzt)
                     iz = iz_idx(kxkyzt_lo, ikxkyzt)
                     do iv = 1, nvpa
                        do imu = 1, nmu
                           ghrs(nmu * (iv - 1) + imu) = -code_dt * modmw(iv, imu, iz, is) / modmwnorm(iz)
                        end do
                     end do
            call solve_test_particle_operator(ghrs, ikxkyzt, is)
                     do ix = 1, nvpa
                        gvmu_kxkyzt(ix, :, ikxkyzt, is) = ghrs(nmu * (ix - 1) + 1:nmu * (ix - 1) + nmu)
                     end do
                  end do
                  gvmu = 0.
                  call scatter(kxkyz2kxkyzt, gvmu_kxkyzt, gvmu)
                  deallocate (gvmu_kxkyzt)
               else
                  gvmu = 0.
                  call get_psi_response(ll2, mm2, jj2, isb, gvmu)
               end if

               if ((jj1 == 0) .and. (ll1 == 0) .and. (field_particle_density_from_test_particle)) then
                  ! get the fields C_testpart[isa,isa+isb+isc...][dh_a/dpsi^ab], C_testpart[isa,isa+isb+isc...][dh_a/dpsi^ac], C_testpart[isa,isa+isb+isc...][dh_a/dpsi^ad] ...
                  call get_testpart_density(isa, isb, gvmu, field)

                  do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                     iky = iky_idx(kxkyz_lo, ikxkyz)
                     ikx = ikx_idx(kxkyz_lo, ikxkyz)
                     iz = iz_idx(kxkyz_lo, ikxkyz)
                     it = it_idx(kxkyz_lo, ikxkyz)
                     ! AVB: check - is index
                     fp_response(2 + (idx1 - 1) * nspec + (isb - 1), 2 + (idx2 - 1) * nspec:1 + idx2 * nspec, ikxkyz) = -field(iky, ikx, iz, it, :)
                  end do
               else
                  ! get fields Q_(isa,isb)[dh_{is}/dpsi^{is,isa}], Q_(isa,isb)[dh_{is}/dpsi^{is,isa}], ...
                  call get_psi(gvmu, field, isa, isb, ll1, mm1, jj1)
                  do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
                     iky = iky_idx(kxkyz_lo, ikxkyz)
                     ikx = ikx_idx(kxkyz_lo, ikxkyz)
                     iz = iz_idx(kxkyz_lo, ikxkyz)
                     it = it_idx(kxkyz_lo, ikxkyz)
                     ! AVB: check - is index
                     fp_response(2 + (idx1 - 1) * nspec + (isb - 1), 2 + (idx2 - 1) * nspec:1 + idx2 * nspec, ikxkyz) = -field(iky, ikx, iz, it, :)
                  end do
               end if
            end do
         end do

         ! add 1 to diagonal
         do idx1 = 2, 1 + (jmax + 1) * (lmax + 1)**2 * nspec**2
            fp_response(idx1, idx1, :) = fp_response(idx1, idx1, :) + 1
         end do

      end if

      ! LU decomposition for response
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         call lu_decomposition(fp_response(:, :, ikxkyz), diff_idx(:, ikxkyz), dum2)
      end do

      conservative_wgts = .false.
      call set_vpa_weights(conservative_wgts)

      deallocate (dum1, dum3, field)

   end subroutine init_fp_conserve

!###############################################################################
!############################ FINISH COLLISION MODEL ###########################
!###############################################################################

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine reset_dt_collisions_fp(collisions_implicit)

      implicit none

      ! Arguments
      logical, intent(in) :: collisions_implicit

      !-------------------------------------------------------------------------

      ! The explicit algorithm multiplies by code_dt at every time step, but the block matrices,
      ! the field-particle terms and the response matrices of the implicit algorithm contain
      ! code_dt, so rebuild them for the new time step (otherwise the collision frequency is
      ! effectively multiplied by dt_old/dt_new). The velocity grids, Legendre polynomials,
      ! Bessel functions and collision frequencies do not depend on code_dt and are kept.
      if (.not. collisions_implicit) return
      if (.not. initialised_fokker_planck) return
      call finish_fp_diffmatrix
      call finish_fp_response
      call finish_deltaj
      call init_fp_diffmatrix
      call init_deltaj_vmu
      call init_fp_conserve

   end subroutine reset_dt_collisions_fp

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine finish_collisions_fp

      implicit none

      call finish_nusDpa
      call finish_fp_diffmatrix
      call finish_fp_response
      call finish_deltaj
      if (allocated(legendre_vpamu)) deallocate (legendre_vpamu)
      if (allocated(jm)) deallocate (jm)
      if (allocated(izl)) deallocate (izl, iz_of_izl, local_z)

      initialised_fokker_planck = .false.

   end subroutine finish_collisions_fp

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine finish_deltaj

      implicit none
      if (allocated(deltaj)) deallocate (deltaj)
      if (allocated(deltaj_tp)) deallocate (deltaj_tp)
      if (allocated(psijnorm)) deallocate (psijnorm)
      if (allocated(mwnorm)) deallocate (mwnorm)
      if (allocated(modmwnorm)) deallocate (modmwnorm)

   end subroutine finish_deltaj

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine finish_fp_diffmatrix

      implicit none

      if (allocated(aa_blcs)) deallocate (aa_blcs)
      if (allocated(bb_blcs)) deallocate (bb_blcs)
      if (allocated(cc_blcs)) deallocate (cc_blcs)
      if (allocated(cdiffmat_band)) deallocate (cdiffmat_band)
      if (allocated(blockmatrix)) deallocate (blockmatrix)
      if (allocated(blockmatrix_sum)) deallocate (blockmatrix_sum)
      if (allocated(tp_f0moments)) deallocate (tp_f0moments)

   end subroutine finish_fp_diffmatrix

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine finish_fp_response

      implicit none

      if (allocated(fp_response)) deallocate (fp_response)
      if (allocated(diff_idx)) deallocate (diff_idx)

   end subroutine finish_fp_response
   
   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine finish_nusDpa

      implicit none

      if (allocated(nuD)) deallocate (nuD)
      if (allocated(nupa)) deallocate (nupa)
      if (allocated(nux)) deallocate (nux)
      if (allocated(mw)) deallocate (mw)
      if (allocated(modmw)) deallocate (modmw)

      if (allocated(velvpamu)) deallocate (velvpamu)

   end subroutine finish_nusDpa
   
!###############################################################################
!################################# CALCULATIONS ################################
!###############################################################################

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine get_psi_response(ll, mm, jj, isa, response)

      ! solve for responses dh_a / dpsi^{aa}, dh_a / dpsi^{ab}, dh_a / dpsi^{ac} ... for all species b, c, ...

      use calculations_finite_differences, only: tridag
      use linear_solve, only: lu_decomposition
      use grids_time, only: code_dt
      use grids_species, only: spec, nspec
      use grids_velocity, only: nmu, nvpa
      use grids_velocity, only: set_vpa_weights
      use parallelisation_layouts, only: kxkyz_lo, kxkyzt_lo
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, it_idx
      use initialise_redistribute, only: kxkyz2kxkyzt
      use redistribute, only: scatter
      use constants, only: pi

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(out) :: response
      integer, intent(in) :: ll, mm, jj, isa
      complex, dimension(:), allocatable :: ghrs
      complex, dimension(:, :, :, :), allocatable :: response_kxkyzt
      integer :: ikxkyzt, iky, ikx, iz, is, ia, iv, imu
      real :: clm

      !-------------------------------------------------------------------------

      allocate (ghrs(nmu * nvpa))

      ia = 1

      clm = sqrt(((2 * ll + 1) * gamma(ll - mm + 1.)) / (4 * pi * gamma(ll + mm + 1.)))

      ! calculate response dh/dpsi_jlm, for unit impulse to psi_jlm.
      ! Solve runs on kxkyzt_lo (matrix is for species `isa`, which lives on
      ! whichever rank owns the kxkyzt point at (iky,ikx,iz,it)). Inner do is
      ! loop covers what was kxkyz_lo's parallelised species axis -- it is the
      ! "isb" index and may differ from isa.
      allocate (response_kxkyzt(nvpa, nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec))
      response_kxkyzt = 0.

      do ikxkyzt = kxkyzt_lo%llim_proc, kxkyzt_lo%ulim_proc
         iky = iky_idx(kxkyzt_lo, ikxkyzt)
         ikx = ikx_idx(kxkyzt_lo, ikxkyzt)
         iz = iz_idx(kxkyzt_lo, ikxkyzt)

         do is = 1, nspec  ! isb (was kxkyz_lo's species index)

            ! supply unit impulse to psi_j^(lm)^{isa,isb}
            do iv = 1, nvpa
               do imu = 1, nmu
                  if (mm == 0) then
                     ghrs(nmu * (iv - 1) + imu) = code_dt * spec(isa)%vnew(is) * clm * legendre_vpamu(ll, mm, iv, imu, izl(iz)) &
                                   * jm(imu, mm, iky, ikx, izl(iz), isa) * (spec(isa)%mass / spec(is)%mass)**(-1.5) * deltaj(ll, jj, isa, is, iv, imu, ia, izl(iz))
                  else if (mm > 0) then
                     ghrs(nmu * (iv - 1) + imu) = code_dt * spec(isa)%vnew(is) * clm * legendre_vpamu(ll, mm, iv, imu, izl(iz)) &
                                   * jm(imu, mm, iky, ikx, izl(iz), isa) * (spec(isa)%mass / spec(is)%mass)**(-1.5) * deltaj(ll, jj, isa, is, iv, imu, ia, izl(iz))
                  else if (mm < 0) then
                     ghrs(nmu * (iv - 1) + imu) = (-1)**mm * code_dt * spec(isa)%vnew(is) * clm * legendre_vpamu(ll, mm, iv, imu, izl(iz)) &
                              * jm(imu, abs(mm), iky, ikx, izl(iz), isa) * (spec(isa)%mass / spec(is)%mass)**(-1.5) * deltaj(ll, jj, isa, is, iv, imu, ia, izl(iz))
                  end if
               end do
            end do

            ! solve for response
            ! need to solve [1 - Deltat C_{test}] dh_a/dpsi^ab = delta_{ab} for dh_a/dpsi^ab. Here, C_{test} includes self collisions and a-b collisions;
            ! note that the matrix is the one of species `isa` (cross-species access), not of the loop's `is`
            call solve_test_particle_operator(ghrs, ikxkyzt, isa)

            do iv = 1, nvpa
               response_kxkyzt(iv, :, ikxkyzt, is) = ghrs(nmu * (iv - 1) + 1:nmu * (iv - 1) + nmu)
            end do

            ! to zero l=1,j=1 term:
            if (no_j1l1) then
               if ((ll == 1) .and. (jj == 1)) then
                  response_kxkyzt(:, :, ikxkyzt, is) = 0.
               end if
            end if

            if (no_j1l2) then
               if ((ll == 2) .and. (jj == 1)) then
                  response_kxkyzt(:, :, ikxkyzt, is) = 0.
               end if
            end if

            if (no_j0l2) then
               if ((ll == 2) .and. (jj == 0)) then
                  response_kxkyzt(:, :, ikxkyzt, is) = 0.
               end if
            end if

            if (spitzer_problem) then
               if (.not. ((isa == 2) .and. (is == 2))) then
                  response_kxkyzt(:, :, ikxkyzt, is) = 0.
               end if
            end if

         end do
      end do

      ! scatter back to kxkyz_lo so the caller sees the same shape as before
      call scatter(kxkyz2kxkyzt, response_kxkyzt, response)

      deallocate (response_kxkyzt)
      deallocate (ghrs)

   end subroutine get_psi_response

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine get_psi(g, fld, isa, isb, ll, mm, jj)

      ! if isb==0:
      ! get the fields psi_aa^lmj(g_a), psi_ab^lmj(g_b), psi_ac^lmj(g_c) ..., for species b, c, ...
      ! if isb/=0, fix the species indices of the operator psi_ab, ie
      ! get the fields psi_ab^lmj(g_a), psi_ab^lmj(g_b), psi_ab^lmj(g_c) ...

      use mp, only: sum_allreduce
      use grids_z, only: nzgrid
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_velocity, only: set_vpa_weights, nvpa, nmu, vpa
      use parallelisation_layouts, only: kxkyz_lo, kxkyzt_lo
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx, idx
      use constants, only: pi
      use grids_species, only: spec, nspec
      use grids_time, only: code_dt

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld
      integer, intent(in) :: isa, isb, ll, mm, jj

      integer :: ikxkyz, ikxkyzt, iky, ikx, iz, it, is, ia, iv
      complex, dimension(:, :), allocatable :: g0
      complex, dimension(:), allocatable :: ghrs
      real :: clm
      logical :: conservative_wgts

      !-------------------------------------------------------------------------

      allocate (g0(nvpa, nmu))
      allocate (ghrs(nmu * nvpa))

      ! Each processor only fills the (ky,kx,z,tube,species) points it owns, and <fld> is summed
      ! over all processors at the end, so the other points must be zero. Otherwise the values that
      ! <fld> holds on entry (e.g. the previous psi moment) are added once for every other
      ! processor, which made the field-particle part depend on the number of processors and blow up
      fld = 0.

      if (conservative_test_particle_density) then
         conservative_wgts = .true.
         call set_vpa_weights(conservative_wgts)
      else
         conservative_wgts = .false.
         call set_vpa_weights(conservative_wgts)
      end if

      clm = (-1)**mm * sqrt(((2 * ll + 1) * gamma(ll + mm + 1.)) / (4 * pi * gamma(ll - mm + 1.)))

      ia = 1

      if (isb == 0) then
         ! Every point writes only its own fld(iky, ikx, iz, it, is), so the points are shared among the threads
         !$omp parallel do default(none) &
         !$omp   firstprivate(kxkyz_lo, isa, isb, ll, mm, jj, clm, ia, nvpa, nmu) &
         !$omp   private(ikxkyz, ikxkyzt, iky, ikx, iz, it, is, iv, g0, ghrs) &
         !$omp   shared(g, fld, mw, jm, legendre_vpamu, deltaj, deltaj_tp, spec, blockmatrix, code_dt, vpa, velvpamu, &
         !$omp          kxkyzt_lo, exact_conservation_block_products, exact_conservation_test_particle_matrix, izl)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            if (mm == 0) then

               if ((exact_conservation_block_products) .and. (((ll == 0) .and. (jj == 1)) .or. ((ll == 1) .and. (jj == 0)))) then
           g0 = clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,-mm,iky,ikx,izl(iz),is),1,nvpa)*legendre_vpamu(ll,-mm,:,:,izl(iz))*deltaj_tp(ll,jj,is,isa,:,:,ia,izl(iz))

               else if ((exact_conservation_test_particle_matrix) .and. ((ll == 1) .and. (jj == 0))) then
                  do iv = 1, nvpa
                     ghrs(nmu * (iv - 1) + 1:nmu * iv) = clm * (spec(is)%mass / spec(isa)%mass)**2 * jm(:, -mm, iky, ikx, izl(iz), is) * g(iv, :, ikxkyz)
                  end do
                  ! blockmatrix is on kxkyzt_lo with explicit (target_is, source_isb).
                  ! NOTE: this branch assumes this rank owns the (iky,ikx,iz,it) point
                  ! in kxkyzt_lo. With exact_conservation_test_particle_matrix enabled, get_psi may need
                  ! to be restructured to loop over kxkyzt_lo with a kxkyz<->kxkyzt
                  ! redistribute, similar to the zgbtrs sites elsewhere in this module.
                  ikxkyzt = idx(kxkyzt_lo, iky, ikx, iz, it)
                  ghrs = matmul(-blockmatrix(:, :, ikxkyzt, is, isa) / code_dt / spec(is)%vnew(isa), ghrs)
                  do iv = 1, nvpa
                     g0(iv, :) = -vpa(iv) * ghrs(nmu * (iv - 1) + 1:nmu * iv)
                  end do
                  if (spec(is)%vnew(isa) == 0) then
                     g0 = 0.
                  end if

               else if ((exact_conservation_test_particle_matrix) .and. ((ll == 0) .and. (jj == 1))) then
                  do iv = 1, nvpa
                     ghrs(nmu * (iv - 1) + 1:nmu * iv) = clm * (spec(is)%mass / spec(isa)%mass)**1.5 * jm(:, -mm, iky, ikx, izl(iz), is) * g(iv, :, ikxkyz)
                  end do
                  ikxkyzt = idx(kxkyzt_lo, iky, ikx, iz, it)
                  ghrs = matmul(-blockmatrix(:, :, ikxkyzt, is, isa) / code_dt / spec(is)%vnew(isa), ghrs)
                  do iv = 1, nvpa
                     g0(iv, :) = velvpamu(iv, :, iz)**2 * ghrs(nmu * (iv - 1) + 1:nmu * iv)
                  end do
                  if (spec(is)%vnew(isa) == 0) then
                     g0 = 0.
                  end if

               else
              g0 = clm*g(:,:,ikxkyz)/mw(:,:,iz,is)*spread(jm(:,-mm,iky,ikx,izl(iz),is),1,nvpa)*legendre_vpamu(ll,-mm,:,:,izl(iz))*deltaj(ll,jj,is,isa,:,:,ia,izl(iz))

               end if
            else if (mm < 0) then
               if ((exact_conservation_block_products) .and. (((ll == 0) .and. (jj == 1)) .or. ((ll == 1) .and. (jj == 0)))) then
                  g0 = (-1)**mm * clm * g(:, :, ikxkyz) / mw(:, :, iz, is) * spread(jm(:, abs(mm), iky, ikx, izl(iz), is), 1, nvpa) &
                       * legendre_vpamu(ll, -mm, :, :, izl(iz)) * deltaj_tp(ll, jj, is, isa, :, :, ia, izl(iz))
               else
                  g0 = (-1)**mm * clm * g(:, :, ikxkyz) / mw(:, :, iz, is) * spread(jm(:, abs(mm), iky, ikx, izl(iz), is), 1, nvpa) &
                       * legendre_vpamu(ll, -mm, :, :, izl(iz)) * deltaj(ll, jj, is, isa, :, :, ia, izl(iz))
               end if

            else if (mm > 0) then
               if ((exact_conservation_block_products) .and. (((ll == 0) .and. (jj == 1)) .or. ((ll == 1) .and. (jj == 0)))) then
                  g0 = clm * g(:, :, ikxkyz) / mw(:, :, iz, is) * spread(jm(:, mm, iky, ikx, izl(iz), is), 1, nvpa) &
                       * legendre_vpamu(ll, -mm, :, :, izl(iz)) * deltaj_tp(ll, jj, is, isa, :, :, ia, izl(iz))
               else
                  g0 = clm * g(:, :, ikxkyz) / mw(:, :, iz, is) * spread(jm(:, mm, iky, ikx, izl(iz), is), 1, nvpa) &
                       * legendre_vpamu(ll, -mm, :, :, izl(iz)) * deltaj(ll, jj, is, isa, :, :, ia, izl(iz))
               end if
            end if

            call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))

            ! to zero l=1,j=1 term:

         end do
         !$omp end parallel do

      else ! isb /= 0:
         ! apply the operator \psi_ab^jlm[] to every species index of g; only used in the calculation of the response matrix
         ! where the species index, is, of g contains the response \delta h_a / \delta \psi^{a,is}
         ! Every point writes only its own fld(iky, ikx, iz, it, is), so the points are shared among the threads
         !$omp parallel do default(none) &
         !$omp   firstprivate(kxkyz_lo, isa, isb, ll, mm, jj, clm, ia, nvpa, nmu) &
         !$omp   private(ikxkyz, ikxkyzt, iky, ikx, iz, it, is, iv, g0, ghrs) &
         !$omp   shared(g, fld, mw, jm, legendre_vpamu, deltaj, deltaj_tp, spec, blockmatrix, code_dt, vpa, velvpamu, &
         !$omp          kxkyzt_lo, exact_conservation_block_products, exact_conservation_test_particle_matrix, izl)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            if (mm == 0) then
               if ((exact_conservation_block_products) .and. (((ll == 0) .and. (jj == 1)) .or. ((ll == 1) .and. (jj == 0)))) then
                  g0 = clm * g(:, :, ikxkyz) / mw(:, :, iz, is) * spread(jm(:, mm, iky, ikx, izl(iz), isb), 1, nvpa) &
                       * legendre_vpamu(ll, -mm, :, :, izl(iz)) * deltaj_tp(ll, jj, isb, isa, :, :, ia, izl(iz))

               else if ((exact_conservation_test_particle_matrix) .and. ((ll == 1) .and. (jj == 0))) then
                  do iv = 1, nvpa
                     ghrs(nmu * (iv - 1) + 1:nmu * iv) = clm * (spec(isb)%mass / spec(isa)%mass)**2 &
                                                         * jm(:, mm, iky, ikx, izl(iz), isb) * g(iv, :, ikxkyz)
                  end do
                  ! Under kxkyzt_lo, all species are local at any (iky,ikx,iz,it):
                  ! the previous "search for ikxkyz_isb at (iky,ikx,iz,it,isb)" nested
                  ! loop is no longer needed -- a single direct lookup suffices.
                  ! Requires this rank to own the kxkyzt point; see TODO note above.
                  ikxkyzt = idx(kxkyzt_lo, iky, ikx, iz, it)
                  ghrs = matmul(-blockmatrix(:, :, ikxkyzt, isb, isa) / code_dt / spec(isb)%vnew(isa), ghrs)
                  do iv = 1, nvpa
                     g0(iv, :) = -vpa(iv) * ghrs(nmu * (iv - 1) + 1:nmu * iv)
                  end do
                  if (spec(isb)%vnew(isa) == 0) then
                     g0 = 0.
                  end if

               else if ((exact_conservation_test_particle_matrix) .and. ((ll == 0) .and. (jj == 1))) then
                  do iv = 1, nvpa
                     ghrs(nmu * (iv - 1) + 1:nmu * iv) = clm * (spec(isb)%mass / spec(isa)%mass)**1.5 &
                                                         * jm(:, mm, iky, ikx, izl(iz), isb) * g(iv, :, ikxkyz)
                  end do
                  ikxkyzt = idx(kxkyzt_lo, iky, ikx, iz, it)
                  ghrs = matmul(-blockmatrix(:, :, ikxkyzt, isb, isa) / code_dt / spec(isb)%vnew(isa), ghrs)
                  do iv = 1, nvpa
                     g0(iv, :) = velvpamu(iv, :, iz)**2 * ghrs(nmu * (iv - 1) + 1:nmu * iv)
                  end do
                  if (spec(isb)%vnew(isa) == 0) then
                     g0 = 0.
                  end if
               else
                  g0 = clm * g(:, :, ikxkyz) / mw(:, :, iz, is) * spread(jm(:, mm, iky, ikx, izl(iz), isb), 1, nvpa) &
                       * legendre_vpamu(ll, -mm, :, :, izl(iz)) * deltaj(ll, jj, isb, isa, :, :, ia, izl(iz))
               end if
            else if (mm < 0) then
               if ((exact_conservation_block_products) .and. (((ll == 0) .and. (jj == 1)) .or. ((ll == 1) .and. (jj == 0)))) then
                  g0 = (-1)**mm * clm * g(:, :, ikxkyz) / mw(:, :, iz, is) * spread(jm(:, abs(mm), iky, ikx, izl(iz), isb), 1, nvpa) &
                       * legendre_vpamu(ll, -mm, :, :, izl(iz)) * deltaj_tp(ll, jj, isb, isa, :, :, ia, izl(iz))
               else
                  g0 = (-1)**mm * clm * g(:, :, ikxkyz) / mw(:, :, iz, is) * spread(jm(:, abs(mm), iky, ikx, izl(iz), isb), 1, nvpa) &
                       * legendre_vpamu(ll, -mm, :, :, izl(iz)) * deltaj(ll, jj, isb, isa, :, :, ia, izl(iz))
               end if

            else if (mm > 0) then
               if ((exact_conservation_block_products) .and. (((ll == 0) .and. (jj == 1)) .or. ((ll == 1) .and. (jj == 0)))) then
                  g0 = clm * g(:, :, ikxkyz) / mw(:, :, iz, is) * spread(jm(:, mm, iky, ikx, izl(iz), isb), 1, nvpa) &
                       * legendre_vpamu(ll, -mm, :, :, izl(iz)) * deltaj_tp(ll, jj, isb, isa, :, :, ia, izl(iz))
               else
                  g0 = clm * g(:, :, ikxkyz) / mw(:, :, iz, is) * spread(jm(:, mm, iky, ikx, izl(iz), isb), 1, nvpa) &
                       * legendre_vpamu(ll, -mm, :, :, izl(iz)) * deltaj(ll, jj, isb, isa, :, :, ia, izl(iz))
               end if
            end if

            call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))

            ! to zero l=1,j=1 term:

         end do
         !$omp end parallel do
      end if

      ! Leave out some of the moments (and, for the Spitzer problem, all but the electron-electron ones); this
      ! zeroes the whole species slices of fld, as before point by point in the loops above
      if ((no_j1l1 .and. ll == 1 .and. jj == 1) .or. (no_j1l2 .and. ll == 2 .and. jj == 1) &
          .or. (no_j0l2 .and. ll == 2 .and. jj == 0)) fld = 0.
      if (spitzer_problem) then
         if (isb == 0) then
            do is = 1, nspec
               if (.not. ((isa == 2) .and. (is == 2))) fld(:, :, :, :, is) = 0.
            end do
         else
            if (.not. ((isa == 2) .and. (isb == 2))) fld = 0.
         end if
      end if

      ! normalise psijs
      if (isb == 0) then
         do is = 1, nspec
            do iz = -nzgrid, nzgrid
               fld(:, :, iz, :, is) = fld(:, :, iz, :, is) / psijnorm(ll, jj, isa, is, iz)
            end do
         end do
      else
         do iz = -nzgrid, nzgrid
            fld(:, :, iz, :, :) = fld(:, :, iz, :, :) / psijnorm(ll, jj, isa, isb, iz)
         end do
      end if

      deallocate (g0)
      deallocate (ghrs)

      conservative_wgts = .false.
      call set_vpa_weights(conservative_wgts)

      call sum_allreduce(fld)

   end subroutine get_psi

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine get_testpart_density(isa, isb, g, fld)

      ! if isb==0:
      ! get the field tp_den_isa(g_isa), store it in fld(:,:,:,:,isa)

      ! if isa=0, fix the species indices of the operator tp_den, ie
      ! get the fields tp_den_isb(g_a), tp_den_isb(g_b), tp_den_isb(g_c) ...

      use mp, only: sum_allreduce
      use grids_z, only: nzgrid
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_velocity, only: set_vpa_weights, nvpa, nmu
      use parallelisation_layouts, only: kxkyz_lo, kxkyzt_lo
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx, idx
      use constants, only: pi
      use grids_time, only: code_dt

      implicit none

      integer, intent(in) :: isa, isb
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld

      integer :: ikxkyz, ikxkyzt, iky, ikx, iz, it, is, ia, iv
      complex, dimension(:, :), allocatable :: g0
      complex, dimension(:), allocatable :: ghrs

      !-------------------------------------------------------------------------

      allocate (g0(nvpa, nmu))
      allocate (ghrs(nmu * nvpa))

      ia = 1

      fld = 0.

      if (isb == 0) then
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            if (is /= isa) cycle
            do iv = 1, nvpa
               ghrs(nmu * (iv - 1) + 1:nmu * iv) = g(iv, :, ikxkyz)
            end do
            ! blockmatrix_sum contains sum of interspecies test particle operators.
            ! Now indexed by (kxkyzt, target_is); requires this rank to own the
            ! kxkyzt point at (iky,ikx,iz,it) (see TODO note in get_psi above).
            ikxkyzt = idx(kxkyzt_lo, iky, ikx, iz, it)
            ghrs = matmul(-blockmatrix_sum(:, :, ikxkyzt, is) / code_dt, ghrs)
            do iv = 1, nvpa
               g0(iv, :) = ghrs(nmu * (iv - 1) + 1:nmu * iv)
            end do
            call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
         end do
      else if (isa == isb) then
         ! apply the operator tp_den_isb[] to every species index of g; only used in the calculation of the response matrix
         ! where the species index, is, of g contains the response \delta h_a / \delta \psi^{a,is}; always given on an x_a grid
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            do iv = 1, nvpa
               ghrs(nmu * (iv - 1) + 1:nmu * iv) = g(iv, :, ikxkyz)
            end do
            ! Under kxkyzt_lo, all species are local at any (iky,ikx,iz,it):
            ! the previous "search for ikxkyz_isb at (iky,ikx,iz,it,isb)" nested
            ! loop is no longer needed -- a single direct lookup suffices.
            ikxkyzt = idx(kxkyzt_lo, iky, ikx, iz, it)
            ghrs = matmul(-blockmatrix_sum(:, :, ikxkyzt, isb) / code_dt, ghrs)
            do iv = 1, nvpa
               g0(iv, :) = ghrs(nmu * (iv - 1) + 1:nmu * iv)
            end do
            call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
         end do
      else
         fld = 0.
      end if

      deallocate (g0)
      deallocate (ghrs)

      call sum_allreduce(fld)

   end subroutine get_testpart_density

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   recursive subroutine gamlow(a, x, gl)

      ! recursive calculation of lower incomplete gamma function for half-integer a > 0

      use constants, only: pi
      use spfunc, only: erf => erf_ext

      implicit none
      
      real, intent(in) :: a
      real, intent(in) :: x
      real, intent(out) :: gl
      real :: glm1

      !-------------------------------------------------------------------------

      if (a == 0.5) then
         gl = sqrt(pi) * erf(sqrt(x))
      else
         call gamlow(a - 1., x, glm1)
         gl = (a - 1.) * glm1 - x**(a - 1.) * exp(-x)
      end if

   end subroutine gamlow

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   recursive subroutine gamup(a, x, gu)

      ! recursive calculation of the upper incomplete gamma function for integer a > 0

      use constants, only: pi
      use spfunc, only: erf => erf_ext

      implicit none
      
      real, intent(in) :: a
      real, intent(in) :: x
      real, intent(out) :: gu
      real :: gum1

      !-------------------------------------------------------------------------

      if (a == 1.0) then
         gu = exp(-x)
      else
         call gamup(a - 1., x, gum1)
         gu = (a - 1.) * gum1 + x**(a - 1.) * exp(-x)
      end if

   end subroutine gamup

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine calc_delta0(xa, jj, ll, isa, isb, delt0)

      ! calculate Delta0^{j,l,ab}(xa) (on xa grid)
      ! j and l denote the degree and index of the associated laguerre polynomial

      use grids_species, only: spec
      use constants, only: pi

      implicit none

      integer, intent(in) :: jj, ll, isa, isb
      real, intent(in) :: xa
      real, intent(out) :: delt0

      real :: massr, ckjl, xb, gaml1, gaml2, gamu1, gamu2
      integer :: kk

      !-------------------------------------------------------------------------

      massr = spec(isa)%mass / spec(isb)%mass
      xb = xa / sqrt(massr)

      delt0 = 0
      do kk = 0, jj
         ckjl = (-1)**kk * gamma(jj + ll + 0.5 + 1) / (gamma(jj - kk + 1.) * gamma(ll + kk + 0.5 + 1) * gamma(kk + 1.))
         call gamlow(1.5 + ll + kk, xb**2, gaml1)
         call gamlow(2.5 + ll + kk, xb**2, gaml2)
         call gamup(1.+kk, xb**2, gamu1)
         call gamup(2.+kk, xb**2, gamu2)
         delt0 = delt0 + ckjl * ((2 * ll + 1.) * xb**(ll + 2.*kk) * exp(-xb**2) &
                                 - xb * (1.-massr) * (-(ll + 1.) / xb**(ll + 2.) * gaml1 + ll * xb**(ll - 1.) * gamu1) &
                                 - (1./xb**(ll + 1.) * gaml1 + xb**ll * gamu1) &
                                 + massr * xb**2 * ((ll + 1.) * (ll + 2.) / (2 * ll + 3.) * (xb**(-ll - 3.) * gaml2 + xb**ll * gamu1) &
                                                    - ll * (ll - 1.) / (2 * ll - 1.) * (xb**(-ll - 1.) * gaml1 + xb**(ll - 2.) * gamu2)))
      end do
      delt0 = delt0 * 4 * pi / (pi**1.5) * exp(-xa**2) * massr / (2 * ll + 1.)

   end subroutine calc_delta0

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   recursive subroutine calc_deltaj_vmu(jj, nn, ll, isa, isb, deltj)

      ! calculate Delta_j^{ab,l}[x_b^l L_n^{l+0.5}(x_b^2) exp(-x_b^2)](xa) (on x_a grid)
      ! these are normalised, and calculated without the collision frequency
      ! in contrast to Hirshman & Sigmar 1976

      use grids_velocity, only: mu, nmu, vpa, nvpa
      use grids_z, only: nzgrid
      use geometry, only: bmag

      implicit none

      integer, intent(in) :: jj, nn, ll, isa, isb
      real, dimension(:, :, -nzgrid:), intent(out) :: deltj
      real, dimension(nvpa, nmu, -nzgrid:nzgrid) :: deltajm1_n, deltajm1_j
      integer :: iv, imu, iz, ia
      real :: v
      real, dimension(-nzgrid:nzgrid) :: psijm1_n

      !-------------------------------------------------------------------------

      ia = 1
      if (jj == 0) then
         ! Only on the z points of this processor (local_z, see init_local_z), zero elsewhere
         deltj = 0.
         !$omp parallel do default(none) collapse(3) &
         !$omp   firstprivate(nvpa, nmu, nn, ll, isa, isb, ia) &
         !$omp   private(iv, imu, iz, v) &
         !$omp   shared(nzgrid, local_z, vpa, mu, bmag, deltj)
         do iv = 1, nvpa
            do imu = 1, nmu
               do iz = -nzgrid, nzgrid
                  if (.not. local_z(iz)) cycle
                  v = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mu(imu))
                  call calc_delta0(v, nn, ll, isa, isb, deltj(iv, imu, iz))
               end do
            end do
         end do
         !$omp end parallel do
      else
         ! get Delta_[j-1]^{ab,l}[x_b^l L_n^{l+0.5}(x_b^2) exp(-x_b^2)](xa)
         call calc_deltaj_vmu(jj - 1, nn, ll, isa, isb, deltajm1_n)
         ! get Delta_[j-1]^{ab,l}[x_b^l L_[j-1]^{l+0.5}(x_b^2) exp(-x_b^2)](xa)
         call calc_deltaj_vmu(jj - 1, jj - 1, ll, isa, isb, deltajm1_j)
         ! get psi_[j-1]^{ab,l}[x_b^l L_n^{l+0.5}(x_b^2) exp(-x_b^2)](xa)
         call calc_psi_vmu(jj - 1, nn, ll, isa, isb, psijm1_n)
         deltj = deltajm1_n - spread(spread(psijm1_n, 1, nvpa), 2, nmu) * deltajm1_j
      end if

   end subroutine calc_deltaj_vmu

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine vLj_vmu(jj, ll, vLj)

      use grids_velocity, only: mu, nmu, vpa, nvpa
      use grids_z, only: nzgrid
      use geometry, only: bmag

      implicit none

      integer, intent(in) :: jj, ll
      real, dimension(:, :, -nzgrid:), intent(out) :: vLj
      integer :: iv, imu, iz, ia
      real :: v

      !-------------------------------------------------------------------------

      ia = 1
      vLj = 0.
      do iv = 1, nvpa
         do imu = 1, nmu
            do iz = -nzgrid, nzgrid
               if (.not. local_z(iz)) cycle
               v = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mu(imu))
               vLj(iv, imu, iz) = v**ll * associated_laguerre(jj, ll + 1./2., v**2)
            end do
         end do
      end do

   end subroutine vLj_vmu

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   recursive subroutine calc_psi_vmu(jj, nn, ll, isa, isb, psij)

      ! calculate psi_j^{ab,l}[x_b^l L_n^{l+0.5}(x_b^2) exp(-x_b^2)]

      ! have defined deltaj without collision frequency
      ! and normalised everthing to species thermal speeds

      use grids_species, only: spec
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_velocity, only: nmu, nvpa
      use grids_z, only: nzgrid

      implicit none

      integer, intent(in) :: jj, nn, ll, isa, isb
      real, dimension(nvpa, nmu, -nzgrid:nzgrid) :: deltj_j, vLj, vLn !deltj_n
      real, dimension(-nzgrid:nzgrid), intent(out) :: psij
      integer :: iz
      real, dimension(-nzgrid:nzgrid) :: num, den

      !-------------------------------------------------------------------------

      if (jj == 0) then
         if ((ll == 0) .and. (nn == 0)) then ! never used
            num = 1
            den = 1
         else
            ! get delta_j^{ba}(x_a^l L_j^{l+0.5}(x_a^2) exp(-x_a^2)), on x_b grid
            call calc_deltaj_vmu(jj, jj, ll, isb, isa, deltj_j)
            ! multiply by x_b^l L_n(x_b^2) and integrate
            call vLj_vmu(nn, ll, vLn)
            do iz = -nzgrid, nzgrid
               if (.not. local_z(iz)) cycle
               call integrate_vmu(vLn(:, :, iz) * deltj_j(:, :, iz), iz, num(iz)) ! numerator in psijl
            end do

            ! if mb / ma < 1 use self-adjointness to avoid resolution problems
            if (spec(isb)%mass / spec(isa)%mass < 1.) then
               call calc_deltaj_vmu(jj, jj, ll, isa, isb, deltj_j)
               ! need to account for mass ratio normalisation
               ! and collision frequency
               ! \int x_a^l L_j^{l+1/2}(x_a^2) \Delta_j'^{ab}[\tilde{f}_b/F0b * F0b] x_a^2 dx_a
               !   = ma^0.5/mb^0.5 * ma^3/mb^3 \int f_b/F0b \Delta_j^{ba}'[x_a^l L_j^{l+1/2}(x_a^2) \tilde{F0a}] x_b^2 dx_b
               call vLj_vmu(jj, ll, vLj)
               do iz = -nzgrid, nzgrid
                  if (.not. local_z(iz)) cycle
                  call integrate_vmu(spec(isb)%mass**3.5 / spec(isa)%mass**3.5 * vLj(:, :, iz) * deltj_j(:, :, iz), iz, den(iz)) ! denominator in psijl
               end do
            else
               call vLj_vmu(jj, ll, vLj)
               do iz = -nzgrid, nzgrid
                  if (.not. local_z(iz)) cycle
                  call integrate_vmu(vLj(:, :, iz) * deltj_j(:, :, iz), iz, den(iz)) ! denominator in psijl
               end do
            end if

         end if
      else
         ! get delta_j^{ba}(x_a^l L_j^{l+0.5}(x_a^2) exp(-x_a^2)), on x_b grid
         call calc_deltaj_vmu(jj, jj, ll, isb, isa, deltj_j)
         ! multiply by x_b^l L_n(x_b^2) and integrate
         call vLj_vmu(nn, ll, vLn)
         do iz = -nzgrid, nzgrid
            if (.not. local_z(iz)) cycle
            call integrate_vmu(vLn(:, :, iz) * deltj_j(:, :, iz), iz, num(iz)) ! numerator in psijl
         end do

         ! if mb / ma < 1 use self-adjointness to avoid resolution problems
         if (spec(isb)%mass / spec(isa)%mass < 1.) then
            call calc_deltaj_vmu(jj, jj, ll, isa, isb, deltj_j)
            ! need to account for mass ratio normalisation
            ! and collision frequency
            ! \int x_a^l L_j^{l+1/2}(x_a^2) \Delta_j'^{ab}[\tilde{f}_b/F0b * F0b] x_a^2 dx_a
            !   = ma^0.5/mb^0.5 * ma^3/mb^3 \int f_b/F0b \Delta_j^{ba}'[x_a^l L_j^{l+1/2}(x_a^2) \tilde{F0a}] x_b^2 dx_b
            call vLj_vmu(jj, ll, vLj)
            do iz = -nzgrid, nzgrid
               if (.not. local_z(iz)) cycle
               call integrate_vmu(spec(isb)%mass**3.5 / spec(isa)%mass**3.5 * vLj(:, :, iz) * deltj_j(:, :, iz), iz, den(iz)) ! denominator in psijl
            end do
         else
            call vLj_vmu(jj, ll, vLj)
            do iz = -nzgrid, nzgrid
               if (.not. local_z(iz)) cycle
               call integrate_vmu(vLj(:, :, iz) * deltj_j(:, :, iz), iz, den(iz)) ! denominator in psijl
            end do
         end if

      end if

      ! Only on the z points of this processor (local_z, see init_local_z), zero elsewhere
      psij = 0.
      do iz = -nzgrid, nzgrid
         if (.not. local_z(iz)) cycle
         psij(iz) = num(iz) / den(iz)
      end do

   end subroutine calc_psi_vmu

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   elemental function associated_laguerre(n, alpha, x)

      integer, intent(in) :: n
      real, intent(in) :: x
      real, intent(in) :: alpha
      integer :: k
      real :: associated_laguerre, p, p1, p2

      !-------------------------------------------------------------------------

      p1 = dble(1.0)
      p2 = dble(1.0) + alpha - x

      if (n == 0) then
         associated_laguerre = p1
         return
      else if (n == 1) then
         associated_laguerre = p2
         return
      end if

      do k = 2, n
         p = ((dble(2.0) * k - dble(1.0) + alpha - x) * p2 - (k - dble(1.0) + alpha) * p1) / k
         p1 = p2
         p2 = p
      end do

      associated_laguerre = p

   end function associated_laguerre

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   elemental function associated_legendre(l, m, x)

      integer, intent(in) :: l, m
      double precision, intent(in) :: x
      integer :: k
      double precision :: associated_legendre, p, p1, p2, fac
      double precision :: pi

      !-------------------------------------------------------------------------

      pi = 3.14159265359

      ! to start the recursion, use that P_l^m = 0 for l < abs(m)
      ! and P_l^l = (-1)^l*(2l-1)!!(1-x^2)^(l/2)
      ! where (2l-1)!! = 2**l * Gamma(l+0.5) / sqrt(pi)
      p1 = 0.
      p2 = (-1)**abs(m) * 2**abs(m) * gamma(abs(m) + 0.5) / sqrt(pi) * (1.-x**2)**(abs(m) / 2.)

      if (abs(m) > l) then
         associated_legendre = 0.
         return
      end if

      if (l == 0) then
         associated_legendre = 1.
         return
      end if

      if (l == m) then
         associated_legendre = p2
         return
      end if

      if (l == -m) then
         fac = (-1)**abs(m) * gamma(l - abs(m) + 1.) / gamma(l + abs(m) + 1.)
         associated_legendre = p2 * fac
         return
      end if

      do k = abs(m) + 1, l
         p = ((dble(2.0) * k - dble(1.0)) * x * p2 - (k - dble(1.0) + abs(m)) * p1) / (k - abs(m))
         p1 = p2
         p2 = p
      end do

      if (m < 0) then
         fac = (-1)**abs(m) * gamma(l - abs(m) + 1.) / gamma(l + abs(m) + 1.)
         p = p * fac
      end if

      associated_legendre = p

   end function associated_legendre
   
!###############################################################################
!######################### ADVANCE COLLISIONS: EXPLICIT ########################
!###############################################################################

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine advance_collisions_fp_explicit(g, phi, bpar, gke_rhs)

      use mp, only: mp_abort
      use timers, only: region_start, region_end
      use redistribute, only: scatter, gather
      use grids_time, only: code_dt
      use grids_z, only: nzgrid, ntubes
      use parameters_physics, only: fphi
      use parameters_physics, only: full_flux_surface
      use grids_kxky, only: naky, nakx
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: set_vpa_weights
      use parallelisation_layouts, only: vmu_lo, kxkyz_lo
      use parallelisation_layouts, only: is_idx, iky_idx, ikx_idx, iz_idx
      use initialise_redistribute, only: kxkyz2vmu
      use arrays_distribution_function, only: gvmu
      use calculations_tofrom_ghf, only: g_to_h

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gke_rhs

      integer :: is, ikxkyz, imu, iv, ivmu, ikx, iky, iz, ia
      logical :: conservative_wgts

      complex, dimension(:, :, :, :, :), allocatable :: tmp_vmulo

      complex, dimension(:, :, :), allocatable :: mucoll_fp
      complex, dimension(:, :, :), allocatable :: coll_fp

      !-------------------------------------------------------------------------

      ia = 1

      if (full_flux_surface) then
         call mp_abort("collisions not currently supported for full_flux_surface=T.  Aborting.")
      end if

      call region_start('collisions')

      allocate (tmp_vmulo(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      ! want exact conservation properties for collision operator
      conservative_wgts = .true.
      call set_vpa_weights(conservative_wgts)

      ! switch from g = <f> to h = f + Z*e*phi/T * F0
      tmp_vmulo = g
      call g_to_h(tmp_vmulo, phi, bpar, fphi)

      ! remap so that (vpa,mu) local
      call region_start('collisions_redistribute')
      call scatter(kxkyz2vmu, tmp_vmulo, gvmu)
      call region_end('collisions_redistribute')

      ia = 1

      ! take vpa derivatives
      allocate (coll_fp(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc)); coll_fp = 0.0
      allocate (mucoll_fp(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc)); mucoll_fp = 0.0

      if (conservative_test_particle_density) then
         !$omp parallel default(none) &
         !$omp firstprivate(kxkyz_lo, nmu, nvpa, ia, include_vpa_derivatives, include_mu_derivatives, scale_flr_terms) &
         !$omp private(ikxkyz, iky, ikx, iz, is, imu, iv) &
         !$omp shared(gvmu, coll_fp, mucoll_fp)
         !$omp do
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            if (include_vpa_derivatives) then
               do imu = 1, nmu
                  call vpa_differential_operator_fp_conservative(gvmu(:, :, ikxkyz), coll_fp(:, :, ikxkyz), imu, iz, is, ia)
               end do
            end if
            if (include_mu_derivatives) then
               do iv = 1, nvpa
                  call mu_differential_operator_fp_conservative(gvmu(:, :, ikxkyz), mucoll_fp(:, :, ikxkyz), iv, iz, is, ia, iky, ikx, scale_flr_terms)
               end do
            end if
            gvmu(:, :, ikxkyz) = coll_fp(:, :, ikxkyz) + mucoll_fp(:, :, ikxkyz)
         end do
         !$omp end do
         !$omp end parallel
      else
         !$omp parallel default(none) &
         !$omp firstprivate(kxkyz_lo, nmu, nvpa, ia, include_vpa_derivatives, include_mu_derivatives, scale_flr_terms) &
         !$omp private(ikxkyz, iky, ikx, iz, is, imu, iv) &
         !$omp shared(gvmu, coll_fp, mucoll_fp)
         !$omp do
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            if (include_vpa_derivatives) then
               do imu = 1, nmu
                  call vpa_differential_operator_fp(gvmu(:, :, ikxkyz), coll_fp(:, :, ikxkyz), imu, iz, is, ia)
               end do
            end if
            if (include_mu_derivatives) then
               do iv = 1, nvpa
                  call mu_differential_operator_fp(gvmu(:, :, ikxkyz), mucoll_fp(:, :, ikxkyz), iv, iz, is, ia, iky, ikx, scale_flr_terms)
               end do
            end if
            gvmu(:, :, ikxkyz) = coll_fp(:, :, ikxkyz) + mucoll_fp(:, :, ikxkyz)
         end do
         !$omp end do
         !$omp end parallel
      end if

      deallocate (coll_fp, mucoll_fp)

      call gather(kxkyz2vmu, gvmu, tmp_vmulo)

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         gke_rhs(:, :, :, :, ivmu) = gke_rhs(:, :, :, :, ivmu) + code_dt * tmp_vmulo(:, :, :, :, ivmu)
      end do

      deallocate (tmp_vmulo)

      ! reset to default integration wgts
      conservative_wgts = .false.
      call set_vpa_weights(conservative_wgts)

      call region_end('collisions')

   end subroutine advance_collisions_fp_explicit

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine vpa_differential_operator_fp(h, Dh, imu, iz, is, ia)

      use grids_velocity, only: nvpa, vpa, dvpa, mu, dmu, nmu, maxwell_mu
      use geometry, only: bmag
      use constants, only: pi
      use grids_species, only: spec

      implicit none

      complex, dimension(:, :), intent(out) :: Dh
      complex, dimension(:, :), intent(in) :: h
      integer, intent(in) :: imu, iz, ia, is
      integer :: iv
      complex :: dmuhp, dmuhm, dvpahp, dvpahm
      real :: xp, xm, vpap, vpam, nupap, nupam, nuDp, nuDm, mwp, mwm

      !-------------------------------------------------------------------------

      iv = 1
      vpap = 0.5 * (vpa(iv) + vpa(iv + 1))
      xp = sqrt(vpap**2 + 2 * bmag(ia, iz) * mu(imu))
      nuDp = spec(is)%vnew(is) * (erf(xp) - (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2)) / xp**3
      nupap = spec(is)%vnew(is) * 2 * (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2) / xp**3
      mwp = exp(-vpap**2) * maxwell_mu(1, iz, imu, is)
      dvpahp = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dvpa

      if (imu == 1) then
         ! first-order accurate, as in implicit routine:
         dmuhp = (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv + 1, imu) / mw(iv + 1, imu, iz, is)) / dmu(imu)
      else if (imu == nmu) then
         ! first-order accurate, as in implicit routine:
         dmuhp = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is)) / dmu(imu - 1)
      else
         dmuhp = ((h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is)) * dmu(imu) / dmu(imu - 1) &
+ (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv + 1, imu) / mw(iv + 1, imu, iz, is)) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu))
      end if
      Dh(iv,imu) = (2*0.5*(nupap*vpap**2 + 2*nuDp*bmag(ia,iz)*mu(imu))*mwp*dvpahp + vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,is)*mw(iv+1,imu,iz,is)*dmuhp)/(2*dvpa)

      iv = nvpa
      vpam = 0.5 * (vpa(iv) + vpa(iv - 1))
      xm = sqrt(vpam**2 + 2 * bmag(ia, iz) * mu(imu))
      nuDm = spec(is)%vnew(is) * (erf(xm) - (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2)) / xm**3
      nupam = spec(is)%vnew(is) * 2 * (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2) / xm**3
      mwm = exp(-vpam**2) * maxwell_mu(1, iz, imu, is)
      dvpahm = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dvpa

      if (imu == 1) then
         ! first-order accurate, as in implicit routine:
         dmuhm = (h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dmu(imu)
      else if (imu == nmu) then
         ! first-order accurate, as in implicit routine:
         dmuhm = (h(iv - 1, imu) / mw(iv - 1, imu, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) / dmu(imu - 1)
      else
         dmuhm = ((h(iv - 1, imu) / mw(iv - 1, imu, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) * dmu(imu) / dmu(imu - 1) &
+ (h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu))
      end if
      Dh(iv,imu) = (-2*0.5*(nupam*vpam**2 + 2*nuDm*bmag(ia,iz)*mu(imu))*mwm*dvpahm - vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,is)*mw(iv-1,imu,iz,is)*dmuhm) / (2*dvpa)

      do iv = 2, nvpa - 1
         ! AVB: interior nodes:
         ! quantities at half-grid-points:
         vpap = 0.5 * (vpa(iv) + vpa(iv + 1))
         vpam = 0.5 * (vpa(iv) + vpa(iv - 1))
         xp = sqrt(vpap**2 + 2 * bmag(ia, iz) * mu(imu))
         xm = sqrt(vpam**2 + 2 * bmag(ia, iz) * mu(imu))
         nuDp = spec(is)%vnew(is) * (erf(xp) - (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2)) / xp**3
         nuDm = spec(is)%vnew(is) * (erf(xm) - (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2)) / xm**3
         nupap = spec(is)%vnew(is) * 2 * (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2) / xp**3
         nupam = spec(is)%vnew(is) * 2 * (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2) / xm**3
         mwp = exp(-vpap**2) * maxwell_mu(1, iz, imu, is)
         mwm = exp(-vpam**2) * maxwell_mu(1, iz, imu, is)
         dvpahp = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dvpa
         dvpahm = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dvpa

         if (imu == 1) then
            ! or first-order accurate, as in implicit routine:
            dmuhp = (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv + 1, imu) / mw(iv + 1, imu, iz, is)) / dmu(imu)
            dmuhm = (h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dmu(imu)
         else if (imu == nmu) then
            ! or first-order accurate, as in implicit routine:
            dmuhp = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is)) / dmu(imu - 1)
            dmuhm = (h(iv - 1, imu) / mw(iv - 1, imu, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) / dmu(imu - 1)
         else
            dmuhp = ((h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is)) * dmu(imu) / dmu(imu - 1) &
+ (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv + 1, imu) / mw(iv + 1, imu, iz, is)) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu))
            dmuhm = ((h(iv - 1, imu) / mw(iv - 1, imu, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) * dmu(imu) / dmu(imu - 1) &
+ (h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu))
         end if
         Dh(iv, imu) = (2 * 0.5 * (nupap * vpap**2 + 2 * nuDp * bmag(ia, iz) * mu(imu)) * mwp * dvpahp &
                        + vpa(iv + 1) * mu(imu) * nux(iv + 1, imu, iz, is, is) * mw(iv + 1, imu, iz, is) * dmuhp &
                        - 2 * 0.5 * (nupam * vpam**2 + 2 * nuDm * bmag(ia, iz) * mu(imu)) * mwm * dvpahm &
                        - vpa(iv - 1) * mu(imu) * nux(iv - 1, imu, iz, is, is) * mw(iv - 1, imu, iz, is) * dmuhm) / (2 * dvpa)
      end do

   end subroutine vpa_differential_operator_fp

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine mu_differential_operator_fp(h, Dh, iv, iz, is, ia, iky, ikx, scale_flr_terms)

      use grids_velocity, only: nmu, mu, dmu, vpa, dvpa, nvpa, maxwell_vpa
      use geometry, only: bmag
      use grids_species, only: spec
      use arrays, only: kperp2
      use constants, only: pi

      implicit none

      complex, dimension(:, :), intent(in) :: h
      complex, dimension(:, :), intent(out) :: Dh
      integer, intent(in) :: iv, iz, is, ia, iky, ikx
      real, intent(in) :: scale_flr_terms
      complex ::  Dvpah, Dvpah_p, Dvpah_m, Dmuh, Dmuh_m, Dmuh_p, Dmuh1, Dmuh2
      real :: nuDp, nuDm, nupap, nupam, mup, mum, xp, xm, mwm, mwp
      integer :: imu

      !-------------------------------------------------------------------------

      imu = 1
      ! vpa-differential terms:
      if (iv == 1) then
         ! AVB: first order accurate:
         Dvpah = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dvpa
         Dvpah_p = (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv, imu + 1) / mw(iv, imu + 1, iz, is)) / dvpa
      else if (iv == nvpa) then
         ! AVB: first order accurate:
         Dvpah = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dvpa
         Dvpah_p = (h(iv, imu + 1) / mw(iv, imu + 1, iz, is) - h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is)) / dvpa
      else
         Dvpah = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / (2 * dvpa)
         Dvpah_p = (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is)) / (2 * dvpa)
      end if

      ! first mu-derivative term at mu_{i}:
      ! use ghost cell at mu_{0} = 0, where mu*vpa*nux(vpa,mu)*F0 vanishes, dmu(0) = mu(1).
      Dmuh1 = ((vpa(iv) * mu(imu) * nux(iv, imu, iz, is, is) * mw(iv, imu, iz, is) * Dvpah) * dmu(imu) / mu(imu) &
              +(vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,is)*mw(iv,imu+1,iz,is)*Dvpah_p - vpa(iv)*mu(imu)*nux(iv,imu,iz,is,is)*mw(iv,imu,iz,is)*Dvpah)*mu(imu)/dmu(imu)) / (mu(imu)+dmu(imu))
      ! first derivative of h, at mu_{i+1/2}, and at mu_i:
      Dmuh = (h(iv, imu + 1) / mw(iv, imu + 1, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dmu(imu) ! first-order accurate, as used in implicit routine
      ! for second-order accuracy:
      Dmuh_p = (h(iv, imu + 1) / mw(iv, imu + 1, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dmu(imu) ! second-order accurate
      ! quantities at mu_{i+1/2}:
      mup = 0.5 * (mu(imu) + mu(imu + 1))
      mwp = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mup)
      xp = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mup)
      nuDp = spec(is)%vnew(is) * (erf(xp) - (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2)) / xp**3
      nupap = spec(is)%vnew(is) * 2 * (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2) / xp**3
      ! second mu-derivative term at mu_{i}:
      ! use d/dmu[...]_{1} = ([...]_{1+1/2} - [...]_{0})/(dmu_{1}/2+mu(1)), where [...]_{0} is a ghost cell at mu_{0} = 0, with [...]_{0} = 0.
 Dmuh2 = ( (2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh )*dmu(imu)/2./mu(imu) &
               + (2 * (nupap * mup**2 + nuDp * vpa(iv)**2 / (2 * bmag(ia, iz)) * mup) * mwp * Dmuh_p &
               - 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh)*mu(imu)/(dmu(imu)/2.) )/(mu(imu)+dmu(imu)/2.)
      ! add differential terms and gyro-diffusive term:
      Dh(iv, imu) = Dmuh2 + Dmuh1 - scale_flr_terms * 0.5 * kperp2(iky, ikx, ia, iz) * (spec(is)%smz / bmag(ia, iz))**2 &
                * (nupa(iv, imu, iz, is, is) * bmag(ia, iz) * mu(imu) + nuD(iv, imu, iz, is, is) * (vpa(iv)**2 + bmag(ia, iz) * mu(imu))) * h(iv, imu)

      imu = nmu
      ! vpa-differential terms:
      if (iv == 1) then
         ! AVB: first order accurate:
         Dvpah = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dvpa
         Dvpah_m = (h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is) - h(iv, imu - 1) / mw(iv, imu - 1, iz, is)) / dvpa
      else if (iv == nvpa) then
         ! AVB: first order accurate:
         Dvpah = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dvpa
         Dvpah_m = (h(iv, imu - 1) / mw(iv, imu - 1, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) / dvpa
      else
         Dvpah = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / (2 * dvpa)
         Dvpah_m = (h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) / (2 * dvpa)
      end if

      ! first mu-derivative term at mu_{nmu}:
      Dmuh1 = ((vpa(iv) * mu(imu) * nux(iv, imu, iz, is, is) * mw(iv, imu, iz, is) * Dvpah &
                - vpa(iv) * mu(imu - 1) * nux(iv, imu - 1, iz, is, is) * mw(iv, imu - 1, iz, is) * Dvpah_m) * dmu(imu - 1) / dmu(imu - 1) &
        + (-vpa(iv) * mu(imu) * nux(iv, imu, iz, is, is) * mw(iv, imu, iz, is) * Dvpah) * dmu(imu - 1) / dmu(imu - 1)) / (dmu(imu - 1) + dmu(imu - 1))
      ! first derivative of h, at mu_{nmu} and mu_{nmu-1/2}:
      Dmuh_m = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv, imu - 1) / mw(iv, imu - 1, iz, is)) / dmu(imu - 1)
      Dmuh = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv, imu - 1) / mw(iv, imu - 1, iz, is)) / dmu(imu - 1) ! first-order accurate
      ! quantities at mu_{nmu-1/2}:
      mum = 0.5 * (mu(imu) + mu(imu - 1))
      mwm = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mum)
      xm = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mum)
      nuDm = spec(is)%vnew(is) * (erf(xm) - (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2)) / xm**3
      nupam = spec(is)%vnew(is) * 2 * (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2) / xm**3
      ! second mu-derivative term at mu_{nmu}:
      ! use d/dmu[...]_{nmu} = ([...]_{nmu+1} - [...]_{nmu-1/2})/(dmu_{nmu-1}/2+dmu(nmu-1)), where [...]_{nmu+1} is a ghost cell at mu = mu_{nmu} + dmu(nmu-1), with [...]_{nmu+1} = 0.
      Dmuh2 = ((2 * (nupa(iv, imu, iz, is, is) * mu(imu)**2 + nuD(iv, imu, iz, is, is) * vpa(iv)**2 / (2 * bmag(ia, iz)) &
                     *mu(imu))*mw(iv,imu,iz,is)*Dmuh - 2*(nupam*mum**2+nuDm*vpa(iv)**2/(2*bmag(ia,iz))*mum)*mwm*Dmuh_m)*dmu(imu-1)/(dmu(imu-1)/2) &
               + (-2 * (nupa(iv, imu, iz, is, is) * mu(imu)**2 + nuD(iv, imu, iz, is, is) * vpa(iv)**2 / (2 * bmag(ia, iz)) &
                        * mu(imu)) * mw(iv, imu, iz, is) * Dmuh) * dmu(imu - 1) / 2./dmu(imu - 1)) / (dmu(imu - 1) / 2.+dmu(imu - 1))
      ! add differential terms and gyro-diffusive term:
      Dh(iv, imu) = Dmuh2 + Dmuh1 - scale_flr_terms * 0.5 * kperp2(iky, ikx, ia, iz) * (spec(is)%smz / bmag(ia, iz))**2 &
                * (nupa(iv, imu, iz, is, is) * bmag(ia, iz) * mu(imu) + nuD(iv, imu, iz, is, is) * (vpa(iv)**2 + bmag(ia, iz) * mu(imu))) * h(iv, imu)

      do imu = 2, nmu - 1
         ! vpa-differential terms:
         if (iv == 1) then
            ! AVB: first order accurate:
            Dvpah = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dvpa
            Dvpah_p = (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv, imu + 1) / mw(iv, imu + 1, iz, is)) / dvpa
            Dvpah_m = (h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is) - h(iv, imu - 1) / mw(iv, imu - 1, iz, is)) / dvpa
         else if (iv == nvpa) then
            ! AVB: first order accurate:
            Dvpah = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dvpa
            Dvpah_p = (h(iv, imu + 1) / mw(iv, imu + 1, iz, is) - h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is)) / dvpa
            Dvpah_m = (h(iv, imu - 1) / mw(iv, imu - 1, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) / dvpa
         else
            Dvpah = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / (2 * dvpa)
            Dvpah_p = (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is)) / (2 * dvpa)
            Dvpah_m = (h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) / (2 * dvpa)
         end if
         ! first mu-derivative of vpa-derivative term, at mu_{i}:
         Dmuh1 = ((vpa(iv) * mu(imu) * nux(iv, imu, iz, is, is) * mw(iv, imu, iz, is) * Dvpah &
                   - vpa(iv) * mu(imu - 1) * nux(iv, imu - 1, iz, is, is) * mw(iv, imu - 1, iz, is) * Dvpah_m) * dmu(imu) / dmu(imu - 1) &
                  + (vpa(iv) * mu(imu + 1) * nux(iv, imu + 1, iz, is, is) * mw(iv, imu + 1, iz, is) * Dvpah_p &
                  - vpa(iv) * mu(imu) * nux(iv, imu, iz, is, is) * mw(iv, imu, iz, is) * Dvpah) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu))
         ! first mu-derivatives of h, at mu_i, mu_{i+1/2} and mu_{i-1/2}:
         Dmuh = ((h(iv, imu) / mw(iv, imu, iz, is) - h(iv, imu - 1) / mw(iv, imu - 1, iz, is)) * dmu(imu) / dmu(imu - 1) &
                + (h(iv, imu + 1) / mw(iv, imu + 1, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu))
         Dmuh_m = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv, imu - 1) / mw(iv, imu - 1, iz, is)) / dmu(imu - 1)
         Dmuh_p = (h(iv, imu + 1) / mw(iv, imu + 1, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dmu(imu)
         ! quantities at mu_{i+1/2} and mu_{i-1/2}:
         mup = 0.5 * (mu(imu) + mu(imu + 1))
         mum = 0.5 * (mu(imu) + mu(imu - 1))
         mwp = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mup)
         mwm = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mum)
         xp = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mup)
         xm = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mum)
         nuDp = spec(is)%vnew(is) * (erf(xp) - (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2)) / xp**3
         nuDm = spec(is)%vnew(is) * (erf(xm) - (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2)) / xm**3
         nupap = spec(is)%vnew(is) * 2 * (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2) / xp**3
         nupam = spec(is)%vnew(is) * 2 * (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2) / xm**3
         ! second mu-derivative term at mu_{i}:
         Dmuh2 = ( ( 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh &
                    - 2 * (nupam * mum**2 + nuDm * vpa(iv)**2 / (2 * bmag(ia, iz)) * mum) * mwm * Dmuh_m) * dmu(imu) / dmu(imu - 1) &
                 + (2 * (nupap * mup**2 + nuDp * vpa(iv)**2 / (2 * bmag(ia, iz)) * mup) * mwp * Dmuh_p - 2 * (nupa(iv, imu, iz, is, is) * mu(imu)**2 &
                       +nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh )*dmu(imu-1)/dmu(imu) )*2/(dmu(imu-1)+dmu(imu))
         ! add differential terms and gyro-diffusive term:
         Dh(iv,imu) = Dmuh2 + Dmuh1 - scale_flr_terms*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*(nupa(iv,imu,iz,is,is)*bmag(ia,iz)*mu(imu) &
                                                                      + nuD(iv, imu, iz, is, is) * (vpa(iv)**2 + bmag(ia, iz) * mu(imu))) * h(iv, imu)
      end do

   end subroutine mu_differential_operator_fp

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine vpa_differential_operator_fp_conservative(h, Dh, imu, iz, is, ia)

      use grids_velocity, only: nvpa, vpa, dvpa, mu, dmu, nmu, maxwell_mu
      use geometry, only: bmag
      use constants, only: pi
      use grids_species, only: spec

      implicit none

      complex, dimension(:, :), intent(out) :: Dh
      complex, dimension(:, :), intent(in) :: h
      integer, intent(in) :: imu, iz, ia, is
      integer :: iv
      complex :: dmuhp, dmuhm, dvpahp, dvpahm
      real :: xp, xm, vpap, vpam, nupap, nupam, nuDp, nuDm, mwp, mwm

      !-------------------------------------------------------------------------

      iv = 1
      vpap = 0.5 * (vpa(iv) + vpa(iv + 1))
      xp = sqrt(vpap**2 + 2 * bmag(ia, iz) * mu(imu))
      nuDp = spec(is)%vnew(is) * (erf(xp) - (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2)) / xp**3
      nupap = spec(is)%vnew(is) * 2 * (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2) / xp**3
      mwp = exp(-vpap**2) * maxwell_mu(1, iz, imu, is)
      dvpahp = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dvpa

      if (imu == 1) then
         ! first-order accurate, as in implicit routine:
         dmuhp = (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv + 1, imu) / mw(iv + 1, imu, iz, is)) / dmu(imu)
      else if (imu == nmu) then
         ! first-order accurate, as in implicit routine:
         dmuhp = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is)) / dmu(imu - 1)
      else
         dmuhp = ((h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is)) * dmu(imu) / dmu(imu - 1) &
+ (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv + 1, imu) / mw(iv + 1, imu, iz, is)) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu))
      end if
      Dh(iv,imu) = (2*0.5*(nupap*vpap**2 + 2*nuDp*bmag(ia,iz)*mu(imu))*mwp*dvpahp + vpa(iv+1)*mu(imu)*nux(iv+1,imu,iz,is,is)*mw(iv+1,imu,iz,is)*dmuhp)/(2*dvpa)

      iv = nvpa
      vpam = 0.5 * (vpa(iv) + vpa(iv - 1))
      xm = sqrt(vpam**2 + 2 * bmag(ia, iz) * mu(imu))
      nuDm = spec(is)%vnew(is) * (erf(xm) - (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2)) / xm**3
      nupam = spec(is)%vnew(is) * 2 * (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2) / xm**3
      mwm = exp(-vpam**2) * maxwell_mu(1, iz, imu, is)
      dvpahm = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dvpa

      if (imu == 1) then
         ! first-order accurate, as in implicit routine:
         dmuhm = (h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dmu(imu)
      else if (imu == nmu) then
         ! first-order accurate, as in implicit routine:
         dmuhm = (h(iv - 1, imu) / mw(iv - 1, imu, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) / dmu(imu - 1)
      else
         dmuhm = ((h(iv - 1, imu) / mw(iv - 1, imu, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) * dmu(imu) / dmu(imu - 1) &
+ (h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu))
      end if
      Dh(iv,imu) = (-2*0.5*(nupam*vpam**2 + 2*nuDm*bmag(ia,iz)*mu(imu))*mwm*dvpahm - vpa(iv-1)*mu(imu)*nux(iv-1,imu,iz,is,is)*mw(iv-1,imu,iz,is)*dmuhm) / (2*dvpa)

      do iv = 2, nvpa - 1
         ! AVB: interior nodes:
         ! quantities at half-grid-points:
         vpap = 0.5 * (vpa(iv) + vpa(iv + 1))
         vpam = 0.5 * (vpa(iv) + vpa(iv - 1))
         xp = sqrt(vpap**2 + 2 * bmag(ia, iz) * mu(imu))
         xm = sqrt(vpam**2 + 2 * bmag(ia, iz) * mu(imu))
         nuDp = spec(is)%vnew(is) * (erf(xp) - (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2)) / xp**3
         nuDm = spec(is)%vnew(is) * (erf(xm) - (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2)) / xm**3
         nupap = spec(is)%vnew(is) * 2 * (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2) / xp**3
         nupam = spec(is)%vnew(is) * 2 * (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2) / xm**3
         mwp = exp(-vpap**2) * maxwell_mu(1, iz, imu, is)
         mwm = exp(-vpam**2) * maxwell_mu(1, iz, imu, is)
         dvpahp = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dvpa
         dvpahm = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dvpa

         if (imu == 1) then
            ! or first-order accurate, as in implicit routine:
            dmuhp = (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv + 1, imu) / mw(iv + 1, imu, iz, is)) / dmu(imu)
            dmuhm = (h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dmu(imu)
         else if (imu == nmu) then
            ! or first-order accurate, as in implicit routine:
            dmuhp = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is)) / dmu(imu - 1)
            dmuhm = (h(iv - 1, imu) / mw(iv - 1, imu, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) / dmu(imu - 1)
         else
            dmuhp = ((h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is)) * dmu(imu) / dmu(imu - 1) &
+ (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv + 1, imu) / mw(iv + 1, imu, iz, is)) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu))
            dmuhm = ((h(iv - 1, imu) / mw(iv - 1, imu, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) * dmu(imu) / dmu(imu - 1) &
+ (h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu))
         end if

         if (iv == 2) then
            ! assume vpa*mu*nux*F0*dh/dmu vanishes at iv=1, to ensure density conservation
            Dh(iv, imu) = (2 * 0.5 * (nupap * vpap**2 + 2 * nuDp * bmag(ia, iz) * mu(imu)) * mwp * dvpahp &
                           + vpa(iv + 1) * mu(imu) * nux(iv + 1, imu, iz, is, is) * mw(iv + 1, imu, iz, is) * dmuhp &
                           - 2 * 0.5 * (nupam * vpam**2 + 2 * nuDm * bmag(ia, iz) * mu(imu)) * mwm * dvpahm &
                           - 0 * vpa(iv - 1) * mu(imu) * nux(iv - 1, imu, iz, is, is) * mw(iv - 1, imu, iz, is) * dmuhm) / (2 * dvpa)
         else if (iv == nvpa - 1) then
            ! assume vpa*mu*nux*F0*dh/dmu vanishes at iv=nvpa, to ensure density conservation
            Dh(iv, imu) = (2 * 0.5 * (nupap * vpap**2 + 2 * nuDp * bmag(ia, iz) * mu(imu)) * mwp * dvpahp &
                           + 0 * vpa(iv + 1) * mu(imu) * nux(iv + 1, imu, iz, is, is) * mw(iv + 1, imu, iz, is) * dmuhp &
                           - 2 * 0.5 * (nupam * vpam**2 + 2 * nuDm * bmag(ia, iz) * mu(imu)) * mwm * dvpahm &
                           - vpa(iv - 1) * mu(imu) * nux(iv - 1, imu, iz, is, is) * mw(iv - 1, imu, iz, is) * dmuhm) / (2 * dvpa)
         else
            Dh(iv, imu) = (2 * 0.5 * (nupap * vpap**2 + 2 * nuDp * bmag(ia, iz) * mu(imu)) * mwp * dvpahp &
                           + vpa(iv + 1) * mu(imu) * nux(iv + 1, imu, iz, is, is) * mw(iv + 1, imu, iz, is) * dmuhp &
                           - 2 * 0.5 * (nupam * vpam**2 + 2 * nuDm * bmag(ia, iz) * mu(imu)) * mwm * dvpahm &
                           - vpa(iv - 1) * mu(imu) * nux(iv - 1, imu, iz, is, is) * mw(iv - 1, imu, iz, is) * dmuhm) / (2 * dvpa)
         end if

      end do

   end subroutine vpa_differential_operator_fp_conservative

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine mu_differential_operator_fp_conservative(h, Dh, iv, iz, is, ia, iky, ikx, scale_flr_terms)

      use grids_velocity, only: nmu, mu, dmu, vpa, dvpa, nvpa, maxwell_vpa
      use geometry, only: bmag
      use grids_species, only: spec
      use arrays, only: kperp2
      use constants, only: pi

      implicit none

      complex, dimension(:, :), intent(in) :: h
      complex, dimension(:, :), intent(out) :: Dh
      integer, intent(in) :: iv, iz, is, ia, iky, ikx
      real, intent(in) :: scale_flr_terms
      complex ::  Dvpah, Dvpah_p, Dvpah_m, Dmuh, Dmuh_m, Dmuh_p, Dmuh1, Dmuh2
      real :: nuDp, nuDm, nupap, nupam, mup, mum, xp, xm, mwm, mwp
      integer :: imu

      !-------------------------------------------------------------------------

      imu = 1
      ! vpa-differential terms:
      if (iv == 1) then
         ! AVB: first order accurate:
         Dvpah = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dvpa
         Dvpah_p = (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv, imu + 1) / mw(iv, imu + 1, iz, is)) / dvpa
      else if (iv == nvpa) then
         ! AVB: first order accurate:
         Dvpah = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dvpa
         Dvpah_p = (h(iv, imu + 1) / mw(iv, imu + 1, iz, is) - h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is)) / dvpa
      else
         Dvpah = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / (2 * dvpa)
         Dvpah_p = (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is)) / (2 * dvpa)
      end if
      ! first mu-derivative term at mu_{i}:
      ! use ghost cell at mu_{0} = 0, where mu*vpa*nux(vpa,mu)*F0 vanishes, dmu(0) = mu(1)
      ! to ensure conservation of density we approximate as follows:
      Dmuh1 = (vpa(iv)*mu(imu)*nux(iv,imu,iz,is,is)*mw(iv,imu,iz,is)*Dvpah + vpa(iv)*mu(imu+1)*nux(iv,imu+1,iz,is,is)*mw(iv,imu+1,iz,is)*Dvpah_p) / (2.*dmu(imu))

      ! first derivative of h, at mu_{i+1/2}, and at mu_i:
      Dmuh = (h(iv, imu + 1) / mw(iv, imu + 1, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dmu(imu) ! first-order accurate, as used in implicit routine
      ! for second-order accuracy:
      Dmuh_p = (h(iv, imu + 1) / mw(iv, imu + 1, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dmu(imu) ! second-order accurate
      ! quantities at mu_{i+1/2}:
      mup = 0.5 * (mu(imu) + mu(imu + 1))
      mwp = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mup)
      xp = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mup)
      nuDp = spec(is)%vnew(is) * (erf(xp) - (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2)) / xp**3
      nupap = spec(is)%vnew(is) * 2 * (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2) / xp**3
      ! second mu-derivative term at mu_{i}:
      ! use d/dmu[...]_{1} = ([...]_{1+1/2} - [...]_{0})/(dmu_{1}/2+mu(1)), where [...]_{0} is a ghost cell at mu_{0} = 0, with [...]_{0} = 0.
 Dmuh2 = ( (2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh )*dmu(imu)/2./mu(imu) &
               + (2 * (nupap * mup**2 + nuDp * vpa(iv)**2 / (2 * bmag(ia, iz)) * mup) * mwp * Dmuh_p &
        - 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh)*mu(imu)/(dmu(imu)/2.) )/(mu(imu)+dmu(imu)/2.)
      ! add differential terms and gyro-diffusive term:
      Dh(iv, imu) = Dmuh2 + Dmuh1 - scale_flr_terms * 0.5 * kperp2(iky, ikx, ia, iz) * (spec(is)%smz / bmag(ia, iz))**2 &
                * (nupa(iv, imu, iz, is, is) * bmag(ia, iz) * mu(imu) + nuD(iv, imu, iz, is, is) * (vpa(iv)**2 + bmag(ia, iz) * mu(imu))) * h(iv, imu)

      do imu = 2, nmu - 1
         ! vpa-differential terms:
         if (iv == 1) then
            ! AVB: first order accurate:
            Dvpah = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dvpa
            Dvpah_p = (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv, imu + 1) / mw(iv, imu + 1, iz, is)) / dvpa
            Dvpah_m = (h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is) - h(iv, imu - 1) / mw(iv, imu - 1, iz, is)) / dvpa
         else if (iv == nvpa) then
            ! AVB: first order accurate:
            Dvpah = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dvpa
            Dvpah_p = (h(iv, imu + 1) / mw(iv, imu + 1, iz, is) - h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is)) / dvpa
            Dvpah_m = (h(iv, imu - 1) / mw(iv, imu - 1, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) / dvpa
         else
            Dvpah = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / (2 * dvpa)
            Dvpah_p = (h(iv + 1, imu + 1) / mw(iv + 1, imu + 1, iz, is) - h(iv - 1, imu + 1) / mw(iv - 1, imu + 1, iz, is)) / (2 * dvpa)
            Dvpah_m = (h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) / (2 * dvpa)
         end if

         ! first mu-derivative of vpa-derivative term, at mu_{i}:
         if (imu == nmu - 1) then
            ! to ensure conservation of density we assume mu*nux*F0*d(h/F0)/dvpa vanishes at nmu
            Dmuh1 = (-vpa(iv) * mu(imu - 1) * nux(iv, imu - 1, iz, is, is) * mw(iv, imu - 1, iz, is) * Dvpah_m) / (dmu(imu - 1) + dmu(imu - 1))
         else
            Dmuh1 = (-vpa(iv) * mu(imu - 1) * nux(iv, imu - 1, iz, is, is) * mw(iv, imu - 1, iz, is) * Dvpah_m &
                     + vpa(iv) * mu(imu + 1) * nux(iv, imu + 1, iz, is, is) * mw(iv, imu + 1, iz, is) * Dvpah_p) / (dmu(imu - 1) + dmu(imu - 1))
         end if

         ! first mu-derivatives of h, at mu_i, mu_{i+1/2} and mu_{i-1/2}:
          Dmuh   = ( (h(iv,imu)/mw(iv,imu,iz,is)-h(iv,imu-1)/mw(iv,imu-1,iz,is))*dmu(imu)/dmu(imu-1) + (h(iv,imu+1)/mw(iv,imu+1,iz,is)-h(iv,imu)/mw(iv,imu,iz,is))&
                   * dmu(imu - 1) / dmu(imu)) / (dmu(imu - 1) + dmu(imu))
         Dmuh_m = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv, imu - 1) / mw(iv, imu - 1, iz, is)) / dmu(imu - 1)
         Dmuh_p = (h(iv, imu + 1) / mw(iv, imu + 1, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dmu(imu)
         ! quantities at mu_{i+1/2} and mu_{i-1/2}:
         mup = 0.5 * (mu(imu) + mu(imu + 1))
         mum = 0.5 * (mu(imu) + mu(imu - 1))
         mwp = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mup)
         mwm = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mum)
         xp = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mup)
         xm = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mum)
         nuDp = spec(is)%vnew(is) * (erf(xp) - (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2)) / xp**3
         nuDm = spec(is)%vnew(is) * (erf(xm) - (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2)) / xm**3
         nupap = spec(is)%vnew(is) * 2 * (erf(xp) - xp * (2 / sqrt(pi)) * exp(-xp**2)) / (2 * xp**2) / xp**3
         nupam = spec(is)%vnew(is) * 2 * (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2) / xm**3
         ! second mu-derivative term at mu_{i}:
         Dmuh2 = ( ( 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh &
                    - 2 * (nupam * mum**2 + nuDm * vpa(iv)**2 / (2 * bmag(ia, iz)) * mum) * mwm * Dmuh_m) * dmu(imu) / dmu(imu - 1) &
                  + (2 * (nupap * mup**2 + nuDp * vpa(iv)**2 / (2 * bmag(ia, iz)) * mup) * mwp * Dmuh_p &
            - 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh )*dmu(imu-1)/dmu(imu) )*2/(dmu(imu-1)+dmu(imu))
         ! add differential terms and gyro-diffusive term:
         Dh(iv,imu) = Dmuh2 + Dmuh1 - scale_flr_terms*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*(nupa(iv,imu,iz,is,is)*bmag(ia,iz)*mu(imu) &
                                                                      + nuD(iv, imu, iz, is, is) * (vpa(iv)**2 + bmag(ia, iz) * mu(imu))) * h(iv, imu)
      end do

      imu = nmu
      ! vpa-differential terms:
      if (iv == 1) then
         ! AVB: first order accurate:
         Dvpah = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv, imu) / mw(iv, imu, iz, is)) / dvpa
         Dvpah_m = (h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is) - h(iv, imu - 1) / mw(iv, imu - 1, iz, is)) / dvpa
      else if (iv == nvpa) then
         ! AVB: first order accurate:
         Dvpah = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / dvpa
         Dvpah_m = (h(iv, imu - 1) / mw(iv, imu - 1, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) / dvpa
      else
         Dvpah = (h(iv + 1, imu) / mw(iv + 1, imu, iz, is) - h(iv - 1, imu) / mw(iv - 1, imu, iz, is)) / (2 * dvpa)
         Dvpah_m = (h(iv + 1, imu - 1) / mw(iv + 1, imu - 1, iz, is) - h(iv - 1, imu - 1) / mw(iv - 1, imu - 1, iz, is)) / (2 * dvpa)
      end if
      ! first mu-derivative term at mu_{nmu}:
      ! to ensure conservation of density, assume that term is zero beyond nmu
      Dmuh1 = (-vpa(iv) * mu(imu - 1) * nux(iv, imu - 1, iz, is, is) * mw(iv, imu - 1, iz, is) * Dvpah_m) / (dmu(imu - 1) + dmu(imu - 1))

      ! first derivative of h, at mu_{nmu} and mu_{nmu-1/2}:
      Dmuh_m = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv, imu - 1) / mw(iv, imu - 1, iz, is)) / dmu(imu - 1)
      Dmuh = (h(iv, imu) / mw(iv, imu, iz, is) - h(iv, imu - 1) / mw(iv, imu - 1, iz, is)) / dmu(imu - 1) ! first-order accurate
      ! quantities at mu_{nmu-1/2}:
      mum = 0.5 * (mu(imu) + mu(imu - 1))
      mwm = maxwell_vpa(iv, is) * exp(-2 * bmag(ia, iz) * mum)
      xm = sqrt(vpa(iv)**2 + 2 * bmag(ia, iz) * mum)
      nuDm = spec(is)%vnew(is) * (erf(xm) - (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2)) / xm**3
      nupam = spec(is)%vnew(is) * 2 * (erf(xm) - xm * (2 / sqrt(pi)) * exp(-xm**2)) / (2 * xm**2) / xm**3
      ! second mu-derivative term at mu_{nmu}:
      ! to ensure density conservation
      ! use d/dmu[...]_{nmu} = ([...]_{nmu+1/2} - [...]_{nmu-1/2})/(dmu_{nmu-1}/2+dmu(nmu-1)/2), where [...]_{nmu+1/2} is a ghost cell at mu = mu_{nmu} + dmu(nmu-1)/2, with [...]_{nmu+1/2} = 0.
      Dmuh2 = ( ( 2*(nupa(iv,imu,iz,is,is)*mu(imu)**2+nuD(iv,imu,iz,is,is)*vpa(iv)**2/(2*bmag(ia,iz))*mu(imu))*mw(iv,imu,iz,is)*Dmuh &
                 - 2 * (nupam * mum**2 + nuDm * vpa(iv)**2 / (2 * bmag(ia, iz)) * mum) * mwm * Dmuh_m) * dmu(imu - 1) / dmu(imu - 1) &
               + (-2 * (nupa(iv, imu, iz, is, is) * mu(imu)**2 + nuD(iv, imu, iz, is, is) * vpa(iv)**2 / (2 * bmag(ia, iz)) * mu(imu)) &
                  * mw(iv, imu, iz, is) * Dmuh) * dmu(imu - 1) / dmu(imu - 1)) / (dmu(imu - 1) / 2.+dmu(imu - 1) / 2.)
      ! add differential terms and gyro-diffusive term:
      Dh(iv,imu) = Dmuh2 + Dmuh1 - scale_flr_terms*0.5*kperp2(iky,ikx,ia,iz)*(spec(is)%smz/bmag(ia,iz))**2*(nupa(iv,imu,iz,is,is)*bmag(ia,iz)*mu(imu) &
                                                                      + nuD(iv, imu, iz, is, is) * (vpa(iv)**2 + bmag(ia, iz) * mu(imu))) * h(iv, imu)

   end subroutine mu_differential_operator_fp_conservative

!###############################################################################
!######################### ADVANCE COLLISIONS: IMPLICIT ########################
!###############################################################################

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine advance_collisions_fp_implicit(phi, apar, bpar)

      use grids_z, only: nzgrid
      use grids_velocity, only: set_vpa_weights
      use arrays_distribution_function, only: gvmu

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar

      logical :: conservative_wgts

      !-------------------------------------------------------------------------

      if (conservative_test_particle_density) then
         conservative_wgts = .true.
         call set_vpa_weights(conservative_wgts)
      end if
      if (exact_conservation_test_particle_matrix) then
         conservative_wgts = .false.
         call set_vpa_weights(conservative_wgts)
      end if

      call advance_implicit_fp(phi, apar, bpar, gvmu)

   end subroutine advance_collisions_fp_implicit

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine advance_implicit_fp(phi, apar, bpar, g)

      use mp, only: sum_allreduce
      use calculations_finite_differences, only: tridag
      use linear_solve, only: lu_back_substitution
      use grids_time, only: code_dt
      use parameters_physics, only: fphi
      use parameters_physics, only: include_apar, include_bpar
      use grids_species, only: nspec, spec
      use grids_z, only: nzgrid, ntubes
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_velocity, only: nmu, nvpa
      use grids_velocity, only: vpa
      use grids_velocity, only: set_vpa_weights
      use grids_kxky, only: naky, nakx
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: kxkyzt_lo
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx, it_idx
      use initialise_redistribute, only: kxkyz2kxkyzt
      use redistribute, only: gather, scatter
      use calculations_tofrom_ghf, only: g_to_h
      use field_equations_fluxtube, only: advance_fields_fluxtube_using_field_equations
      use constants, only: pi
      use grids_time, only: code_dt

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in out) :: g

      complex, dimension(:, :, :, :, :), allocatable :: flds, flds_solved
      complex, dimension(:, :, :), allocatable :: g_in
      complex, dimension(:), allocatable :: ghrs
      complex, dimension(:, :, :, :, :), allocatable :: field
      ! Temporary kxkyzt-shaped buffer for the implicit zgbtrs solves; the
      ! solves run on whichever rank owns the kxkyzt point at (iky,ikx,iz,it),
      ! after which the result is scattered back to g (kxkyz_lo).
      complex, dimension(:, :, :, :), allocatable :: g_kxkyzt


      integer :: ikxkyz, ikxkyzt, iky, ikx, iz, is, iv, it, ia
      integer :: idx1, ij, il, im, jj, ll, mm, ll1, mm1, jj1, isa, isb
      real :: clm

      real :: spitzer_i1, spitzer_i2, applied_Epar, b_dot_gradz_lnp0, b_dot_gradz_lnT0

      complex, dimension(:, :, :, :), allocatable :: apar_old, bpar_old, phi_zero
      logical :: em_fields

      !-------------------------------------------------------------------------

      ! store input g for use later, as g will be overwritten below
      allocate (g_in(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

      ia = 1

      if (spitzer_problem) then
         ! to solve the Spitzer problem, we add a source term associated with a constant
         ! applied electric field, and pressure and temperature gradients here
         ! all other non-collisional terms are disabled, and Delta t --> \infty

         applied_Epar = 0.01
         b_dot_gradz_lnp0 = 0
         b_dot_gradz_lnT0 = 0.01

         spitzer_i1 = (applied_Epar - b_dot_gradz_lnp0) * spitzer_scale_epar
         spitzer_i2 = b_dot_gradz_lnT0 * spitzer_scale_grad_t

         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            if (is == 1) cycle
            g(:, :, ikxkyz) = g(:, :, ikxkyz) - 1./sqrt(spec(is)%mass) * code_dt * (spread(vpa, 2, nmu) * spitzer_i1 &
                                          + (spread(vpa, 2, nmu) * velvpamu(:, :, iz)**2 - 5./2.*spread(vpa, 2, nmu)) * spitzer_i2) * mw(:, :, iz, is)
         end do
      end if

      g_in = g

      ! h = g + (Z/T) J0 phi F0 + 4 mu (J1/b) bpar F0 uses the B_par of the start of the step, in the inhomogeneous
      ! solution as well as in the final solve. Before, the inhomogeneous solution lacked the B_par part of h, and
      ! its field solve below overwrote B_par (and A_par) with the fields of h_inh, which lacks the Maxwellian
      ! (phi and B_par) part of h; that B_par was then used to convert between g and h.
      em_fields = include_apar .or. include_bpar
      if (em_fields) then
         allocate (apar_old, source=apar); allocate (bpar_old, source=bpar)
         allocate (phi_zero, mold=phi); phi_zero = 0.
         call g_to_h(g, phi_zero, bpar_old, fphi)
      end if

      allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      allocate (flds(naky, nakx, -nzgrid:nzgrid, ntubes, nresponse))

      ! since backwards difference in time, (I-dt*D)h_inh^{n+1} = g^{***}
      ! invert above equation to get h_inh^{n+1}
      ! cdiffmat_band lives on kxkyzt_lo, so the solve is done on kxkyzt and
      ! the RHS / result are redistributed via kxkyz2kxkyzt.
      allocate (g_kxkyzt(nvpa, nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec))
      call gather(kxkyz2kxkyzt, g, g_kxkyzt)
      !$omp parallel default(none) &
      !$omp firstprivate(kxkyzt_lo, nvpa, nmu, nspec) &
      !$omp private(ikxkyzt, is, iv, ghrs) &
      !$omp shared(g_kxkyzt)
      allocate (ghrs(nmu * nvpa))
      !$omp do collapse(2)
      do is = 1, nspec
         do ikxkyzt = kxkyzt_lo%llim_proc, kxkyzt_lo%ulim_proc
            do iv = 1, nvpa
               ghrs(nmu * (iv - 1) + 1:nmu * iv) = g_kxkyzt(iv, :, ikxkyzt, is)
            end do
            call solve_test_particle_operator(ghrs, ikxkyzt, is)
            do iv = 1, nvpa
               g_kxkyzt(iv, :, ikxkyzt, is) = ghrs(nmu * (iv - 1) + 1:nmu * iv)
            end do
         end do
      end do
      !$omp end do
      deallocate (ghrs)
      !$omp end parallel
      call scatter(kxkyz2kxkyzt, g_kxkyzt, g)
      deallocate (g_kxkyzt)

      ! obtain phi^{n+1} and conservation terms using response matrix approach

      ! first get phi_inh^{n+1}
      call advance_fields_fluxtube_using_field_equations(g, phi, apar, bpar, dist='h')
      flds(:, :, :, :, 1) = phi

      if (em_fields) then
         apar = apar_old; bpar = bpar_old
      end if

      ! next get the psi^{s1s2,jlm}_inh^{n+1}
      ! note g contains h_s1_inh, h_s2_inh, ... ie all species
      if (include_field_particle_part) then
         ! layout of field(,,,,:) is phi; jlm0 psi_aa, jlm0 psi_ab, jlm0 psi_ba,  jlm0 psi_bb; jlm1 psi_aa, etc, because we want species to be contiguous
         do idx1 = 1, (jmax + 1) * (lmax + 1)**2 * nspec

            isa = 1 + int((idx1 - 1) / ((jmax + 1) * (lmax + 1)**2)) !1 + mod(idx1-1,(jmax+1)*(lmax+1)**2)
            ij = 1 + mod(1 + mod(idx1 - 1, (jmax + 1) * (lmax + 1)**2) - 1, jmax + 1)
            il = 1 + int(sqrt(1.0 * (1 + mod(idx1 - 1, (jmax + 1) * (lmax + 1)**2) - ij) / (jmax + 1)))
            im = (1 + mod(idx1 - 1, (jmax + 1) * (lmax + 1)**2) - ij) / (jmax + 1) - (il - 1)**2 + 1
            ll = il - 1
            mm = -ll + im - 1
            jj = ij - 1

            if (field_particle_density_from_test_particle .and. (jj == 0) .and. (ll == 0)) then
               ! get the density produced by the combined test particle operator C_{test,isa} = C_{test,isa,isa} + C_{test,isa,isb} + ...
               ! this is stored in field(:,:,:,:,isa), all other entries are zero
               field = 0.
               call get_testpart_density(isa, 0, g, field)
            else
               ! get psi^{isa[is1....isN],jlm}_inh^{n+1}
               call get_psi(g, field, isa, 0, ll, mm, jj)
            end if

            ! add to rhs vector
            flds(:, :, :, :, 2 + (idx1 - 1) * nspec:1 + idx1 * nspec) = field(:, :, :, :, :)

         end do
      end if

      ! AVB: obtain phi^{n+1} and psijlm^{n+1} from response matrix
      ! The response matrix is only applied on the processor that owns is = 1 for each
      ! (ky,kx,z,tube) point, so store the solution in <flds_solved>, which is summed over
      ! all processors below. Note that <flds> contains the inhomogeneous fields on every
      ! processor, so summing <flds> itself would add these for every other processor.
      allocate (flds_solved(naky, nakx, -nzgrid:nzgrid, ntubes, nresponse)); flds_solved = 0.
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         ! all is indices inside ikxkyz super-index have same info
         ! no need to compute multiple times
         is = is_idx(kxkyz_lo, ikxkyz); if (is /= 1) cycle
         call lu_back_substitution(fp_response(:, :, ikxkyz), diff_idx(:, ikxkyz), flds(iky, ikx, iz, it, :))
         flds_solved(iky, ikx, iz, it, :) = flds(iky, ikx, iz, it, :)
      end do

      ! Collect the solution of the response matrix from all processors
      call sum_allreduce(flds_solved)
      flds = flds_solved
      deallocate (flds_solved)

      phi(:, :, :, :) = flds(:, :, :, :, 1)

      g = g_in

      ! RHS is g^{***} + Ze/T*<phi^{n+1}>*F0 + sum_jlm psi_jlm^{n+1}*delta_jl
      ! first two terms added via g_to_h subroutine
      call g_to_h(g, phi, bpar, fphi)

      ! add field particle contribution to RHS:
      if (include_field_particle_part) then
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            do idx1 = 1, (jmax + 1) * (lmax + 1)**2
               ij = 1 + mod(idx1 - 1, jmax + 1)
               il = 1 + int(sqrt(1.0 * (idx1 - ij) / (jmax + 1)))
               im = (idx1 - ij) / (jmax + 1) - (il - 1)**2 + 1
               ll1 = il - 1
               mm1 = -ll1 + im - 1
               jj1 = ij - 1

               if (field_particle_density_from_test_particle .and. (jj1 == 0) .and. (ll1 == 0)) then
                  isb = is
                  g(:,:,ikxkyz) = g(:,:,ikxkyz) - code_dt*flds(iky,ikx,iz,it, 2 + (is-1)*((jmax+1)*(lmax+1)**2*nspec) + (idx1-1)*nspec + (isb-1))&
                                  * modmw(:, :, iz, is) / modmwnorm(iz)
               else
                  clm = sqrt(((2 * ll1 + 1) * gamma(ll1 - mm1 + 1.)) / (4 * pi * gamma(ll1 + mm1 + 1.)))
                  do isb = 1, nspec
                     if (mm1 == 0) then
                g(:,:,ikxkyz) = g(:,:,ikxkyz) + code_dt*spec(is)%vnew(isb)*clm*legendre_vpamu(ll1,mm1,:,:,izl(iz))*spread(jm(:,mm1,iky,ikx,izl(iz),is),1,nvpa) &
                                      * flds(iky, ikx, iz, it, 2 + (is - 1) * ((jmax + 1) * (lmax + 1)**2 * nspec) + (idx1 - 1) * nspec + (isb - 1)) &
                                        * (spec(is)%mass / spec(isb)%mass)**(-1.5) * deltaj(ll1, jj1, is, isb, :, :, ia, izl(iz))
                     else if (mm1 > 0) then
                g(:,:,ikxkyz) = g(:,:,ikxkyz) + code_dt*spec(is)%vnew(isb)*clm*legendre_vpamu(ll1,mm1,:,:,izl(iz))*spread(jm(:,mm1,iky,ikx,izl(iz),is),1,nvpa) &
                                      * flds(iky, ikx, iz, it, 2 + (is - 1) * ((jmax + 1) * (lmax + 1)**2 * nspec) + (idx1 - 1) * nspec + (isb - 1)) &
                                        * (spec(is)%mass / spec(isb)%mass)**(-1.5) * deltaj(ll1, jj1, is, isb, :, :, ia, izl(iz))
                     else if (mm1 < 0) then
 g(:,:,ikxkyz) = g(:,:,ikxkyz) + (-1)**mm1*code_dt*spec(is)%vnew(isb)*clm*legendre_vpamu(ll1,mm1,:,:,izl(iz))*spread(jm(:,abs(mm1),iky,ikx,izl(iz),is),1,nvpa) &
                                      * flds(iky, ikx, iz, it, 2 + (is - 1) * ((jmax + 1) * (lmax + 1)**2 * nspec) + (idx1 - 1) * nspec + (isb - 1)) &
                                        * (spec(is)%mass / spec(isb)%mass)**(-1.5) * deltaj(ll1, jj1, is, isb, :, :, ia, izl(iz))
                     end if
                  end do

               end if

            end do
         end do

      end if

      deallocate (flds)

      ! invert system to get h^{n+1}.
      ! cdiffmat_band lives on kxkyzt_lo, so the solve is done on kxkyzt and
      ! the RHS / result are redistributed via kxkyz2kxkyzt.
      allocate (g_kxkyzt(nvpa, nmu, kxkyzt_lo%llim_proc:kxkyzt_lo%ulim_alloc, nspec))
      call gather(kxkyz2kxkyzt, g, g_kxkyzt)
      !$omp parallel default(none) &
      !$omp firstprivate(kxkyzt_lo, nvpa, nmu, nspec) &
      !$omp private(ikxkyzt, is, iv, ghrs) &
      !$omp shared(g_kxkyzt)
      allocate (ghrs(nmu * nvpa))
      !$omp do collapse(2)
      do is = 1, nspec
         do ikxkyzt = kxkyzt_lo%llim_proc, kxkyzt_lo%ulim_proc
            do iv = 1, nvpa
               ghrs(nmu * (iv - 1) + 1:nmu * iv) = g_kxkyzt(iv, :, ikxkyzt, is)
            end do
            call solve_test_particle_operator(ghrs, ikxkyzt, is)
            do iv = 1, nvpa
               g_kxkyzt(iv, :, ikxkyzt, is) = ghrs(nmu * (iv - 1) + 1:nmu * iv)
            end do
         end do
      end do
      !$omp end do
      deallocate (ghrs)
      !$omp end parallel
      call scatter(kxkyz2kxkyzt, g_kxkyzt, g)
      deallocate (g_kxkyzt)

      ! get g^{n+1} from h^{n+1} and phi^{n+1}
      call g_to_h(g, phi, bpar, -fphi)

      deallocate (g_in)
      if (em_fields) deallocate (apar_old, bpar_old, phi_zero)
      deallocate (field)

   end subroutine advance_implicit_fp

   !****************************************************************************
   !                   Coefficient of the gyro-diffusive term
   !****************************************************************************
   ! The finite-Larmor-radius (classical diffusion) term of the test-particle operator is diagonal in velocity space:
   ! C_test[h] contains -scale_flr_terms * kperp2 * gyro_diffusion_coefficient * h. It is used both in the matrix of the implicit
   ! solve (init_fp_diffmatrix) and in the density correction after the solve (solve_test_particle_operator).
   !****************************************************************************
   real function gyro_diffusion_coefficient(iv, imu, iz, is)

      use grids_velocity, only: vpa, mu
      use grids_species, only: spec, nspec
      use geometry, only: bmag

      implicit none

      integer, intent(in) :: iv, imu, iz, is

      real :: eiediff, eidefl, nusum
      integer :: isb, ia

      ia = 1
      ! intra-species
      nusum = nupa(iv, imu, iz, is, is) * bmag(ia, iz) * mu(imu) + scale_pitch_angle_scattering * nuD(iv, imu, iz, is, is) * (vpa(iv)**2 + bmag(ia, iz) * mu(imu))
      ! inter-species
      do isb = 1, nspec
         if (isb == is) cycle
         if ((is == 2) .and. (isb == 1)) then
            eiediff = scale_energy_diffusion_ei
            eidefl = scale_pitch_angle_scattering_ei
         else
            eiediff = 1
            eidefl = 1
         end if
         nusum = nusum + eiediff * nupa(iv, imu, iz, is, isb) * bmag(ia, iz) * mu(imu) &
                 + eidefl * scale_pitch_angle_scattering * nuD(iv, imu, iz, is, isb) * (vpa(iv)**2 + bmag(ia, iz) * mu(imu))
      end do
      ! (k_perp rho)^2 / k_perp^2 of species is, with bess_fac as in the argument of the Bessel functions
      gyro_diffusion_coefficient = 0.5 * (spec(is)%bess_fac * spec(is)%smz / bmag(ia, iz))**2 * nusum

   end function gyro_diffusion_coefficient

   !****************************************************************************
   !                       Solve [1 - Deltat C_test] h = ghrs
   !****************************************************************************
   ! The velocity-space part of the test-particle operator conserves density analytically, but on the non-uniform mu
   ! grid (without conservative_test_particle_density) the implicit solve changed the density by about 1e-5 of h per
   ! time step (benchmark C2, Nov 2026). With kinetic electrons, quasineutrality amplifies this by 1/(k_perp rho_i)^2,
   ! so that collisions alone were unstable. With density_correction_solve, this density error is removed after the
   ! solve, by adding a multiple of modmw, which has no momentum and no energy moment. The finite-Larmor-radius
   ! (classical diffusion) term of the operator (gyro_diffusion_coefficient) is diagonal in velocity space and changes
   ! the gyrocentre density physically, so only the remainder is removed: with flr_diag = 1 + Deltat (the diagonal
   ! FLR term), the solution x of [flr_diag - Deltat C_v] x = ghrs, where C_v is the velocity-space part, must
   ! satisfy n[flr_diag x] = n[ghrs], where n[.] is the density moment on the grid. (Conserving n[x] instead made the
   ! ion-ion collisions unstable in benchmark C2.) The correction is linear in ghrs, so it is applied in the same
   ! way to the response matrices.
   !****************************************************************************
   subroutine solve_test_particle_operator(ghrs, ikxkyzt, is)

      use parallelisation_layouts, only: kxkyzt_lo, iky_idx, ikx_idx, iz_idx
      use grids_velocity, only: nvpa, nmu
      use grids_time, only: code_dt
      use arrays, only: kperp2
      use calculations_velocity_integrals, only: integrate_vmu
      use banded_lu, only: band_lu_solve

      implicit none

      complex, dimension(:), intent(in out) :: ghrs
      integer, intent(in) :: ikxkyzt, is

      complex, dimension(:, :), allocatable :: hvmu
      real, dimension(:, :), allocatable :: flr_diag
      complex :: density_rhs, density_error
      real :: modmw_density
      integer :: iv, imu, ia, iky, ikx, iz

      iky = iky_idx(kxkyzt_lo, ikxkyzt)
      ikx = ikx_idx(kxkyzt_lo, ikxkyzt)
      iz = iz_idx(kxkyzt_lo, ikxkyzt)

      if (density_correction_solve) then
         ia = 1
         allocate (hvmu(nvpa, nmu), flr_diag(nvpa, nmu))
         hvmu = transpose(reshape(ghrs, shape=(/nmu, nvpa/)))
         call integrate_vmu(hvmu, iz, density_rhs)
         ! 1 plus the diagonal finite-Larmor-radius term of cdiffmat_band (see init_fp_diffmatrix)
         do iv = 1, nvpa
            do imu = 1, nmu
               flr_diag(iv, imu) = 1.
               if (include_test_particle_part) flr_diag(iv, imu) = flr_diag(iv, imu) + code_dt * scale_flr_terms * kperp2(iky, ikx, ia, iz) * gyro_diffusion_coefficient(iv, imu, iz, is)
            end do
         end do
      end if

      call band_lu_solve(nvpa * nmu, nmu + 1, nmu + 1, cdiffmat_band(:, :, ikxkyzt, is), ipiv(:, ikxkyzt, is), ghrs)

      if (density_correction_solve) then
         ! density error: n[flr_diag x] - n[ghrs], removed with modmw
         hvmu = flr_diag * transpose(reshape(ghrs, shape=(/nmu, nvpa/)))
         call integrate_vmu(hvmu, iz, density_error)
         density_error = density_error - density_rhs
         call integrate_vmu(flr_diag * modmw(:, :, iz, is), iz, modmw_density)
         do iv = 1, nvpa
            ghrs(nmu * (iv - 1) + 1:nmu * iv) = ghrs(nmu * (iv - 1) + 1:nmu * iv) &
                                                - modmw(iv, :, iz, is) / modmw_density * density_error
         end do
         deallocate (hvmu, flr_diag)
      end if

   end subroutine solve_test_particle_operator

end module collisions_fokkerplanck
