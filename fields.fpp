module fields

   use common_types, only: eigen_type

   use mpi

   use common_types, only: coupled_alpha_type, gam0_ffs_type

   implicit none

   public :: init_fields, finish_fields
   public :: advance_fields, get_fields
   public :: advance_apar
   public :: get_radial_correction
   public :: enforce_reality_field
   public :: rescale_fields
   public :: get_fields_by_spec, get_fields_by_spec_idx
   public :: gamtot_h, gamtot3_h
   public :: apar_denom
   public :: time_field_solve
   public :: fields_updated
   public :: get_dchidy, get_dchidx
   public :: efac, efacp
   public :: nfields

   private

   real, dimension(:, :, :), allocatable ::  apar_denom
   real :: gamtot_h, gamtot3_h, efac, efacp

   !> arrays allocated/used if simulating a full flux surface
   type(coupled_alpha_type), dimension(:, :, :), allocatable :: gam0_ffs
   type(gam0_ffs_type), dimension(:, :), allocatable :: lu_gam0_ffs
   complex, dimension(:), allocatable :: adiabatic_response_factor

   logical :: fields_updated = .false.
   logical :: fields_initialized = .false.
#ifdef ISO_C_BINDING
   integer :: phi_shared_window = MPI_WIN_NULL
#endif
   logical :: debug = .false.

   integer :: zm
   integer :: nfields

   real, dimension(2, 5) :: time_field_solve = 0.

   interface get_dchidy
      module procedure get_dchidy_4d
      module procedure get_dchidy_2d
   end interface get_dchidy

   interface advance_fields
      module procedure advance_fields_vmu_lo
      module procedure advance_fields_kxkyz_lo
   end interface advance_fields

contains

   subroutine init_fields

      use mp, only: proc0
      use linear_solve, only: lu_decomposition
      use physics_flags, only: full_flux_surface
      use physics_flags, only: include_apar

      implicit none

      debug = debug .and. proc0

      if (full_flux_surface) then
         call init_fields_ffs
      else
         call init_fields_fluxtube
      end if

      nfields = 1
      if (include_apar) nfields = nfields + 1

   end subroutine init_fields

   !> @todo would be tidier if the code related to radial profile variation
   !> were gathered into a separate subroutine or subroutines

   !> init_fields_fluxtube allocates and fills arrays needed during main time advance
   !> loop for the field solve for flux tube simulations
   subroutine init_fields_fluxtube

      use mp, only: sum_allreduce
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, onlY: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use dist_fn_arrays, only: kperp2, dkperp2dr
      use gyro_averages, only: aj0v, aj1v
      use run_parameters, only: fphi
      use run_parameters, only: ky_solve_radial
      use run_parameters, only: maxwellian_normalization
      use physics_parameters, only: tite, nine, beta
      use physics_flags, only: radial_variation
      use species, only: spec, has_electron_species, ion_species
      use stella_geometry, only: dl_over_b, dBdrho, bmag
      use zgrid, only: nzgrid
      use vpamu_grids, only: nvpa, nmu, mu
      use vpamu_grids, only: vpa, vperp2
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: integrate_vmu
      use species, only: spec
      use kt_grids, only: naky, nakx, akx
      use kt_grids, only: zonal_mode
      use physics_flags, only: include_apar, include_bpar
      use physics_flags, only: adiabatic_option_switch
      use physics_flags, only: adiabatic_option_fieldlineavg
      use fields_arrays, only: gamtot, dgamtotdr, gamtot3
      use fields_arrays, only: gamtot13, gamtot31, gamtot33

      implicit none

      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      real :: tmp, wgt
      real, dimension(:, :), allocatable :: g0
      real, dimension(:), allocatable :: g1

      ia = 1
      zm = 0

      ! do not see why this is before fields_initialized check below
      call allocate_arrays

      if (fields_initialized) return
      fields_initialized = .true.

      ! could move these array allocations to allocate_arrays to clean up code
      if (.not. allocated(gamtot)) allocate (gamtot(naky, nakx, -nzgrid:nzgrid)); gamtot = 0.
      if (.not. allocated(gamtot3)) then
         if (.not. has_electron_species(spec) &
             .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            allocate (gamtot3(nakx, -nzgrid:nzgrid)); gamtot3 = 0.
         else
            allocate (gamtot3(1, 1)); gamtot3 = 0.
         end if
      end if
      if (.not. allocated(apar_denom)) then
         if (include_apar) then
            allocate (apar_denom(naky, nakx, -nzgrid:nzgrid)); apar_denom = 0.
         else
            allocate (apar_denom(1, 1, 1)); apar_denom = 0.
         end if
      end if

      if (.not. allocated(gamtot33)) then
         if (include_bpar) then
            allocate (gamtot33(naky, nakx, -nzgrid:nzgrid)); gamtot33 = 0.
         else
            allocate (gamtot33(1, 1, 1)); gamtot33 = 0.
         end if
      end if

      ! gamtot13 and gamtot31 required if fphi!=0 and fbpar!=0
      if (.not. allocated(gamtot13)) then
         if (include_bpar) then
            allocate (gamtot13(naky, nakx, -nzgrid:nzgrid)); gamtot13 = 0.
         else
            allocate (gamtot13(1, 1, 1)); gamtot13 = 0.
         end if
      end if

      if (.not. allocated(gamtot31)) then
         if (include_bpar) then
            allocate (gamtot31(naky, nakx, -nzgrid:nzgrid)); gamtot31 = 0.
         else
            allocate (gamtot31(1, 1, 1)); gamtot31 = 0.
         end if
      end if

      if (radial_variation) then
         if (.not. allocated(dgamtotdr)) allocate (dgamtotdr(naky, nakx, -nzgrid:nzgrid)); dgamtotdr = 0.
      else
         if (.not. allocated(dgamtotdr)) allocate (dgamtotdr(1, 1, 1)); dgamtotdr = 0.
      end if

      if (fphi > epsilon(0.0)) then
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            ! gamtot does not depend on flux tube index,
            ! so only compute for one flux tube index
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            g0 = spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa)
            if (.not. maxwellian_normalization) then
               g0 = g0 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
            end if
            wgt = spec(is)%z * spec(is)%z * spec(is)%dens_psi0 / spec(is)%temp
            call integrate_vmu(g0, iz, tmp)
            gamtot(iky, ikx, iz) = gamtot(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(gamtot)

         gamtot_h = sum(spec%z * spec%z * spec%dens / spec%temp)

         if (radial_variation) then
            allocate (g1(nmu))
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               it = it_idx(kxkyz_lo, ikxkyz)
               ! gamtot does not depend on flux tube index,
               ! so only compute for one flux tube index
               if (it /= 1) cycle
               iky = iky_idx(kxkyz_lo, ikxkyz)
               ikx = ikx_idx(kxkyz_lo, ikxkyz)
               iz = iz_idx(kxkyz_lo, ikxkyz)
               is = is_idx(kxkyz_lo, ikxkyz)
               g1 = aj0v(:, ikxkyz) * aj1v(:, ikxkyz) * (spec(is)%smz)**2 &
                    * (kperp2(iky, ikx, ia, iz) * vperp2(ia, iz, :) / bmag(ia, iz)**2) &
                    * (dkperp2dr(iky, ikx, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                    / (1.0 - aj0v(:, ikxkyz)**2 + 100.*epsilon(0.0))

               g0 = spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa) &
                    * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is) &
                    * (-spec(is)%tprim * (spread(vpa**2, 2, nmu) + spread(vperp2(ia, iz, :), 1, nvpa) - 2.5) &
                       - spec(is)%fprim + (dBdrho(iz) / bmag(ia, iz)) * (1.0 - 2.0 * spread(mu, 1, nvpa) * bmag(ia, iz)) &
                       + spread(g1, 1, nvpa))
               wgt = spec(is)%z * spec(is)%z * spec(is)%dens / spec(is)%temp
               call integrate_vmu(g0, iz, tmp)
               dgamtotdr(iky, ikx, iz) = dgamtotdr(iky, ikx, iz) + tmp * wgt
            end do
            call sum_allreduce(dgamtotdr)

            deallocate (g1)

         end if
         ! avoid divide by zero when kx=ky=0
         ! do not evolve this mode, so value is irrelevant
         if (zonal_mode(1) .and. akx(1) < epsilon(0.) .and. has_electron_species(spec)) then
            gamtot(1, 1, :) = 0.0
            dgamtotdr(1, 1, :) = 0.0
            zm = 1
         end if

         if (.not. has_electron_species(spec)) then
            efac = tite / nine * (spec(ion_species)%dens / spec(ion_species)%temp)
            efacp = efac * (spec(ion_species)%tprim - spec(ion_species)%fprim)
            gamtot = gamtot + efac
            gamtot_h = gamtot_h + efac
            if (radial_variation) dgamtotdr = dgamtotdr + efacp
            if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
               if (zonal_mode(1)) then
                  gamtot3_h = efac / (sum(spec%zt * spec%z * spec%dens))
                  do ikx = 1, nakx
                     ! avoid divide by zero for kx=ky=0 mode,
                     ! which we do not need anyway
                     !if (abs(akx(ikx)) < epsilon(0.)) cycle
                     tmp = 1./efac - sum(dl_over_b(ia, :) / gamtot(1, ikx, :))
                     gamtot3(ikx, :) = 1./(gamtot(1, ikx, :) * tmp)
                  end do
                  if (akx(1) < epsilon(0.)) then
                     gamtot3(1, :) = 0.0
                  end if
               end if
            end if
         end if

         deallocate (g0)

         if (radial_variation .and. ky_solve_radial > 0) call init_radial_field_solve

      end if

      if (include_apar) then
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            ! apar_denom does not depend on flux tube index,
            ! so only compute for one flux tube index
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            g0 = spread(maxwell_vpa(:, is) * vpa**2, 2, nmu) * maxwell_fac(is) &
                 * spread(maxwell_mu(ia, iz, :, is) * aj0v(:, ikxkyz)**2, 1, nvpa)
            wgt = 2.0 * beta * spec(is)%z * spec(is)%z * spec(is)%dens / spec(is)%mass
            call integrate_vmu(g0, iz, tmp)
            apar_denom(iky, ikx, iz) = apar_denom(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(apar_denom)
         apar_denom = apar_denom + kperp2(:, :, ia, :)

         deallocate (g0)
      end if

   end subroutine init_fields_fluxtube

   subroutine init_radial_field_solve
      use mp, only: job
#ifdef ISO_C_BINDING
      use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_intptr_t
      use fields_arrays, only: qn_window, phi_shared
      use mp, only: sgproc0, curr_focus, sharedsubprocs
      use mp, only: scope, real_size, nbytes_real
      use mp, only: split_n_tasks, create_shared_memory_window
      use mpi
#endif
      use run_parameters, only: ky_solve_radial, ky_solve_real
      use species, only: spec, has_electron_species
      use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use zgrid, only: nzgrid, ntubes, nztot
      use species, only: spec
      use kt_grids, only: naky, nakx
      use kt_grids, only: zonal_mode, rho_d_clamped
      use physics_flags, only: adiabatic_option_switch
      use physics_flags, only: adiabatic_option_fieldlineavg
      use linear_solve, only: lu_decomposition, lu_inverse
      use multibox, only: init_mb_get_phi
      use fields_arrays, only: gamtot, dgamtotdr
      use fields_arrays, only: phi_solve, c_mat, theta
      use file_utils, only: runtype_option_switch, runtype_multibox

      implicit none

      integer :: iz, ikx, iky, ia, zmi, naky_r
      real :: dum
      logical :: has_elec, adia_elec
#ifdef ISO_C_BINDING
      integer :: prior_focus, ierr
      integer :: counter, c_lo, c_hi
      integer(c_intptr_t):: cur_pos
      integer(kind=MPI_ADDRESS_KIND) :: win_size
      complex, dimension(:), pointer :: phi_shared_temp
      type(c_ptr) :: cptr
#endif

      complex, dimension(:, :), allocatable :: g0k, g0x

      ia = 1

      naky_r = min(naky, ky_solve_radial)

      has_elec = has_electron_species(spec)
      adia_elec = .not. has_elec .and. zonal_mode(1) &
                  .and. adiabatic_option_switch == adiabatic_option_fieldlineavg

      if (runtype_option_switch == runtype_multibox .and. job == 1 .and. ky_solve_real) then
         call init_mb_get_phi(has_elec, adia_elec, efac, efacp)
      elseif (runtype_option_switch /= runtype_multibox .or. (job == 1 .and. .not. ky_solve_real)) then
         allocate (g0k(1, nakx))
         allocate (g0x(1, nakx))

         if (.not. allocated(phi_solve)) allocate (phi_solve(naky_r, -nzgrid:nzgrid))
#ifdef ISO_C_BINDING
         prior_focus = curr_focus
         call scope(sharedsubprocs)
         !the following is to parallelize the calculation of QN for radial variation sims
         if (debug) write (*, *) 'fields::init_fields::phi_shared_init'
         if (phi_shared_window == MPI_WIN_NULL) then
            win_size = 0
            if (sgproc0) then
               win_size = int(naky * nakx * nztot * ntubes, MPI_ADDRESS_KIND) * 2 * real_size !complex size
            end if

            call create_shared_memory_window(win_size, phi_shared_window, cur_pos)

            cptr = transfer(cur_pos, cptr)

            if (.not. associated(phi_shared)) then
               ! associate array with lower bounds of 1
               call c_f_pointer(cptr, phi_shared_temp, (/naky * nakx * nztot * ntubes/))
               ! now get the correct bounds
               phi_shared(1:naky, 1:nakx, -nzgrid:nzgrid, 1:ntubes) => phi_shared_temp
            end if
            call mpi_win_fence(0, phi_shared_window, ierr)
         end if

         if (debug) write (*, *) 'fields::init_fields::qn_window_init'
         if (qn_window == MPI_WIN_NULL) then
            win_size = 0
            if (sgproc0) then
               win_size = int(nakx * nztot * naky_r, MPI_ADDRESS_KIND) * 4_MPI_ADDRESS_KIND &
                          + int(nakx**2 * nztot * naky_r, MPI_ADDRESS_KIND) * 2 * real_size !complex size
            end if

            call create_shared_memory_window(win_size, qn_window, cur_pos)

            !allocate the memory
            do iky = 1, naky_r
               zmi = 0
               if (iky == 1) zmi = zm !zero mode may or may not be included in matrix
               do iz = -nzgrid, nzgrid
                  if (.not. associated(phi_solve(iky, iz)%zloc)) then
                     allocate (phi_solve(iky, iz)%zloc(nakx - zmi, nakx - zmi))
                     cptr = transfer(cur_pos, cptr)
                     call c_f_pointer(cptr, phi_solve(iky, iz)%zloc, (/nakx - zmi, nakx - zmi/))
                  end if
                  cur_pos = cur_pos + (nakx - zmi)**2 * 2 * nbytes_real
                  if (.not. associated(phi_solve(iky, iz)%idx)) then
                     cptr = transfer(cur_pos, cptr)
                     call c_f_pointer(cptr, phi_solve(iky, iz)%idx, (/nakx - zmi/))
                  end if
                  cur_pos = cur_pos + (nakx - zmi) * 4
               end do
            end do

            call mpi_win_fence(0, qn_window, ierr)
         end if

         call split_n_tasks(nztot * naky_r, c_lo, c_hi)

         call scope(prior_focus)
         counter = 0
#else
         do iky = 1, naky_r
            zmi = 0
            if (iky == 1) zmi = zm !zero mode may or may not be included in matrix
            do iz = -nzgrid, nzgrid
               if (.not. associated(phi_solve(iky, iz)%zloc)) &
                  allocate (phi_solve(iky, iz)%zloc(nakx - zmi, nakx - zmi))
               if (.not. associated(phi_solve(iky, iz)%idx)) &
                  allocate (phi_solve(iky, iz)%idx(nakx - zmi))
            end do
         end do
#endif

         do iky = 1, naky_r
            zmi = 0
            if (iky == 1) zmi = zm !zero mode may or may not be included in matrix
            do iz = -nzgrid, nzgrid
#ifdef ISO_C_BINDING
               counter = counter + 1
               if ((counter >= c_lo) .and. (counter <= c_hi)) then
#endif
                  phi_solve(iky, iz)%zloc = 0.0
                  phi_solve(iky, iz)%idx = 0
                  do ikx = 1 + zmi, nakx
                     g0k(1, :) = 0.0
                     g0k(1, ikx) = dgamtotdr(iky, ikx, iz)

                     call transform_kx2x_unpadded(g0k, g0x)
                     g0x(1, :) = rho_d_clamped * g0x(1, :)
                     call transform_x2kx_unpadded(g0x, g0k)

                     !row column
                     phi_solve(iky, iz)%zloc(:, ikx - zmi) = g0k(1, (1 + zmi):)
                     phi_solve(iky, iz)%zloc(ikx - zmi, ikx - zmi) = phi_solve(iky, iz)%zloc(ikx - zmi, ikx - zmi) &
                                                                     + gamtot(iky, ikx, iz)
                  end do

                  call lu_decomposition(phi_solve(iky, iz)%zloc, phi_solve(iky, iz)%idx, dum)
#ifdef ISO_C_BINDING
               end if
#endif
            end do
         end do

         if (adia_elec) then
            if (.not. allocated(c_mat)) allocate (c_mat(nakx, nakx)); 
            if (.not. allocated(theta)) allocate (theta(nakx, nakx, -nzgrid:nzgrid)); 
            !get C
            do ikx = 1, nakx
               g0k(1, :) = 0.0
               g0k(1, ikx) = 1.0

               call transform_kx2x_unpadded(g0k, g0x)
               g0x(1, :) = (efac + efacp * rho_d_clamped) * g0x(1, :)
               call transform_x2kx_unpadded(g0x, g0k)

               !row column
               c_mat(:, ikx) = g0k(1, :)
            end do

            !get Theta
            do iz = -nzgrid, nzgrid

               !get Theta
               do ikx = 1, nakx
                  g0k(1, :) = 0.0
                  g0k(1, ikx) = dgamtotdr(1, ikx, iz) - efacp

                  call transform_kx2x_unpadded(g0k, g0x)
                  g0x(1, :) = rho_d_clamped * g0x(1, :)
                  call transform_x2kx_unpadded(g0x, g0k)

                  !row column
                  theta(:, ikx, iz) = g0k(1, :)
                  theta(ikx, ikx, iz) = theta(ikx, ikx, iz) + gamtot(1, ikx, iz) - efac
               end do
            end do
         end if
         deallocate (g0k, g0x)
      end if

   end subroutine init_radial_field_solve

   !> init_fields_ffs allocates and fills arrays needed during main time advance
   !> loop for the field solve for full_flux_surface simulations
   subroutine init_fields_ffs

      use species, only: modified_adiabatic_electrons

      implicit none

      if (fields_initialized) return
      fields_initialized = .true.

      !> allocate arrays such as phi that are needed
      !> throughout the simulation
      call allocate_arrays

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
      use vpamu_grids, only: vperp2
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

   subroutine allocate_arrays

      use fields_arrays, only: phi, phi_old
      use fields_arrays, only: apar, apar_old
      use fields_arrays, only: bpar, bpar_old
      use fields_arrays, only: phi_corr_QN, phi_corr_GA
      use fields_arrays, only: apar_corr_QN, apar_corr_GA
      use zgrid, only: nzgrid, ntubes
      use stella_layouts, only: vmu_lo
      use physics_flags, only: radial_variation
      use kt_grids, only: naky, nakx

      implicit none

      if (.not. allocated(phi)) then
         allocate (phi(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi = 0.
      end if
      if (.not. allocated(phi_old)) then
         allocate (phi_old(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi_old = 0.
      end if
      if (.not. allocated(apar)) then
         allocate (apar(naky, nakx, -nzgrid:nzgrid, ntubes))
         apar = 0.
      end if
      if (.not. allocated(apar_old)) then
         allocate (apar_old(naky, nakx, -nzgrid:nzgrid, ntubes))
         apar_old = 0.
      end if
      if (.not. allocated(bpar)) then
         allocate (bpar(naky, nakx, -nzgrid:nzgrid, ntubes))
         bpar = 0.
      end if
      if (.not. allocated(bpar_old)) then
         allocate (bpar_old(naky, nakx, -nzgrid:nzgrid, ntubes))
         bpar_old = 0.
      end if
      if (.not. allocated(phi_corr_QN) .and. radial_variation) then
         allocate (phi_corr_QN(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi_corr_QN = 0.
      end if
      if (.not. allocated(phi_corr_GA) .and. radial_variation) then
         allocate (phi_corr_GA(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         phi_corr_GA = 0.
      end if
      if (.not. allocated(apar_corr_QN) .and. radial_variation) then
         !allocate (apar_corr(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (apar_corr_QN(1, 1, 1, 1))
         apar_corr_QN = 0.
      end if
      if (.not. allocated(apar_corr_GA) .and. radial_variation) then
         !allocate (apar_corr(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (apar_corr_GA(1, 1, 1, 1, 1))
         apar_corr_GA = 0.
      end if

   end subroutine allocate_arrays

   subroutine enforce_reality_field(fin)

!DSO> while most of the modes in the box have reality built in (as we
!     throw out half the kx-ky plane, modes with ky=0 do not have
!     this enforcement built in. In theory this should not be a problem
!     as these modes should be stable, but I made this function (and
!     its relative in the dist file) just in case

      use kt_grids, only: nakx
      use zgrid, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(inout) :: fin

      integer ikx

      fin(1, 1, :, :) = real(fin(1, 1, :, :))
      do ikx = 2, nakx / 2 + 1
         fin(1, ikx, :, :) = 0.5 * (fin(1, ikx, :, :) + conjg(fin(1, nakx - ikx + 2, :, :)))
         fin(1, nakx - ikx + 2, :, :) = conjg(fin(1, ikx, :, :))
      end do

   end subroutine enforce_reality_field

   subroutine advance_fields_vmu_lo(g, phi, apar, dist)

      use mp, only: proc0
      use stella_layouts, only: vmu_lo
      use job_manage, only: time_message
      use redistribute, only: scatter
      use dist_fn_arrays, only: gvmu
      use zgrid, only: nzgrid
      use dist_redistribute, only: kxkyz2vmu
      use run_parameters, only: fields_kxkyz
      use physics_flags, only: full_flux_surface

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar
      character(*), intent(in) :: dist

      if (fields_updated) return

      !> time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')
      !> fields_kxkyz = F is the default
      if (fields_kxkyz) then
         !> first gather (vpa,mu) onto processor for v-space operations
         !> v-space operations are field solve, dg/dvpa, and collisions
         if (debug) write (*, *) 'fields::advance_fields_vmu_lo::scatter'
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' fields_redist')
         call scatter(kxkyz2vmu, g, gvmu)
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' fields_redist')
         !> given gvmu with vpa and mu local, calculate the corresponding fields
         if (debug) write (*, *) 'fields::advance_fields_vmu_lo::get_fields'
         call get_fields(gvmu, phi, apar, dist)
      else
         if (full_flux_surface) then
            if (debug) write (*, *) 'fields::advance_fields_vmu_lo::get_fields_ffs'
            call get_fields_ffs(g, phi, apar)
         else
            call get_fields_vmulo(g, phi, apar, dist)
         end if
      end if

      !> set a flag to indicate that the fields have been updated
      !> this helps avoid unnecessary field solves
      fields_updated = .true.
      !> time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

   end subroutine advance_fields_vmu_lo

   subroutine advance_fields_kxkyz_lo(gvmu, phi, apar, dist)

      use mp, only: proc0
      use stella_layouts, only: kxkyz_lo
      use job_manage, only: time_message
      use redistribute, only: scatter
      use zgrid, only: nzgrid
      use run_parameters, only: fields_kxkyz
      use physics_flags, only: full_flux_surface

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: gvmu
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar
      character(*), intent(in) :: dist

      if (fields_updated) return

      !> time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

      !> given gvmu with vpa and mu local, calculate the corresponding fields
      if (debug) write (*, *) 'dist_fn::advance_stella::advance_fields_kxkyz_lo::get_fields'
      call get_fields(gvmu, phi, apar, dist)

      !> set a flag to indicate that the fields have been updated
      !> this helps avoid unnecessary field solves
      fields_updated = .true.
      !> time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

   end subroutine advance_fields_kxkyz_lo

   subroutine get_fields(g, phi, apar, dist, skip_fsa)

      use mp, only: proc0
      use mp, only: sum_allreduce, mp_abort
      use job_manage, only: time_message
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use dist_fn_arrays, only: kperp2
      use gyro_averages, only: gyro_average
      use physics_flags, only: include_apar
      use run_parameters, only: fphi
      use physics_parameters, only: beta
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: vpa
      use vpamu_grids, only: integrate_vmu
      use species, only: spec

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar
      logical, optional, intent(in) :: skip_fsa
      character(*), intent(in) :: dist
      complex :: tmp

      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      logical :: skip_fsa_local

      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      if (debug) write (*, *) 'dist_fn::advance_stella::get_fields_kxkyzlo'

      ia = 1

      phi = 0.
      if (fphi > epsilon(0.0)) then
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
            wgt = spec(is)%z * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp
         end do
         deallocate (g0)
         call sum_allreduce(phi)
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')

         call get_phi(phi, dist, skip_fsa_local)

      end if

      apar = 0.
      if (include_apar) then
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            call gyro_average(spread(vpa, 2, nmu) * g(:, :, ikxkyz), ikxkyz, g0)
            wgt = 2.0 * beta * spec(is)%z * spec(is)%dens * spec(is)%stm
            call integrate_vmu(g0, iz, tmp)
            apar(iky, ikx, iz, it) = apar(iky, ikx, iz, it) + tmp * wgt
         end do
         call sum_allreduce(apar)
         if (dist == 'g' .or. dist == 'h') then
            where (spread(kperp2(:, :, ia, :), 4, ntubes) > epsilon(0.0))
               apar = apar / spread(kperp2(:, :, ia, :), 4, ntubes)
            elsewhere
               apar = 0.0
            end where
         else if (dist == 'gbar') then
            apar = apar / spread(apar_denom, 4, ntubes)
         else if (dist == 'gstar') then
            write (*, *) 'APAR NOT SETUP FOR GSTAR YET. aborting.'
            call mp_abort('APAR NOT SETUP FOR GSTAR YET. aborting.')
         else
            if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
            call mp_abort('unknown dist option in get_fields. aborting')
         end if
         deallocate (g0)
      end if

   end subroutine get_fields

   subroutine get_fields_vmulo(g, phi, apar, dist, skip_fsa)

      use mp, only: mp_abort, proc0
      use job_manage, only: time_message
      use stella_layouts, only: vmu_lo, iv_idx
      use gyro_averages, only: gyro_average
      use run_parameters, only: fphi
      use physics_parameters, only: beta
      use physics_flags, only: include_apar
      use physics_flags, only: radial_variation
      use dist_fn_arrays, only: g_scratch
      use zgrid, only: nzgrid
      use vpamu_grids, only: integrate_species
      use vpamu_grids, only: vpa
      use species, only: spec

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar
      logical, optional, intent(in) :: skip_fsa
      character(*), intent(in) :: dist

      integer :: iv, ivmu
      logical :: skip_fsa_local

      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      if (debug) write (*, *) 'dist_fn::advance_stella::get_fields_vmulo'

      phi = 0.
      if (fphi > epsilon(0.0)) then

         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')

         ! gyroaverage the distribution function g at each phase space location
         call gyro_average(g, g_scratch)

         ! <g> requires modification if radial profile variation is included
         if (radial_variation) call add_radial_correction_int_species(g_scratch)

         ! integrate <g> over velocity space and sum over species
         !> store result in phi, which will be further modified below to account for polarization term
         if (debug) write (*, *) 'dist_fn::advance_stella::get_fields_vmulo::integrate_species_phi'
         call integrate_species(g_scratch, spec%z * spec%dens_psi0, phi)

         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')

         call get_phi(phi, dist, skip_fsa_local)

      end if

      apar = 0.
      if (include_apar) then

         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')

         ! if fphi > 0, then g_scratch = <g> already calculated above
         if (fphi < epsilon(0.0)) call gyro_average(g, g_scratch)

         ! for parallel Ampere's Law, need to calculate parallel current rather than density,
         ! so multiply <g> by vpa before integrating
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            ! get the vpa index
            iv = iv_idx(vmu_lo, ivmu)
            ! multiply by vpa
            g_scratch(:, :, :, :, ivmu) = g_scratch(:, :, :, :, ivmu) * vpa(iv)
         end do

         ! integrate vpa*<g> over velocity space and sum over species
         !> store result in apar, which will be further modified below to account for apar pre-factor
         if (debug) write (*, *) 'dist_fn::advance_stella::get_fields_vmulo::integrate_species_apar'
         call integrate_species(g_scratch, spec%z * spec%dens_psi0 * spec%stm_psi0 * beta, apar)

         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')

         ! divide the apar obtained above by the appropriate Apar pre-factor;
         ! this is just kperp2 if g = <f> is used or apar_denom = (kperp2 + ...)
         ! if gbar = g + <vpa*apar/c> * Ze/T * F_0 is used
         call get_apar(apar, dist)

      end if

   end subroutine get_fields_vmulo

   !> get_fields_ffs accepts as input the guiding centre distribution function g
   !> and calculates/returns the electronstatic potential phi for full_flux_surface simulations
   subroutine get_fields_ffs(g, phi, apar)

      use mp, only: mp_abort
      use physics_parameters, only: nine, tite
      use physics_flags, only: include_apar
      use stella_layouts, only: vmu_lo
      use run_parameters, only: fphi
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
      if (include_apar) then
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

   subroutine get_fields_by_spec(g, fld, skip_fsa)

      use mp, only: sum_allreduce
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use gyro_averages, only: gyro_average
      use run_parameters, only: fphi
      use stella_geometry, only: dl_over_b
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: integrate_vmu
      use kt_grids, only: nakx
      use kt_grids, only: zonal_mode
      use species, only: spec, nspec, has_electron_species
      use physics_flags, only: adiabatic_option_switch
      use physics_flags, only: adiabatic_option_fieldlineavg

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld
      logical, optional, intent(in) :: skip_fsa

      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      logical :: skip_fsa_local
      complex, dimension(nspec) :: tmp

      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      if (debug) write (*, *) 'dist_fn::advance_stella::get_fields_by_spec'

      ia = 1

      fld = 0.
      if (fphi > epsilon(0.0)) then
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            wgt = spec(is)%z * spec(is)%dens_psi0
            call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
            g0 = g0 * wgt
            call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
         end do
         call sum_allreduce(fld)

         fld = fld / gamtot_h

         if (.not. has_electron_species(spec) .and. (.not. skip_fsa_local) .and. &
             adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            if (zonal_mode(1)) then
               do ikx = 1, nakx
                  do it = 1, ntubes
                     do is = 1, nspec
                        tmp(is) = sum(dl_over_b(ia, :) * fld(1, ikx, :, it, is))
                        fld(1, ikx, :, it, is) = fld(1, ikx, :, it, is) + tmp(is) * gamtot3_h
                     end do
                  end do
               end do
            end if
         end if

         deallocate (g0)
      end if

   end subroutine get_fields_by_spec

   subroutine get_fields_by_spec_idx(isa, g, fld)

      ! apply phi_isa[ ] to all species indices contained in g
      ! ie get phi_isa[g_is1], phi_isa[g_is2], phi_isa[g_is3] ...

      use mp, only: sum_allreduce
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use gyro_averages, only: gyro_average
      use run_parameters, only: fphi
      use stella_geometry, only: dl_over_b, bmag
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: vperp2, nvpa, nmu
      use vpamu_grids, only: integrate_vmu
      use kt_grids, only: nakx
      use kt_grids, only: zonal_mode
      use species, only: spec, nspec, has_electron_species
      use physics_flags, only: adiabatic_option_switch
      use physics_flags, only: adiabatic_option_fieldlineavg
      use dist_fn_arrays, only: kperp2
      use spfunc, only: j0

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld
      integer, intent(in) :: isa

      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia, imu
      complex, dimension(nspec) :: tmp
      real :: wgt
      real :: arg

      ia = 1

      fld = 0.
      if (fphi > epsilon(0.0)) then
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            wgt = spec(isa)%z * spec(isa)%dens
            do imu = 1, nmu
               ! AVB: changed this for use of j0, check
               arg = spec(isa)%bess_fac * spec(isa)%smz_psi0 * sqrt(vperp2(ia, iz, imu) * kperp2(iky, ikx, ia, iz)) / bmag(ia, iz)
               g0(:, imu) = g(:, imu, ikxkyz) * j0(arg) ! AVB: gyroaverage
            end do
            g0 = g0 * wgt
            call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
         end do
         call sum_allreduce(fld)

         fld = fld / gamtot_h

         if (.not. has_electron_species(spec) .and. &
             adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            if (zonal_mode(1)) then
               do ikx = 1, nakx
                  do it = 1, ntubes
                     do is = 1, nspec
                        tmp(is) = sum(dl_over_b(ia, :) * fld(1, ikx, :, it, is))
                        fld(1, ikx, :, it, is) = fld(1, ikx, :, it, is) + tmp(is) * gamtot3_h
                     end do
                  end do
               end do
            end if
         end if

         deallocate (g0)
      end if

   end subroutine get_fields_by_spec_idx

   subroutine get_phi(phi, dist, skip_fsa)

      use mp, only: proc0, mp_abort, job
      use job_manage, only: time_message
      use physics_flags, only: radial_variation
      use run_parameters, only: ky_solve_radial, ky_solve_real
      use zgrid, only: nzgrid, ntubes
      use stella_geometry, only: dl_over_b
      use kt_grids, only: nakx, naky, zonal_mode
      use physics_flags, only: adiabatic_option_switch
      use physics_flags, only: adiabatic_option_fieldlineavg
      use species, only: spec, has_electron_species
      use multibox, only: mb_get_phi
      use fields_arrays, only: gamtot, gamtot3
      use file_utils, only: runtype_option_switch, runtype_multibox

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
      logical, optional, intent(in) :: skip_fsa
      real, dimension(:, :, :, :), allocatable :: gamtot_t
      integer :: ia, it, ikx
      complex :: tmp
      logical :: skip_fsa_local
      logical :: has_elec, adia_elec
      logical :: global_quasineutrality, center_cell
      logical :: multibox_mode

      character(*), intent(in) :: dist

      if (debug) write (*, *) 'dist_fn::advance_stella::get_phi'

      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ia = 1
      has_elec = has_electron_species(spec)
      adia_elec = .not. has_elec &
                  .and. adiabatic_option_switch == adiabatic_option_fieldlineavg

      global_quasineutrality = radial_variation .and. ky_solve_radial > 0
      multibox_mode = runtype_option_switch == runtype_multibox
      center_cell = multibox_mode .and. job == 1 .and. .not. ky_solve_real

      if (proc0) call time_message(.false., time_field_solve(:, 4), ' get_phi')
      if (dist == 'h') then
         phi = phi / gamtot_h
      else if (dist == 'g' .or. dist == 'gbar') then
         if (global_quasineutrality .and. (center_cell .or. .not. multibox_mode) .and. .not. ky_solve_real) then
            call get_phi_radial(phi)
         else if (global_quasineutrality .and. center_cell .and. ky_solve_real) then
            call mb_get_phi(phi, has_elec, adia_elec)
         else
            ! divide <g> by sum_s (\Gamma_0s-1) Zs^2*e*ns/Ts to get phi
            allocate (gamtot_t(naky, nakx, -nzgrid:nzgrid, ntubes))
            gamtot_t = spread(gamtot, 4, ntubes)
            where (gamtot_t < epsilon(0.0))
               phi = 0.0
            elsewhere
               phi = phi / gamtot_t
            end where
            deallocate (gamtot_t)
         end if
      else
         if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
         call mp_abort('unknown dist option in get_fields. aborting')
         return
      end if

      if (any(gamtot(1, 1, :) < epsilon(0.))) phi(1, 1, :, :) = 0.0
      if (proc0) call time_message(.false., time_field_solve(:, 4), ' get_phi')

      ! now handle adiabatic electrons if needed
      if (proc0) call time_message(.false., time_field_solve(:, 5), 'get_phi_adia_elec')
      if (adia_elec .and. zonal_mode(1) .and. .not. skip_fsa_local) then
         if (debug) write (*, *) 'dist_fn::advance_stella::adiabatic_electrons'

         if (dist == 'h') then
            do it = 1, ntubes
               do ikx = 1, nakx
                  tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                  phi(1, ikx, :, it) = phi(1, ikx, :, it) + tmp * gamtot3_h
               end do
            end do
         else if (dist == 'g' .or. dist == 'gbar') then
            if (global_quasineutrality .and. center_cell .and. ky_solve_real) then
               !this is already taken care of in mb_get_phi
            elseif (global_quasineutrality .and. (center_cell .or. .not. multibox_mode) &
                    .and. .not. ky_solve_real) then
               call add_adiabatic_response_radial(phi)
            else
               do ikx = 1, nakx
                  do it = 1, ntubes
                     tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                     phi(1, ikx, :, it) = phi(1, ikx, :, it) + tmp * gamtot3(ikx, :)
                  end do
               end do
            end if
         else
            if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
            call mp_abort('unknown dist option in get_fields. aborting')
         end if
      end if
      if (proc0) call time_message(.false., time_field_solve(:, 5), 'get_phi_adia_elec')

   end subroutine get_phi

   !> Non-perturbative approach to solving quasineutrality for radially
   !> global simulations
   subroutine get_phi_radial(phi)

#ifdef ISO_C_BINDING
      use mpi
      use mp, only: curr_focus, sharedsubprocs, scope
      use mp, only: split_n_tasks, sgproc0
      use zgrid, only: nztot
      use fields_arrays, only: phi_shared
      use mp_lu_decomposition, only: lu_matrix_multiply_local
#endif
      use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use physics_flags, only: adiabatic_option_switch
      use physics_flags, only: adiabatic_option_fieldlineavg
      use run_parameters, only: ky_solve_radial
      use zgrid, only: nzgrid, ntubes
      use species, only: spec, has_electron_species
      use kt_grids, only: nakx, naky, zonal_mode
      use linear_solve, only: lu_back_substitution
      use fields_arrays, only: gamtot, phi_solve

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
      integer :: it, iz, iky, zmi
      integer :: naky_r
      complex, dimension(:, :), allocatable :: g0k, g0x
      logical :: has_elec, adia_elec
#ifdef ISO_C_BINDING
      integer :: counter, c_lo, c_hi
      integer :: prior_focus, ierr
#endif

      allocate (g0k(1, nakx))
      allocate (g0x(1, nakx))

      has_elec = has_electron_species(spec)
      adia_elec = .not. has_elec &
                  .and. adiabatic_option_switch == adiabatic_option_fieldlineavg

      naky_r = min(naky, ky_solve_radial)
#ifdef ISO_C_BINDING
      prior_focus = curr_focus
      call scope(sharedsubprocs)

      call split_n_tasks(nztot * ntubes * naky_r, c_lo, c_hi)

      call scope(prior_focus)
      counter = 0
      if (sgproc0) phi_shared = phi
      call mpi_win_fence(0, phi_shared_window, ierr)
#endif
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do iky = 1, naky_r
#ifdef ISO_C_BINDING
               counter = counter + 1
               if ((counter >= c_lo) .and. (counter <= c_hi)) then
                  if (.not. (adia_elec .and. zonal_mode(iky))) then
                     zmi = 0
                     if (iky == 1) zmi = zm !zero mode may or may not be included in matrix
                     call lu_back_substitution(phi_solve(iky, iz)%zloc, &
                                               phi_solve(iky, iz)%idx, phi_shared(iky, (1 + zmi):, iz, it))
                     if (zmi > 0) phi(iky, zmi, iz, it) = 0.0
                  end if
               end if
#else
               if (.not. (adia_elec .and. zonal_mode(iky))) then
                  zmi = 0
                  if (iky == 1) zmi = zm !zero mode may or may not be included in matrix
                  call lu_back_substitution(phi_solve(iky, iz)%zloc, &
                                            phi_solve(iky, iz)%idx, phi(iky, (1 + zmi):, iz, it))
                  if (zmi > 0) phi(iky, zmi, iz, it) = 0.0
               end if
#endif
            end do
         end do
      end do
#ifdef ISO_C_BINDING
      call mpi_win_fence(0, phi_shared_window, ierr)
      phi = phi_shared
#endif

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do iky = naky_r + 1, naky
               phi(iky, :, iz, it) = phi(iky, :, iz, it) / gamtot(iky, :, iz)
            end do
         end do
      end do

      if (ky_solve_radial == 0 .and. any(gamtot(1, 1, :) < epsilon(0.))) &
         phi(1, 1, :, :) = 0.0

      deallocate (g0k, g0x)

   end subroutine get_phi_radial

   ! get_apar solves pre-factor * Apar = beta_ref * sum_s Z_s n_s vth_s int d3v vpa * J0 * pdf
   ! for apar, with pdf being either g or gbar (specified by dist input).
   ! the input apar is the RHS of the above equation and is overwritten by the true apar
   ! the pre-factor depends on whether g or gbar is used (kperp2 in former case, with additional
   ! term appearing in latter case)
   subroutine get_apar(apar, dist)

      use mp, only: proc0, mp_abort
      use zgrid, only: nzgrid, ntubes
      use dist_fn_arrays, only: kperp2

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: apar
      character(*), intent(in) :: dist

      integer :: ia

      ! this subroutine only considers flux tubes, so set ia = 1
      ia = 1
      if (dist == 'g') then
         where (spread(kperp2(:, :, ia, :), 4, ntubes) > epsilon(0.0))
            apar = apar / spread(kperp2(:, :, ia, :), 4, ntubes)
         elsewhere
            apar = 0.0
         end where
      else if (dist == 'gbar') then
         apar = apar / spread(apar_denom, 4, ntubes)
      else
         if (proc0) write (*, *) 'unknown dist option in get_apar. aborting'
         call mp_abort('unkown dist option in get_apar. aborting')
      end if

   end subroutine get_apar

   subroutine advance_apar(g, dist, apar)

      use mp, only: mp_abort, sum_allreduce
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use physics_flags, only: include_apar
      use physics_parameters, only: beta
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: nvpa, nmu, vpa
      use vpamu_grids, only: integrate_vmu
      use gyro_averages, only: gyro_average

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      character(*), intent(in) :: dist
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: apar

      integer :: ikxkyz, iky, ikx, iz, it, is
      real :: wgt
      complex :: tmp
      complex, dimension(:, :), allocatable :: scratch

      apar = 0.
      if (include_apar) then
         allocate (scratch(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            call gyro_average(spread(vpa, 2, nmu) * g(:, :, ikxkyz), ikxkyz, scratch)
            wgt = beta * spec(is)%z * spec(is)%dens_psi0 * spec(is)%stm_psi0
            call integrate_vmu(scratch, iz, tmp)
            apar(iky, ikx, iz, it) = apar(iky, ikx, iz, it) + tmp * wgt
         end do
         ! apar for different species may be spread over processors at this point, so
         ! broadcast to all procs and sum over species
         call sum_allreduce(apar)
         ! divide by the appropriate apar pre-factor to get apar
         call get_apar(apar, dist)
         deallocate (scratch)
      end if

   end subroutine advance_apar

   !> Add the adiabatic eletron contribution for globally radial simulations.
   !> This actually entails solving for the whole ky = 0 slice of phi at once (not really adding!)
   subroutine add_adiabatic_response_radial(phi)

#ifdef ISO_C_BINDING
      use mpi
      use mp, only: sgproc0, comm_sgroup
      use fields_arrays, only: qn_zf_window
      use mp_lu_decomposition, only: lu_matrix_multiply_local
#else
      use linear_solve, only: lu_back_substitution
#endif
      use zgrid, only: nzgrid, ntubes
      use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use stella_geometry, only: dl_over_b, d_dl_over_b_drho
      use kt_grids, only: nakx, boundary_size, rho_d_clamped
      use fields_arrays, only: phizf_solve, phi_ext
      use fields_arrays, only: phi_proj, phi_proj_stage, theta
      use fields_arrays, only: exclude_boundary_regions_qn, exp_fac_qn, tcorr_source_qn

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
      integer :: ia, it, iz, ikx
      integer :: inmat
      complex, dimension(:, :), allocatable :: g0k, g1k, g0x
#ifdef ISO_C_BINDING
      integer :: ierr
#endif

      allocate (g0k(1, nakx))
      allocate (g1k(1, nakx))
      allocate (g0x(1, nakx))

      ia = 1

      do it = 1, ntubes
         ! calculate <<g>_psi>_T
         g1k = 0.0
         do iz = -nzgrid, nzgrid - 1
            g0k(1, :) = phi(1, :, iz, it)
            call transform_kx2x_unpadded(g0k, g0x)
            g0x(1, :) = (dl_over_b(ia, iz) + d_dl_over_b_drho(ia, iz) * rho_d_clamped) * g0x(1, :)
            if (exclude_boundary_regions_qn) then
               g0x(1, :) = sum(g0x(1, (boundary_size + 1):(nakx - boundary_size))) &
                           / (nakx - 2 * boundary_size)
               g0x(1, 1:boundary_size) = 0.0
               g0x(1, (nakx - boundary_size + 1):) = 0.0
            else
               g0x(1, :) = sum(g0x(1, :)) / nakx
            end if

            call transform_x2kx_unpadded(g0x, g0k)

            g1k = g1k + g0k
         end do

         phi_proj_stage(:, 1, it) = g1k(1, :)
         if (tcorr_source_qn < epsilon(0.0)) then
            do iz = -nzgrid, nzgrid - 1
               phi(1, :, iz, it) = phi(1, :, iz, it) - g1k(1, :)
            end do
         else
            do iz = -nzgrid, nzgrid - 1
               phi(1, :, iz, it) = phi(1, :, iz, it) &
                                   - (1.-exp_fac_qn) * g1k(1, :) - exp_fac_qn * phi_proj(:, 1, it)
            end do
         end if

#ifdef ISO_C_BINDING
         if (sgproc0) then
#endif
            do iz = -nzgrid, nzgrid - 1
               do ikx = 1, nakx
                  inmat = ikx + nakx * (iz + nzgrid)
                  phi_ext(inmat) = phi(1, ikx, iz, it)
               end do
            end do
#ifdef ISO_C_BINDING
         end if
         call mpi_win_fence(0, qn_zf_window, ierr)
#endif

#ifdef ISO_C_BINDING
         call lu_matrix_multiply_local(comm_sgroup, qn_zf_window, phizf_solve%zloc, phi_ext)
         call mpi_win_fence(0, qn_zf_window, ierr)
#else
         call lu_back_substitution(phizf_solve%zloc, phizf_solve%idx, phi_ext)
#endif

         do iz = -nzgrid, nzgrid - 1
            do ikx = 1, nakx
               inmat = ikx + nakx * (iz + nzgrid)
               phi(1, ikx, iz, it) = phi_ext(inmat)
            end do
         end do

         !enforce periodicity
         phi(1, :, nzgrid, it) = phi(1, :, -nzgrid, it)

         ! calculate Theta.phi
         g1k = 0.0
         do iz = -nzgrid, nzgrid - 1
            do ikx = 1, nakx
               g0k(1, ikx) = sum(theta(ikx, :, iz) * phi(1, :, iz, it))
            end do

            call transform_kx2x_unpadded(g0k, g0x)

            g0x(1, :) = (dl_over_b(ia, iz) + d_dl_over_b_drho(ia, iz) * rho_d_clamped) * g0x(1, :)
            if (exclude_boundary_regions_qn) then
               g0x(1, :) = sum(g0x(1, (boundary_size + 1):(nakx - boundary_size))) &
                           / (nakx - 2 * boundary_size)
               g0x(1, 1:boundary_size) = 0.0
               g0x(1, (nakx - boundary_size + 1):) = 0.0
            else
               g0x(1, :) = sum(g0x(1, :)) / nakx
            end if

            call transform_x2kx_unpadded(g0x, g0k)
            g1k = g1k + g0k
         end do

         phi_proj_stage(:, 1, it) = phi_proj_stage(:, 1, it) - g1k(1, :)
      end do
      deallocate (g0k, g1k, g0x)

   end subroutine add_adiabatic_response_radial

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

   !> Add radial variation of the Jacobian and gyroaveraing in the velocity integration of
   !> <g>, needed for radially global simulations
   subroutine add_radial_correction_int_species(g_in)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: imu_idx, is_idx
      use gyro_averages, only: aj0x, aj1x
      use stella_geometry, only: dBdrho, bmag
      use dist_fn_arrays, only: kperp2, dkperp2dr
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: vperp2
      use kt_grids, only: nakx, naky, multiply_by_rho
      use run_parameters, only: ky_solve_radial
      use species, only: spec

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(inout) :: g_in

      integer :: ivmu, iz, it, ia, imu, is, iky
      complex, dimension(:, :), allocatable :: g0k

      if (ky_solve_radial <= 0) return

      allocate (g0k(naky, nakx))

      ia = 1

      ! loop over super-index ivmu, which include vpa, mu and spec
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! is = species index
         is = is_idx(vmu_lo, ivmu)
         ! imu = mu index
         imu = imu_idx(vmu_lo, ivmu)

         ! loop over flux tubes in flux tube train
         do it = 1, ntubes
            ! loop over zed location within flux tube
            do iz = -nzgrid, nzgrid
               g0k = 0.0
               do iky = 1, min(ky_solve_radial, naky)
                  g0k(iky, :) = g_in(iky, :, iz, it, ivmu) &
                                * (-0.5 * aj1x(iky, :, iz, ivmu) / aj0x(iky, :, iz, ivmu) * (spec(is)%smz)**2 &
                                   * (kperp2(iky, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                   * (dkperp2dr(iky, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                   + dBdrho(iz) / bmag(ia, iz))

               end do
               !g0k(1,1) = 0.
               call multiply_by_rho(g0k)
               g_in(:, :, iz, it, ivmu) = g_in(:, :, iz, it, ivmu) + g0k
            end do
         end do
      end do

      deallocate (g0k)

   end subroutine add_radial_correction_int_species

   !> the following routine gets the correction in phi both from gyroaveraging and quasineutrality
   subroutine get_radial_correction(g, phi0, dist)

      use mp, only: proc0, mp_abort, sum_allreduce
      use stella_layouts, only: vmu_lo
      use gyro_averages, only: gyro_average, gyro_average_j1
      use gyro_averages, only: aj0x, aj1x
      use run_parameters, only: fphi, ky_solve_radial
      use stella_geometry, only: dl_over_b, d_dl_over_b_drho, bmag, dBdrho
      use stella_layouts, only: imu_idx, is_idx
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: integrate_species, vperp2
      use kt_grids, only: nakx, nx, naky, rho_d_clamped
      use kt_grids, only: zonal_mode, multiply_by_rho
      use species, only: spec, has_electron_species
      use fields_arrays, only: phi_corr_QN, phi_corr_GA
      use fields_arrays, only: gamtot, dgamtotdr
      use fields_arrays, only: gamtot3
      use dist_fn_arrays, only: kperp2, dkperp2dr
      use physics_flags, only: adiabatic_option_switch
      use physics_flags, only: adiabatic_option_fieldlineavg
      use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi0
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      character(*), intent(in) :: dist

      integer :: ikx, iky, ivmu, iz, it, ia, is, imu
      complex :: tmp
      real, dimension(:, :, :, :), allocatable :: gamtot_t
      complex, dimension(:, :, :, :), allocatable :: phi1
      complex, dimension(:, :, :), allocatable :: gyro_g
      complex, dimension(:, :), allocatable :: g0k, g1k, g1x

      ia = 1

      if (fphi > epsilon(0.0)) then
         allocate (gyro_g(naky, nakx, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (g0k(naky, nakx))
         allocate (phi1(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi1 = 0.
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                  is = is_idx(vmu_lo, ivmu)
                  imu = imu_idx(vmu_lo, ivmu)

                  g0k = g(:, :, iz, it, ivmu) &
                        * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) &
                           * (spec(is)%smz)**2 &
                           * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                           + dBdrho(iz) / bmag(ia, iz))

                  call gyro_average(g0k, iz, ivmu, gyro_g(:, :, ivmu))
               end do
               call integrate_species(gyro_g, iz, spec%z * spec%dens_psi0, phi1(:, :, iz, it), reduce_in=.false.)
            end do
         end do
         call sum_allreduce(phi1)

         !apply radial operator Xhat
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0k = phi1(:, :, iz, it) - dgamtotdr(:, :, iz) * phi0(:, :, iz, it)
               call multiply_by_rho(g0k)
               phi1(:, :, iz, it) = g0k
            end do
         end do

         if (dist == 'gbar' .or. dist == 'g') then
            allocate (gamtot_t(naky, nakx, -nzgrid:nzgrid, ntubes))
            gamtot_t = spread(gamtot, 4, ntubes)
            where (gamtot_t < epsilon(0.0))
               phi1 = 0.0
            elsewhere
               phi1 = phi1 / gamtot_t
            end where
            deallocate (gamtot_t)
         else if (dist == 'h') then
            if (proc0) write (*, *) 'dist option "h" not implemented in radial_correction. aborting'
            call mp_abort('dist option "h" in radial_correction. aborting')
         else
            if (proc0) write (*, *) 'unknown dist option in radial_correction. aborting'
            call mp_abort('unknown dist option in radial_correction. aborting')
         end if

         if (.not. has_electron_species(spec) .and. &
             adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            if (zonal_mode(1)) then
               if (dist == 'g' .or. dist == 'gbar') then
                  allocate (g1k(1, nakx))
                  allocate (g1x(1, nakx))
                  do it = 1, ntubes
                     do ikx = 1, nakx
                        g1k(1, ikx) = sum(phi0(1, ikx, :, it) &
                                          * (efacp * dl_over_b(ia, :) + efac * d_dl_over_b_drho(ia, :)))
                     end do
                     call transform_kx2x_unpadded(g1k, g1x)
                     g1x(1, :) = rho_d_clamped * g1x(1, :)
                     call transform_x2kx_unpadded(g1x, g1k)

                     do ikx = 1, nakx
                        phi1(1, ikx, :, it) = phi1(1, ikx, :, it) + g1k(1, ikx) / gamtot(1, ikx, :)
                        tmp = sum(dl_over_b(ia, :) * phi1(1, ikx, :, it))
                        phi1(1, ikx, :, it) = phi1(1, ikx, :, it) + gamtot3(ikx, :) * tmp
                     end do
                  end do
                  deallocate (g1k, g1x)
               else
                  if (proc0) write (*, *) 'unknown dist option in radial_correction. aborting'
                  call mp_abort('unknown dist option in radial_correction. aborting')
               end if
            end if
         end if

         !> collect quasineutrality corrections in wavenumber space
         phi_corr_QN = phi1

         !> zero out the ones we have already solved for using the full method
         do iky = 1, min(ky_solve_radial, naky)
            phi_corr_QN(iky, :, :, :) = 0.0
         end do

         deallocate (phi1)

         !> collect gyroaveraging corrections in wavenumber space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call gyro_average_j1(phi0(:, :, iz, it), iz, ivmu, g0k)
                  g0k = -g0k * (spec(is)%smz)**2 &
                        * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                        * 0.5 * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz))

                  call multiply_by_rho(g0k)
                  phi_corr_GA(:, :, iz, it, ivmu) = g0k
               end do
            end do
         end do

         deallocate (g0k)
         deallocate (gyro_g)

      end if

   end subroutine get_radial_correction

   !> rescale fields, including the distribution function
   subroutine rescale_fields(target_amplitude)

      use mp, only: scope, subprocs, crossdomprocs, sum_allreduce
      use fields_arrays, only: phi, apar
      use dist_fn_arrays, only: gnew, gvmu
      use volume_averages, only: volume_average
      use job_manage, only: njobs
      use file_utils, only: runtype_option_switch, runtype_multibox

      implicit none

      real, intent(in) :: target_amplitude
      real :: phi2, rescale

      call volume_average(phi, phi2)

      if (runtype_option_switch == runtype_multibox) then
         call scope(crossdomprocs)
         call sum_allreduce(phi2)
         call scope(subprocs)
         phi2 = phi2 / njobs
      end if

      rescale = target_amplitude / sqrt(phi2)

      phi = rescale * phi
      apar = rescale * apar
      gnew = rescale * gnew
      gvmu = rescale * gvmu

   end subroutine rescale_fields

   !> compute d<chi>/dy in (ky,kx,z,tube) space
   subroutine get_dchidy_4d(phi, apar, bpar, dchidy)

      use constants, only: zi
      use gyro_averages, only: gyro_average
      use gyro_averages, only: gyro_average_j1
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx, imu_idx
      use physics_flags, only: include_apar
      use physics_flags, only: include_bpar
      use run_parameters, only: fphi
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: vpa, mu
      use kt_grids, only: aky, naky, nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dchidy

      integer :: ivmu, iv, is, iky, imu
      complex, dimension(:, :, :, :), allocatable :: field, gyro_tmp

      allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (gyro_tmp(naky, nakx, -nzgrid:nzgrid, ntubes))

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         is = is_idx(vmu_lo, ivmu)
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         ! intermediate calculation to get factor involving phi contribution
         field = fphi * phi
         ! add apar contribution if including it
         if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
         ! take spectral y-derivative
         do iky = 1, naky
            field(iky, :, :, :) = zi * aky(iky) * field(iky, :, :, :)
         end do
         call gyro_average(field, ivmu, dchidy(:, :, :, :, ivmu))
         if (include_bpar) then
            field = 4 * mu(imu) * (spec(is)%tz) * bpar
            do iky = 1, naky
               field(iky, :, :, :) = zi * aky(iky) * field(iky, :, :, :)
            end do
            call gyro_average_j1(field, ivmu, gyro_tmp)
            !> include bpar contribution
            dchidy(:, :, :, :, ivmu) = dchidy(:, :, :, :, ivmu) + gyro_tmp
         end if
      end do

      deallocate (field)
      deallocate (gyro_tmp)

   end subroutine get_dchidy_4d

   !> compute d<chi>/dy in (ky,kx) space
   subroutine get_dchidy_2d(iz, ivmu, phi, apar, bpar, dchidy)

      use constants, only: zi
      use gyro_averages, only: gyro_average
      use gyro_averages, only: gyro_average_j1
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx, imu_idx
      use physics_flags, only: include_apar
      use physics_flags, only: include_bpar
      use run_parameters, only: fphi
      use species, only: spec
      use vpamu_grids, only: vpa, mu
      use kt_grids, only: nakx, aky, naky

      implicit none

      integer, intent(in) :: ivmu, iz
      complex, dimension(:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :), intent(out) :: dchidy

      integer :: iv, is, imu
      complex, dimension(:, :), allocatable :: field, gyro_tmp

      allocate (field(naky, nakx))
      allocate (gyro_tmp(naky, nakx))

      is = is_idx(vmu_lo, ivmu)
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      field = fphi * phi
      if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
      field = zi * spread(aky, 2, nakx) * field
      call gyro_average(field, iz, ivmu, dchidy)
      if (include_bpar) then
         field = 4 * mu(imu) * (spec(is)%tz) * bpar
         field = zi * spread(aky, 2, nakx) * field
         call gyro_average_j1(field, iz, ivmu, gyro_tmp)
         !> include bpar contribution
         dchidy = dchidy + gyro_tmp
      end if
      deallocate (field)
      deallocate (gyro_tmp)

   end subroutine get_dchidy_2d

   !> compute d<chi>/dx in (ky,kx) space
   subroutine get_dchidx(iz, ivmu, phi, apar, bpar, dchidx)

      use constants, only: zi
      use gyro_averages, only: gyro_average
      use gyro_averages, only: gyro_average_j1
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx, imu_idx
      use physics_flags, only: include_apar
      use physics_flags, only: include_bpar
      use run_parameters, only: fphi
      use species, only: spec
      use vpamu_grids, only: vpa, mu
      use kt_grids, only: akx, naky, nakx

      implicit none

      integer, intent(in) :: ivmu, iz
      complex, dimension(:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :), intent(out) :: dchidx

      integer :: iv, is, imu
      complex, dimension(:, :), allocatable :: field, gyro_tmp

      allocate (field(naky, nakx))
      allocate (gyro_tmp(naky, nakx))

      is = is_idx(vmu_lo, ivmu)
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      field = fphi * phi
      if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
      field = zi * spread(akx, 1, naky) * field
      call gyro_average(field, iz, ivmu, dchidx)
      if (include_bpar) then
         field = 4 * mu(imu) * (spec(is)%tz) * bpar
         field = zi * spread(akx, 1, naky) * field
         call gyro_average_j1(field, iz, ivmu, gyro_tmp)
         !> include bpar contribution
         dchidx = dchidx + gyro_tmp
      end if
      deallocate (field)
      deallocate (gyro_tmp)

   end subroutine get_dchidx

   subroutine finish_fields

      use fields_arrays, only: phi, phi_old
      use fields_arrays, only: apar, apar_old
      use fields_arrays, only: bpar, bpar_old
      use fields_arrays, only: phi_corr_QN, phi_corr_GA
      use fields_arrays, only: apar_corr_QN, apar_corr_GA
      use fields_arrays, only: gamtot, dgamtotdr, gamtot3
      use fields_arrays, only: gamtot13, gamtot33, gamtot31
      use fields_arrays, only: c_mat, theta
#ifdef ISO_C_BINDING
      use fields_arrays, only: qn_window
      use mpi
#else
      use fields_arrays, only: phi_solve
#endif
      implicit none

#ifdef ISO_C_BINDING
      integer ierr
#endif

      if (allocated(phi)) deallocate (phi)
      if (allocated(phi_old)) deallocate (phi_old)
      if (allocated(apar)) deallocate (apar)
      if (allocated(apar_old)) deallocate (apar_old)
      if (allocated(bpar)) deallocate (bpar)
      if (allocated(bpar_old)) deallocate (bpar_old)
      if (allocated(phi_corr_QN)) deallocate (phi_corr_QN)
      if (allocated(phi_corr_GA)) deallocate (phi_corr_GA)
      if (allocated(apar_corr_QN)) deallocate (apar_corr_QN)
      if (allocated(apar_corr_GA)) deallocate (apar_corr_GA)
      if (allocated(gamtot)) deallocate (gamtot)
      if (allocated(gamtot3)) deallocate (gamtot3)
      if (allocated(dgamtotdr)) deallocate (dgamtotdr)
      if (allocated(apar_denom)) deallocate (apar_denom)
      if (allocated(gamtot33)) deallocate (gamtot33)
      if (allocated(gamtot13)) deallocate (gamtot13)
      if (allocated(gamtot31)) deallocate (gamtot31)
#ifdef ISO_C_BINDING
      if (phi_shared_window /= MPI_WIN_NULL) call mpi_win_free(phi_shared_window, ierr)
      if (qn_window /= MPI_WIN_NULL) then
         call mpi_win_free(qn_window, ierr)
      end if
#else
      if (allocated(phi_solve)) deallocate (phi_solve)
#endif
      if (allocated(c_mat)) deallocate (c_mat)
      if (allocated(theta)) deallocate (theta)

      !> arrays only allocated/used if simulating a full flux surface
      if (allocated(gam0_ffs)) deallocate (gam0_ffs)
      if (allocated(lu_gam0_ffs)) deallocate (lu_gam0_ffs)
      if (allocated(adiabatic_response_factor)) deallocate (adiabatic_response_factor)

      fields_initialized = .false.

   end subroutine finish_fields

end module fields
