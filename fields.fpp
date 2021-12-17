module fields

   use common_types, only: eigen_type
#if defined MPI && defined ISO_C_BINDING
   use mpi
#endif

   implicit none

   public :: init_fields, finish_fields
   public :: advance_fields, get_fields
   public :: get_fields_vmulo_single
   public :: get_radial_correction
   public :: enforce_reality_field
   public :: get_fields_by_spec, get_fields_by_spec_idx
   public :: gamtot_h, gamtot3_h, gamtot3, dgamtot3dr
   public :: time_field_solve
   public :: fields_updated
   public :: get_dchidy, get_dchidx
   public :: efac, efacp
   public :: get_gyroaverage_chi
   public :: get_chi

   private

   real, dimension(:, :, :), allocatable ::  apar_denom
   real, dimension(:, :), allocatable :: gamtot3
   real :: gamtot_h, gamtot3_h, efac, efacp
   complex, dimension(:, :), allocatable :: b_mat

   real, dimension(:, :), allocatable :: dgamtot3dr

   complex, dimension(:, :), allocatable :: save1, save2

   logical :: fields_updated = .false.
   logical :: fields_initialized = .false.
#if defined MPI && defined ISO_C_BINDING
   logical :: qn_window_initialized = .false.
   integer :: phi_shared_window = MPI_WIN_NULL
#endif
   logical :: debug = .false.

   integer :: zm

   real, dimension(2, 5) :: time_field_solve

   interface get_dchidy
      module procedure get_dchidy_4d
      module procedure get_dchidy_2d
   end interface

   interface get_dchidx
      module procedure get_dchidx_4d
      module procedure get_dchidx_2d
   end interface

   interface get_gyroaverage_chi
      module procedure get_gyroaverage_chi_4d
      module procedure get_gyroaverage_chi_2d
   end interface

contains

   subroutine init_fields

      use mp, only: sum_allreduce, job, proc0
#if defined MPI && defined ISO_C_BINDING
      use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
      use fields_arrays, only: qn_window, phi_shared
      use mp, only: sgproc0, curr_focus, mp_comm, sharedsubprocs
      use mp, only: scope, real_size, nbytes_real, iproc, nproc
      use mpi
#else
#endif
      use linear_solve, only: lu_decomposition
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, onlY: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use dist_fn_arrays, only: kperp2, dkperp2dr
      use gyro_averages, only: aj0v, aj1v
      use run_parameters, only: fphi, fapar, fbpar
      use run_parameters, only: ky_solve_radial, ky_solve_real
      use physics_parameters, only: tite, nine, beta
      use physics_flags, only: radial_variation
      use species, only: spec, has_electron_species, ion_species
      use stella_geometry, only: dl_over_b, d_dl_over_b_drho, dBdrho, bmag
      use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
      use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use zgrid, only: nzgrid, ntubes, nztot
      use vpamu_grids, only: nvpa, nmu, mu
      use vpamu_grids, only: vpa, vperp2
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: integrate_vmu
      use species, only: spec, nspec
      use kt_grids, only: naky, nakx, akx
      use kt_grids, only: zonal_mode, rho_d_clamped
      use dist_fn, only: adiabatic_option_switch
      use dist_fn, only: adiabatic_option_fieldlineavg
      use multibox, only: init_mb_get_phi
      use fields_arrays, only: gamtot, dgamtotdr, gamtot3, dgamtot3dr
      use fields_arrays, only: phi_solve, c_mat, theta
      use fields_arrays, only: gamtot13, gamtot31, gamtot33, apar_denom
      use file_utils, only: runtype_option_switch, runtype_multibox

      implicit none

      integer :: ikxkyz, iz, it, ikx, iky, is, ia, zmi, naky_r
      real :: tmp, tmp2, wgt, dum
      real, dimension(:, :), allocatable :: g0
      real, dimension(:), allocatable :: g1
      logical :: has_elec, adia_elec
#if defined MPI && ISO_C_BINDING
      integer :: prior_focus, ierr
      integer :: counter, c_lo, c_hi, c_max, c_div, c_mod
      integer :: disp_unit = 1
      integer*8 :: cur_pos
      integer(kind=MPI_ADDRESS_KIND) :: win_size
      complex, dimension(:), pointer :: phi_shared_temp
      type(c_ptr) :: cptr
#endif

      complex, dimension(:, :), allocatable :: g0k, g0x

      debug = debug .and. proc0

      ia = 1
      zm = 0

      ! do not see why this is before fields_initialized check below
      call allocate_arrays

      if (fields_initialized) return
      fields_initialized = .true.

      ! could move these array allocations to allocate_arrays to clean up code
      if (.not. allocated(gamtot)) allocate (gamtot(naky, nakx, -nzgrid:nzgrid)); gamtot = 0.
      if (.not. allocated(gamtot13)) allocate (gamtot13(naky, nakx, -nzgrid:nzgrid)); gamtot13 = 0.
      if (.not. allocated(gamtot31)) allocate (gamtot31(naky, nakx, -nzgrid:nzgrid)); gamtot31 = 0.
      if (.not. allocated(gamtot33)) allocate (gamtot33(naky, nakx, -nzgrid:nzgrid)); gamtot33 = 0.
      ! Bob: gamone needed for calculation of bpar using distribution function g,
      ! but not currently implemented.
      !if (.not.allocated(gamone)) allocate (gamone(naky,nakx,-nzgrid:nzgrid)) ; gamone = 0.
      if (.not. allocated(gamtot3)) then
         if (.not. has_electron_species(spec) &
             .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            allocate (gamtot3(nakx, -nzgrid:nzgrid)); gamtot3 = 0.
         else
            allocate (gamtot3(1, 1)); gamtot3 = 0.
         end if
      end if
      if (.not. allocated(apar_denom)) then
         if (fapar > epsilon(0.0)) then
            allocate (apar_denom(naky, nakx, -nzgrid:nzgrid)); apar_denom = 0.
         else
            allocate (apar_denom(1, 1, 1)); apar_denom = 0.
         end if
      end if

      if (radial_variation) then
         if (.not. allocated(dgamtotdr)) allocate (dgamtotdr(naky, nakx, -nzgrid:nzgrid)); dgamtotdr = 0.
         if (.not. allocated(dgamtot3dr)) then
            if (.not. has_electron_species(spec) &
                .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
               allocate (dgamtot3dr(nakx, -nzgrid:nzgrid)); dgamtot3dr = 0.
               allocate (save1(nakx, ntubes)); save1 = 0.
               allocate (save2(nakx, ntubes)); save2 = 0.
            else
               allocate (dgamtot3dr(1, 1)); dgamtot3dr = 0.
            end if
         end if
      else
         if (.not. allocated(dgamtotdr)) allocate (dgamtotdr(1, 1, 1)); dgamtotdr = 0.
         if (.not. allocated(dgamtot3dr)) allocate (dgamtot3dr(1, 1)); dgamtot3dr = 0.
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
            g0 = spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa) &
                 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
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
                     if (radial_variation) then
                        tmp2 = (spec(ion_species)%tprim - spec(ion_species)%fprim) / efac &
                               + sum(d_dl_over_b_drho(ia, :) / gamtot(1, ikx, :)) &
                               - sum(dl_over_b(ia, :) * dgamtotdr(1, ikx, :) &
                                     / gamtot(1, ikx, :)**2)
                        dgamtot3dr(ikx, :) = gamtot3(ikx, :) &
                                             * (-dgamtotdr(1, ikx, :) / gamtot(1, ikx, :) + tmp2 / tmp)
                     end if
                  end do
                  if (akx(1) < epsilon(0.)) then
                     gamtot3(1, :) = 0.0
                     dgamtot3dr(1, :) = 0.0
                     zm = 1
                  end if
               end if
            end if
         end if

         if (radial_variation .and. ky_solve_radial > 0) then

            naky_r = min(naky, ky_solve_radial)

            has_elec = has_electron_species(spec)
            adia_elec = .not. has_elec .and. zonal_mode(1) &
                        .and. adiabatic_option_switch == adiabatic_option_fieldlineavg

            if (runtype_option_switch == runtype_multibox .and. job == 1 .and. ky_solve_real) then
               call init_mb_get_phi(has_elec, adia_elec, efac, efacp)
            elseif (runtype_option_switch /= runtype_multibox .or. &
                    (job == 1 .and. .not. ky_solve_real)) then
               allocate (g0k(1, nakx))
               allocate (g0x(1, nakx))

               if (.not. allocated(phi_solve)) allocate (phi_solve(naky_r, -nzgrid:nzgrid))
#if defined MPI && ISO_C_BINDING
               prior_focus = curr_focus
               call scope(sharedsubprocs)
               !the following is to parallelize the calculation of QN for radial variation sims
               if (debug) write (*, *) 'fields::init_fields::phi_shared_init'
               if (phi_shared_window == MPI_WIN_NULL) then
                  win_size = 0
                  if (sgproc0) then
                     win_size = int(naky * nakx * nztot * ntubes, MPI_ADDRESS_KIND) * 2 * real_size !complex size
                  end if

                  call mpi_win_allocate_shared(win_size, disp_unit, MPI_INFO_NULL, &
                                               mp_comm, cptr, phi_shared_window, ierr)

                  if (.not. sgproc0) then
                     !make sure all the procs have the right memory address
                     call mpi_win_shared_query(phi_shared_window, 0, win_size, disp_unit, cptr, ierr)
                  end if
                  call mpi_win_fence(0, phi_shared_window, ierr)

                  if (.not. associated(phi_shared)) then
                     ! associate array with lower bounds of 1
                     call c_f_pointer(cptr, phi_shared_temp, (/naky * nakx * nztot * ntubes/))
                     ! now get the correct bounds
                     phi_shared(1:naky, 1:nakx, -nzgrid:nzgrid, 1:ntubes) => phi_shared_temp
                  end if
                  call mpi_win_fence(0, phi_shared_window, ierr)
               end if

               if (debug) write (*, *) 'fields::init_fields::qn_window_init'
               if ((.not. qn_window_initialized) .or. (qn_window == MPI_WIN_NULL)) then
                  win_size = 0
                  if (sgproc0) then
                     win_size = int(nakx * nztot * naky_r, MPI_ADDRESS_KIND) * 4_MPI_ADDRESS_KIND &
                                + int(nakx**2 * nztot * naky_r, MPI_ADDRESS_KIND) * 2 * real_size !complex size
                  end if

                  call mpi_win_allocate_shared(win_size, disp_unit, MPI_INFO_NULL, &
                                               mp_comm, cptr, qn_window, ierr)

                  if (.not. sgproc0) then
                     !make sure all the procs have the right memory address
                     call mpi_win_shared_query(qn_window, 0, win_size, disp_unit, cptr, ierr)
                  end if
                  call mpi_win_fence(0, qn_window, ierr)
                  cur_pos = transfer(cptr, cur_pos)

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

                  qn_window_initialized = .true.
               end if
               c_max = nztot * naky_r
               c_div = c_max / nproc
               c_mod = mod(c_max, nproc)

               c_lo = iproc * c_div + 1 + min(iproc, c_mod)
               c_hi = c_lo + c_div - 1
               if (iproc < c_mod) c_hi = c_hi + 1

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
#if defined MPI && ISO_C_BINDING
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
#if defined MPI && ISO_C_BINDING
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
         end if
         deallocate (g0)
      end if

      if (fapar > epsilon(0.)) then
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
            ! apar_denom = kperp^2 + 2 beta * sum(Z^2  * n / m * integrate_vmu (vpa*vpa*exp(-v^2) J0^2) )
            g0 = spread(maxwell_vpa(:, is) * vpa**2, 2, nmu) * maxwell_fac(is) &
                 * spread(maxwell_mu(ia, iz, :, is) * aj0v(:, ikxkyz) * aj0v(:, ikxkyz), 1, nvpa)
            wgt = 2.0 * beta * spec(is)%z * spec(is)%z * spec(is)%dens / spec(is)%mass
            call integrate_vmu(g0, iz, tmp)
            apar_denom(iky, ikx, iz) = apar_denom(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(apar_denom)
         apar_denom = apar_denom + kperp2(:, :, ia, :)

         deallocate (g0)
      end if

      ! Bob: Calculate gamtot13, gamtot31, gamtot33. These are needed provided
      ! fbpar != 0
      if (fbpar > epsilon(0.0)) then
         ! gamtot13
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            ! gamtot13 does not depend on flux tube index,
            ! so only compute for one flux tube index
            ! gamtot13 = -4 * sum_s (Z*n* integrate_vmu(mu*exp(-v^2) * J0 *J1/gamma))
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            g0 = spread((mu(:) * aj0v(:, ikxkyz) * aj1v(:, ikxkyz)), 1, nvpa) &
                 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
            wgt = -4 * spec(is)%z * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            gamtot13(iky, ikx, iz) = gamtot13(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(gamtot13)
         g0 = 0

         ! gamtot31 = 2 * beta * sum_s (Z*n* integrate_vmu(mu*exp(-v^2) * J0 *J1/gamma))
         !          = -gamtot13/2 * beta
         gamtot31 = -gamtot13 / 2 * beta

         ! gamtot33
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            ! gamtot33 does not depend on flux tube index,
            ! so only compute for one flux tube index
            ! gamtot33 = 1 + 8 * beta * sum_s (n*T* integrate_vmu(mu*mu*exp(-v^2) *(J1/gamma)*(J1/gamma)))
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            g0 = spread((mu(:) * mu(:) * aj1v(:, ikxkyz) * aj1v(:, ikxkyz)), 1, nvpa) &
                 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
            wgt = 8 * spec(is)%temp * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            gamtot33(iky, ikx, iz) = gamtot33(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(gamtot33)

         gamtot33 = 1.0 + beta * gamtot33
         deallocate (g0)
         ! gamtot_h = sum(spec%z*spec%z*spec%dens/spec%temp)
         !
         ! if (radial_variation) then
         !   allocate (g1(nmu))
         !   do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         !     it = it_idx(kxkyz_lo,ikxkyz)
         !     ! gamtot does not depend on flux tube index,
         !     ! so only compute for one flux tube index
         !     if (it /= 1) cycle
         !     iky = iky_idx(kxkyz_lo,ikxkyz)
         !     ikx = ikx_idx(kxkyz_lo,ikxkyz)
         !     iz = iz_idx(kxkyz_lo,ikxkyz)
         !     is = is_idx(kxkyz_lo,ikxkyz)
         !     g1 = aj0v(:,ikxkyz)*aj1v(:,ikxkyz)*(spec(is)%smz)**2 &
         !        * (kperp2(iky,ikx,ia,iz)*vperp2(ia,iz,:)/bmag(ia,iz)**2) &
         !        * (dkperp2dr(iky,ikx,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
         !        / (1.0 - aj0v(:,ikxkyz)**2 + 100.*epsilon(0.0))
         !
         !     g0 = spread((1.0 - aj0v(:,ikxkyz)**2),1,nvpa) &
         !        * spread(maxwell_vpa(:,is),2,nmu)*spread(maxwell_mu(ia,iz,:,is),1,nvpa)*maxwell_fac(is) &
         !        * (-spec(is)%tprim*(spread(vpa**2,2,nmu)+spread(vperp2(ia,iz,:),1,nvpa)-2.5) &
         !           -spec(is)%fprim + (dBdrho(iz)/bmag(ia,iz))*(1.0 - 2.0*spread(mu,1,nvpa)*bmag(ia,iz)) &
         !           + spread(g1,1,nvpa))
         !     wgt = spec(is)%z*spec(is)%z*spec(is)%dens/spec(is)%temp
         !     call integrate_vmu (g0, iz, tmp)
         !     dgamtotdr(iky,ikx,iz) = dgamtotdr(iky,ikx,iz) + tmp*wgt
         !   end do
         !   call sum_allreduce (dgamtotdr)
         !
         !   deallocate (g1)
         !
         ! endif
         ! ! avoid divide by zero when kx=ky=0
         ! ! do not evolve this mode, so value is irrelevant
         ! if (zonal_mode(1).and.akx(1)<epsilon(0.).and.has_electron_species(spec)) then
         !   gamtot(1,1,:)    = 0.0
         !   dgamtotdr(1,1,:) = 0.0
         !   zm = 1
         ! endif
         !
         ! if (.not.has_electron_species(spec)) then
         !    efac = tite/nine * (spec(ion_species)%dens/spec(ion_species)%temp)
         !    efacp = efac*(spec(ion_species)%tprim - spec(ion_species)%fprim)
         !    gamtot   = gamtot   + efac
         !    gamtot_h = gamtot_h + efac
         !    if(radial_variation) dgamtotdr = dgamtotdr + efacp
         !    if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         !       if (zonal_mode(1)) then
         !          gamtot3_h = tite/(nine*sum(spec%zt*spec%z*spec%dens))
         !          do ikx = 1, nakx
         !             ! avoid divide by zero for kx=ky=0 mode,
         !             ! which we do not need anyway
         !             !if (abs(akx(ikx)) < epsilon(0.)) cycle
         !             tmp = 1./efac - sum(dl_over_b(ia,:)/gamtot(1,ikx,:))
         !             gamtot3(ikx,:) = 1./(gamtot(1,ikx,:)*tmp)
         !             if (radial_variation) then
         !               tmp2 = (spec(ion_species)%tprim - spec(ion_species)%fprim)/efac &
         !                      + sum(d_dl_over_b_drho(ia,:)/gamtot(1,ikx,:)) &
         !                      - sum(dl_over_b(ia,:)*dgamtotdr(1,ikx,:) &
         !                          / gamtot(1,ikx,:)**2)
         !               dgamtot3dr(ikx,:)  = gamtot3(ikx,:) &
         !                                  * (-dgamtotdr(1,ikx,:)/gamtot(1,ikx,:) + tmp2/tmp)
         !             endif
         !          end do
         !          if(akx(1)<epsilon(0.)) then
         !             gamtot3(1,:)    = 0.0
         !             dgamtot3dr(1,:) = 0.0
         !             zm = 1
         !          endif
         !       end if
         !    end if
      end if

!    filename=trim(run_name)//".gamtot"
!    open(3636,file=trim(filename),status='unknown')
!    do iky = 1, naky
!      do ikx = 1, nakx
!        do iz = -nzgrid,nzgrid
!          write(3636,'(4e25.8)') gamtot(iky,ikx,iz), dgamtotdr(iky,ikx,iz), &
!                                 gamtot3(ikx,iz), dgamtot3dr(ikx,iz)
!        enddo
!      enddo
!    enddo
!    close(3636)

   end subroutine init_fields

   subroutine allocate_arrays

      use fields_arrays, only: phi, apar, bpar, phi_old
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
      if (.not. allocated(apar)) then
         allocate (apar(naky, nakx, -nzgrid:nzgrid, ntubes))
         apar = 0.
      end if
      if (.not. allocated(bpar)) then
         allocate (bpar(naky, nakx, -nzgrid:nzgrid, ntubes))
         bpar = 0.
      end if
      if (.not. allocated(phi_old)) then
         allocate (phi_old(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi_old = 0.
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

   subroutine advance_fields(g, phi, apar, bpar, dist)

      use mp, only: proc0
      use stella_layouts, only: vmu_lo
      use job_manage, only: time_message
      use redistribute, only: scatter
      use dist_fn_arrays, only: gvmu
      use zgrid, only: nzgrid
      use dist_redistribute, only: kxkyz2vmu
      use run_parameters, only: fields_kxkyz

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(inout) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar, bpar
      character(*), intent(in) :: dist

      if (fields_updated) return

      ! time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')
      if (fields_kxkyz) then
         ! first gather (vpa,mu) onto processor for v-space operations
         ! v-space operations are field solve, dg/dvpa, and collisions
         if (debug) write (*, *) 'dist_fn::advance_stella::scatter'
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' fields_redist')
         call scatter(kxkyz2vmu, g, gvmu)
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' fields_redist')
         ! given gvmu with vpa and mu local, calculate the corresponding fields
         if (debug) write (*, *) 'dist_fn::advance_stella::get_fields'
         call get_fields(gvmu, phi, apar, bpar, dist)
      else
         call get_fields_vmulo(g, phi, apar, bpar, dist)
      end if

      ! set a flag to indicate that the fields have been updated
      ! this helps avoid unnecessary field solves
      fields_updated = .true.
      ! time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

   end subroutine advance_fields

   subroutine get_fields(g, phi, apar, bpar, dist, skip_fsa)

      use mp, only: proc0
      use mp, only: sum_allreduce, mp_abort
      use job_manage, only: time_message
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use dist_fn_arrays, only: kperp2
      use gyro_averages, only: gyro_average, gyro_average_j1
      use run_parameters, only: fphi, fapar, fbpar
      use physics_parameters, only: beta
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: vpa, mu
      use vpamu_grids, only: integrate_vmu
      use kt_grids, only: nakx, naky
      use species, only: spec
      use constants, only: pi
      use fields_arrays, only: gamtot, gamtot13, gamtot31, gamtot33, apar_denom
      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(inout) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar, bpar
      logical, optional, intent(in) :: skip_fsa
      character(*), intent(in) :: dist
      complex :: tmp

      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      complex, dimension(:, :, :, :), allocatable :: antot1, antot2, antot3
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      logical :: skip_fsa_local

      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Bob: In the electromagnetic stella form, the field equations are
    !!! coupled, so each field calculation requires multiple integrals
    !!! of the distribution function. Thus it makes sense to calculate the
    !!! 3 integrals first, before calculating the fields:
    !!!
    !!! antot1 = \sum_s Z_s * dens_s * integrate_vmu(gyro_average(gbar) )
    !!! antot2 = beta * \sum_s Z_s * dens_s * v_{th,s,norm} * integrate_vmu(vpa * gyro_average(gbar) )
    !!! antot3 = -2 * beta * \sum_s dens_s * temp_s * integrate_vmu(mu * gyro_average1(gbar))
    !!!
    !!! "antot" variable name used because similar integrals with this name appear in GS2.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (dist == 'gbar') then
         allocate (g0(nvpa, nmu))
         allocate (antot1(naky, nakx, -nzgrid:nzgrid, ntubes))
         allocate (antot2(naky, nakx, -nzgrid:nzgrid, ntubes))
         allocate (antot3(naky, nakx, -nzgrid:nzgrid, ntubes))
         ! Initialise fields
         phi = 0.
         apar = 0.
         bpar = 0.
         if (fphi > epsilon(0.0)) then
            antot1 = 0.

            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               iz = iz_idx(kxkyz_lo, ikxkyz)
               it = it_idx(kxkyz_lo, ikxkyz)
               ikx = ikx_idx(kxkyz_lo, ikxkyz)
               iky = iky_idx(kxkyz_lo, ikxkyz)
               is = is_idx(kxkyz_lo, ikxkyz)
               call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)

               !antot1 = \sum_s Z_s * dens_s * integrate_vmu(gyro_average(ghat) )
               wgt = spec(is)%z * spec(is)%dens_psi0
               call integrate_vmu(g0, iz, tmp)
               antot1(iky, ikx, iz, it) = antot1(iky, ikx, iz, it) + wgt * tmp

            end do
            call sum_allreduce(antot1)
         end if

         if (fbpar > epsilon(0.0)) then
            antot3 = 0.

            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               iz = iz_idx(kxkyz_lo, ikxkyz)
               it = it_idx(kxkyz_lo, ikxkyz)
               ikx = ikx_idx(kxkyz_lo, ikxkyz)
               iky = iky_idx(kxkyz_lo, ikxkyz)
               is = is_idx(kxkyz_lo, ikxkyz)

               ! antot3 = -2 * beta * \sum_s dens_s * temp_s * integrate_vmu(mu * gyro_average1(gbar))
               call gyro_average_j1(g(:, :, ikxkyz), ikxkyz, g0)
               wgt = -2 * beta * spec(is)%dens_psi0 * spec(is)%temp_psi0
               call integrate_vmu((g0 * spread(mu, 1, nvpa)), iz, tmp)
               antot3(iky, ikx, iz, it) = antot3(iky, ikx, iz, it) + wgt * tmp

            end do
            call sum_allreduce(antot3)
         end if

         ! The fields are a funcion of (ky, kx, z, itube).
         ! antot variables are functions of (iky,ikx,iz,it)
         ! gamtot variables are functions of (iky, ikx, iz)
         ! ia = 1

         if (fphi > epsilon(0.0)) then
            if (fbpar > epsilon(0.0)) then
               phi = (antot1 - (spread(gamtot13, 4, ntubes) / spread(gamtot33, 4, ntubes)) * antot3) &
                     / (spread(gamtot, 4, ntubes) - (spread(gamtot13, 4, ntubes) * spread(gamtot31, 4, ntubes) / spread(gamtot33, 4, ntubes)))
               bpar = (antot3 - (spread(gamtot31, 4, ntubes) / spread(gamtot, 4, ntubes)) * antot1) &
                      / (spread(gamtot33, 4, ntubes) - (spread(gamtot13, 4, ntubes) * spread(gamtot31, 4, ntubes)) / spread(gamtot, 4, ntubes))
            else
               phi = (antot1 / (spread(gamtot, 4, ntubes)))
            end if

         else
            if (fbpar > epsilon(0.0)) then
               bpar = antot3 / (spread(gamtot33, 4, ntubes))
            end if
         end if

         if (fapar > epsilon(0.0)) then
            antot2 = 0.

            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               iz = iz_idx(kxkyz_lo, ikxkyz)
               it = it_idx(kxkyz_lo, ikxkyz)
               ikx = ikx_idx(kxkyz_lo, ikxkyz)
               iky = iky_idx(kxkyz_lo, ikxkyz)
               is = is_idx(kxkyz_lo, ikxkyz)
               call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)

               ! antot2 = beta * \sum_s Z_s * dens_s * v_{th,s,norm} * integrate_vmu(vpa * gyro_average(ghat) )
               wgt = beta * spec(is)%z * spec(is)%dens_psi0 * spec(is)%stm_psi0
               ! Bob: g0 has dimensions (ivpa, imu)
               ! So want to multiply g0 by vpa, spread out nmu times over axis 1
               call integrate_vmu((g0 * spread(vpa, 2, nmu)), iz, tmp)
               antot2(iky, ikx, iz, it) = antot2(iky, ikx, iz, it) + wgt * tmp

            end do

            call sum_allreduce(antot2)
            apar = antot2 / spread(apar_denom, 4, ntubes)

         end if

         deallocate (antot1)
         deallocate (antot2)
         deallocate (antot3)
         deallocate (g0)

      else
         ia = 1

         phi = 0.
         if (fphi > epsilon(0.0)) then
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

            call get_phi(phi, dist)

         end if

         apar = 0.
         if (fapar > epsilon(0.0)) then
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
            if (dist == 'h') then
               apar = apar / spread(kperp2(:, :, ia, :), 4, ntubes)
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

      !! Bob: This may be handy (calculating bpar for dist. fn g), but not for now.
      !! Calculate bpar, using g as the distribution function.
         ! bpar(kx, ky, z, t) = -beta/sqrt(pi) sum_s dens * { (8*bmag*temp_s integral(dvpa * dmu_s * g_sk * mu_s * aj1) )
         !                              + Z_s/bmag * phi * gamone  }
         ! bpar = 0.
         ! if (fbpar > epsilon(0.0)) then
         !   allocate (g0(nvpa,nmu))
         !   do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         !      iz = iz_idx(kxkyz_lo,ikxkyz)
         !      it = it_idx(kxkyz_lo,ikxkyz)
         !      ikx = ikx_idx(kxkyz_lo,ikxkyz)
         !      iky = iky_idx(kxkyz_lo,ikxkyz)
         !      is = is_idx(kxkyz_lo,ikxkyz)
         !      !!!!! Bob: Calculate bpar here
         !      ! integral = integral(dvpa * dmu_s * g_sk * mu_s * aj1)
         !      ! tmp = dens * ( 8*bmag*temp_s*integral + Z_s/bmag * phi * gamone)
         !      ! bpar(iky,ikx,iz,it) = bpar(iky,ikx,iz,it) + tmp
         !      ! g has dimensions (vpa, mu, ikxkyz), where ikxkyz contains kx, ky, z, ispecies, itube
         !      ! We're looping over ikxkyz, so everything should be a fn of vpa, mu or a single quantity
         !      ! integral should be over vpa, mu, so the integrand (let's call it g0) should be a function of vpa, mu only
         !      !!!!!!!!
         !      g0 = g(:,:,ikxkyz) * spread(mu, 1, nvpa) * spread(aj1v(:,ikxkyz),1,nvpa)
         !      ! write(*,*) "g0 = ", g0
         !      call integrate_vmu (g0, iz, tmp)
         !      ! write(*,*) "tmp = ", tmp
         !      !tmp = dens * (8*bmag*temp_s*integral + Z_s/bmag * phi * gamone)
         !      tmp = spec(is)%dens * ( 8*bmag(ia,iz)*spec(is)%temp*tmp + spec(is)%z/bmag(ia,iz) * phi(iky,ikx,iz,it) * gamone(iky,ikx,iz))
         !      ! write(*,*) "tmp = ", tmp
         !      ! gamone(iky,ikx,iz)
         !      bpar(iky,ikx,iz,it) = bpar(iky,ikx,iz,it) + tmp
         !   end do
         !   ! write(*,*) "bpar = ", bpar
         !   ! stop "Aborting now"
         !   !call sum_allreduce (apar)
         !   if (dist == 'h') then
         !      !apar = apar/spread(kperp2(:,:,ia,:),4,ntubes)
         !      write (*,*) 'BPAR NOT SETUP FOR H YET. aborting.'
         !      call mp_abort('BPAR NOT SETUP FOR H YET. aborting.')
         !   else if (dist == 'gbar') then
         !      !apar = apar/spread(apar_denom,4,ntubes)
         !   else if (dist == 'gstar') then
         !      write (*,*) 'BPAR NOT SETUP FOR GSTAR YET. aborting.'
         !      call mp_abort('BPAR NOT SETUP FOR GSTAR YET. aborting.')
         !   else
         !      if (proc0) write (*,*) 'unknown dist option in get_fields. aborting'
         !      call mp_abort ('unknown dist option in get_fields. aborting')
         !   end if
         !   deallocate (g0)
         ! end if

         !write(*,*) "bpar = ", bpar
      end if

   end subroutine get_fields

   subroutine get_fields_vmulo(g, phi, apar, bpar, dist, skip_fsa)

      use mp, only: mp_abort, sum_allreduce, proc0
      use job_manage, only: time_message
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: imu_idx, is_idx, iv_idx
      use gyro_averages, only: gyro_average, gyro_average_j1, aj0x, aj1x
      use run_parameters, only: fphi, fapar, fbpar
      use physics_parameters, only: beta
      use stella_geometry, only: dBdrho, bmag
      use physics_flags, only: radial_variation
      use dist_fn_arrays, only: kperp2, dkperp2dr
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: integrate_species, vperp2
      use vpamu_grids, only: vpa, mu, wgts_mu, wgts_vpa
      use kt_grids, only: nakx, naky, multiply_by_rho
      use run_parameters, only: ky_solve_radial
      use species, only: spec
      use fields_arrays, only: gamtot, gamtot13, gamtot31, gamtot33, apar_denom

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar, bpar
      logical, optional, intent(in) :: skip_fsa
      character(*), intent(in) :: dist
      real, dimension(:), allocatable :: wgts
      integer :: ivmu, iz, it, ia, imu, is, iv, iky
      logical :: skip_fsa_local
      complex, dimension(:, :, :), allocatable :: gyro_g
      complex, dimension(:, :), allocatable :: g0k
      complex, dimension(:, :, :, :), allocatable :: antot1, antot2, antot3

      ! Implement electromagnetic field calculations. Currently incompatible
      ! with global stella, so if this occurs, abort.
      if (dist == 'gbar') then
         if (radial_variation) then
            call mp_abort('radial_variation not supported for dist=gbar. Aborting.')
         end if

         skip_fsa_local = .false.
         if (present(skip_fsa)) skip_fsa_local = skip_fsa

         if (debug) write (*, *) 'dist_fn::advance_stella::get_fields_vmulo'

         ia = 1

         allocate (gyro_g(naky, nakx, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

         ! Initialise fields
         phi = 0.
         apar = 0.
         bpar = 0.
         if (fphi > epsilon(0.0)) then
            allocate (antot1(naky, nakx, -nzgrid:nzgrid, ntubes))
            antot1 = 0.
            wgts = spec%z * spec%dens_psi0
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               is = is_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               iv = iv_idx(vmu_lo, ivmu)
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
           antot1(:, :, iz, it) = antot1(:, :, iz, it) + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * wgts(is) * g(:, :, iz, it, ivmu) * aj0x(:, :, iz, ivmu)
                     !antot1 = \sum_s Z_s * dens_s * integrate_vmu(gyro_average(ghat) )
                  end do
               end do
            end do
            call sum_allreduce(antot1)
         end if

         if (fbpar > epsilon(0.0)) then
            allocate (antot3(naky, nakx, -nzgrid:nzgrid, ntubes))

            antot3 = 0.
            wgts = -2 * beta * spec%dens_psi0 * spec%temp_psi0
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               is = is_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               iv = iv_idx(vmu_lo, ivmu)
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
                     ! antot3 = -2 * beta * \sum_s dens_s * temp_s * integrate_vmu(mu * gyro_average1(ghat))
 antot3(:, :, iz, it) = antot3(:, :, iz, it) + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * wgts(is) * g(:, :, iz, it, ivmu) * aj0x(:, :, iz, ivmu) * mu(imu)
                  end do
               end do
            end do
            call sum_allreduce(antot3)
         end if

         ! The fields are a funcion of (ky, kx, z, itube).
         ! antot variables are functions of (iky,ikx,iz,it)
         ! gamtot variables are functions of (iky, ikx, iz)
         ! ia = 1

         if (fphi > epsilon(0.0)) then
            if (fbpar > epsilon(0.0)) then
               phi = (antot1 - (spread(gamtot13, 4, ntubes) / spread(gamtot33, 4, ntubes)) * antot3) &
                     / (spread(gamtot, 4, ntubes) - (spread(gamtot13, 4, ntubes) * spread(gamtot31, 4, ntubes) / spread(gamtot33, 4, ntubes)))
               bpar = (antot3 - (spread(gamtot31, 4, ntubes) / spread(gamtot, 4, ntubes)) * antot1) &
                      / (spread(gamtot33, 4, ntubes) - (spread(gamtot13, 4, ntubes) * spread(gamtot31, 4, ntubes)) / spread(gamtot, 4, ntubes))
            else
               phi = (antot1 / (spread(gamtot, 4, ntubes)))
            end if

         else
            if (fbpar > epsilon(0.0)) then
               bpar = antot3 / (spread(gamtot33, 4, ntubes))
            end if
         end if

         if (fapar > epsilon(0.0)) then
            allocate (antot2(naky, nakx, -nzgrid:nzgrid, ntubes))
            antot2 = 0.

            wgts = beta * spec%z * spec%dens_psi0 * spec%stm_psi0
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               is = is_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               iv = iv_idx(vmu_lo, ivmu)
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
                     ! antot2 = beta * \sum_s Z_s * dens_s * v_{th,s,norm} * integrate_vmu(vpa * gyro_average(ghat) )
 antot2(:, :, iz, it) = antot2(:, :, iz, it) + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * wgts(is) * g(:, :, iz, it, ivmu) * aj0x(:, :, iz, ivmu) * vpa(iv)
                  end do
               end do
            end do
            call sum_allreduce(antot2)
            apar = antot2 / spread(apar_denom, 4, ntubes)
         end if

         if (allocated(antot1)) deallocate (antot1)
         if (allocated(antot2)) deallocate (antot2)
         if (allocated(antot3)) deallocate (antot3)
         deallocate (gyro_g)

      else

         ia = 1
         phi = 0.
         if (fphi > epsilon(0.0)) then ! Normally, fphi = 1.0
            allocate (g0k(naky, nakx))
            allocate (gyro_g(naky, nakx, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                     is = is_idx(vmu_lo, ivmu)
                     imu = imu_idx(vmu_lo, ivmu)
                     call gyro_average(g(:, :, iz, it, ivmu), iz, ivmu, gyro_g(:, :, ivmu))
                     g0k = 0.0
                     if (radial_variation) then
                        do iky = 1, min(ky_solve_radial, naky)
                           g0k(iky, :) = gyro_g(iky, :, ivmu) &
                                         * (-0.5 * aj1x(iky, :, iz, ivmu) / aj0x(iky, :, iz, ivmu) * (spec(is)%smz)**2 &
                                            * (kperp2(iky, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                            * (dkperp2dr(iky, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                            + dBdrho(iz) / bmag(ia, iz))

                        end do
                        call multiply_by_rho(g0k)

                     end if

                     gyro_g(:, :, ivmu) = gyro_g(:, :, ivmu) + g0k

                  end do
                  call integrate_species(gyro_g, iz, spec%z * spec%dens_psi0, phi(:, :, iz, it), reduce_in=.false.)
               end do
            end do
            deallocate (gyro_g)
            call sum_allreduce(phi)

            call get_phi(phi, dist, skip_fsa_local)

         end if

         apar = 0.
         if (fapar > epsilon(0.0)) then
            ! FLAG -- NEW LAYOUT NOT YET SUPPORTED !!
            call mp_abort('APAR NOT YET SUPPORTED FOR NEW FIELD SOLVE. ABORTING.')
            !        allocate (g0(-nvgrid:nvgrid,nmu))
            !        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            !           iz = iz_idx(kxkyz_lo,ikxkyz)
            !           ikx = ikx_idx(kxkyz_lo,ikxkyz)
            !           iky = iky_idx(kxkyz_lo,ikxkyz)
            !           is = is_idx(kxkyz_lo,ikxkyz)
            !           g0 = spread(aj0v(:,ikxkyz),1,nvpa)*spread(vpa,2,nmu)*g(:,:,ikxkyz)
            !           wgt = 2.0*beta*spec(is)%z*spec(is)%dens*spec(is)%stm
            !           call integrate_vmu (g0, iz, tmp)
            !           apar(iky,ikx,iz) = apar(iky,ikx,iz) + tmp*wgt
            !        end do
            !        call sum_allreduce (apar)
            !        if (dist == 'h') then
            !           apar = apar/kperp2
            !        else if (dist == 'gbar') then
            !           apar = apar/apar_denom
            !        else if (dist == 'gstar') then
            !           write (*,*) 'APAR NOT SETUP FOR GSTAR YET. aborting.'
            !           call mp_abort('APAR NOT SETUP FOR GSTAR YET. aborting.')
            !        else
            !           if (proc0) write (*,*) 'unknown dist option in get_fields. aborting'
            !           call mp_abort ('unknown dist option in get_fields. aborting')
            !        end if
            !        deallocate (g0)
         end if
      end if

   end subroutine get_fields_vmulo

   subroutine get_fields_vmulo_single(g, iky, ikx, iz, phi, apar, bpar, dist)

      use mp, only: mp_abort, sum_allreduce
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: imu_idx, is_idx, iv_idx
      use gyro_averages, only: gyro_average, gyro_average_j1
      use run_parameters, only: fphi, fapar, fbpar
      use physics_parameters, only: beta
      use physics_flags, only: radial_variation
      use zgrid, only: ntubes
      use vpamu_grids, only: integrate_species
      use vpamu_grids, only: vpa, mu
      use kt_grids, only: nakx, naky, multiply_by_rho
      use run_parameters, only: ky_solve_radial
      use species, only: spec
      use fields_arrays, only: gamtot, gamtot13, gamtot31, gamtot33, apar_denom

      implicit none

      complex, dimension(vmu_lo%llim_proc:), intent(in) :: g
      integer, intent(in) :: iky, ikx, iz
      complex, intent(out) :: phi, apar, bpar
      character(*), intent(in) :: dist

      integer :: ivmu, iv, imu
      complex, dimension(:), allocatable :: gyro_g
      complex :: antot1, antot2, antot3

      ! Implement electromagnetic field calculations. Currently incompatible
      ! with global stella, so if this occurs, abort.
      if (dist == 'gbar') then
         if (radial_variation) then
            call mp_abort('radial_variation not supported for dist=gbar. Aborting.')
         end if

         allocate (gyro_g(vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         phi = 0.
         apar = 0.
         bpar = 0.
         if (fphi > epsilon(0.0)) then
            antot1 = 0.

            call gyro_average(g, iky, ikx, iz, gyro_g)
            !antot1 = \sum_s Z_s * dens_s * integrate_vmu(gyro_average(ghat) )
            call integrate_species(gyro_g, iz, spec%z * spec%dens_psi0, antot1, reduce_in=.false.)
            call sum_allreduce(antot1)
         end if

         if (fbpar > epsilon(0.0)) then
            antot3 = 0.
            call gyro_average_j1(g, iky, ikx, iz, gyro_g)
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               imu = imu_idx(vmu_lo, ivmu)
               gyro_g(ivmu) = gyro_g(ivmu) * mu(imu)
            end do
            ! antot3 = -2 * beta * \sum_s dens_s * temp_s * integrate_vmu(mu * gyro_average1(ghat))
            call integrate_species(gyro_g, iz, (-2 * beta * spec%dens_psi0 * spec%temp_psi0), antot3, reduce_in=.false.)
            call sum_allreduce(antot3)
         end if

         ! The fields are a funcion of (ky, kx, z, itube).
         ! antot variables are functions of (iky,ikx,iz,it)
         ! gamtot variables are functions of (iky, ikx, iz)
         ! ia = 1

         if (fphi > epsilon(0.0)) then
            if (fbpar > epsilon(0.0)) then
               phi = (antot1 - gamtot13(iky, ikx, iz) / gamtot33(iky, ikx, iz) * antot3) &
                     / (gamtot(iky, ikx, iz) - (gamtot13(iky, ikx, iz) * gamtot31(iky, ikx, iz)) / gamtot33(iky, ikx, iz))
               bpar = (antot3 - gamtot31(iky, ikx, iz) / gamtot(iky, ikx, iz) * antot1) &
                      / (gamtot33(iky, ikx, iz) - (gamtot13(iky, ikx, iz) * gamtot31(iky, ikx, iz)) / gamtot(iky, ikx, iz))
            else
               phi = antot1 / gamtot(iky, ikx, iz)
            end if

         else
            bpar = antot3 / gamtot33(iky, ikx, iz)
         end if

         if (fapar > epsilon(0.0)) then
            antot2 = 0.

            call gyro_average(g, iky, ikx, iz, gyro_g)
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               gyro_g(ivmu) = gyro_g(ivmu) * vpa(iv)
            end do

            ! antot2 = beta * \sum_s Z_s * dens_s * v_{th,s,norm} * integrate_vmu(vpa * gyro_average(ghat) )
            call integrate_species(gyro_g, iz, (beta * spec%z * spec%dens_psi0 * spec%stm_psi0), antot2, reduce_in=.false.)
            call sum_allreduce(antot2)
            apar = antot2 / apar_denom(iky, ikx, iz)
         end if
         deallocate (gyro_g)

      else
         ! No reason to be here - call abort
         call mp_abort('get_fields_vmulo_single_kxkyz called with dist!=gbar (unexpected behaviour). Aborting.')

      end if
   end subroutine get_fields_vmulo_single

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
      use dist_fn, only: adiabatic_option_switch
      use dist_fn, only: adiabatic_option_fieldlineavg
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

   subroutine get_fields_by_spec(g, fld, skip_fsa)

      use mp, only: sum_allreduce, mp_abort
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use gyro_averages, only: gyro_average
      use run_parameters, only: fphi, fapar, fbpar
      use stella_geometry, only: dl_over_b
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: integrate_vmu
      use kt_grids, only: nakx
      use kt_grids, only: zonal_mode
      use species, only: spec, nspec, has_electron_species
      use dist_fn, only: adiabatic_option_switch
      use dist_fn, only: adiabatic_option_fieldlineavg

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
         if (fapar > epsilon(0.0)) then
            ! Abort
            call mp_abort('get_fields_by_spec currently doesn`t support apar. aborting')
         end if
         if (fbpar > epsilon(0.0)) then
            ! Abort
            call mp_abort('get_fields_by_spec currently doesn`t support bpar. aborting')
         end if
         deallocate (g0)
      end if

   end subroutine get_fields_by_spec
   ! Bob: The following subroutine takes the fields and returns chi
   subroutine get_chi(phi, apar, bpar, ivmu, chi)

      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: chi
      integer, intent(in) :: ivmu

      write (*, *) "Not currently supported!"
      stop "Stopping"

   end subroutine get_chi

   ! The following subroutine takes the fields(ky,kx,z,tube) and returns
   ! gyroaverage(chi)(ky,kx,z,tube) = (J0*phi - 2*vpa*vths*J0*apar + 4*mu*(T/Z)*(J1/gamma) * bpar)
   !
   subroutine get_gyroaverage_chi_4d(ivmu, phi, apar, bpar, gyro_chi)

      use gyro_averages, only: gyro_average, gyro_average_j1
      use stella_layouts, only: vmu_lo
      use vpamu_grids, only: vpa, mu
      use stella_layouts, only: imu_idx, is_idx, iv_idx
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use run_parameters, only: fphi, fapar, fbpar
      implicit none
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
      integer, intent(in) :: ivmu
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: gyro_chi
      integer :: is, imu, iv
      complex, dimension(:, :, :, :), allocatable :: gyro_field

      ! Get vpa, mu and species from ivmu
      is = is_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      iv = iv_idx(vmu_lo, ivmu)

      allocate (gyro_field(naky, nakx, -nzgrid:nzgrid, ntubes))

      if (fphi > 0.0) then
         call gyro_average(phi, ivmu, gyro_field)
         gyro_chi = fphi * gyro_field
      else
         gyro_chi = 0.
      end if

      if (fapar > 0.0) then
         call gyro_average(apar, ivmu, gyro_field)
         gyro_chi = gyro_chi - fapar * 2 * vpa(iv) * spec(is)%stm * gyro_field
      end if

      if (fbpar > 0.0) then
         call gyro_average_j1(bpar, ivmu, gyro_field)
         gyro_chi = gyro_chi + fbpar * 4 * mu(imu) * (spec(is)%tz) * gyro_field
      end if

      deallocate (gyro_field)

   end subroutine get_gyroaverage_chi_4d

   ! The following subroutine takes the fields(ky,kx) and returns
   ! gyroaverage(chi)(ky,kx) = (J0*phi - 2*vpa*vths*J0*apar + 4*mu*(T/Z)*(J1/gamma) * bpar)
   subroutine get_gyroaverage_chi_2d(ivmu, phi, apar, bpar, gyro_chi)

      use gyro_averages, only: gyro_average, gyro_average_j1
      use stella_layouts, only: vmu_lo
      use vpamu_grids, only: vpa, mu
      use stella_layouts, only: imu_idx, is_idx, iv_idx
      use species, only: spec
      use kt_grids, only: naky, nakx
      use run_parameters, only: fphi, fapar, fbpar

      implicit none

      complex, dimension(:, :), intent(in) :: phi, apar, bpar
      integer, intent(in) :: ivmu
      complex, dimension(:, :), intent(out) :: gyro_chi
      integer :: is, imu, iv
      complex, dimension(:, :), allocatable :: gyro_field

      ! Get vpa, mu and species from ivmu
      is = is_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      iv = iv_idx(vmu_lo, ivmu)

      allocate (gyro_field(naky, nakx))

      if (fphi > 0.0) then
         call gyro_average(phi, ivmu, gyro_field)
         gyro_chi = fphi * gyro_field
      else
         gyro_chi = 0.0
      end if

      if (fapar > 0.0) then
         call gyro_average(apar, ivmu, gyro_field)
         gyro_chi = gyro_chi - fapar * 2 * vpa(iv) * spec(is)%stm * gyro_field
      end if

      if (fbpar > 0.0) then
         call gyro_average_j1(bpar, ivmu, gyro_field)
         gyro_chi = gyro_chi + fbpar * 4 * mu(imu) * (spec(is)%tz) * gyro_field
      end if

      deallocate (gyro_field)

   end subroutine get_gyroaverage_chi_2d

   subroutine get_phi(phi, dist, skip_fsa)

      use mp, only: proc0, mp_abort, job
#if defined MPI && ISO_C_BINDING
      use mpi
      use mp, only: proc0, iproc, nproc, sgproc0, comm_sgroup
      use mp, only: curr_focus, sharedsubprocs, scope
      use zgrid, only: nztot
      use fields_arrays, only: qn_zf_window, phi_shared
      use mp_lu_decomposition, only: lu_matrix_multiply_local
#endif
      use job_manage, only: time_message
      use physics_flags, only: full_flux_surface, radial_variation
      use run_parameters, only: ky_solve_radial, ky_solve_real
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: swap_kxky_ordered, nakx, naky, rho_d_clamped, zonal_mode, boundary_size
      use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use stella_geometry, only: dl_over_b, d_dl_over_b_drho
      use dist_fn, only: adiabatic_option_switch
      use dist_fn, only: adiabatic_option_fieldlineavg
      use species, only: spec, has_electron_species
      use multibox, only: mb_get_phi
      use fields_arrays, only: gamtot, gamtot3, phi_solve, phizf_solve, phi_ext
      use fields_arrays, only: phi_proj, phi_proj_stage, theta
      use fields_arrays, only: exclude_boundary_regions_qn, exp_fac_qn, tcorr_source_qn
      use file_utils, only: runtype_option_switch, runtype_multibox
      use linear_solve, only: lu_back_substitution

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
      logical, optional, intent(in) :: skip_fsa
      integer :: ia, it, iz, ikx, iky, zmi
      integer :: inmat, naky_r
      complex, dimension(:, :), allocatable :: g0k, g1k, g0x, g0a
      complex :: tmp
      logical :: skip_fsa_local
      logical :: has_elec, adia_elec
#if defined MPI && ISO_C_BINDING
      integer :: counter, c_lo, c_hi, c_max, c_div, c_mod
      integer :: prior_focus, ierr
#endif

      character(*), intent(in) :: dist

      if (debug) write (*, *) 'dist_fn::advance_stella::get_phi'

      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ia = 1
      has_elec = has_electron_species(spec)
      adia_elec = .not. has_elec &
                  .and. adiabatic_option_switch == adiabatic_option_fieldlineavg

      if (proc0) call time_message(.false., time_field_solve(:, 4), ' get_phi')
      if (dist == 'h') then
         phi = phi / gamtot_h
      else if (dist == 'gbar') then
         if (full_flux_surface) then
!        ! need to invert 1-gamma0 operator
!        allocate (phi_swap(naky_all,ikx_max))
!        do it = 1, ntubes
!           do iz = -nzgrid, nzgrid
!              call swap_kxky_ordered (phi, phi_swap)
!              do ikx = 1, ikx_max
!                 ! if ky values uncoupled, then obtain phi by simple divide
!                 if (maxval(ia_max_gam0a(:,ikx,iz)) == 1) then
!                    phi_swap(:,ikx) = phi_swap(:,ikx)/gam0a(1,:,ikx,iz)
!                 else
!                    ! for ky values that require info from other ky values
!                    ! must solve linear system
!                    call
!                 end if
!              end do
!           end do
!        end do
!        deallocate (phi_swap)
         else if ((radial_variation .and. ky_solve_radial > 0 &
                   .and. runtype_option_switch /= runtype_multibox) &
                  .or. &!DSO -> sorry for this if statement
                  (radial_variation .and. ky_solve_radial > 0 .and. job == 1 &
                   .and. runtype_option_switch == runtype_multibox &
                   .and. .not. ky_solve_real)) then
            allocate (g0k(1, nakx))
            allocate (g0x(1, nakx))
            allocate (g0a(1, nakx))

            naky_r = min(naky, ky_solve_radial)
#if defined MPI && ISO_C_BINDING
            prior_focus = curr_focus
            call scope(sharedsubprocs)
            c_max = nztot * ntubes * naky_r
            c_div = c_max / nproc
            c_mod = mod(c_max, nproc)

            c_lo = iproc * c_div + 1 + min(iproc, c_mod)
            c_hi = c_lo + c_div - 1
            if (iproc < c_mod) c_hi = c_hi + 1

            call scope(prior_focus)
            counter = 0
            if (sgproc0) phi_shared = phi
            call mpi_win_fence(0, phi_shared_window, ierr)
#endif
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do iky = 1, naky_r
#if defined MPI && ISO_C_BINDING
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
#if defined MPI && ISO_C_BINDING
            call mpi_win_fence(0, phi_shared_window, ierr)
            phi = phi_shared
#endif

            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do iky = naky_r, naky
                     phi(iky, :, iz, it) = phi(iky, :, iz, it) / gamtot(iky, :, iz)
                  end do
               end do
            end do

            if (ky_solve_radial == 0 .and. any(gamtot(1, 1, :) < epsilon(0.))) &
               phi(1, 1, :, :) = 0.0

            deallocate (g0k, g0x, g0a)
         else if (radial_variation .and. ky_solve_radial > 0 .and. job == 1 &
                  .and. runtype_option_switch == runtype_multibox) then
            call mb_get_phi(phi, has_elec, adia_elec)
         else
            phi = phi / spread(gamtot, 4, ntubes)
            if (any(gamtot(1, 1, :) < epsilon(0.))) phi(1, 1, :, :) = 0.0
         end if
      else
         if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
         call mp_abort('unknown dist option in get_fields. aborting')
         return
      end if

      if (any(gamtot(1, 1, :) < epsilon(0.))) phi(1, 1, :, :) = 0.0
      if (proc0) call time_message(.false., time_field_solve(:, 4), ' get_phi')

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
         else if (dist == 'gbar') then
            if (radial_variation .and. ky_solve_radial > 0 .and. job == 1 &
                .and. runtype_option_switch == runtype_multibox .and. ky_solve_real) then
               !this is already taken care of in mb_get_phi
            elseif ((radial_variation .and. ky_solve_radial > 0 &
                     .and. runtype_option_switch /= runtype_multibox) &
                    .or. &
                    (radial_variation .and. ky_solve_radial > 0 .and. job == 1 &
                     .and. runtype_option_switch == runtype_multibox &
                     .and. .not. ky_solve_real)) then
               allocate (g0k(1, nakx))
               allocate (g1k(1, nakx))
               allocate (g0x(1, nakx))

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

#if defined MPI && ISO_C_BINDING
                  if (sgproc0) then
#endif
                     do iz = -nzgrid, nzgrid - 1
                        do ikx = 1, nakx
                           inmat = ikx + nakx * (iz + nzgrid)
                           phi_ext(inmat) = phi(1, ikx, iz, it)
                        end do
                     end do
#if defined MPI && ISO_C_BINDING
                  end if
                  call mpi_win_fence(0, qn_zf_window, ierr)
#endif

#if defined MPI && ISO_C_BINDING
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
            else
               if (radial_variation) then
                  do it = 1, ntubes
                     do ikx = 1, nakx
                        ! DSO - this is sort of hack in order to avoid extra communications
                        !       However, get_radial_correction should be called immediately
                        !       after advance_fields, so it should be ok...
                        save1(ikx, it) = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                        save2(ikx, it) = sum(d_dl_over_b_drho(ia, :) * phi(1, ikx, :, it))
                     end do
                  end do
               end if
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

   ! the following routine gets the correction in phi both from gyroaveraging and quasineutrality
   ! the output, phi,
   subroutine get_radial_correction(g, phi_in, dist)

      use mp, only: proc0, mp_abort, sum_allreduce
      use stella_layouts, only: vmu_lo
      use gyro_averages, only: gyro_average, gyro_average_j1
      use gyro_averages, only: aj0x, aj1x
      use run_parameters, only: fphi, ky_solve_radial
      use stella_geometry, only: dl_over_b, bmag, dBdrho
      use stella_layouts, only: imu_idx, is_idx
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: integrate_species, vperp2
      use kt_grids, only: nakx, nx, naky
      use kt_grids, only: zonal_mode, multiply_by_rho
      use species, only: spec, has_electron_species
      use fields_arrays, only: phi_corr_QN, phi_corr_GA
      use fields_arrays, only: gamtot, dgamtotdr
      use fields_arrays, only: gamtot3, dgamtot3dr
      use dist_fn_arrays, only: kperp2, dkperp2dr
      use dist_fn, only: adiabatic_option_switch
      use dist_fn, only: adiabatic_option_fieldlineavg

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi_in
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      character(*), intent(in) :: dist

      integer :: ikx, iky, ivmu, iz, it, ia, is, imu
      complex :: tmp
      complex, dimension(:, :, :, :), allocatable :: phi
      complex, dimension(:, :, :), allocatable :: gyro_g
      complex, dimension(:, :), allocatable :: g0k, g0x

      ia = 1

      if (fphi > epsilon(0.0)) then
         allocate (gyro_g(naky, nakx, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (g0k(naky, nakx))
         allocate (g0x(naky, nx))
         allocate (phi(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi = 0.
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
                           + dBdrho(iz) / bmag(ia, iz) - dgamtotdr(:, :, iz) / gamtot(:, :, iz))

                  g0k(1, 1) = 0.

                  call gyro_average(g0k, iz, ivmu, gyro_g(:, :, ivmu))
               end do
               call integrate_species(gyro_g, iz, spec%z * spec%dens_psi0, phi(:, :, iz, it), reduce_in=.false.)
            end do
         end do
         call sum_allreduce(phi)

         if (dist == 'gbar') then
            !call get_phi (phi)
            phi = phi / spread(gamtot, 4, ntubes)
            phi(1, 1, :, :) = 0.0
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
               if (dist == 'gbar') then
                  do it = 1, ntubes
                     do ikx = 1, nakx
                        tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                        phi(1, ikx, :, it) = phi(1, ikx, :, it) &
                                             + tmp * gamtot3(ikx, :) &
                                             + dgamtot3dr(ikx, :) * save1(ikx, it) &
                                             + gamtot3(ikx, :) * save2(ikx, it)
                     end do
                  end do
               else
                  if (proc0) write (*, *) 'unknown dist option in radial_correction. aborting'
                  call mp_abort('unknown dist option in radial_correction. aborting')
               end if
            end if
         end if

         !collect quasineutrality corrections in wavenumber space
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0k = phi(:, :, iz, it)
               call multiply_by_rho(g0k)
               phi_corr_QN(:, :, iz, it) = g0k
            end do
         end do
         !zero out the ones we have already solved for using the full method
         do iky = 1, min(ky_solve_radial, naky)
            phi_corr_QN(iky, :, :, :) = 0.0
         end do

         deallocate (phi)

         !collect gyroaveraging corrections in wavenumber space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call gyro_average_j1(phi_in(:, :, iz, it), iz, ivmu, g0k)
                  g0k = -g0k * (spec(is)%smz)**2 &
                        * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                        * 0.5 * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz))

                  call multiply_by_rho(g0k)
                  phi_corr_GA(:, :, iz, it, ivmu) = g0k
               end do
            end do
         end do

         deallocate (g0x, g0k)
         deallocate (gyro_g)

      end if

   end subroutine get_radial_correction

   ! Take phi, apar, bpar(ky,kx,z,tube) and return
   ! d<chi>/dy (ky,kx,z,tube,vmu)
   subroutine get_dchidy_4d(phi, apar, bpar, dchidy)

      use constants, only: zi
      use stella_layouts, only: vmu_lo
      use run_parameters, only: fphi, fapar
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, aky, naky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dchidy

      integer :: ivmu
      complex, dimension(:, :, :, :), allocatable :: gyro_chi
      complex, dimension(:, :, :, :), allocatable :: coefficient

      allocate (gyro_chi(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (coefficient(naky, nakx, -nzgrid:nzgrid, ntubes))
      coefficient = spread(spread(spread(zi * aky, 2, nakx), 3, 2 * nzgrid + 1), 4, ntubes)
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         call get_gyroaverage_chi(ivmu, phi, apar, bpar, gyro_chi)
         dchidy(:, :, :, :, ivmu) = coefficient * gyro_chi
      end do

      deallocate (gyro_chi, coefficient)

   end subroutine get_dchidy_4d

   ! Take phi, apar, bpar(ky, kx) and return
   ! d<chi>/dy (ky,kx)
   subroutine get_dchidy_2d(ivmu, phi, apar, bpar, dchidy)

      use constants, only: zi
      use kt_grids, only: nakx, aky, naky

      implicit none

      integer, intent(in) :: ivmu
      complex, dimension(:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :), intent(out) :: dchidy

      !integer :: iv, is
      complex, dimension(:, :), allocatable :: gyro_chi

      allocate (gyro_chi(naky, nakx))
      call get_gyroaverage_chi(ivmu, phi, apar, bpar, gyro_chi)
      dchidy = zi * spread(aky, 2, nakx) * gyro_chi
      deallocate (gyro_chi)

   end subroutine get_dchidy_2d

   ! Take phi, apar, bpar(ky,kx,z,tube) and return
   ! d<chi>/dx (ky,kx,z,tube,vmu)
   subroutine get_dchidx_4d(phi, apar, bpar, dchidx)

      use constants, only: zi
      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, akx, naky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dchidx

      integer :: ivmu
      complex, dimension(:, :, :, :), allocatable :: gyro_chi

      allocate (gyro_chi(naky, nakx, -nzgrid:nzgrid, ntubes))

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         call get_gyroaverage_chi(ivmu, phi, apar, bpar, gyro_chi)
         dchidx(:, :, :, :, ivmu) = zi * spread(spread(spread(akx, 1, naky), 3, 2 * nzgrid + 1), 4, ntubes) &
                                    * gyro_chi
      end do

      deallocate (gyro_chi)

   end subroutine get_dchidx_4d

   ! Take phi, apar, bpar(ky, kx) and return
   ! d<chi>/dx (ky,kx)
   subroutine get_dchidx_2d(ivmu, phi, apar, bpar, dchidx)

      use constants, only: zi
      use kt_grids, only: akx, naky, nakx

      implicit none

      integer, intent(in) :: ivmu
      complex, dimension(:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :), intent(out) :: dchidx

      complex, dimension(:, :), allocatable :: gyro_chi

      allocate (gyro_chi(naky, nakx))
      call get_gyroaverage_chi(ivmu, phi, apar, bpar, gyro_chi)
      dchidx = zi * spread(akx, 1, naky) * gyro_chi
      deallocate (gyro_chi)

   end subroutine get_dchidx_2d

   subroutine finish_fields

      use fields_arrays, only: phi, phi_old
      use fields_arrays, only: phi_corr_QN, phi_corr_GA
      use fields_arrays, only: apar, apar_corr_QN, apar_corr_GA, bpar
      use fields_arrays, only: gamtot, dgamtotdr, gamtot3, dgamtot3dr
      use fields_arrays, only: c_mat, theta
      use fields_arrays, only: gamtot, dgamtotdr, gamtot13, gamtot31, gamtot33
#if defined MPI && ISO_C_BINDING
      use fields_arrays, only: qn_window
      use mpi
#else
      use fields_arrays, only: phi_solve
#endif
      implicit none

#if defined MPI && ISO_C_BINDING
      integer ierr
#endif

      if (allocated(phi)) deallocate (phi)
      if (allocated(apar)) deallocate (apar)
      if (allocated(bpar)) deallocate (bpar)
      if (allocated(phi_old)) deallocate (phi_old)
      if (allocated(phi_corr_QN)) deallocate (phi_corr_QN)
      if (allocated(phi_corr_GA)) deallocate (phi_corr_GA)
      if (allocated(apar)) deallocate (apar)
      if (allocated(apar_corr_QN)) deallocate (apar_corr_QN)
      if (allocated(apar_corr_GA)) deallocate (apar_corr_GA)
      if (allocated(gamtot)) deallocate (gamtot)
      if (allocated(gamtot13)) deallocate (gamtot13)
      if (allocated(gamtot31)) deallocate (gamtot31)
      if (allocated(gamtot33)) deallocate (gamtot33)
      if (allocated(gamtot3)) deallocate (gamtot3)
      if (allocated(dgamtotdr)) deallocate (dgamtotdr)
      if (allocated(dgamtot3dr)) deallocate (dgamtot3dr)
      if (allocated(apar_denom)) deallocate (apar_denom)
      if (allocated(save1)) deallocate (save1)
      if (allocated(save2)) deallocate (save2)

#if defined MPI && ISO_C_BINDING
      if (phi_shared_window /= MPI_WIN_NULL) call mpi_win_free(phi_shared_window, ierr)
      if (qn_window_initialized .and. qn_window /= MPI_WIN_NULL) then
         call mpi_win_free(qn_window, ierr)
         qn_window_initialized = .false.
      end if
#else
      if (allocated(phi_solve)) deallocate (phi_solve)
#endif
      if (allocated(c_mat)) deallocate (c_mat)
      if (allocated(theta)) deallocate (theta)

      fields_initialized = .false.

   end subroutine finish_fields

end module fields
