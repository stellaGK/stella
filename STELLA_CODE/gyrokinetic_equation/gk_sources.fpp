!###############################################################################
!###############################################################################
!###############################################################################
! 
! This module ...
! 
!###############################################################################
module gk_sources

#ifdef ISO_C_BINDING
   use mpi
#endif

   implicit none

   public :: init_sources, finish_sources
   public :: init_quasineutrality_source
   public :: init_source_timeaverage
   public :: update_quasineutrality_source
   public :: source_option_switch, source_option_none
   public :: source_option_krook, source_option_projection
   public :: include_qn_source
   public :: update_tcorr_krook
   public :: project_out_zero
   public :: add_krook_operator
   public :: tcorr_source, exclude_boundary_regions, exp_fac
   public :: int_krook, int_proj
   public :: qn_source_initialized
   public :: time_sources

   private

   logical :: krook_odd, exclude_boundary_regions
   logical :: from_zero
   logical :: conserve_momentum, conserve_density
   integer:: ikxmax_source
   real :: nu_krook, tcorr_source, int_krook, int_proj
   real :: exp_fac

   logical :: qn_source_initialized, include_qn_source

   logical :: debug = .false.

   real, dimension(2, 2) :: time_sources = 0.

   integer :: source_option_switch
   integer, parameter :: source_option_none = 1, &
                         source_option_krook = 2, &
                         source_option_projection = 3

contains

   subroutine init_sources

      use mp, only: job, proc0
      use parameters_physics, only: fphi
      use parameters_multibox, only: ky_solve_radial, ky_solve_real
      use parameters_kxky_grid, only: naky, nakx
      use grids_kxky, only: zonal_mode
      use grids_z, only: nzgrid, ntubes
      use stella_layouts, only: vmu_lo
      use arrays_store_distribution_fn, only: g_krook, g_proj, g_symm
      use arrays_store_fields, only: phi_proj, phi_proj_stage
      use parameters_physics, only: radial_variation
      use grids_species, only: spec, has_electron_species
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg
      use file_utils, only: runtype_option_switch, runtype_multibox
      use parameters_kxky_grid, only: ikx_max
      use namelist_sources, only: read_namelist_sources
      use arrays_store_useful, only: tcorr_source_qn, exclude_boundary_regions_qn
      use mp, only: broadcast
      implicit none

      logical :: has_elec, adia_elec
      real :: fac

      debug = debug .and. proc0

      call read_namelist_sources(source_option_switch, nu_krook, tcorr_source, &
                              ikxmax_source, krook_odd, exclude_boundary_regions, &
                              tcorr_source_qn, exclude_boundary_regions_qn, from_zero, &
                              conserve_momentum, conserve_density)

      ikxmax_source = min(ikxmax_source, ikx_max)

      int_proj = -1.
      int_krook = -1.

      call broadcast_arrays

      if (source_option_switch == source_option_krook .and. .not. allocated(g_krook)) then
         allocate (g_krook(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         g_krook = 0.
      end if

      if (source_option_switch == source_option_projection .and. .not. allocated(g_proj)) then
         allocate (g_proj(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         g_proj = 0.
      end if

      if (.not. allocated(phi_proj)) then
         allocate (phi_proj(nakx, -nzgrid:nzgrid, ntubes)); phi_proj = 0.
      end if
      if (.not. allocated(phi_proj_stage)) then
         allocate (phi_proj_stage(nakx, -nzgrid:nzgrid, ntubes)); phi_proj_stage = 0.
      end if

      if ((conserve_momentum .or. conserve_density) .and. .not. allocated(g_symm)) then
         allocate (g_symm(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      end if

      fac = 1.
      if (from_zero) fac = 0.

      if (int_krook < 0.) int_krook = fac * tcorr_source
      if (int_proj < 0.) int_proj = fac * tcorr_source

      include_qn_source = .false.
      if (fphi > epsilon(0.0) .and. radial_variation .and. ky_solve_radial > 0) then
         has_elec = has_electron_species(spec)
         adia_elec = .not. has_elec .and. zonal_mode(1) &
                     .and. adiabatic_option_switch == adiabatic_option_fieldlineavg
         if (adia_elec) then
            if (runtype_option_switch /= runtype_multibox .or. (job == 1 .and. .not. ky_solve_real)) then
               include_qn_source = .true.
            end if
         end if
      end if

   contains

      subroutine broadcast_arrays

         implicit none

         call broadcast(source_option_switch)
         call broadcast(exclude_boundary_regions)
         call broadcast(exclude_boundary_regions_qn)
         call broadcast(nu_krook)
         call broadcast(tcorr_source)
         call broadcast(tcorr_source_qn)
         call broadcast(ikxmax_source)
         call broadcast(krook_odd)
         call broadcast(from_zero)
         call broadcast(conserve_momentum)
         call broadcast(conserve_density)

      end subroutine broadcast_arrays

   end subroutine init_sources

   subroutine init_source_timeaverage

      use stella_time, only: code_dt
      use arrays_store_useful, only: tcorr_source_qn, exp_fac_qn

      implicit none

      if (tcorr_source > 0.0) then
         exp_fac = exp(-code_dt / tcorr_source)
      else
         exp_fac = 0.0
      end if

      if (tcorr_source_qn > 0.0) then
         exp_fac_qn = exp(-code_dt / tcorr_source_qn)
      else
         exp_fac_qn = 0.0
      end if

   end subroutine init_source_timeaverage

   subroutine finish_sources

      use arrays_store_distribution_fn, only: g_krook, g_proj, g_symm
      use arrays_store_fields, only: phi_proj, phi_proj_stage
#ifdef ISO_C_BINDING
      use arrays_store_useful, only: qn_zf_window
#else
      use arrays_store_fields, only: phizf_solve, phi_ext
#endif

      implicit none

      integer :: ierr

      if (allocated(g_krook)) deallocate (g_krook)
      if (allocated(g_proj)) deallocate (g_proj)
      if (allocated(g_symm)) deallocate (g_symm)
      if (allocated(phi_proj)) deallocate (phi_proj)
      if (allocated(phi_proj_stage)) deallocate (phi_proj_stage)

#ifdef ISO_C_BINDING
      if (qn_zf_window /= MPI_WIN_NULL) then
         call mpi_win_free(qn_zf_window, ierr)
      end if
#else
      if (associated(phizf_solve%zloc)) deallocate (phizf_solve%zloc)
      if (associated(phizf_solve%idx)) deallocate (phizf_solve%idx)
      if (associated(phi_ext)) deallocate (phi_ext)
#endif

   end subroutine finish_sources

   subroutine add_krook_operator(g, gke_rhs)

      use mp, only: proc0
      use job_manage, only: time_message
      use grids_z, only: nzgrid, ntubes
      use constants, only: pi, zi
      use grids_kxky, only: akx, zonal_mode 
      use parameters_multibox, only: boundary_size
      use parameters_kxky_grid, only: nakx
      use stella_layouts, only: vmu_lo
      use stella_time, only: code_dt
      use arrays_store_distribution_fn, only: g_krook, g_symm
      use calculations_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

      implicit none

      complex :: tmp
      integer :: ikx, jkx, iz, it, ia, ivmu, npts

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in), target :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gke_rhs

      complex, dimension(:, :, :, :, :), pointer :: g_work

      complex, dimension(:, :), allocatable :: g0k, g0x, g1x
      real, dimension(:), allocatable :: basis_func

      ia = 1
      if (.not. zonal_mode(1)) return

      if (proc0) call time_message(.false., time_sources(:, 1), ' sources')

      if (debug) write (6, *) 'sources::add_krook_operator'

      g_work => g
      if (conserve_momentum .or. conserve_density) then
         g_work => g_symm
         g_work = g
      end if

      if (debug) write (6, *) 'sources::add_krook_operator::conservation'

      if (conserve_momentum) call enforce_momentum_conservation(g_work)
      if (conserve_density) call enforce_density_conservation(g_work(1, :, :, :, :))

      if (exclude_boundary_regions) then
         if (debug) write (6, *) 'sources::add_krook_operator::exclude_boundary_regions'
         npts = nakx - 2 * boundary_size
         allocate (g0k(1, nakx))
         allocate (g0x(1, nakx))
         allocate (g1x(1, nakx))
         allocate (basis_func(npts))
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g0k(1, :) = g_work(1, :, iz, it, ivmu)
                  g1x = 0.
                  call transform_kx2x_unpadded(g0k, g0x)
                  do ikx = 1, ikxmax_source
                     if (ikx == 1) then
                        basis_func = 1.0
                        tmp = sum(g0x(1, (boundary_size + 1):(nakx - boundary_size))) / real(npts)
                     else
                        do jkx = 1, npts
                           basis_func(jkx) = sin(2.0 * pi * (ikx - 1) * jkx / real(npts + 1))
                        end do
                        tmp = 2.0 * sum(basis_func * g0x(1, (boundary_size + 1):(nakx - boundary_size))) / real(npts + 1)
                     end if
                     if (tcorr_source > epsilon(0.0)) then
                        tmp = (code_dt * tmp + exp_fac * int_krook * g_krook(ikx, iz, it, ivmu)) &
                              / (code_dt + exp_fac * int_krook)
                     end if
                     do jkx = 1, npts
                        g1x(1, boundary_size + jkx) = g1x(1, boundary_size + jkx) + tmp * basis_func(jkx)
                     end do
                  end do
                  call transform_x2kx_unpadded(g1x, g0k)
                  gke_rhs(1, :, iz, it, ivmu) = gke_rhs(1, :, iz, it, ivmu) - code_dt * nu_krook * g0k(1, :)
               end do
            end do
         end do
         deallocate (g0k, g0x, g1x, basis_func)
      else
         if (debug) write (6, *) 'sources::add_krook_operator::include_boundary_regions'
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do ikx = 1, nakx
                     if (abs(akx(ikx)) > akx(ikxmax_source)) cycle
                     tmp = g_work(1, ikx, iz, it, ivmu)
                     if (krook_odd .and. abs(akx(ikx)) > epsilon(0.0)) tmp = zi * aimag(tmp)
                     if (tcorr_source <= epsilon(0.0)) then
                        gke_rhs(1, ikx, iz, it, ivmu) = gke_rhs(1, ikx, iz, it, ivmu) - code_dt * nu_krook * tmp
                     else
                        gke_rhs(1, ikx, iz, it, ivmu) = gke_rhs(1, ikx, iz, it, ivmu) - code_dt * nu_krook &
                                                        * (code_dt * tmp + exp_fac * int_krook * g_krook(ikx, iz, it, ivmu)) &
                                                        / (code_dt + exp_fac * int_krook)
                     end if
                  end do
               end do
            end do
         end do
      end if

      if (proc0) call time_message(.false., time_sources(:, 1), ' sources')

   end subroutine add_krook_operator

   subroutine update_tcorr_krook(g)

      use mp, only: proc0
      use job_manage, only: time_message
      use constants, only: pi, zi
      use arrays_store_distribution_fn, only: g_krook, g_symm
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: akx, zonal_mode
      use parameters_multibox, only: boundary_size
      use parameters_kxky_grid, only: nakx
      use stella_layouts, only: vmu_lo
      use stella_time, only: code_dt
      use calculations_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), target, intent(in) :: g
      complex, dimension(:, :), allocatable :: g0k, g0x

      complex, dimension(:, :, :, :, :), pointer :: g_work

      integer :: ivmu, iz, it, ikx, jkx, ia, npts
      real :: int_krook_old
      complex :: tmp

      if (.not. zonal_mode(1)) return

      if (proc0) call time_message(.false., time_sources(:, 1), ' sources')

      ia = 1

      if (debug) write (6, *) 'sources::update_tcorr_krook'

      g_work => g
      if (conserve_momentum .or. conserve_density) then
         g_work => g_symm
         g_work = g
      end if

      if (conserve_momentum) call enforce_momentum_conservation(g_work)
      if (conserve_density) call enforce_density_conservation(g_work(1, :, :, :, :))

      int_krook_old = int_krook
      int_krook = code_dt + exp_fac * int_krook_old

      if (exclude_boundary_regions) then
         npts = nakx - 2 * boundary_size
         allocate (g0k(1, nakx))
         allocate (g0x(1, nakx))
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g0k(1, :) = g_work(1, :, iz, it, ivmu)
                  call transform_kx2x_unpadded(g0k, g0x)
                  do ikx = 1, ikxmax_source
                     if (ikx == 1) then
                        tmp = sum(g0x(1, (boundary_size + 1):(nakx - boundary_size))) / real(npts)
                     else
                        tmp = 0.
                        do jkx = 1, npts
                           tmp = tmp + sin(2.0 * pi * (ikx - 1) * jkx / real(npts + 1)) * g0x(1, boundary_size + jkx)
                        end do
                        tmp = 2.0 * tmp / real(npts + 1)
                     end if
                     g_krook(ikx, iz, it, ivmu) = (code_dt * tmp + exp_fac * int_krook_old * g_krook(ikx, iz, it, ivmu)) &
                                                  / int_krook
                  end do
               end do
            end do
         end do
         deallocate (g0k, g0x)
      else
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do ikx = 1, nakx
                     tmp = g(1, ikx, iz, it, ivmu)
                     if (krook_odd .and. abs(akx(ikx)) > epsilon(0.0)) tmp = zi * aimag(tmp)
                     g_krook(ikx, iz, it, ivmu) = (code_dt * tmp + exp_fac * int_krook_old * g_krook(ikx, iz, it, ivmu)) / int_krook
                  end do
               end do
            end do
         end do
      end if

      if (proc0) call time_message(.false., time_sources(:, 1), ' sources')

   end subroutine update_tcorr_krook

   subroutine enforce_momentum_conservation(g_work)

      use mp, only: proc0
      use job_manage, only: time_message
      use redistribute, only: scatter, gather
      use stella_layouts, only: vmu_lo, kxkyz_lo
      use stella_layouts, only: imu_idx, is_idx, iv_idx
      use grids_velocity, only: nvgrid, nvpa, nmu
      use calculations_redistribute, only: kxkyz2vmu
      use arrays_store_distribution_fn, only: gvmu
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(inout) :: g_work

      integer :: ikxkyz, imu, iv, iv2
      complex :: tmp

      if (proc0) call time_message(.false., time_sources(:, 2), ' source_redist')
      call scatter(kxkyz2vmu, g_work, gvmu)
      if (proc0) call time_message(.false., time_sources(:, 2), ' source_redist')

      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         do imu = 1, nmu
            do iv = 1, nvgrid
               iv2 = nvpa - iv + 1
               tmp = 0.5 * (gvmu(iv, imu, ikxkyz) + gvmu(iv2, imu, ikxkyz))
               gvmu(iv, imu, ikxkyz) = tmp
               gvmu(iv2, imu, ikxkyz) = tmp
            end do
         end do
      end do

      if (proc0) call time_message(.false., time_sources(:, 2), ' source_redist')
      call gather(kxkyz2vmu, gvmu, g_work)
      if (proc0) call time_message(.false., time_sources(:, 2), ' source_redist')

   end subroutine enforce_momentum_conservation

   subroutine enforce_density_conservation(g_work)

      use mp, only: sum_allreduce
      use grids_species, only: spec
      use parameters_physics, only: radial_variation
      use calculations_velocity_integrals, only: integrate_species
      use grids_velocity, only: mu, vpa, vperp2
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use parameters_kxky_grid, only: nakx
      use grids_kxky, only: rho_d_clamped
      use stella_layouts, only: vmu_lo, imu_idx, is_idx, iv_idx
      use geometry, only: bmag, dBdrho, dl_over_b, d_dl_over_b_drho
      use calculations_gyro_averages, only: gyro_average
      use arrays_gyro_averages, only: aj0x, aj1x
      use arrays_store_useful, only: kperp2, dkperp2dr
      use grids_z, only: nzgrid, ntubes
      use calculations_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

      implicit none

      complex, dimension(:, -nzgrid:, :, vmu_lo%llim_proc:), intent(inout) :: g_work
      integer :: ia, ikx, it, iz, imu, iv, ivmu, is

      complex, dimension(:, :), allocatable :: gyro_g, g0k, g0x
      complex, dimension(:), allocatable :: g_fsa
      real :: energy

      ia = 1

      allocate (gyro_g(nakx, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      allocate (g0k(1, nakx))
      allocate (g0x(1, nakx))
      allocate (g_fsa(nakx)); g_fsa = 0.

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               is = is_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               iv = iv_idx(vmu_lo, ivmu)
               gyro_g(:, ivmu) = g_work(:, iz, it, ivmu) * aj0x(1, :, iz, ivmu)
               g0k = 0.0
               if (radial_variation) then
                  g0k(1, :) = gyro_g(:, ivmu) &
                              * (-0.5 * aj1x(1, :, iz, ivmu) / aj0x(1, :, iz, ivmu) * (spec(is)%smz)**2 &
                                 * (kperp2(1, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                 * (dkperp2dr(1, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                 + dBdrho(iz) / bmag(ia, iz))

                  call transform_kx2x_unpadded(g0k, g0x)
                  g0x(1, :) = rho_d_clamped * g0x(1, :)
                  call transform_x2kx_unpadded(g0x, g0k)
               end if
               gyro_g(:, ivmu) = gyro_g(:, ivmu) + g0k(1, :)
            end do

            do ikx = 1, nakx
               call integrate_species(gyro_g(ikx, :), iz, spec%dens_psi0, g0k(1, ikx), reduce_in=.false.)
            end do
            call sum_allreduce(g0k)

            !we now have delta n. Flux surface average
            call transform_kx2x_unpadded(g0k, g0x)
            g_fsa = g_fsa + dl_over_b(ia, iz) * g0x(1, :)
            if (radial_variation) g_fsa = g_fsa + d_dl_over_b_drho(ia, iz) * rho_d_clamped * g0x(1, :)
         end do
      end do
      g_fsa = g_fsa / ntubes

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         is = is_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         iv = iv_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0x(1, :) = g_fsa
               !multiply by f0
               g0x(1, :) = g0x(1, :) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
               if (radial_variation) then !variation in the density cancels
                  energy = (vpa(iv)**2 + vperp2(ia, iz, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
                  g0x(1, :) = g0x(1, :) * (1.0 - rho_d_clamped * (spec(is)%tprim * (energy - 1.5) + 2.*mu(imu) * dBdrho(iz)))
               end if

               call transform_x2kx_unpadded(g0x, g0k)
               g_work(:, iz, it, ivmu) = g_work(:, iz, it, ivmu) - g0k(1, :)
            end do
         end do
      end do
      deallocate (gyro_g, g0k, g0x, g_fsa)

   end subroutine enforce_density_conservation

   subroutine project_out_zero(gold, gnew)

      use mp, only: proc0
      use job_manage, only: time_message
      use grids_z, only: nzgrid, ntubes
      use constants, only: pi, zi
      use grids_kxky, only: zonal_mode, akx
      use parameters_multibox, only: boundary_size
      use parameters_kxky_grid, only: nakx
      use stella_layouts, only: vmu_lo
      use stella_time, only: code_dt
      use arrays_store_distribution_fn, only: g_proj, g_symm
      use calculations_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

      implicit none

      complex :: tmp
      integer :: ikx, jkx, iz, it, ia, ivmu, npts

      complex, dimension(:, :), allocatable :: g0k, g0x, g1x
      real, dimension(:), allocatable :: basis_func
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gold
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(inout) :: gnew
      complex, allocatable, dimension(:, :, :, :) :: g

      ia = 1
      if (.not. zonal_mode(1)) return

      if (debug) write (6, *) 'sources::project_out_zero'

      if (proc0) call time_message(.false., time_sources(:, 1), ' sources')

      allocate (g(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      !divide by code_dt to ensure time averaging is performed correctly
      if (conserve_momentum) then
         g_symm = (gnew - gold) / code_dt
         call enforce_momentum_conservation(g_symm)
         g = g_symm(1, :, :, :, :)

      else
         g = (gnew(1, :, :, :, :) - gold(1, :, :, :, :)) / code_dt
      end if

      if (conserve_density) call enforce_density_conservation(g)

      if (exclude_boundary_regions) then
         npts = nakx - 2 * boundary_size
         allocate (g0k(1, nakx))
         allocate (g0x(1, nakx))
         allocate (g1x(1, nakx))
         allocate (basis_func(npts))
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g0k(1, :) = g(:, iz, it, ivmu)
                  g1x = 0.
                  call transform_kx2x_unpadded(g0k, g0x)
                  do ikx = 1, ikxmax_source
                     !physical region should have an odd number of collocation points
                     if (ikx == 1) then
                        basis_func = 1.0
                        tmp = sum(g0x(1, (boundary_size + 1):(nakx - boundary_size))) / real(npts)
                     else
                        ! here we use a Fourier basis due to periodicity,
                        ! though we could use Legendre polynomials
                        ! NB: Only a constant or linear function (or nearly linear, i.e. first
                        ! sine harmonic) make physical sense as sources, so ikxmax_source <= 2
                        do jkx = 1, npts
                           basis_func(jkx) = sin(2.0 * pi * (ikx - 1) * jkx / real(npts + 1))
                        end do
                        tmp = 2.0 * sum(basis_func * g0x(1, (boundary_size + 1):(nakx - boundary_size))) / real(npts + 1)
                     end if
                     if (tcorr_source > epsilon(0.)) then
                        tmp = (code_dt * tmp + exp_fac * int_proj * g_proj(ikx, iz, it, ivmu)) &
                              / (code_dt + exp_fac * int_proj)
                        g_proj(ikx, iz, it, ivmu) = tmp
                     end if
                     do jkx = 1, npts
                        g1x(1, boundary_size + jkx) = g1x(1, boundary_size + jkx) + tmp * basis_func(jkx)
                     end do
                  end do
                  call transform_x2kx_unpadded(g1x, g0k)
                  g(:, iz, it, ivmu) = g0k(1, :)
               end do
            end do
         end do
         deallocate (g0k, g0x, g1x, basis_func)
      else
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do ikx = 1, nakx
                     if (abs(akx(ikx)) > akx(ikxmax_source)) then
                        g(ikx, iz, it, ivmu) = 0.0
                     else
                        tmp = g(ikx, iz, it, ivmu)
                        if (krook_odd .and. abs(akx(ikx)) > epsilon(0.0)) tmp = zi * aimag(tmp)
                        if (tcorr_source <= epsilon(0.)) then
                           g(ikx, iz, it, ivmu) = tmp
                        else
                           g(ikx, iz, it, ivmu) = (code_dt * tmp + exp_fac * int_proj * g_proj(ikx, iz, it, ivmu)) &
                                                  / (code_dt + exp_fac * int_proj)
                        end if
                     end if
                     if (krook_odd .and. abs(akx(ikx)) > epsilon(0.0)) then
                        g_proj(ikx, iz, it, ivmu) = zi * aimag(g(ikx, iz, it, ivmu))
                     else
                        g_proj(ikx, iz, it, ivmu) = g(ikx, iz, it, ivmu)
                     end if
                  end do
               end do
            end do
         end do
      end if

      int_proj = code_dt + exp_fac * int_proj

      gnew(1, :, :, :, :) = gnew(1, :, :, :, :) - code_dt * g

      deallocate (g)

      if (proc0) call time_message(.false., time_sources(:, 1), ' sources')

   end subroutine project_out_zero

   subroutine init_quasineutrality_source

#ifdef ISO_C_BINDING
      use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_intptr_t
      use arrays_store_useful, only: qn_zf_window
      use mp, only: sgproc0, sharedsubprocs, comm_sgroup
      use mp, only: real_size, nbytes_real, create_shared_memory_window
      use mp_lu_decomposition, only: lu_decomposition_local, lu_inverse_local
      use mpi
#endif
      use geometry, only: dl_over_b, d_dl_over_b_drho
      use calculations_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use grids_z, only: nzgrid, nztot
      use parameters_kxky_grid, only: nakx
      use grids_kxky, only: rho_d_clamped
      use parameters_multibox, only: boundary_size
      use linear_solve, only: lu_decomposition
      use arrays_store_fields, only: phizf_solve, phi_ext
      use arrays_store_useful, only: c_mat, theta
      use arrays_store_useful, only: tcorr_source_qn, exclude_boundary_regions_qn, exp_fac_qn

      implicit none

      integer :: iz, ikx, ia, jkx, jz
      integer :: inmat, jnmat, nmat_zf
      real :: dum
#ifdef ISO_C_BINDING
      integer :: ierr, temp_window
      integer(c_intptr_t):: cur_pos
      integer(kind=MPI_ADDRESS_KIND) :: win_size
      type(c_ptr) :: cptr
      complex, dimension(:, :), pointer :: temp_mat
#endif
      complex, dimension(:, :), allocatable :: g0k, g0x, g1k

      ia = 1

      if (qn_source_initialized) return
      qn_source_initialized = .true.

      if (include_qn_source) then
         nmat_zf = nakx * (nztot - 1)
#ifdef ISO_C_BINDING
         if (qn_zf_window == MPI_WIN_NULL) then
            win_size = 0
            if (sgproc0) then
               win_size = int(nmat_zf, MPI_ADDRESS_KIND) * 4_MPI_ADDRESS_KIND &
                          + int(nmat_zf * (nmat_zf + 1), MPI_ADDRESS_KIND) * 2 * real_size !complex size
            end if

            if (debug) write (6, *) 'sources::init_quasineutrality_source::win_allocate'

            call create_shared_memory_window(win_size, qn_zf_window, cur_pos)

            !allocate the memory
            if (.not. associated(phizf_solve%zloc)) then
               cptr = transfer(cur_pos, cptr)
               call c_f_pointer(cptr, phizf_solve%zloc, (/nmat_zf, nmat_zf/))
            end if
            cur_pos = cur_pos + nmat_zf**2 * 2 * nbytes_real

            if (.not. associated(phi_ext)) then
               cptr = transfer(cur_pos, cptr)
               call c_f_pointer(cptr, phi_ext, (/nmat_zf/))
            end if
            cur_pos = cur_pos + nmat_zf * 2 * nbytes_real

            if (.not. associated(phizf_solve%idx)) then
               cptr = transfer(cur_pos, cptr)
               call c_f_pointer(cptr, phizf_solve%idx, (/nmat_zf/))
            end if

            call mpi_win_fence(0, qn_zf_window, ierr)
         end if
#else
         if (debug) write (6, *) 'sources::init_quasineutrality_source::allocate_phizf'
         if (.not. associated(phizf_solve%zloc)) allocate (phizf_solve%zloc(nmat_zf, nmat_zf))
         if (.not. associated(phizf_solve%idx)) allocate (phizf_solve%idx(nmat_zf))
         if (.not. associated(phi_ext)) allocate (phi_ext(nmat_zf))
#endif

#ifdef ISO_C_BINDING
         if (sgproc0) then
#endif
            allocate (g0k(1, nakx))
            allocate (g0x(1, nakx))
            allocate (g1k(1, nakx))

            phizf_solve%zloc = 0.

            !get the big matrix
            do jz = -nzgrid, nzgrid - 1
               do jkx = 1, nakx
                  jnmat = jkx + nakx * (jz + nzgrid)
                  ! C.phi
                  do ikx = 1, nakx
                     inmat = ikx + nakx * (jz + nzgrid)
                     phizf_solve%zloc(inmat, jnmat) = phizf_solve%zloc(inmat, jnmat) + c_mat(ikx, jkx)
                  end do

                  ! -C.<phi>_\psi
                  g0k = 0.0; g0k(1, jkx) = 1.0
                  call transform_kx2x_unpadded(g0k, g0x)
                  g0x(1, :) = (dl_over_b(ia, jz) + d_dl_over_b_drho(ia, jz) * rho_d_clamped) * g0x(1, :)
                  call transform_x2kx_unpadded(g0x, g0k)

                  !set the gauge potential
                  if (jkx == 1) g0k(1, 1) = 0.

                  do ikx = 1, nakx
                     g1k(1, ikx) = sum(c_mat(ikx, :) * g0k(1, :))
                  end do

                  do iz = -nzgrid, nzgrid - 1
                     do ikx = 1, nakx
                        inmat = ikx + nakx * (iz + nzgrid)
                        phizf_solve%zloc(inmat, jnmat) = phizf_solve%zloc(inmat, jnmat) - g1k(1, ikx)
                     end do
                  end do

                  ! get theta.phi
                  g1k(1, :) = theta(:, jkx, jz)

                  ! +theta.phi
                  do ikx = 1, nakx
                     inmat = ikx + nakx * (jz + nzgrid)
                     phizf_solve%zloc(inmat, jnmat) = phizf_solve%zloc(inmat, jnmat) + g1k(1, ikx)
                  end do

                  ! -<<theta.phi>_psi>_T
                  call transform_kx2x_unpadded(g1k, g0x)
                  g0x(1, :) = (dl_over_b(ia, jz) + d_dl_over_b_drho(ia, jz) * rho_d_clamped) * g0x(1, :)

                  if (exclude_boundary_regions_qn) then
                     g0x(1, :) = sum(g0x(1, (boundary_size + 1):(nakx - boundary_size))) &
                                 / (nakx - 2 * boundary_size)
                     g0x(1, 1:boundary_size) = 0.0
                     g0x(1, (nakx - boundary_size + 1):) = 0.0
                  else
                     g0x(1, :) = sum(g0x(1, :)) / nakx
                  end if

                  call transform_x2kx_unpadded(g0x, g0k)

                  if (tcorr_source_qn > epsilon(0.)) then
                     g0k = (1.-exp_fac_qn) * g0k
                  end if

                  do iz = -nzgrid, nzgrid - 1
                     do ikx = 1, nakx
                        inmat = ikx + nakx * (iz + nzgrid)
                        phizf_solve%zloc(inmat, jnmat) = phizf_solve%zloc(inmat, jnmat) &
                                                         - g0k(1, ikx)
                     end do
                  end do
               end do
            end do
            deallocate (g0k, g1k, g0x)
#ifdef ISO_C_BINDING
         end if

         call mpi_win_fence(0, qn_zf_window, ierr)
         if (debug) write (6, *) 'sources::init_quasineutrality_source::lu_decomposition'
         call lu_decomposition_local(comm_sgroup, 0, qn_zf_window, &
                                     phizf_solve%zloc, phizf_solve%idx, dum)

         call mpi_win_fence(0, qn_zf_window, ierr)

         if (debug) write (6, *) 'sources::init_quasineutrality_source::temp_mat'

         win_size = 0
         if (sgproc0) then
            win_size = int(nmat_zf**2, MPI_ADDRESS_KIND) * 2 * real_size !complex size
         end if

         if (debug) write (6, *) 'sources::init_quasineutrality_source::win_allocate'
         call create_shared_memory_window(win_size, temp_window, cur_pos)

         cptr = transfer(cur_pos, cptr)

         call c_f_pointer(cptr, temp_mat, (/nmat_zf, nmat_zf/))

         if (sgproc0) temp_mat = phizf_solve%zloc

         call mpi_win_fence(0, temp_window, ierr)

         ! inverse is calculated since it is more straightforward to parallelize
         ! inverse calculation/matrix multiplication than the lu back substitution
         if (debug) write (6, *) 'sources::init_quasineutrality_source::lu_inverse'
         call lu_inverse_local(comm_sgroup, qn_zf_window, &
                               temp_mat, phizf_solve%idx, phizf_solve%zloc)
         call mpi_win_free(temp_window, ierr)
#else
         if (debug) write (6, *) 'sources::init_quasineutrality_source::lu_decomposition'
         call lu_decomposition(phizf_solve%zloc, phizf_solve%idx, dum)
#endif
      end if
   end subroutine init_quasineutrality_source

   subroutine update_quasineutrality_source

      use arrays_store_fields, only: phi_proj, phi_proj_stage
      use arrays_store_useful, only: tcorr_source_qn, exp_fac_qn

      implicit none

      if (tcorr_source_qn < epsilon(0.)) then
         phi_proj = phi_proj_stage
      else
         phi_proj = exp_fac_qn * phi_proj + (1.-exp_fac_qn) * phi_proj_stage
      end if

   end subroutine update_quasineutrality_source

end module gk_sources
