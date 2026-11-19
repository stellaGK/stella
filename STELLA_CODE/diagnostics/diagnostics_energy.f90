!###############################################################################
!                           ENERGY DIAGNOSTICS
!###############################################################################
!
! This module implements a free-energy closure diagnostic for the gyrokinetic
! equation (GKE): the free-energy change dH/dt is obtained as the projection of
! the GKE RHS onto the free-energy kernel at the single time t_n, and the
! individual RHS terms are projected separately. The exact two-level free-energy
! change is kept as a SEPARATE column, deliberately not claimed to equal the
! projection.
!
! Closure layout
! --------------
! The total free energy is H = H_p + H_chi. The diagnostic writes, per saved step:
!   col  1  time
!   col  2  particle energy H_p
!   col  3  field energy H_chi (electrostatic; EM deferred)
!   col  4  total energy H = H_p + H_chi
!   col  5  dedt_particle = (H_p(t_n) - H_p(t_{n-1}))/(t_n - t_{n-1})
!   col  6  dedt_field    = (H_chi(t_n) - H_chi(t_{n-1}))/(t_n - t_{n-1})
!   col  7  dedt_total    = dedt_particle + dedt_field
!   col  8  dedt_projected = <c_f | RHS(t_n)>  (projection identity, at t_n)
!   col  9  drive (wstar / background-gradient source)
!   col 10  heat sources (nonlocal) -> 0 in the local flux-tube limit
!   col 11  collisions
!   col 12  hyp_v (parallel-velocity hyperdiffusion = diss_vpa)
!   col 13  hyp_z (z hyperdiffusion = diss_zed)
!   col 14  hyp_x/y (perpendicular hyperdiffusion = diss_perp)
!   col 15  nonlinearity (ExB)
!   col 16  parallel streaming
!   col 17  mirror
!   col 18  drifts (drift/drifts term; makes the projection identity close by
!           construction)
!
! The projection identity col 8 = sum(cols 9..18) holds to round-off BY
! CONSTRUCTION because col 18 is defined as the remainder. The physically
! meaningful outputs are the individual term projections (cols 9..18) and the
! exact change (col 7), which is deliberately different from col 8 for
! fast-evolving fields.
!
! Normalisation
! -------------
! The particle energy is
!     H_p = 0.5 * sum_s < h_s * h_s* / F_{0,s} >
! and the field energy is
!     H_chi = sum_k 0.5 * mode_fac(ky) * dVolume * lambda(k) * |phi_k|^2
! with lambda = denominator_fields (the quasineutrality coefficient). All
! scalars are volume averages over the (kx,zed) grid. The drive term is
! projected onto the same h* kernel as all other terms.
!
! Per-mode array normalisation
! ----------------------------
! This is the single, normative convention for EVERY <*_kxkyz> energy array
! (particle energy, field energy, and each per-term array):
!
!     array  = 0.5 * mode_fac(ky) * <integrand>          (NO dVolume inside)
!     scalar = sum_{kx,ky,zed,tube} array * dVolume(kx,zed) / volume_total
!
! with volume_total = sum_{kx,zed} dVolume in a flux tube. This mirrors
! qflux_vs_kxkyz: the volume weight is applied by the consumer, not stored.
! The vmu arrays are already volume-averaged and reproduce the scalar through
!     scalar = sum_{vpa,mu} wgts_vpa * wgts_mu_bare * array,
! where wgts_vpa is the CONSERVATIVE parallel-velocity weight dvpa/sqrt(pi)
! that this diagnostic switches to internally. That weight is written to the
! NetCDF file as <wgts_vpa_energy>; the stored <wgts_vpa> is the default
! (Simpson) weight and does NOT reproduce the scalars.
! field_energy_vmu / field_energy_em_vmu are constant over (vpa,mu) and equal
! to the scalar. They are a plotting convenience, NOT an integral density: a
! velocity-weighted sum returns the scalar times
! sum(wgts_vpa)*sum(wgts_mu_bare), which is not 1.
!
! Full flux surface: the FFS path is a DIFFERENT representation. There the
! volume element depends on the field-line label alpha, so the weight is not
! diagonal in ky and is currently folded into the per-mode arrays (volfac in
! ffs_real_to_ky); the FFS scalar is then the plain sum of the per-mode array.
! The FFS scalar normalisation is known to be inconsistent with the flux tube
! and is deferred (see TODO.md). The representation in use is recorded in the
! 'energy_kxkyz_convention' NetCDF attribute.
!
! Output representations
! -----------------------
! The diagnostic is written out in two representations:
!   - "kxkyz": on the full (ky,kx,zed,tube,species) grid, and
!   - "vmu":  spatially averaged, on the (vpa,mu,species) grid.
! Both are written to NetCDF when <write_energy_kxkyz> / <write_energy_vmu>.
! The combined (sum-over-species) scalars are written to an ascii file
! (<stem>.energy, 18 columns). Per-species files (<stem>.energy_species_<label>)
! are written when <write_energy_species> = .true. (default off).
!
! Electrostatic only: electromagnetic field energy and Poynting/work terms are
! NOT included (deferred).
!###############################################################################
module diagnostics_energy

   implicit none

   public :: energy_diagnostics_type, energy_diag
   public :: init_energy_diagnostic, finish_energy_diagnostic
   public :: get_free_energy

   type :: energy_diagnostics_type

      !> Closure scalars (combined over all species).
      !> - particle_energy_sum: H_p at t_n (particle energy)
      !> - field_energy_sum:    H_chi at t_n (electrostatic + electromagnetic field energy)
      !> - total_energy_sum:    H_p + H_chi (col 3; used for dedt_total)
      !> - dedt_proj_sum:       <c_f|RHS> at t_n
      !> - drive/coll/diss_vpa/diss_zed/diss_perp/nonlinearity/streaming/mirror/
      !>   drifts:           the per-term projections
      real :: particle_energy_sum
      real :: field_energy_sum
      real :: total_energy_sum
      real :: dedt_proj_sum
      real :: dedt_particle_sum, dedt_field_sum, dedt_total_sum
      real :: drive_sum, coll_sum, diss_vpa_sum, diss_zed_sum, diss_perp_sum
      real :: nonlinearity_sum, streaming_sum, mirror_sum, drifts_sum
      real :: last_particle_energy, last_field_energy, last_total_energy, last_time

      !> Per-species closure scalars (nspec). Used for the per-species ascii
      !> files and the combined sums.
      real, dimension(:), allocatable :: energy_spec, dedt_proj_spec
      real, dimension(:), allocatable :: drive_spec, coll_spec
      real, dimension(:), allocatable :: diss_vpa_spec, diss_zed_spec, diss_perp_spec
      real, dimension(:), allocatable :: nonlinearity_spec, streaming_spec, mirror_spec, drifts_spec

      !> Per-species prefactor (dens*temp, or -D_hyper*dens*temp for diss_perp)
      real, dimension(:), allocatable :: factor_spec

      !> Work arrays for the velocity/spatial integrals in the projections.
      complex, dimension(:, :, :, :, :), allocatable :: velocity_integral1
      complex, dimension(:, :, :), allocatable :: spatial_integral1

      !> Per-mode (kxkyz representation) closure terms,
      !> (ky,kx,z,tube,species). Written to NetCDF when write_energy_kxkyz.
      real, dimension(:, :, :, :, :), allocatable :: particle_energy_kxkyz
      real, dimension(:, :, :, :, :), allocatable :: dedt_proj_kxkyz
      real, dimension(:, :, :, :, :), allocatable :: drive_kxkyz, coll_kxkyz
      real, dimension(:, :, :, :, :), allocatable :: diss_vpa_kxkyz, diss_zed_kxkyz
      real, dimension(:, :, :, :, :), allocatable :: diss_perp_kxkyz
      real, dimension(:, :, :, :, :), allocatable :: nonlinearity_kxkyz
      real, dimension(:, :, :, :, :), allocatable :: streaming_kxkyz, mirror_kxkyz, drifts_kxkyz
      real, dimension(:, :, :, :, :), allocatable :: field_energy_kxkyz

      !> Real-space energy densities (unfolded: no dVolume/volume_total and no
      !> 1/(nx*ny) FFT factor). Written to NetCDF when write_energy_xy.
      !> field: (x,y,z,tube); particle: (x,y,z,tube,species).
      real, dimension(:, :, :, :), allocatable :: field_energy_xy_z
      real, dimension(:, :, :, :, :), allocatable :: particle_energy_xy_z

      !> Per-mode (vmu representation) closure terms,
      !> (vpa,mu,species). Written to NetCDF when write_energy_vmu.
      real, dimension(:, :, :), allocatable :: particle_energy_vmu
      real, dimension(:, :, :), allocatable :: dedt_proj_vmu
      real, dimension(:, :, :), allocatable :: drive_vmu, coll_vmu
      real, dimension(:, :, :), allocatable :: diss_vpa_vmu, diss_zed_vmu
      real, dimension(:, :, :), allocatable :: diss_perp_vmu
      real, dimension(:, :, :), allocatable :: nonlinearity_vmu
      real, dimension(:, :, :), allocatable :: streaming_vmu, mirror_vmu, drifts_vmu
      real, dimension(:, :, :), allocatable :: field_energy_vmu
      real, dimension(:, :, :), allocatable :: field_energy_vmu_dummy

      !> Per-species ascii '.energy_species_<label>' files (optional).
      integer, dimension(:), allocatable :: energy_species_units
      character(len=20), dimension(:), allocatable :: energy_species_labels
   end type energy_diagnostics_type

   type(energy_diagnostics_type), save :: energy_diag

   logical, parameter :: debug = .false.

contains

   !> Initialise the free-energy diagnostic. Allocates all per-species scalar
   !> vectors, the work arrays, and (if requested) the per-mode kxkyz/vmu arrays.
   !> Opens the per-species ascii files when write_energy_species = .true.
   subroutine init_energy_diagnostic(write_energy_vmu, write_energy_species, write_energy_xy, restart)

      use mp, only: proc0
      use grids_species, only: nspec, spec, ion_species, electron_species
      use grids_kxky, only: nakx, naky, nx, ny, grid_option_switch, grid_option_box
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: nmu, nvpa
      use parameters_physics, only: radial_variation
      use file_utils, only: open_output_file, flush_output_file

      implicit none

      logical, intent(in) :: write_energy_vmu
      logical, intent(in) :: write_energy_species
      logical, intent(in) :: write_energy_xy
      logical, intent(in) :: restart

      integer :: is, ion_count, elec_count

      if (write_energy_xy .and. (grid_option_switch /= grid_option_box .or. radial_variation)) then
         if (proc0) write (*, *) 'ERROR: write_energy_xy requires grid_option = box and radial_variation = .false.'
         error stop
      end if

      if (write_energy_species .and. proc0) then
         if (.not. allocated(energy_diag%energy_species_units)) allocate (energy_diag%energy_species_units(nspec))
         if (.not. allocated(energy_diag%energy_species_labels)) allocate (energy_diag%energy_species_labels(nspec))
         ion_count = 0
         elec_count = 0
         do is = 1, nspec
            if (spec(is)%type == electron_species) then
               elec_count = elec_count + 1
               energy_diag%energy_species_labels(is) = 'electron_' // trim(achar(48 + min(elec_count, 9)))
            else
               ion_count = ion_count + 1
               energy_diag%energy_species_labels(is) = 'ion_' // trim(achar(48 + min(ion_count, 9)))
            end if
            call open_output_file(energy_diag%energy_species_units(is), '.energy_species_' // energy_diag%energy_species_labels(is), .not. restart)
            write (energy_diag%energy_species_units(is), '(A)') '# 12 columns: time, energy, dedt_proj,'
            write (energy_diag%energy_species_units(is), '(A)') '#            drive, collisions, diss_vpa, diss_zed, diss_perp,'
            write (energy_diag%energy_species_units(is), '(A)') '#            nonlinearity, streaming, mirror, drifts'
            call flush_output_file(energy_diag%energy_species_units(is))
         end do
      end if

      if (.not. allocated(energy_diag%energy_spec)) allocate (energy_diag%energy_spec(nspec)); energy_diag%energy_spec = 0.
      if (.not. allocated(energy_diag%dedt_proj_spec)) allocate (energy_diag%dedt_proj_spec(nspec)); energy_diag%dedt_proj_spec = 0.
      if (.not. allocated(energy_diag%drive_spec)) allocate (energy_diag%drive_spec(nspec)); energy_diag%drive_spec = 0.
      if (.not. allocated(energy_diag%coll_spec)) allocate (energy_diag%coll_spec(nspec)); energy_diag%coll_spec = 0.
      if (.not. allocated(energy_diag%diss_vpa_spec)) allocate (energy_diag%diss_vpa_spec(nspec)); energy_diag%diss_vpa_spec = 0.
      if (.not. allocated(energy_diag%diss_zed_spec)) allocate (energy_diag%diss_zed_spec(nspec)); energy_diag%diss_zed_spec = 0.
      if (.not. allocated(energy_diag%diss_perp_spec)) allocate (energy_diag%diss_perp_spec(nspec)); energy_diag%diss_perp_spec = 0.
      if (.not. allocated(energy_diag%nonlinearity_spec)) allocate (energy_diag%nonlinearity_spec(nspec)); energy_diag%nonlinearity_spec = 0.
      if (.not. allocated(energy_diag%streaming_spec)) allocate (energy_diag%streaming_spec(nspec)); energy_diag%streaming_spec = 0.
      if (.not. allocated(energy_diag%mirror_spec)) allocate (energy_diag%mirror_spec(nspec)); energy_diag%mirror_spec = 0.
      if (.not. allocated(energy_diag%drifts_spec)) allocate (energy_diag%drifts_spec(nspec)); energy_diag%drifts_spec = 0.

      if (.not. allocated(energy_diag%factor_spec)) allocate (energy_diag%factor_spec(nspec)); energy_diag%factor_spec = 0.

      if (.not. allocated(energy_diag%velocity_integral1)) then
         allocate (energy_diag%velocity_integral1(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
         energy_diag%velocity_integral1 = 0.
      end if
      if (.not. allocated(energy_diag%spatial_integral1)) then
         allocate (energy_diag%spatial_integral1(nvpa, nmu, nspec))
         energy_diag%spatial_integral1 = 0.
      end if

      ! The kxkyz-representation per-mode arrays are always allocated (the scalar
      ! energies are obtained by integrating them); writing them to NetCDF is
      ! controlled by <write_energy_kxkyz> in the diagnostics module.
      if (.not. allocated(energy_diag%particle_energy_kxkyz)) allocate (energy_diag%particle_energy_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      if (.not. allocated(energy_diag%dedt_proj_kxkyz)) allocate (energy_diag%dedt_proj_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      if (.not. allocated(energy_diag%drive_kxkyz)) allocate (energy_diag%drive_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      if (.not. allocated(energy_diag%coll_kxkyz)) allocate (energy_diag%coll_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      if (.not. allocated(energy_diag%diss_vpa_kxkyz)) allocate (energy_diag%diss_vpa_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      if (.not. allocated(energy_diag%diss_zed_kxkyz)) allocate (energy_diag%diss_zed_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      if (.not. allocated(energy_diag%diss_perp_kxkyz)) allocate (energy_diag%diss_perp_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      if (.not. allocated(energy_diag%nonlinearity_kxkyz)) allocate (energy_diag%nonlinearity_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      if (.not. allocated(energy_diag%streaming_kxkyz)) allocate (energy_diag%streaming_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      if (.not. allocated(energy_diag%mirror_kxkyz)) allocate (energy_diag%mirror_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      if (.not. allocated(energy_diag%drifts_kxkyz)) allocate (energy_diag%drifts_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      if (.not. allocated(energy_diag%field_energy_kxkyz)) allocate (energy_diag%field_energy_kxkyz(naky, nakx, -nzgrid:nzgrid, ntubes, 1))

      ! The real-space arrays are only needed when write_energy_xy is on.
      if (write_energy_xy) then
         if (.not. allocated(energy_diag%field_energy_xy_z)) then
            allocate (energy_diag%field_energy_xy_z(nx, ny, -nzgrid:nzgrid, ntubes))
            energy_diag%field_energy_xy_z = 0.
         end if
         if (.not. allocated(energy_diag%particle_energy_xy_z)) then
            allocate (energy_diag%particle_energy_xy_z(nx, ny, -nzgrid:nzgrid, ntubes, nspec))
            energy_diag%particle_energy_xy_z = 0.
         end if
      end if

      if (write_energy_vmu) then
         if (.not. allocated(energy_diag%particle_energy_vmu)) allocate (energy_diag%particle_energy_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%dedt_proj_vmu)) allocate (energy_diag%dedt_proj_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%drive_vmu)) allocate (energy_diag%drive_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%coll_vmu)) allocate (energy_diag%coll_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%diss_vpa_vmu)) allocate (energy_diag%diss_vpa_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%diss_zed_vmu)) allocate (energy_diag%diss_zed_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%diss_perp_vmu)) allocate (energy_diag%diss_perp_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%nonlinearity_vmu)) allocate (energy_diag%nonlinearity_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%streaming_vmu)) allocate (energy_diag%streaming_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%mirror_vmu)) allocate (energy_diag%mirror_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%drifts_vmu)) allocate (energy_diag%drifts_vmu(nvpa, nmu, nspec))
         if (.not. allocated(energy_diag%field_energy_vmu)) allocate (energy_diag%field_energy_vmu(nvpa, nmu, 1))
         if (.not. allocated(energy_diag%field_energy_vmu_dummy)) allocate (energy_diag%field_energy_vmu_dummy(nvpa, nmu, 1))
      end if

      energy_diag%particle_energy_sum = 0.
      energy_diag%field_energy_sum = 0.
      energy_diag%total_energy_sum = 0.
      energy_diag%dedt_proj_sum = 0.
      energy_diag%dedt_particle_sum = 0.
      energy_diag%dedt_field_sum = 0.
      energy_diag%dedt_total_sum = 0.
      energy_diag%drive_sum = 0.
      energy_diag%coll_sum = 0.
      energy_diag%diss_vpa_sum = 0.
      energy_diag%diss_zed_sum = 0.
      energy_diag%diss_perp_sum = 0.
      energy_diag%nonlinearity_sum = 0.
      energy_diag%streaming_sum = 0.
      energy_diag%mirror_sum = 0.
      energy_diag%drifts_sum = 0.
      energy_diag%last_particle_energy = 0.
      energy_diag%last_field_energy = 0.
      energy_diag%last_total_energy = 0.
      energy_diag%last_time = -1.0

   end subroutine init_energy_diagnostic

   !> Compute one GKE RHS energy term on the (ky,kx,zed,tube,species) grid.
   !>
   !> Projects < term * h* / F0 >, integrated over velocity space. The per-mode
   !> array (term_kxkyz) stores 0.5*mode_fac(ky)*Re[factor*velocity_integral]
   !> WITHOUT the dVolume weight; the per-species scalars (sum_spec) apply
   !> dVolume/volume_total. See the module header for the convention.
   subroutine get_one_energy_term_kxkyz(h, term, factor_spec, sum_spec, sum_total, term_kxkyz)
      use mp, only: proc0, sum_allreduce
      use arrays_distribution_function, only: g0
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      use grids_species, only: nspec
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use grids_velocity, only: wgts_mu, wgts_vpa
      use calculations_volume_averages, only: mode_fac
      use grids_kxky, only: naky, nakx
      use geometry, only: dVolume
      use calculations_volume_averages, only: volume_total

      implicit none
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: h, term

      real, dimension(nspec), intent(in) :: factor_spec
      real, dimension(nspec), intent(out) :: sum_spec
      real, intent(out) :: sum_total
      real, dimension(:, :, -nzgrid:, :, :), intent(out) :: term_kxkyz

      integer :: ivmu, imu, iv, iz, it, is, ia, ikx, iky

      sum_spec = 0.
      term_kxkyz = 0.
      energy_diag%velocity_integral1 = 0.
      sum_total = 0.
      g0 = 0.

      ia = 1
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do iz = -nzgrid, nzgrid
            g0(:, :, iz, :, ivmu) = term(:, :, iz, :, ivmu) * conjg(h(:, :, iz, :, ivmu)) &
                                    / (maxwell_fac(is) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is))
            energy_diag%velocity_integral1(:, :, iz, :, is) = energy_diag%velocity_integral1(:, :, iz, :, is) + &
                                                          wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g0(:, :, iz, :, ivmu)
         end do
      end do

      call sum_allreduce(energy_diag%velocity_integral1)

      if (proc0) then
         do is = 1, nspec
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do ikx = 1, nakx
                     do iky = 1, naky
                        ! 0.5*mode_fac(ky): count the +-ky pair once (mode_fac=2
                        ! for ky!=0 gives 0.5*2=1; mode_fac=1 for ky=0 gives 0.5).
                        term_kxkyz(iky, ikx, iz, it, is) = 0.5 * mode_fac(iky) * (real(factor_spec(is) * energy_diag%velocity_integral1(iky, ikx, iz, it, is)))
                        sum_spec(is) = sum_spec(is) + term_kxkyz(iky, ikx, iz, it, is) * dVolume(ia, ikx, iz)
                     end do
                  end do
               end do
            end do
            sum_spec(is) = sum_spec(is) / volume_total
            sum_total = sum_total + sum_spec(is)
         end do
      end if
      g0 = 0.

   end subroutine get_one_energy_term_kxkyz

   !> Compute one GKE RHS energy term on the (vpa,mu,species) grid.
   !> Same as get_one_energy_term_kxkyz but the spatial integral is done in the
   !> vmu representation (scatter to vmu, integrate over (vpa,mu) with the bmag
   !> factor from the 2*pi*B_0/m_nu Jacobian).
   subroutine get_one_energy_term_vmu(h, term, factor_spec, term_vmu)
      use mp, only: proc0, sum_allreduce
      use arrays_distribution_function, only: g0, gvmu0
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: is_idx, ikx_idx, iky_idx, iz_idx, it_idx
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      use grids_species, only: nspec
      use grids_z, only: ntubes, nzgrid
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use grids_velocity, only: nvpa, nmu

      use calculations_volume_averages, only: mode_fac
      use geometry, only: dVolume, bmag
      use calculations_volume_averages, only: volume_total

      use redistribute, only: scatter
      use initialise_redistribute, only: kxkyz2vmu

      implicit none
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: h, term

      real, dimension(nspec), intent(in) :: factor_spec

      real, dimension(:, :, :), intent(out) :: term_vmu

      integer :: ivmu, imu, iv, is
      integer :: ikxkyz, iz, it, ia, ikx, iky

      term_vmu = 0.
      energy_diag%spatial_integral1 = 0.
      g0 = 0.
      gvmu0 = 0.

      ia = 1

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0(:, :, iz, it, ivmu) = term(:, :, iz, it, ivmu) * conjg(h(:, :, iz, it, ivmu))
            end do
         end do
      end do

      call scatter(kxkyz2vmu, g0, gvmu0)

      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         is = is_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iky = iky_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            do iv = 1, nvpa
               ! The vmu integral carries the explicit bmag (= B_0) factor from the
               ! 2*pi*B_0/m_nu Jacobian; the kxkyz/ASCII form instead folds B_0 into
               ! the velocity weight wgts_mu = 2*bmag*wgts_mu_bare. Both give the same
               ! total H_p: 0.5*mode_fac*2*bmag (kxkyz, with the 0.5 from the +-ky
               ! convention) equals mode_fac*bmag (vmu).
               energy_diag%spatial_integral1(iv, imu, is) = energy_diag%spatial_integral1(iv, imu, is) + mode_fac(iky) * bmag(ia, iz) * &
               factor_spec(is) * gvmu0(iv, imu, ikxkyz) * dVolume(ia, ikx, iz) / (maxwell_fac(is) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is))
            end do
         end do
      end do

      energy_diag%spatial_integral1 = energy_diag%spatial_integral1 / volume_total
      call sum_allreduce(energy_diag%spatial_integral1)

      if (proc0) then
         do is = 1, nspec
            do imu = 1, nmu
               do iv = 1, nvpa
                  term_vmu(iv, imu, is) = real(energy_diag%spatial_integral1(iv, imu, is))
               end do
            end do
         end do
      end if

      g0 = 0.
      gvmu0 = 0.

   end subroutine get_one_energy_term_vmu

   !> Transform a (ky,kx) Fourier slice to real space (y,x) using the same
   !> unnormalised FFTW transforms as the rest of the code. The output is
   !> Nx*Ny times the normalised inverse transform; callers that need the
   !> volume-average scalar must apply 1/(nx*ny) themselves.
   subroutine transform_kxky_to_xy(f_kxky, f_xy)

      use calculations_kxky, only: swap_kxky
      use calculations_transforms, only: transform_ky2y, transform_kx2x
      use grids_kxky, only: naky, nakx, naky_all, ikx_max, nx, ny

      implicit none

      complex, dimension(naky, nakx), intent(in) :: f_kxky
      real, dimension(ny, nx), intent(out) :: f_xy

      complex, dimension(naky_all, ikx_max) :: f_swap
      complex, dimension(ny, ikx_max) :: f_ykx

      call swap_kxky(f_kxky, f_swap)
      call transform_ky2y(f_swap, f_ykx)
      call transform_kx2x(f_ykx, f_xy)

   end subroutine transform_kxky_to_xy

   !> Calculate the electrostatic field energy H_chi = sum_k 0.5*mode_fac*dVolume*lambda*|phi_k|^2.
   !> lambda = denominator_fields is the k-dependent quasineutrality coefficient
   !> (1-Gamma_0 for kinetic species; already includes the Boltzmann efac for
   !> adiabatic electrons).
   !>
   !> field_energy_kxkyz stores the per-mode value WITHOUT the dVolume weight;
   !> the scalar applies dVolume/volume_total (see the module header).
   subroutine get_field_energy(phi, field_energy_sum, field_energy_kxkyz, write_vmu, field_energy_vmu, &
                               write_energy_xy, field_energy_xy_z)

      use mp, only: proc0
      use arrays, only: denominator_fields
      use grids_kxky, only: naky, nakx, nx, ny, naky_all, ikx_max
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: nvpa, nmu
      use calculations_volume_averages, only: mode_fac
      use geometry, only: dVolume
      use calculations_volume_averages, only: volume_total

      implicit none
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      real, intent(out) :: field_energy_sum
      real, dimension(:, :, -nzgrid:, :, :), intent(out) :: field_energy_kxkyz
      logical, intent(in) :: write_vmu
      real, dimension(:, :, :), intent(out) :: field_energy_vmu
      logical, intent(in) :: write_energy_xy
      real, dimension(:, :, -nzgrid:, :), intent(out), optional :: field_energy_xy_z

      integer :: ikx, iky, iz, it, ia
      real :: tmp
      complex, dimension(naky, nakx) :: lam_phi
      real, dimension(ny, nx) :: phi_xy, lam_phi_xy

      ia = 1
      field_energy_sum = 0.
      field_energy_kxkyz = 0.
      if (write_vmu) field_energy_vmu = 0.

      do ikx = 1, nakx
         do iky = 1, naky
            do iz = -nzgrid, nzgrid
               do it = 1, ntubes
                  ! Per-mode value is stored WITHOUT the volume weight; the
                  ! weight is applied only to the scalar below.
                  tmp = 0.5 * mode_fac(iky) &
                        * denominator_fields(iky, ikx, iz) &
                        * (abs(phi(iky, ikx, iz, it))**2)
                  field_energy_kxkyz(iky, ikx, iz, it, 1) = tmp
                  field_energy_sum = field_energy_sum + tmp * dVolume(ia, ikx, iz)
               end do
            end do
         end do
      end do

      field_energy_sum = field_energy_sum / volume_total

      if (write_energy_xy .and. proc0 .and. present(field_energy_xy_z)) then
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               lam_phi = denominator_fields(:, :, iz) * phi(:, :, iz, it)
               call transform_kxky_to_xy(phi(:, :, iz, it), phi_xy)
               call transform_kxky_to_xy(lam_phi, lam_phi_xy)
               ! Stored unfolded: no dVolume/volume_total and no 1/(nx*ny).
               field_energy_xy_z(:, :, iz, it) = transpose(0.5 * phi_xy * lam_phi_xy)
            end do
         end do
      end if

      if (write_vmu) field_energy_vmu(:, :, 1) = field_energy_sum

   end subroutine get_field_energy

   !> Calculate the full free-energy closure at the current time step.
   !>
   !> Arguments:
   !>   h            - non-adiabatic response at t_n (ky,kx,z,tube,vmu) [in]
   !>   g            - full distribution at t_n (the pristine pre-conversion copy) [in]
   !>   phi          - electrostatic potential at t_n [in]
   !>   istep        - current time step (for the ExB nonlinearity restart flag)
   !>   energy_unit  - ascii unit number for the combined .energy file
   !>   write_energy_vmu - if .true., also compute the vmu-representation terms
   subroutine get_free_energy(h, g, phi, istep, energy_unit, write_energy_vmu, write_energy_xy)

      use file_utils, only: flush_output_file
      use mp, only: proc0, sum_allreduce
      use arrays_distribution_function, only: g1
      use arrays, only: kperp2
      use arrays_fields, only: mod_bpar => bpar
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      use grids_species, only: spec, nspec
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: nmu, nvpa
      use grids_velocity, only: wgts_vpa, wgts_mu, maxwell_fac, maxwell_vpa, maxwell_mu
      use parameters_physics, only: fphi
      use grids_kxky, only: naky, nakx, nx, ny, naky_all, ikx_max
      use grids_time, only: code_time, code_dt
      use parameters_physics, only: include_nonlinear
      use dissipation_hyper, only: D_hyper, k2max, advance_hyper_zed, advance_hyper_vpa, hyp_vpa, hyp_zed
      use gk_magnetic_drift, only: advance_wdriftx_explicit, advance_wdrifty_explicit
      use gk_drive, only: advance_wstar_explicit
      use gk_parallel_streaming, only: advance_parallel_streaming_explicit
      use gk_mirror, only: advance_mirror_explicit
      use gk_ExB_nonlinearity, only: advance_ExB_nonlinearity
      use dissipation_and_collisions, only: advance_collisions_explicit, include_collisions
      use parameters_physics, only: include_mirror
      use parameters_physics, only: include_parallel_streaming
      use parameters_physics, only: include_drive
      use parameters_physics, only: xdriftknob, ydriftknob
      use job_manage, only: time_message
      use timers, only: time_energy_terms
      use debug_flags, only: energy_diagnostics_debug

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(inout) :: g, h
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi

      integer, intent(in) :: istep
      integer, intent(in) :: energy_unit
      logical, intent(in) :: write_energy_vmu
      logical, intent(in) :: write_energy_xy

      integer :: ivmu, imu, iv, iz, it, is, ia, ikx, iky
      real :: particle_energy
      real :: dedt_proj_sum
      real :: drive_sum, coll_sum, diss_vpa_sum, diss_zed_sum, diss_perp_sum
      real :: nonlinearity_sum, streaming_sum, mirror_sum, drifts_sum
      real :: heat_sources_sum, total_energy
      real, dimension(ny, nx) :: h_xy
      logical :: restart_time_step

      ! Gate the drive and drifts (drifts) terms exactly like the GKE advance:
      ! the drive is on when include_drive is set, the drifts on when the knob is
      ! nonzero; the perpendicular hyperdiffusion on when D_hyper > 0.
      ! Note: advance_wstar_explicit / advance_wdriftx|y_explicit still return 0 in
      ! the off case (coefficient zero), but guarding here keeps the diagnostics
      ! exactly zero (and cheaper) when the term is disabled.

      restart_time_step = .false.
      ia = 1

      !--------------------------------------------------------------------
      ! Zero all outputs
      !--------------------------------------------------------------------
      energy_diag%particle_energy_kxkyz = 0.
      energy_diag%dedt_proj_kxkyz = 0.
      energy_diag%drive_kxkyz = 0.
      energy_diag%coll_kxkyz = 0.
      energy_diag%diss_vpa_kxkyz = 0.
      energy_diag%diss_zed_kxkyz = 0.
      energy_diag%diss_perp_kxkyz = 0.
      energy_diag%nonlinearity_kxkyz = 0.
      energy_diag%streaming_kxkyz = 0.
      energy_diag%mirror_kxkyz = 0.
      energy_diag%drifts_kxkyz = 0.
      if (write_energy_xy) then
         energy_diag%field_energy_xy_z = 0.
         energy_diag%particle_energy_xy_z = 0.
      end if
      if (write_energy_vmu) then
         energy_diag%particle_energy_vmu = 0.
         energy_diag%dedt_proj_vmu = 0.
         energy_diag%drive_vmu = 0.
         energy_diag%coll_vmu = 0.
         energy_diag%diss_vpa_vmu = 0.
         energy_diag%diss_zed_vmu = 0.
         energy_diag%diss_perp_vmu = 0.
         energy_diag%nonlinearity_vmu = 0.
         energy_diag%streaming_vmu = 0.
         energy_diag%mirror_vmu = 0.
         energy_diag%drifts_vmu = 0.
      end if

      particle_energy = 0.
      dedt_proj_sum = 0.
      drive_sum = 0.
      coll_sum = 0.
      diss_vpa_sum = 0.
      diss_zed_sum = 0.
      diss_perp_sum = 0.
      nonlinearity_sum = 0.
      streaming_sum = 0.
      mirror_sum = 0.
      drifts_sum = 0.
      heat_sources_sum = 0.
      total_energy = 0.

      ! FLAG - electrostatic for now
      if (fphi > epsilon(0.0)) then

         !--------------------------------------------------------------------
         ! Particle energy H_p(t_n) = 0.5 * <h* | h>
         !--------------------------------------------------------------------
         energy_diag%factor_spec = spec%dens * spec%temp
         call get_one_energy_term_kxkyz(h, h, energy_diag%factor_spec, energy_diag%energy_spec, particle_energy, energy_diag%particle_energy_kxkyz)
         particle_energy = 0.5 * particle_energy
         energy_diag%energy_spec = 0.5 * energy_diag%energy_spec
         energy_diag%particle_energy_kxkyz = 0.5 * energy_diag%particle_energy_kxkyz

         if (write_energy_vmu) then
            call get_one_energy_term_vmu(h, h, energy_diag%factor_spec, energy_diag%particle_energy_vmu)
            energy_diag%particle_energy_vmu = 0.5 * energy_diag%particle_energy_vmu
         end if

         ! Real-space particle-energy density (unfolded: no dVolume/volume_total
         ! and no 1/(nx*ny)). Each rank handles its local vmu slices, then the
         ! density is summed over all ranks.
         if (write_energy_xy) then
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               is = is_idx(vmu_lo, ivmu)
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
                     call transform_kxky_to_xy(h(:, :, iz, it, ivmu), h_xy)
                     ! 0.25 = 0.5 from the 0.5*<h|h> energy norm times the 0.5
                     ! mode-counting factor used by get_one_energy_term_kxkyz.
                     energy_diag%particle_energy_xy_z(:, :, iz, it, is) = &
                        energy_diag%particle_energy_xy_z(:, :, iz, it, is) + &
                        transpose(0.25 * energy_diag%factor_spec(is) &
                           * wgts_vpa(iv) * wgts_mu(ia, iz, imu) &
                           / (maxwell_fac(is) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is)) &
                           * h_xy**2)
                  end do
               end do
            end do
            call sum_allreduce(energy_diag%particle_energy_xy_z)
         end if

         ! Field energy H_chi (electrostatic) and total energy H = H_p + H_chi
         if (write_energy_xy) then
            call get_field_energy(phi, energy_diag%field_energy_sum, energy_diag%field_energy_kxkyz, &
                                  write_energy_vmu, energy_diag%field_energy_vmu, &
                                  write_energy_xy, energy_diag%field_energy_xy_z)
         else
            call get_field_energy(phi, energy_diag%field_energy_sum, energy_diag%field_energy_kxkyz, &
                                  write_energy_vmu, energy_diag%field_energy_vmu, .false.)
         end if
         total_energy = particle_energy + energy_diag%field_energy_sum
         energy_diag%particle_energy_sum = particle_energy
         energy_diag%total_energy_sum = total_energy

         !--------------------------------------------------------------------
         ! col 4: drive (project the wstar rate onto h)
         !--------------------------------------------------------------------
         energy_diag%factor_spec = spec%dens * spec%temp
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 1), 'energy drive')
         if (include_drive) then
            g1 = 0.
            call advance_wstar_explicit(phi, g1)
            g1 = g1 / code_dt
            call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%drive_spec, drive_sum, energy_diag%drive_kxkyz)
            if (write_energy_vmu) call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%drive_vmu)
         else
            drive_sum = 0.
            energy_diag%drive_spec = 0.
         end if
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 1), 'energy drive')

         !--------------------------------------------------------------------
         ! col 6: collisions
         !--------------------------------------------------------------------
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 2), 'energy collisions')
         if (include_collisions) then
            energy_diag%factor_spec = spec%dens * spec%temp
            g1 = 0.
            call advance_collisions_explicit(g, phi, mod_bpar, g1)
            g1 = g1 / code_dt
            call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%coll_spec, coll_sum, energy_diag%coll_kxkyz)
            if (write_energy_vmu) call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%coll_vmu)
         else
            coll_sum = 0.
            energy_diag%coll_spec = 0.
         end if
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 2), 'energy collisions')

         !--------------------------------------------------------------------
         ! col 7: hyp_v (parallel-velocity hyperdiffusion)
         !--------------------------------------------------------------------
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 3), 'energy diss_vpa')
         if (hyp_vpa) then
            energy_diag%factor_spec = spec%dens * spec%temp
            g1 = 0.
            call advance_hyper_vpa(g, g1)
            g1 = g1 / code_dt
            call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%diss_vpa_spec, diss_vpa_sum, energy_diag%diss_vpa_kxkyz)
            if (write_energy_vmu) call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%diss_vpa_vmu)
         else
            diss_vpa_sum = 0.
            energy_diag%diss_vpa_spec = 0.
         end if
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 3), 'energy diss_vpa')

         !--------------------------------------------------------------------
         ! col 8: hyp_z (z hyperdiffusion)
         !--------------------------------------------------------------------
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 4), 'energy diss_zed')
         if (hyp_zed) then
            energy_diag%factor_spec = spec%dens * spec%temp
            g1 = 0.
            call advance_hyper_zed(g, g1)
            g1 = g1 / code_dt
            call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%diss_zed_spec, diss_zed_sum, energy_diag%diss_zed_kxkyz)
            if (write_energy_vmu) call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%diss_zed_vmu)
         else
            diss_zed_sum = 0.
            energy_diag%diss_zed_spec = 0.
         end if
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 4), 'energy diss_zed')

         !--------------------------------------------------------------------
         ! col 9: hyp_x/y (perpendicular hyperdiffusion)
         !--------------------------------------------------------------------
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 5), 'energy diss_perp')
         if (D_hyper > 0.0) then
            energy_diag%factor_spec = -D_hyper * spec%dens * spec%temp
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               is = is_idx(vmu_lo, ivmu)
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
                     g1(:, :, iz, it, ivmu) = g(:, :, iz, it, ivmu) * (kperp2(:, :, ia, iz) / k2max)**2
                  end do
               end do
            end do
            call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%diss_perp_spec, diss_perp_sum, energy_diag%diss_perp_kxkyz)
            if (write_energy_vmu) call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%diss_perp_vmu)
         else
            diss_perp_sum = 0.
            energy_diag%diss_perp_spec = 0.
         end if
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 5), 'energy diss_perp')

         !--------------------------------------------------------------------
         ! col 10: nonlinearity
         !--------------------------------------------------------------------
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 6), 'energy nonlinearity')
         if (include_nonlinear) then
            energy_diag%factor_spec = spec%dens * spec%temp
            g1 = 0.
            call advance_ExB_nonlinearity(g, g1, restart_time_step, istep)
            g1 = g1 / code_dt
            call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%nonlinearity_spec, nonlinearity_sum, energy_diag%nonlinearity_kxkyz)
            if (write_energy_vmu) call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%nonlinearity_vmu)
         else
            nonlinearity_sum = 0.
            energy_diag%nonlinearity_spec = 0.
         end if
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 6), 'energy nonlinearity')

         !--------------------------------------------------------------------
         ! col 11: parallel streaming term
         !--------------------------------------------------------------------
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 7), 'energy streaming')
         if (include_parallel_streaming) then
            energy_diag%factor_spec = spec%dens * spec%temp
            g1 = 0.
            call advance_parallel_streaming_explicit(g, phi, mod_bpar, g1)
            g1 = g1 / code_dt
            call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%streaming_spec, streaming_sum, energy_diag%streaming_kxkyz)
            if (write_energy_vmu) call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%streaming_vmu)
         else
            streaming_sum = 0.
            energy_diag%streaming_spec = 0.
         end if
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 7), 'energy streaming')

         !--------------------------------------------------------------------
         ! col 12: mirror term
         !--------------------------------------------------------------------
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 8), 'energy mirror')
         if (include_mirror) then
            energy_diag%factor_spec = spec%dens * spec%temp
            g1 = 0.
            call advance_mirror_explicit(g, g1)
            g1 = g1 / code_dt
            call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%mirror_spec, mirror_sum, energy_diag%mirror_kxkyz)
            if (write_energy_vmu) call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%mirror_vmu)
         else
            mirror_sum = 0.
            energy_diag%mirror_spec = 0.
         end if
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 8), 'energy mirror')

         !--------------------------------------------------------------------
         ! col 12: drifts + remaining = the drift/drifts term, defined as
         ! the remainder of the projection identity. Here we
         ! compute the drift projection directly and assign it to drifts, which
         ! makes the projection identity close exactly by construction. Guarded on
         ! the drift knobs (matching the GKE advance), not just the master flags.
         !--------------------------------------------------------------------
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 9), 'energy drifts')
         if (abs(xdriftknob) > epsilon(0.0) .or. abs(ydriftknob) > epsilon(0.0)) then
            energy_diag%factor_spec = spec%dens * spec%temp
            g1 = 0.
            call advance_wdriftx_explicit(g, phi, mod_bpar, g1)
            call advance_wdrifty_explicit(g, phi, mod_bpar, g1)
            g1 = g1 / code_dt
            call get_one_energy_term_kxkyz(h, g1, energy_diag%factor_spec, energy_diag%drifts_spec, drifts_sum, energy_diag%drifts_kxkyz)
            if (write_energy_vmu) call get_one_energy_term_vmu(h, g1, energy_diag%factor_spec, energy_diag%drifts_vmu)
         else
            drifts_sum = 0.
            energy_diag%drifts_spec = 0.
         end if
         if (energy_diagnostics_debug) call time_message(.false., time_energy_terms(:, 9), 'energy drifts')

         !--------------------------------------------------------------------
         ! col 3: dedt_projected = <c_f | RHS(t_n)> = sum of all per-term projections
         ! (col 5 heat_sources = 0 in the local flux-tube limit)
         !--------------------------------------------------------------------
         ! dedt_projected = <c_f | RHS(t_n)> = sum of all per-term projections
         !--------------------------------------------------------------------
         dedt_proj_sum = drive_sum + heat_sources_sum + coll_sum + diss_vpa_sum + diss_zed_sum + diss_perp_sum + &
                          nonlinearity_sum + streaming_sum + mirror_sum + drifts_sum
         energy_diag%dedt_proj_sum = dedt_proj_sum

         ! Per-species dedt_projected = sum of per-species term projections
         energy_diag%dedt_proj_spec = energy_diag%drive_spec + energy_diag%coll_spec + energy_diag%diss_vpa_spec + &
                                      energy_diag%diss_zed_spec + energy_diag%diss_perp_spec + energy_diag%nonlinearity_spec + &
                                      energy_diag%streaming_spec + energy_diag%mirror_spec + energy_diag%drifts_spec

         ! Per-mode dedt_projected = sum of per-mode term projections
         energy_diag%dedt_proj_kxkyz = energy_diag%drive_kxkyz + energy_diag%coll_kxkyz + energy_diag%diss_vpa_kxkyz + &
                                       energy_diag%diss_zed_kxkyz + energy_diag%diss_perp_kxkyz + energy_diag%nonlinearity_kxkyz + &
                                       energy_diag%streaming_kxkyz + energy_diag%mirror_kxkyz + energy_diag%drifts_kxkyz
         if (write_energy_vmu) then
            energy_diag%dedt_proj_vmu = energy_diag%drive_vmu + energy_diag%coll_vmu + energy_diag%diss_vpa_vmu + &
                                        energy_diag%diss_zed_vmu + energy_diag%diss_perp_vmu + energy_diag%nonlinearity_vmu + &
                                        energy_diag%streaming_vmu + energy_diag%mirror_vmu + energy_diag%drifts_vmu
         end if

      end if

      !--------------------------------------------------------------------
      ! exact change in energy, kept SEPARATE from dedt_projected: the two
      ! differ for fast-evolving fields, so dedt_total is reported but never
      ! equated to the projection. We track the previous diagnostic call's
      ! particle/field energy and time over the whole run; on the first call
      ! last_*_energy = 0 and last_time = -1.
      !--------------------------------------------------------------------
      energy_diag%dedt_particle_sum = (particle_energy - energy_diag%last_particle_energy) / (code_time - energy_diag%last_time)
      energy_diag%dedt_field_sum = (energy_diag%field_energy_sum - energy_diag%last_field_energy) / (code_time - energy_diag%last_time)
      energy_diag%dedt_total_sum = energy_diag%dedt_particle_sum + energy_diag%dedt_field_sum
      energy_diag%last_particle_energy = particle_energy
      energy_diag%last_field_energy = energy_diag%field_energy_sum
      energy_diag%last_total_energy = total_energy
      energy_diag%last_time = code_time

      !--------------------------------------------------------------------
      ! ASCII output: 18 columns
      !   time, particle_energy, field_energy, total_energy, dedt_particle,
      !   dedt_field, dedt_total, dedt_projected, drive, heat_sources,
      !   collisions, diss_vpa, diss_zed, diss_perp, drifts, streaming,
      !   mirror, nonlinearity
      !--------------------------------------------------------------------
      if (proc0) then
         write (energy_unit, '(18e20.8E3)') code_time, particle_energy, energy_diag%field_energy_sum, total_energy, &
            energy_diag%dedt_particle_sum, energy_diag%dedt_field_sum, energy_diag%dedt_total_sum, dedt_proj_sum, &
            drive_sum, heat_sources_sum, coll_sum, diss_vpa_sum, diss_zed_sum, &
            diss_perp_sum, drifts_sum, streaming_sum, mirror_sum, nonlinearity_sum
         call flush_output_file(energy_unit)
      end if

      ! Per-species ascii files (optional)
      if (proc0 .and. allocated(energy_diag%energy_species_units)) then
         do is = 1, nspec
            write (energy_diag%energy_species_units(is), '(12e15.6E3)') &
               code_time, energy_diag%energy_spec(is), energy_diag%dedt_proj_spec(is), &
               energy_diag%drive_spec(is), energy_diag%coll_spec(is), &
               energy_diag%diss_vpa_spec(is), energy_diag%diss_zed_spec(is), &
               energy_diag%diss_perp_spec(is), energy_diag%nonlinearity_spec(is), &
               energy_diag%streaming_spec(is), energy_diag%mirror_spec(is), &
               energy_diag%drifts_spec(is)
         end do
         do is = 1, nspec
            call flush_output_file(energy_diag%energy_species_units(is))
         end do
      end if

   end subroutine get_free_energy

   !> Finalise the free-energy diagnostic: deallocate arrays and close files.
   subroutine finish_energy_diagnostic
      implicit none

      if (allocated(energy_diag%energy_spec)) deallocate (energy_diag%energy_spec)
      if (allocated(energy_diag%dedt_proj_spec)) deallocate (energy_diag%dedt_proj_spec)
      if (allocated(energy_diag%drive_spec)) deallocate (energy_diag%drive_spec)
      if (allocated(energy_diag%coll_spec)) deallocate (energy_diag%coll_spec)
      if (allocated(energy_diag%diss_vpa_spec)) deallocate (energy_diag%diss_vpa_spec)
      if (allocated(energy_diag%diss_zed_spec)) deallocate (energy_diag%diss_zed_spec)
      if (allocated(energy_diag%diss_perp_spec)) deallocate (energy_diag%diss_perp_spec)
      if (allocated(energy_diag%nonlinearity_spec)) deallocate (energy_diag%nonlinearity_spec)
      if (allocated(energy_diag%streaming_spec)) deallocate (energy_diag%streaming_spec)
      if (allocated(energy_diag%mirror_spec)) deallocate (energy_diag%mirror_spec)
      if (allocated(energy_diag%drifts_spec)) deallocate (energy_diag%drifts_spec)
      if (allocated(energy_diag%factor_spec)) deallocate (energy_diag%factor_spec)

      if (allocated(energy_diag%particle_energy_kxkyz)) deallocate (energy_diag%particle_energy_kxkyz)
      if (allocated(energy_diag%dedt_proj_kxkyz)) deallocate (energy_diag%dedt_proj_kxkyz)
      if (allocated(energy_diag%drive_kxkyz)) deallocate (energy_diag%drive_kxkyz)
      if (allocated(energy_diag%coll_kxkyz)) deallocate (energy_diag%coll_kxkyz)
      if (allocated(energy_diag%diss_vpa_kxkyz)) deallocate (energy_diag%diss_vpa_kxkyz)
      if (allocated(energy_diag%diss_zed_kxkyz)) deallocate (energy_diag%diss_zed_kxkyz)
      if (allocated(energy_diag%diss_perp_kxkyz)) deallocate (energy_diag%diss_perp_kxkyz)
      if (allocated(energy_diag%nonlinearity_kxkyz)) deallocate (energy_diag%nonlinearity_kxkyz)
      if (allocated(energy_diag%streaming_kxkyz)) deallocate (energy_diag%streaming_kxkyz)
      if (allocated(energy_diag%mirror_kxkyz)) deallocate (energy_diag%mirror_kxkyz)
      if (allocated(energy_diag%drifts_kxkyz)) deallocate (energy_diag%drifts_kxkyz)
      if (allocated(energy_diag%field_energy_kxkyz)) deallocate (energy_diag%field_energy_kxkyz)
      if (allocated(energy_diag%field_energy_xy_z)) deallocate (energy_diag%field_energy_xy_z)
      if (allocated(energy_diag%particle_energy_xy_z)) deallocate (energy_diag%particle_energy_xy_z)

      if (allocated(energy_diag%particle_energy_vmu)) deallocate (energy_diag%particle_energy_vmu)
      if (allocated(energy_diag%dedt_proj_vmu)) deallocate (energy_diag%dedt_proj_vmu)
      if (allocated(energy_diag%drive_vmu)) deallocate (energy_diag%drive_vmu)
      if (allocated(energy_diag%coll_vmu)) deallocate (energy_diag%coll_vmu)
      if (allocated(energy_diag%diss_vpa_vmu)) deallocate (energy_diag%diss_vpa_vmu)
      if (allocated(energy_diag%diss_zed_vmu)) deallocate (energy_diag%diss_zed_vmu)
      if (allocated(energy_diag%diss_perp_vmu)) deallocate (energy_diag%diss_perp_vmu)
      if (allocated(energy_diag%nonlinearity_vmu)) deallocate (energy_diag%nonlinearity_vmu)
      if (allocated(energy_diag%streaming_vmu)) deallocate (energy_diag%streaming_vmu)
      if (allocated(energy_diag%mirror_vmu)) deallocate (energy_diag%mirror_vmu)
      if (allocated(energy_diag%drifts_vmu)) deallocate (energy_diag%drifts_vmu)
      if (allocated(energy_diag%field_energy_vmu)) deallocate (energy_diag%field_energy_vmu)
      if (allocated(energy_diag%field_energy_vmu_dummy)) deallocate (energy_diag%field_energy_vmu_dummy)

      if (allocated(energy_diag%velocity_integral1)) deallocate (energy_diag%velocity_integral1)
      if (allocated(energy_diag%spatial_integral1)) deallocate (energy_diag%spatial_integral1)

      if (allocated(energy_diag%energy_species_units)) deallocate (energy_diag%energy_species_units)
      if (allocated(energy_diag%energy_species_labels)) deallocate (energy_diag%energy_species_labels)

   end subroutine finish_energy_diagnostic

end module diagnostics_energy
