module dist_fn

   implicit none

   public :: init_gxyz
   public :: init_dist_fn, finish_dist_fn
   public :: checksum

   private

   interface checksum
      module procedure checksum_field
      module procedure checksum_dist
   end interface

   logical :: dist_fn_initialized = .false.
   logical :: gxyz_initialized = .false.
   logical :: kp2init = .false.
   logical :: dkp2drinit = .false.
   logical :: vp2init = .false.

   logical :: debug = .false.

contains

   subroutine init_gxyz(restarted)

      use dist_fn_arrays, only: gvmu, gold, gnew
      use redistribute, only: gather, scatter
      use dist_redistribute, only: kxkyz2vmu
      use physics_flags, only: radial_variation
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
      use kt_grids, only: nalpha, nakx, naky, multiply_by_rho
      use vpamu_grids, only: mu, vpa, vperp2
      use zgrid, only: nzgrid, ntubes
      use species, only: spec, pfac
      use stella_geometry, only: dBdrho, gfac

      implicit none

      real :: corr
      integer :: ivmu, is, imu, iv, it, iz, ia
      real, dimension(:, :), allocatable :: energy
      complex, dimension(:, :), allocatable :: g0k
      logical, intent(in) :: restarted

      if (gxyz_initialized) return
      gxyz_initialized = .false.

      ! get version of g that has ky,kx,z local
      call gather(kxkyz2vmu, gvmu, gnew)

      ia = 1

      !calculate radial corrections to F0 for use in Krook operator, as well as g1 from initialization
      if (radial_variation) then
         !init_g uses maxwellians, so account for variation in temperature, density, and B

         allocate (energy(nalpha, -nzgrid:nzgrid))
         allocate (g0k(naky, nakx))

         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            energy = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid

                  corr = -(pfac * (spec(is)%fprim + spec(is)%tprim * (energy(ia, iz) - 1.5)) &
                           + 2 * gfac * mu(imu) * dBdrho(iz))

                  if (.not. restarted) then
                     g0k = corr * gnew(:, :, iz, it, ivmu)
                     call multiply_by_rho(g0k)
                     gnew(:, :, iz, it, ivmu) = gnew(:, :, iz, it, ivmu) + g0k
                  end if
               end do
            end do
         end do
         deallocate (energy, g0k)

         if (.not. restarted) call scatter(kxkyz2vmu, gnew, gvmu)
      end if

      gold = gnew

   end subroutine init_gxyz

   subroutine init_dist_fn

      use mp, only: proc0
      use physics_flags, only: radial_variation
      use stella_layouts, only: init_dist_fn_layouts
      use gyro_averages, only: init_bessel

      implicit none

      if (dist_fn_initialized) return
      dist_fn_initialized = .true.

      debug = debug .and. proc0

      if (debug) write (*, *) 'dist_fn::init_dist_fn::allocate_arrays'
      call allocate_arrays

      !> allocate and initialise kperp2 and dkperp2dr
      if (debug) write (*, *) 'dist_fn::init_dist_fn::init_kperp2'
      call init_kperp2
      if (radial_variation) call init_dkperp2dr

      !> allocate and initialise vperp2
      if (debug) write (*, *) 'dist_fn::init_dist_fn::init_vperp2'
      call init_vperp2

      !> init_bessel sets up arrays needed for gyro-averaging;
      !> for a flux tube simulation, this is j0 and j1;
      !> for a flux annulus simulation, gyro-averaging is non-local in ky
      !> and so more effort is required
      if (debug) write (*, *) 'dist_fn::init_dist_fn::init_bessel'
      call init_bessel

   end subroutine init_dist_fn

   !> init_kperp2 allocates and initialises the kperp2 array
   subroutine init_kperp2

      use dist_fn_arrays, only: kperp2
      use stella_geometry, only: gds2, gds21, gds22
      use stella_geometry, only: geo_surf, q_as_x
      use zgrid, only: nzgrid
      use kt_grids, only: naky, nakx, theta0
      use kt_grids, only: akx, aky
      use kt_grids, only: zonal_mode
      use kt_grids, only: nalpha

      implicit none

      integer :: iky, ikx

      if (kp2init) return
      kp2init = .true.

      !> allocate the kperp2 array to contain |k_perp|^2
      allocate (kperp2(naky, nakx, nalpha, -nzgrid:nzgrid))

      do iky = 1, naky
         if (zonal_mode(iky)) then
            do ikx = 1, nakx
               if (q_as_x) then
                  kperp2(iky, ikx, :, :) = akx(ikx) * akx(ikx) * gds22
               else
                  kperp2(iky, ikx, :, :) = akx(ikx) * akx(ikx) * gds22 / (geo_surf%shat**2)
               end if
            end do
         else
            do ikx = 1, nakx
               kperp2(iky, ikx, :, :) = aky(iky) * aky(iky) &
                                        * (gds2 + 2.0 * theta0(iky, ikx) * gds21 &
                                           + theta0(iky, ikx) * theta0(iky, ikx) * gds22)
            end do
         end if
      end do

      ! NB: should really avoid this by using higher resolution when reading in VMEC geometry and then
      ! NB: course-graining if necessary to map onto lower-resolution stella grid
      ! ensure kperp2 is positive everywhere (only might go negative if using full-flux-surface due to interpolation)
      where (kperp2 < 0.0)
         kperp2 = 0.0
      end where

      call enforce_single_valued_kperp2

   end subroutine init_kperp2

   !> init_dkperp2dr allocates and initialises the dkperp2dr array, needed for radial variation
   subroutine init_dkperp2dr

      use dist_fn_arrays, only: kperp2, dkperp2dr
      use stella_geometry, only: dgds2dr, dgds21dr, dgds22dr
      use stella_geometry, only: geo_surf, q_as_x
      use zgrid, only: nzgrid
      use kt_grids, only: naky, nakx, theta0
      use kt_grids, only: akx, aky
      use kt_grids, only: zonal_mode
      use kt_grids, only: nalpha

      implicit none

      integer :: iky, ikx

      if (dkp2drinit) return
      dkp2drinit = .true.

      allocate (dkperp2dr(naky, nakx, nalpha, -nzgrid:nzgrid))
      do iky = 1, naky
         if (zonal_mode(iky)) then
            do ikx = 1, nakx
               if (q_as_x) then
                  where (kperp2(iky, ikx, :, :) > epsilon(0.0))
                     dkperp2dr(iky, ikx, :, :) = akx(ikx) * akx(ikx) * dgds22dr / kperp2(iky, ikx, :, :)
                  elsewhere
                     dkperp2dr(iky, ikx, :, :) = 0.0
                  end where
               else
                  where (kperp2(iky, ikx, :, :) > epsilon(0.0))
                     dkperp2dr(iky, ikx, :, :) = akx(ikx) * akx(ikx) * dgds22dr / (geo_surf%shat**2 * kperp2(iky, ikx, :, :))
                  elsewhere
                     dkperp2dr(iky, ikx, :, :) = 0.0
                  end where
               end if
            end do
         else
            do ikx = 1, nakx
               dkperp2dr(iky, ikx, :, :) = aky(iky) * aky(iky) &
                                           * (dgds2dr + 2.0 * theta0(iky, ikx) * dgds21dr &
                                              + theta0(iky, ikx) * theta0(iky, ikx) * dgds22dr)
               dkperp2dr(iky, ikx, :, :) = dkperp2dr(iky, ikx, :, :) / kperp2(iky, ikx, :, :)
               if (any(kperp2(iky, ikx, :, :) < epsilon(0.))) dkperp2dr(iky, ikx, :, :) = 0.
            end do
         end if
      end do

   end subroutine init_dkperp2dr

   subroutine enforce_single_valued_kperp2

      use dist_fn_arrays, only: kperp2
      use kt_grids, only: naky, nalpha
      use zgrid, only: nzgrid
      use extended_zgrid, only: neigen, nsegments, ikxmod

      implicit none

      integer :: iky, ie, iseg
      real, dimension(:), allocatable :: tmp

      allocate (tmp(nalpha)); tmp = 0.0

      do iky = 1, naky
         do ie = 1, neigen(iky)
            if (nsegments(ie, iky) > 1) then
               do iseg = 2, nsegments(ie, iky)
                  tmp = 0.5 * (kperp2(iky, ikxmod(iseg - 1, ie, iky), :, nzgrid) + kperp2(iky, ikxmod(iseg, ie, iky), :, -nzgrid))
                  kperp2(iky, ikxmod(iseg, ie, iky), :, -nzgrid) = tmp
                  kperp2(iky, ikxmod(iseg - 1, ie, iky), :, nzgrid) = tmp
               end do
            end if
         end do
      end do

      deallocate (tmp)

   end subroutine enforce_single_valued_kperp2

   subroutine allocate_arrays

      use stella_layouts, only: kxkyz_lo, vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use vpamu_grids, only: nvpa, nmu
      use dist_fn_arrays, only: gnew, gold, g_gyro, g_secondary_source
      use dist_fn_arrays, only: gvmu

      implicit none

      if (.not. allocated(gnew)) &
         allocate (gnew(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      gnew = 0.
      if (.not. allocated(gold)) &
         allocate (gold(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      gold = 0.
      if (.not. allocated(g_secondary_source)) &
         allocate (g_secondary_source(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g_secondary_source = 0.
      if (.not. allocated(g_gyro)) &
         allocate (g_gyro(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g_gyro = 0.
      if (.not. allocated(gvmu)) &
         allocate (gvmu(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      gvmu = 0.

   end subroutine allocate_arrays

   subroutine init_vperp2

      use stella_geometry, only: bmag
      use zgrid, only: nzgrid
      use vpamu_grids, only: vperp2
      use vpamu_grids, only: nmu, mu
      use kt_grids, only: nalpha

      implicit none

      integer :: imu

      if (vp2init) return
      vp2init = .true.

      if (.not. allocated(vperp2)) allocate (vperp2(nalpha, -nzgrid:nzgrid, nmu)); vperp2 = 0.

      do imu = 1, nmu
         vperp2(:, :, imu) = 2.0 * mu(imu) * bmag
      end do

   end subroutine init_vperp2

   subroutine finish_dist_fn

      use gyro_averages, only: finish_bessel

      implicit none

      call finish_bessel
      call finish_kperp2
      call finish_vperp2
      call deallocate_arrays

      dist_fn_initialized = .false.
      gxyz_initialized = .false.

   end subroutine finish_dist_fn

   subroutine deallocate_arrays

      use dist_fn_arrays, only: gnew, gold, g_secondary_source, g_gyro, gvmu

      implicit none

      if (allocated(gnew)) deallocate (gnew)
      if (allocated(gold)) deallocate (gold)
      if (allocated(g_secondary_source)) deallocate (g_secondary_source)
      if (allocated(g_gyro)) deallocate (g_gyro)
      if (allocated(gvmu)) deallocate (gvmu)

   end subroutine deallocate_arrays

   subroutine finish_kperp2

      use dist_fn_arrays, only: kperp2, dkperp2dr

      implicit none

      if (allocated(kperp2)) deallocate (kperp2)
      if (allocated(dkperp2dr)) deallocate (dkperp2dr)

      kp2init = .false.
      dkp2drinit = .false.

   end subroutine finish_kperp2

   subroutine finish_vperp2

      use vpamu_grids, only: vperp2

      implicit none

      if (allocated(vperp2)) deallocate (vperp2)

      vp2init = .false.

   end subroutine finish_vperp2

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

end module dist_fn
