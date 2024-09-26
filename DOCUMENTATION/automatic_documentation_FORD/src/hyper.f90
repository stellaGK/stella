module hyper

   implicit none

   public :: read_parameters_hyper
   public :: init_hyper
   public :: advance_hyper_dissipation
   public :: advance_hyper_vpa
   public :: advance_hyper_zed
   public :: D_hyper, D_zed, D_vpa
   public :: hyp_vpa, hyp_zed
   public :: k2max

   private

   logical :: use_physical_ksqr, scale_to_outboard
   real :: D_hyper, D_zed, D_vpa
   logical :: hyp_vpa, hyp_zed
   real :: tfac
   real :: k2max

contains

   subroutine read_parameters_hyper

      use file_utils, only: input_unit_exist
      use parameters_physics, only: full_flux_surface, radial_variation
      use mp, only: proc0, broadcast

      implicit none

      namelist /hyper/ D_hyper, D_zed, D_vpa, hyp_zed, hyp_vpa, use_physical_ksqr, scale_to_outboard

      integer :: in_file
      logical :: dexist

      if (proc0) then
         use_physical_ksqr = .not. (full_flux_surface .or. radial_variation)  ! use kperp2, instead of akx^2 + aky^2
         scale_to_outboard = .false.                                          ! scales hyperdissipation to zed = 0
         D_hyper = 0.05
         D_zed = 0.05
         D_vpa = 0.05
         hyp_vpa = .false.
         hyp_zed = .false.

         in_file = input_unit_exist("hyper", dexist)
         if (dexist) read (unit=in_file, nml=hyper)
      end if

      call broadcast(use_physical_ksqr)
      call broadcast(scale_to_outboard)
      call broadcast(D_hyper)
      call broadcast(D_zed)
      call broadcast(D_vpa)
      call broadcast(hyp_vpa)
      call broadcast(hyp_zed)

   end subroutine read_parameters_hyper

   subroutine init_hyper

      use parameters_kxky_grids, only: ikx_max, nakx, naky
      use grids_kxky, only: aky, akx, theta0
      use zgrid, only: nzgrid, zed
      use geometry, only: geo_surf, q_as_x
      use arrays_dist_fn, only: kperp2

      implicit none

      integer :: iky, ikx, iz, ia
      real :: temp

      ia = 1

      if (.not. use_physical_ksqr) then
         !> avoid spatially dependent kperp (through the geometric coefficients)
         !> still allowed to vary along zed with global magnetic shear
         !> useful for full_flux_surface and radial_variation runs
         tfac = geo_surf%shat**2

         !> q_as_x uses a different definition of theta0
         if (q_as_x) tfac = 1.0

         if (scale_to_outboard) then
            !get k2max at outboard midplane
            k2max = akx(ikx_max)**2 + aky(naky)**2
         else
            k2max = -1.0
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  do iky = 1, naky
                     temp = aky(iky)**2 * (1.0 + tfac * (zed(iz) - theta0(iky, ikx))**2)
                     if (temp > k2max) k2max = temp
                  end do
               end do
            end do
         end if
      else
         if (scale_to_outboard) then
            !> get k2max at outboard midplane
            k2max = maxval(kperp2(:, :, ia, 0))
         else
            k2max = maxval(kperp2)
         end if
      end if
      if (k2max < epsilon(0.0)) k2max = 1.0

   end subroutine init_hyper

   subroutine advance_hyper_dissipation(g)

      use stella_time, only: code_dt
      use zgrid, only: nzgrid, ntubes, zed
      use stella_layouts, only: vmu_lo
      use arrays_dist_fn, only: kperp2
      use parameters_kxky_grids, only: naky
      use grids_kxky, only: aky, akx, theta0, zonal_mode

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g

      integer :: ia, ivmu, iz, it, iky

      ia = 1

      if (.not. use_physical_ksqr) then
         !> avoid spatially dependent kperp
         !> add in hyper-dissipation of form dg/dt = -D*(k/kmax)^4*g
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do iky = 1, naky
                     if (zonal_mode(iky)) then
                        g(iky, :, iz, it, ivmu) = g(iky, :, iz, it, ivmu) / (1.+code_dt * (akx(:)**2 / k2max)**2 * D_hyper)
                     else
                        g(iky, :, iz, it, ivmu) = g(iky, :, iz, it, ivmu) / (1.+code_dt * (aky(iky)**2 &
                                                                                 * (1.0 + tfac * (zed(iz) - theta0(iky, :))**2) / k2max)**2 * D_hyper)
                     end if
                  end do
               end do
            end do
         end do
      else
         !> add in hyper-dissipation of form dg/dt = -D*(k/kmax)^4*g
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            g(:, :, :, :, ivmu) = g(:, :, :, :, ivmu) / (1.+code_dt * (spread(kperp2(:, :, ia, :), 4, ntubes) / k2max)**2 * D_hyper)
         end do
      end if

   end subroutine advance_hyper_dissipation
 
   subroutine advance_hyper_vpa(g, dgdvpa)
   !> computes the fourth derivative of g in vpa and returns this in dgdvpa
   !> multiplied by the vpa diffusion coefficient
      use stella_time, only: code_dt
      use zgrid, only: nzgrid
      use stella_layouts, only: vmu_lo, kxkyz_lo
      use redistribute, only: gather, scatter
      use dist_redistribute, only: kxkyz2vmu
      use vpamu_grids, only: nmu, nvpa, dvpa

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dgdvpa

      complex, dimension(:, :, :), allocatable :: g0v, g1v

      allocate (g0v(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      allocate (g1v(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

      call scatter(kxkyz2vmu, g, g0v)
      call get_dgdvpa_fourth_order(g0v, g1v)
      call gather(kxkyz2vmu, g1v, dgdvpa)
      dgdvpa = -code_dt * D_vpa * dvpa**4 / 16 * dgdvpa
      
      deallocate (g0v)
      deallocate (g1v)

   end subroutine advance_hyper_vpa

   subroutine get_dgdvpa_fourth_order(g, gout)

      use finite_differences, only: fourth_derivate_second_centered_vpa
      use stella_layouts, only: kxkyz_lo, iz_idx, is_idx
      use vpamu_grids, only: nvpa, nmu, dvpa

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(inout) :: gout

      integer :: ikxkyz, imu, iz, is
      complex, dimension(:), allocatable :: tmp

      allocate (tmp(nvpa))
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            call fourth_derivate_second_centered_vpa(1, g(:, imu, ikxkyz), dvpa, tmp)
            gout(:, imu, ikxkyz) = tmp
         end do
      end do

      deallocate (tmp)
   end subroutine get_dgdvpa_fourth_order

   subroutine advance_hyper_zed(g, dgdz)
   !> computes the fourth derivative of g in z and returns this in
   !> dgdz multiplied by the z hyper diffusion coefficient
      use stella_time, only: code_dt
      use zgrid, only: nzgrid, delzed
      use stella_layouts, only: vmu_lo

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dgdz

      call get_dgdz_fourth_order(g, dgdz)
      dgdz = -code_dt * D_zed * delzed(0)**4 / 16 * dgdz

   end subroutine advance_hyper_zed

   subroutine get_dgdz_fourth_order(g, dgdz)

      use finite_differences, only: fourth_derivative_second_centered_zed
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx
      use zgrid, only: nzgrid, delzed, ntubes
      use extended_zgrid, only: neigen, nsegments
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: fill_zed_ghost_zones
      use extended_zgrid, only: periodic
      use parameters_kxky_grids, only: naky

      use stella_layouts, only: iv_idx, imu_idx, is_idx

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(inout) :: dgdz

      integer :: iseg, ie, iky, iv, it, ivmu, imu
      complex, dimension(2) :: gleft, gright
      ! FLAG -- assuming delta zed is equally spaced below!
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         do iky = 1, naky
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  do iseg = 1, nsegments(ie, iky)
                     ! first fill in ghost zones at boundaries in g(z)
                     call fill_zed_ghost_zones(it, iseg, ie, iky, g(:, :, :, :, ivmu), gleft, gright)
                     ! now get dg/dz
                     call fourth_derivative_second_centered_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
                                                                g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it, ivmu), &
                                                                delzed(0), gleft, gright, periodic(iky), &
                                                                dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it, ivmu))
                  end do
               end do
            end do
         end do
      end do
   end subroutine get_dgdz_fourth_order

end module hyper
