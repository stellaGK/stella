module hyper

   implicit none

   public :: read_parameters_hyper
   public :: init_hyper
   public :: advance_hyper_dissipation

   private

   logical :: use_physical_ksqr
   real :: D_hyper
   real :: tfac
   real :: k2max

contains

   subroutine read_parameters_hyper

      use file_utils, only: input_unit_exist
      use physics_flags, only: full_flux_surface, radial_variation
      use mp, only: proc0, broadcast

      implicit none

      namelist /hyper/ D_hyper, use_physical_ksqr

      integer :: in_file
      logical :: dexist

      if (proc0) then
         use_physical_ksqr = .not. (full_flux_surface .or. radial_variation)
         D_hyper = 0.05

         in_file = input_unit_exist("hyper", dexist)
         if (dexist) read (unit=in_file, nml=hyper)
      end if

      call broadcast(use_physical_ksqr)
      call broadcast(D_hyper)

   end subroutine read_parameters_hyper

   subroutine init_hyper

      use kt_grids, only: ikx_max, naky
      use kt_grids, only: aky, akx
      use stella_geometry, only: geo_surf, q_as_x
      use dist_fn_arrays, only: kperp2

      implicit none

      integer :: ia

      ia = 1

      if (.not. use_physical_ksqr) then
         ! avoid spatially dependent kperp

         !get k2max at outboard midplane
         k2max = akx(ikx_max)**2 + aky(naky)**2
         if (k2max < epsilon(0.0)) k2max = 1.0

         tfac = geo_surf%shat**2
         if (q_as_x) tfac = 1.0
      else
         !get k2max at outboard midplane
         k2max = maxval(kperp2(:,:,ia,0)) 
         if (k2max < epsilon(0.0)) k2max = 1.0
      end if

   end subroutine init_hyper

   subroutine advance_hyper_dissipation(g)

      use stella_time, only: code_dt
      use zgrid, only: nzgrid, ntubes, zed
      use stella_layouts, only: vmu_lo
      use dist_fn_arrays, only: kperp2
      use kt_grids, only: naky
      use kt_grids, only: aky, akx, theta0, zonal_mode

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g

      integer :: ia, ivmu, iz, it, iky

      ia = 1

      if (.not. use_physical_ksqr) then
         ! avoid spatially dependent kperp
         ! add in hyper-dissipation of form dg/dt = -D*(k/kmax)^4*g
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
         ! add in hyper-dissipation of form dg/dt = -D*(k/kmax)^4*g
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            g(:, :, :, :, ivmu) = g(:, :, :, :, ivmu) / (1.+code_dt * (spread(kperp2(:, :, ia, :), 4, ntubes) / k2max)**2 * D_hyper)
         end do
      end if

   end subroutine advance_hyper_dissipation

end module hyper
