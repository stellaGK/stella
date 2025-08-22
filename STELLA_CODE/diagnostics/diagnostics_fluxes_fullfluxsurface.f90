
!###############################################################################
!#################### CALCULATE FLUXES FOR FULL FLUX SURFACE ###################
!###############################################################################
! TODO-GA-CLEANUP Fluxes and moments for FFS
 
module diagnostics_fluxes_fullfluxsurface

   implicit none
 
   public :: calculate_fluxes_fullfluxsurface
   public :: calculate_moments_fullfluxsurface

   private     

contains

!###############################################################################
!############################## CALCULATE FLUXES ###############################
!###############################################################################

   !============================================================================
   !============================ CALCULATE FLUXES ==============================
   !============================================================================
   !> Calculate the total particle, momentum and heat fluxes (pflux_vs_s, vflux_vs_s, qflux_vs_s)
   !> and the contributions from a given (kx,ky,z) location (pflux_kxkyz, vflux_kxkyz, qflux_kxkyz)
   !> inputs are the particle density (dens), parallel flow (upar) and pressure (pres)
   subroutine calculate_fluxes_fullfluxsurface(dens, upar, pres, pflux_vs_s, &
         vflux_vs_s, qflux_vs_s, pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts)

      use constants, only: zi
      use z_grid, only: nzgrid, delzed
      use parameters_kxky_grid, only: naky, nakx, ny
      use grids_kxky, only: aky, dy
      use arrays_fields, only: phi
      use geometry, only: grad_x, jacob

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: dens, upar, pres
      real, dimension(:), intent(out) :: pflux_vs_s, vflux_vs_s, qflux_vs_s
      real, dimension(:, :, -nzgrid:, :, :), intent(out) :: pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts

      integer :: iky, it
      real, dimension(:, :), allocatable :: flxfac 

      complex, dimension(:, :, :), allocatable :: dphidy

      !> assume a single flux annulus
      it = 1
      
      !> Set all fluxes to zero just in case
      pflux_vs_s = 0.; vflux_vs_s = 0.; qflux_vs_s = 0.
      pflux_kxkyzts = 0.; vflux_kxkyzts = 0.; qflux_kxkyzts = 0.

      allocate (dphidy(naky, nakx, -nzgrid:nzgrid))

      !> Obtain the y-component of the electric field that appears as a factor
      !> in the flux expression due to the radial component of the ExB velocity
      do iky = 1, naky
	 	!dphidy(iky, :, :) = zi * aky(iky) * conjg(phi(iky, :, :, it))
         dphidy(iky, :, :) = phi(iky, :, :, it) * aky(iky)
      end do
      
      !> Calculate Jacobian for spacial integral 
      allocate (flxfac(ny, -nzgrid:nzgrid))

      flxfac = spread(delzed * dy, 1, ny) * jacob
      flxfac(:, -nzgrid) = 0.5 * flxfac(:, -nzgrid)
      flxfac(:, nzgrid) = 0.5 * flxfac(:, nzgrid)
      
      flxfac = flxfac/(sum(flxfac * grad_x) )

      call get_one_flux_ffs(dens, dphidy, flxfac, pflux_vs_s, pflux_kxkyzts(:, :, :, it, :))
      call get_one_flux_ffs(pres, dphidy, flxfac, qflux_vs_s, qflux_kxkyzts(:, :, :, it, :))
      call get_one_flux_ffs(upar, dphidy, flxfac, vflux_vs_s, vflux_kxkyzts(:, :, :, it, :))

      deallocate (dphidy, flxfac)

   end subroutine calculate_fluxes_fullfluxsurface

   !============================================================================
   !============================ CALCULATE FLUXES ==============================
   !============================================================================
   subroutine get_one_flux_ffs(mom, dphidy, flxfac, flx, flx_vs_kxkyz)

      use species, only: nspec
      use z_grid, only: nzgrid
      use parameters_kxky_grid, only: naky, nakx
      use volume_averages, only: mode_fac

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: mom
      complex, dimension(:, :, -nzgrid:), intent(in) :: dphidy
      real, dimension(:, -nzgrid:), intent(in) :: flxfac
      real, dimension(:), intent(out) :: flx
      real, dimension(:, :, -nzgrid:, :), intent(out) :: flx_vs_kxkyz

      integer :: iky, ikx, iz, is
      complex, dimension(:, :, :, :), allocatable :: mom_ky

      allocate (mom_ky(naky, nakx, -nzgrid:nzgrid, nspec))

      flx = 0.0
      flx_vs_kxkyz = 0.0

      !> Multiply input moment (e.g. density, momentum, or energy) but Jacobian (called flxfac)
      !> Need to do this because Jacobian has y-dependence 
      !> then Fourier transform back to (ky,kx)-space
      
      call get_modified_fourier_coefficient(mom, mom_ky, flxfac)
      do is = 1, nspec
         !> pflx_vs_kxkyz is the particle flux before summing over (kx,ky) and integrating over z
!flx_vs_kxkyz(:, :, :, is) = flxfac * aimag(mom_ky(:, :, :, is) * dphidy)
         flx_vs_kxkyz(:, :, :, is) = aimag(mom_ky(:, :, :, is) * conjg(dphidy(:, :, :)))
         !> calculate the volume average of the particle flux
         !> note that the factor of 1/B that appears in the Jacobian has already been taken into account
         !> in the numerator of the flux surface average
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               do iky = 1, naky
                  !if(iky == ikx == 1) then
                  !   cycle 
                  !else
                  flx(is) = flx(is) + 0.5 * mode_fac(iky) * flx_vs_kxkyz(iky, ikx, iz, is)
                  !end if
               end do
            end do
         end do
      end do

      deallocate (mom_ky)

   end subroutine get_one_flux_ffs

   !============================================================================
   !============================ CALCULATE FLUXES ==============================
   !============================================================================
   subroutine get_modified_fourier_coefficient(moment, moment_ky, flxfac)

      use species, only: nspec
      use z_grid, only: nzgrid
      use parameters_kxky_grid, only: ikx_max, naky_all, ny
      use calculations_kxky, only: swap_kxky_back
      use stella_transforms, only: transform_y2ky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: moment
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: moment_ky

      integer :: ikx, iz, is
      complex, dimension(:, :), allocatable :: tmp_kykx
      complex, dimension(:, :), allocatable :: tmp_ykx
      real, dimension(:, -nzgrid:), intent(in) :: flxfac

      integer :: is_end

      allocate (tmp_kykx(naky_all, ikx_max))
      allocate (tmp_ykx(ny, ikx_max))
      !> is_end is either nspec for fluxes, or ntubes when calcualting |phi|^2
      is_end = size(moment, dim=4)

      do is = 1, is_end
         do iz = -nzgrid, nzgrid
            do ikx = 1, ikx_max
               !> divide the input moment by the magnetic field strength
               !> to account for Jacobian in flux-surface average
               tmp_ykx(:, ikx) = moment(:, ikx, iz, is) * flxfac(:, iz)
            end do
            !> transform the B-modified input moment from y to ky space
            call transform_y2ky(tmp_ykx, tmp_kykx)
            !> swap from all ky and kx >= 0 to all kx and ky >= 0
            call swap_kxky_back(tmp_kykx, moment_ky(:, :, iz, is))
         end do
      end do
      deallocate (tmp_kykx, tmp_ykx)

   end subroutine get_modified_fourier_coefficient

   subroutine calculate_moments_fullfluxsurface(g, dens, upar, pres)

      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use species, only: spec, nspec
      use z_grid, only: nzgrid
      use velocity_grids, only: integrate_vmu_ffs
      use velocity_grids, only: vpa, vperp2
      use parameters_kxky_grid, only: naky_all, ikx_max, ny
      use calculations_kxky, only: swap_kxky
      use arrays_dist_fn, only: g0, g1, g2
      use gyro_averages, only: gyro_average, j0_ffs
      use arrays_fields, only: phi
      use stella_transforms, only: transform_ky2y
      use grids_kxky, only: aky, theta0
      use parameters_kxky_grid, only: nakx
      use constants, only: zi, pi
      use z_grid, only: ntubes
      use g_tofrom_h, only: g_to_h

      !> For momentum flux
      use gyro_averages, only: j1_ffs
      use geometry, only: gds21, gds22, gds2
      use geometry, only: geo_surf
      use geometry, only: gradzeta_grady_R2overB2, gradzeta_gradx_R2overB2, b_dot_grad_zeta_RR

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dens, upar, pres

      real, dimension(:), allocatable :: dens_wgts, upar_wgts, pres_wgts
      !> f_swap will contain delta f(ky,kx) on a grid with all kys and kx >= 0
      complex, dimension(:, :), allocatable :: f_swap
      !> fy will contain delta f(y,kx) on a grid with kx >= 0
      complex, dimension(:, :, :), allocatable :: fy, f2y, f3y 
      !> integrand will contain the integrand in the velocity moment integrals
      complex, dimension(:), allocatable :: integrand

      integer :: iy, ikx, iz, it
      integer :: ivmu, iv, imu, is
      real :: fac1, fac2

      !> species-dependent factor by which velocity moments must be multiplied
      !> to get density, pressure, etc.
      allocate (dens_wgts(nspec))
      allocate (upar_wgts(nspec))
      allocate (pres_wgts(nspec))

      dens = 0.; upar = 0.; pres = 0.
      !> set species-dependent factors needed for density, parallel flow and pressure
      dens_wgts = spec%dens
      pres_wgts = spec%dens * spec%temp
      upar_wgts = spec%stm !spec%dens * sqrt(spec%mass * spec%temp)

      !> Already allocated arrays (allocated in time_advance.f90). Set to zero just in case
      g0 = 0. ; g1 = 0. ; g2 = 0.

      g0 = g
	   !> This returns f = g + Z/T * (<phi> - phi)
      !> NB: g_to_f0 TODO-GA: check 
      call g_to_f1(g, phi, g0)

      !> g1 = J0 * f
      call gyro_average(g0, g1, j0_ffs)
      !> g2 = J1 * f
      call gyro_average(g0, g2, j1_ffs)
      
      allocate (f_swap(naky_all, ikx_max))
      allocate (integrand(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      allocate (fy(ny, ikx_max, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); fy = 0.
      allocate (f2y(ny, ikx_max, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); f2y = 0.
      allocate (f3y(ny, ikx_max, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); f3y = 0.

      !> assume only a single flux annulus
      it = 1
      do iz = -nzgrid, nzgrid
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            !> for every (z,vpa,mu,spec) point
            !> switch from ky >= 0 and kx = [-kxmax, kxmax]
            !> to ky = [-kymax, kymax] and kx >= 0
            call swap_kxky(g1(:, :, iz, it, ivmu), f_swap)
            !> for every (z,vpa,mu,spec) point, Fourier tranform from ky to y space to get
            !> the kx component of <f(y,x)>_r
            call transform_ky2y(f_swap, fy(:, :, ivmu))

            !! The following are only needed for the momentum flux because we need J1*f 
            !! J1* zi * ky * f
            g2(:, :, iz, it, ivmu) = zi * g2(:, :, iz, it, ivmu) * spread(aky, 2, nakx)
            call swap_kxky(g2(:, :, iz, it, ivmu), f_swap)
            call transform_ky2y(f_swap, f2y(:, :, ivmu))
            !! J1 * zi * kx * f
            g2(:, :, iz, it, ivmu) = g2(:, :, iz, it, ivmu) * theta0(:, :)
            call swap_kxky(g2(:, :, iz, it, ivmu), f_swap)
            call transform_ky2y(f_swap, f3y(:, :, ivmu))
         end do

         do ikx = 1, ikx_max
            do iy = 1, ny
               !> the integrand for the density moment is the distribution function
               integrand = fy(iy, ikx, :)
               !> integrate over v-space to get the density, normalised by the reference density.
               call integrate_vmu_ffs(integrand, dens_wgts, iy, iz, dens(iy, ikx, iz, :))

               !> the integrand for the pressure moment is the energy-weighted distribution function
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                  iv = iv_idx(vmu_lo, ivmu)
                  imu = imu_idx(vmu_lo, ivmu)
                  integrand(ivmu) = fy(iy, ikx, ivmu) * (vpa(iv)**2 + vperp2(iy, iz, imu))
               end do
               !> integrate over v-space to get the pressure, normalised by the reference pressure.
               call integrate_vmu_ffs(integrand, pres_wgts, iy, iz, pres(iy, ikx, iz, :))

               fac1 = gradzeta_grady_R2overB2(iy, iz) * gds21(iy, iz) / geo_surf%shat - gradzeta_gradx_R2overB2(iy, iz) * gds2(iy, iz)
               fac2 = gradzeta_grady_R2overB2(iy, iz) * gds22(iy, iz) / geo_surf%shat - gradzeta_gradx_R2overB2(iy, iz) * gds21(iy, iz)
               !> the integrand for the parallel flow moment is the parallel velocity
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                  iv = iv_idx(vmu_lo, ivmu)
                  is = is_idx(vmu_lo, ivmu)
                  imu = imu_idx(vmu_lo, ivmu)
                  integrand(ivmu) = fy(iy, ikx, ivmu) * b_dot_grad_zeta_RR(iy, iz) &
                                    - vperp2(iy, iz, imu) * spec(is)%smz * (f2y(iy, ikx, ivmu) * fac1 + f3y(iy, ikx, ivmu) * fac2)
               end do
               !> integrate over v-space to get the parallel flow, normalised by the reference thermal speed.
               call integrate_vmu_ffs(integrand, upar_wgts, iy, iz, upar(iy, ikx, iz, :))
            end do
         end do
      end do

      deallocate (dens_wgts, upar_wgts, pres_wgts)
      deallocate (f_swap)
      deallocate (integrand)
      deallocate (fy, f2y, f3y)

   end subroutine calculate_moments_fullfluxsurface

   !> the Fourier components of the guiding centre distribution function
   !> normalized by the equilibrium Maxwellian is passed in as g,
   !> along with the Fourier components of the electrostatic potential, phi.
   !> g_to_f calculates the Maxwellian-normalized distribution function f,
   !> which is related to g via
   !> f = g + (Ze/T)*(<phi>_R - phi)
   subroutine g_to_f0(g, phi, f)

      use stella_layouts, only: vmu_lo, is_idx
      use species, only: spec
      use z_grid, only: nzgrid, ntubes
      use gyro_averages, only: gyro_average, j0_ffs

      use stella_transforms, only: transform_ky2y, transform_y2ky
      use calculations_kxky, only: swap_kxky, swap_kxky_back
      use velocity_grids, only: maxwell_vpa, maxwell_mu
      use parameters_kxky_grid, only: naky, naky_all, nakx, ikx_max, ny
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: f
      complex, dimension(:, :), allocatable :: phi_swap
      complex, dimension(:, :), allocatable :: phiy
      complex, dimension(:, :, :), allocatable :: adjust

      complex, dimension (:,:,:,:), allocatable :: g_store
      integer :: ivmu, is, it
      integer :: iz, iv, imu, ia

      allocate (phi_swap(naky_all, ikx_max))
      allocate (phiy(ny, ikx_max))
      allocate (adjust(naky, nakx, -nzgrid:nzgrid))

      allocate(g_store(naky, nakx, -nzgrid:nzgrid, ntubes)) ; g_store = 0.0

      it = 1

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
	      !> compute phi * maxwellian in real space and transform back to k-space
         do iz = -nzgrid, nzgrid
!!!            call swap_kxky(adjust(:, :, iz), phi_swap)
            call swap_kxky(phi(:, :, iz,1), phi_swap)
            call transform_ky2y(phi_swap, phiy(:, :))

            !> phiy = maxwellian * (<phi> - phi) 
            do ia = 1, ny
               phiy(ia, :) = phiy(ia, :) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is)
            end do

            !> Put this back into k-space for future calculations
            call transform_y2ky(phiy(:, :), phi_swap)
            call swap_kxky_back(phi_swap, adjust(:, :, iz))
         end do

         !> calculate the normalized f, given phi and <phi>_R (temporarily stored in f)
         !! f = g + Z/T * (<phi> - phi) * maxwellian 
         is = is_idx(vmu_lo, ivmu)
         f(:, :, :, :, ivmu) = g(:, :, :, :, ivmu) + spec(is)%zt * spread(adjust(:, :, :), 4, ntubes)
         call gyro_average(f(:, :, :, :, ivmu), g_store(:,:,:,:), j0_ffs(:, :, :, ivmu)) 

         f(:, :, :, :, ivmu) = g_store - spec(is)%zt * spread(adjust(:, :, :), 4, ntubes)
      end do

      deallocate(adjust, phi_swap, phiy)
      deallocate(g_store) 
      
   end subroutine g_to_f0

   subroutine g_to_f1(g, phi, f)

      use stella_layouts, only: vmu_lo, is_idx
      use species, only: spec
      use z_grid, only: nzgrid, ntubes
      use gyro_averages, only: gyro_average, j0_ffs

      use stella_transforms, only: transform_ky2y, transform_y2ky
      use calculations_kxky, only: swap_kxky, swap_kxky_back
      use velocity_grids, only: maxwell_vpa, maxwell_mu
      use parameters_kxky_grid, only: naky, naky_all, nakx, ikx_max, ny
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: f
      complex, dimension(:, :), allocatable :: phi_swap
      complex, dimension(:, :), allocatable :: phiy
      complex, dimension(:, :, :), allocatable :: adjust

      complex, dimension (:,:,:,:), allocatable :: g_store
      integer :: ivmu, is, it
      integer :: iz, iv, imu, ia

      allocate (phi_swap(naky_all, ikx_max))
      allocate (phiy(ny, ikx_max))
      allocate (adjust(naky, nakx, -nzgrid:nzgrid))

      allocate(g_store(naky, nakx, -nzgrid:nzgrid, ntubes)) ; g_store = 0.0

      it = 1

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         call gyro_average(phi, g_store, j0_ffs(:, :, :, ivmu)) 
         adjust = g_store(:, :, :, 1) - phi(:, :, :, 1)
	      !> compute phi * maxwellian in real space and transform back to k-space
         do iz = -nzgrid, nzgrid
            call swap_kxky(adjust(:, :, iz), phi_swap)
            !! call swap_kxky(phi(:, :, iz,1), phi_swap)
            call transform_ky2y(phi_swap, phiy(:, :))

            !> phiy = maxwellian * phi 
            do ia = 1, ny
               phiy(ia, :) = phiy(ia, :) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is)
            end do

            !> Put this back into k-space for future calculations
            call transform_y2ky(phiy(:, :), phi_swap)
            call swap_kxky_back(phi_swap, adjust(:, :, iz))
         end do
         !> Note that adjust here is (<phi>_R - phi) * maxwellian 
         !> calculate the normalized f, given (<phi>_R - phi) 
         !> f = g + Z/T * (<phi> - phi) * maxwellian 
         f(:, :, :, :, ivmu) = g(:, :, :, :, ivmu) + spec(is)%zt * spread(adjust(:, :, :), 4, ntubes)
         is = is_idx(vmu_lo, ivmu)
         !> f = g + J0 * phi * maxwellian
      end do

      deallocate(adjust, phi_swap, phiy)
      deallocate(g_store) 
   end subroutine g_to_f1

end module diagnostics_fluxes_fullfluxsurface
