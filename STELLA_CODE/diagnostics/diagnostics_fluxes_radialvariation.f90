
!###############################################################################
!#################### CALCULATE FLUXES FOR RADIAL VARIATION ####################
!###############################################################################
 
module diagnostics_fluxes_radialvariation

   implicit none
 
   public :: calculate_fluxes_radialvariation

   private     

   ! Debugging
   logical :: debug = .false.

   interface get_one_flux_vmulo
      module procedure get_one_flux_vmulo_int
      module procedure get_one_flux_vmulo_kxkyz
   end interface

contains

!###############################################################################
!############################### CALCULATE FLUXES ##############################
!###############################################################################

   !==============================================
   !============ GET FLUXES VMULO ================
   !==============================================
   subroutine calculate_fluxes_radialvariation(g, phi, pflux_vs_s, vflux_vs_s, qflux_vs_s, pflux_vs_kxs,  &
         vflux_vs_kxs, qflux_vs_kxs, pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts)
 
      use mp, only: sum_reduce
      use constants, only: zi
      use arrays_dist_fn, only: g1, g2, kperp2, dkperp2dr
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use species, only: spec
      use geometry, only: grho_norm, bmag, btor
      use geometry, only: drhodpsi
      use geometry, only: gds21, gds22
      use geometry, only: dgds21dr, dgds22dr
      use geometry, only: geo_surf
      use geometry, only: dBdrho, dIdrho
      use z_grid, only: nzgrid, ntubes
      use velocity_grids, only: vperp2, vpa, mu
      use velocity_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use numerical_parameters, only: fphi
      use numerical_parameters, only: maxwellian_normalization
      use grids_kxky, only: aky, theta0
      use kxky_grid_parameters, only: naky, nakx
      use calculations_kxky, only: multiply_by_rho
      use physics_parameters, only: radial_variation
      use gyro_averages, only: gyro_average, gyro_average_j1, aj0x, aj1x
      
      ! Flags 
      use parameters_diagnostics, only: write_radial_fluxes 
      use parameters_diagnostics, only: flux_norm

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      real, dimension(:), intent(out) :: pflux_vs_s, vflux_vs_s, qflux_vs_s
      real, dimension(:, :), intent(out) :: pflux_vs_kxs, vflux_vs_kxs, qflux_vs_kxs
      real, dimension(:, :, -nzgrid:, :, :), intent(out) :: pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts 

      integer :: ivmu, imu, iv, iz, it, is, ia
      real :: flx_norm
      complex, dimension(:, :), allocatable :: g0k, g1k
      
      ! Track the code
      if (debug) write (*, *) 'diagnostics::calculate_fluxes_radialvariation'

      pflux_vs_s = 0.; vflux_vs_s = 0.; qflux_vs_s = 0.
      pflux_vs_kxs = 0.; vflux_vs_kxs = 0.; qflux_vs_kxs = 0.
      pflux_kxkyzts = 0.; vflux_kxkyzts = 0.; qflux_kxkyzts = 0.

      ia = 1
      if (flux_norm) then
         flx_norm = 1./grho_norm
      else
         flx_norm = 1.
      end if

      allocate (g0k(naky, nakx))
      allocate (g1k(naky, nakx))

      ! FLAG - electrostatic for now
      ! get electrostatic contributions to fluxes

      if (fphi > epsilon(0.0)) then
         ia = 1

         !get particle flux
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            call gyro_average(g(:, :, :, :, ivmu), ivmu, g1(:, :, :, :, ivmu))

            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  if (radial_variation) then
                     g0k = g1(:, :, iz, it, ivmu) &
                           * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                              * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                              * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                              + dBdrho(iz) / bmag(ia, iz))

                     call multiply_by_rho(g0k)
                     g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) + g0k
                  end if

                  !subtract adiabatic contribution part of g
                  g0k = spec(is)%zt * fphi * phi(:, :, iz, it) * aj0x(:, :, iz, ivmu)**2
                  if (.not. maxwellian_normalization) then
                     g0k = g0k * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                  end if
                  if (radial_variation) then
                     g1k = g0k * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                                  - spec(is)%fprim - 2.0 * dBdrho(iz) * mu(imu) &
                                  - aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                                  * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                  * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                  + dBdrho(iz) / bmag(ia, iz))

                     call multiply_by_rho(g1k)

                     g0k = g0k + g1k
                  end if
                  g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) + g0k

               end do
            end do
         end do
         call get_one_flux_vmulo(flx_norm * spec%dens_psi0, g1, phi, pflux_vs_s)
         call get_one_flux_vmulo(flx_norm * spec%dens_psi0, g1, phi, pflux_kxkyzts)

         if (write_radial_fluxes) then
            call get_one_flux_radial(flx_norm * spec%dens_psi0, g1, phi, pflux_vs_kxs)
         end if

         !get heat flux
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            call gyro_average(g(:, :, :, :, ivmu), ivmu, g1(:, :, :, :, ivmu))
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid

                  g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) * (vpa(iv)**2 + vperp2(ia, iz, imu))

                  if (radial_variation) then
                     g0k = g1(:, :, iz, it, ivmu) &
                           * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                              * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                              * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                              + dBdrho(iz) / bmag(ia, iz) &
                              + 2.0 * mu(imu) * dBdrho(iz) / (vpa(iv)**2 + vperp2(ia, iz, imu)))

                     call multiply_by_rho(g0k)

                     g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) + g0k

                  end if

                  !subtract adiabatic contribution part of g
                  g0k = spec(is)%zt * fphi * phi(:, :, iz, it) * aj0x(:, :, iz, ivmu)**2 &
                        * (vpa(iv)**2 + vperp2(ia, iz, imu))
                  if (.not. maxwellian_normalization) then
                     g0k = g0k * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                  end if
                  if (radial_variation) then
                     g1k = g0k * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                                  - spec(is)%fprim - 2.0 * dBdrho(iz) * mu(imu) &
                                  - aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                                  * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                  * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                  + dBdrho(iz) / bmag(ia, iz) &
                                  + 2.0 * mu(imu) * dBdrho(iz) / (vpa(iv)**2 + vperp2(ia, iz, imu)))

                     call multiply_by_rho(g1k)

                     g0k = g0k + g1k
                  end if
                  g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) + g0k
               end do
            end do
         end do
         call get_one_flux_vmulo(flx_norm * spec%dens_psi0 * spec%temp_psi0, g1, phi, qflux_vs_s)
         call get_one_flux_vmulo(flx_norm * spec%dens_psi0 * spec%temp_psi0, g1, phi, qflux_kxkyzts)

         if (write_radial_fluxes) then
            call get_one_flux_radial(flx_norm * spec%dens_psi0 * spec%temp_psi0, g1, phi, qflux_vs_kxs)
         end if

         ! get momentum flux
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  ! parallel component
                  g0k = g(:, :, iz, it, ivmu) * vpa(iv) * geo_surf%rmaj * btor(iz) / bmag(ia, iz)
                  call gyro_average(g0k, iz, ivmu, g1(:, :, iz, it, ivmu))

                  if (radial_variation) then
                     g0k = g1(:, :, iz, it, ivmu) &
                           * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                              * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                              * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                              + dIdrho / (geo_surf%rmaj * btor(iz)))

                     call multiply_by_rho(g0k)

                     g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) + g0k

                  end if
                  !subtract adiabatic contribution part of g
                  g0k = spec(is)%zt * fphi * phi(:, :, iz, it) * aj0x(:, :, iz, ivmu)**2 &
                        * vpa(iv) * geo_surf%rmaj * btor(iz) / bmag(ia, iz)
                  if (.not. maxwellian_normalization) then
                     g0k = g0k * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                  end if
                  if (radial_variation) then
                     g1k = g0k * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                                  - spec(is)%fprim - 2.0 * dBdrho(iz) * mu(imu) &
                                  - aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                                  * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                  * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                  + dIdrho / (geo_surf%rmaj * btor(iz)))

                     call multiply_by_rho(g1k)

                     g0k = g0k + g1k
                  end if
                  g1(:, :, iz, it, ivmu) = g1(:, :, iz, it, ivmu) + g0k

                  ! perpendicular component
                  g0k = -g(:, :, iz, it, ivmu) * zi * spread(aky, 2, nakx) * vperp2(ia, iz, imu) * geo_surf%rhoc &
                        * (gds21(ia, iz) + theta0 * gds22(ia, iz)) * spec(is)%smz &
                        / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2)

                  call gyro_average_j1(g0k, iz, ivmu, g2(:, :, iz, it, ivmu))
                  if (radial_variation) then
                     g0k = -g(:, :, iz, it, ivmu) * zi * spread(aky, 2, nakx) * vperp2(ia, iz, imu) * geo_surf%rhoc &
                           * (dgds21dr(ia, iz) + theta0 * dgds22dr(ia, iz)) * aj1x(:, :, iz, ivmu) * spec(is)%smz &
                           / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2)

                     g0k = g0k - g(:, :, iz, it, ivmu) * zi * spread(aky, 2, nakx) * vperp2(ia, iz, imu) * geo_surf%rhoc &
                           * (gds21(ia, iz) + theta0 * gds22(ia, iz)) * spec(is)%smz &
                           / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2) &
                           * (0.5 * aj0x(:, :, iz, ivmu) - aj1x(:, :, iz, ivmu)) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz))

                     g0k = g0k + g2(:, :, iz, it, ivmu) &
                           * (-geo_surf%d2qdr2 * geo_surf%rhoc / (geo_surf%shat * geo_surf%qinp) &
                              - geo_surf%d2psidr2 * drhodpsi)

                     call multiply_by_rho(g0k)

                     g2(:, :, iz, it, ivmu) = g2(:, :, iz, it, ivmu) + g0k
                  end if

                  !subtract adiabatic contribution part of g
                  g0k = -spec(is)%zt * fphi * phi(:, :, iz, it) * aj0x(:, :, iz, ivmu) * aj1x(:, :, iz, ivmu) &
                        * zi * spread(aky, 2, nakx) * vperp2(ia, iz, imu) * geo_surf%rhoc &
                        * (gds21(ia, iz) + theta0 * gds22(ia, iz)) * spec(is)%smz &
                        / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2)
                  if (.not. maxwellian_normalization) then
                     g0k = g0k * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                  end if

                  if (radial_variation) then
                     g1k = -spec(is)%zt * fphi * phi(:, :, iz, it) * aj0x(:, :, iz, ivmu) * aj1x(:, :, iz, ivmu) &
                           * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is) &
                           * zi * spread(aky, 2, nakx) * vperp2(ia, iz, imu) * geo_surf%rhoc &
                           * (dgds21dr(ia, iz) + theta0 * dgds22dr(ia, iz)) * spec(is)%smz &
                           / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2)

                     g1k = g1k - spec(is)%zt * fphi * phi(:, :, iz, it) * aj0x(:, :, iz, ivmu) &
                           * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is) &
                           * zi * spread(aky, 2, nakx) * vperp2(ia, iz, imu) * geo_surf%rhoc &
                           * (gds21(ia, iz) + theta0 * gds22(ia, iz)) * spec(is)%smz &
                           / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2) &
                           * (0.5 * aj0x(:, :, iz, ivmu) - aj1x(:, :, iz, ivmu)) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz))

                     g1k = g1k + &
                           g0k * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                                  - spec(is)%fprim - 2.0 * dBdrho(iz) * mu(imu) &
                                  - 0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                                  * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                  * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                  - geo_surf%d2qdr2 * geo_surf%rhoc / (geo_surf%shat * geo_surf%qinp) &
                                  - geo_surf%d2psidr2 * drhodpsi)

                     call multiply_by_rho(g1k)

                     g0k = g0k + g1k
                  end if
                  g2(:, :, iz, it, ivmu) = g2(:, :, iz, it, ivmu) + g0k
               end do
            end do
         end do

         g1 = g1 + g2
         call get_one_flux_vmulo(flx_norm * spec%dens_psi0 * sqrt(spec%mass * spec%temp_psi0), g1, phi, vflux_vs_s)
         call get_one_flux_vmulo(flx_norm * spec%dens_psi0 * sqrt(spec%mass * spec%temp_psi0), g1, phi, vflux_kxkyzts)

         if (write_radial_fluxes) then
            call get_one_flux_radial(flx_norm * spec%dens_psi0 * sqrt(spec%mass * spec%temp_psi0), g1, phi, vflux_vs_kxs)
         end if

      end if

      if (allocated(g0k)) deallocate (g0k)
      if (allocated(g1k)) deallocate (g1k)

   end subroutine calculate_fluxes_radialvariation

   !==============================================
   !============ GET ONE FLUX VMULO ==============
   !==============================================
   subroutine get_one_flux_vmulo_int(weights, gin, fld, flxout)

      use velocity_grids, only: integrate_vmu
      use stella_layouts, only: vmu_lo
      use grids_kxky, only: aky
      use parameters_multibox, only: boundary_size
      use kxky_grid_parameters, only: nakx, naky
      use z_grid, only: nzgrid, ntubes
      use species, only: nspec
      use volume_averages, only: mode_fac
      use geometry, only: dVolume
      use stella_transforms, only: transform_kx2x_unpadded
      use physics_parameters, only: radial_variation

      implicit none

      real, dimension(:), intent(in) :: weights
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: fld
      real, dimension(:), intent(in out) :: flxout

      complex, dimension(:, :, :, :, :), allocatable :: totals
      complex, dimension(:, :), allocatable :: g0x, g1x

      integer :: ia, is, it, iz, ikx
      real, dimension(nspec) :: flux_sum
      real :: volume, factor

      allocate (totals(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))

      ia = 1
      flux_sum = 0.

      ! The factor in front of all the flux definitions is <factor> = -sgn(psi_t)/2
      ! The constants which differ for each flux are gathered in <weights> and added in <integrate_vmu()>
      ! For the particle flux <weights> = ñ_s/<|\tilde{∇}ρ>_ζ = <flx_norm> * <spec%dens_psi0>
      ! For the heat flux <weights> = ñ_s*T̃_s/<|\tilde{∇}ρ>_ζ = <flx_norm> * <spec%dens_psi0> * <spec%temp_psi0>
      ! For the momentum flux <weights> = ñ_s*sqrt(m̃*T̃_s)/<|\tilde{∇}ρ>_ζ = <flx_norm> * <spec%dens_psi0> * sqrt(<spec%mass> * <spec%temp_psi0>)
      factor = - 0.5 ! The clebsch_factor is inside <dVolume>, since it is inside <jacob>

      call integrate_vmu(gin, weights, totals)
      if (radial_variation) then !do it in real-space
         allocate (g0x(naky, nakx))
         allocate (g1x(naky, nakx))
         do is = 1, nspec
            volume = 0.
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call transform_kx2x_unpadded(totals(:, :, iz, it, is), g0x)
                  call transform_kx2x_unpadded(fld(:, :, iz, it), g1x)
                  do ikx = boundary_size + 1, nakx - boundary_size
                     flux_sum(is) = flux_sum(is) + &
                                    sum(factor * mode_fac * aky * aimag(g0x(:, ikx) * conjg(g1x(:, ikx))) * dVolume(ia, ikx, iz))
                     volume = volume + dVolume(ia, ikx, iz)
                  end do
               end do
            end do
         end do
         deallocate (g0x, g1x)
      else
         do is = 1, nspec
            volume = 0.
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do ikx = 1, nakx
                     flux_sum(is) = flux_sum(is) + &
                                    sum(factor * mode_fac * aky * aimag(totals(:, ikx, iz, it, is) * conjg(fld(:, ikx, iz, it))) * dVolume(ia, 1, iz))
                  end do
                  volume = volume + dVolume(ia, 1, iz)
               end do
            end do
         end do
      end if

      flxout = flxout + flux_sum / volume

      deallocate (totals)

   end subroutine get_one_flux_vmulo_int

   subroutine get_one_flux_vmulo_kxkyz(weights, gin, fld, flxout)

      use velocity_grids, only: integrate_vmu
      use stella_layouts, only: vmu_lo
      use grids_kxky, only: aky
      use kxky_grid_parameters, only: nakx, naky
      use z_grid, only: nzgrid, ntubes
      use species, only: nspec
      use volume_averages, only: mode_fac

      implicit none

      real, dimension(:), intent(in) :: weights
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: fld
      real, dimension(:, :, -nzgrid:, :, :), intent(in out) :: flxout

      complex, dimension(:, :, :, :, :), allocatable :: totals

      integer :: ia, is, it, iz, ikx

      allocate (totals(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))

      ia = 1 
      call integrate_vmu(gin, weights, totals)
      do is = 1, nspec
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  flxout(:, ikx, iz, it, is) = 0.5 * mode_fac * aky * aimag(totals(:, ikx, iz, it, is) * conjg(fld(:, ikx, iz, it)))
               end do
            end do
         end do
      end do

      deallocate (totals)

   end subroutine get_one_flux_vmulo_kxkyz

   !==============================================
   !=========== GET ONE FLUX RADIAL ==============
   !==============================================
   subroutine get_one_flux_radial(weights, gin, fld, flxout)

      use velocity_grids, only: integrate_vmu
      use geometry, only: dVolume
      use stella_layouts, only: vmu_lo
      use grids_kxky, only: aky
      use kxky_grid_parameters, only: nakx, naky
      use z_grid, only: nzgrid, ntubes
      use species, only: nspec
      use volume_averages, only: mode_fac
      use stella_transforms, only: transform_kx2x_unpadded 

      implicit none

      real, dimension(:), intent(in) :: weights
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: fld
      real, dimension(:, :), intent(in out) :: flxout

      real, dimension(:), allocatable :: dV_rad
      complex, dimension(:, :, :, :, :), allocatable :: totals

      complex, dimension(:, :), allocatable :: g0x, g1x

      integer :: ia, is, it, iz, ikx

      allocate (dV_rad(nakx))
      allocate (g0x(naky, nakx))
      allocate (g1x(naky, nakx))
      allocate (totals(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))

      ia = 1

      dV_rad = sum(sum(dVolume, 3), 1) * ntubes

      ! NB: this returns the flux-surface-averaged radial fluxes. Keep in mind that the
      !     volume element in a flux-surface, dV, may not be uniform across surfaces, so
      !     one cannot simply sum across the radius here to get the total flux; rather, one
      !     would have to multiply by dV/V across the radius first
      call integrate_vmu(gin, weights, totals)
      do is = 1, nspec
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               call transform_kx2x_unpadded(totals(:, :, iz, it, is), g0x)
               call transform_kx2x_unpadded(fld(:, :, iz, it), g1x)
               do ikx = 1, nakx
                  flxout(ikx, is) = flxout(ikx, is) + sum(0.5 * mode_fac * aky * aimag(g0x(:, ikx) * conjg(g1x(:, ikx))) * dVolume(ia, ikx, iz) / dV_rad(ikx))
               end do
            end do
         end do
      end do

      deallocate (dV_rad, g0x, g1x, totals)

   end subroutine get_one_flux_radial


end module diagnostics_fluxes_radialvariation
