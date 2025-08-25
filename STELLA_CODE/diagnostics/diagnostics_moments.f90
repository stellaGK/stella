
!###############################################################################
!############################### DIAGNOSE MOMENTS ##############################
!###############################################################################
! 
! Routines for calculating and writing the moments. 
! 
! The dens_vs_kykxzts is denoted by dens.
! The temp_vs_kykxzts is denoted by temp.
! The parallel velocity is denoted by upar.
! 
!###############################################################################
 
module diagnostics_moments

   implicit none
 
   public :: init_diagnostics_moments 
   public :: write_moments_to_netcdf_file 

   private 

   ! Debugging
   logical :: debug = .false.

contains

!###############################################################################
!################################ WRITE MOMENTS ################################
!###############################################################################

   !============================================================================
   !================= CALCULATE AND WRITE MOMENTS TO NETCDF FILE ===============
   !============================================================================
   subroutine write_moments_to_netcdf_file(nout, timer)

      ! Data
      use arrays_store_distribution_fn, only: gnew

      ! Dimensions
      use parameters_kxky_grid, only: naky, nakx
      use grids_z, only: nztot, ntubes
      use grids_species, only: nspec
      
      ! Flags 
      use parameters_physics, only: radial_variation
      use parameters_physics, only: full_flux_surface

      ! Write to netcdf file 
      use stella_io, only: write_radial_moments_nc
      use stella_io, only: write_moments_nc
      
      ! Routines
      use job_manage, only: time_message
      use mp, only: proc0
      
      ! Input file
      use parameters_diagnostics, only: write_radial_moments
      use parameters_diagnostics, only: write_moments

      implicit none 

      ! The pointer in the netcdf file and a timer
      real, dimension(:), intent(in out) :: timer   
      integer, intent(in) :: nout    

      ! Variables needed to write and calculate diagnostics 
      complex, dimension(:, :, :, :, :), allocatable :: dens_vs_kykxzts, upar_vs_kykxzts, temp_vs_kykxzts, spitzer2_vs_kykxzts 
      real, dimension(:, :), allocatable :: dens_kxs, upar_kxs, temp_kxs 

      !---------------------------------------------------------------------- 

      ! Only continue if the moments have to be written
      if ((.not. write_moments) .and. (.not. write_radial_moments)) return  

      ! Start timer
      if (proc0) call time_message(.false., timer(:), 'Write moments')
      
      ! Allocate the arrays for the moments
      allocate (dens_vs_kykxzts(naky, nakx, nztot, ntubes, nspec))
      allocate (upar_vs_kykxzts(naky, nakx, nztot, ntubes, nspec))
      allocate (temp_vs_kykxzts(naky, nakx, nztot, ntubes, nspec))
      allocate (spitzer2_vs_kykxzts(naky, nakx, nztot, ntubes, nspec))
      if (write_radial_moments) allocate (dens_kxs(nakx, nspec))
      if (write_radial_moments) allocate (upar_kxs(nakx, nspec))
      if (write_radial_moments) allocate (temp_kxs(nakx, nspec)) 

      ! Calculate the moments delta n(ky,kx,z,tube,s); delta T(ky,kx,z,tube,s); delta u_par(ky,kx,z,tube,s)
      if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_moments'

      ! Calculate the moments if <radial_variation> = True
      if (radial_variation .or. write_radial_moments) then 
         call get_moments_radial_variation(gnew, dens_vs_kykxzts, upar_vs_kykxzts, temp_vs_kykxzts, dens_kxs, upar_kxs, temp_kxs, spitzer2_vs_kykxzts)
      end if
      
      ! Calculate the moments if <full_flux_surface> = True
      if (full_flux_surface .and. write_moments) then  
         
         ! TODO-GA The moments for FFS are calculated in the fluxes routine
         ! Since the fluxes rely on the moments

      ! Calculate the moments for a flux tube simulation
      else if (write_moments) then
         call get_moments_fluxtube(gnew, dens_vs_kykxzts, upar_vs_kykxzts, temp_vs_kykxzts, spitzer2_vs_kykxzts)
      end if

      ! Write the moments to the netcdf file
      if (proc0 .and. write_moments) call write_moments_nc(nout, dens_vs_kykxzts, upar_vs_kykxzts, temp_vs_kykxzts, spitzer2_vs_kykxzts)
      if (proc0 .and. write_radial_moments) call write_radial_moments_nc(nout, dens_kxs, upar_kxs, temp_kxs) 

      ! Deallocate the arrays for the moments
      deallocate (dens_vs_kykxzts, upar_vs_kykxzts, temp_vs_kykxzts, spitzer2_vs_kykxzts)
      if (allocated(dens_kxs)) deallocate (dens_kxs)
      if (allocated(upar_kxs)) deallocate (upar_kxs)
      if (allocated(temp_kxs)) deallocate (temp_kxs)

       ! End timer
       if (proc0) call time_message(.false., timer(:), 'Write moments')
 
   end subroutine write_moments_to_netcdf_file
   

   !============================================================================
   !====================== GET MOMENTS FOR THE FLUX TUBE =======================
   !============================================================================
   ! h is gyrophase independent, but is in gyrocenter coordinates,
   ! so requires a J_0 to get to particle coordinates
   ! <f>_r = h J_0 - Ze*phi/T * F0
   ! g     = h     - Ze*<phi>_R/T * F0
   ! <f>_r = g J_0 + Ze*(J_0<phi>_R-phi)/T * F0
   !============================================================================
   subroutine get_moments_fluxtube(g, density, upar_vs_kykxzts, temperature, spitzer2_vs_kykxzts)

      use grids_z, only: nzgrid, ntubes
      use grids_species, only: spec, nspec
      use grids_velocity, only: vpa, vperp2, integrate_vmu
      use grids_velocity, only: maxwell_mu, ztmax, maxwell_fac, maxwell_vpa
      use parameters_kxky_grid, only: naky, nakx
      use calculations_kxky, only: multiply_by_rho
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use calculations_gyro_averages, only: gyro_average
      use arrays_gyro_averages, only: aj0x
      use arrays_store_fields, only: phi
      use parameters_physics, only: fphi
      use parameters_numerical, only: maxwellian_normalization
      use parameters_physics, only: radial_variation
      use calculations_transforms, only: transform_kx2x_unpadded
      
      ! Import temp arrays g1 and g2 with dimensions (nky, nkx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
      use arrays_store_distribution_fn, only: g_gyro => g1 
      use arrays_store_distribution_fn, only: integrand => g2 

      implicit none

		! The distribution function enters with dimensions (ky, kx, z, tube, ivmus)
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      
		! The moments are returned with dimensions (ky, kx, z, tube, s)
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: density, temperature, upar_vs_kykxzts, spitzer2_vs_kykxzts

		! Local variables
      integer :: ivmu, iv, imu, is, ia
      
      ! We only have one field line because <full_flux_surface> = .false.
      ia = 1
 
      !=========================================================================
      !                     TURBULENT DENSITY FLUCTUATIONS                     !
      !=========================================================================
		! The normalized turbulent density fluctuations are calculated as:
		!		<dens> = tilde{δn} = (δn_s / n_r) (a / rho_r) 
		! 		<dens> = tilde{n_s} * velocity_integral( g*J0 + (Zs/Ts)*phi*(J0^2 - 1) * exp(-E_s/T_s) )
		! We do this in the following steps
		! 		<g_gyro> = g*J0 = <g> * <aj0x(iky, ikx, iz, ivmu)>
		! 		<integrand> = (Zs/Ts)*phi*(J0^2 - 1) = <spec(is)%zt> * <phi> * (<aj0x>**2 - 1)
		! 		<integrand> = (Zs/Ts)*phi*(J0^2 - 1) * exp(-E_s/T_s) = <integrand> * <maxwell_vpa> * <maxwell_mu> * <maxwell_fac> 
		! 		<integrand> = g*J0 + (Zs/Ts)*phi*(J0^2 - 1) * exp(-E_s/T_s) = <g_gyro> + <integrand> 
      !=========================================================================
      
      ! Calculate <integrand> = g*J0 + (Zs/Ts)*phi*(J0^2 - 1) * exp(-E_s/T_s)
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         call gyro_average(g(:, :, :, :, ivmu), ivmu, g_gyro(:, :, :, :, ivmu))
         integrand(:, :, :, :, ivmu) = spread(aj0x(:, :, :, ivmu)**2 - 1.0, 4, ntubes) * spec(is)%zt * fphi * phi
         if (.not. maxwellian_normalization) then
            integrand(:, :, :, :, ivmu) = integrand(:, :, :, :, ivmu) * maxwell_vpa(iv, is) * &
                  spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx) * maxwell_fac(is), 4, ntubes)
         end if
         integrand(:, :, :, :, ivmu) = integrand(:, :, :, :, ivmu) + g_gyro(:, :, :, :, ivmu)
      end do
      
      ! Calculate <dens> = tilde{n_s} * velocity_integral( g*J0 + (Zs/Ts)*phi*(J0^2 - 1) * exp(-E_s/T_s) )
      call integrate_vmu(integrand, spec%dens_psi0, density)

      !=========================================================================
      !                   TURBULENT TEMPERATURE FLUCTUATIONS                   !
      !=========================================================================
		! Calculate the turbulent temperature fluctuations
		!		<temp> = tilde{δT} = (δT_s / T_r) (a / rho_r) 
		! 		<temp> = tilde{T_s} * velocity_integral( g*J0 + (Zs/Ts)*phi*(J0^2 - 1) * ... * exp(-E_s/T_s) )
      !=========================================================================
         
      ! Calculate <integrand> = g*J0 + (Zs/Ts)*phi*(J0^2 - 1) * ... * exp(-E_s/T_s)
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         if (maxwellian_normalization) then
            integrand(:, :, :, :, ivmu) = (g_gyro(:, :, :, :, ivmu) + spec(is)%zt &
                  * spread(aj0x(:, :, :, ivmu)**2 - 1.0, 4, ntubes) * phi * fphi) &
                  * (vpa(iv)**2 + spread(spread(spread(vperp2(1, :, imu), 1, naky), 2, nakx), 4, ntubes) - 1.5) / 1.5
         else
            integrand(:, :, :, :, ivmu) = (g_gyro(:, :, :, :, ivmu) + ztmax(iv, is) &
                  * spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx) &
                  * maxwell_fac(is) * (aj0x(:, :, :, ivmu)**2 - 1.0), 4, ntubes) * phi * fphi) &
                  * (vpa(iv)**2 + spread(spread(spread(vperp2(1, :, imu), 1, naky), 2, nakx), 4, ntubes) - 1.5) / 1.5
         end if
      end do
      
      ! Calculate <temp> = tilde{T_s}/tilde{n_s} * velocity_integral( g*J0 + (Zs/Ts)*phi*(J0^2 - 1) * ... * exp(-E_s/T_s) )
      call integrate_vmu(integrand, spec%temp_psi0 * spec%dens_psi0, temperature)
      
      !=========================================================================
      !                                 SPITZER                                !      
      !=========================================================================
      ! For Spitzer problem tests of the collision operator 
      !========================================================================= 
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         integrand(:,:,:,:,ivmu) = g(:,:,:,:,ivmu) * ( vpa(iv) * (vpa(iv)**2 + &
                  spread(spread(spread(vperp2(1,:,imu),1,naky),2,nakx),4,ntubes)) - 5./2. * vpa(iv) )
      end do
      call integrate_vmu(integrand, spec%stm, spitzer2_vs_kykxzts) ! AVB: stm is the thermal speed

      !=========================================================================
      !                                  UPAR                                  !
      !=========================================================================
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         integrand(:, :, :, :, ivmu) = vpa(iv) * g_gyro(:, :, :, :, ivmu)
      end do
      call integrate_vmu(integrand, spec%stm_psi0, upar_vs_kykxzts)

   end subroutine get_moments_fluxtube
 
   !============================================================================
   !==================== GET MOMENTS FOR RADIAL VARIATION =====================
   !============================================================================
   subroutine get_moments_radial_variation(g, dens, upar_vs_kykxzts, temp, dens_kxs, upar_kxs, temp_kxs, spitzer2_vs_kykxzts)

      use grids_z, only: nzgrid, ntubes
      use grids_species, only: spec, nspec
      use grids_velocity, only: integrate_vmu
      use grids_velocity, only: vpa, vperp2, mu
      use grids_velocity, only: maxwell_mu, ztmax, maxwell_fac, maxwell_vpa
      use parameters_kxky_grid, only: naky, nakx
      use grids_kxky, only: rho_d_clamped
      use calculations_kxky, only: multiply_by_rho
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use arrays_store_distribution_fn, only: g1, g2
      use arrays_store_useful, only: kperp2, dkperp2dr
      use geometry, only: bmag, dBdrho
      use geometry, only: dl_over_b, d_dl_over_b_drho
      use calculations_gyro_averages, only: gyro_average
      use arrays_gyro_averages, only: aj0x, aj1x
      use arrays_store_fields, only: phi, phi_corr_QN, phi_proj
      use parameters_physics, only: fphi
      use parameters_numerical, only: maxwellian_normalization
      use parameters_physics, only: radial_variation
      use calculations_transforms, only: transform_kx2x_unpadded
      
      ! Input file
      use parameters_diagnostics, only: write_radial_moments
      use parameters_diagnostics, only: write_moments

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: dens, upar_vs_kykxzts, temp, spitzer2_vs_kykxzts
      real, dimension(:, :), intent(out) :: dens_kxs, upar_kxs, temp_kxs

      complex, dimension(:, :), allocatable :: g0k, g1k, g1x
      real :: zero

      integer :: ivmu, iv, imu, is, ia
      integer :: iz, it

      if (radial_variation) then
         allocate (g0k(naky, nakx))
      end if
      if (write_radial_moments) then
         allocate (g1k(1, nakx))
         allocate (g1x(1, nakx))
      end if

      ! Hack below. Works since J0^2 - 1 and its derivative are zero at the origin
      zero = 100.*epsilon(0.)

      ! h is gyrophase independent, but is in gyrocenter coordinates,
      ! so requires a J_0 to get to particle coordinates
      ! <f>_r = h J_0 - Ze*phi/T * F0
      ! g     = h     - Ze*<phi>_R/T * F0
      ! <f>_r = g J_0 + Ze*(J_0<phi>_R-phi)/T * F0

      ! calculate the integrand appearing in the integral for the dens_vs_kykxzts
      ia = 1
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         ! obtain the gyro-average of g that appears in the dens_vs_kykxzts integral
         call gyro_average(g(:, :, :, :, ivmu), ivmu, g1(:, :, :, :, ivmu))
         ! FLAG -- AJ0X NEEDS DEALING WITH BELOW
         g2(:, :, :, :, ivmu) = spread(aj0x(:, :, :, ivmu)**2 - 1.0, 4, ntubes) * spec(is)%zt * fphi * phi
         if (.not. maxwellian_normalization) then
            g2(:, :, :, :, ivmu) = g2(:, :, :, :, ivmu) * maxwell_vpa(iv, is) * &
                                   spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx) * maxwell_fac(is), 4, ntubes)
         end if
         g2(:, :, :, :, ivmu) = g2(:, :, :, :, ivmu) + g1(:, :, :, :, ivmu)
         ! g2(:, :, :, :, ivmu) = g1(:, :, :, :, ivmu) + ztmax(iv, is) &
         !                        * spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx) &
         !                                 * maxwell_fac(is) * (aj0x(:, :, :, ivmu)**2 - 1.0), 4, ntubes) * fphi * phi

         if (radial_variation) then
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  !phi
                  g0k = ztmax(iv, is) * maxwell_mu(ia, iz, imu, is) &
                        * maxwell_fac(is) * (aj0x(:, :, iz, ivmu)**2 - 1.0) * fphi * phi(:, :, iz, it) &
                        * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                           - spec(is)%fprim + (dBdrho(iz) / bmag(ia, iz)) * (1.0 - 2.0 * mu(imu) * bmag(ia, iz)) &
                           - aj1x(:, :, iz, ivmu) * aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                           * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                           / (aj0x(:, :, iz, ivmu)**2 - 1.0 + zero))

                  !g
                  g0k = g0k + g1(:, :, iz, it, ivmu) &
                        * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                           * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                           + dBdrho(iz) / bmag(ia, iz))

                  call multiply_by_rho(g0k)

                  ! g0k(1,1) = 0.0

                  !phi QN
                  g0k = g0k + ztmax(iv, is) * maxwell_mu(ia, iz, imu, is) * fphi &
                        * maxwell_fac(is) * (aj0x(:, :, iz, ivmu)**2 - 1.0) * phi_corr_QN(:, :, iz, it)

                  g2(:, :, iz, it, ivmu) = g2(:, :, iz, it, ivmu) + g0k
               end do
            end do
         end if
      end do
      call integrate_vmu(g2, spec%dens_psi0, dens)

      if (write_radial_moments) then
         dens_kxs = 0.0
         do is = 1, nspec
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g1k(1, :) = dens(1, :, iz, it, is) - phi_proj(:, 1, it)
                  call transform_kx2x_unpadded(g1k, g1x)
                  dens_kxs(:, is) = dens_kxs(:, is) &
                                  + real(g1x(1, :) * (dl_over_b(ia, iz) + rho_d_clamped * d_dl_over_b_drho(ia, iz)))
               end do
            end do
         end do
         dens_kxs = dens_kxs / ntubes
      end if

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         if (maxwellian_normalization) then
            g2(:, :, :, :, ivmu) = (g1(:, :, :, :, ivmu) + spec(is)%zt &
                                    * spread(aj0x(:, :, :, ivmu)**2 - 1.0, 4, ntubes) * phi * fphi) &
                                   * (vpa(iv)**2 + spread(spread(spread(vperp2(1, :, imu), 1, naky), 2, nakx), 4, ntubes) - 1.5) / 1.5
         else
            g2(:, :, :, :, ivmu) = (g1(:, :, :, :, ivmu) + ztmax(iv, is) &
                                    * spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx) &
                                             * maxwell_fac(is) * (aj0x(:, :, :, ivmu)**2 - 1.0), 4, ntubes) * phi * fphi) &
                                   * (vpa(iv)**2 + spread(spread(spread(vperp2(1, :, imu), 1, naky), 2, nakx), 4, ntubes) - 1.5) / 1.5
         end if
         if (radial_variation) then
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  !phi
                  g0k = ztmax(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is) &
                        * (vpa(iv)**2 + vperp2(ia, iz, imu) - 1.5) / 1.5 &
                        * (aj0x(:, :, iz, ivmu)**2 - 1.0) * phi(:, :, iz, it) * fphi &
                        * (-spec(is)%tprim * (vpa(iv)**2 + vperp2(ia, iz, imu) - 2.5) &
                           - spec(is)%fprim + (dBdrho(iz) / bmag(ia, iz)) * (1.0 - 2.0 * mu(imu) * bmag(ia, iz)) &
                           - aj1x(:, :, iz, ivmu) * aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                           * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                           / (aj0x(:, :, iz, ivmu)**2 - 1.0 + zero) &
                           + 2.0 * mu(imu) * dBdrho(iz) / (vpa(iv)**2 + vperp2(ia, iz, imu) - 1.5))

                  !g
                  g0k = g0k + g1(:, :, iz, it, ivmu) * (vpa(iv)**2 + vperp2(ia, iz, imu) - 1.5) / 1.5 &
                        * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                           * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                           + dBdrho(iz) / bmag(ia, iz) &
                           + 2.0 * mu(imu) * dBdrho(iz) / (vpa(iv)**2 + vperp2(ia, iz, imu) - 1.5))

                  call multiply_by_rho(g0k)

                  !phi QN
                  g0k = g0k + fphi * ztmax(iv, is) * maxwell_mu(ia, iz, imu, is) &
                        * (vpa(iv)**2 + vperp2(ia, iz, imu) - 1.5) / 1.5 &
                        * maxwell_fac(is) * (aj0x(:, :, iz, ivmu)**2 - 1.0) * phi_corr_QN(:, :, iz, it)

                  g2(:, :, iz, it, ivmu) = g2(:, :, iz, it, ivmu) + g0k
               end do
            end do
         end if
      end do
      ! integrate to get dTs/Tr
      !    call integrate_vmu (g2, spec%temp, temp)
      call integrate_vmu(g2, spec%temp_psi0 * spec%dens_psi0, temp)

      if (write_radial_moments) then
         temp_kxs = 0.0
         do is = 1, nspec
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g1k(1, :) = temp(1, :, iz, it, is)
                  call transform_kx2x_unpadded(g1k, g1x)
                  temp_kxs(:, is) = temp_kxs(:, is) &
                                  + real(g1x(1, :) * (dl_over_b(ia, iz) + rho_d_clamped * d_dl_over_b_drho(ia, iz)))
               end do
            end do
         end do
         temp_kxs = temp_kxs / ntubes
      end if

      ! for Spitzer problem tests of the collision operator
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         g2(:,:,:,:,ivmu) = g(:,:,:,:,ivmu)*( vpa(iv)*(vpa(iv)**2+spread(spread(spread(vperp2(1,:,imu),1,naky),2,nakx),4,ntubes)) - 5./2.*vpa(iv) )
      end do
      call integrate_vmu(g2, spec%stm, spitzer2_vs_kykxzts) ! AVB: stm is the thermal speed

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         g2(:, :, :, :, ivmu) = vpa(iv) * g1(:, :, :, :, ivmu)
         if (radial_variation) then
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  !g
                  g0k = vpa(iv) * g1(:, :, iz, it, ivmu) &
                        * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                           * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                           + dBdrho(iz) / bmag(ia, iz))

                  call multiply_by_rho(g0k)

                  g2(:, :, iz, it, ivmu) = g2(:, :, iz, it, ivmu) + g0k
               end do
            end do
         end if
      end do
      call integrate_vmu(g2, spec%stm_psi0, upar_vs_kykxzts)
      if (write_radial_moments) then
         upar_kxs = 0.0
         do is = 1, nspec
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g1k(1, :) = upar_vs_kykxzts(1, :, iz, it, is)
                  call transform_kx2x_unpadded(g1k, g1x)
                  upar_kxs(:, is) = upar_kxs(:, is) &
                                  + real(g1x(1, :) * (dl_over_b(ia, iz) + rho_d_clamped * d_dl_over_b_drho(ia, iz)))
               end do
            end do
         end do
         upar_kxs = upar_kxs / ntubes
      end if

      if (allocated(g0k)) deallocate (g0k)
      if (allocated(g1k)) deallocate (g1k)
      if (allocated(g1x)) deallocate (g1x)

   end subroutine get_moments_radial_variation
   

!###############################################################################
!############################ INITALIZE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !======================== INITALIZE THE DIAGNOSTICS =========================
   !============================================================================  
   subroutine init_diagnostics_moments()  

      use mp, only: proc0

      implicit none 

      !----------------------------------------------------------------------
      
      ! Only debug on the first processor
      debug = debug .and. proc0

   end subroutine init_diagnostics_moments 

end module diagnostics_moments

