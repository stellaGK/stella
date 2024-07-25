
!###############################################################################
!####################### CALCULATE FLUXES IN A FLUXTUBE ########################
!###############################################################################
! 
! Routines for calculating and writing the turbulent fluxes. 
! 
! 
! DEFINITION DISTRIBUTION
! -----------------------
! δf = h - Zs*e/Ts * phi * F0
! h = <h>_gyroaverage
! g = <delta f>_gyroaverage 
! g = h - Zs*e/Ts * <phi>_gyroaverage * F0 
! δf = g + Zs*e/Ts * (<phi>_gyroaverage - phi) * F0   
! 
! 
! DEFINITIONS FLUXES
! ------------------ 
! The particle flux is the first order moment of the fluctuating distribution function, 
!     Gamma_s = δn_s * vec{upar}_s = int (δf) * vec{v} dvec{v}
! 
! The radial component of the particle flux, along ∇ψ, is defined as 
!     Gamma_s,radial = 1/<|∇ψ|>_ψ < int (vec{v_Es} . ∇ψ ) δf  dvec{v} >_ψ
! 
! In stella we calculate the flux surface averaged radial component of the particle flux
!     pflux(kx,ky,s) = <Re[Gamma_s,radial]>_ψ 
!        = -(1/C)(ñ_s/2)(k̃_y/<∇̃ρ>_zeta) (2*B̃/√π) int_0^inf dμ̃ int_-inf^inf dṽ int J*(Im[conj(phi)*δf*J0])*dz / int(J*dz) 
!        = -(1/C)(ñ_s/2)(k̃_y/<∇̃ρ>_zeta) * int(Im[conj(phi) * (2*B̃/√π) int_0^inf dμ̃ int_-inf^inf dṽ (δf*J0) ])*J*dz ) / int(J*dz) 
!        = sum_z (-(1/C)(k̃_y/2) * ñ_s * Im[conj(phi)*integrate_vmu(δf*J0)]*J*dz) / (sum_z J*dz) * 1/<|tilde{∇}ρ>_ζ
! 
! The expression for the heat flux and the momentum flux are similar and we can write,
!     flux(kx,ky,s) = sum_z ( -(1/C)*(k̃_y/2) * CONSTANT * Im[conj(phi)*integrate_vmu(VELOCITY_INTEGRAND)] * NORM ) 
!                   = CONSTANT * get_one_flux(norm, velocityintegrand_vs_vpamu, phi)
! 
! Here NORM is a function of <z> and takes care of the flux surface average (or field line average) + a factor 1/<|tilde{∇}ρ>_ζ
! The extra factor 1/<|tilde{∇}ρ>_ζ is included if <flux_norm> = True, otherwise it is not included.
!     NORM(iz) = J[iz]*delzed[iz]/(sum_z J*dz) * 1/<|tilde{∇}ρ>_ζ = < . >_ψ[iz] / <|tilde{∇}ρ>_ζ
! 
!            CONSTANT       VELOCITY_INTEGRAND
! pflux:       ñ_s                δf*J0 
! qflux:     ñ_s*T̃_s            δf*J0*v^2 
! vflux:   ñ_s*√(m̃_s*T̃_s)    δf*J0*vpa*R*Btor/B - i*δf*J1*ky*vperp^2*rho*(gds21+theta0*gds22)*(sqrt(T*m)/Z)/(q*shat*B^2))
!  
! Note that for pflux and qflux, we can replace δf by g or h in the integrand.
! This has been tested numerically, and gives the same fluxes.
!###############################################################################
 
module diagnose_fluxes_fluxtube

   implicit none
 
   public :: calculate_fluxes_fluxtube  

   private      

contains

!###############################################################################
!############################### CALCULATE FLUXES ##############################
!###############################################################################
 
   !============================================================================
   !====================== CALCULATE ONE FLUX(KX,KY,Z,S) =======================
   !============================================================================
   ! In stella we calculate the flux surface averaged radial component of the particle flux
   !     pflux(kx,ky,s) = <Re[Gamma_s,radial]>_ψ 
   !        = -sgn(psi_t)(ñ_s/2)(k̃_y/<∇̃ρ>_zeta) (2*B̃/√π) int_0^inf dμ̃ int_-inf^inf dṽ int J*(Im[conj(phi)*δf*J0])*dz / int(J*dz) 
   !        = -sgn(psi_t)(ñ_s/2)(k̃_y/<∇̃ρ>_zeta) * int(Im[conj(phi) * (2*B̃/√π) int_0^inf dμ̃ int_-inf^inf dṽ (δf*J0) ])*J*dz ) / int(J*dz) 
   !        = sum_z (-sgn(psi_t)(k̃_y/2) * ñ_s * Im[conj(phi)*integrate_vmu(δf*J0)]*J*dz) / (sum_z J*dz) * 1/<|tilde{∇}ρ>_ζ
   ! 
   ! The expression for the heat flux and the momentum flux are similar and we can write,
   !     flux(kx,ky,s) = sum_z ( -sgn(psi_t)*(k̃_y/2) * CONSTANT * Im[conj(phi)*integrate_vmu(VELOCITY_INTEGRAND)] * NORM ) 
   !                   = CONSTANT * get_one_flux(norm, velocityintegrand_vs_vpamu, phi)
   ! 
   ! Here NORM is a function of <z> and takes care of the flux surface average (or field line average) + a factor 1/<|tilde{∇}ρ>_ζ
   ! The extra factor 1/<|tilde{∇}ρ>_ζ is included if <flux_norm> = True, otherwise it is not included.
   !     NORM(iz) = J[iz]*delzed[iz]/(sum_z J*dz) * 1/<|tilde{∇}ρ>_ζ = < . >_ψ[iz] / <|tilde{∇}ρ>_ζ
   ! 
   !            CONSTANT       VELOCITY_INTEGRAND
   ! pflux:       ñ_s                δf*J0 
   ! qflux:     ñ_s*T̃_s            δf*J0*v^2 
   ! vflux:   ñ_s*√(m̃_s*T̃_s)    δf*J0*vpa*R*Btor/B - i*δf*J1*ky*vperp^2*rho*(gds21+theta0*gds22)*(sqrt(T*m)/Z)/(q*shat*B^2))
   ! 
   ! For <flux>(s) we integrate over all the z-points so we need to add the factor <fluxnorm_vs_z> = < . >_ψ / <∇̃ρ>_ψ
   ! For <flux>(ky,kx,z,t,s) we do not integrate over z so we only add the factor 1/<∇̃ρ>_ψ 
   !  
   ! This subroutine assumes that the trubulent component of the distribution function, δf, is passed in.
   ! We know that δf = h - Zs*e/Ts * phi * F0, assuming now that h = <h>_gyroaverage and defining g = <delta f>_gyroaverage 
   ! we have g = h - Zs*e/Ts * <phi>_gyroaverage * F0 and δf = g + Zs*e/Ts * (<phi>_gyroaverage - phi) * F0   
   ! It's been tested numerically, and whether we give <g>, <h> or <δf> does not make a difference for 
   ! <qflux> or <pflux>, but it does matter for <vflux>! Only <δf> is the correct options for <vflux> TODO is it? 
   !============================================================================ 
   subroutine calculate_fluxes_fluxtube(df_vs_vpamuikxkyzs, pflux_vs_s, vflux_vs_s, &
      qflux_vs_s, pflux_vs_kxkyzts, vflux_vs_kxkyzts, qflux_vs_kxkyzts, flux_norm, write_fluxes_kxkyz)

      ! Flags
      use physics_flags, only: include_apar, include_bpar

      ! Data 
      use fields_arrays, only: phi, apar, bpar
      use run_parameters, only: fphi, fapar
      
      ! Geometry 
      use stella_geometry, only: bmag, btor, gds2, gds21, gds22, geo_surf 
      use stella_geometry, only: gradzeta_gradx_RRoverBB
      use stella_geometry, only: gradzeta_grady_RRoverBB
      use stella_geometry, only: b_dot_grad_zeta_RR
      use stella_geometry, only: jacob, grho, btor 
      
      ! Dimensions
      use vpamu_grids, only: vperp2, vpa, mu
      use vpamu_grids, only: nvpa, nmu
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: aky, theta0
      use species, only: spec, nspec

      ! Calculations
      use gyro_averages, only: gyro_average, gyro_average_j1
      use constants, only: zi

      ! Routines
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx, kxkyz_lo
      use mp, only: sum_reduce

      implicit none

      ! We receive the distribution function <δf>(vpa,mu,i[kx,ky,z,s]) 
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: df_vs_vpamuikxkyzs

      ! Return the fluxes in the flux tube 
      real, dimension(:, :, -nzgrid:, :, :), intent(out) :: pflux_vs_kxkyzts, vflux_vs_kxkyzts, qflux_vs_kxkyzts
      real, dimension(:), intent(out) :: pflux_vs_s, vflux_vs_s, qflux_vs_s

      ! Toggle to include or not include the factor 1/<∇̃ρ>_ψ in the flux definition
      logical, intent(in) :: flux_norm, write_fluxes_kxkyz

      ! Local variables
      complex, dimension(:, :), allocatable :: velocityintegrand_vs_vpamu, temp1_vs_vpamu, temp2_vs_vpamu
      real, dimension(:), allocatable :: fluxnorm_vs_z
      integer :: ikxkyz, iky, ikx, iz, it, is, ia
      real :: one_over_nablarho

      ! Allocate arrays
      allocate (fluxnorm_vs_z(-nzgrid:nzgrid))
      allocate (velocityintegrand_vs_vpamu(nvpa, nmu), temp1_vs_vpamu(nvpa, nmu), temp2_vs_vpamu(nvpa, nmu))

      ! Make sure that the ararys we will fill are empty
      pflux_vs_s = 0.; vflux_vs_s = 0.; qflux_vs_s = 0.
      pflux_vs_kxkyzts = 0.; vflux_vs_kxkyzts = 0.; qflux_vs_kxkyzts = 0.

      ! Get <fluxnorm_vs_z> which takes care of the flux surface average, 
      ! and the factor <one_over_nablarho> which is used if <flux_norm> = True, otherwise its set to 1.
      call get_factor_for_fluxsurfaceaverage(fluxnorm_vs_z, one_over_nablarho, flux_norm) 

      ! We only have one flux tube since <radial_variation> = False
      ia = 1

      !===================================
      !           ELECTROSTATIC           
      !===================================
      ! Get electrostatic contributions to the fluxes
      ! TODO heat flux and particle flux should be the same with g or h, test both options!
      if (fphi > epsilon(0.0)) then
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc

            ! Get the (kx,ky,z,tube,s) indices on this processor
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! For the particle flux we have <VELOCITY_INTEGRAND>(vpa,mu) = h(mu,vpa)*J0 
            ! Add the contribution (-sgn(psi_t)*(k̃_y/2)*Im[conj(phi)*integrate_vmu(VELOCITY_INTEGRAND)]*NORM) to the particle flux 
            call gyro_average(df_vs_vpamuikxkyzs(:, :, ikxkyz), ikxkyz, velocityintegrand_vs_vpamu)
            call get_one_flux(iky, iz, fluxnorm_vs_z(iz), velocityintegrand_vs_vpamu, phi(iky, ikx, iz, it), pflux_vs_s(is))
            if (write_fluxes_kxkyz) call get_one_flux(iky, iz, one_over_nablarho, velocityintegrand_vs_vpamu, & 
                  phi(iky, ikx, iz, it), pflux_vs_kxkyzts(iky, ikx, iz, it, is))

            ! For the heat flux we have <VELOCITY_INTEGRAND>(vpa,mu) = h(mu,vpa)*J0*v^2
            ! Add the contribution (-sgn(psi_t)*(k̃_y/2)*Im[conj(phi)*integrate_vmu(VELOCITY_INTEGRAND)]*NORM) to the heat flux 
            velocityintegrand_vs_vpamu = velocityintegrand_vs_vpamu * (spread(vpa**2, 2, nmu) + spread(vperp2(1, iz, :), 1, nvpa))
            call get_one_flux(iky, iz, fluxnorm_vs_z(iz), velocityintegrand_vs_vpamu, phi(iky, ikx, iz, it), qflux_vs_s(is))
            if (write_fluxes_kxkyz) call get_one_flux(iky, iz, one_over_nablarho, velocityintegrand_vs_vpamu, &
                  phi(iky, ikx, iz, it), qflux_vs_kxkyzts(iky, ikx, iz, it, is))

            ! For the parallel (first term) and perpendicular (second term) component of the momentum flux we have,
            ! <VELOCITY_INTEGRAND>(vpa,mu) = δf*J0*vpa*(b.∇ζ)*R^2 
            !			- i*δf*J1*ky*vperp^2*rho*(gds21+theta0*gds22)*(∇ζ.∇y)*(sqrt(T*m)/Z)/(q*shat*B^2))
            ! 			+ i*δf*J1*ky*vperp^2*rho*(theta0*gds21+gds2)*(∇ζ.∇x)*(sqrt(T*m)/Z)/(q*B^2)) 
            
            ! First add the parallel component of the momentum flux: δf*J0*vpa*(b.∇ζ)*R^2 
            velocityintegrand_vs_vpamu = df_vs_vpamuikxkyzs(:, :, ikxkyz) * spread(vpa, 2, nmu) * b_dot_grad_zeta_RR(ia, iz)
            call gyro_average(velocityintegrand_vs_vpamu, ikxkyz, temp1_vs_vpamu)
            
            ! Next add the perpendicular component
            velocityintegrand_vs_vpamu = -df_vs_vpamuikxkyzs(:, :, ikxkyz) * spread(vperp2(ia, iz, :), 1, nvpa) * spec(is)%smz &
                    * zi * aky(iky) * (gradzeta_grady_RRoverBB(ia, iz) * (gds21(ia, iz) &
                    + theta0(iky, ikx) * gds22(ia, iz)) / geo_surf%shat &
                    - gradzeta_gradx_RRoverBB(ia, iz) * (theta0(iky, ikx) * gds21(ia, iz) + gds2(ia, iz)))
            call gyro_average_j1(velocityintegrand_vs_vpamu, ikxkyz, temp2_vs_vpamu)
            
            ! Sum parallel and perpendicular components together
            velocityintegrand_vs_vpamu = temp1_vs_vpamu + temp2_vs_vpamu
	
				! Add the contribution (-sgn(psi_t)*(k̃_y/2)*Im[conj(phi)*integrate_vmu(VELOCITY_INTEGRAND)]*NORM) to the momentum flux
            call get_one_flux(iky, iz, fluxnorm_vs_z(iz), velocityintegrand_vs_vpamu, phi(iky, ikx, iz, it), vflux_vs_s(is))
            if (write_fluxes_kxkyz) call get_one_flux(iky, iz, one_over_nablarho, velocityintegrand_vs_vpamu, phi(iky, ikx, iz, it), vflux_vs_kxkyzts(iky, ikx, iz, it, is))

         end do
      end if

      !===================================
      !          ELECTROMAGNETIC                    
      !===================================
      if (include_apar) then
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz) 

            ! Apar contribution to particle flux
            temp1_vs_vpamu = -2.0*df_vs_vpamuikxkyzs(:, :, ikxkyz) * spec(is)%stm * spread(vpa, 2, nmu)
            call gyro_average(temp1_vs_vpamu, ikxkyz, velocityintegrand_vs_vpamu)
            call get_one_flux(iky, iz, fluxnorm_vs_z(iz), velocityintegrand_vs_vpamu, apar(iky, ikx, iz, it), pflux_vs_s(is))
            if (write_fluxes_kxkyz) call get_one_flux(iky, iz, one_over_nablarho, velocityintegrand_vs_vpamu, & 
               apar(iky, ikx, iz, it), pflux_vs_kxkyzts(iky, ikx, iz, it, is))

            ! Apar contribution to heat flux
            velocityintegrand_vs_vpamu = velocityintegrand_vs_vpamu * (spread(vpa**2, 2, nmu) + spread(vperp2(ia, iz, :), 1, nvpa))
            call get_one_flux(iky, iz, fluxnorm_vs_z(iz), velocityintegrand_vs_vpamu, apar(iky, ikx, iz, it), qflux_vs_s(is))
            if (write_fluxes_kxkyz) call get_one_flux(iky, iz, one_over_nablarho, velocityintegrand_vs_vpamu, & 
               apar(iky, ikx, iz, it), qflux_vs_kxkyzts(iky, ikx, iz, it, is))

            ! TODO -- NEED TO ADD IN CONTRIBUTION FROM BOLTZMANN PIECE !! Only valid for axis-symmetric devices now
            ! Apar contribution to the parallel (first term) and perpendicular (second term) component of the momentum flux
            velocityintegrand_vs_vpamu = -2.0*spread(vpa**2, 2, nmu) * spec(is)%stm * df_vs_vpamuikxkyzs(:, :, ikxkyz) &
                    * geo_surf%rmaj * btor(iz) / bmag(1, iz)
            call gyro_average(velocityintegrand_vs_vpamu, ikxkyz, temp1_vs_vpamu)
            velocityintegrand_vs_vpamu = 2.0*spread(vpa, 2, nmu) * spec(is)%stm * df_vs_vpamuikxkyzs(:, :, ikxkyz) &
                    * zi * aky(iky) * spread(vperp2(ia, iz, :), 1, nvpa) * geo_surf%rhoc &
                    * (gds21(ia, iz) + theta0(iky, ikx) * gds22(ia, iz)) * spec(is)%smz &
                    / (geo_surf%qinp * geo_surf%shat * bmag(ia, iz)**2)
            call gyro_average_j1(velocityintegrand_vs_vpamu, ikxkyz, temp2_vs_vpamu)
            velocityintegrand_vs_vpamu = temp1_vs_vpamu + temp2_vs_vpamu
            call get_one_flux(iky, iz, fluxnorm_vs_z(iz), velocityintegrand_vs_vpamu, apar(iky, ikx, iz, it), vflux_vs_s(is))
            if (write_fluxes_kxkyz) call get_one_flux(iky, iz, one_over_nablarho, velocityintegrand_vs_vpamu, & 
               apar(iky, ikx, iz, it), vflux_vs_kxkyzts(iky, ikx, iz, it, is))

         end do
      end if 
      
      
      if (include_bpar) then
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! Bpar contribution to particle flux
            temp1_vs_vpamu = 4.0 * df_vs_vpamuikxkyzs(:, :, ikxkyz) * spec(is)%tz * spread(mu, 1, nvpa)
            call gyro_average_j1(temp1_vs_vpamu, ikxkyz, velocityintegrand_vs_vpamu)
            call get_one_flux(iky, iz, fluxnorm_vs_z(iz), velocityintegrand_vs_vpamu, bpar(iky, ikx, iz, it), pflux_vs_s(is))
            call get_one_flux(iky, iz, one_over_nablarho, velocityintegrand_vs_vpamu, bpar(iky, ikx, iz, it), pflux_vs_kxkyzts(iky, ikx, iz, it, is))

            ! Bpar contribution to heat flux
            velocityintegrand_vs_vpamu = velocityintegrand_vs_vpamu * (spread(vpa**2, 2, nmu) + spread(vperp2(ia, iz, :), 1, nvpa))
            call get_one_flux(iky, iz, fluxnorm_vs_z(iz), velocityintegrand_vs_vpamu, bpar(iky, ikx, iz, it), qflux_vs_s(is))
            call get_one_flux(iky, iz, one_over_nablarho, velocityintegrand_vs_vpamu, bpar(iky, ikx, iz, it), qflux_vs_kxkyzts(iky, ikx, iz, it, is))

            ! Bpar contribution to momentum flux
            ! NOT SUPPORTED, REQUIRES d J1(x)/ d x 
            velocityintegrand_vs_vpamu = 4.0 * spec(is)%tz * spread(mu, 1, nvpa) * df_vs_vpamuikxkyzs(:, :, ikxkyz) &
                    * spread(vpa, 2, nmu) * geo_surf%rmaj * btor(iz) / bmag(ia, iz)
            call gyro_average_j1(velocityintegrand_vs_vpamu, ikxkyz, temp1_vs_vpamu)
            temp2_vs_vpamu = 0.0
            velocityintegrand_vs_vpamu = temp1_vs_vpamu + temp2_vs_vpamu
            call get_one_flux(iky, iz, fluxnorm_vs_z(iz), velocityintegrand_vs_vpamu, bpar(iky, ikx, iz, it), vflux_vs_s(is)) 
            call get_one_flux(iky, iz, one_over_nablarho, velocityintegrand_vs_vpamu, bpar(iky, ikx, iz, it), vflux_vs_kxkyzts(iky, ikx, iz, it, is))
         end do
      end if

      ! Sum the values on all proccesors and send them to <proc0>
      call sum_reduce(pflux_vs_s, 0)
      call sum_reduce(qflux_vs_s, 0)
      call sum_reduce(vflux_vs_s, 0)

      ! Add the CONSTANT to the fluxes (ñ_s for pflux; ñ_s*T̃_s for qflux and ñ_s*√(m̃_s*T̃_s) for vflux) 
      pflux_vs_s = pflux_vs_s * spec%dens_psi0
      qflux_vs_s = qflux_vs_s * spec%dens_psi0 * spec%temp_psi0
      vflux_vs_s = vflux_vs_s * spec%dens_psi0 * sqrt(spec%mass * spec%temp_psi0)

      ! Normalise to account for contributions from multiple flux tubes in a flux tube train
      pflux_vs_s = pflux_vs_s / real(ntubes)
      qflux_vs_s = qflux_vs_s / real(ntubes)
      vflux_vs_s = vflux_vs_s / real(ntubes)

      ! Sum the values on all proccesors and send them to <proc0>
      call sum_reduce(pflux_vs_kxkyzts, 0)
      call sum_reduce(qflux_vs_kxkyzts, 0)
      call sum_reduce(vflux_vs_kxkyzts, 0)

      ! Add the CONSTANT to the fluxes (ñ_s for pflux; ñ_s*T̃_s for qflux and ñ_s*√(m̃_s*T̃_s) for vflux) 
      do is = 1, nspec
         pflux_vs_kxkyzts(:, :, :, :, is) = pflux_vs_kxkyzts(:, :, :, :, is) * spec(is)%dens_psi0
         qflux_vs_kxkyzts(:, :, :, :, is) = qflux_vs_kxkyzts(:, :, :, :, is) * spec(is)%dens_psi0 * spec(is)%temp_psi0
         vflux_vs_kxkyzts(:, :, :, :, is) = vflux_vs_kxkyzts(:, :, :, :, is) * spec(is)%dens_psi0 * sqrt(spec(is)%mass * spec(is)%temp_psi0)
      end do

      ! Deallocate the arrays
      deallocate (velocityintegrand_vs_vpamu, temp1_vs_vpamu, temp2_vs_vpamu)
      deallocate (fluxnorm_vs_z)

   end subroutine calculate_fluxes_fluxtube

   !============================================================================
   !==================== CALCULATE ONE FLUX(iKX,iKY,iZ,iS) =====================
   !============================================================================
   ! For example, the flux surface averaged particle flux in the radial direction is defined as,
   ! 
   !     pflux(kx,ky,s) = <Re[Gamma_s,radial]>_ψ 
   !        = -sgn(psi_t)(ñ_s/2)(k̃_y/<∇̃ρ>_zeta) (2*B̃/√π) int_0^inf dμ̃ int_-inf^inf dṽ int J*(Im[conj(phi)*δf*J0])*dz / int(J*dz) 
   !        = -sgn(psi_t)(ñ_s/2)(k̃_y/<∇̃ρ>_zeta) * int(Im[conj(phi) * (2*B̃/√π) int_0^inf dμ̃ int_-inf^inf dṽ (δf*J0) ])*J*dz ) / int(J*dz) 
   !        = sum_z (-sgn(psi_t)(k̃_y/2) * ñ_s * Im[conj(phi)*integrate_vmu(δf*J0)]*J*dz) / (sum_z J*dz) * 1/<|tilde{∇}ρ>_ζ
   ! 
   ! The expression for the heat flux and the momentum flux are similar and we can write,
   !     flux(kx,ky,s) = sum_z ( -sgn(psi_t)*(k̃_y/2) * CONSTANT * Im[conj(phi)*integrate_vmu(VELOCITY_INTEGRAND)] * NORM )
   ! 
   ! This routine will calculate flux(kx,ky,z) without the CONSTANT, as it can be added later
   !    flux(kx,ky,s) = CONSTANT * get_one_flux(norm, velocityintegrand_vs_vpamu, phi)
   ! 
   ! Here NORM is a function of <z> and takes care of the flux surface average (or field line average) + a factor 1/<|tilde{∇}ρ>_ζ
   ! The extra factor 1/<|tilde{∇}ρ>_ζ is included if <flux_norm> = True, otherwise it is not included.
   !     NORM(iz) = J[iz]*delzed[iz]/(sum_z J*dz) * 1/<|tilde{∇}ρ>_ζ = < . >_ψ[iz] / <|tilde{∇}ρ>_ζ
   ! 
   !            CONSTANT       VELOCITY_INTEGRAND
   ! pflux:       ñ_s                δf*J0 
   ! qflux:     ñ_s*T̃_s            δf*J0*v^2 
   ! vflux:   ñ_s*√(m̃_s*T̃_s)    δf*J0*vpa*R*Btor/B - i*δf*J1*ky*vperp^2*rho*(gds21+theta0*gds22)*(sqrt(T*m)/Z)/(q*shat*B^2))
   ! 
   ! If the goal is to integrate over z, then NORM = < . >_ψ[iz] / <∇̃ρ>_ψ if
   ! If we want to calculate the flux at every z-point, then NORM = 1 / <∇̃ρ>_ψ and the integral can be performed manually later
   !============================================================================ 
   subroutine get_one_flux(iky, iz, norm, velocityintegrand_vs_vpamu, phi, flux_out)

      use vpamu_grids, only: integrate_vmu
      use volume_averages, only: mode_fac
      use stella_geometry, only: flux_fac
      use kt_grids, only: aky

      implicit none
      
      ! We receive <velocityintegrand>(vpa,mu) and <phi> for a specific point i[kx,ky,z,tube,s]
      complex, dimension(:, :), intent(in) :: velocityintegrand_vs_vpamu
      complex, intent(in) :: phi

      ! We add each contribution (-sgn(psi_t)*(k̃_y/2)*Im[conj(phi)*integrate_vmu(VELOCITY_INTEGRAND)]*NORM) to <flux_out>  
      real, intent(in out) :: flux_out
      real, intent(in) :: norm

      ! The specific ky-point is needed to multiple with <ky>
      ! The specific z-point is needed since the integration weights for the velocity depend on <z>
      integer, intent(in) :: iky, iz

      ! We first calculate the velocity integral, and then the contribution to the flux
      complex :: velocity_integral

      ! Calculate the velocity integral <velocity_integral> = (2B̃/√π) int( int( (velocityintegrand(vpa,mu) ) dṽparallel) dμ̃)
      call integrate_vmu(velocityintegrand_vs_vpamu, iz, velocity_integral)

      ! Add the contribution (-(1/C)*(k̃_y/2)*Im[conj(phi)*VELOCITY_INTEGRAL]*NORM) to the flux 
      flux_out = flux_out + flux_fac * mode_fac(iky) * aky(iky) * aimag(velocity_integral * conjg(phi)) * norm

   end subroutine get_one_flux

   !============================================================================
   !=========================== FLUX SURFACE AVERAGE ===========================
   !============================================================================
   ! 
   ! The flux surface average reduces to the field line average in the flux tube approximation
   !     <Q>_ψ = <Q>_fluxsurface = int Q dV / int dV = int Q J dz / int J dz
   !     <Q>_ψ[iz] = Q[iz]*jacob[iz]*delzed[iz] / sum(jacob[iz]*delzed[iz])
   ! 
   ! Therefore the integration weights or the normalization factor for the fluxes is defined as 
   !     fluxnorm_vs_z[iz] = < . >_ψ = jacob[iz]*delzed[iz] / sum(jacob[iz]*delzed[iz])
   ! 
   ! In the definition of the fluxes we have a factor 1/<∇̃ρ>_fluxsurface which is calculated as 
   !      1/<∇̃ρ>_ψ = int dV / int ∇̃ρ dV = sum(jacob[iz]*delzed[iz]) / sum(grho[iz]*jacob[iz]*delzed[iz])
   !      grho = a|∇ρ0| = a (dρ0/dψ) ∇ψ = 1/ρ0 * ∇ψ/(a*Bref) = sqrt(|grad_psi_grad_psi|)/ρ0
   ! 
   ! If <flux_norm> = True then we absorb the factor 1/<∇̃ρ>_fluxsurface in the definition of <fluxnorm_vs_z> 
   !     fluxnorm_vs_z[iz] = < . >_ψ / <∇̃ρ>_ψ = jacob[iz]*delzed[iz] / sum(grho[iz]*jacob[iz]*delzed[iz]) 
   !============================================================================
   subroutine get_factor_for_fluxsurfaceaverage(fluxnorm_vs_z, one_over_nablarho, flux_norm) 

      use stella_geometry, only: jacob, grho
      use zgrid, only: delzed, nzgrid  

      implicit none

      ! Construct the factor in front of the flux definition 
      real, dimension(:), allocatable, intent(out) :: fluxnorm_vs_z 
      real, intent(out) :: one_over_nablarho

      ! Toggle to include or not include the factor 1/<∇̃ρ>_ψ in the flux definition
      logical, intent(in) :: flux_norm

      ! Local variables
      real, dimension(:), allocatable :: jacob_times_delzed_vs_z 

      ! Allocate array
      allocate (fluxnorm_vs_z(-nzgrid:nzgrid))
      allocate (jacob_times_delzed_vs_z(-nzgrid:nzgrid))

      ! Multiply the Jacobian with the step size in z and make each end count towards half the flux
      jacob_times_delzed_vs_z = jacob(1, :) * delzed
      jacob_times_delzed_vs_z(-nzgrid) = 0.5 * jacob_times_delzed_vs_z(-nzgrid)
      jacob_times_delzed_vs_z(nzgrid) = 0.5 * jacob_times_delzed_vs_z(nzgrid)

      ! Define the factor in front of the defintion of the fluxes 
      ! If <flux_norm> = False:  fluxnorm_vs_z[iz] = < . >_ψ = jacob[iz]*delzed[iz] / sum(jacob[iz]*delzed[iz])
      ! If <flux_norm> = True:   fluxnorm_vs_z[iz] = < . >_ψ / <∇̃ρ>_ψ = jacob[iz]*delzed[iz] / sum(grho[iz]*jacob[iz]*delzed[iz]) 
      if (.not. flux_norm) fluxnorm_vs_z = jacob_times_delzed_vs_z / sum(jacob_times_delzed_vs_z) 
      if (flux_norm) fluxnorm_vs_z = jacob_times_delzed_vs_z / sum(jacob_times_delzed_vs_z * grho(1, :))

      ! If <flux_norm> = True, we include the factor <∇̃ρ>_ψ in the flux definition, otherwise we ignore it 
      ! 1/<∇̃ρ>_ψ = int dV / int ∇̃ρ dV = sum(jacob[iz]*delzed[iz]) / sum(grho[iz]*jacob[iz]*delzed[iz])
      if (flux_norm) one_over_nablarho = sum(jacob_times_delzed_vs_z) / sum(jacob_times_delzed_vs_z * grho(1, :))
      if (.not. flux_norm) one_over_nablarho = 1.0

      ! Deallocate arrays
      deallocate (jacob_times_delzed_vs_z)

   end subroutine get_factor_for_fluxsurfaceaverage

end module diagnose_fluxes_fluxtube
