!###############################################################################
!############### MAGNETIC DRIFT TERM OF THE GYROKINETIC EQUATION ###############
!###############################################################################
! 
! This module evolves the magnetic drift term:
!     - i omega_{d,k,s} (g_{k,s} + Z_s/T_s J_0 ϕ_k exp(-v²)
! 
! 
! Mathematics
! -----------
! 
! The normalized magnetic drift frequency omega_{d,k,s} is given by equation (24) in [2019 - Barnes]:
!     omega_{d,k,s} = Ts/Zs 1/B (v_parallel² v_kappa + mu_s v_∇B) (ky ∇y + kx ∇x)
! 
! Split up the magnetic drift frequencies in the contributions from curvature and ∇B
!     omega_{curvature,d,k,s} = Ts/Zs 1/B (v_parallel² v_kappa) (ky ∇y + kx ∇x)
!                             = v_parallel² 0.5 Ts/Zs (kx/shat*<cvdrift0> + ky*<cvdrift>)
!                             = -1/code_dt (kx*<wcvdriftx> + ky*<wcvdrifty>)
!            omega_{∇B,d,k,s} = Ts/Zs 1/B (mu_s v_∇B) (ky ∇y + kx ∇x)
!                             = v_perp² 0.25 Ts/Zs (kx/shat*<gbvdrift0> + ky*<gbdrift>)
!                             = -1/code_dt (kx*<wgbdriftx> + ky*<wgbdrifty>)
! 
! Variables defined in this routine
! ---------------------------------
! Curvature (cv) components of the magnetic drift along (kx,ky):
!    kx*<wcvdriftx> + ky*<wcvdrifty> = - code_dt * omega_{curvature,d,k,s}
!    <wcvdriftx> = -0.5 * code_dt * v_parallel² * Ts/Zs * <cvdrift0>/shat
!    <wcvdrifty> = -0.5 * code_dt * v_parallel² * Ts/Zs * <cvdrift>
! 
! Gradient B (gb) components of the magnetic drift along (kx,ky):
!    kx*<wgbdriftx> + ky*<wgbdrifty> = - code_dt * omega_{∇B,d,k,s}
!    <wgbdriftx> = -0.25 * code_dt * v_perp² * Ts/Zs * <gbvdrift0>/shat
!    <wgbdrifty> = -0.25 * code_dt * v_perp² * Ts/Zs * <gbvdrift>
! 
!###############################################################################
module gk_magnetic_drift

   ! Load debug flags
   use debug_flags, only: debug => time_advance_debug

   implicit none
   
   ! Make these routine available to gk_time_advance()
   public :: init_wdrift
   public :: finish_wdrift
   public :: advance_wdriftx_explicit
   public :: advance_wdrifty_explicit

   private

contains

   !****************************************************************************
   !                           Initialise explicit drifts
   !****************************************************************************
   subroutine init_wdrift

      use neoclassical_terms, only: include_neoclassical_terms
      use arrays, only: initialised_wdrift
      
      implicit none

      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_wdrift) return
      initialised_wdrift = .true.
      
      ! Allocate arrays that will be used thoughout the time advance
      call allocate_arrays_wdrift

      ! Add the neoclassical terms in a seperate subroutine since they 
      ! require a lot of care with regards to vpa-space
      if (include_neoclassical_terms) then
         call init_wdrift_with_neoclassical_terms
      else
         call init_wdrift_without_neoclassical_terms
      end if

   end subroutine init_wdrift
   
   !------------------------- Without neoclassical terms -----------------------
   subroutine init_wdrift_without_neoclassical_terms
      
      ! Parallelisation
      use mp, only: mp_abort
      
      ! Numerics
      use parameters_numerical, only: maxwellian_normalization
      
      ! Geometry
      use geometry, only: cvdrift, gbdrift
      use geometry, only: cvdrift0, gbdrift0
      use geometry, only: geo_surf, q_as_x
      
      ! Grids
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use grids_time, only: code_dt
      use grids_species, only: spec
      use grids_z, only: nzgrid
      use grids_kxky, only: nalpha
      use grids_velocity, only: vpa, vperp2, mu
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
      ! Scale the drift term in the gyrokinetic equation
      use parameters_physics, only: xdriftknob, ydriftknob
      
      ! This routine fills the following arrays, with dimensions [ialpha, iz, ivmu]
      use arrays, only: wdriftx_g, wdrifty_g
      use arrays, only: wdriftx_phi, wdrifty_phi
      use arrays, only: wdriftx_bpar, wdrifty_bpar
      
      implicit none
      
      ! Indices
      integer :: ivmu, iv, imu, is

      ! Gather calculations
      real :: fac
      
      ! Temporary arrays
      real, dimension(:, :), allocatable :: wcvdrifty, wgbdrifty
      real, dimension(:, :), allocatable :: wcvdriftx, wgbdriftx
      
      !-------------------------------------------------------------------------
      ! In this routine we calculate the following arrays with dimensions [ialpha, iz, ivmu]
      !    <wcvdriftx> = 0.5 * code_dt * v_parallel² * Ts/Zs * <cvdrift0>/shat
      !    <wcvdrifty> = 0.5 * code_dt * v_parallel² * Ts/Zs * <cvdrift>
      !    <wgbdriftx> = 0.25 * code_dt * v_perp² * Ts/Zs * <gbvdrift0>/shat
      !    <wgbdrifty> = 0.25 * code_dt * v_perp² * Ts/Zs * <gbvdrift>
      ! 
      ! The magnetic drift term in the gyrokinetic equation is given by,
      !  - i omega_{d,k,s} (g_{k,s} + Z_s/T_s J_0 ϕ_k exp(-v²) )
      ! 
      ! Therefore, we add the fixed terms for g_{k,s} and ϕ_k
      !    <wdriftx_g>[ialpha, iz, ivmu] = <wcvdriftx> + <wgbdriftx>
      !    <wdrifty_g>[ialpha, iz, ivmu] = <wcvdrifty> + <wgbdrifty>
      !    <wdriftx_phi>[ialpha, iz, ivmu] = Z_s/T_s * exp(-v²) * (<wcvdriftx> + <wgbdriftx>)
      !    <wdrifty_phi>[ialpha, iz, ivmu] = Z_s/T_s * exp(-v²) * (<wcvdrifty> + <wgbdrifty>)
      !-------------------------------------------------------------------------

      ! Allocate temporary arrays
      allocate (wcvdrifty(nalpha, -nzgrid:nzgrid))
      allocate (wgbdrifty(nalpha, -nzgrid:nzgrid))
      allocate (wcvdriftx(nalpha, -nzgrid:nzgrid))
      allocate (wgbdriftx(nalpha, -nzgrid:nzgrid))

      ! Iterate over velocity space
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         !----------------------------------------------------------------------
         !----------------------- Calculate y-components -----------------------
         !----------------------------------------------------------------------
         
         ! Calculate <wcvdrifty> = 0.5 * code_dt * Ts/Zs * v_parallel² * <cvdrift>
         ! and <wgbdrifty> = 0.25 * code_dt * Ts/Zs * v_perp² * <gbvdrift>
         ! We also add the input parameter <ydriftknob> to rescale the y-drifts
         fac = -ydriftknob * 0.5 * code_dt * spec(is)%tz_psi0
         wcvdrifty = fac * cvdrift * vpa(iv) * vpa(iv)
         wgbdrifty = fac * gbdrift * 0.5 * vperp2(:, :, imu)
         
         ! Calculate <wdrifty_g>[ialpha, iz, ivmu] = <wcvdrifty> + <wgbdrifty>
         wdrifty_g(:, :, ivmu) = wcvdrifty + wgbdrifty
         
         ! Calculate <wdrifty_phi>[ialpha, iz, ivmu] = Z_s/T_s * exp(-v²) * (<wcvdrifty> + <wgbdrifty>)
         wdrifty_phi(:, :, ivmu) = spec(is)%zt * (wcvdrifty + wgbdrifty)
         if (.not. maxwellian_normalization) then
            wdrifty_phi(:, :, ivmu) = wdrifty_phi(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
         end if
         
         ! TODO - write documentation and neoclassical terms not supported
         wdrifty_bpar(:,:,ivmu) = 4.0 * mu(imu) * wdrifty_phi(:, :, ivmu) * spec(is)%tz
         
         !----------------------------------------------------------------------
         !----------------------- Calculate x-components -----------------------
         !----------------------------------------------------------------------

         ! Calculate <wcvdriftx> = 0.5 * code_dt * v_parallel² * Ts/Zs * <cvdrift0>/shat
         ! and <wgbdriftx> = 0.25 * code_dt * v_perp² * Ts/Zs * <gbvdrift0>/shat
         ! Note that if x = q instead of x = r, we do not need the shat term here
         if (q_as_x) then
            fac = -xdriftknob * 0.5 * code_dt * spec(is)%tz_psi0
         else
            fac = -xdriftknob * 0.5 * code_dt * spec(is)%tz_psi0 / geo_surf%shat
         end if
         wcvdriftx = fac * cvdrift0 * vpa(iv) * vpa(iv)
         wgbdriftx = fac * gbdrift0 * 0.5 * vperp2(:, :, imu)
         
         ! Calculate <wdriftx_g>[ialpha, iz, ivmu] = <wcvdriftx> + <wgbdriftx>
         wdriftx_g(:, :, ivmu) = wcvdriftx + wgbdriftx
         
         ! Calculate <wdriftx_phi>[ialpha, iz, ivmu] = Z_s/T_s * exp(-v²) * (<wcvdriftx> + <wgbdriftx>)
         wdriftx_phi(:, :, ivmu) = spec(is)%zt * (wcvdriftx + wgbdriftx)
         if (.not. maxwellian_normalization) then
            wdriftx_phi(:, :, ivmu) = wdriftx_phi(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
         end if
         
         ! TODO - write documentation and neoclassical terms not supported
         wdriftx_bpar(:,:,ivmu) = 4.0 * mu(imu) * wdriftx_phi(:, :, ivmu) * spec(is)%tz

      end do

      deallocate (wcvdriftx, wgbdriftx, wcvdrifty, wgbdrifty) 

   end subroutine init_wdrift_without_neoclassical_terms
   
   !-------------------------- With neoclassical terms -------------------------
   subroutine init_wdrift_with_neoclassical_terms

      use mp, only: mp_abort
      use geometry, only: cvdrift, gbdrift
      use geometry, only: cvdrift0, gbdrift0
      use geometry, only: gds23, gds24
      use geometry, only: geo_surf, q_as_x
      use geometry, only: dxdpsi, drhodpsi, dydalpha
      use neoclassical_terms, only: include_neoclassical_terms
      use neoclassical_terms, only: dphineo_dzed, dphineo_drho, dphineo_dalpha
      use neoclassical_terms, only: dfneo_dvpa, dfneo_dzed, dfneo_dalpha
      use parameters_numerical, only: maxwellian_normalization
      use parameters_physics, only: xdriftknob, ydriftknob
      
      ! Grids
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use grids_time, only: code_dt
      use grids_species, only: spec
      use grids_z, only: nzgrid
      use grids_kxky, only: nalpha
      use grids_velocity, only: vpa, vperp2, mu
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
      ! This routine fills the following arrays, with dimensions [ialpha, iz, ivmu]
      use arrays, only: wdriftx_g, wdrifty_g
      use arrays, only: wdriftx_phi, wdrifty_phi
      use arrays, only: wdriftx_bpar, wdrifty_bpar
      
      implicit none
      
      ! Indices
      integer :: ivmu, iv, imu, is

      ! Gather calculations
      real :: fac
      
      ! Temporary arrays
      real, dimension(:, :), allocatable :: wcvdrifty, wgbdrifty
      real, dimension(:, :), allocatable :: wcvdriftx, wgbdriftx
      
      !-------------------------------------------------------------------------
      ! FLAG -- need to deal with shat=0 case.  ideally move away from q as x-coordinate
      !-------------------------------------------------------------------------

      ! Allocate temporary arrays
      allocate (wcvdrifty(nalpha, -nzgrid:nzgrid))
      allocate (wgbdrifty(nalpha, -nzgrid:nzgrid))
      allocate (wcvdriftx(nalpha, -nzgrid:nzgrid))
      allocate (wgbdriftx(nalpha, -nzgrid:nzgrid))

      ! Iterate over velocity space
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         
         fac = -ydriftknob * 0.5 * code_dt * spec(is)%tz_psi0
         
         ! Calculate <wcvdrifty> = 0.5 * code_dt * Ts/Zs * v_parallel² * <cvdrift>
         ! This is the curvature drift piece of wdrifty with missing factor of vpa
         ! vpa factor is missing to avoid singularity when including
         ! non-Maxwellian corrections to equilibrium
         wcvdrifty = fac * cvdrift * vpa(iv)
         
         ! This is the grad-B drift piece of wdrifty
         wgbdrifty = fac * gbdrift * 0.5 * vperp2(:, :, imu)
         wdrifty_g(:, :, ivmu) = wcvdrifty * vpa(iv) + wgbdrifty
         
         ! if including neoclassical correction to equilibrium Maxwellian,
         ! then add in v_E^{nc} . grad y dg/dy coefficient here
         if (include_neoclassical_terms) then
               wdrifty_g(:, :, ivmu) = wdrifty_g(:, :, ivmu) + code_dt * 0.5 * (gds23 * dphineo_dzed + drhodpsi * dydalpha * dphineo_drho)
         end if

         wdrifty_phi(:, :, ivmu) = spec(is)%zt * (wgbdrifty + wcvdrifty * vpa(iv))

         ! if maxwwellian_normalization = .true., evolved distribution function is normalised by a Maxwellian
         ! otherwise, it is not; a Maxwellian weighting factor must thus be included
         if (.not. maxwellian_normalization) then
               wdrifty_phi(:, :, ivmu) = wdrifty_phi(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
         end if
         
         ! assign wdrifty_bpar, neoclassical terms not supported
         wdrifty_bpar(:,:,ivmu) = 4.0 * mu(imu) * wdrifty_phi(:, :, ivmu) * spec(is)%tz
         
         ! if including neoclassical corrections to equilibrium,
         ! add in -(Ze/m) * v_curv/vpa . grad y d<phi>/dy * dF^{nc}/dvpa term
         ! and v_E . grad z dF^{nc}/dz (here get the dphi/dy part of v_E)
         if (include_neoclassical_terms) then
               ! NB: the below neoclassical correction needs to be divided by an equilibrium Maxwellian
               ! if maxwellian_normalization = .true.
               if (maxwellian_normalization) then
               call mp_abort("include_neoclassical_terms=T not currently supported for maxwellian_normalization=T.  aborting")
               end if
               wdrifty_phi(:, :, ivmu) = wdrifty_phi(:, :, ivmu) &
                  - 0.5 * spec(is)%zt * dfneo_dvpa(:, :, ivmu) * wcvdrifty &
                  - code_dt * 0.5 * dfneo_dzed(:, :, ivmu) * gds23
         end if

         if (q_as_x) then
               fac = -xdriftknob * 0.5 * code_dt * spec(is)%tz_psi0
         else
               fac = -xdriftknob * 0.5 * code_dt * spec(is)%tz_psi0 / geo_surf%shat
         end if
         
         ! This is the curvature drift piece of wdriftx with missing factor of vpa
         ! vpa factor is missing to avoid singularity when including
         ! non-Maxwellian corrections to equilibrium
         wcvdriftx = fac * cvdrift0 * vpa(iv)
         
         ! This is the grad-B drift piece of wdriftx
         wgbdriftx = fac * gbdrift0 * 0.5 * vperp2(:, :, imu)
         wdriftx_g(:, :, ivmu) = wcvdriftx * vpa(iv) + wgbdriftx
         
         ! if including neoclassical correction to equilibrium Maxwellian,
         ! then add in v_E^{nc} . grad x dg/dx coefficient here
         if (include_neoclassical_terms) then
               wdriftx_g(:, :, ivmu) = wdriftx_g(:, :, ivmu) + code_dt * 0.5 * (gds24 * dphineo_dzed - dxdpsi * dphineo_dalpha)
         end if
         
         wdriftx_phi(:, :, ivmu) = spec(is)%zt * (wgbdriftx + wcvdriftx * vpa(iv))
         
         ! if maxwellian_normalizatiion = .true., evolved distribution function is normalised by a Maxwellian
         ! otherwise, it is not; a Maxwellian weighting factor must thus be included
         if (.not. maxwellian_normalization) then
               wdriftx_phi(:, :, ivmu) = wdriftx_phi(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
         end if
         
         ! assign wdriftx_bpar, neoclassical terms not supported
         wdriftx_bpar(:,:,ivmu) = 4.0 * mu(imu) * wdriftx_phi(:, :, ivmu) * spec(is)%tz
         
         ! if including neoclassical corrections to equilibrium,
         ! add in (Ze/m) * v_curv/vpa . grad x d<phi>/dx * dF^{nc}/dvpa term
         ! and v_E . grad z dF^{nc}/dz (here get the dphi/dx part of v_E)
         ! and v_E . grad alpha dF^{nc}/dalpha (dphi/dx part of v_E)
         if (include_neoclassical_terms) then
               ! NB: the below neoclassical correction needs to be divided by an equilibrium Maxwellian
               ! if running with maxwellian_normalzation = .true.
               if (maxwellian_normalization) then
               call mp_abort("include_neoclassical_terms=T not currently supported for maxwellian_normalization=T.  aborting")
               end if
               wdriftx_phi(:, :, ivmu) = wdriftx_phi(:, :, ivmu) &
                                       - 0.5 * spec(is)%zt * dfneo_dvpa(:, :, ivmu) * wcvdriftx &
                                       + code_dt * 0.5 * (dfneo_dalpha(:, :, ivmu) * dxdpsi - dfneo_dzed(:, :, ivmu) * gds24)
         end if

      end do

      deallocate (wcvdriftx, wgbdriftx, wcvdrifty, wgbdrifty) 

   end subroutine init_wdrift_with_neoclassical_terms
      
   !---------------------------- Allocate arrays ---------------------------
   subroutine allocate_arrays_wdrift
      
      ! Grids
      use grids_z, only: nzgrid
      use grids_kxky, only: nalpha
      use stella_layouts, only: vmu_lo
    
      ! Allocate the following arrays with dimensions [ialpha, iz, ivmu]
      use arrays, only: wdriftx_g, wdrifty_g
      use arrays, only: wdriftx_phi, wdrifty_phi
      use arrays, only: wdriftx_bpar, wdrifty_bpar

      implicit none

      ! Allocate wdriftx_phi, the factor multiplying dphi/dx in the magnetic drift term
      if (.not. allocated(wdriftx_phi)) then
         allocate (wdriftx_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         wdriftx_phi = 0.0
      end if
      
      ! Allocate wdrifty_phi, the factor multiplying dphi/dy in the magnetic drift term
      if (.not. allocated(wdrifty_phi)) then
         allocate (wdrifty_phi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         wdrifty_phi = 0.0
      end if
      
      ! Allocate wdriftx_bpar, the factor multiplying dbpar/dx in the magnetic drift term
      if (.not. allocated(wdriftx_bpar)) then
         allocate (wdriftx_bpar(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         wdriftx_bpar = 0.0
      end if
      
      ! Allocate wdrifty_bpar, the factor multiplying dbpar/dy in the magnetic drift term
      if (.not. allocated(wdrifty_bpar)) then
         allocate (wdrifty_bpar(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         wdrifty_bpar = 0.0
      end if
      
      ! Allocate wdriftx_g, the factor multiplying dg/dx in the magnetic drift term
      if (.not. allocated(wdriftx_g)) then
         allocate (wdriftx_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         wdriftx_g = 0.0
      end if
      
      ! Allocate wdrifty_g, the factor multiplying dg/dy in the magnetic drift term
      if (.not. allocated(wdrifty_g)) then
         allocate (wdrifty_g(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         wdrifty_g = 0.0
      end if
      
   end subroutine allocate_arrays_wdrift

   !*****************************************************************************
   !                           Advance explicit drifts
   !*****************************************************************************
   ! advance_wdrifty_explicit subroutine calculates and adds the y-component of the
   ! magnetic drift term to the RHS of the GK equation
   !*****************************************************************************
   subroutine advance_wdrifty_explicit(g, phi, bpar, gout)

      use mp, only: proc0
      use stella_layouts, only: vmu_lo
      use job_manage, only: time_message
      use calculations_transforms, only: transform_ky2y
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: nakx, ikx_max, naky, naky_all, ny
      use calculations_kxky, only: swap_kxky
      use parameters_physics, only: full_flux_surface, include_bpar
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      use arrays, only: wdrifty_g, wdrifty_phi, wdrifty_bpar
      use arrays_distribution_function, only: g_scratch
      use calculations_add_explicit_terms, only: add_explicit_term, add_explicit_term_ffs
      use calculations_kxky_derivatives, only: get_dgdy
      use arrays, only: time_gke

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      integer :: ivmu, iz, it
      complex, dimension(:, :, :, :), allocatable :: dphidy, dbpardy
      complex, dimension(:, :, :, :, :), allocatable :: g0k, g0y
      complex, dimension(:, :), allocatable :: g0k_swap

      !-------------------------------------------------------------------------

      ! start the timing of the y component of the magnetic drift advance
      if (proc0) call time_message(.false., time_gke(:, 4), ' dgdy advance')

      allocate (dphidy(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dbpardy(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g0k(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gke::advance_wdrifty_explicit::get_dgdy'
      ! calculate dg/dy in (ky,kx) space
      call get_dgdy(g, g0k)
      ! calculate dbpar/dy in (ky,kx) space
      if (include_bpar) call get_dgdy(bpar, dbpardy)

      if (full_flux_surface) then
         ! assume only a single flux surface simulated
         it = 1
         allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (g0k_swap(naky_all, ikx_max))
         ! transform dg/dy from k-space to y-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               do iz = -nzgrid, nzgrid
               call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
               call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
               end do
         end do

         ! add vM . grad y dg/dy term to equation
         call add_explicit_term_ffs(g0y, wdrifty_g, gout)

         ! > calculate dphi/dy in (ky,kx) space
         ! Here g_scratch is <phi> in k-space that has been pre-calculated and stored
         call get_dgdy(g_scratch, g0k)

         ! transform d<phi>/dy from k-space to y-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               do iz = -nzgrid, nzgrid
               call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
               call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
               end do
         end do

         ! add vM . grad y d<phi>/dy term to equation
         call add_explicit_term_ffs(g0y, wdrifty_phi, gout)

         deallocate (g0y, g0k_swap)
      else
         if (debug) write (*, *) 'time_advance::solve_gke::add_dgdy_term'
         ! add vM . grad y dg/dy term to equation
         call add_explicit_term(g0k, wdrifty_g(1, :, :), gout)

         ! Note that this is here because for FFS te gyro-average is calculated once outside this routine
         ! TODO-GA: can we do something similar for fluxtube to save cpu time?
         ! calculate dphi/dy in (ky,kx) space
         call get_dgdy(phi, dphidy)

         ! get <dphi/dy> in k-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               call gyro_average(dphidy, ivmu, g0k(:, :, :, :, ivmu))
         end do

         ! add vM . grad y d<phi>/dy term to equation
         call add_explicit_term(g0k, wdrifty_phi(1, :, :), gout)
         
         if (include_bpar) then
               ! get <dbpar/dy> in k-space
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               call gyro_average_j1(dbpardy, ivmu, g0k(:, :, :, :, ivmu))
               end do
               ! add vM . grad y (4 mu d<bpar>/dy) term to equation
               call add_explicit_term(g0k, wdrifty_bpar(1, :, :), gout)            
         end if
      end if
      deallocate (g0k, dphidy, dbpardy)

      ! stop the timing of the y component of the magnetic drift advance
      if (proc0) call time_message(.false., time_gke(:, 4), ' dgdy advance')

   end subroutine advance_wdrifty_explicit

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   ! advance_wdriftx_explicit subroutine calculates and adds the x-component of the
   ! magnetic drift term to the RHS of the GK equation
   !****************************************************************************
   subroutine advance_wdriftx_explicit(g, phi, bpar, gout)

      use mp, only: proc0
      use stella_layouts, only: vmu_lo
      use job_manage, only: time_message
      use calculations_transforms, only: transform_ky2y
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: nakx, ikx_max, naky, naky_all, ny
      use grids_kxky, only: akx
      use calculations_kxky, only: swap_kxky
      use parameters_physics, only: full_flux_surface, include_bpar
      use calculations_gyro_averages, only: gyro_average
      use arrays, only: wdriftx_g, wdriftx_phi, wdriftx_bpar
      use arrays_distribution_function, only: g_scratch
      use calculations_kxky_derivatives, only: get_dgdx
      use calculations_add_explicit_terms, only: add_explicit_term, add_explicit_term_ffs
      use arrays, only: time_gke

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      integer :: ivmu, iz, it
      complex, dimension(:, :, :, :), allocatable :: dphidx, dbpardx
      complex, dimension(:, :, :, :, :), allocatable :: g0k, g0y
      complex, dimension(:, :), allocatable :: g0k_swap

      !-------------------------------------------------------------------------

      ! start the timing of the x component of the magnetic drift advance
      if (proc0) call time_message(.false., time_gke(:, 5), ' dgdx advance')

      ! do not calculate if wdriftx terms are all zero
      if (maxval(abs(akx)) < epsilon(0.)) then
         if (proc0) call time_message(.false., time_gke(:, 5), ' dgdx advance')
         return
      end if

      allocate (dphidx(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dbpardx(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g0k(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      if (debug) write (*, *) 'time_advance::solve_gke::get_dgdx'
      ! calculate dg/dx in (ky,kx) space
      call get_dgdx(g, g0k)

      ! calculate dbpar/dx in (ky,kx) space
      if (include_bpar) call get_dgdx(bpar, dbpardx)

      if (full_flux_surface) then
         ! assume a single flux surface is simulated
         it = 1
         allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (g0k_swap(naky_all, ikx_max))
         ! transform dg/dx from k-space to y-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               do iz = -nzgrid, nzgrid
               call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
               call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
               end do
         end do
         ! add vM . grad x dg/dx term to equation
         call add_explicit_term_ffs(g0y, wdriftx_g, gout)

         ! Here g_scratch is <phi> in k-space that has been pre-calculated and stored
         ! get <dphi/dx> in k-space
         call get_dgdx(g_scratch, g0k)

         ! transform d<phi>/dx from k-space to y-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               do iz = -nzgrid, nzgrid
               call swap_kxky(g0k(:, :, iz, it, ivmu), g0k_swap)
               call transform_ky2y(g0k_swap, g0y(:, :, iz, it, ivmu))
               end do
         end do
         ! add vM . grad x d<phi>/dx term to equation
         call add_explicit_term_ffs(g0y, wdriftx_phi, gout)
         deallocate (g0y, g0k_swap)
      else
         if (debug) write (*, *) 'time_advance::solve_gke::add_dgdx_term'
         ! add vM . grad x dg/dx term to equation
         call add_explicit_term(g0k, wdriftx_g(1, :, :), gout)
         ! calculate dphi/dx in (ky,kx) space
         call get_dgdx(phi, dphidx)
         ! get <dphi/dx> in k-space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               call gyro_average(dphidx, ivmu, g0k(:, :, :, :, ivmu))
         end do
         ! add vM . grad x d<phi>/dx term to equation
         call add_explicit_term(g0k, wdriftx_phi(1, :, :), gout)
         if (include_bpar) then
               ! get <dbpar/dx> in k-space
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               call gyro_average(dbpardx, ivmu, g0k(:, :, :, :, ivmu))
               end do
               ! add vM . grad x ( 4 mu d<bpar>/dx ) term to equation
               call add_explicit_term(g0k, wdriftx_bpar(1, :, :), gout)
         end if
      end if
      deallocate (g0k, dphidx, dbpardx)

      ! stop the timing of the x component of the magnetic drift advance
      if (proc0) call time_message(.false., time_gke(:, 5), ' dgdx advance')
      
   end subroutine advance_wdriftx_explicit

   !*****************************************************************************
   !                           Finalise explicit drifts
   !*****************************************************************************
   subroutine finish_wdrift

      use arrays, only: wdriftx_g, wdrifty_g
      use arrays, only: wdriftx_phi, wdrifty_phi
      use arrays, only: initialised_wdrift

      implicit none

      if (allocated(wdriftx_g)) deallocate (wdriftx_g)
      if (allocated(wdrifty_g)) deallocate (wdrifty_g)
      if (allocated(wdriftx_phi)) deallocate (wdriftx_phi)
      if (allocated(wdrifty_phi)) deallocate (wdrifty_phi)

      initialised_wdrift = .false.

   end subroutine finish_wdrift
   
end module gk_magnetic_drift
