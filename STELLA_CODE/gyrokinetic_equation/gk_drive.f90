!###############################################################################
!################### DRIVE TERM OF THE GYROKINETIC EQUATION ####################
!###############################################################################
! 
! This module evolves the drive term:
!     - i omega_{*,k,s} J_0 ϕ_k
! 
! 
! Mathematics
! -----------
! 
! The normalized diamagnetic frequency omega_{*,k,s} is given by equation (22) in [2019 - Barnes]:
!     omega_{*,k,s} = (1/C) (ky_SI*rho_ref/2) a*Bref (dy_SI/dalpha_SI) exp(-v²) d ln F_s / d psi_SI
!     
! The equations are derived using a generalized Clebsch form
!     B = <clebsch_factor> ∇ψ x ∇α
! 
! Consider the following normalizations
!     rho = r_SI/a
!     psi = psi_SI/(a²*Bref)
! 
! The normalized perpendicular gradients are rescaled with (rho_ref/a) to be of order one
!     <dydalpha> = (rho_ref/a) (d(y_SI/rho_ref)/dalpha) = (1/a)(dy_SI/dalpha)
!     <drhodpsi> = drho/dpsi = d(r/a)/d(psi_SI/(a^2*Br)) = (a*Bref) * dr/dpsi_SI
! 
! The normalized gradient scale lengths are calculated with respect to d/dx = d/dr
!     <fprim> = - a/n_s (dn_s/dx) = - a/n_s (dn_s/dr)
!     <tprim> = - a/T_s (dT_s/dx) = - a/T_s (dT_s/dr)
! 
! The radial derivative of the Maxwellian distribution function F_s with respect to rho = r/a is given by,
!     d F_s / drho = [a/n_s dn_s/dr + a/T_s dT_s/dr ( v_parallel² + 2 mu B - 3/2 )] F_s

! Rewriting the normalized diamagnetic frequency omega_{*,k,s} in function of code variables gives,
!     omega_{*,k,s} = (1/C) (ky_SI*rho_ref/2) a*Bref (dy_SI/dalpha_SI) exp(-v²) d ln F_s / d psi_SI
!                   = (1/C) (ky/2) a*Bref a*<dydalpha> exp(-v²) d ln F_s / d psi * (1/a²Bref)
!                   = ky * 1/(2C) <dydalpha> exp(-v²) (drho/dpsi) d ln F_s / d rho
!                   = ky * 0.5/<clebsch_factor> <dydalpha> exp(-v²) <drhodpsi> [<fprim> + <tprim> (v_parallel² + 2 mu B - 1.5)] 
! 
! Gather the term containing the kinetic energy as
!     <energy>[ialpha,iz] = v_parallel² + 2 mu B = vpa(iv)**2 + vperp2(ialpha, iz, imu)
! 
! The exponential coming from the Maxwellian is calculated as
!     exp(-v²) = exp(-v²_parallel -v²_perp) = exp(-<vpa>*<vpa>) * exp (-2*<mu>*<bmag>)
!              = maxwell_vpa(iv, is) * maxwell_mu(ialpha, iz, imu, is)
! 
! Note that the radial variation terms are equal to one if radial_variation = False, i.e.,
!      spec(is)%temp_psi0 / spec(is)%temp = 1
!      maxwell_fac = spec%dens / spec%dens_psi0 * (spec%temp_psi0 / spec%temp)**1.5 = 1
! 
! 
! Variables defined in this routine
! ---------------------------------
!     wstar = - code_dt*omega_{*,k,s}/ky = - code_dt/(2C) <dydalpha> exp(-v²) (drho/dpsi) d ln F_s / d rho
! 
!###############################################################################
module gk_drive

   use debug_flags, only: debug => time_advance_debug

   implicit none
   
   ! Make these routine available to gk_time_advance()
   public :: init_wstar
   public :: finish_wstar
   public :: advance_wstar_explicit

   private

contains

   !****************************************************************************
   !                   INITIALISE THE DIAMAGNETIC FREQUENCY                    !
   !****************************************************************************
   subroutine init_wstar

      use mp, only: mp_abort
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use stella_time, only: code_dt
      use geometry, only: dydalpha, drhodpsi, clebsch_factor
      use neoclassical_terms, only: include_neoclassical_terms, dfneo_drho
      use parameters_numerical, only: maxwellian_normalization
      use arrays_store_useful, only: wstar, wstarinit
      
      ! Grids
      use parameters_kxky_grid, only: nalpha
      use grids_z, only: nzgrid
      use grids_species, only: spec
      use grids_velocity, only: vperp2, vpa
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
      ! Rescale the drive term with <wstarknob>
      use parameters_physics, only: wstarknob

      implicit none

      ! Indices
      integer :: is, imu, iv, ivmu
      
      ! To make the calculations easier to follow, we
      ! calculate <energy> = v_parallel² + 2 mu B for each velocity point
      real, dimension(:, :), allocatable :: energy
         
      !-------------------------------------------------------------------------

      ! Only intialise omega_{*,k,s} once
      if (wstarinit) return
      wstarinit = .true.

      ! Allocate omega_{*,k,s} = wstar[ialpha, iz, i[mu,vpa,s]]
      if (.not. allocated(wstar)) then
         allocate (wstar(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstar = 0.0
      end if

      ! Allocate <energy>[ialpha,iz]
      allocate (energy(nalpha, -nzgrid:nzgrid))
      
      ! Iterate over velocity space
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         is = is_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         iv = iv_idx(vmu_lo, ivmu)
         
         ! Calculate <energy>[ialpha,iz] = v_parallel² + 2 mu B = vpa(iv)**2 + vperp2(ialpha, iz, imu)
         energy = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
         
         ! Calculate wstar = - code_dt*omega_{*,k,s}/ky = - code_dt * 0.5/C <dydalpha> exp(-v²) (drho/dpsi) d ln F_s / d rho 
         !  = - code_dt * 0.5/<clebsch_factor> <dydalpha> <drhodpsi> [<fprim> + <tprim> (v_parallel² + 2 mu B - 1.5)] exp(-v²)
         if (include_neoclassical_terms) then
            if (maxwellian_normalization) then
               call mp_abort("include_neoclassical_terms = T not yet supported for maxwellian_normalization = T. Aborting.")
            else
            wstar(:, :, ivmu) = - (1/clebsch_factor) * dydalpha * drhodpsi * wstarknob * 0.5 * code_dt &
                * (maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is) &
                * (spec(is)%fprim + spec(is)%tprim * (energy - 1.5)) - dfneo_drho(:, :, ivmu))
            end if
         else
             wstar(:, :, ivmu) = - (1/clebsch_factor) * dydalpha * drhodpsi * wstarknob * 0.5 * code_dt &
                  * (spec(is)%fprim + spec(is)%tprim * (energy - 1.5))
         end if
         
         ! Add exp(-v²) = exp(-v²_parallel -v²_perp) = maxwell_vpa(iv, is) * maxwell_mu(ialpha, iz, imu, is) to wstar
         if (.not. maxwellian_normalization) then
               wstar(:, :, ivmu) = wstar(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)
         end if
         
      end do

      deallocate (energy)

   end subroutine init_wstar

   !*****************************************************************************
   !                           ADVANCE DRIVE TERM
   !*****************************************************************************
   subroutine advance_wstar_explicit(phi, gout)
   
      ! Grids
      use grids_z, only: nzgrid
      use stella_layouts, only: vmu_lo
      
      ! Physics
      use parameters_physics, only: full_flux_surface
   
      implicit none
   
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
   
      if (full_flux_surface) then
         call advance_wstar_explicit_ffs(gout)
      else 
         call advance_wstar_explicit_flux_tube(phi, gout)
      end if
   
   end subroutine advance_wstar_explicit
   
   !-------------------------------- Flux tube ---------------------------------
   subroutine advance_wstar_explicit_flux_tube(phi, gout)

      ! Parallelisation
      use mp, only: proc0
      
      ! Data arrays
      use arrays_store_useful, only: wstar
      use arrays_store_fields, only: apar, bpar
      
      ! Grids
      use stella_layouts, only: vmu_lo
      use grids_z, only: nzgrid, ntubes
      use parameters_kxky_grid, only: naky, nakx
      
      ! Calculations
      use add_explicit_terms, only: add_explicit_term
      use calculations_kxky_derivatives, only: get_dchidy
      
      ! Time this routine
      use arrays_store_useful, only: time_gke
      use job_manage, only: time_message

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
      complex, dimension(:, :, :, :, :), allocatable :: g0
         
      !-------------------------------------------------------------------------
      ! Add the drive term to the gyrokinetic equation:
      !      - i omega_{*,k,s} J_0 ϕ_k = <wstar> * i ky J_0 ϕ_k
      ! 
      ! First we calculate i ky J_0 ϕ_k which corresponds to d<chi>_theta/dy
      ! with chi = ϕ − v · δA/c = fphi ∗ phi − fapar ∗ vpa(iv) ∗ sqrt(T/m) ∗ apar 
      !     <g0> = i ky J_0 ϕ_k = get_dchidy(phi, apar, bpar, g0)
      ! 
      ! Then multiply with <wstar> and add it to the right-hand-side of the gyrokinetic equation
      !     add_explicit_term(g0, wstar(1, :, :), gout)
      !-------------------------------------------------------------------------

      ! Start timing the time advance due to the driving gradients
      if (proc0) call time_message(.false., time_gke(:, 6), ' wstar advance')

      ! Allocate temporary array for <g0> = i ky J_0 ϕ_k
      allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      
      ! Calculate <g0> = i ky J_0 ϕ_k = d<chi>_theta/dy
      if (debug) write (*, *) 'time_advance::solve_gke::get_dchidy'
      call get_dchidy(phi, apar, bpar, g0)
      
      ! Add the drive term to the right-hand-side of the gyrokinetic equation
      if (debug) write (*, *) 'time_advance::solve_gke::add_wstar_term'
      call add_explicit_term(g0, wstar(1, :, :), gout)

      ! Deallocate <g0> = i ky J_0 ϕ_k
      deallocate (g0)

      ! Stop timing the time advance due to the driving gradients
      if (proc0) call time_message(.false., time_gke(:, 6), ' wstar advance')

   end subroutine advance_wstar_explicit_flux_tube
   
   !---------------------------- Full flux surface -----------------------------
   subroutine advance_wstar_explicit_ffs(gout)

      use mp, only: proc0, mp_abort
      use job_manage, only: time_message
      use stella_layouts, only: vmu_lo
      use calculations_transforms, only: transform_ky2y
      use grids_z, only: nzgrid, ntubes
      use parameters_kxky_grid, only: naky, naky_all, nakx, ikx_max, ny
      use calculations_kxky, only: swap_kxky
      use arrays_store_useful, only: wstar
      use arrays_store_distribution_fn, only: g_scratch
      use calculations_gyro_averages, only: gyro_average
      use add_explicit_terms, only: add_explicit_term_ffs
      use calculations_kxky_derivatives, only: get_dgdy, get_dchidy
      use arrays_store_useful, only: time_gke

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
      complex, dimension(:, :, :, :, :), allocatable :: g0, g0y
      complex, dimension(:, :), allocatable :: g0_swap
      integer :: iz, it, ivmu
         
      !-------------------------------------------------------------------------

      ! Start timing the time advance due to the driving gradients
      if (proc0) call time_message(.false., time_gke(:, 6), ' wstar advance')
      
      ! Assume only a single flux surface simulated
      it = 1

      ! Allocate temporary arrays
      allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      allocate (g0_swap(naky_all, ikx_max))

      ! Calculate d<phi>/dy in k-space, i.e., calculate i*ky*J_0*chi, and save it as <g0>
      ! Here g_scratch is <phi> in k-space that has been pre-calculated and stored
      call get_dgdy(g_scratch, g0)
      
      ! Transform d<phi>/dy from ky-space to y-space
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do iz = -nzgrid, nzgrid
            call swap_kxky(g0(:, :, iz, it, ivmu), g0_swap)
            call transform_ky2y(g0_swap, g0y(:, :, iz, it, ivmu))
         end do
      end do
      
      ! multiply d<chi>/dy with omega_* coefficient and add to source (RHS of GK eqn)
      call add_explicit_term_ffs(g0y, wstar, gout)
      
      ! Deallocate
      deallocate (g0y, g0_swap)
      deallocate (g0)

      ! Stop timing the time advance due to the driving gradients
      if (proc0) call time_message(.false., time_gke(:, 6), ' wstar advance')

   end subroutine advance_wstar_explicit_ffs

   !*****************************************************************************
   !                           FINALISE DRIVE TERM
   !*****************************************************************************
   subroutine finish_wstar

      use arrays_store_useful, only: wstar, wstarinit

      implicit none

      if (allocated(wstar)) deallocate (wstar)
      wstarinit = .false.

   end subroutine finish_wstar
   
end module gk_drive
