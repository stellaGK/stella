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

! ================================================================================================================================================================================= !
! -------------------------------------------------------------------- Add documentation for NEO's corrections! ------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

module gk_drive

   ! Load debug flags
   use debug_flags, only: debug => time_advance_debug
   use neoclassical_terms_neo, only: neoclassical_is_enabled

   implicit none
   
   ! Make these routine available to gk_time_advance()
   public :: init_wstar, init_wpol
   public :: finish_wstar, finish_wpol
   public :: advance_wstar_explicit, advance_wpol_explicit

   integer :: neo_option_switch
   integer, parameter :: neo_option_sfincs = 1
   integer, parameter :: neo_option_NEO = 2

   private

contains

   !****************************************************************************
   !                   INITIALISE THE DIAMAGNETIC FREQUENCY                    !
   !****************************************************************************
   subroutine init_wstar

      ! Parallelisation
      use mp, only: mp_abort
      use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      
      ! Grids
      use grids_time, only: code_dt
      use grids_kxky, only: nalpha
      use grids_z, only: nzgrid
      use grids_species, only: spec
      use grids_velocity, only: vperp2, vpa
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
      use geometry, only: dydalpha, drhodpsi, clebsch_factor

      use neoclassical_terms, only: include_neoclassical_terms, dfneo_drho

      use neoclassical_terms_neo, only:  neo_h, neo_phi             
      use neoclassical_terms_neo, only: dneo_h_dpsi, dneo_phi_dpsi   

      use arrays, only: wstar, initialised_wstar

      ! Rescale the drive term with <wstarknob>
      use parameters_physics, only: wstarknob

      use neoclassical_diagnostics, only: write_wpol_diagnostic

      implicit none

      ! Indices
      integer :: is, imu, iv, ivmu, iz
      
      ! To make the calculations easier to follow, we
      ! calculate <energy> = v_parallel² + 2 mu B for each velocity point
      real, dimension(:, :), allocatable :: energy
         
      !-------------------------------------------------------------------------

      ! Only intialise omega_{*,k,s} once
      if (initialised_wstar) return
      initialised_wstar = .true.

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
         !  = - code_dt * 0.5/<clebsch_factor> <dydalpha> <drhodpsi> [<fprim> + <tprim> (v_parallel² + 2 mu B - 1.5)] exp(-v²). 
         ! This block only computes when sfincs is chosen for the neoclassical option.
         if (include_neoclassical_terms .and. neo_option_switch == neo_option_sfincs) then
             wstar(:, :, ivmu) = - (1/clebsch_factor) * dydalpha * drhodpsi * wstarknob * 0.5 * code_dt &
                 * (maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is) &
                 * (spec(is)%fprim + spec(is)%tprim * (energy - 1.5)) - dfneo_drho(:, :, ivmu))
         else
             wstar(:, :, ivmu) = - (1/clebsch_factor) * dydalpha * drhodpsi * wstarknob * 0.5 * code_dt &
             * (spec(is)%fprim + spec(is)%tprim * (energy - 1.5))            
         end if
        
         if (neoclassical_is_enabled()) then
	     do iz = -nzgrid, nzgrid
                 wstar(:, iz, ivmu) = wstar(:, iz, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) &
                 * (1 + neo_h(iz, ivmu, 1) - spec(is)%z * neo_phi(iz, 1)) - (1/clebsch_factor) * dydalpha * drhodpsi * wstarknob * 0.5 * code_dt &
                 * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) * (dneo_h_dpsi(iz, ivmu, 1) - spec(is)%z * dneo_phi_dpsi(iz, 1)) 
             end do  
         else
             wstar(:, :, ivmu) = wstar(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is) 
         end if
      end do

      deallocate (energy)

      ! Read out wstar using the wpol diagnostic. 
      call write_wpol_diagnostic(wstar)

   end subroutine init_wstar

! ================================================================================================================================================================================== !
! --------------------------------------------------------------- If NEO's corrections are enabled, initialise wpol. --------------------------------------------------------------- !
! ================================================================================================================================================================================== !

    subroutine init_wpol
        ! Parallelisation.
        use mp, only: mp_abort, proc0
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      
        ! Grids.
        use grids_time, only: code_dt
        use grids_kxky, only: nalpha
        use grids_z, only: nzgrid
        use grids_species, only: spec
        use grids_velocity, only: vperp2, vpa
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
        use geometry, only: clebsch_factor, dxdpsi            
        use geometry_miller, only: local

        use neoclassical_terms_neo, only: neoclassical_is_enabled, dneo_h_dz, dneo_phi_dz
        use neoclassical_diagnostics, only: write_wpol_diagnostic

        use arrays, only: wpol, initialised_wpol

        ! Rescale the drive term with <wstarknob>.
        use parameters_physics, only: wstarknob

        implicit none

        ! Indices.
        integer :: is, imu, iv, ivmu, iz
        
        ! If neoclassical terms are NOT enabled, exit the subroutine immediately. 
        if (.not. neoclassical_is_enabled()) return
 
        ! Only intialise omega_{pol,k,s} once.
        if (initialised_wpol) return
        initialised_wpol = .true.

        ! Allocate omega_{pol,k,s} = wpol[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(wpol)) then
            allocate (wpol(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wpol = 0.0
        end if
      
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
         
            do iz = -nzgrid, nzgrid
                wpol(:, iz, ivmu) = 0 
                ! (1/clebsch_factor) * (1/local%qinp) * dxdpsi * wstarknob * 0.5 * code_dt &
                ! * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) * (dneo_h_dz(iz, ivmu, 1) - spec(is)%z * dneo_phi_dz(iz, 1))
            end do
        end do 
       
        ! Read out wpol array for inspection. 
        ! call write_wpol_diagnostic(wpol)

    end subroutine init_wpol

! ================================================================================================================================================================================== !
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- !
! ================================================================================================================================================================================== !


   !*****************************************************************************
   !                           ADVANCE DRIVE TERM(S)
   !*****************************************************************************

   ! Advances the wstar term. 
   subroutine advance_wstar_explicit(phi, gout)
   
      ! Grids
      use grids_z, only: nzgrid
      use parallelisation_layouts, only: vmu_lo
      
      ! Physics
      use parameters_physics, only: full_flux_surface
   
      implicit none
   
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      !-------------------------------------------------------------------------
   
      if (full_flux_surface) then
         call advance_wstar_explicit_ffs(gout)
      else 
         call advance_wstar_explicit_flux_tube(phi, gout)
      end if
   
   end subroutine advance_wstar_explicit


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !


   !-------------------------------- Flux tube ---------------------------------
   subroutine advance_wstar_explicit_flux_tube(phi, gout)

      ! Parallelisation
      use mp, only: proc0
      
      ! Data arrays
      use arrays, only: wstar
      use arrays_fields, only: apar, bpar
      
      ! Grids
      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, nakx
      
      ! Calculations
      use calculations_add_explicit_terms, only: add_explicit_term
      use calculations_kxky_derivatives, only: get_dchidy
      
      ! Time this routine
      use timers, only: time_gke
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

      ! Add the drive term to the right-hand-side of the gyrokinetic equation. 
      if (debug) write (*, *) 'time_advance::solve_gke::add_wstar_term'
      call add_explicit_term(g0, wstar(1, :, :), gout)

      ! Deallocate <g0> = i ky J_0 ϕ_k.
      deallocate (g0)

      

      ! Stop timing the time advance due to the driving gradients
      if (proc0) call time_message(.false., time_gke(:, 6), ' wstar advance')

   end subroutine advance_wstar_explicit_flux_tube
   
   !---------------------------- Full flux surface -----------------------------
   subroutine advance_wstar_explicit_ffs(gout)

      use mp, only: proc0, mp_abort
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo
      use calculations_transforms, only: transform_ky2y
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, naky_all, nakx, ikx_max, ny
      use calculations_kxky, only: swap_kxky
      use arrays, only: wstar
      use arrays_distribution_function, only: phi_gyro
      use calculations_gyro_averages, only: gyro_average
      use calculations_add_explicit_terms, only: add_explicit_term_ffs
      use calculations_kxky_derivatives, only: get_dgdy, get_dchidy
      use timers, only: time_gke

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
      ! Here phi_gyro is <phi> in k-space that has been pre-calculated and stored
      call get_dgdy(phi_gyro, g0)
      
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


! =============================================================================================================================================================================== !
! ---------------------------------------------------- Advance the wpol contribution explicitly if NEO's corrections are included ----------------------------------------------- !
! =============================================================================================================================================================================== !

   subroutine advance_wpol_explicit(phi, gout)
      ! Parallelisation.
      use mp, only: proc0
      
      ! Data arrays.
      use arrays, only: wpol
      use arrays_fields, only: apar, bpar
      
      ! Grids.
      use parallelisation_layouts, only: vmu_lo, is_idx, iv_idx, imu_idx
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, nakx
      
      ! Calculations.
      use calculations_add_explicit_terms, only: add_explicit_term
      use calculations_kxky_derivatives, only: get_dchidx

      use neoclassical_diagnostics, only: write_dchidx_diagnostic_in_advance_wpol_routine

      ! Currently this subroutine is not timed. 
      ! use timers, only: time_gke
      ! use job_manage, only: time_message

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
      complex, dimension(:, :, :, :, :), allocatable :: g0
   
      ! Start timing the time advance due to the poloidal drive.
      ! if (proc0) call time_message(.false., time_gke(:, 6), ' wpol advance')

      if (proc0) then
          print *, "DEBUG: Entering advance_wpol_explicit"
          print *, "DEBUG: phi size: ", size(phi)
          print *, "DEBUG: phi shape: ", shape(phi)
          ! Check if phi contains any NaNs or actual values
          print *, "DEBUG: phi max abs val: ", maxval(abs(phi))
          print *, "DEBUG: phi sum: ", sum(phi)
      end if

      ! Allocate temporary array for <g0> = i ky J_0 ϕ_k.
      allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (proc0) print *, "g0 allocated successfully in the wpol time advance."

      call get_dchidx(phi, apar, bpar, g0)

      if (proc0) print *, "dchidx_4d called successfully."
 
      ! DCHIDX DIAGNOSTIC. 

      call write_dchidx_diagnostic_in_advance_wpol_routine(g0)
     
      print *, "DEBUG Step 1: gout (RHS) max BEFORE: ", maxval(abs(gout))

      ! Add the poloidal drive term to the right-hand-side of the gyrokinetic equation.
      if (debug) write (*, *) 'time_advance::solve_gke::add_wpol_term'
      call add_explicit_term(g0, wpol(1, :, :), gout)
      if (proc0) print *, "gout added to the GKE equation."

      print *, "DEBUG Step 1: gout (RHS) max AFTER: ", maxval(abs(gout))

      ! Deallocate <g0> = i kx J_0 ϕ_k.
      deallocate (g0)

      ! Stop timing the time advance due to the driving gradients.
      ! if (proc0) call time_message(.false., time_gke(:, 6), ' wpol advance')

   end subroutine advance_wpol_explicit

! =============================================================================================================================================================================== !
! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- !
! =============================================================================================================================================================================== !


   !*****************************************************************************
   !                           FINALISE DRIVE TERM(S)
   !*****************************************************************************
   subroutine finish_wstar

      use arrays, only: wstar, initialised_wstar

      implicit none

      if (allocated(wstar)) deallocate (wstar)
      initialised_wstar = .false.

   end subroutine finish_wstar
   
! =============================================================================================================================================================================== !
! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- !
! =============================================================================================================================================================================== !

    subroutine finish_wpol

      use neoclassical_terms_neo, only: neoclassical_is_enabled
      use arrays, only: wpol, initialised_wpol

      implicit none

      ! If neoclassical terms are NOT enabled, exit the subroutine immediately. 
      if (.not. neoclassical_is_enabled()) return

      if (allocated(wpol)) deallocate (wpol)
      initialised_wpol = .false.

   end subroutine finish_wpol


! =============================================================================================================================================================================== !
! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- !
! =============================================================================================================================================================================== !

end module gk_drive
