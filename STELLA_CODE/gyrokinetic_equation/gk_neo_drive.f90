! ================================================================================================================================================================================= !
! ------------------------------------------------------------------- Evolves the neoclassical gradient drive. -------------------------------------------------------------------- !​
! ================================================================================================================================================================================= !
! 
! This module evolves the neoclassical gradient drive. In the conventional theory, this has only a ψ derivative component, which is captured by wstar defined in gk_drive.f90.
! When F_1 is included, the gradient drive has both a ψ derivative component and a θ (poloidal angle) derivative component. When using a miller geometry, this is equivalent to a 
! z derivative component.     
!         
! In this routine we define two coeffecients, wstar1 and wpol. wstar1 is the higher-order counterpart to wstar. 
! 
! 
!
! wstar1 must be multiplied by ik_y * <Χ_k> and then added to the RHS of the GKE. wpol must be multiplied by ik_x * <Χ_k> and added to the RHS of the GKE. 
! 
! ================================================================================================================================================================================= !

module gk_neo_drive

   ! Load debug flags.
   ! use debug_flags, only: debug => neo_drive_debug
   
   implicit none

   ! Make routines available to other modules. 
   public :: initialised_wstar1, initialised_wpol
   public :: init_wstar1, init_wpol
   public :: finish_wstar1, finish_wpol 
   public :: advance_wstar1_explicit, advance_wpol_explicit 

   private
   
   ! Only initialise once.
   logical :: initialised_wstar1 = .false.
   logical :: initialised_wpol   = .false.

contains

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------- Initialise wstar1. ------------------------------------------------------------------------------ ! 
! ================================================================================================================================================================================= !

    subroutine init_wstar1
        ! Parallelisation
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      
        ! Grids
        use grids_time, only: code_dt
        use grids_kxky, only: nalpha
        use grids_z, only: nzgrid, zed
        use grids_species, only: spec
        use grids_velocity, only: vperp2, vpa
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
        use geometry, only: dydalpha, drhodpsi, clebsch_factor, dxdpsi
        use geometry, only: bmag, gradx_dot_grady
        use geometry_miller, only: local 

        use neoclassical_terms_neo, only:  neo_h, neo_phi             
        use neoclassical_terms_neo, only: dneo_h_dpsi, dneo_phi_dpsi   
        use neoclassical_terms_neo, only: dneo_h_dz, dneo_phi_dz

        use arrays, only: wstar1, initialised_wstar1

        ! Rescale the drive term with <wstar1knob>.
        ! use parameters_physics, only: wstarknob

        implicit none

        ! Indices.
        integer :: is, imu, iv, ivmu, iz
      
        ! To make the calculations easier to follow, we calculate <energy> = v_parallel² + 2 mu B for each velocity point.
        real, dimension(:, :), allocatable :: energy
         
        ! Only intialise omega_{*,k,s,1} once.
        if (initialised_wstar1) return
        initialised_wstar1 = .true.

        ! Allocate omega_{*,k,s,1} = wstar1[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(wstar1)) then
            allocate (wstar1(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstar1 = 0.0
        end if

        ! Allocate <energy>[ialpha,iz].
        allocate (energy(nalpha, -nzgrid:nzgrid))
      
        ! Iterate over velocity space
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
         
            ! Calculate <energy>[ialpha,iz] = v_parallel² + 2 mu B = vpa(iv)**2 + vperp2(ialpha, iz, imu).
            energy(:, :) = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
         
            ! Calculate wstar = - code_dt*omega_{*,k,s}/ky = - code_dt * 0.5/C <dydalpha> exp(-v²) (drho/dpsi) d ln F_s / d rho 
            !  = - code_dt * 0.5/<clebsch_factor> <dydalpha> <drhodpsi> [<fprim> + <tprim> (v_parallel² + 2 mu B - 1.5)] exp(-v²). 
            ! This block only computes when sfincs is chosen for the neoclassical option.

            wstar1(:, :, ivmu) = - (1/clebsch_factor) * dydalpha * drhodpsi * 0.5 * code_dt * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)            

            do iz = -nzgrid, nzgrid
                wstar1(:, iz, ivmu) = wstar1(:, iz, ivmu) * ( ( spec(is)%fprim + spec(is)%tprim * ( energy(:, iz) - 1.5 ) ) * ( neo_h(iz, ivmu, 1) - spec(is)%z * neo_phi(iz, 1) ) &
                - ( dneo_h_dpsi(iz, ivmu, 1) - spec(is)%z * dneo_phi_dpsi(iz, 1) ) )
            end do  

            ! We must add the remaining corrections to wstar1, namely those proportinal to the z derivative of F_1. 

            ! do iz = -nzgrid, nzgrid
                ! wstar1(:, iz, ivmu) = 
            ! end do
      end do

      deallocate (energy)

    end subroutine init_wstar1


! ================================================================================================================================================================================= !
! -------------------------------------------------------------------------------- Initialise wpol. ------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_wpol
        implicit none
    end subroutine init_wpol               


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------- Advance wstar1 explicitly. -------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_wstar1_explicit(phi, gout)
        ! Parallelisation.
        use mp, only: proc0
      
        ! Data arrays.
        use arrays, only: wstar1
        use arrays_fields, only: apar, bpar
      
        ! Grids.
        use parallelisation_layouts, only: vmu_lo
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx
      
        ! Calculations.
        use calculations_add_explicit_terms, only: add_explicit_term
        use calculations_kxky_derivatives, only: get_dchidy
      
        ! Time this routine.
        use timers, only: time_gke
        use job_manage, only: time_message

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout
        complex, dimension(:, :, :, :, :), allocatable :: g0
         
        ! ========================================================================================== !
        ! ------------------------------------------------------------------------------------------ !
        ! ========================================================================================== !
        !                                                                                            ! 
        ! Add the neoclassical  drive term to the GKE:                                               !
        !                                                                                            ! 
        ! -i omega_{*,k,s, 1} _k = <wstar1> * i ky <Χ_k>                                             !
        !                                                                                            !
        ! First we calculate i ky * <Χ_k> which corresponds to d<chi>_theta/dy:                      !
        !                                                                                            ! 
        ! Then multiply with <wstar> and add it to the right-hand-side of the GKE:                   !
        !                                                                                            !
        ! add_explicit_term(g0, wstar1(1, :, :), gout)                                               !
        !                                                                                            ! 
        ! ========================================================================================== !
        ! ------------------------------------------------------------------------------------------ !
        ! ========================================================================================== !

        ! Start timing the time advance due to the driving gradient.
        if (proc0) call time_message(.false., time_gke(:, 6), ' wstar1 advance')

        ! Allocate temporary array for <g0> = i ky J_0 ϕ_k.
        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      
        ! Calculate <g0> = i ky J_0 ϕ_k = d<chi>_theta/dy.
        ! if (debug) write (*, *) 'time_advance::solve_gke::get_dchidy'
        call get_dchidy(phi, apar, bpar, g0)

        ! Add the drive term to the right-hand-side of the gyrokinetic equation. 
        ! if (debug) write (*, *) 'time_advance::solve_gke::add_wstar1_term'
        call add_explicit_term(g0, wstar1(1, :, :), gout)

        ! Deallocate <g0> = i ky J_0 ϕ_k.
        deallocate (g0)

        ! Stop timing the time advance due to the driving gradient.
        if (proc0) call time_message(.false., time_gke(:, 6), ' wstar1 advance')
 
    end subroutine advance_wstar1_explicit


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------- Advance wpol explicitly. --------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_wpol_explicit
        implicit none
    end subroutine advance_wpol_explicit


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------- Finish wstar1. -------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_wstar1
        use arrays, only: wstar1

        implicit none

        if (allocated(wstar1)) deallocate (wstar1)
        initialised_wstar1 = .false.

    end subroutine finish_wstar1


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------------- Finish wpol. --------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_wpol   
        use arrays, only: wpol       

        implicit none

        if (allocated(wpol)) deallocate (wpol)
        initialised_wpol = .false.

    end subroutine finish_wpol   


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_drive
