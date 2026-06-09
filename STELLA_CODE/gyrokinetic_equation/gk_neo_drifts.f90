! ================================================================================================================================================================================= !
! ------------------------------------------------------------ Evolves the neoclassical grad-B and curvature drift terms. --------------------------------------------------------- !​
! ================================================================================================================================================================================= !
!
! This module evolves the neoclassical drift terms on the RHS of the GKE. There will be a magnetic drift (ω_d) term and a curvature drift (ω_κ) term. 
! There will be a contribution proportional to kx and another one proportional to ky. The dimensionless magnetic drift terms are given by
!
! ================================================================================================================================================================================= !

module gk_neo_drifts

   implicit none

   ! Make routines available to other modules. 
   public :: initialised_neo_curv_drift
   public :: init_neo_curv_drift
   public :: finish_neo_curv_drift
   public :: advance_neo_curv_drift_explicit

   private

   ! Only initialise once.
   logical :: initialised_neo_curv_drift = .false.

contains

! ================================================================================================================================================================================= !
! ----------------------------------------------------------------------- Initialise the curvature drift. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_neo_curv_drift
        ! Parallelisation.
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        ! Grids. 
        use grids_time, only: code_dt
        use grids_species, only: spec
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
        use grids_velocity, only: vpa, mu
        use grids_z, only: nzgrid
        use grids_kxky, only: nalpha

        ! Geometry. 
        use geometry, only: bmag
        use geometry, only: B_times_kappa_dot_gradx, B_times_kappa_dot_grady
        use geometry, only: B_times_gradB_dot_gradx, B_times_gradB_dot_grady

        ! Neoclassical. 
        use neoclassical_terms_neo, only: neo_vpa_fac

        ! Arrays. 
        use arrays, only: neocurvx, neocurvy, initialised_neo_curv_drift

        implicit none

        integer :: iz, iv, is, imu, ivmu

        ! Only intialise once.
        if (initialised_neo_curv_drift) return
        initialised_neo_curv_drift = .true.

        ! Allocate neocurvx = neocurvx[ialpha, iz, i[mu,vpa,s]] and neocurvy = neocurvy[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neocurvx)) then
            allocate (neocurvx(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neocurvx = 0.0
        end if
        if (.not. allocated(neocurvy)) then
            allocate (neocurvy(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neocurvy = 0.0
        end if

        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
           
            do iz = -nzgrid, nzgrid            
                neocurvx(:, iz, ivmu) = vpa(iv) * B_times_kappa_dot_gradx(:, iz) + mu(imu) * B_times_gradB_dot_gradx(:, iz) / vpa(iv)
                
                neocurvy(:, iz, ivmu) = vpa(iv) * B_times_kappa_dot_grady(:, iz) + mu(imu) * B_times_gradB_dot_grady(:, iz) / vpa(iv)

                ! Multiply by the neoclassical distribution factor. 
                neocurvx(:, iz, ivmu) = neocurvx(:, iz, ivmu) * 0.5 * code_dt * neo_vpa_fac(iz, ivmu, 1) & 
                * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) / ( bmag(:, iz) ** 2 )

                neocurvy(:, iz, ivmu) = neocurvy(:, iz, ivmu) * 0.5 * code_dt * neo_vpa_fac(iz, ivmu, 1) &
                * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) / ( bmag(:, iz) ** 2 )
            end do
        end do

    end subroutine init_neo_curv_drift         


! ================================================================================================================================================================================= !
! -------------------------------------------------------------------- Advance the curvature drift explicitly. -------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_curv_drift_explicit(phi, gout)
        ! Parallelisation.
        use mp, only: proc0
        use parallelisation_layouts, only: vmu_lo
      
        ! Data arrays.
        use arrays, only: neocurvx, neocurvy
        use arrays_fields, only: apar, bpar      

        ! Grids. 
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx
      
        ! Calculations.
        use calculations_add_explicit_terms, only: add_explicit_term
        use calculations_kxky_derivatives, only: get_dchidx, get_dchidy

        ! Time this routine.
        use timers, only: time_gke
        use job_manage, only: time_message

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout        
        complex, dimension(:, :, :, :, :), allocatable :: g0x, g0y

        ! ======================================================================================= ! 
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= !
        !                                                                                         ! 
        ! Here we define two temporary arrays, g0x and g0y. These will hold                       ! 
        ! <g0x> = ∂<Χ_k>/∂x = i * kx * <Χ_k> and <g0y> = ∂<Χ_k>/∂y = i * ky * <Χ_k> respectively. ! 
        !                                                                                         ! 
        ! get_dchidx(phi, apar, bpar, g0x)                                                        !     
        ! get_dchidy(phi, apar, bpar, g0y)                                                        !
        !                                                                                         ! 
        ! We then multiply g0x by neocurvx and g0y by neocurvy and add the results seperately to  ! 
        ! the RHS of the GKE:                                                                     !
        !                                                                                         ! 
        ! add_explicit_term(g0x, neocurvx(1, :, :), gout)                                         !
        ! add_explicit_term(g0y, neocyrvy(1, :, :), gout)                                         !
        !                                                                                         !
        ! ======================================================================================= !
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= !

        ! Start timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neocurvx and neocurvy advance')

        ! Allocate temporary array for <g0x> and <g0y>.
        allocate (g0x(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        allocate (g0y(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
 
        ! Construct the derivative of the generalised potential. 
        call get_dchidx(phi, apar, bpar, g0x)
        call get_dchidy(phi, apar, bpar, g0y)        
        
        ! Add the terms to the right-hand-side of the GKE by multiplying by the appropriate coeffecients. 
        call add_explicit_term(g0x, neocurvx(1, :, :), gout)
        call add_explicit_term(g0y, neocurvy(1, :, :), gout)

        ! Deallocate <g0x> and <g0y>.
        deallocate (g0x)
        deallocate (g0y)

        ! Stop timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neocurvx and neocurvy advance')
    end subroutine advance_neo_curv_drift_explicit


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------- Finish the curvature drift. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_curv_drift
        use arrays, only: neocurvx, neocurvy

        implicit none

        if (allocated(neocurvx)) deallocate (neocurvx)
        if (allocated(neocurvy)) deallocate (neocurvy)
    end subroutine finish_neo_curv_drift


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_drifts
