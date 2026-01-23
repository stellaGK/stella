! ================================================================================================================================================================================= !
! --------------------------------- Evolves neoclassical corrections proportional to the z derivative of the gyroaveraged generalised, ∂<Χ_k>/∂z. --------------------------------- !​
! ================================================================================================================================================================================= !
! 
! This module evolves the following higher order neoclassical corrections: 
!
! = 1/2 * Z/T * v_{th,s} * b.∇z * exp(-v²) * ( ∂H_1/∂v∥|_μ - 2v∥ * ( H_1 - e * ϕ₀¹ ) ) * ∂<Χ_k>/∂z            
!         
! Define the neoclassical ∂<Χ_k>/∂z coefficient as: 
! 
! <neo_dchidz_coeff> = 1/2 * Z/T * v_{th,s} * b.∇z * exp(-v²) * ( ∂H_1/∂v∥|_μ - 2v∥ * ( H_1 - e * ϕ₀¹ ) ) * code_dt
!
! This must be multiplied by ∂<Χ_k>/∂z and then added to the RHS of the GKE.
! 
! ================================================================================================================================================================================= !

module gk_neo_dchidz_terms

   ! Load debug flags.
   ! use debug_flags, only: debug => neoclassical_dchidz_terms_debug
   
   implicit none

   ! Make routines available to other modules. 
   public :: initialised_neo_dchidz_terms
   public :: init_neo_dchidz_terms, finish_neo_dchidz_terms
   public :: advance_neo_dchidz_terms_explicit, advance_neo_dchidz_terms_implicit 

   private
   
   ! Only initialise once.
   logical :: initialised_neo_dchidz_terms = .false.

contains

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------ Initialise the neoclassical ∂<Χ_k>/∂z terms. ----------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_neo_dchidz_terms
        ! Parallelisation.
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        use grids_time, only: code_dt
        use grids_species, only: spec, nspec
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
        use grids_velocity, only: mu, vperp2, vpa
        use grids_z, only: nzgrid, nztot
        use grids_kxky, only: nalpha

        use geometry, only: bmag, dbdzed, b_dot_gradz

        use neoclassical_terms_neo, only: dneo_h_dvpa, neo_h, neo_phi

        use parameters_physics, only: include_apar

        use arrays, only: neo_dchidz_coeff, initialised_neo_dchidz_terms

        implicit none

        integer :: iz, iv, is, imu, ivmu

        ! Only intialise once.
        if (initialised_neo_dchidz_terms) return
        initialised_neo_dchidz_terms = .true.

        ! Allocate neo_dchidz_coeff = neo_dchidz_coeff[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_dchidz_coeff)) then
            allocate (neo_dchidz_coeff(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_dchidz_coeff = 0.0
        end if
        
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
         
            ! Calculate neo_dchidz_coeff at each grid point. Calculation is broken up for ease of reading.  
            do iz = -nzgrid, nzgrid
                ! First compute the magnetic geometry prefactor. 
                neo_dchidz_coeff(:, iz, ivmu) = 0.5 * b_dot_gradz(:, iz) 

                ! Multiply by the species dependent prefactor. 
                neo_dchidz_coeff(:, iz, ivmu) = neo_dchidz_coeff(:, iz, ivmu) * (spec(is)%z / spec(is)%temp) * maxwell_vpa(iv, is) &
                * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) * spec(is)%stm_psi0

                ! Multiply by the neoclassical distribution prefactor. 
                neo_dchidz_coeff(:, iz, ivmu) = neo_dchidz_coeff(:, iz, ivmu) * ( dneo_h_dvpa(iz, ivmu, 1) - 2 * vpa(iv) &
                * ( neo_h(iz, ivmu, 1) - spec(is)%z * neo_phi(iz, 1) ) ) 

                ! Finally, multiply by code_dt. 
                neo_dchidz_coeff(:, iz, ivmu) = neo_dchidz_coeff(:, iz, ivmu) * code_dt 
            end do 
        end do

    end subroutine init_neo_dchidz_terms

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------- Advance the terms explicitly. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_dchidz_terms_explicit(phi, gout)
        ! Parallelisation.
        use mp, only: proc0
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      
        ! Data arrays.
        use arrays, only: neo_dchidz_coeff
        use arrays_fields, only: apar, bpar      

        ! Grids. 
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx
      
        ! For calculating <Χ_k> and ∂<Χ_k>/∂z. 

        use gk_neo_chi_terms, only: get_chi
        use gk_parallel_streaming, only: get_dgdz_centered

        ! Calculations.
        use calculations_add_explicit_terms, only: add_explicit_term

        ! Time this routine.
        use timers, only: time_gke
        use job_manage, only: time_message

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout        
        complex, dimension(:, :, :, :, :), allocatable :: g0, dg0dz

        integer :: iv, is, imu, ivmu

        ! ======================================================================================= ! 
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= ! 
        !                                                                                         !
        ! Calculate the gyrokinetic potential:                                                    ! 
        !                                                                                         !
        ! <g0> = <Χ_k>                                                                            !
        !                                                                                         !
        ! Calculate the parallel derivative:                                                      !
        !                                                                                         !
        ! <dg0dz> = ∂<Χ_k>/∂z                                                                     !
        !                                                                                         !
        ! Mutlipy this by neo_dchidz_coeff and add to the right-hand-side of the GKE:             !
        !                                                                                         ! 
        ! add_explicit_term(g0, neo_dchidz_coeff(1, :, :), gout)                                  !
        !                                                                                         ! 
        ! ======================================================================================= !
        ! --------------------------------------------------------------------------------------- ! 
        ! ======================================================================================= !

        ! Start timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_dchidz_coeff advance')

        ! Allocate temporary array for <g0> = <Χ_k> and <dg0dz> = ∂<Χ_k>/∂z.
        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        allocate (dg0dz(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
 
        ! Construct the generalised potential. 
        call get_chi(phi, apar, bpar, g0)        

        ! Calculate the parallel derivative. 
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
         
            ! Note that this should be a centered difference to avoid numerical unpleasantness to do with inexact cancellations in later velocity integration?
            ! See Appendix of the stella JCP 2019 for details?
            call get_dgdz_centered(g0(:, :, :, :, ivmu), ivmu, dg0dz(:, :, :, :, ivmu))
        end do

        ! Add the term to the right-hand-side of the GKE. 
        call add_explicit_term(dg0dz, neo_dchidz_coeff(1, :, :), gout)

        ! Deallocate <g0> and <dg0dz>.
        deallocate (g0)
        deallocate (dg0dz)

        ! Stop timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_dchidz_coeff advance')

    end subroutine advance_neo_dchidz_terms_explicit


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------- Advance the terms implicitly. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_dchidz_terms_implicit
        implicit none
    end subroutine advance_neo_dchidz_terms_implicit


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------- Finish the terms. ------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_dchidz_terms
        use arrays, only: neo_dchidz_coeff, initialised_neo_dchidz_terms

        implicit none

        if (allocated(neo_dchidz_coeff)) deallocate (neo_dchidz_coeff)
        initialised_neo_dchidz_terms = .false.

    end subroutine finish_neo_dchidz_terms


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------------- Utilities. ----------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_dchidz_terms
