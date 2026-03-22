! ================================================================================================================================================================================= !
! -------------------------------------------- Evolves neoclassical parallel streaming corrections proportional to the gyroaveraged generalised, <Χ_k>. --------------------------- !​
! ================================================================================================================================================================================= !
! 
! This module evolves the following higher order neoclassical corrections: 
!
! = 
!      
! Define the neoclassical stremaing coefficient as: 
! 
! <neoclassical_stream_coeff> = 
!
! This must be multiplied by <Χ_k> and then added to the RHS of the GKE.
!
! ================================================================================================================================================================================= !

module gk_neo_parallel_streaming

   ! Load debug flags.
   ! use debug_flags, only: debug => neoclassical_stream_terms_debug
   
   implicit none

   ! Make routines available to other modules. 
   public :: initialised_neo_stream
   public :: init_neo_stream, finish_neo_stream
   public :: advance_neo_stream_explicit

   private
   
   ! Only initialise once.
   logical :: initialised_neo_stream = .false.

contains

! ================================================================================================================================================================================= !
! -------------------------------------------------------------------- Initialise the neoclassical Χ_k terms. --------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_neo_stream
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

        use neoclassical_terms_neo, only: neo_phi, neo_h, dneo_phi_dz, dneo_h_dz, d2neo_h_dzdmu

        use arrays, only: neo_stream_coeff, initialised_neo_stream

        implicit none

        integer :: iz, iv, is, imu, ivmu

        ! Allocate neo_stream_coeff = neo_stream_coeff[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_stream_coeff)) then
            allocate (neo_stream_coeff(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_stream_coeff = 0.0
        end if

        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
         
            ! Calcualte the species dependent factor. 
            neo_stream_coeff(:, :, ivmu) = spec(is)%stm * spec(is)%zt 

            ! Multiply by the z-dependent factor. 
            do iz = -nzgrid, nzgrid
                neo_stream_coeff(:, iz, ivmu) = - neo_stream_coeff(:, iz, ivmu) * ( 0.5 * vpa(iv) / bmag(:, iz) ) &
                * b_dot_gradz(:, iz) * ( d2neo_h_dzdmu(iz, ivmu, 1) - 2 * dbdzed(:, iz) * ( neo_h(iz, ivmu, 1) - spec(is)%z * neo_phi(iz) ) & 
                - 2 * bmag(:, iz) * ( dneo_h_dz(iz, ivmu, 1) - spec(is)%z * dneo_phi_dz(iz) ) ) &
                * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is)
            end do 
        end do

        neo_stream_coeff = neo_stream_coeff * code_dt

    end subroutine init_neo_stream


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------- Advance the terms explicitly. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_stream_explicit(phi, gout)
        ! Parallelisation.
        use mp, only: proc0
        use parallelisation_layouts, only: vmu_lo
      
        ! Data arrays.
        use arrays, only: neo_stream_coeff
        use arrays_fields, only: apar, bpar      

        ! Grids. 
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx
      
        ! Calculations.
        use calculations_add_explicit_terms, only: add_explicit_term
        use gk_neo_chi_terms, only: get_chi

        ! Time this routine.
        use timers, only: time_gke
        use job_manage, only: time_message

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout        
        complex, dimension(:, :, :, :, :), allocatable :: g0 

        ! ======================================================================================= ! 
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= !
        !                                                                                         ! 
        ! Calculate the gyrokinetic potential:                                                    ! 
        !                                                                                         !
        ! <g0> = Χ_k                                                                              !
        !                                                                                         !
        ! Mutlipy this by neo_stream_coeff and add to the right-hand-side of the GKE:             !
        !                                                                                         ! 
        ! add_explicit_term(g0, neo_stream_coeff(1, :, :), gout)                                  !
        !                                                                                         !
        ! ======================================================================================= !
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= !

        ! Start timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_stream_coeff advance')

        ! Allocate temporary array for <g0> = J_0 Χ_k.
        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
 
        ! Construct the generalised potential. 
        call get_chi(phi, apar, bpar, g0)        
        
        ! Add the term to the right-hand-side of the GKE. 
        call add_explicit_term(g0, neo_stream_coeff(1, :, :), gout)

        ! Deallocate <g0>.
        deallocate (g0)

        ! Stop timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_stream_coeff advance')

    end subroutine advance_neo_stream_explicit


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------- Finish the terms. ------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_stream
        use arrays, only: neo_stream_coeff, initialised_neo_stream

        implicit none

        if (allocated(neo_stream_coeff)) deallocate (neo_stream_coeff)
        initialised_neo_stream = .false.

    end subroutine finish_neo_stream


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_parallel_streaming
