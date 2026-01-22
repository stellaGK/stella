 ! ================================================================================================================================================================================= !
! -------------------------------------------------------------- Evolves neoclassical corrections proportional to <A∥_k>. --------------------------------------------------------- !​
! ================================================================================================================================================================================= !
! 
! This module evolves the following higher order neoclassical corrections: 
!
! = (Z * μ)/(m * v∥) * b.∇B * J₀ * A∥_k * exp(-v²) * ( ∂H_1/∂v∥|_μ - 2v∥ * ( H_1 - e * ϕ₀¹ ) )            
!         
! Define the neoclassical apar coefficient as: 
! 
! <neo_apar_coeff> = (Z * μ)/(m * v∥) * b.∇B * exp(-v²) * ( ∂H_1/∂v∥|_μ - 2v∥ * ( H_1 - e * ϕ₀¹ ) )  * code_dt
!
! This must be multiplied by <A∥_k> = J₀ * A∥_k and then added to the RHS of the GKE.
! 
! ================================================================================================================================================================================= !

module gk_neo_apar_terms

   ! Load debug flags.
   ! use debug_flags, only: debug => neo_apar_terms_debug
   
   implicit none

   ! Make routines available to other modules. 
   public :: initialised_neo_apar_terms
   public :: init_neo_apar_terms, finish_neo_apar_terms
   public :: advance_neo_apar_terms_explicit

   private
   
   ! Only initialise once.
   logical :: initialised_neo_apar_terms = .false.

contains

! ================================================================================================================================================================================= !
! -------------------------------------------------------------------- Initialise the neoclassical A∥_k terms. -------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_neo_apar_terms
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

        use arrays, only: neo_apar_coeff, initialised_neo_apar_terms

        implicit none

        integer :: iz, iv, is, imu, ivmu

        ! Only intialise once.
        if (initialised_neo_apar_terms) return
        initialised_neo_apar_terms = .true.

        ! Allocate neo_apar_coeff = neo_apar_coeff[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_apar_coeff)) then
            allocate (neo_apar_coeff(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_apar_coeff = 0.0
        end if
        
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
         
            ! Calculate neo_apar_coeff at each grid point. Calculation is broken up for ease of reading.  
            do iz = -nzgrid, nzgrid
                ! First compute the magnetic geometry prefactor. 
                neo_apar_coeff(:, iz, ivmu) = b_dot_gradz(:, iz) * dbdzed(:, iz)

                ! Multiply by the species dependent prefactor. 
                neo_apar_coeff(:, iz, ivmu) = neo_apar_coeff(:, iz, ivmu) * spec(is)%z/spec(is)%mass &
                * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is)

                ! Multiply by the μ grid point and divide by the v∥ grid point.
                neo_apar_coeff(:, iz, ivmu) = neo_apar_coeff(:, iz, ivmu) * ( mu(imu)/vpa(iv) )

                ! Multiply by the neoclassical distribution prefactor. 
                neo_apar_coeff(:, iz, ivmu) = neo_apar_coeff(:, iz, ivmu) &
                * ( dneo_h_dvpa(iz, ivmu, 1) - 2 * vpa(iv) * ( neo_h(iz, ivmu, 1) - spec(is)%z * neo_phi(iz, 1) ) )

                ! Finally, multiply by code_dt. 
                neo_apar_coeff(:, iz, ivmu) = neo_apar_coeff(:, iz, ivmu) * code_dt
            end do 
        end do

    end subroutine init_neo_apar_terms

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------- Advance the terms explicitly. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_apar_terms_explicit(gout)
        ! Parallelisation.
        use mp, only: proc0
        use parallelisation_layouts, only: vmu_lo
      
        ! Data arrays.
        use arrays, only: neo_apar_coeff
        use arrays_fields, only: apar

        ! Grids. 
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx
      
        ! Calculations.
        use calculations_add_explicit_terms, only: add_explicit_term

        ! Time this routine.
        use timers, only: time_gke
        use job_manage, only: time_message

        implicit none

        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout        
        complex, dimension(:, :, :, :, :), allocatable :: g0 

        ! ======================================================================================= ! 
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= !
        !                                                                                         !
        ! Calculate the gyroaveraged <A∥_k>:                                                      ! 
        !                                                                                         !
        ! <g0> = <A∥_k>                                                                           !
        !                                                                                         !
        ! Mutlipy this by neo_apar_coeff and add to the right-hand-side of the GKE:               !
        !                                                                                         ! 
        ! add_explicit_term(g0, neo_apar_coeff(1, :, :), gout)                                    !
        !                                                                                         !
        ! ======================================================================================= !
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= !

        ! Start timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_apar_coeff advance')

        ! Allocate temporary array for <g0> = <A∥_k>.
        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
 
        ! Construct <A∥_k> potential. 
        call get_apar(apar, g0)        
        
        ! Add the term to the right-hand-side of the GKE. 
        call add_explicit_term(g0, neo_apar_coeff(1, :, :), gout)

        ! Deallocate <g0>.
        deallocate (g0)

        ! Stop timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_apar_coeff advance')

    end subroutine advance_neo_apar_terms_explicit


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------- Finish the A∥_k terms. ----------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_apar_terms
        use arrays, only: neo_apar_coeff, initialised_neo_apar_terms

        implicit none

        if (allocated(neo_apar_coeff)) deallocate (neo_apar_coeff)
        initialised_neo_apar_terms = .false.

    end subroutine finish_neo_apar_terms


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------------- Utilities. ----------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------ Calculate the gyroaveraged A∥_k. ----------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine get_apar(apar, gyro_apar)      
        ! Parallelisation.
        use parallelisation_layouts, only: vmu_lo
        use parallelisation_layouts, only: is_idx, iv_idx, imu_idx
      
        ! Flags.
        use parameters_physics, only: include_apar
      
        ! Grids.
        use grids_species, only: spec
        use grids_z, only: nzgrid, ntubes
        use grids_velocity, only: vpa, mu
        use grids_kxky, only: nakx, naky, aky
      
        ! Calculations.
        use calculations_gyro_averages, only: gyro_average

        implicit none

        ! Arguments.
        complex, dimension(:, :, -nzgrid:, :), intent(in)                     :: apar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: gyro_apar

        ! Local variables.
        integer :: ivmu, iv, is, imu
        complex, dimension(:, :, :, :), allocatable :: field

        ! Allocate temporary array for <gyro_apar> = A∥_k. 
        allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes))
 
        ! Iterate over the (mu,vpa,s) points. 

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)

            ! Calculate the apar field. 
            field = apar

            ! Gyroaverage.
            call gyro_average(field, ivmu, gyro_apar(:, :, :, :, ivmu))
        end do

        ! Deallocate temporary array.
        deallocate (field)

    end subroutine get_apar

! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_apar_terms
