! ================================================================================================================================================================================= !
! -------------------------------------------- Evolves neoclassical corrections proportional to the gyroaveraged generalised, <Χ_k>. ---------------------------------------------- !​
! ================================================================================================================================================================================= !
! 
! This module evolves the following higher order neoclassical corrections: 
!
! = 1/2B * Z/T * v_{th,s} * b.∇B * exp(-v²) * ( v∥/B * ∂F_1/∂μ|_v∥ - ∂F_1/∂v∥|_μ) * <Χ_k>           
!
! Derivatives coming from the Maxwellian normalisation cancel one another. 
!         
! Define the neoclassical chi coefficient as: 
! 
! <neoclassical_chi_coeff> = 1/2B Z/T * v_{th,s} * b.∇B * exp(-v²) * ( v∥/B * ∂F_1/∂μ|_v∥ - ∂F_1/∂v∥|_μ ) * code_dt
!
! This must be multiplied by <Χ_k> and then added to the RHS of the GKE.
!
! Neoclassical corrections proportional to <Χ_k> and the magnetic curvature drift will probably be handled in a seperate module for easier interpretation.   
! 
! ================================================================================================================================================================================= !

module gk_neo_chi_terms

   ! Load debug flags.
   ! use debug_flags, only: debug => neoclassical_chi_terms_debug
   
   implicit none

   ! Make routines available to other modules. 
   public :: initialised_neo_chi_terms
   public :: init_neo_chi_terms, finish_neo_chi_terms
   public :: advance_neo_chi_terms_explicit, advance_neo_chi_terms_implicit 
   public :: get_chi

   private
   
   ! Only initialise once.
   logical :: initialised_neo_chi_terms = .false.

contains

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------ Initialise the neoclassical Χ_k terms. ----------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_neo_chi_terms
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

        use neoclassical_terms_neo, only: dneo_h_dvpa, dneo_h_dmu

        use parameters_physics, only: include_apar

        use arrays, only: neo_chi_coeff, initialised_neo_chi_terms

        implicit none

        integer :: iz, iv, is, imu, ivmu

        ! Only intialise once.
        if (initialised_neo_chi_terms) return
        initialised_neo_chi_terms = .true.

        ! Allocate neo chi_coeff = neo_chi_coeff[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_chi_coeff)) then
            allocate (neo_chi_coeff(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_chi_coeff = 0.0
        end if
        
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
         
            ! Calculate neo_chi_coeff at each grid point. Calculation is broken up for ease of reading.  
            do iz = -nzgrid, nzgrid
                ! First compute the magnetic geometry prefactor. 
                neo_chi_coeff(:, iz, ivmu) = (0.5/bmag(:, iz)) * b_dot_gradz(:, iz) * dbdzed(:, iz)

                ! Multiply by the species dependent prefactor. 
                neo_chi_coeff(:, iz, ivmu) = neo_chi_coeff(:, iz, ivmu) * (spec(is)%z / spec(is)%temp) * maxwell_vpa(iv, is) &
                * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) * spec(is)%stm_psi0

                ! Multiply by the neoclassical distribution prefactor. 
                neo_chi_coeff(:, iz, ivmu) = neo_chi_coeff(:, iz, ivmu) *  ( ( vpa(iv)/bmag(:, iz) ) * dneo_h_dmu(iz, ivmu, 1) - dneo_h_dvpa(iz, ivmu, 1) )

                ! Finally, multiply by code_dt. 
                neo_chi_coeff(:, iz, ivmu) = neo_chi_coeff(:, iz, ivmu) * code_dt 
            end do 
        end do

    end subroutine init_neo_chi_terms

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------- Advance the terms explicitly. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_chi_terms_explicit(phi, gout)
        ! Parallelisation.
        use mp, only: proc0
        use parallelisation_layouts, only: vmu_lo
      
        ! Data arrays.
        use arrays, only: neo_chi_coeff
        use arrays_fields, only: apar, bpar      

        ! Grids. 
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx
      
        ! Calculations.
        use calculations_add_explicit_terms, only: add_explicit_term

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
        ! Mutlipy this by neo_chi_coeff and add to the right-hand-side of the GKE:                !
        !                                                                                         ! 
        ! add_explicit_term(g0, neo_chi_coeff(1, :, :), gout)                                     !
        !                                                                                         !
        ! ======================================================================================= !
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= !

        ! Start timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_chi_coeff advance')

        ! Allocate temporary array for <g0> = J_0 Χ_k.
        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
 
        ! Construct the generalised potential. 
        call get_chi(phi, apar, bpar, g0)        
        
        ! Add the term to the right-hand-side of the GKE. 
        call add_explicit_term(g0, neo_chi_coeff(1, :, :), gout)

        ! Deallocate <g0>.
        deallocate (g0)

        ! Stop timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_chi_coeff advance')

    end subroutine advance_neo_chi_terms_explicit


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------- Advance the terms implicitly. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_chi_terms_implicit
        implicit none
    end subroutine advance_neo_chi_terms_implicit


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------- Finish the terms. ------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_chi_terms
        use arrays, only: neo_chi_coeff, initialised_neo_chi_terms

        implicit none

        if (allocated(neo_chi_coeff)) deallocate (neo_chi_coeff)
        initialised_neo_chi_terms = .false.

    end subroutine finish_neo_chi_terms


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------------- Utilities. ----------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------- Calculate the gyroaveraged potential. -------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine get_chi(phi, apar, bpar, chi)      
        ! Parallelisation.
        use parallelisation_layouts, only: vmu_lo
        use parallelisation_layouts, only: is_idx, iv_idx, imu_idx
      
        ! Flags.
        use parameters_physics, only: include_apar, include_bpar
        use parameters_physics, only: fphi
      
        ! Grids.
        use grids_species, only: spec
        use grids_z, only: nzgrid, ntubes
        use grids_velocity, only: vpa, mu
        use grids_kxky, only: nakx, naky, aky
      
        ! Calculations.
        use calculations_gyro_averages, only: gyro_average
        use calculations_gyro_averages, only: gyro_average_j1

        implicit none

        ! Arguments.
        complex, dimension(:, :, -nzgrid:, :), intent(in)                     :: phi, apar, bpar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: chi

        ! Local variables.
        integer :: ivmu, iv, is, imu
        complex, dimension(:, :, :, :), allocatable :: field, gyro_tmp

        ! Allocate temporary array for <g0> = Χ_k. 
        allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes))
        allocate (gyro_tmp(naky, nakx, -nzgrid:nzgrid, ntubes))

        ! Construct the generalised potential.
        ! Iterate over the (mu,vpa,s) points.

        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)

            ! Calculate phi.
            field = fphi * phi

            ! If apar is present, we must account for this.
            if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar

            ! Gyroaverage the J_0 contribution.
            call gyro_average(field, ivmu, chi(:, :, :, :, ivmu))

            ! If bpar is present, we must account for this too.
            if (include_bpar) then
                field = 4.0 * mu(imu) * (spec(is)%tz) * bpar
               
                ! Gyroaverage the J_1 contribution.
                call gyro_average_j1(field, ivmu, gyro_tmp)
              
                chi(:, :, :, :, ivmu) = chi(:, :, :, :, ivmu) + gyro_tmp
            end if
        end do

        ! Deallocate temporary arrays.
        deallocate (field)
        deallocate (gyro_tmp)

    end subroutine get_chi

! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_chi_terms
