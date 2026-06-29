 ! ================================================================================================================================================================================= !
! -------------------------------------------------------------- Evolves neoclassical corrections proportional to <A∥_k>. --------------------------------------------------------- !​
! ================================================================================================================================================================================= !
! 
! This module evolves the following higher order neoclassical corrections: 
!
! = 
!         
! Define the neoclassical apar coefficient as: 
! 
! <neo_apar_coeff> =
!
! This must be multiplied by <A∥_k> = J₀ * A∥_k and then added to the RHS of the GKE.
! 
! ================================================================================================================================================================================= !

module gk_neo_mirror
 
   implicit none

   ! Make routines available to other modules. 
   public :: initialised_neo_mirror
   public :: init_neo_mirror, finish_neo_mirror
   public :: advance_neo_mirror_explicit

   private
   
   ! Only initialise once.
   logical :: initialised_neo_mirror = .false.

contains

! ================================================================================================================================================================================= !
! -------------------------------------------------------------------- Initialise the neoclassical A∥_k terms. -------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_neo_mirror
        ! Parallelisation.
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        ! Grids.
        use grids_time, only: code_dt
        use grids_species, only: spec, nspec
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
        use grids_velocity, only: mu, vpa, nmu, nvpa, vperp2
        use grids_z, only: nzgrid, nztot
        use grids_kxky, only: nalpha

        ! Geometry.
        use geometry, only: bmag, dbdzed, b_dot_gradz

        ! Neoclassical.
        use neoclassical_terms_neo, only: neo_vpa_fac, neo_mu_fac
        use neoclassical_terms_neo, only: neo_mu_fac_global

        use arrays, only: neo_mirror, initialised_neo_mirror

        ! For switching streaming on and off.
        use parameters_physics, only: neomirrorknob

        implicit none

        ! Local variables. 
        integer :: iz, iv, is, imu, ivmu
        real, dimension(:, :, :), allocatable :: d2dF1_dvpadmu

        ! Only intialise once.
        if (initialised_neo_mirror) return
        initialised_neo_mirror = .true.

        ! Allocate neo_mirror = neo_mirror[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_mirror)) then
            allocate (neo_mirror(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_mirror = 0.0
        end if
 
        ! Allocate the mixed derivative in F_1.
        allocate(d2dF1_dvpadmu(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc, 1)); d2dF1_dvpadmu = 0.0

        ! Get the mixed derivative of F_1, parallelised over the velocity space. 
        call get_vpa_derivative_explicit(neo_mu_fac_global, d2dF1_dvpadmu)

        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid  
                neo_mirror(:, iz, ivmu) = neomirrorknob * code_dt * spec(is)%z * mu(imu) * b_dot_gradz(:, iz) * dbdzed(:, iz) &
                * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) / spec(is)%mass 

                neo_mirror(:, iz, ivmu) = neo_mirror(:, iz, ivmu) &
                * ( neo_vpa_fac(iz, ivmu, 1) / vpa(iv) - neo_mu_fac(iz, ivmu, 1) / bmag(:, iz) - vpa(iv) * d2dF1_dvpadmu(iz, ivmu, 1) / bmag(:, iz) )
            end do 
        end do

        ! Deallocate temporary array. 
        deallocate(d2dF1_dvpadmu)

    end subroutine init_neo_mirror

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------- Advance the terms explicitly. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_mirror_explicit(apar, gout)
        ! Parallelisation.
        use mp, only: proc0
        use parallelisation_layouts, only: vmu_lo
        use parallelisation_layouts, only: is_idx, iv_idx, imu_idx

      
        ! Data arrays.
        use arrays, only: neo_mirror

        ! Grids. 
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx
      
        ! Calculations.
        use calculations_add_explicit_terms, only: add_explicit_term
        use calculations_gyro_averages, only: gyro_average

        ! Time this routine.
        use timers, only: time_gke
        use job_manage, only: time_message

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: apar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout        

        ! Local variables.
        integer :: ivmu, iv, is, imu
        complex, dimension(:, :, :, :), allocatable :: field
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
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_mirror advance')

        ! Allocate temporary array for <g0> = <A∥_k>.
        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes))
 
        ! Construct <A∥_k>.
        ! Iterate over the (mu,vpa,s) points. 
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)

            ! Calculate the apar field. 
            field = apar

            ! Gyroaverage.
            call gyro_average(field, ivmu, g0(:, :, :, :, ivmu))
        end do

        ! Add the term to the right-hand-side of the GKE. 
        call add_explicit_term(g0, neo_mirror(1, :, :), gout)

        ! Deallocate <g0>.
        deallocate (g0)
        deallocate (field)

        ! Stop timing the time advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_mirror advance')

    end subroutine advance_neo_mirror_explicit


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------- Finish the mirror correction. ----------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_mirror
        use arrays, only: neo_mirror, initialised_neo_mirror

        implicit none

        if (allocated(neo_mirror)) deallocate (neo_mirror)
        initialised_neo_mirror = .false.

    end subroutine finish_neo_mirror

! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------- Utilities. ------------------------------------------------------------------------------------ ! 
! ================================================================================================================================================================================= !

! ================================================================================================================================================================================= !
! -- Get the vpa derivative of an array that is local in the velocity data - needed for the explicit mirror advance when apar is included. Uses a third order upwind scheme based - !
! --------------------------------------------------------------------- on the sign of the mirror coeffecient. -------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine get_vpa_derivative_explicit(g, dgdvpa)
        ! Parallelisation. 
        use parallelisation_layouts, only: vmu_lo
        use neoclassical_terms_neo, only: distribute_vmus_over_procs
    
        ! Caclculations. 
        use calculations_finite_differences, only: third_order_upwind

        ! Grids.
        use grids_z, only: nzgrid
        use grids_velocity, only: dvpa, nmu, nvpa
        use grids_species, only: nspec

        ! Conventional mirror term. 
        use gk_mirror, only: mirror_sign

        implicit none

        real, dimension(-nzgrid:, :, :, :, :), intent(in) :: g
        real, dimension(-nzgrid:, vmu_lo%llim_proc:, :), intent(in out) :: dgdvpa

        integer :: iv, imu, is, iz, ia
        real, dimension(:, :, :, :, :), allocatable :: dgdvpa_global
        real, dimension(:), allocatable :: tmp

        ia = 1
        
        ! =================================================================== !

        allocate (tmp(nvpa))
        allocate (dgdvpa_global(-nzgrid:nzgrid, nvpa, nmu, nspec, 1))

        do iz = -nzgrid, nzgrid
            do imu = 1, nmu
                do is = 1, nspec
                    call third_order_upwind(1, g(iz, :, imu, is, 1), dvpa, mirror_sign(1, iz), tmp)

                    dgdvpa_global(iz, :, imu, is, 1) = tmp
                end do
            end do
        end do

        ! Parallelise over the velocity space.
        do iz = -nzgrid, nzgrid
            call distribute_vmus_over_procs(dgdvpa_global(iz, :, :, :, 1), dgdvpa(iz, :, 1))
        end do

        deallocate(tmp)
        deallocate(dgdvpa_global)

   end subroutine get_vpa_derivative_explicit


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_mirror
