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
        use neoclassical_terms_neo, only: dneo_h_dmu_global
        use neoclassical_terms_neo, only: dneo_h_dvpa, dneo_h_dmu

        use arrays, only: neo_mirror_apar_1, neo_mirror_apar_2
        use arrays, only: initialised_neo_mirror

        ! For switching streaming on and off.
        use parameters_physics, only: neomirrorknob

        implicit none

        ! Local variables. 
        integer :: iz, iv, is, imu, ivmu
        real, dimension(:, :, :), allocatable :: d2neo_h_dvpadmu

        ! Only intialise once.
        if (initialised_neo_mirror) return
        initialised_neo_mirror = .true.

        ! Allocate neo_mirror_apar_1 = neo_mirror_apar_1[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_mirror_apar_1)) then
            allocate (neo_mirror_apar_1(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_mirror_apar_1 = 0.0
        end if

        ! Allocate neo_mirror_apar_2 = neo_mirror_apar_2[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_mirror_apar_2)) then
            allocate (neo_mirror_apar_2(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_mirror_apar_2 = 0.0
        end if
 
        ! Allocate the mixed derivative in F_1.
        allocate(d2neo_h_dvpadmu(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc, 1)); d2neo_h_dvpadmu = 0.0

        ! Get the mixed derivative of F_1, parallelised over the velocity space. 
        call get_vpa_derivative_explicit(dneo_h_dmu_global, d2neo_h_dvpadmu)

        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid  
                ! This is the term multipling apar. 
                neo_mirror_apar_1(:, iz, ivmu) = neomirrorknob * code_dt * spec(is)%zt * spec(is)%stm * mu(imu) * b_dot_gradz(:, iz) * dbdzed(:, iz) &
                * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is)  

                neo_mirror_apar_1(:, iz, ivmu) = neo_mirror_apar_1(:, iz, ivmu) &
                * ( vpa(iv) * dneo_h_dmu(iz, ivmu, 1) / bmag(:, iz) - 0.5 * d2neo_h_dvpadmu(iz, ivmu, 1) / bmag(:, iz) + 2.0 * vpa(iv) * neo_vpa_fac(iz, ivmu, 1) )
            end do 

            do iz = -nzgrid, nzgrid
                ! This is the term multipling the vpa derivative of apar. 
                neo_mirror_apar_2(:, iz, ivmu) = neomirrorknob * code_dt * 0.5 * spec(is)%zt * spec(is)%stm * mu(imu) * b_dot_gradz(:, iz) * dbdzed(:, iz) &
                * ( dneo_h_dvpa(iz, ivmu, 1) / vpa(iv) - dneo_h_dmu(iz, ivmu, 1) / bmag(:, iz) ) * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is)
            end do
        end do

        ! Deallocate temporary array. 
        deallocate(d2neo_h_dvpadmu)

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
        use arrays, only: neo_mirror_apar_1, neo_mirror_apar_2

        ! Grids. 
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx
        use grids_species, only: spec
        use grids_velocity, only: vpa
      
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
 
        ! Construct 2 v_th vpa <A∥_k>.
        ! Iterate over the (mu,vpa,s) points. 
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)

            ! Calculate the apar field. 
            field = 2.0 * spec(is)%stm * vpa(iv) * apar 

            ! Gyroaverage.
            call gyro_average(field, ivmu, g0(:, :, :, :, ivmu))
        end do

        ! Add this term to the right-hand-side of the GKE. 
        call add_explicit_term(g0, neo_mirror_apar_1(1, :, :), gout)

        ! We now need to add the second term, proportional to the vpa derivative of apar. 
        ! Construct 2 v_th <A∥_k>.
        ! Iterate over the (mu,vpa,s) points. 
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)

            ! Calculate the apar field. 
            field = 2.0 * spec(is)%stm * apar

            ! Gyroaverage.
            call gyro_average(field, ivmu, g0(:, :, :, :, ivmu))
        end do

        ! Add this term to the right-hand-side of the GKE. 
        call add_explicit_term(g0, neo_mirror_apar_2(1, :, :), gout)

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
        use arrays, only: neo_mirror_apar_1, neo_mirror_apar_2
        use arrays, only: initialised_neo_mirror

        implicit none

        if (allocated(neo_mirror_apar_1)) deallocate (neo_mirror_apar_1)
        if (allocated(neo_mirror_apar_2)) deallocate (neo_mirror_apar_2)
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
