! ================================================================================================================================================================================= !
! ------------------------------------------------------------ Evolves the neoclassical grad-B and curvature drift terms. --------------------------------------------------------- !​
! ================================================================================================================================================================================= !
!
! This module evolves the neoclassical drift terms on the RHS of the GKE. There will be a magnetic drift (ω_d) term and a curvature drift (ω_κ) term. 
! There will be a contribution proportional to kx and another one proportional to ky. The dimensionless magnetic drift terms are given by:
!
! ================================================================================================================================================================================= !

module gk_neo_drifts

   implicit none

   ! Make routines available to other modules. 
   public :: initialised_neo_wdrifty
   public :: init_neo_wdrifty
   public :: finish_neo_wdrifty
   public :: advance_neo_wdrifty_explicit
   
   public :: initialised_neo_wdriftx
   public :: init_neo_wdriftx
   public :: finish_neo_wdriftx
   public :: advance_neo_wdriftx_explicit

   private

   ! Only initialise once.
   logical :: initialised_neo_wdrifty = .false.
   logical :: initialised_neo_wdriftx = .false.

contains

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------ Initialise the y-component of the drifts. -------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_neo_wdrifty
        ! Parallelisation.
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        ! Grids. 
        use grids_time, only: code_dt
        use grids_species, only: spec
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac, vpa, mu
        use grids_z, only: nzgrid
        use grids_kxky, only: nalpha

        ! Geometry. 
        use geometry, only: bmag, B_times_kappa_dot_grady, B_times_gradB_dot_grady

        ! Neoclassical. 
        use neoclassical_terms_neo, only: neo_vpa_fac 

        ! Arrays. 
        use arrays, only: neo_wdrifty, initialised_neo_wdrifty

        ! For switching the ydrift on and off.
        use parameters_physics, only: neoydriftknob

        implicit none

        integer :: iz, iv, is, imu, ivmu

        ! Only intialise once.
        if (initialised_neo_wdrifty) return
        initialised_neo_wdrifty = .true.

        ! Allocate neo_wdrifty = neo_wdrifty[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_wdrifty)) then
            allocate (neo_wdrifty(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_wdrifty = 0.0
        end if

        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
           
            do iz = -nzgrid, nzgrid            
                neo_wdrifty(:, iz, ivmu) = vpa(iv) * B_times_kappa_dot_grady(:, iz) + mu(imu) * B_times_gradB_dot_grady(:, iz) / vpa(iv) 

                neo_wdrifty(:, iz, ivmu) = neo_wdrifty(:, iz, ivmu) * neoydriftknob * code_dt * 0.5 &
                * neo_vpa_fac(iz, ivmu, 1) * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) / ( bmag(:, iz) ** 2 )
            end do
        end do

    end subroutine init_neo_wdrifty


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------ Initialise the x-component of the drifts. -------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_neo_wdriftx
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
        use geometry, only: bmag, B_times_kappa_dot_gradx, B_times_gradB_dot_gradx

        ! Neoclassical. 
        use neoclassical_terms_neo, only: neo_vpa_fac

        ! Arrays. 
        use arrays, only: neo_wdriftx, initialised_neo_wdriftx

        ! Parameters. 
        use parameters_physics, only: neoxdriftknob

        implicit none

        integer :: iz, iv, is, imu, ivmu

        ! Only intialise once.
        if (initialised_neo_wdriftx) return
        initialised_neo_wdriftx = .true.

        ! Allocate neo_wdriftx = neo_wdriftx[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(neo_wdriftx)) then
            allocate (neo_wdriftx(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); neo_wdriftx = 0.0
        end if

        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
           
            do iz = -nzgrid, nzgrid                          
                neo_wdriftx(:, iz, ivmu) = vpa(iv) * B_times_kappa_dot_gradx(:, iz) + mu(imu) * B_times_gradB_dot_gradx(:, iz) / vpa(iv)

                neo_wdriftx(:, iz, ivmu) = neo_wdriftx(:, iz, ivmu) * neoxdriftknob * code_dt * 0.5 & 
                * neo_vpa_fac(iz, ivmu, 1) * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is) / ( bmag(:, iz) ** 2 )
            end do
        end do

    end subroutine init_neo_wdriftx                  


! ================================================================================================================================================================================= !
! -------------------------------------------------------------- Advance the y component of the drifts explicitly. ---------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_wdrifty_explicit(phi, apar, bpar, gout)
        ! Parallelisation.
        use mp, only: proc0
        use parallelisation_layouts, only: vmu_lo
        use parallelisation_layouts, only: iv_idx, imu_idx, is_idx

        ! Constants
        use constants, only: zi

        ! Data arrays.
        use arrays, only: neo_wdrifty   

        ! Parameters.
        use parameters_physics, only: fphi, include_apar, include_bpar

        ! Grids. 
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx, aky 
        use grids_species, only: spec
        use grids_velocity, only: vpa, mu
      
        ! Calculations.
        use calculations_add_explicit_terms, only: add_explicit_term
        use calculations_gyro_averages, only: gyro_average, gyro_average_j1

        ! Time this routine.
        use timers, only: time_gke
        use job_manage, only: time_message

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout        
        complex, dimension(:, :, :, :, :), allocatable :: g0y

        ! Local variables.
        integer :: ivmu, iv, is, iky, imu
        complex, dimension(:, :, :, :), allocatable :: field

        ! ======================================================================================= ! 
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= !
        !                                                                                         ! 
        ! Here we define the temporary array, g0y. These will hold                                ! 
        ! <g0y> = ∂<Χ_k>/∂y = i * ky * <Χ_k>.                                                     !
        !                                                                                         ! 
        ! We then multiply g0y by neo_wdrifty and add the result to the RHS of the GKE:           !
        !                                                                                         ! 
        ! add_explicit_term(g0y, neo_wdrifty(1, :, :), gout)                                      !
        !                                                                                         !
        ! ======================================================================================= !
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= !

        ! Start timing the advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_wdrifty advance')

        ! Allocate temporary arrays.
        allocate(g0y(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        allocate(field(naky, nakx, -nzgrid:nzgrid, ntubes))
       
        ! Calculate and add the phi contribution.
        ! Iterate over the (mu,vpa,s) points.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
         
            ! Calculate phi factor. 
            field = fphi * phi
       
            do iky = 1, naky
                field(iky, :, :, :) = zi * aky(iky) * field(iky, :, :, :)
            end do
         
            ! Gyroaverage.
            call gyro_average(field, ivmu, g0y(:, :, :, :, ivmu))
        end do

        ! Add the term to the RHS of the GKE by multiplying by the appropriate coeffecient. 
        call add_explicit_term(g0y, neo_wdrifty(1, :, :), gout)


        ! Calculate and add the apar contribution if needed. 
        if (include_apar) then
            ! Iterate over the (mu,vpa,s) points.
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                is = is_idx(vmu_lo, ivmu)
                iv = iv_idx(vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)

                ! Calculate apar factor. 
                field = 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar

                do iky = 1, naky
                    field(iky, :, :, :) = zi * aky(iky) * field(iky, :, :, :)
                end do

                ! Gyroaverage.
                call gyro_average(field, ivmu, g0y(:, :, :, :, ivmu))
            end do

            ! Add the terms to the RHS of the GKE by multiplying by the appropriate coeffecient. 
            call add_explicit_term(g0y, neo_wdrifty(1, :, :), gout)
        end if

        
        ! Calculate and add the apar contribution if needed. 
        if (include_bpar) then
            ! Iterate over the (mu,vpa,s) points.
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                is = is_idx(vmu_lo, ivmu)
                iv = iv_idx(vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)

                ! Calculate bpar factor. 
                field = 4.0 * mu(imu) * spec(is)%tz * bpar

                do iky = 1, naky
                    field(iky, :, :, :) = zi * aky(iky) * field(iky, :, :, :)
                end do

                ! Gyroaverage.
                call gyro_average_j1(field, ivmu, g0y(:, :, :, :, ivmu))
            end do

            ! Add the terms to the RHS of the GKE by multiplying by the appropriate coeffecient. 
            call add_explicit_term(g0y, neo_wdrifty(1, :, :), gout)
        end if

        ! Deallocate temporary arrays.
        deallocate(g0y)
        deallocate(field)

        ! Stop timing the advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_wdrifty advance')

    end subroutine advance_neo_wdrifty_explicit


! ================================================================================================================================================================================= !
! -------------------------------------------------------------- Advance the x component of the drifts explicitly. ---------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_neo_wdriftx_explicit(phi, apar, bpar, gout)
        ! Parallelisation.
        use mp, only: proc0
        use parallelisation_layouts, only: vmu_lo
        use parallelisation_layouts, only: iv_idx, imu_idx, is_idx

        ! Constants
        use constants, only: zi

        ! Data arrays.
        use arrays, only: neo_wdriftx

        ! Parameters.
        use parameters_physics, only: fphi, include_apar, include_bpar

        ! Grids. 
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx, akx 
        use grids_species, only: spec
        use grids_velocity, only: vpa, mu
      
        ! Calculations.
        use calculations_add_explicit_terms, only: add_explicit_term
        use calculations_gyro_averages, only: gyro_average, gyro_average_j1

        ! Time this routine.
        use timers, only: time_gke
        use job_manage, only: time_message

        implicit none

        complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
        complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout        

        ! Local variables.
        integer :: ivmu, iv, is, ikx, imu
        complex, dimension(:, :, :, :, :), allocatable :: g0x
        complex, dimension(:, :, :, :), allocatable :: field

        ! ======================================================================================= ! 
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= !
        !                                                                                         ! 
        ! Here we define the temporary array, g0x. This will hold                                 ! 
        ! <g0x> = ∂<Χ_k>/∂x = i * kx * <Χ_k>.                                                     ! 
        !                                                                                         ! 
        ! We then multiply g0x by neo_wdriftx and add the result to the RHS of the GKE:           !
        !                                                                                         ! 
        ! add_explicit_term(g0x, neo_wdriftx(1, :, :), gout)                                      !
        !                                                                                         !
        ! ======================================================================================= !
        ! --------------------------------------------------------------------------------------- !
        ! ======================================================================================= !

        ! Start timing the advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_wdriftx advance')

        ! Allocate temporary arrays.
        allocate (g0x(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
        allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes))
       
        ! Iterate over the (mu,vpa,s) points.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
         
            ! Calculate phi factor. 
            field = fphi * phi
       
            do ikx = 1, nakx
                field(:, ikx, :, :) = zi * akx(ikx) * field(:, ikx, :, :)
            end do
         
            ! Gyroaverage.
            call gyro_average(field, ivmu, g0x(:, :, :, :, ivmu))
        end do

        ! Add the terms to the RHS of the GKE. 
        call add_explicit_term(g0x, neo_wdriftx(1, :, :), gout)


        ! Calculate and add the apar contribution if needed. 
        if (include_apar) then
            ! Iterate over the (mu,vpa,s) points.
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                is = is_idx(vmu_lo, ivmu)
                iv = iv_idx(vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)

                ! Calculate apar factor. 
                field = 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar

                do ikx = 1, nakx
                    field(:, ikx, :, :) = zi * akx(ikx) * field(:, ikx, :, :)
                end do

                ! Gyroaverage.
                call gyro_average(field, ivmu, g0x(:, :, :, :, ivmu))
            end do

            ! Add the terms to the RHS of the GKE. 
            call add_explicit_term(g0x, neo_wdriftx(1, :, :), gout)
        end if

        
        ! Add bpar contribution. 
        if (include_bpar) then
            ! Iterate over the (mu,vpa,s) points.
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                is = is_idx(vmu_lo, ivmu)
                iv = iv_idx(vmu_lo, ivmu)
                imu = imu_idx(vmu_lo, ivmu)

                ! Calculate bpar factor. 
                field = 4.0 * mu(imu) * spec(is)%tz * bpar

                do ikx = 1, nakx
                    field(:, ikx, :, :) = zi * akx(ikx) * field(:, ikx, :, :)
                end do

                ! Gyroaverage.
                call gyro_average_j1(field, ivmu, g0x(:, :, :, :, ivmu))
            end do

            ! Add the terms to the RHS of the GKE. 
            call add_explicit_term(g0x, neo_wdriftx(1, :, :), gout)
        end if

        ! Deallocate temporary arrays.
        deallocate (g0x)
        deallocate (field)

        ! Stop timing the advance.
        if (proc0) call time_message(.false., time_gke(:, 6), 'neo_wdriftx advance')

    end subroutine advance_neo_wdriftx_explicit


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------- Finish the y component of the drifts. -------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_neo_wdrifty
        use arrays, only: neo_wdrifty, initialised_neo_wdrifty

        implicit none

        if (allocated(neo_wdrifty)) deallocate (neo_wdrifty)
        initialised_neo_wdrifty = .false.

    end subroutine finish_neo_wdrifty


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------- Finish the x component of the drifts. -------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine finish_neo_wdriftx
        use arrays, only: neo_wdriftx, initialised_neo_wdriftx

        implicit none

        if (allocated(neo_wdriftx)) deallocate (neo_wdriftx)
        initialised_neo_wdriftx = .false.

    end subroutine finish_neo_wdriftx


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_drifts
