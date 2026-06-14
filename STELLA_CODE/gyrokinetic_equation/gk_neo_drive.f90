! ================================================================================================================================================================================= !
! ------------------------------------------------------------------- Evolves the neoclassical gradient drive. -------------------------------------------------------------------- !​
! ================================================================================================================================================================================= !
! 
! This module evolves the HO neoclassical gradient drive:
!
! (k⟂ x B₀) ⋅∇F₁ = ∂F₁/∂ψ (k⟂ x B₀) ⋅∇ψ + ∂F₁/∂θ (k⟂ x B₀) ⋅∇θ 
!
! In the conventional theory this has only a ∂F₀/∂ψ component, which is captured by wstar defined in gk_drive.f90. When F₁ is included, the gradient drive has both a ∂F₁/∂ψ  
! component and a ∂F₁/∂θ component. When using a miller geometry, this is equivalent to ∂F₁/∂z. In this routine we define two coeffecients, wstar1 and wpol. wstar1 is the 
! higher-order counterpart to wstar. The dimensionless expression is  
! 
! wstar1 = wstar1psi + wstar1z
!
! where 
!
! wstar1psi =
!
! and 
!
! wstar1z = 
!
! wstar1 must then be multiplied by i * ky * <Χ_k> and then added to the RHS of the GKE. Similarly wpol is given by: 
! 
! wpol = 
!
! This must be multiplied by i * kx * <Χ_k> and added to the RHS of the GKE.
!
! ================================================================================================================================================================================= !

module gk_neo_drive

    ! Load debug flags.
    ! use debug_flags, only: debug => neo_drive_debug
   
    implicit none

    ! Make routines available to other modules. 
    public :: initialised_wstar1y, initialised_wstar1x
    public :: init_wstar1y, init_wstar1x
    public :: finish_wstar1y, finish_wstar1x 
    public :: advance_wstar1y_explicit, advance_wstar1x_explicit 

    private
   
    ! Only initialise once.
    logical :: initialised_wstar1y = .false.
    logical :: initialised_wstar1x = .false.

contains

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------ Initialise wstar1y. ------------------------------------------------------------------------------ ! 
! ================================================================================================================================================================================= !

    subroutine init_wstar1y
        ! Parallelisation.
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      
        ! Grids.
        use grids_time, only: code_dt
        use grids_kxky, only: nalpha
        use grids_z, only: nzgrid, zed
        use grids_species, only: spec
        use grids_velocity, only: vperp2, vpa, mu
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
        ! Geometry. 
        use geometry, only: dydalpha, drhodpsi, clebsch_factor, dxdpsi 
        use geometry, only: bmag, gradx_dot_grady, B_times_gradB_dot_grady
        use geometry, only: geo_surf, Rmajor

        ! NEO data.
        use neoclassical_terms_neo, only: neo_vpa_fac
        use neoclassical_terms_neo, only: neo_h, neo_phi
        use neoclassical_terms_neo, only: dneo_h_dpsi, dneo_phi_dpsi   
        use neoclassical_terms_neo, only: dneo_h_dz, dneo_phi_dz

        ! Arrays. 
        use arrays, only: wstar1y, initialised_wstar1y

        implicit none

        ! Indices.
        integer :: is, imu, iv, ivmu, iz
     
        ! wstar1y has a component which is proportional to ∂F_1/∂ψ, we will call this wstar1ypsi. 
        ! Similarly the component proportional to ∂F_1/∂z will be called wstar1yz. 
        ! The component proportional to the parallel velocity derivative is called wstar1yvpa. 
        ! Splitting up the calculation this way should make the maths easier to follow.
        real, dimension(:, :, :), allocatable :: wstar1ypsi, wstar1yz, wstar1yvpa         

        ! Only intialise omega_{*,k,s,1,y} once.
        if (initialised_wstar1y) return
        initialised_wstar1y = .true.

        ! Allocate omega_{*,k,s,1,y} = wstar1y[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(wstar1y)) then
            allocate (wstar1y(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstar1y = 0.0
        end if

        ! Allocate the temporary arrays. 
        allocate (wstar1ypsi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstar1ypsi = 0.0      
        allocate (wstar1yz(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstar1yz = 0.0
        allocate (wstar1yvpa(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstar1yvpa = 0.0

        ! First calculate the component proportional to the psi derivative of F_1. 
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid
                wstar1ypsi(:, iz, ivmu) = dydalpha * drhodpsi * ( dneo_h_dpsi(iz, ivmu, 1) - spec(is)%z * dneo_phi_dpsi(iz) ) / clebsch_factor
            end do
        end do

        ! Now calculate the component proportional to the z derivative of F_1. 
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid 
                wstar1yz(:, iz, ivmu) = - dydalpha * drhodpsi * geo_surf%shat * zed(iz) / ( geo_surf%rhoc * clebsch_factor ) &
                - clebsch_factor * gradx_dot_grady(:, iz) / ( Rmajor(iz) * Rmajor(iz) * bmag(:, iz) * bmag(:, iz) * geo_surf%qinp * dxdpsi )
               
                ! Multiply by the F_1 factor.
                wstar1yz(:, iz, ivmu) = wstar1yz(:, iz, ivmu) * ( dneo_h_dz(iz, ivmu, 1) - spec(is)%z * dneo_phi_dz(iz) )
            end do
        end do

        ! Finally add the remaining corrections to wstar1y, namely the one proportinal to the vpa derivative of F₁.
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid
                wstar1yvpa(:, iz, ivmu) = - mu(imu) * B_times_gradB_dot_grady(:, iz)  / ( vpa(iv) * bmag(:, iz) * bmag(:, iz) ) 
 
                ! Multiply by the F_1 factor.
                wstar1yvpa(:, iz, ivmu) = wstar1yvpa(:, iz, ivmu) * ( neo_vpa_fac(iz, ivmu, 1) + 2.0 * vpa(iv) * ( neo_h(iz, ivmu, 1) - spec(is)%z * neo_phi(iz) ) ) 
            end do
        end do

        ! Finally, calculate wstar1y = wstar1ypsi + wstar1yz + wstar1yvpa, including the common factor.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
          
            do iz = -nzgrid, nzgrid
                wstar1y(:, iz, ivmu) = 0.5 * code_dt * ( wstar1ypsi(:, iz, ivmu) + wstar1yz(:, iz, ivmu) + wstar1yvpa(:, iz, ivmu) ) &
                * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is)
            end do
        end do

        ! Deallocate temporary arrays. 
        deallocate(wstar1ypsi)
        deallocate(wstar1yz)
        deallocate(wstar1yvpa)

    end subroutine init_wstar1y


! ================================================================================================================================================================================= !
! ----------------------------------------------------------------------------- Initialise wstar1x. ------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_wstar1x
        ! Parallelisation.
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      
        ! Grids.
        use grids_time, only: code_dt
        use grids_kxky, only: nalpha
        use grids_z, only: nzgrid, zed
        use grids_species, only: spec
        use grids_velocity, only: vpa, mu 
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
        ! Geometry.
        use geometry, only: dydalpha, dxdpsi, drhodpsi, clebsch_factor
        use geometry, only: bmag, b_dot_gradz, gradx_dot_gradx, B_times_gradB_dot_gradx
        use geometry, only: geo_surf, Rmajor
   
        ! NEO data.
        use neoclassical_terms_neo, only: neo_vpa_fac
        use neoclassical_terms_neo, only: dneo_h_dz, dneo_phi_dz
        use neoclassical_terms_neo, only: neo_h, neo_phi

        ! Arrays. 
        use arrays, only: wstar1x, initialised_wstar1x

        implicit none

        ! Local variables.
        integer :: is, imu, iv, ivmu, iz
        real, dimension(:, :, :), allocatable :: wstar1xz, wstar1xvpa         

        ! Only intialise once.
        if (initialised_wstar1x) return
        initialised_wstar1x = .true.

        ! Allocate temporary arrays. 
        allocate (wstar1xz(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstar1xz = 0.0
        allocate (wstar1xvpa(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstar1xvpa = 0.0

        ! Allocate wstar1x = wstar1x[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(wstar1x)) then
            allocate (wstar1x(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstar1x = 0.0
        end if
  
        ! First calculate the component proportional to the z derivative of F_1. 
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid
                wstar1xz(:, iz, ivmu) = ( (dxdpsi / clebsch_factor) - clebsch_factor * gradx_dot_gradx(:, iz) / (Rmajor(iz) * Rmajor(iz) * bmag(:, iz) * bmag(:, iz) * dxdpsi ) ) &
                / geo_surf%qinp

                ! Multiply by the F_1 factor.
                wstar1xz(:, iz, ivmu) = wstar1xz(:, iz, ivmu) * ( dneo_h_dz(iz, ivmu, 1) - spec(is)%z * dneo_phi_dz(iz) ) 
            end do  
        end do

        ! Now calculate the component proportional to the vpa derivative of F_1.
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid
                wstar1xvpa(:, iz, ivmu) = - mu(imu) * B_times_gradB_dot_gradx(:, iz)  / ( vpa(iv) * bmag(:, iz) * bmag(:, iz) )

                ! Multiply by the F_1 factor. 
                wstar1xvpa(:, iz, ivmu) = wstar1xvpa(:, iz, ivmu) * ( neo_vpa_fac(iz, ivmu, 1) + 2.0 * vpa(iv) * ( neo_h(iz, ivmu, 1) - spec(is)%z * neo_phi(iz) ) )
            end do
        end do

        ! Finally, calculate wstar1x = wstar1xz + wstar1xvpa, including the common factor.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            do iz = -nzgrid, nzgrid
                wstar1x(:, iz, ivmu) = 0.5 * code_dt * ( wstar1xz(:, iz, ivmu) + wstar1xvpa(:, iz, ivmu) ) &
                * maxwell_vpa(iv, is) * maxwell_mu(:, iz, imu, is) * maxwell_fac(is)
            end do
        end do

        ! Deallocate temporary arrays. 
        deallocate (wstar1xz)
        deallocate (wstar1xvpa)

    end subroutine init_wstar1x               


! ================================================================================================================================================================================= !
! -------------------------------------------------------------------------- Advance wstar1y explicitly. -------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_wstar1y_explicit(phi, gout)
        ! Parallelisation.
        use mp, only: proc0
      
        ! Data arrays.
        use arrays, only: wstar1y
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
        ! Add the k_y drive term to the GKE:                                                         !
        !                                                                                            ! 
        ! wstar1y * i * ky * <Χ_k>                                                                   !
        !                                                                                            !
        ! First we calculate i * ky * <Χ_k> which corresponds to ∂<Χ_k>/∂y:                          !
        !                                                                                            ! 
        ! Then multiply with <wstar1y> and add it to the RHS of the GKE:                             !
        !                                                                                            !
        ! add_explicit_term(g0, wstar1y(1, :, :), gout)                                              !
        !                                                                                            ! 
        ! ========================================================================================== !
        ! ------------------------------------------------------------------------------------------ !
        ! ========================================================================================== !

        ! Start timing the time advance due to the driving gradient.
        if (proc0) call time_message(.false., time_gke(:, 6), ' wstar1y advance')

        ! Allocate temporary array for <g0> = ∂Χ_k/∂y = i *ky * <Χ_k>.
        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      
        ! Calculate <g0>.
        ! if (debug) write (*, *) 'time_advance::solve_gke::get_dchidy'
        call get_dchidy(phi, apar, bpar, g0)

        ! Add the drive term to the RHS of the gyrokinetic equation. 
        ! if (debug) write (*, *) 'time_advance::solve_gke::add_wstar1y_term'
        call add_explicit_term(g0, wstar1y(1, :, :), gout)

        ! Deallocate <g0>.
        deallocate (g0)

        ! Stop timing the time advance due to the driving gradient.
        if (proc0) call time_message(.false., time_gke(:, 6), ' wstar1y advance')
 
    end subroutine advance_wstar1y_explicit


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------- Advance wstar1x explicitly. --------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_wstar1x_explicit(phi, gout)
        ! Parallelisation.
        use mp, only: proc0
      
        ! Data arrays.
        use arrays, only: wstar1x
        use arrays_fields, only: apar, bpar
      
        ! Grids.
        use parallelisation_layouts, only: vmu_lo
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx
      
        ! Calculations.
        use calculations_add_explicit_terms, only: add_explicit_term
        use calculations_kxky_derivatives, only: get_dchidx
      
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
        ! Add the k_x drive term to the GKE:                                                         !
        !                                                                                            ! 
        ! wstar1x * i * kx * <Χ_k>                                                                   !
        !                                                                                            !
        ! First we calculate i * kx * <Χ_k> which corresponds to d<Χ_k>/dx:                          !
        !                                                                                            ! 
        ! Then multiply with wstar1x and add it to the right-hand-side of the GKE:                   !
        !                                                                                            !
        ! add_explicit_term(g0, wstar1x(1, :, :), gout)                                              !
        !                                                                                            ! 
        ! ========================================================================================== !
        ! ------------------------------------------------------------------------------------------ !
        ! ========================================================================================== !

        ! Start timing the time advance due to the driving gradient.
        if (proc0) call time_message(.false., time_gke(:, 6), ' wstar1x advance')

        ! Allocate temporary array for <g0> = d<Χ_k>/dx  = i * kx * <Χ_k>.
        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      
        ! Calculate <g0>.
        ! if (debug) write (*, *) 'time_advance::solve_gke::get_dchidx'
        call get_dchidx(phi, apar, bpar, g0)

        ! Add the drive term to the RHS of the GKE. 
        ! if (debug) write (*, *) 'time_advance::solve_gke::add_wstar1x_term'
        call add_explicit_term(g0, wstar1x(1, :, :), gout)

        ! Deallocate <g0>.
        deallocate (g0)

        ! Stop timing the time advance due to the driving gradient.
        if (proc0) call time_message(.false., time_gke(:, 6), ' wstar1x advance')

    end subroutine advance_wstar1x_explicit


! ================================================================================================================================================================================= !
! -------------------------------------------------------------------------------- Finish wstar1y. -------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_wstar1y
        use arrays, only: wstar1y

        implicit none

        if (allocated(wstar1y)) deallocate (wstar1y)
        initialised_wstar1y = .false.

    end subroutine finish_wstar1y


! ================================================================================================================================================================================= !
! -------------------------------------------------------------------------------- Finish wstar1x. -------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine finish_wstar1x   
        use arrays, only: wstar1x       

        implicit none

        if (allocated(wstar1x)) deallocate (wstar1x)
        initialised_wstar1x = .false.

    end subroutine finish_wstar1x   


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

end module gk_neo_drive
