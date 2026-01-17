! ================================================================================================================================================================================= !
! ------------------------------------------------------------------- Evolves the neoclassical gradient drive. -------------------------------------------------------------------- !​
! ================================================================================================================================================================================= !
! 
! This module evolves the neoclassical gradient drive:
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
! wstar1psi = - (1/2C) * exp(-v²) * dy/dα * dρ/dψ * ( ∂F₁/∂ρ + F₁ * dlnF₀/dρ )
!
! and 
!
! wstar1z =  (1/2C) * exp(-v²) * ( dψ/dx * ( |∇x ⋅∇y|/( q * R₀² * B₀² ) ) + (ŝ/ρ) * dρ/dψ * dy/dα * z ) * ∂F₁/∂z
!
! wstar1 must then be multiplied by i * ky * <Χ_k> and then added to the RHS of the GKE. Similarly wpol is given by: 
! 
! wpol = (1/2Cq) * exp(-v²) * ( dψ/dx * ( |∇x ⋅∇x|/( R₀² * B₀² ) ) - dx/dψ ) * ∂F₁/∂z
!
! This must be multiplied by i * kx * <Χ_k> and added to the RHS of the GKE.
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
        ! Parallelisation.
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      
        ! Grids.
        use grids_time, only: code_dt
        use grids_kxky, only: nalpha
        use grids_z, only: nzgrid, zed
        use grids_species, only: spec
        use grids_velocity, only: vperp2, vpa
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
        ! Geometry. 
        use geometry, only: dydalpha, drhodpsi, clebsch_factor, dxdpsi
        use geometry, only: bmag, gradx_dot_grady
        use geometry, only: geo_surf

        ! NEO data.
        use neoclassical_terms_neo, only: neo_h, neo_phi             
        use neoclassical_terms_neo, only: dneo_h_dpsi, dneo_phi_dpsi   
        use neoclassical_terms_neo, only: dneo_h_dz, dneo_phi_dz

        ! Arrays. 
        use arrays, only: wstar1, initialised_wstar1

        ! Rescale the drive term with <wstar1knob>.
        ! use parameters_physics, only: wstarknob

        implicit none

        ! Indices.
        integer :: is, imu, iv, ivmu, iz
     
        ! wstar1 has a component which is proportional to ∂F_1/∂ψ, we will call this wstar1psi. 
        ! Similarly the component proportional to ∂F_1/∂z will be called wstarz. 
        ! We need to declare temporary arrays for both of these; this should make the maths easier to follow.
        real, dimension(:, :, :), allocatable :: wstar1psi, wstar1z         
        
        ! To make the calculations easier to follow, we also calculate energy = v_parallel² + 2 mu B for each velocity point.
        real, dimension(:, :), allocatable :: energy

        ! Only intialise omega_{*,k,s,1} once.
        if (initialised_wstar1) return
        initialised_wstar1 = .true.

        ! Allocate omega_{*,k,s,1} = wstar1[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(wstar1)) then
            allocate (wstar1(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wstar1 = 0.0
        end if

        ! Allocate the temporary arrays. 

        allocate (energy(nalpha, -nzgrid:nzgrid))        
        allocate (wstar1psi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))      
        allocate (wstar1z(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)

            ! Calculate <energy>[ialpha,iz] = v_parallel² + 2 mu B = vpa(iv)**2 + vperp2(ialpha, iz, imu). 
            energy(:, :) = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
 
            ! Mutliply by the magnetic geometry prefactor.
            wstar1psi(:, :, ivmu) = - (0.5/clebsch_factor) * dydalpha * drhodpsi 

            ! Mutliply by the species dependent prefactor.
            wstar1psi(:, :, ivmu) = wstar1psi(:, :, ivmu) * maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)

            do iz = -nzgrid, nzgrid
                ! Multiply by the neolcassical distribution factor. 
                wstar1psi(:, iz, ivmu) = wstar1psi(:, iz, ivmu) * ( dneo_h_dpsi(iz, ivmu, 1) - spec(is)%z * dneo_phi_dpsi(iz, 1) &
                + ( neo_h(iz, ivmu, 1) - spec(is)%z * neo_phi(iz, 1)) * ( spec(is)%fprim + spec(is)%tprim * ( energy(:, iz) - 1.5 ) ) )
            end do
      end do

      ! The energy array is no longer needed. 
      deallocate(energy)

      ! We must add the remaining corrections to wstar1, namely those proportinal to the z derivative of F₁.
      ! Iterate over velocity space.
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          is = is_idx(vmu_lo, ivmu)
          imu = imu_idx(vmu_lo, ivmu)
          iv = iv_idx(vmu_lo, ivmu)

          ! Mutliply by the species dependent prefactor.
          wstar1z(:, :, ivmu) = maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)

          do iz = -nzgrid, nzgrid
              ! Multiply by the neolcassical distribution factor. 
              wstar1z(:, iz, ivmu) = wstar1z(:, iz, ivmu) * ( dneo_h_dz(iz, ivmu, 1) - spec(is)%z * dneo_phi_dz(iz, 1) )

              ! Mutliply by the magnetic geometry prefactor.
              wstar1z(:, iz, ivmu) = wstar1z(:, iz, ivmu) * (0.5/clebsch_factor) * ( gradx_dot_grady(:, iz)/( dxdpsi * geo_surf%qinp &
              * ( bmag(:, iz) * geo_surf%rmaj )**2 ) + (geo_surf%shat/geo_surf%rhoc) * drhodpsi * dydalpha * zed(iz) )
          end do
      end do

      ! Finally, calculate wstar1 = wstar1psi + wstar1z.
      wstar1 = ( wstar1psi + wstar1z ) * code_dt

      ! Deallocate the remaining temporary arrays. 
      deallocate(wstar1psi)
      deallocate(wstar1z)

    end subroutine init_wstar1


! ================================================================================================================================================================================= !
! -------------------------------------------------------------------------------- Initialise wpol. ------------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine init_wpol
        ! Parallelisation.
        use mp, only: mp_abort
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      
        ! Grids.
        use grids_time, only: code_dt
        use grids_kxky, only: nalpha
        use grids_z, only: nzgrid, zed
        use grids_species, only: spec
        use grids_velocity, only: vpa
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
        ! Geometry.
        use geometry, only: dydalpha, drhodpsi, clebsch_factor, dxdpsi
        use geometry, only: bmag, gradx_dot_gradx
        use geometry, only: geo_surf
   
        ! NEO data.
        use neoclassical_terms_neo, only: dneo_h_dz, dneo_phi_dz

        ! Arrays. 
        use arrays, only: wpol, initialised_wpol

        ! Rescale the drive term with <wpolknob>.
        ! use parameters_physics, only: wpolknob

        implicit none

        ! Indices.
        integer :: is, imu, iv, ivmu, iz
         
        ! Only intialise omega_pol once.
        if (initialised_wpol) return
        initialised_wpol = .true.

        ! Allocate omega_pol = wpol[ialpha, iz, i[mu,vpa,s]].
        if (.not. allocated(wpol)) then
            allocate (wpol(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); wpol = 0.0
        end if
  
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
          
            ! Multiply by the species dependent factor.
            wpol(:, :, ivmu) = maxwell_vpa(iv, is) * maxwell_mu(:, :, imu, is) * maxwell_fac(is)

            ! Multiply by the magnetic geometry factor. 
            wpol(:, :, ivmu) = ( 0.5 / ( clebsch_factor * geo_surf%qinp ) ) * ( gradx_dot_gradx(:, :)/( dxdpsi * ( bmag(:, :) * geo_surf%rmaj )**2 )  - dxdpsi )

            ! Multiply by the neoclassical coeffecient. 
            do iz = -nzgrid, nzgrid
                wpol(:, iz, ivmu) = wpol(:, iz, ivmu) * (dneo_h_dz(iz, ivmu, 1) - spec(is)%z * dneo_phi_dz(iz, 1))
            end do  

            ! Finally multipy by code_dt.
            wpol(:, :, ivmu) =  wpol(:, :, ivmu) * code_dt
        end do

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
        ! Add the neoclassical ψ drive term to the GKE:                                              !
        !                                                                                            ! 
        ! wstar1 * i * ky * <Χ_k>                                                                    !
        !                                                                                            !
        ! First we calculate i * ky * <Χ_k> which corresponds to ∂<Χ_k>/∂y:                          !
        !                                                                                            ! 
        ! Then multiply with <wstar1> and add it to the RHS of the GKE:                              !
        !                                                                                            !
        ! add_explicit_term(g0, wstar1(1, :, :), gout)                                               !
        !                                                                                            ! 
        ! ========================================================================================== !
        ! ------------------------------------------------------------------------------------------ !
        ! ========================================================================================== !

        ! Start timing the time advance due to the driving gradient.
        if (proc0) call time_message(.false., time_gke(:, 6), ' wstar1 advance')

        ! Allocate temporary array for <g0> = ∂Χ_k/∂y = i *ky * <Χ_k>.
        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      
        ! Calculate <g0>.
        ! if (debug) write (*, *) 'time_advance::solve_gke::get_dchidy'
        call get_dchidy(phi, apar, bpar, g0)

        ! Add the drive term to the RHS of the gyrokinetic equation. 
        ! if (debug) write (*, *) 'time_advance::solve_gke::add_wstar1_term'
        call add_explicit_term(g0, wstar1(1, :, :), gout)

        ! Deallocate <g0>.
        deallocate (g0)

        ! Stop timing the time advance due to the driving gradient.
        if (proc0) call time_message(.false., time_gke(:, 6), ' wstar1 advance')
 
    end subroutine advance_wstar1_explicit


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------- Advance wpol explicitly. --------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================= !

    subroutine advance_wpol_explicit(phi, gout)
        ! Parallelisation.
        use mp, only: proc0
      
        ! Data arrays.
        use arrays, only: wpol
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
        ! Add the neoclassical z drive term to the GKE:                                              !
        !                                                                                            ! 
        ! wpol * i * kx * <Χ_k>                                                                      !
        !                                                                                            !
        ! First we calculate i * kx * <Χ_k> which corresponds to d<Χ_k>/dx:                          !
        !                                                                                            ! 
        ! Then multiply with wpol and add it to the right-hand-side of the GKE:                      !
        !                                                                                            !
        ! add_explicit_term(g0, wpol(1, :, :), gout)                                                 !
        !                                                                                            ! 
        ! ========================================================================================== !
        ! ------------------------------------------------------------------------------------------ !
        ! ========================================================================================== !

        ! Start timing the time advance due to the driving gradient.
        if (proc0) call time_message(.false., time_gke(:, 6), ' wpol advance')

        ! Allocate temporary array for <g0> = d<Χ_k>/dx  = i * kx * <Χ_k>.
        allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      
        ! Calculate <g0>.
        ! if (debug) write (*, *) 'time_advance::solve_gke::get_dchidx'
        call get_dchidx(phi, apar, bpar, g0)

        ! Add the drive term to the RHS of the GKE. 
        ! if (debug) write (*, *) 'time_advance::solve_gke::add_wpol_term'
        call add_explicit_term(g0, wpol(1, :, :), gout)

        ! Deallocate <g0>.
        deallocate (g0)

        ! Stop timing the time advance due to the driving gradient.
        if (proc0) call time_message(.false., time_gke(:, 6), ' wpol advance')

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
