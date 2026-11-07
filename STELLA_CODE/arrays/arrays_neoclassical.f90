! ================================================================================================================================================================================= !
! ------------------------------------------------ Arrays holding all neoclassical quantities needed for HO GK simulations in stella. ----------------------------------------------- !
! ================================================================================================================================================================================= !
!
! This module owns every neoclassical array used in the GK-Maxwell equations. 
! It contains no physics and no knowledge of where the data comes from; a drift kinetic provider
! (either NEO or SFINCS) is responsible for populating these arrays after they have been allocated here.
!
! Arrays defined here are:
!
!   - Local distributions and their derivatives (not distributed over ranks), mostly for use in the 
!   field equations and some distribution transforming routines: 
!   (nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec)
!
!   - Velocity-distributed distributions, mostly for use in the GKE:
!   (nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc)
!  
!   - ϕ^1_0 and its derivatives : (nalpha, -nzgrid:nzgrid)
!
! Usage:
!   call allocate_neoclassical_arrays    ! Before the provider populates the arrays.
!   call deallocate_neoclassical_arrays  ! At finish.
!
! ================================================================================================================================================================================= !

module arrays_neoclassical
    implicit none

    public :: neoclassical_is_enabled
    public :: include_neoclassical_terms

    public :: allocate_neoclassical_arrays
    public :: deallocate_neoclassical_arrays
    public :: initialised_neoclassical_arrays

    public :: neo_h, dneo_h_dpsi, dneo_h_dz          
    public :: neo_phi, dneo_phi_dpsi, dneo_phi_dz    

    public :: dneo_h_dmu, dneo_h_dvpa
    public :: dneo_h_dmu_global, dneo_h_dvpa_global

    public :: neo_vpa_fac, neo_mu_fac                
    public :: neo_vpa_fac_global, neo_mu_fac_global  
    
    public :: neo_zed_fac
    
    private

    ! (nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec).
    real, dimension(:, :, :, :, :), allocatable :: dneo_h_dvpa_global, dneo_h_dmu_global
    real, dimension(:, :, :, :, :), allocatable :: neo_vpa_fac_global, neo_mu_fac_global

    ! (nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc).
    real, dimension(:, :, :), allocatable :: neo_h, dneo_h_dpsi, dneo_h_dz
    real, dimension(:, :, :), allocatable :: dneo_h_dmu, dneo_h_dvpa
    real, dimension(:, :, :), allocatable :: neo_vpa_fac, neo_mu_fac
    real, dimension(:, :, :), allocatable :: neo_zed_fac

    ! (nalpha, -nzgrid:nzgrid).
    real, dimension(:, :), allocatable :: neo_phi, dneo_phi_dpsi, dneo_phi_dz

    logical :: include_neoclassical_terms = .false.
    logical :: initialised_neoclassical_arrays = .false.

contains

! ================================================================================================================================================================================= !
! ------------------------------------------------------------------ Allocate all neoclassical arrays. ---------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    logical function neoclassical_is_enabled()
        neoclassical_is_enabled = include_neoclassical_terms
    end function neoclassical_is_enabled


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------ Allocate all neoclassical arrays. ---------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine allocate_neoclassical_arrays
        ! Parallelisation.
        use parallelisation_layouts, only: vmu_lo

        ! Grids.
        use grids_z, only: nzgrid
        use grids_kxky, only: nalpha
        use grids_velocity, only: nvpa, nmu
        use grids_species, only: nspec

        implicit none

        if (initialised_neoclassical_arrays) return
        initialised_neoclassical_arrays = .true.

        ! Global arrays: (nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec).
        if (.not. allocated(dneo_h_dmu_global)) allocate(dneo_h_dmu_global(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec)); dneo_h_dmu_global = 0.0
        if (.not. allocated(dneo_h_dvpa_global)) allocate(dneo_h_dvpa_global(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec)); dneo_h_dvpa_global = 0.0
        if (.not. allocated(neo_vpa_fac_global)) allocate(neo_vpa_fac_global(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec)); neo_vpa_fac_global = 0.0
        if (.not. allocated(neo_mu_fac_global)) allocate(neo_mu_fac_global(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec)); neo_mu_fac_global = 0.0

        ! Velocity-distributed arrays: (nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc).
        if (.not. allocated(neo_h)) allocate(neo_h(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc)); neo_h = 0.0 
        if (.not. allocated(dneo_h_dpsi)) allocate(dneo_h_dpsi(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc)); dneo_h_dpsi = 0.0 
        if (.not. allocated(dneo_h_dz)) allocate(dneo_h_dz(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc)); dneo_h_dz = 0.0
        if (.not. allocated(dneo_h_dmu)) allocate(dneo_h_dmu(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc)); dneo_h_dmu = 0.0
        if (.not. allocated(dneo_h_dvpa)) allocate(dneo_h_dvpa(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc)); dneo_h_dvpa = 0.0
        if (.not. allocated(neo_vpa_fac)) allocate(neo_vpa_fac(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc)); neo_vpa_fac = 0.0
        if (.not. allocated(neo_mu_fac)) allocate(neo_mu_fac(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc)); neo_mu_fac = 0.0
        if (.not. allocated(neo_zed_fac)) allocate(neo_zed_fac(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc)); neo_zed_fac = 0.0

        ! ϕ^1_0 arrays: (nalpha, -nzgrid:nzgrid).
        if (.not. allocated(neo_phi)) allocate(neo_phi(nalpha, -nzgrid:nzgrid)); neo_phi = 0.0
        if (.not. allocated(dneo_phi_dpsi)) allocate(dneo_phi_dpsi(nalpha, -nzgrid:nzgrid)); dneo_phi_dpsi = 0.0
        if (.not. allocated(dneo_phi_dz)) allocate(dneo_phi_dz(nalpha, -nzgrid:nzgrid)); dneo_phi_dz = 0.0

    end subroutine allocate_neoclassical_arrays


! ================================================================================================================================================================================= !
! ----------------------------------------------------------------- Deallocate all neoclassical arrays. --------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine deallocate_neoclassical_arrays
        implicit none

        if (allocated(neo_h)) deallocate(neo_h)
        if (allocated(dneo_h_dpsi)) deallocate(dneo_h_dpsi)
        if (allocated(dneo_h_dz)) deallocate(dneo_h_dz)
        if (allocated(dneo_h_dmu)) deallocate(dneo_h_dmu)
        if (allocated(dneo_h_dvpa)) deallocate(dneo_h_dvpa)
        if (allocated(neo_vpa_fac)) deallocate(neo_vpa_fac)
        if (allocated(neo_mu_fac)) deallocate(neo_mu_fac)
        if (allocated(neo_zed_fac)) deallocate(neo_zed_fac)

        if (allocated(neo_phi)) deallocate(neo_phi)
        if (allocated(dneo_phi_dpsi)) deallocate(dneo_phi_dpsi)
        if (allocated(dneo_phi_dz)) deallocate(dneo_phi_dz)

        if (allocated(dneo_h_dmu_global)) deallocate(dneo_h_dmu_global)
        if (allocated(dneo_h_dvpa_global)) deallocate(dneo_h_dvpa_global)
        if (allocated(neo_vpa_fac_global)) deallocate(neo_vpa_fac_global)
        if (allocated(neo_mu_fac_global)) deallocate(neo_mu_fac_global)

        initialised_neoclassical_arrays = .false.

    end subroutine deallocate_neoclassical_arrays


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------ End Module. -------------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

end module arrays_neoclassical