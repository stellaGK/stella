! ================================================================================================================================================================================= !
! ------------------------ Routines for remapping NEO data on to stella grids and calculating all quantites needed for higher order GK calculations. ------------------------------ !
! ================================================================================================================================================================================= !
! 
! NEO uses pitch angle cosine, ξ = v∥​/v, and normalised energy, E, for velocity coordinates. stella uses v∥​ and μ. A remapping of the NEO data on to the stella grids is required. 
! 
! We must also reconstruct the NEO H_1 from the amplitudes provided by out.neo.f: see https://gacode.io/neo/outputs.html#neo-out-neo-f for details on the Legendre/Laguerre
! representation of H_1 in terms of the amplitudes. 
!
! This representation can also be used to calculate the derivatives of H_1 in v∥​ and μ. The spatial gradients of H_1 (and F_1) may be calculated by finite difference methods. 
!
! ================================================================================================================================================================================= !

module neoclassical_terms_neo
    implicit none
 
    ! Load debug flags?

    ! Make routines available to other modules.
    public :: read_parameters_neoclassical
    public :: neoclassical_is_enabled
    public :: init_neoclassical_terms_neo
    public :: finish_neoclassical_terms_neo

    public :: neo_h, dneo_h_dpsi, dneo_h_dz, dneo_h_dvpa, dneo_h_dmu           ! Will represent NEO's distribution and its derivatives in real space and velocity space.  
    public :: neo_phi, dneo_phi_dpsi, dneo_phi_dz                              ! Will represents NEO's ϕ^1_0 and its derivatives in real space. 
    public :: d2neo_h_dzdmu, d2neo_h_dzdvpa

    public :: neo_dens                                                         ! Holds moments of the NEO H_1 distribution. neo_dens is the zeroeth order moment. 
    public :: neo_fac

    public :: initialised_neoclassical_terms_neo

    private

    logical :: include_neoclassical_terms                                           
    integer :: neo_option_switch                                                   ! Should be = 2 for NEO.
    integer :: nradii                                                              ! Can only be 3.
    real    :: drho                                                                ! Typically taken as 0.01.

    integer :: iz, unit
    character(len=128) :: filename

    real, dimension(:, :, :), allocatable :: neo_h, dneo_h_dpsi, dneo_h_dz
    real, dimension(:), allocatable :: neo_phi, dneo_phi_dpsi, dneo_phi_dz
    real, dimension(:, :, :), allocatable :: dneo_h_dvpa, dneo_h_dmu, d2neo_h_dzdmu, d2neo_h_dzdvpa
    real, dimension(:, :), allocatable :: neo_dens, neo_dens_right, neo_dens_left
    real, dimension(:, :), allocatable :: neo_u_par, neo_u_par_right, neo_u_par_left
    real, dimension(:, :), allocatable :: neo_fac
    
    ! For testing velocity derivatives.
    real, dimension(:, :), allocatable :: neo_dens_vpa_deriv
    real, dimension(:, :), allocatable :: neo_dens_mu_deriv

    logical :: initialised_neoclassical_terms_neo = .false.

contains

! ================================================================================================================================================================================= !
! ----------------------------------------------------------------------- Read the neoclassical namelist. ------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine read_parameters_neoclassical
        use namelist_neoclassical_input, only: read_namelist_neoclassical_input
              
        implicit none

        call read_namelist_neoclassical_input(include_neoclassical_terms, neo_option_switch, nradii, drho) ! Read the neoclassical namelist.
        call broadcast_neoclassical_namelist_options                                                       ! Broadcast this information to all other processors. 

    contains 

        subroutine broadcast_neoclassical_namelist_options
            use mp, only: broadcast

            implicit none

            call broadcast(include_neoclassical_terms)
            call broadcast(neo_option_switch)
            call broadcast(nradii)
            call broadcast(drho)
       
        end subroutine broadcast_neoclassical_namelist_options
    end subroutine read_parameters_neoclassical    


! ================================================================================================================================================================================= !
! -------------------------------------------------- A logical function telling stella when to include NEO corrections and routines. ---------------------------------------------- !
! ================================================================================================================================================================================= !

    logical function neoclassical_is_enabled()
        neoclassical_is_enabled = include_neoclassical_terms .and. neo_option_switch == 2
    end function neoclassical_is_enabled


! ================================================================================================================================================================================= !
! ----------------------------------------------------------------- Initiliase the neoclassical terms. ---------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine init_neoclassical_terms_neo
        !MP.
        use mp, only: proc0, broadcast

        ! Parallelisation. 
        use parallelisation_layouts, only: vmu_lo

        ! Grids. 
        use grids_z, only: nzgrid
        use grids_velocity, only: nvpa, nmu 
        use grids_species, only: nspec
        
        ! For NEO's data. 
        use NEO_interface, only: read_basic_neo_files, read_neo_f_and_phi, neo_grid_data, neo_version_data    
        use neoclassical_diagnostics, only: write_neo_phi_on_stella_z_grid_diagnostic, write_neo_distribution_on_stella_grids_diagnostic
        use neoclassical_diagnostics, only: write_distribution_moment_diagnostic

        implicit none

        real, dimension(:, :, :, :, :), allocatable :: neo_h_hat_in, neo_h_hat_right_in, neo_h_hat_left_in                ! Holds vectors for reconstructing NEO H_1 on 3 flux surfaces.
        real, dimension(:, :), allocatable :: neo_phi_in, neo_phi_right_in, neo_phi_left_in                               ! Holds NEO ϕ^1_0 on 3 flux surfaces.

        ! Intermediate arrays hold NEO h_hat data evaluated on the stella z grid, but on the NEO velocity grids.
        real, dimension(:, :, :, :, :), allocatable :: neo_h_hat_z_grid, neo_h_hat_right_z_grid, neo_h_hat_left_z_grid 
        real, dimension(:), allocatable :: neo_phi_right, neo_phi_left                                                 ! Holds NEO ϕ^1_0 data evaluated on the stella z grid.

        ! Intermediate arrays for holding data associated with the NEO H_1 z derivative.
        real, dimension(:, :, :, :, :), allocatable :: dneo_h_hat_dz_z_grid, dneo_h_dz_global

        ! Holds NEO H_1 data evaluated on the stella z, v∥​ and μ grids. Since ϕ^1_0 is independent of velocity variables, there are no accompanying arrays for ϕ^1_0 here.
        real, dimension(:, :, :, :, :), allocatable :: neo_h_global, neo_h_global_right, neo_h_global_left

        ! Holds NEO H_1 derivative data evaluated on the stella z, v∥​ and μ grids.
        real, dimension(:, :, :, :, :), allocatable :: dneo_h_dvpa_global, dneo_h_dmu_global, d2neo_h_dzdmu_global, d2neo_h_dzdvpa_global

        ! Holds NEO H_1 data evaluated on the stella z, v∥​ and μ grids, compacted into 3 indicdes.
        real, dimension(:, :, :), allocatable :: neo_h_right, neo_h_left                 

        integer :: iz
        integer :: surface_index
        integer :: output_unit 

        type(neo_grid_data) :: neo_grid
        type(neo_version_data) :: neo_version

        if (initialised_neoclassical_terms_neo) return ! Initialise only once.
        initialised_neoclassical_terms_neo = .true.

        if (proc0) then
            call read_basic_neo_files(neo_grid, neo_version)
            write(*, '(A)') '############################################################'
            write(*, '(A)') '               O(Δ¹ δ¹ v_th/L F₀) CORRECTIONS               ' 
            write(*, '(A)') '############################################################'
            write(*, '("Reading neo files created on system ",A," at ",A," (commit : ",A,")")') neo_version%system, neo_version%date, neo_version%commit
            write(*, '(A)') '                                                            '
        end if
       
        call broadcast(neo_grid%n_species)     ! Broadcast the NEO grid structure to other processors.
        call broadcast(neo_grid%n_energy)
        call broadcast(neo_grid%n_xi)
        call broadcast(neo_grid%n_theta)
        call broadcast(neo_grid%n_radial)
    
        if (.not. allocated(neo_grid%theta)) allocate(neo_grid%theta(neo_grid%n_theta))
        if (.not. allocated(neo_grid%radius)) allocate(neo_grid%radius(neo_grid%n_radial))
    
        call broadcast(neo_grid%theta)
        call broadcast(neo_grid%radius)
       
        ! Allocate all of the arrays that will be used in the higher order corrections. 
        if (.not. allocated(neo_h)) allocate(neo_h(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))
        if (.not. allocated(neo_phi)) allocate(neo_phi(-nzgrid:nzgrid))
        if (.not. allocated(dneo_h_dpsi)) allocate(dneo_h_dpsi(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))
        if (.not. allocated(dneo_phi_dpsi)) allocate(dneo_phi_dpsi(-nzgrid:nzgrid))
        if (.not. allocated(dneo_h_dz)) allocate(dneo_h_dz(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))
        if (.not. allocated(dneo_phi_dz)) allocate(dneo_phi_dz(-nzgrid:nzgrid))
        if (.not. allocated(dneo_h_dvpa)) allocate(dneo_h_dvpa(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))
        if (.not. allocated(dneo_h_dmu)) allocate(dneo_h_dmu(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))
        if (.not. allocated(d2neo_h_dzdmu)) allocate(d2neo_h_dzdmu(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))
        if (.not. allocated(d2neo_h_dzdvpa)) allocate(d2neo_h_dzdvpa(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))
        if (.not. allocated(neo_dens)) allocate(neo_dens(-nzgrid:nzgrid, neo_grid%n_species))
        if (.not. allocated(neo_fac)) allocate(neo_fac(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc))

        ! Allocate all temporary arrays needed for initilization. 
        call allocate_temp_arrays

        if (proc0) then
            call read_neo_f_and_phi(neo_h_hat_in, neo_phi_in, neo_grid)                                 
            call read_neo_f_and_phi(neo_h_hat_right_in, neo_phi_right_in, neo_grid, suffix = '.right')   
            call read_neo_f_and_phi(neo_h_hat_left_in, neo_phi_left_in, neo_grid, suffix = '.left')    
        end if        
        
        ! Broadcast the read in data. 
        call broadcast(neo_h_hat_in) ; call broadcast(neo_phi_in)             
        call broadcast(neo_h_hat_right_in) ; call broadcast(neo_phi_right_in)
        call broadcast(neo_h_hat_left_in) ; call broadcast(neo_phi_left_in)
   
        ! Interpolates the NEO h_hat data on to the stella z grid for all three flux surfaces. h_hat data is still on the NEO velocity grids at this stage.
        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_in, neo_grid, 1, neo_h_hat_z_grid, .false.)             
        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_right_in, neo_grid, 1, neo_h_hat_right_z_grid, .false.)  
        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_left_in, neo_grid, 1, neo_h_hat_left_z_grid, .false.)    

        ! Interpolates the NEO potential data on to the stella z grid for all three flux surfaces.
        call get_neo_phi_on_stella_z_grid(neo_phi_in, neo_grid, 1, neo_phi, .false.)             
        call get_neo_phi_on_stella_z_grid(neo_phi_right_in, neo_grid, 1, neo_phi_right, .false.)       
        call get_neo_phi_on_stella_z_grid(neo_phi_left_in, neo_grid, 1, neo_phi_left, .false.)         

        ! PHI DIAGNOSTIC.
        ! if (proc0) then
            ! call write_neo_phi_on_stella_z_grid_diagnostic(neo_grid, neo_phi, "neo_phi_on_stella_z_grid_central_surface")
            ! call write_neo_phi_on_stella_z_grid_diagnostic(neo_grid, neo_phi_left, "neo_phi_on_stella_z_grid_left_surface")
            ! call write_neo_phi_on_stella_z_grid_diagnostic(neo_grid, neo_phi_right, "neo_phi_on_stella_z_grid_right_surface")
        ! end if

        ! Now reconstruct H_1 on stellas v∥​ and μ grids for central surface.
        call get_neo_h_on_stella_grids(neo_h_hat_right_z_grid, neo_grid, 1, neo_h_global_right)                                             
        call get_neo_h_on_stella_grids(neo_h_hat_left_z_grid, neo_grid, 1, neo_h_global_left)      
        call get_neo_h_on_stella_grids(neo_h_hat_z_grid, neo_grid, 1, neo_h_global)
 
        ! DISTRIBUTION DIAGNOSTIC.
        ! if (proc0) then
            ! call write_neo_distribution_on_stella_grids_diagnostic(neo_grid, neo_h_global, 1, "neo_h_on_stella_z_grid_central_surface")
        ! end if

        ! We also require the zeroeth order moment for the distribution for testing purposes. 
        ! The central surface data is also used in the QN condition. Construct the moments from the global data rather than the local data. 
        call get_neo_moment(neo_h_global, neo_phi, neo_grid, "zeroeth", neo_dens)
        ! call get_neo_moment(neo_h_global_right, neo_phi_right, neo_grid, "zeroeth", neo_dens_right)
        ! call get_neo_moment(neo_h_global_left, neo_phi_left, neo_grid, "zeroeth", neo_dens_left)

        ! Write out the density arrays to output files.  
        ! if (proc0) then
            ! call write_distribution_moment_diagnostic(neo_grid, neo_dens, "zeroeth_moment_on_stella_z_grid_central_surface")
            ! call write_distribution_moment_diagnostic(neo_grid, neo_dens_right, "zeroeth_moment_on_stella_z_grid_right_surface")
            ! call write_distribution_moment_diagnostic(neo_grid, neo_dens_left, "zeroeth_moment_on_stella_z_grid_left_surface")
        ! end if

        ! We also require the first order moment for the distribution for testing purposes. Construct the moments from the global data rather than the local data. 
        ! call get_neo_moment(neo_h_global, neo_phi, neo_grid, "first", neo_u_par)
        ! call get_neo_moment(neo_h_global_right, neo_phi_right, neo_grid, "first", neo_u_par_right)
        ! call get_neo_moment(neo_h_global_left, neo_phi_left, neo_grid, "first", neo_u_par_left)

        ! Write out the parallel flow arrays to output files.  
        ! if (proc0) then
            ! call write_distribution_moment_diagnostic(neo_grid, neo_u_par, "first_moment_on_stella_z_grid_central_surface")
            ! call write_distribution_moment_diagnostic(neo_grid, neo_u_par_right, "first_moment_on_stella_z_grid_right_surface")
            ! call write_distribution_moment_diagnostic(neo_grid, neo_u_par_left, "first_moment_on_stella_z_grid_left_surface")
        ! end if

        ! From the first order moments (parallel flows), we can calculate the bootstrap current by taking the field line average. 
        ! call get_bootstrap_current(neo_u_par, "bootstrap_current_central_surface")
        ! call get_bootstrap_current(neo_u_par_right, "bootstrap_current_right_surface")
        ! call get_bootstrap_current(neo_u_par_left, "bootstrap_current_left_surface")

        ! Now compact distribution into 3 indices for use in the GK equation and also for calculating the derivatives.  
        do iz = -nzgrid, nzgrid
            call distribute_vmus_over_procs(neo_h_global(iz, :, :, :, 1), neo_h(iz, :, 1))      
            call distribute_vmus_over_procs(neo_h_global_right(iz, :, :, :, 1), neo_h_right(iz, :, 1))
            call distribute_vmus_over_procs(neo_h_global_left(iz, :, :, :, 1), neo_h_left(iz, :, 1))
        end do        

        ! Now that we have H_1 (not normalised to the Maxwellian here) and ϕ^1_0 for the three flux surfaces, the radial, z, v∥​ and μ derivatives are needed.
        ! Calculate the psi derivative arrays via a central difference method. 
        dneo_h_dpsi = (neo_h_right - neo_h_left) / (2 * drho)
        dneo_phi_dpsi = (neo_phi_right - neo_phi_left) / (2 * drho) 

        ! PHI DIAGNOSTIC.
        ! if (proc0) then
            ! call write_neo_phi_on_stella_z_grid_diagnostic(neo_grid, dneo_phi_dpsi, "dneo_phi_drho_on_stella_z_grid_central_surface")
        ! end if

        ! z derivatives can also be obtained via splines.
        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_in, neo_grid, 1, dneo_h_hat_dz_z_grid, .true.)
        call get_neo_h_on_stella_grids(dneo_h_hat_dz_z_grid, neo_grid, 1, dneo_h_dz_global)
         
        ! Now compact distribution into 3 indices for use in the GK equation and also for calculating the derivatives.  
        do iz = -nzgrid, nzgrid
            call distribute_vmus_over_procs(dneo_h_dz_global(iz, :, :, :, 1), dneo_h_dz(iz, :, 1))
        end do

        call get_neo_phi_on_stella_z_grid(neo_phi_in, neo_grid, 1, dneo_phi_dz, .true.)

        ! PHI DIAGNOSTIC.
        ! if (proc0) then
            ! call write_neo_phi_on_stella_z_grid_diagnostic(neo_grid, dneo_phi_dz, "dneo_phi_dz_on_stella_z_grid_central_surface")
        ! end if

        ! We will need the derivatives of the distribution with respect to stellas velocity variables.
        call get_neo_h_velocity_derivs_on_stella_grids(neo_h_hat_z_grid, neo_grid, 1, dneo_h_dvpa_global, dneo_h_dmu_global)

        ! We also need the mixed z and velocity derivatives. 
        call get_neo_h_velocity_derivs_on_stella_grids(dneo_h_hat_dz_z_grid, neo_grid, 1, d2neo_h_dzdvpa_global, d2neo_h_dzdmu_global)

        ! DISTRIBUTION DIAGNOSTIC.
        ! if (proc0) then
            ! call write_neo_distribution_on_stella_grids_diagnostic(neo_grid, dneo_h_dvpa_global, 1, "neo_h_vpa_deriv_on_stella_z_grid_central_surface")
            ! call write_neo_distribution_on_stella_grids_diagnostic(neo_grid, dneo_h_dmu_global, 1, "neo_h_mu_deriv_on_stella_z_grid_central_surface")
        ! end if

       ! Collapse derivative arrays into 3 dimensions for use in GKE. 
        do iz = -nzgrid, nzgrid
            call distribute_vmus_over_procs(dneo_h_dvpa_global(iz, :, :, :, 1), dneo_h_dvpa(iz, :, 1))
            call distribute_vmus_over_procs(dneo_h_dmu_global(iz, :, :, :, 1), dneo_h_dmu(iz, :, 1))
            call distribute_vmus_over_procs(d2neo_h_dzdvpa_global(iz, :, :, :, 1), d2neo_h_dzdvpa(iz, :, 1))
            call distribute_vmus_over_procs(d2neo_h_dzdmu_global(iz, :, :, :, 1), d2neo_h_dzdmu(iz, :, 1))
        end do
             
        ! For testing purposes, we also take the moments of the vpa and mu derivatives of neo_h. This is done with the global neo_h data. 
        ! call get_neo_h_velocity_derivative_moment(dneo_h_dvpa_global, neo_h_global, neo_phi, neo_grid, "vpa", neo_dens_vpa_deriv)
        ! call get_neo_h_velocity_derivative_moment(dneo_h_dmu_global, neo_h_global, neo_phi, neo_grid, "mu", neo_dens_mu_deriv)        

        ! if (proc0) then
            ! call write_distribution_moment_diagnostic(neo_grid, neo_dens_vpa_deriv, "neo_distribution_vpa_deriv_first_moment_on_stella_z_grid_central_surface") 
            ! call write_distribution_moment_diagnostic(neo_grid, neo_dens_mu_deriv, "neo_distribution_mu_deriv_first_moment_on_stella_z_grid_central_surface")
        ! end if

        ! Check for Nan's in NEO data being passed to the GKE. 
        ! call check_real_nans_r3(neo_h, "neo_h")
        ! call check_real_nans_r3(dneo_h_dpsi, "dneo_h_dpsi")
        ! call check_real_nans_r3(dneo_h_dz, "dneo_h_dz")
        ! call check_real_nans_r3(dneo_h_dvpa, "dneo_h_dvpa")
        ! call check_real_nans_r3(dneo_h_dmu, "dneo_h_dmu")
        
        ! Finally, deallocate all temporary arrays.
        call deallocate_temp_arrays

    contains


    ! ========================================================================================================================================================================== !
    ! ------------------------------------ Allocates temporary arrays for H_1, ϕ^1_0 and their derivatives, used in init_neoclassical_terms_neo. ------------------------------- !
    ! ========================================================================================================================================================================== !

    subroutine allocate_temp_arrays
        implicit none

        ! Allocate the h_hat arrays on NEO grids.
        if (.not. allocated(neo_h_hat_in)) allocate(neo_h_hat_in(neo_grid%n_theta, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial))  
        if (.not. allocated(neo_h_hat_right_in)) allocate(neo_h_hat_right_in(neo_grid%n_theta, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial))
        if (.not. allocated(neo_h_hat_left_in)) allocate(neo_h_hat_left_in(neo_grid%n_theta, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial))
      
        ! Allocate the ϕ^1_0 arrays on NEOs z grid.
        if (.not. allocated(neo_phi_in)) allocate(neo_phi_in(neo_grid%n_theta, neo_grid%n_radial))   
        if (.not. allocated(neo_phi_right_in)) allocate(neo_phi_right_in(neo_grid%n_theta, neo_grid%n_radial))
        if (.not. allocated(neo_phi_left_in)) allocate(neo_phi_left_in(neo_grid%n_theta, neo_grid%n_radial))

         ! Allocate the NEO h_hat arrays on stellas z grid. Data will still be on NEOs velocity grids at this point. 
        if (.not. allocated(neo_h_hat_z_grid)) allocate(neo_h_hat_z_grid(-nzgrid:nzgrid, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial))
        if (.not. allocated(neo_h_hat_right_z_grid)) allocate(neo_h_hat_right_z_grid(-nzgrid:nzgrid, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial)) 
        if (.not. allocated(neo_h_hat_left_z_grid)) allocate(neo_h_hat_left_z_grid(-nzgrid:nzgrid, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial))   
        
        ! Allocate the z derivative of the NEO h_hat arrays on stellas z grid. Data will still be on NEOs velocity grids at this point.
        if (.not. allocated(dneo_h_hat_dz_z_grid)) allocate(dneo_h_hat_dz_z_grid(-nzgrid:nzgrid, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial))
        if (.not. allocated(dneo_h_dz_global)) allocate(dneo_h_dz_global(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))
        if (.not. allocated(d2neo_h_dzdmu_global)) allocate(d2neo_h_dzdmu_global(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))
        if (.not. allocated(d2neo_h_dzdvpa_global)) allocate(d2neo_h_dzdvpa_global(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))      

        ! Allocate the ϕ^1_0 arrays (for left and right flux surfaces) on stellas z grid.                                                                                                    
        if (.not. allocated(neo_phi_right)) allocate(neo_phi_right(-nzgrid:nzgrid)) 
        if (.not. allocated(neo_phi_left)) allocate(neo_phi_left(-nzgrid:nzgrid))    

        ! Allocate the NEO H_1 5D arrays on stellas z and velocity grids.
        if (.not. allocated(neo_h_global)) allocate(neo_h_global(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))
        if (.not. allocated(neo_h_global_right)) allocate(neo_h_global_right(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))
        if (.not. allocated(neo_h_global_left)) allocate(neo_h_global_left(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))

        ! Allocate the NEO H_1 3D arrays (for left and right flux surfaces) on stellas z and velocity grids.
        if (.not. allocated(neo_h_right)) allocate(neo_h_right(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial)) 
        if (.not. allocated(neo_h_left)) allocate(neo_h_left(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))   

        ! Allocate the NEO zeroeth order moments, for testing purposes.
        if (.not. allocated(neo_dens_right)) allocate(neo_dens_right(-nzgrid:nzgrid, neo_grid%n_species))
        if (.not. allocated(neo_dens_left)) allocate(neo_dens_left(-nzgrid:nzgrid, neo_grid%n_species))

        ! Allocate the NEO first order moments, for testing purposes.
        if (.not. allocated(neo_u_par)) allocate(neo_u_par(-nzgrid:nzgrid, neo_grid%n_species))
        if (.not. allocated(neo_u_par_right)) allocate(neo_u_par_right(-nzgrid:nzgrid, neo_grid%n_species))
        if (.not. allocated(neo_u_par_left)) allocate(neo_u_par_left(-nzgrid:nzgrid, neo_grid%n_species))

        ! FOR TESTING PURPOSES.
        if (.not. allocated(dneo_h_dmu_global)) allocate(dneo_h_dmu_global(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))
        if (.not. allocated(dneo_h_dvpa_global)) allocate(dneo_h_dvpa_global(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))

        if (.not. allocated(neo_dens_vpa_deriv)) allocate(neo_dens_vpa_deriv(-nzgrid:nzgrid, neo_grid%n_species))
        if (.not. allocated(neo_dens_mu_deriv)) allocate(neo_dens_mu_deriv(-nzgrid:nzgrid, neo_grid%n_species))

    end subroutine allocate_temp_arrays
    

    ! ========================================================================================================================================================================== !
    ! ----------------------------------- Deallocates temporary arrays for H_1, ϕ^1_0 and their derivatives, used in init_neoclassical_terms_neo. ------------------------------ !
    ! ========================================================================================================================================================================== !

    subroutine deallocate_temp_arrays
        implicit none

        if (allocated(neo_h_hat_in)) deallocate(neo_h_hat_in) 
        if (allocated(neo_h_hat_right_in)) deallocate(neo_h_hat_right_in)
        if (allocated(neo_h_hat_left_in)) deallocate(neo_h_hat_left_in)
        if (allocated(neo_phi_in)) deallocate(neo_phi_in)  
        if (allocated(neo_phi_right_in)) deallocate(neo_phi_right_in)
        if (allocated(neo_phi_left_in)) deallocate(neo_phi_left_in)
        if (allocated(neo_h_hat_z_grid)) deallocate(neo_h_hat_z_grid)     
        if (allocated(neo_h_hat_right_z_grid)) deallocate(neo_h_hat_right_z_grid)    
        if (allocated(neo_h_hat_left_z_grid)) deallocate(neo_h_hat_left_z_grid)
        if (allocated(dneo_h_hat_dz_z_grid)) deallocate(dneo_h_hat_dz_z_grid)
        if (allocated(dneo_h_dz_global)) deallocate(dneo_h_dz_global)
        if (allocated(d2neo_h_dzdvpa_global)) deallocate(d2neo_h_dzdvpa_global)
        if (allocated(d2neo_h_dzdmu_global)) deallocate(d2neo_h_dzdmu_global)   
        if (allocated(neo_phi_right)) deallocate(neo_phi_right) 
        if (allocated(neo_phi_left)) deallocate(neo_phi_left)
        if (allocated(neo_h_global)) deallocate(neo_h_global)       
        if (allocated(neo_h_global_right)) deallocate(neo_h_global_right)
        if (allocated(neo_h_global_left)) deallocate(neo_h_global_left)
        if (allocated(neo_h_right)) deallocate(neo_h_right)   
        if (allocated(neo_h_left)) deallocate(neo_h_left)  
        if (allocated(neo_dens_right)) deallocate(neo_dens_right)
        if (allocated(neo_dens_left)) deallocate(neo_dens_left)
        if (allocated(neo_u_par)) deallocate(neo_u_par)
        if (allocated(neo_u_par_right)) deallocate(neo_u_par_right)
        if (allocated(neo_u_par_left)) deallocate(neo_u_par_left)
        if (allocated(dneo_h_dvpa_global)) deallocate(dneo_h_dvpa_global)
        if (allocated(neo_dens_vpa_deriv)) deallocate(neo_dens_vpa_deriv)
        if (allocated(neo_dens_mu_deriv)) deallocate(neo_dens_mu_deriv) 
 
    end subroutine deallocate_temp_arrays

    
    end subroutine init_neoclassical_terms_neo


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------- Finish the neoclassical terms. ---------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine finish_neoclassical_terms_neo
        implicit none

        if (allocated(neo_h)) deallocate(neo_h)
        if (allocated(neo_phi)) deallocate(neo_phi)
        if (allocated(dneo_h_dpsi)) deallocate(dneo_h_dpsi)
        if (allocated(dneo_phi_dpsi)) deallocate(dneo_phi_dpsi)
        if (allocated(dneo_h_dz)) deallocate(dneo_h_dz)
        if (allocated(dneo_phi_dz)) deallocate(dneo_phi_dz)
        if (allocated(dneo_h_dvpa)) deallocate(dneo_h_dvpa)
        if (allocated(dneo_h_dmu)) deallocate(dneo_h_dmu)
        if (allocated(d2neo_h_dzdmu)) deallocate(d2neo_h_dzdmu)
        if (allocated(d2neo_h_dzdvpa)) deallocate(d2neo_h_dzdvpa)
        if (allocated(neo_dens)) deallocate(neo_dens)
        if (allocated(neo_fac)) deallocate(neo_fac)

        initialised_neoclassical_terms_neo = .false.
        
    end subroutine finish_neoclassical_terms_neo


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------- Utilities. -------------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

! ================================================================================================================================================================================= !
! ------------------------------------------------------------- Collapse distributuion arrays from 5 dimensions to 3. ------------------------------------------------------------- !
! ================================================================================================================================================================================= !

   subroutine distribute_vmus_over_procs(local, distributed)
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      use grids_z, only: nzgrid

      implicit none

      real, dimension(:, :, :), intent(in) :: local
      real, dimension(vmu_lo%llim_proc:), intent(out) :: distributed

      integer :: ivmu, iv, imu, is

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         distributed(ivmu) = local(iv, imu, is)
      end do
   end subroutine distribute_vmus_over_procs


! ================================================================================================================================================================================= !
! ---------------------------------- Interpolate the NEO h_hat data from the NEO θ grid to the stella z grid for a specified flux surface. ---------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine get_neo_h_hat_on_stella_z_grid(neo_h_hat_in, neo_grid, surface_index, neo_h_hat_z_grid, derivative)
        ! Grids. 
        use grids_z, only: nzgrid, zed
        
        ! NEO.
        use NEO_interface, only: neo_grid_data

        ! Splines. 
        use periodic_splines, only: periodic_spline, new_periodic_spline, delete_periodic_spline
        
        ! Constants.
        use constants, only: twopi

        implicit none

        real, intent(in)  :: neo_h_hat_in(:, :, :, :, :)    
        real, intent(out) :: neo_h_hat_z_grid(-nzgrid:, :, :, :, :)
        type(neo_grid_data), intent(in) :: neo_grid
        integer, intent(in) :: surface_index
        logical, intent(in), optional :: derivative
            
        ! Local variables.
        logical :: local_derivative
        type(periodic_spline) :: the_spline
        integer :: ix, ie, is, iz                  

        if (present(derivative)) then
            local_derivative = derivative
        else
            local_derivative = .false.
        end if

        do ix = 1, neo_grid%n_xi + 1
            do ie = 1, neo_grid%n_energy + 1
                do is = 1, neo_grid%n_species
                    the_spline = new_periodic_spline(neo_grid%theta, neo_h_hat_in(:, ix, ie, is, surface_index), twopi)

                    if (local_derivative) then
                        do iz = -nzgrid, nzgrid
                            neo_h_hat_z_grid(iz, ix, ie, is, surface_index) = the_spline%derivative(zed(iz))
                        end do
                    else
                        do iz = -nzgrid, nzgrid
                            neo_h_hat_z_grid(iz, ix, ie, is, surface_index) = the_spline%interpolate(zed(iz))
                        end do
                    end if

                    call delete_periodic_spline(the_spline)
                end do
            end do
        end do
    end subroutine get_neo_h_hat_on_stella_z_grid


! ================================================================================================================================================================================= !
! --------------------------------- Interpolate the NEO ϕ^1_0 data from the NEO θ grid to the stella z grid for a specified flux surface. ----------------------------------------- !
! ------------------------------------------------------------ This can also handle the z derivative of ϕ^1_0. -------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine get_neo_phi_on_stella_z_grid(neo_phi_in, neo_grid, surface_index, neo_phi, derivative)
        ! Grids.
        use grids_z, only: nzgrid, zed

        ! NEO. 
        use NEO_interface, only: neo_grid_data
  
        ! Splines.
        use periodic_splines, only: periodic_spline, new_periodic_spline, delete_periodic_spline

        ! Utilities. 
        use constants, only: twopi
    
        implicit none

        real, intent(in)  :: neo_phi_in(:, :)    
        type(neo_grid_data), intent(in) :: neo_grid
        integer, intent(in) :: surface_index
        real, intent(out) :: neo_phi(-nzgrid:)                 
        logical, intent(in), optional :: derivative
        
        logical :: local_derivative
        type(periodic_spline) :: the_spline
        integer :: iz

        the_spline = new_periodic_spline(neo_grid%theta, neo_phi_in(:, surface_index), twopi)

        if (present(derivative)) then
            local_derivative = derivative
        else
            local_derivative = .false.
        end if

        do iz = -nzgrid, nzgrid
            if (local_derivative) then
                neo_phi(iz) = the_spline%derivative(zed(iz))
            else
                neo_phi(iz) = the_spline%interpolate(zed(iz))
            end if
        end do
 
        call delete_periodic_spline(the_spline)

    end subroutine get_neo_phi_on_stella_z_grid


! ================================================================================================================================================================================= !
! --------------------------------------------- Reconstruct NEO H_1 on stella z, v∥​ and μ grids for a selected flux surface. ------------------------------------------------------ !
! ================================================================================================================================================================================= !
!
! We have neo_h_hat which are the amplitudes in the Laguerre-Legendre basis. The full neo_h(r, θ, x_a, ξ) (for one particle species) is given by:
!
! neo_h(r, θ, x_a, ξ)  = F_MB(r, θ, x_a) * [Sum_ie Sum_ix{L_ie^(k(ix)+1/2)(x_a^2) P_ix(ξ) neo_h_hat}].
!
! Here E = v/(sqrt(2)*v_th,a), ξ = v∥​/v, k(ix) = 0, ix = 0 and k(ix) = 1, ix > 0. L are the associated Laguerre polynomials and P are the Legendre 
! polynomials. See: https://gacode.io/neo/outputs.html#neo-out-neo-f for more details.
!
! We need to evaluate this on stellas v-space grids for v∥​ and μ. For each point on these grids, we can evaluate x_a and ξ, then evaluate the polynomials.  
! The result will be neo_h on the stellas grids for a selected flux surface. 
!
! ================================================================================================================================================================================= !

    subroutine get_neo_h_on_stella_grids(neo_h_hat_z_grid, neo_grid, surface_index, neo_h, suffix)
	! Grids. 
        use grids_z, only: nzgrid, zed
        use grids_velocity, only: vpa, nvpa, mu, nmu
        use grids_species, only: spec  

        ! Geometry.
        use geometry, only: bmag

        ! NEO.                  
        use NEO_interface, only: neo_grid_data

        ! bmag is calculated and stored in stella as a two-dimensional array depending on alpha (for dealing with non-axisymmetric stellarator configurations) and z. 
        ! Since we are dealing with an axisymmetric tokamak, bmag should be independent of alpha for a given z coordinate.
        ! Here, loops over ia are therefore not needed and all calculations are based on ia = 1. 
 
        implicit none

        real, intent(in)  :: neo_h_hat_z_grid(-nzgrid:, :, :, :, :)
        type(neo_grid_data), intent(in) :: neo_grid
        integer, intent(in) :: surface_index
        real, intent(out) :: neo_h(-nzgrid:, :, :, :, :)         

        real :: xi_in, E_in
        integer :: iz, iv, imu, is

        integer :: unit
        character(len=256) :: filename
        character(len=*), intent(in), optional :: suffix           
        
        do is = 1, neo_grid%n_species
            do iz = -nzgrid, nzgrid
                do iv = 1, nvpa
                    do imu = 1, nmu
                        ! Calculate (ξ_in, E_in) from the (v∥​, μ) stella grid point.
			xi_in = vpa(iv) / sqrt( vpa(iv)**2 + 2 * bmag(1, iz) * mu(imu) ) 
			E_in = vpa(iv)**2 + 2 * bmag(1, iz) * mu(imu)

                        ! Construct neo_h at the given (ξ_in, E_in) point. 
                        neo_h(iz, iv, imu, is, surface_index) = get_neo_h_at_xi_energy(neo_h_hat_z_grid(iz, :, :, is, surface_index), E_in, xi_in, neo_grid)
                    end do
                end do
            end do
        end do

    end subroutine get_neo_h_on_stella_grids


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------- Constructs H_1 for a single (ξ_in, E_in) pair. ----------------------------------------------------------------- !
! ================================================================================================================================================================================= !
!
! Performs the sum over the passed distribution function weighted by Legendre and associated Laguerre basis functions at a single (ξ_in, E_in).
! The result is the neo distribution function, H_1 (not normalised to the Maxwellian, F_0), on the passed (ξ_in, E_in) values.
!
! ================================================================================================================================================================================= !

    pure real function get_neo_h_at_xi_energy(neo_h_hat_z_grid, E_in, xi_in, neo_grid) result(the_sum)
        use polynomials, only: get_legendre, get_laguerre
        use NEO_interface, only: neo_grid_data
    
        implicit none
    
        real, dimension(:, :), intent(in) :: neo_h_hat_z_grid ! A 2-D slice of the 5-dimensional neo_h_hat_z_grid data, at fixed z, species and surface index. 
        real, intent(in) :: xi_in, E_in
        type(neo_grid_data), intent(in) :: neo_grid
    
        ! Local variables.
        integer :: ie, ix
        real, dimension(:, :), allocatable :: L ! To holds two sets of Laguerre polynomials.
        real, dimension(:), allocatable    :: P ! To hold Legendre Polynomials.
       
        allocate(P(neo_grid%n_xi + 1))
        allocate(L(neo_grid%n_energy + 1, 2))
        
        the_sum = 0.0

        ! Evaluate the Legendre polynomials.
        call get_legendre(xi_in, P)

        ! We need to evaluate the Laguerre function with two different values of k, the first for the first ξ index and the other for all other ξ indices.
        call get_laguerre(E_in, 0.5, L(:, 1))
        call get_laguerre(E_in, 1.5, L(:, 2))

        ! Now iterate over data and perform the sum. First add in the piece for the first ξ index. 
        do ie = 1, neo_grid%n_energy + 1
       	    the_sum = the_sum + P(1) * L(ie, 1) * neo_h_hat_z_grid(1, ie)

            do ix = 2, neo_grid%n_xi + 1
                the_sum = the_sum + P(ix) * L(ie, 2) * neo_h_hat_z_grid(ix, ie)
            end do
        end do

        deallocate(P)
        deallocate(L)
    end function get_neo_h_at_xi_energy


! ================================================================================================================================================================================= !
! -------------------------------------------------------------- Constructs dH_1/dξ for a single (ξ_in, E_in) pair. --------------------------------------------------------------- !
! ================================================================================================================================================================================= !
!
! Performs the sum over the passed distribution function weighted by Legendre derivative and associated Laguerre basis functions at a single (ξ_in, E_in).
! The result is the derivative of the neo distribution function with respect to ξ, dH_1/dξ (not normalised to the Maxwellian, F_0), on the passed (ξ_in, E_in) values.
!
! ================================================================================================================================================================================= !

    pure real function get_neo_h_at_xi_deriv_energy(neo_h_hat_z_grid, E_in, xi_in, neo_grid) result(the_sum)
        use polynomials, only: get_legendre_deriv, get_laguerre
        use NEO_interface, only: neo_grid_data
    
        implicit none
    
        real, dimension(:, :), intent(in) :: neo_h_hat_z_grid ! A 2-D slice of the 5-dimensional neo_h_hat_z_grid data, at fixed z, species and surface index. 
        real, intent(in) :: xi_in, E_in
        type(neo_grid_data), intent(in) :: neo_grid
    
        ! Local variables. 
        integer :: ie, ix
        real, dimension(:, :), allocatable :: L ! To holds two sets of Laguerre polynomials.
        real, dimension(:), allocatable    :: dP ! To hold Legendre derivative polynomials.
                
        allocate(dP(neo_grid%n_xi + 1))
        allocate(L(neo_grid%n_energy + 1, 2))

        the_sum = 0.0

        ! Evaluate the Legendre polynomials.
        call get_legendre_deriv(xi_in, dP)

        ! We need to evaluate the Laguerre function with two different values of k, the first for the first ξ index and the other for all other ξ indices.
        call get_laguerre(E_in, 0.5, L(:, 1))
        call get_laguerre(E_in, 1.5, L(:, 2))
        
        ! Now iterate over data and perform the sum. First add in the piece for the first ξ index. 
        do ie = 1, neo_grid%n_energy + 1
       	    the_sum = the_sum + dP(1) * L(ie, 1) * neo_h_hat_z_grid(1, ie)

            ! Now add in all the other xi indicies.
            do ix = 2, neo_grid%n_xi + 1
                the_sum = the_sum + dP(ix) * L(ie, 2) * neo_h_hat_z_grid(ix, ie)
            end do
        end do

        deallocate(dP)
        deallocate(L)
    end function get_neo_h_at_xi_deriv_energy


! ================================================================================================================================================================================= !
! -------------------------------------------------------------- Constructs dH_1/dE for a single (ξ_in, E_in) pair. --------------------------------------------------------------- !
! ================================================================================================================================================================================= !
!
! Performs the sum over the passed distribution function weighted by Legendre and associated Laguerre derivative basis functions at a single (ξ_in, E_in).
! The result is the derivative of the neo distribution function with respect to E, dH_1/dE (not normalised to the Maxwellian, F_0), on the passed (ξ_in, E_in) values.
!
! ================================================================================================================================================================================= !

    pure real function get_neo_h_at_xi_energy_deriv(neo_h_hat_z_grid, E_in, xi_in, neo_grid) result(the_sum)
        use polynomials, only: get_legendre, get_laguerre_deriv
        use NEO_interface, only: neo_grid_data

        implicit none

        real, dimension(:, :), intent(in) :: neo_h_hat_z_grid ! A 2-D slice of the 5-dimensional neo_h_hat_z_grid data, at fixed z, species and surface index. 
        real, intent(in) :: xi_in, E_in
        type(neo_grid_data), intent(in) :: neo_grid

        ! Local variables.
        integer :: ie, ix
        real, dimension(:), allocatable    :: P ! To hold Legendre Polynomials.
        real, dimension(:, :), allocatable :: dL ! To holds two sets of Laguerre derivative polynomials.

        allocate(P(neo_grid%n_xi + 1))
        allocate(dL(neo_grid%n_energy + 1, 2))

        the_sum = 0.0

        ! Evaluate the Legendre polynomials.
        call get_legendre(xi_in, P)

        ! We need to evaluate the Laguerre function with two different values of k, the first for the first ξ index and the other for all other ξ indices.

        call get_laguerre_deriv(E_in, 0.5, dL(:, 1))
        call get_laguerre_deriv(E_in, 1.5, dL(:, 2))

        ! Now iterate over data and perform the sum. 
        do ie = 1, neo_grid%n_energy + 1
            ! First add in the piece for the first ξ index at a given E index.
            the_sum = the_sum + P(1) * dL(ie, 1) * neo_h_hat_z_grid(1, ie)

            ! Now add in all the other ξ indicies.
            do ix = 2, neo_grid%n_xi + 1
                the_sum = the_sum + P(ix) * dL(ie, 2) * neo_h_hat_z_grid(ix, ie)
            end do
        end do

        deallocate(P)
        deallocate(dL)
    end function get_neo_h_at_xi_energy_deriv


! ================================================================================================================================================================================= !
! -------------------------------------- Reconstruct NEO dH_1/dv∥|_μ and dH_1/dμ|_v∥ on stella z, v∥​ and μ grids for a selected flux surface. ------------------------------------- !
! ================================================================================================================================================================================= !
!
! We can calculate via the chain rule: dH_1/dv∥​|_μ = v∥​ dH/dE|_ξ + (1 - ξ^2)/sqrt{2E} dH/dξ|_E 
!
! We can also calculate: dH_1/dμ|_v∥​ = B_0 dH/dE|_ξ - (ξ B_0)/2E dH/dξ}|_E
!
! ================================================================================================================================================================================= !

    subroutine get_neo_h_velocity_derivs_on_stella_grids(neo_h_hat_z_grid, neo_grid, surface_index, dneo_h_dvpa_global, dneo_h_dmu_global)
        use grids_z, only: nzgrid, zed
        use grids_velocity, only: vpa, nvpa, mu, nmu
        use geometry, only: bmag
        
        ! NEO.
        use NEO_interface, only: neo_grid_data
        use neoclassical_diagnostics, only: write_neo_distribution_on_stella_grids_diagnostic 
 
        implicit none

        real, intent(in)                :: neo_h_hat_z_grid(-nzgrid:, :, :, :, :)
        type(neo_grid_data), intent(in) :: neo_grid
        integer, intent(in)             :: surface_index
        real, intent(out)               :: dneo_h_dvpa_global(-nzgrid:, :, :, :, :)
        real, intent(out)               :: dneo_h_dmu_global(-nzgrid:, :, :, :, :)

        real    :: xi_in, E_in
        integer :: iz, iv, imu, is
        real    :: dneo_h_dxi, dneo_h_dE

        do is = 1, neo_grid%n_species
            do iz = -nzgrid, nzgrid
                do iv = 1, nvpa
                    do imu = 1, nmu
                        ! Calculate (ξ_in, E_in) from the (v∥​, μ) stella grid point.
			xi_in = vpa(iv) / sqrt(vpa(iv)**2 + 2 * bmag(1, iz) * mu(imu)) 
			E_in = vpa(iv)**2 + 2 * bmag(1, iz) * mu(imu)

                        ! Construct dH_1/dξ|_E at the given (ξ_in, E_in) point. 
                        dneo_h_dxi = get_neo_h_at_xi_deriv_energy(neo_h_hat_z_grid(iz, :, :, is, surface_index), E_in, xi_in, neo_grid)

                        ! Construct dH_1/dE​|_ξ at the given (ξ_in, E_in) point.
                        dneo_h_dE = get_neo_h_at_xi_energy_deriv(neo_h_hat_z_grid(iz, :, :, is, surface_index), E_in, xi_in, neo_grid)

                        ! Construct dH_1/dv∥|_μ at the given (ξ_in, E_in) point. 
                        dneo_h_dvpa_global(iz, iv, imu, is, surface_index) = 2 * vpa(iv) * dneo_h_dE + 2 * bmag(1, iz) * mu(imu)  * dneo_h_dxi / ( E_in**(1.5) )                           
 
                        ! Construct dH_1/dμ|_v∥ at the given (ξ_in, E_in) point.
                        dneo_h_dmu_global(iz, iv, imu, is, surface_index) = 2 * bmag(1, iz) * dneo_h_dE - bmag(1, iz) * vpa(iv) * dneo_h_dxi / ( E_in**(1.5) ) 
                    end do
                end do
            end do
        end do

    end subroutine get_neo_h_velocity_derivs_on_stella_grids


! ================================================================================================================================================================================= !
! ---------------------------------- Calculate the zeroeth order moment of H_1. This is used in the field equations and also for testing purposes. -------------------------------- !
! ================================================================================================================================================================================= !

    subroutine get_neo_moment(neo_h_global, neo_phi, neo_grid, moment_type, neo_mom)
        ! Grids.
        use grids_z, only: nzgrid
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac, nvpa, vpa, nmu
        use grids_species, only: spec

        ! NEO data.
        use NEO_interface, only: neo_grid_data 

        ! Calculations. 
        use calculations_velocity_integrals, only: integrate_vmu

        implicit none

        real, dimension(-nzgrid:, :, :, :, :), intent(in) :: neo_h_global
        real, dimension(-nzgrid:)                         :: neo_phi 
        type(neo_grid_data), intent(in)                   :: neo_grid
        character(len=*), intent(in)                      :: moment_type
        real, dimension(-nzgrid:, :), intent(out)         :: neo_mom        
        
        ! Local variables.
        integer           :: is, imu, iv, iz 
        real, allocatable :: tmp(:, :, :, :)
        real :: vel_weight

        allocate(tmp(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species))
        tmp = 0.0
        neo_mom = 0.0

        ! Calculate the physical F_1. 
        do is = 1, neo_grid%n_species
            do iv = 1, nvpa
                do imu = 1, nmu
                    do iz = -nzgrid, nzgrid           
                        select case (trim(moment_type))

                        case ("zeroeth")
                            tmp(iz, iv, imu, is) = ( neo_h_global(iz, iv, imu, is, 1) + spec(1)%z * neo_phi(iz) ) * maxwell_vpa(iv, 1) * maxwell_mu(1, iz, imu, 1) * maxwell_fac(1)

                        case ("first")
                            tmp(iz, iv, imu, is) = vpa(iv) * sqrt(2.0) * spec(is)%stm * neo_h_global(iz, iv, imu, is, 1) &
                            * maxwell_vpa(iv, 1) * maxwell_mu(1, iz, imu, 1) * maxwell_fac(1)
                        end select
                    end do
                end do
            end do
        end do 

        ! Integrate over velocity space.
        do is = 1, neo_grid%n_species
            do iz = -nzgrid, nzgrid
                call integrate_vmu(tmp(iz, :, :, is), iz, neo_mom(iz, is))
            end do
        end do

        ! Deallocate the temporary array.
        deallocate(tmp)
    end subroutine get_neo_moment


! ================================================================================================================================================================================= !
! - Calculate neo_fac, a combination of NEO terms that is found in several of the higher order corrections. Calculate this once here rather than calculating individually in each - ! 
! ----------------------------------------- of the gyrokinetic terms modules. This is given by: neo_fac = v∥/B * ∂H_1/∂μ|_v∥ - ∂H_1/∂v∥|_μ. --------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine get_neo_fac(dneo_h_dvpa, dneo_h_dmu, neo_fac)
        ! Parallelisation.
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        ! Grids.
        use grids_z, only: nzgrid
        use grids_velocity, only: maxwell_vpa, maxwell_mu
        use grids_velocity, only: vpa

        ! Geometry. 
        use geometry, only: bmag

        implicit none 

        real, dimension(-nzgrid:, vmu_lo%llim_proc:, :), intent(in) :: dneo_h_dvpa
        real, dimension(-nzgrid:, vmu_lo%llim_proc:, :), intent(in) :: dneo_h_dmu
        real, dimension(-nzgrid:, vmu_lo%llim_proc:), intent(out)   :: neo_fac

        ! Local variables.
        integer :: is, imu, iv, ivmu, iz
        
        ! Calculate neo_fac.  
        ! Iterate over velocity space.
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
                        
            do iz = -nzgrid, nzgrid
                neo_fac(iz, ivmu) = ( vpa(iv) * dneo_h_dmu(iz, ivmu, 1) / bmag(1, iz) - dneo_h_dvpa(iz, ivmu, 1) ) 
            end do
        end do

    end subroutine get_neo_fac 


! ================================================================================================================================================================================= !
! ----------------------------- Calculate the H_1/F_1 velocity derivatives via finite difference methods on the stella grids, for testing purposes. ------------------------------- !
! ================================================================================================================================================================================= !

    subroutine get_neo_h_velocity_derivative_moment(neo_h_deriv_global, neo_h_global, neo_phi, neo_grid, vel_type, neo_mom)
        ! Grids.
        use grids_z, only: nzgrid
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac, nvpa, vpa, nmu, mu
        use grids_species, only: spec

        ! Geometry. 
        use geometry, only: bmag

        ! NEO data.
        use NEO_interface, only: neo_grid_data 

        ! Calculations. 
        use calculations_velocity_integrals, only: integrate_vmu

        ! Constants.
        use constants, only: pi
 
        implicit none

        real, intent(in)                :: neo_h_deriv_global(-nzgrid:, :, :, :, :)
        real, intent(in)                :: neo_h_global(-nzgrid:, :, :, :, :)
        real, intent(in)                :: neo_phi(-nzgrid:)
        type(neo_grid_data), intent(in) :: neo_grid
        character(len=*), intent(in)    :: vel_type
        real, intent(out)               :: neo_mom(-nzgrid:, :)         

        ! Local variables.
        integer           :: is, imu, iv, iz
        real, allocatable :: tmp(:, :, :, :)

        allocate(tmp(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species))

        tmp = 0.0
        neo_mom = 0.0        
 
        do is = 1, neo_grid%n_species
            do iv = 1, nvpa
                do imu = 1, nmu
                    do iz = -nzgrid, nzgrid           
                        select case (trim(vel_type))
                        case ("mu")
                            tmp(iz, iv, imu, is) = mu(imu) * ( neo_h_deriv_global(iz, iv, imu, is, 1) &
                            - 2 * bmag(1, iz) * ( neo_h_global(iz, iv, imu, is, 1) - spec(is)%z * neo_phi(iz) ) ) * maxwell_vpa(iv, is) * maxwell_mu(1, iz, imu, is) * maxwell_fac(is)

                        case ("vpa")
                            tmp(iz, iv, imu, is) = vpa(iv) * ( neo_h_deriv_global(iz, iv, imu, is, 1) &
                            - 2 * vpa(iv) * ( neo_h_global(iz, iv, imu, is, 1) - spec(is)%z * neo_phi(iz) ) ) * maxwell_vpa(iv, is) * maxwell_mu(1, iz, imu, is) * maxwell_fac(is)
              
                        end select
                    end do
                end do
            end do
        end do

        ! Integrate over velocity space.
        do is = 1, neo_grid%n_species
            do iz = -nzgrid, nzgrid
                call integrate_vmu(tmp(iz, :, :, is), iz, neo_mom(iz, is))
            end do
        end do

        ! Deallocate the temporary array.
        deallocate(tmp)
    end subroutine get_neo_h_velocity_derivative_moment


! ================================================================================================================================================================================= !
! -------------------------------------------- Calculate the bootstrap current from the first order moment on a given flux surface. ----------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine get_bootstrap_current(neo_mom, file_tag)
        ! Grids. 
        use grids_z, only: nzgrid
        use grids_species, only: spec
      
        ! Geometry.
        use geometry, only: bmag, dl_over_b

        implicit none

        real, dimension(-nzgrid:, :), intent(in) :: neo_mom
        character(len=*), intent(in)             :: file_tag

        ! Local variables. 
        integer :: it, iz, unit
        real    :: avg, bootstrap
        character(len=256) :: filename

        bootstrap = 0.0
        avg = 0.0 
      
        do iz = -nzgrid, nzgrid 
            bootstrap = bootstrap + ( neo_mom(iz, 1) / ( sqrt(2.0) * spec(1)%stm ) - neo_mom(iz, 2) / ( sqrt(2.0) * spec(2)%stm ) ) * bmag(1, iz) * dl_over_b(1, iz)
            avg = avg + dl_over_b(1, iz)
        end do
      
        bootstrap = bootstrap / avg

        ! Write out the bootstrap current to the output file. 
        if (len_trim(file_tag) > 0) then
            filename = trim(file_tag) // ".dat"
        end if

        unit = 99

        open(newunit=unit, file=trim(filename), status='replace', action='write')

        write(unit,'(A)') '# ================================================= #'
        write(unit,'(A)') '#       Diagnostic output: Bootstrap Current        #'
        write(unit,'(A)') '# ================================================= #'
        write(unit,'(ES16.8)') bootstrap
                 
        close(unit)
   end subroutine get_bootstrap_current


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------ Check for Nan's. --------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine check_real_nans_r3(array, name)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      
        implicit none
    
        ! Strictly Rank 3 real array
        real, dimension(:, :, :), intent(in) :: array 
        character(len=*), intent(in)         :: name
    
        logical :: is_bad
    
        ! Modern compilers can now easily map this to the elemental intrinsic
        is_bad = any(.not. ieee_is_finite(array))
    
        if (is_bad) then
            write(*,'(A,A,A)') ">>> MATH ERROR: NaN or Inf detected in [", name, "] <<<"
            ! Optional: Print the location of the first NaN to help debugging
            ! error stop "Execution halted: Invalid real number found in Rank 3 array."
        else
            write(*,'(A,A,A)') "Check Passed: All values are finite in [", name, "]"
        end if
    
end subroutine check_real_nans_r3


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------ End Module. -------------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

end module neoclassical_terms_neo
