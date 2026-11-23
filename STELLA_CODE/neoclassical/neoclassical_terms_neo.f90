! ================================================================================================================================================================================= !
! ---------------------------------- Routines for mapping NEO data on to stella grids and calculating all quantites needed for HO GK simulations. --------------------------------- !
! ================================================================================================================================================================================= !
! 
! NEO data is provided as a function of the poloidal angle, ϑ, with a domain of -π ≤ ϑ < π. stella uses a field following coordinate, z = ϑ, for Miller geometries, with a domain of
! -nπ ≤ ϑ < nπ, where n is a positive integer defined by the user. NEO data must be interpolated from the NEO ϑ grid to stellas z grid. This is done using perioidic splines. z
! derivatives of NEO quantities are also calculated with perioidic splines.  
!
! NEO uses pitch angle cosine, ξ = v∥​/v, and normalised energy, E, for velocity coordinates. stella uses v∥​ and μ. H_1 data must be correctly mapped from one velocity space to the 
! other using the Legendre/Laguerre representation: https://gacode.io/neo/outputs.html#neo-out-neo-f. This representation can also be used to calculate the derivatives of H_1 
! in v∥​ and μ. 
!
! The radial gradient of H_1 and F_1 may be calculated by central difference methods. 
!
! The neoclassical arrays used in the GK-Maxwell equations are owned by arrays_neoclassical.f90. This module populates them with NEO data. They carry a field-line (α) index as 
! their first dimension, in preparation for general (non-axisymmetric) fields:
!
!   - Global distributions and their derivatives : (nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec)
!   - ϕ^1_0 and its derivatives                  : (nalpha, -nzgrid:nzgrid)
!   - Velocity-distributed distributions         : (nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc)
!
! NEO is axisymmetric, so every quantity is computed once for ia = 1 and then copied to ia = 2, ..., nalpha (see copy_alpha_one_to_all). 
! The NEO radial index is not carried in the final arrays, since n_radial = 1 in every NEO run.
!
! ================================================================================================================================================================================= !

module neoclassical_terms_neo
    use arrays_neoclassical, only: include_neoclassical_terms
    use arrays_neoclassical, only: neo_h, dneo_h_drho, dneo_h_dz
    use arrays_neoclassical, only: neo_phi, dneo_phi_drho, dneo_phi_dz
    use arrays_neoclassical, only: dneo_h_dmu, dneo_h_dvpa
    use arrays_neoclassical, only: dneo_h_dmu_global, dneo_h_dvpa_global
    use arrays_neoclassical, only: neo_vpa_fac, neo_mu_fac
    use arrays_neoclassical, only: neo_vpa_fac_global, neo_mu_fac_global
    use arrays_neoclassical, only: neo_zed_fac

    implicit none

    public :: read_parameters_neoclassical
    ! True when NEO is the selected drift kinetic provider.
    public :: neo_is_enabled
    
    public :: init_neoclassical_terms_neo
    public :: finish_neoclassical_terms_neo

    public :: distribute_vmus_over_procs

    public :: initialised_neoclassical_terms_neo

    public :: get_neo_moment, get_neo_h_velocity_derivative_moment

    ! Transport verification, also used by neoclassical_terms_sfincs to check SFINCS data on stella grids against NEO.
    public :: stella_to_neo_units_data
    public :: get_stella_to_neo_units
    public :: verify_neo_transport_on_surface
    public :: open_neo_check_file, close_neo_check_file_and_summarise
    public :: add_neo_check_notice

    private

    ! Copies the ia = 1 slice of an array (α is the first dimension) on to all other ia.
    interface copy_alpha_one_to_all
        module procedure copy_alpha_one_to_all_2d
        module procedure copy_alpha_one_to_all_3d
        module procedure copy_alpha_one_to_all_5d
    end interface copy_alpha_one_to_all
                                           
    integer :: neo_option_switch           ! Should be = 2 for NEO.
    integer :: nradii                      ! Currently can only be 3.
    real    :: drho                        ! Typically taken as 0.01.

    integer :: iz, unit
    character(len=128) :: filename

    ! NEO-specific intermediate used to build neo_zed_fac; not used in the GKE, so it is not part of arrays_neoclassical.
    ! (nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec).
    real, dimension(:, :, :, :, :), allocatable :: neo_zed_fac_global

    ! DIAGNOSTICS. 
    real, dimension(:, :), allocatable :: neo_dens, neo_dens_right, neo_dens_left
    real, dimension(:, :), allocatable :: neo_u_par, neo_u_par_right, neo_u_par_left
    real, dimension(:, :), allocatable :: neo_dens_vpa_deriv
    real, dimension(:, :), allocatable :: neo_dens_mu_deriv

    ! Ratios of stella reference units to NEO normalisation units (stella / NEO), used to verify the NEO data on stella grids against out.neo.transport:
    !   n   = n_ref / n_norm,  T = T_ref / T_norm,  m = m_ref / m_norm,  vth = v_th,ref / v_norm = sqrt(2 T / m),  B = B_ref / B_unit.
    ! The type must be defined here, in the specification part, before any variable of this type is declared.
    type stella_to_neo_units_data
        real :: n = 1.0, T = 1.0, m = 1.0, vth = sqrt(2.0), B = 1.0
    end type stella_to_neo_units_data

    type(stella_to_neo_units_data) :: stella_to_neo_units

    ! VERIFICATION against out.neo.transport. The full check is written to <run_name>//neo_check_ext; only a one-line summary goes to the screen.
    ! check_unit is the unit every verification routine writes to (the screen until the file is opened).
    character(len=*), parameter :: neo_check_ext = '.neo_transport_check.dat'
    ! Relative tolerance for every stella/NEO comparison.
    real, parameter :: neo_check_rel_tol = 0.01
    integer, parameter :: check_pass = 0, check_warn = 1, check_fail = 2
    integer :: check_unit = 6
    ! File extension, summary title and name of the code whose data are on stella's grids (set by open_neo_check_file).
    character(len=64) :: check_ext = neo_check_ext
    character(len=64) :: check_title = 'NEO transport check'
    character(len=16) :: check_code_label = 'stella'
    ! Counts over all surfaces: diagnostics compared, passed, warned (outside tolerance but resolution-limited by cancellation), failed,
    ! and setup notices (missing files, rotation in NEO, inconsistent gradients or species ratios) that are not diagnostics themselves.
    integer :: n_checks = 0, n_pass = 0, n_warn = 0, n_fail = 0, n_notices = 0


    logical :: initialised_neoclassical_terms_neo = .false.

contains

! ================================================================================================================================================================================= !
! ----------------------------------------------------------------------- Read the neoclassical namelist. ------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine read_parameters_neoclassical
        use namelist_neoclassical_input, only: read_namelist_neoclassical_input
              
        implicit none

        call read_namelist_neoclassical_input(include_neoclassical_terms, neo_option_switch, nradii, drho)
        call broadcast_neoclassical_namelist_options                                                       

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
! ---------------------------------------------------- Logical function telling stella when to use NEO as the neoclassical provider. ---------------------------------------------- !
! ================================================================================================================================================================================= !

    logical function neo_is_enabled()
        neo_is_enabled = include_neoclassical_terms .and. neo_option_switch == 2
    end function neo_is_enabled


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------- Initiliase the neoclassical terms. ----------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine init_neoclassical_terms_neo
        !MP.
        use mp, only: proc0, broadcast

        ! Parallelisation. 
        use parallelisation_layouts, only: vmu_lo
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        ! Grids. 
        use grids_z, only: nzgrid, delzed
        use grids_kxky, only: nalpha
        use grids_velocity, only: nvpa, nmu, mu    
        use grids_species, only: nspec, spec
        
        ! Geometry. 
        use geometry, only: get_dzed, bmag, dbdzed

        ! For NEO's data. 
        use NEO_interface, only: read_basic_neo_files, read_neo_f_and_phi, read_neo_transport, read_neo_equil
        use NEO_interface, only: neo_grid_data, neo_version_data, neo_transport_data
        use NEO_interface, only: neo_equil_data, neo_species_data    
        use neoclassical_diagnostics, only: write_neo_phi_on_stella_z_grid_diagnostic
        use neoclassical_diagnostics, only: write_neo_distribution_on_stella_grids_diagnostic
        use neoclassical_diagnostics, only: write_distribution_moment_diagnostic

        ! Geometry.
        use geometry, only: get_dzed, bmag

        ! Neoclassical arrays (allocation is a no-op if already done).
        use arrays_neoclassical, only: allocate_neoclassical_arrays

        implicit none

        ! Temporaray variables. 
        ! Holds basis vectors for reconstructing NEO H_1 on 3 flux surfaces.
        real, dimension(:, :, :, :, :), allocatable :: neo_h_hat_in, neo_h_hat_right_in, neo_h_hat_left_in                
        ! Holds NEO ϕ^1_0 on 3 flux surfaces.
        real, dimension(:, :), allocatable :: neo_phi_in, neo_phi_right_in, neo_phi_left_in                               

        ! Intermediate arrays hold NEO h_hat data evaluated on the stella z grid, but on the NEO velocity grids.
        real, dimension(:, :, :, :, :), allocatable :: neo_h_hat_z_grid, neo_h_hat_right_z_grid, neo_h_hat_left_z_grid
        ! Holds NEO ϕ^1_0 data evaluated on the stella z grid (for ia = 1).
        real, dimension(:), allocatable :: neo_phi_right, neo_phi_left                                                 

        ! Intermediate arrays for holding data associated with the NEO H_1 z derivative.
        ! dneo_h_dz_global is on the stella z, v∥​ and μ grids for ia = 1: (-nzgrid:nzgrid, nvpa, nmu, nspec).
        real, dimension(:, :, :, :, :), allocatable :: dneo_h_hat_dz_z_grid
        real, dimension(:, :, :, :), allocatable :: dneo_h_dz_global

        ! Holds NEO H_1 data evaluated on the stella z, v∥​ and μ grids for ia = 1: (-nzgrid:nzgrid, nvpa, nmu, nspec).
        ! Since ϕ^1_0 is independent of velocity variables, there are no accompanying arrays for ϕ^1_0 here.
        real, dimension(:, :, :, :), allocatable :: neo_h_global, neo_h_global_right, neo_h_global_left

        ! Holds NEO H_1 data on the left and right flux surfaces for ia = 1, compacted into 2 indices: (-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc).
        real, dimension(:, :), allocatable :: neo_h_right, neo_h_left                 


        integer :: iz, ivmu, iv, imu, is, ia
        integer :: surface_index
 
        type(neo_grid_data) :: neo_grid
        type(neo_transport_data) :: neo_transport_central, neo_transport_left, neo_transport_right
        type(neo_version_data) :: neo_version
        type(neo_equil_data)   :: neo_equil
        type(neo_species_data) :: neo_species

        ! Profiles on the left and right surfaces (used only to verify the NEO data against out.neo.transport on those surfaces).
        ! Species masses and charges are the same on every surface, so the central neo_species is used throughout.
        type(neo_equil_data)           :: neo_equil_left, neo_equil_right
        type(stella_to_neo_units_data) :: stella_to_neo_units_left, stella_to_neo_units_right

        ! ===================================================================================== !

        ! All quantities are computed on the first field line, then copied to all others.
        ia = 1

        if (initialised_neoclassical_terms_neo) return
        initialised_neoclassical_terms_neo = .true.

        if (proc0) then
            call read_basic_neo_files(neo_grid, neo_version, neo_equil, neo_species)
            write(*, '(A)') '############################################################'
            write(*, '(A)') '               O(Δ¹ δ¹ v_th/L F₀) CORRECTIONS               ' 
            write(*, '(A)') '############################################################'
            write(*, '("Reading neo files created on system ",A," at ",A," (commit : ",A,")")') neo_version%system, neo_version%date, neo_version%commit
            write(*, '(A)') '                                                            '
        end if
    
        ! Broadcasst grid sizes to all processors. 
        call broadcast(neo_grid%n_species)     
        call broadcast(neo_grid%n_energy)
        call broadcast(neo_grid%n_xi)
        call broadcast(neo_grid%n_theta)
        call broadcast(neo_grid%n_radial)

        if (.not. allocated(neo_grid%theta)) allocate(neo_grid%theta(neo_grid%n_theta))
        if (.not. allocated(neo_grid%radius)) allocate(neo_grid%radius(neo_grid%n_radial))

        call broadcast(neo_grid%theta)
        call broadcast(neo_grid%radius)
       
        ! The neoclassical arrays used in the GKE are owned by arrays_neoclassical. They should already be allocated 
        ! during initialisation; this call only guarantees it (it returns immediately if they are).
        call allocate_neoclassical_arrays

        ! NEO-specific intermediate: (nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec).
        if (.not. allocated(neo_zed_fac_global)) allocate(neo_zed_fac_global(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec))

        ! Allocate all temporary arrays needed for initilization. 
        call allocate_temp_arrays

        ! Read in raw NEO H_1 and ϕ^1_0 data for each surface and broadcast to all processprs.
        if (proc0) then
            call read_neo_f_and_phi(neo_h_hat_in, neo_phi_in, neo_grid)                                 
            call read_neo_f_and_phi(neo_h_hat_right_in, neo_phi_right_in, neo_grid, suffix = '.right')   
            call read_neo_f_and_phi(neo_h_hat_left_in, neo_phi_left_in, neo_grid, suffix = '.left')    
        end if        

        ! Read in transport data on all three surfaces. Not verbose: the raw values and any missing files are reported in the verification file instead.
        if (proc0) then
            call read_neo_transport(neo_transport_central, neo_grid, verbose = .false.)
            call read_neo_transport(neo_transport_left, neo_grid, suffix = ".left",  verbose = .false.)
            call read_neo_transport(neo_transport_right, neo_grid, suffix = ".right", verbose = .false.)
        end if 

        call broadcast(neo_h_hat_in) ; call broadcast(neo_phi_in)             
        call broadcast(neo_h_hat_right_in) ; call broadcast(neo_phi_right_in)
        call broadcast(neo_h_hat_left_in) ; call broadcast(neo_phi_left_in)
   
        ! Interpolates the NEO h_hat data on to the stella z grid for all three flux surfaces. h_hat data is still on the NEO velocity grids at this stage.
        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_in, neo_grid, 1, neo_h_hat_z_grid, .false.)             
        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_right_in, neo_grid, 1, neo_h_hat_right_z_grid, .false.)  
        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_left_in, neo_grid, 1, neo_h_hat_left_z_grid, .false.)    

        ! Interpolates the NEO potential data on to the stella z grid for all three flux surfaces.
        call get_neo_phi_on_stella_z_grid(neo_phi_in, neo_grid, 1, neo_phi(ia, :), .false.)             
        call get_neo_phi_on_stella_z_grid(neo_phi_right_in, neo_grid, 1, neo_phi_right, .false.)       
        call get_neo_phi_on_stella_z_grid(neo_phi_left_in, neo_grid, 1, neo_phi_left, .false.)         

        ! DIAGNOSTIC.
        if (proc0) then
            call write_neo_phi_on_stella_z_grid_diagnostic(neo_phi(ia, :), "neo_phi_on_stella_z_grid_central_surface")
            call write_neo_phi_on_stella_z_grid_diagnostic(neo_phi_left, "neo_phi_on_stella_z_grid_left_surface")
            call write_neo_phi_on_stella_z_grid_diagnostic(neo_phi_right, "neo_phi_on_stella_z_grid_right_surface")
        end if

        ! Now reconstruct H_1 on stellas v∥​ and μ grids for all three surfaces (ia = 1).
        call get_neo_h_on_stella_grids(neo_h_hat_z_grid, neo_grid, 1, neo_h_global)
        ! The left and right surfaces are needed for the ρ derivatives (via neo_h_right, neo_h_left) and for the transport verification.
        call get_neo_h_on_stella_grids(neo_h_hat_right_z_grid, neo_grid, 1, neo_h_global_right)
        call get_neo_h_on_stella_grids(neo_h_hat_left_z_grid, neo_grid, 1, neo_h_global_left)

        ! DIAGNOSTIC.
        ! call write_neo_distribution_on_stella_grids_diagnostic(neo_grid, neo_h_global, surface_index, 'neo_h_on_stella_grids_central')
 
        ! DIAGNOSTIC.
        call get_neo_moment(neo_h_global, neo_phi(ia, :), "zeroeth", neo_dens)
        ! call get_neo_moment(neo_h_global_right, neo_phi_right, "zeroeth", neo_dens_right)
        ! call get_neo_moment(neo_h_global_left, neo_phi_left, "zeroeth", neo_dens_left)

        ! DIAGNOSTIC.
        ! Write out the density arrays to output files.  
        if (proc0) then
            call write_distribution_moment_diagnostic(neo_dens, "zeroeth_moment_on_stella_z_grid_central_surface")
            ! call write_distribution_moment_diagnostic(neo_dens_right, "zeroeth_moment_on_stella_z_grid_right_surface")
            ! call write_distribution_moment_diagnostic(neo_dens_left, "zeroeth_moment_on_stella_z_grid_left_surface")
        end if

        ! DIAGNOSTIC.
        call get_neo_moment(neo_h_global, neo_phi(ia, :), "first", neo_u_par)
        ! call get_neo_moment(neo_h_global_right, neo_phi_right, "first", neo_u_par_right)
        ! call get_neo_moment(neo_h_global_left, neo_phi_left, "first", neo_u_par_left)

        ! DIAGNOSTIC.
        ! Write out the parallel flow arrays to output files.  
        if (proc0) then
            call write_distribution_moment_diagnostic(neo_u_par, "first_moment_on_stella_z_grid_central_surface")
            ! call write_distribution_moment_diagnostic(neo_u_par_right, "first_moment_on_stella_z_grid_right_surface")
            ! call write_distribution_moment_diagnostic(neo_u_par_left, "first_moment_on_stella_z_grid_left_surface")
        end if

        ! VERIFICATION. Recompute the out.neo.transport quantities from the NEO data on stella grids, convert to NEO units and compare, on each surface.
        ! Only proc0 holds the NEO profiles and the transport data.
        !
        ! The NEO data on each surface are functions of v / v_th,s on that surface, but are evaluated on stella's (v∥, μ) grid in units of the central v_th,s.
        ! Moments of the left/right H_1 are therefore in units of the local n_s and v_th,s, so the flows and the current are converted with that surface's
        ! profiles (stella_to_neo_units_left/right). Φ_1 is normalised to the same T_norm on every surface, so <Φ_1^2> always uses the central conversion.
        ! The left/right surfaces still use the central stella geometry (bmag, dl_over_b, B_ref/B_unit).
        if (proc0) then
            call read_neo_equil(neo_equil_left, neo_grid, suffix = '.left')
            call read_neo_equil(neo_equil_right, neo_grid, suffix = '.right')

            ! Everything below is written to <run_name>.neo_transport_check.dat; a one-line summary is printed at the end.
            call open_neo_check_file

            call get_stella_to_neo_units(stella_to_neo_units, neo_equil, neo_species, neo_grid%n_species, 'central')
            call get_stella_to_neo_units(stella_to_neo_units_left, neo_equil_left, neo_species, neo_grid%n_species, 'left')
            call get_stella_to_neo_units(stella_to_neo_units_right, neo_equil_right, neo_species, neo_grid%n_species, 'right')

            call verify_neo_transport_on_surface('central', neo_h_global, neo_phi(ia, :), neo_transport_central, neo_equil, &
                                                 stella_to_neo_units, stella_to_neo_units)
            call verify_neo_transport_on_surface('left', neo_h_global_left, neo_phi_left, neo_transport_left, neo_equil_left, &
                                                 stella_to_neo_units, stella_to_neo_units_left)
            call verify_neo_transport_on_surface('right', neo_h_global_right, neo_phi_right, neo_transport_right, neo_equil_right, &
                                                 stella_to_neo_units, stella_to_neo_units_right)

            call close_neo_check_file_and_summarise
        end if

        ! Now compact distribution into 3 indices for use in the GKE and also for calculating the derivatives.  
        do iz = -nzgrid, nzgrid
            call distribute_vmus_over_procs(neo_h_global(iz, :, :, :), neo_h(ia, iz, :))      
            call distribute_vmus_over_procs(neo_h_global_right(iz, :, :, :), neo_h_right(iz, :))
            call distribute_vmus_over_procs(neo_h_global_left(iz, :, :, :), neo_h_left(iz, :))
        end do        

        ! z derivatives can also be obtained via splines.
        call get_neo_phi_on_stella_z_grid(neo_phi_in, neo_grid, 1, dneo_phi_dz(ia, :), .true.)

        ! We will need the derivatives of the distribution with respect to stellas velocity variables.
        call get_neo_h_velocity_derivs_on_stella_grids(neo_h_hat_z_grid, neo_grid, 1, dneo_h_dvpa_global(ia, :, :, :, :), dneo_h_dmu_global(ia, :, :, :, :))

        ! Collapse derivative arrays into 3 dimensions for use in GKE. 
        do iz = -nzgrid, nzgrid
            call distribute_vmus_over_procs(dneo_h_dvpa_global(ia, iz, :, :, :), dneo_h_dvpa(ia, iz, :))
            call distribute_vmus_over_procs(dneo_h_dmu_global(ia, iz, :, :, :), dneo_h_dmu(ia, iz, :))
        end do

        ! Calculate the velocity factors appearing in the GKE. 
        call get_velocity_factors(dneo_h_dvpa_global(ia, :, :, :, :), dneo_h_dmu_global(ia, :, :, :, :), neo_h_global, neo_phi(ia, :), &
                                  neo_vpa_fac_global(ia, :, :, :, :), neo_mu_fac_global(ia, :, :, :, :), neo_grid)

        ! Collapse derivative arrays into 3 dimensions for use in GKE. 
        do iz = -nzgrid, nzgrid     
            call distribute_vmus_over_procs(neo_vpa_fac_global(ia, iz, :, :, :), neo_vpa_fac(ia, iz, :))
            call distribute_vmus_over_procs(neo_mu_fac_global(ia, iz, :, :, :), neo_mu_fac(ia, iz, :))
        end do

        ! Calculate the neo_h rho derivative on the central surface at vpa and mu, and the neo_phi rho derivative on the central surface.  
        call get_rho_derivatives(neo_h_right, neo_h_left, neo_phi_right, neo_phi_left, dneo_h_dvpa(ia, :, :), dneo_h_dmu(ia, :, :), drho, dneo_h_drho(ia, :, :), dneo_phi_drho(ia, :))
        
        ! DIAGNOSTIC.
        if (proc0) then
            call write_neo_phi_on_stella_z_grid_diagnostic(dneo_phi_drho(ia, :), "dneo_phi_drho_on_stella_z_grid_central_surface")
        end if

        ! DIAGNOSTIC.
        if (proc0) then
            call write_neo_distribution_on_stella_grids_diagnostic(neo_grid, dneo_h_dvpa_global(ia, :, :, :, :), 1, "neo_h_vpa_deriv_on_stella_z_grid_central_surface")
            call write_neo_distribution_on_stella_grids_diagnostic(neo_grid, dneo_h_dmu_global(ia, :, :, :, :), 1, "neo_h_mu_deriv_on_stella_z_grid_central_surface")
        end if
             
        ! DIAGNOSTIC. 
        call get_neo_h_velocity_derivative_moment(dneo_h_dvpa_global(ia, :, :, :, :), neo_h_global, neo_phi(ia, :), "vpa", neo_dens_vpa_deriv)
        call get_neo_h_velocity_derivative_moment(dneo_h_dmu_global(ia, :, :, :, :), neo_h_global, neo_phi(ia, :), "mu", neo_dens_mu_deriv)        

        ! DIAGNOSTIC.
        if (proc0) then
            call write_distribution_moment_diagnostic(neo_dens_vpa_deriv, "neo_distribution_vpa_deriv_first_moment_on_stella_z_grid_central_surface") 
            call write_distribution_moment_diagnostic(neo_dens_mu_deriv, "neo_distribution_mu_deriv_first_moment_on_stella_z_grid_central_surface")
        end if

        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_in, neo_grid, 1, dneo_h_hat_dz_z_grid, .true.)
        call get_neo_h_on_stella_grids(dneo_h_hat_dz_z_grid, neo_grid, 1, dneo_h_dz_global)
         
        ! The derivative resulting from this is at fixed NEO velocity variables.
        ! Transform this to fixed stella velocity variables. 
        do iz = -nzgrid, nzgrid
            do iv = 1, nvpa
                do imu = 1, nmu
                    do is = 1, nspec
                        dneo_h_dz_global(iz, iv, imu, is) = dneo_h_dz_global(iz, iv, imu, is) &
                        + mu(imu) * dbdzed(ia, iz) * dneo_h_dmu_global(ia, iz, iv, imu, is) &
                        / bmag(ia, iz) 
                    end do 
                end do
            end do
        end do

        ! Now compact distribution into 3 indices for use in the GK equation and also for calculating the derivatives.  
        do iz = -nzgrid, nzgrid
            call distribute_vmus_over_procs(dneo_h_dz_global(iz, :, :, :), dneo_h_dz(ia, iz, :))
        end do


        call get_zed_factor(neo_zed_fac_global(ia, :, :, :, :), dneo_h_dz_global, dneo_phi_dz(ia, :), neo_h_global, neo_phi(ia, :), neo_grid)

        ! Now compact distribution into 3 indices for use in the GKE.  
        do iz = -nzgrid, nzgrid
            call distribute_vmus_over_procs(neo_zed_fac_global(ia, iz, :, :, :), neo_zed_fac(ia, iz, :))
        end do

        ! DIAGNOSTIC.
        if (proc0) then
            call write_neo_phi_on_stella_z_grid_diagnostic(dneo_phi_dz(ia, :), "dneo_phi_dz_on_stella_z_grid_central_surface")
        end if

        ! NEO is axisymmetric, so every field line is identical: copy the ia = 1 data on to ia = 2, ..., nalpha.
        call copy_alpha_one_to_all(neo_h)
        call copy_alpha_one_to_all(dneo_h_drho)
        call copy_alpha_one_to_all(dneo_h_dz)
        call copy_alpha_one_to_all(dneo_h_dmu)
        call copy_alpha_one_to_all(dneo_h_dvpa)
        call copy_alpha_one_to_all(neo_vpa_fac)
        call copy_alpha_one_to_all(neo_mu_fac)
        call copy_alpha_one_to_all(neo_zed_fac)

        call copy_alpha_one_to_all(neo_phi)
        call copy_alpha_one_to_all(dneo_phi_drho)
        call copy_alpha_one_to_all(dneo_phi_dz)

        call copy_alpha_one_to_all(dneo_h_dmu_global)
        call copy_alpha_one_to_all(dneo_h_dvpa_global)
        call copy_alpha_one_to_all(neo_vpa_fac_global)
        call copy_alpha_one_to_all(neo_mu_fac_global)
        call copy_alpha_one_to_all(neo_zed_fac_global)
        
        ! OR MMS TESTING (overwrites all field lines).
        ! call mms_overwrite(neo_h, neo_vpa_fac, neo_mu_fac, neo_mu_fac_global, neo_phi, neo_grid)

        ! Deallocate all temporary arrays.
        call deallocate_temp_arrays

    contains


    ! ========================================================================================================================================================================== !
    ! ------------------------------------------------------- Allocates temporary arrays for H_1, ϕ^1_0 and their derivatives. ------------------------------------------------- !
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
        if (.not. allocated(dneo_h_dz_global)) allocate(dneo_h_dz_global(-nzgrid:nzgrid, nvpa, nmu, nspec))

        ! Allocate the ϕ^1_0 arrays (for left and right flux surfaces) on stellas z grid.                                                                                                    
        if (.not. allocated(neo_phi_right)) allocate(neo_phi_right(-nzgrid:nzgrid)) 
        if (.not. allocated(neo_phi_left)) allocate(neo_phi_left(-nzgrid:nzgrid))    

        ! Allocate the NEO H_1 4D arrays (ia = 1) on stellas z and velocity grids.
        if (.not. allocated(neo_h_global)) allocate(neo_h_global(-nzgrid:nzgrid, nvpa, nmu, nspec))
        if (.not. allocated(neo_h_global_right)) allocate(neo_h_global_right(-nzgrid:nzgrid, nvpa, nmu, nspec))
        if (.not. allocated(neo_h_global_left)) allocate(neo_h_global_left(-nzgrid:nzgrid, nvpa, nmu, nspec))

        ! Allocate the NEO H_1 2D arrays (for left and right flux surfaces, ia = 1) on stellas z and velocity grids.
        if (.not. allocated(neo_h_right)) allocate(neo_h_right(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc)) 
        if (.not. allocated(neo_h_left)) allocate(neo_h_left(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc))   

        ! DIAGNOSTICS.
        ! Allocate the NEO zeroeth order moments.
        if (.not. allocated(neo_dens)) allocate(neo_dens(-nzgrid:nzgrid, neo_grid%n_species))
        if (.not. allocated(neo_dens_right)) allocate(neo_dens_right(-nzgrid:nzgrid, neo_grid%n_species))
        if (.not. allocated(neo_dens_left)) allocate(neo_dens_left(-nzgrid:nzgrid, neo_grid%n_species))

        ! DIAGNOSTICS.
        ! Allocate the NEO first order moments, FOR TESTING PURPOSES.
        if (.not. allocated(neo_u_par)) allocate(neo_u_par(-nzgrid:nzgrid, neo_grid%n_species))
        if (.not. allocated(neo_u_par_right)) allocate(neo_u_par_right(-nzgrid:nzgrid, neo_grid%n_species))
        if (.not. allocated(neo_u_par_left)) allocate(neo_u_par_left(-nzgrid:nzgrid, neo_grid%n_species))

        ! DIAGNOSTICS.
        if (.not. allocated(neo_dens_vpa_deriv)) allocate(neo_dens_vpa_deriv(-nzgrid:nzgrid, neo_grid%n_species))
        if (.not. allocated(neo_dens_mu_deriv)) allocate(neo_dens_mu_deriv(-nzgrid:nzgrid, neo_grid%n_species))

    end subroutine allocate_temp_arrays
    

    ! ========================================================================================================================================================================== !
    ! -------------------------------------------------- Deallocates temporary arrays for H_1, ϕ^1_0 and their derivatives. ---------------------------------------------------- !
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
        if (allocated(neo_phi_right)) deallocate(neo_phi_right) 
        if (allocated(neo_phi_left)) deallocate(neo_phi_left)
        if (allocated(neo_h_global)) deallocate(neo_h_global)
        if (allocated(neo_h_global_right)) deallocate(neo_h_global_right)
        if (allocated(neo_h_global_left)) deallocate(neo_h_global_left)
        if (allocated(neo_h_right)) deallocate(neo_h_right)   
        if (allocated(neo_h_left)) deallocate(neo_h_left)  

        ! DIAGNOSTICS.          
        if (allocated(neo_dens)) deallocate(neo_dens)
        if (allocated(neo_dens_left)) deallocate(neo_dens_left)
        if (allocated(neo_dens_right)) deallocate(neo_dens_right)

        ! DIAGNOSTICS.
        if (allocated(neo_u_par)) deallocate(neo_u_par)
        if (allocated(neo_u_par_left)) deallocate(neo_u_par_left)
        if (allocated(neo_u_par_right)) deallocate(neo_u_par_right)
        
        ! DIAGNOSTICS.
        if (allocated(neo_dens_vpa_deriv)) deallocate(neo_dens_vpa_deriv)
        if (allocated(neo_dens_mu_deriv)) deallocate(neo_dens_mu_deriv)
    
    end subroutine deallocate_temp_arrays

    end subroutine init_neoclassical_terms_neo


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------- Finish the neoclassical terms. ---------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !
!
! Only NEO-specific arrays are deallocated here. The shared neoclassical arrays are deallocated by deallocate_neoclassical_arrays in arrays_neoclassical.
!
! ================================================================================================================================================================================= !

    subroutine finish_neoclassical_terms_neo
        implicit none
 
        if (allocated(neo_zed_fac_global)) deallocate(neo_zed_fac_global)

        initialised_neoclassical_terms_neo = .false.
        
    end subroutine finish_neoclassical_terms_neo


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------- Utilities. -------------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

! ================================================================================================================================================================================= !
! ------------------------------------------------- Copy the ia = 1 slice of an array on to all other field lines, ia = 2, ..., nalpha. ------------------------------------------ !
! ================================================================================================================================================================================= !

    subroutine copy_alpha_one_to_all_2d(arr)
        implicit none

        real, dimension(:, :), intent(in out) :: arr
        integer :: ia

        do ia = 2, size(arr, 1)
            arr(ia, :) = arr(1, :)
        end do
    end subroutine copy_alpha_one_to_all_2d

    subroutine copy_alpha_one_to_all_3d(arr)
        implicit none

        real, dimension(:, :, :), intent(in out) :: arr
        integer :: ia

        do ia = 2, size(arr, 1)
            arr(ia, :, :) = arr(1, :, :)
        end do
    end subroutine copy_alpha_one_to_all_3d

    subroutine copy_alpha_one_to_all_5d(arr)
        implicit none

        real, dimension(:, :, :, :, :), intent(in out) :: arr
        integer :: ia

        do ia = 2, size(arr, 1)
            arr(ia, :, :, :, :) = arr(1, :, :, :, :)
        end do
    end subroutine copy_alpha_one_to_all_5d


! ================================================================================================================================================================================= !
! ------------------------------------------------------------- Collapse distributuion arrays from 5 dimensions to 3. ------------------------------------------------------------- !
! ================================================================================================================================================================================= !
!
! Takes a (nvpa, nmu, nspec) slice at fixed (ia, iz) and packs it into the local vmu_lo range at the same (ia, iz).
!
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
        real, intent(out) :: neo_phi(-nzgrid:)                 ! A single field line: pass e.g. neo_phi(ia, :).
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
! The result will be neo_h on the stellas grids for a selected flux surface, on a single field line (ia = 1): (-nzgrid:nzgrid, nvpa, nmu, nspec). 
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
        real, intent(out) :: neo_h(-nzgrid:, :, :, :)         

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
                        neo_h(iz, iv, imu, is) = get_neo_h_at_xi_energy(neo_h_hat_z_grid(iz, :, :, is, surface_index), E_in, xi_in, neo_grid)
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
                the_sum = the_sum + P(ix) * sqrt(E_in) * L(ie, 2) * neo_h_hat_z_grid(ix, ie)
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
                the_sum = the_sum + dP(ix) * sqrt(E_in) * L(ie, 2) * neo_h_hat_z_grid(ix, ie)
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
    
        use polynomials, only: get_legendre, get_laguerre, get_laguerre_deriv
        use NEO_interface, only: neo_grid_data

        implicit none

        real, dimension(:, :), intent(in) :: neo_h_hat_z_grid
        real, intent(in) :: xi_in, E_in
        type(neo_grid_data), intent(in) :: neo_grid

        integer :: ie, ix
        real :: sqrtE
        real, dimension(:), allocatable    :: P   ! Legendre polynomials.
        real, dimension(:, :), allocatable :: dL  ! (:,1) dL^{1/2}/dE, (:,2) dL^{3/2}/dE, (:,3) L^{3/2}.

        allocate(P(neo_grid%n_xi + 1))
        allocate(dL(neo_grid%n_energy + 1, 3))

        the_sum = 0.0
        sqrtE = sqrt(max(E_in, tiny(E_in)))

        call get_legendre(xi_in, P)
        call get_laguerre_deriv(E_in, 0.5, dL(:, 1))
        call get_laguerre_deriv(E_in, 1.5, dL(:, 2))
        call get_laguerre(E_in, 1.5, dL(:, 3))

        do ie = 1, neo_grid%n_energy + 1
            ! xi index 0: basis L_e^{1/2}(E), unchanged.
            the_sum = the_sum + P(1) * dL(ie, 1) * neo_h_hat_z_grid(1, ie)

            ! xi index > 0: basis sqrt(E) L_e^{3/2}(E) (NEO laguerre_method = 1).
            do ix = 2, neo_grid%n_xi + 1
                the_sum = the_sum + P(ix) * (0.5 * dL(ie, 3) / sqrtE + sqrtE * dL(ie, 2)) &
                                  * neo_h_hat_z_grid(ix, ie)
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
! Outputs are on a single field line (ia = 1): (-nzgrid:nzgrid, nvpa, nmu, nspec). Pass e.g. dneo_h_dvpa_global(ia, :, :, :, :).
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
        real, intent(out)               :: dneo_h_dvpa_global(-nzgrid:, :, :, :)
        real, intent(out)               :: dneo_h_dmu_global(-nzgrid:, :, :, :)

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
                        dneo_h_dvpa_global(iz, iv, imu, is) = 2 * vpa(iv) * dneo_h_dE + 2 * bmag(1, iz) * mu(imu)  * dneo_h_dxi / ( E_in**(1.5) )                           
 
                        ! Construct dH_1/dμ|_v∥ at the given (ξ_in, E_in) point.
                        dneo_h_dmu_global(iz, iv, imu, is) = 2 * bmag(1, iz) * dneo_h_dE - bmag(1, iz) * vpa(iv) * dneo_h_dxi / ( E_in**(1.5) ) 
                    end do
                end do
            end do
        end do

    end subroutine get_neo_h_velocity_derivs_on_stella_grids


! ================================================================================================================================================================================= !
! ---------------------------------------------------------- Calculate the velocity factors appearing in the GKE. ----------------------------------------------------------------- !
! ================================================================================================================================================================================= !
!                                                                                                                                                                                   ! 
! All arrays are on a single field line (ia = 1). Distributions: (-nzgrid:nzgrid, nvpa, nmu, nspec). ϕ^1_0: (-nzgrid:nzgrid).                                                      ! 
!                                                                                                                                                                                   !
! ================================================================================================================================================================================= !

    subroutine get_velocity_factors(dneo_h_dvpa_global, dneo_h_dmu_global, neo_h_global, neo_phi, neo_vpa_fac_global, neo_mu_fac_global, neo_grid)
        ! Grids.
        use grids_z, only: nzgrid
        use grids_velocity, only: vpa, nvpa, nmu
        use grids_species, only: spec

        ! Geometry. 
        use geometry, only: bmag 

        ! NEO data.
        use NEO_interface, only: neo_grid_data

        implicit none    

        real, dimension(-nzgrid:, :, :, :), intent(in)  :: dneo_h_dvpa_global, dneo_h_dmu_global, neo_h_global
        real, dimension(-nzgrid:), intent(in)           :: neo_phi
        real, dimension(-nzgrid:, :, :, :), intent(out) :: neo_vpa_fac_global, neo_mu_fac_global
        type(neo_grid_data), intent(in)                 :: neo_grid

        ! Local variables. 
        integer :: is, imu, iv, iz, ia

        ! ================================================================================ ! 

        ! Compute for ia = 1 only since each field line on our flux surface is identical in an axisymmetric geometry. 
        ia = 1

        do is = 1, neo_grid%n_species
            do iv = 1, nvpa
                do imu = 1, nmu
                    do iz = -nzgrid, nzgrid
                        ! Calculate the parallel velocity derivative factor.
                        neo_vpa_fac_global(iz, iv, imu, is) = dneo_h_dvpa_global(iz, iv, imu, is) - 2.0 * vpa(iv) * neo_h_global(iz, iv, imu, is) &
                        + 2.0 * vpa(iv) * spec(is)%zt * neo_phi(iz)

                        ! Calculate the magnetic moment derivative factor.
                        neo_mu_fac_global(iz, iv, imu, is) = dneo_h_dmu_global(iz, iv, imu, is) - 2.0 * bmag(ia, iz) * neo_h_global(iz, iv, imu, is) &
                        + 2.0 * bmag(ia, iz) * spec(is)%zt * neo_phi(iz)
                    end do
                end do
            end do
        end do
    end subroutine get_velocity_factors


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------- Calculate the zed factor appearing in the GKE. ----------------------------------------------------------------- !
! ================================================================================================================================================================================= !
!                                                                                                                                                                                   ! 
! All arrays are on a single field line (ia = 1). Distributions: (-nzgrid:nzgrid, nvpa, nmu, nspec). ϕ^1_0: (-nzgrid:nzgrid).                                                      ! 
!                                                                                                                                                                                   !
! ================================================================================================================================================================================= !

    subroutine get_zed_factor(neo_zed_fac_global, dneo_h_dz_global, dneo_phi_dz, neo_h_global, neo_phi, neo_grid)
        ! Grids.
        use grids_z, only: nzgrid
        use grids_velocity, only: vpa, mu, nvpa, nmu
        use grids_species, only: spec

        ! Geometry. 
        use geometry, only: bmag, dbdzed 

        ! NEO data.
        use NEO_interface, only: neo_grid_data

        implicit none    

        real, dimension(-nzgrid:, :, :, :), intent(in)  :: dneo_h_dz_global, neo_h_global
        real, dimension(-nzgrid:), intent(in)           :: neo_phi, dneo_phi_dz
        real, dimension(-nzgrid:, :, :, :), intent(out) :: neo_zed_fac_global
        type(neo_grid_data), intent(in)                 :: neo_grid

        ! Local variables. 
        integer :: is, imu, iv, iz, ia

        ! ================================================================================ ! 

        ! Compute for ia = 1 only since each field line on our flux surface is identical in an axisymmetric geometry. 
        ia = 1

        do is = 1, neo_grid%n_species
            do iv = 1, nvpa
                do imu = 1, nmu
                    do iz = -nzgrid, nzgrid
                        neo_zed_fac_global(iz, iv, imu, is) = dneo_h_dz_global(iz, iv, imu, is) - spec(is)%zt * dneo_phi_dz(iz) &
                        - 2.0 * mu(imu) * dbdzed(ia, iz) * ( neo_h_global(iz, iv, imu, is) - spec(is)%zt * neo_phi(iz) )
                    end do
                end do
            end do
        end do

    end subroutine get_zed_factor

! ================================================================================================================================================================================= !
! ------------------------------------------------ Calculate the rho derivative of neo_h and neo_phi on the central flux surface. ----------------------------------------------- !
! ================================================================================================================================================================================= !
!
! All arrays are on a single field line (ia = 1). Distributions: (-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc). ϕ^1_0: (-nzgrid:nzgrid).
!
! ================================================================================================================================================================================= !

    subroutine get_rho_derivatives(neo_h_right, neo_h_left, neo_phi_right, neo_phi_left, dh_dvpa, dh_dmu, drho, dneo_h_drho, dneo_phi_drho)
        ! Parallelisation.
        use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

        ! Grids.
        use grids_z, only: nzgrid
        use grids_velocity, only: vpa, mu
        use grids_species, only: spec

        ! Geometry.
        use geometry, only: bmag, dBdrho

        implicit none

        real, intent(in) :: drho
        real, dimension(-nzgrid:, vmu_lo%llim_proc:), intent(in) :: neo_h_right, neo_h_left
        real, dimension(-nzgrid:, vmu_lo%llim_proc:), intent(in) :: dh_dvpa, dh_dmu
        real, dimension(-nzgrid:), intent(in) :: neo_phi_right, neo_phi_left

        real, dimension(-nzgrid:, vmu_lo%llim_proc:), intent(out) :: dneo_h_drho
        real, dimension(-nzgrid:), intent(out) :: dneo_phi_drho

        ! Local variables.
        integer :: ivmu, is, imu, iv, ia

        ! Assume one magnetic field line.
        ia = 1

        ! ================================================================================ !

        ! The neo_phi derivative is trivial to calculate via a central difference algorithm.
        dneo_phi_drho = (neo_phi_right - neo_phi_left) / (2.0 * drho)

        ! Central difference of neo_h: this is the derivative at fixed NEO variables (E, ξ).
        dneo_h_drho = (neo_h_right - neo_h_left) / (2.0 * drho)

        ! Transform to fixed stella variables (v∥, μ). 
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            dneo_h_drho(:, ivmu) = dneo_h_drho(:, ivmu) &
                + mu(imu) * ( dBdrho(:) / bmag(ia, :) + spec(is)%tprim ) * dh_dmu(:, ivmu) &
                + 0.5 * spec(is)%tprim * vpa(iv) * dh_dvpa(:, ivmu)
        end do

    end subroutine get_rho_derivatives

! ================================================================================================================================================================================= !
! ---------------------------------- Calculate the zeroeth order moment of H_1. This is used in the field equations and also for testing purposes. -------------------------------- !
! ================================================================================================================================================================================= !
!
! Inputs are on a single field line (ia = 1): neo_h_global (-nzgrid:nzgrid, nvpa, nmu, nspec), neo_phi (-nzgrid:nzgrid).
!
! ================================================================================================================================================================================= !

    subroutine get_neo_moment(neo_h_global, neo_phi, moment_type, neo_mom)
        ! Grids.
        use grids_z, only: nzgrid
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac, nvpa, vpa, nmu
        use grids_species, only: spec, nspec

        ! Geometry (flux-surface average weights and |B|).
        use geometry, only: dl_over_b, bmag, geo_surf, drhodpsip

        ! Calculations.
        use calculations_velocity_integrals, only: integrate_vmu

        ! Output.
        use mp, only: proc0

        implicit none

        real, dimension(-nzgrid:, :, :, :), intent(in)    :: neo_h_global
        real, dimension(-nzgrid:)                         :: neo_phi
        character(len=*), intent(in)                      :: moment_type
        real, dimension(-nzgrid:, :), intent(out)         :: neo_mom
        
        ! Local variables.
        integer           :: is, imu, iv, iz, ia
        real, allocatable :: tmp(:, :, :, :)

        allocate(tmp(-nzgrid:nzgrid, nvpa, nmu, nspec))
        tmp = 0.0
        neo_mom = 0.0

        ia = 1

        ! Calculate the physical F_1.
        do is = 1, nspec
            do iv = 1, nvpa
                do imu = 1, nmu
                    do iz = -nzgrid, nzgrid
                        select case (trim(moment_type))

                        case ("zeroeth")
                            tmp(iz, iv, imu, is) = ( neo_h_global(iz, iv, imu, is) &
                            - spec(is)%zt * neo_phi(iz) ) &
                            * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)

                        case ("first")
                            tmp(iz, iv, imu, is) = vpa(iv) * neo_h_global(iz, iv, imu, is) &
                            * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                        end select
                    end do
                end do
            end do
        end do

        ! Integrate over velocity space.
        do is = 1, nspec
            do iz = -nzgrid, nzgrid
                call integrate_vmu(tmp(iz, :, :, is), iz, neo_mom(iz, is))
            end do
        end do

        ! Deallocate the temporary array.
        deallocate(tmp)
    end subroutine get_neo_moment


! ================================================================================================================================================================================= !
! ----------------------------- Calculate the H_1/F_1 velocity derivatives via finite difference methods on the stella grids, for testing purposes. ------------------------------- !
! ================================================================================================================================================================================= !
!
! Inputs are on a single field line (ia = 1): distributions (-nzgrid:nzgrid, nvpa, nmu, nspec), neo_phi (-nzgrid:nzgrid).
!
! ================================================================================================================================================================================= !

    subroutine get_neo_h_velocity_derivative_moment(neo_h_deriv_global, neo_h_global, neo_phi, vel_type, neo_mom)
        ! Grids.
        use grids_z, only: nzgrid
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac, nvpa, vpa, nmu, mu
        use grids_species, only: spec, nspec

        ! Geometry. 
        use geometry, only: bmag

        ! Calculations. 
        use calculations_velocity_integrals, only: integrate_vmu

        ! Constants.
        use constants, only: pi
 
        implicit none

        real, intent(in)                :: neo_h_deriv_global(-nzgrid:, :, :, :)
        real, intent(in)                :: neo_h_global(-nzgrid:, :, :, :)
        real, intent(in)                :: neo_phi(-nzgrid:)
        character(len=*), intent(in)    :: vel_type
        real, intent(out)               :: neo_mom(-nzgrid:, :)         

        ! Local variables.
        integer           :: is, imu, iv, iz, ia
        real, allocatable :: tmp(:, :, :, :)

        ia = 1

        allocate(tmp(-nzgrid:nzgrid, nvpa, nmu, nspec))

        tmp = 0.0
        neo_mom = 0.0        
 
        do is = 1, nspec
            do iv = 1, nvpa
                do imu = 1, nmu
                    do iz = -nzgrid, nzgrid           
                        select case (trim(vel_type))
                        case ("mu")
                            tmp(iz, iv, imu, is) = mu(imu) * ( neo_h_deriv_global(iz, iv, imu, is) &
                            - 2 * bmag(ia, iz) * ( neo_h_global(iz, iv, imu, is) - spec(is)%zt * neo_phi(iz) ) ) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) &
                            * maxwell_fac(is)

                        case ("vpa")
                            tmp(iz, iv, imu, is) = vpa(iv) * ( neo_h_deriv_global(iz, iv, imu, is) &
                            - 2 * vpa(iv) * ( neo_h_global(iz, iv, imu, is) - spec(is)%zt * neo_phi(iz) ) ) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) &
                            * maxwell_fac(is)
              
                        end select
                    end do
                end do
            end do
        end do

        ! Integrate over velocity space.
        do is = 1, nspec
            do iz = -nzgrid, nzgrid
                call integrate_vmu(tmp(iz, :, :, is), iz, neo_mom(iz, is))
            end do
        end do

        ! Deallocate the temporary array.
        deallocate(tmp)
    end subroutine get_neo_h_velocity_derivative_moment


! ================================================================================================================================================================================= !
! ---------------------------------- Convert NEO transport calculations in stella into NEO units. --------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine get_stella_to_neo_units(conv, neo_equil, neo_species, nspec_neo, surface_label)
        ! Grids.
        use grids_species, only: spec, nspec

        ! Geometry.
        use geometry, only: geo_surf, drhodpsip

        ! NEO data types.
        use NEO_interface, only: neo_equil_data, neo_species_data

        implicit none

        type(stella_to_neo_units_data), intent(out) :: conv
        type(neo_equil_data), intent(in)            :: neo_equil
        type(neo_species_data), intent(in)          :: neo_species
        integer, intent(in)                         :: nspec_neo
        character(len=*), intent(in)                :: surface_label

        ! Local variables.
        real, dimension(:), allocatable :: rn, rT, rm
        integer :: is, ns
        real, parameter :: tol = 1.0e-6

        ! Missing out.neo.equil (or out.neo.species): keep the default ratios (n = T = m = B = 1, vth = sqrt(2)).
        if (.not. (allocated(neo_equil%density) .and. allocated(neo_equil%temperature) .and. allocated(neo_species%mass))) then
            write (check_unit, '(a)') ' NOTICE: no NEO profiles for the '//trim(surface_label)//' surface; using default stella -> NEO ratios.'
            n_notices = n_notices + 1
            return
        end if

        ns = min(nspec, nspec_neo)
        if (nspec /= nspec_neo) then
            write (check_unit, '(a, i3, a, i3)') ' NOTICE: stella has nspec = ', nspec, ' but NEO has n_species = ', nspec_neo
            n_notices = n_notices + 1
        end if

        ! Allocate temporaray arrays. 
        allocate (rn(ns), rT(ns), rm(ns))

        ! Ratios (NEO profile on this surface) / (stella profile on the central surface).
        ! Central surface: n = n_ref/n_norm, T = T_ref/T_norm, m = m_ref/m_norm.
        ! Left/right surfaces: n = (n_ref/n_norm), n_s(surface)/n_s(central);
        ! T = (T_ref/T_norm) T_s(surface)/T_s(central),
        ! Masses are obviously the same on all surfaces.
        do is = 1, ns
            rn(is) = neo_equil%density(is, 1) / spec(is)%dens
            rT(is) = neo_equil%temperature(is, 1) / spec(is)%temp
            rm(is) = neo_species%mass(is) / spec(is)%mass
        end do

        conv%n = rn(1); conv%T = rT(1); conv%m = rm(1)
        if (any(abs(rn - conv%n) > tol * abs(conv%n)) .or. any(abs(rT - conv%T) > tol * abs(conv%T)) &
            .or. any(abs(rm - conv%m) > tol * abs(conv%m))) then
            write (check_unit, '(a)') ' NOTICE: stella and NEO species ratios differ between species on the '//trim(surface_label)// &
            ' surface (check species ordering, or unequal gradients across species).'
            n_notices = n_notices + 1
        end if

        conv%vth = sqrt(2.0 * conv%T / conv%m)                      ! v_th,ref / v_norm
        conv%B   = abs(geo_surf%rhoc * drhodpsip / geo_surf%qinp)   ! B_ref / B_unit

        write (check_unit, '(a, a, a, 5es13.5)') ' stella -> NEO units (', trim(surface_label), ' surface): n, T, m, v_th, B ratios = ', &
        conv%n, conv%T, conv%m, conv%vth, conv%B

        ! Deallocate temporaray arrays. 
        deallocate (rn, rT, rm)

    end subroutine get_stella_to_neo_units


! ================================================================================================================================================================================= !
! ------------------------------------------- Compare a quantity recomputed on stella grids (in NEO units) with the raw value from out.neo.transport. ----------------------------- !
! ------------------------------------------------------- Here we compute all of the quantities found in out.neo.transport on each surface ---------------------------------------- !
! ================================================================================================================================================================================= !

    ! Classification of each stella/NEO comparison (one diagnostic):
    !   SUCCESS : within the tolerance.
    !   WARNING : outside the tolerance, but the quantity is a small residual of much larger, cancelling contributions (value_abs given, |value| / value_abs <
    !             cancellation_threshold) and the difference relative to value_abs is within the tolerance, i.e. limited by stella's velocity/z resolution.
    !   FAILED  : outside the tolerance otherwise.
    ! value_abs is the same calculation with |.| taken everywhere (no cancellation); pass it for the radial fluxes.
    subroutine compare_with_neo(label, stella_value, neo_value, rel_tol, abs_tol, value_abs)

        implicit none

        character(len=*), intent(in) :: label
        real, intent(in)             :: stella_value, neo_value, rel_tol
        ! Optional: when |neo_value| < abs_tol (e.g. a quantity that is zero in this run), compare absolute differences instead.
        real, intent(in), optional   :: abs_tol
        real, intent(in), optional   :: value_abs

        real :: rel_diff, abs_diff, ratio, diff_over_abs
        logical :: cancelling
        real, parameter :: tiny_value = 1.0e-30
        real, parameter :: cancellation_threshold = 1.0e-3

        if (present(abs_tol)) then
            if (abs(neo_value) < abs_tol) then
                abs_diff = abs(stella_value - neo_value)
                if (abs_diff <= abs_tol) then
                    write (check_unit, '(a, a, a, es14.6, a, es14.6, a, es10.3, a)') ' SUCCESS: ', label, ':  '//trim(check_code_label)//' = ', stella_value, &
                        '   NEO = ', neo_value, '   (absolute difference ', abs_diff, ')'
                    call record_check(check_pass)
                else
                    write (check_unit, '(a, a, a, es14.6, a, es14.6, a, es10.3, a, es10.3, a)') ' FAILED:  ', label, ':  '//trim(check_code_label)//' = ', stella_value, &
                        '   NEO = ', neo_value, '   (absolute difference ', abs_diff, ', tolerance ', abs_tol, ')'
                    call record_check(check_fail)
                end if
                return
            end if
        end if

        rel_diff = abs(stella_value - neo_value) / max(abs(neo_value), tiny_value)

        ! Degree of cancellation, if the no-cancellation value is available.
        cancelling = .false.
        if (present(value_abs)) then
            if (value_abs > 0.0) then
                ratio = abs(stella_value) / value_abs
                diff_over_abs = abs(stella_value - neo_value) / value_abs
                cancelling = ratio < cancellation_threshold
            end if
        end if

        if (rel_diff <= rel_tol) then
            write (check_unit, '(a, a, a, es14.6, a, es14.6, a, f8.3, a)') ' SUCCESS: ', label, ':  '//trim(check_code_label)//' = ', stella_value, &
                '   NEO = ', neo_value, '   (', 100.0 * rel_diff, ' % difference)'
            call record_check(check_pass)
        else if (cancelling .and. diff_over_abs <= rel_tol) then
            write (check_unit, '(a, a, a, es14.6, a, es14.6, a, f8.3, a, f6.2, a)') ' WARNING: ', label, ':  '//trim(check_code_label)//' = ', stella_value, &
                '   NEO = ', neo_value, '   (', 100.0 * rel_diff, ' % difference, tolerance ', 100.0 * rel_tol, ' %)'
            write (check_unit, '(a, es10.3, a, es10.3, a)') '          strongly cancelling: |value| / |value|_no-cancellation = ', ratio, &
                ', |stella - NEO| / |value|_no-cancellation = ', diff_over_abs, ';'
            write (check_unit, '(a)') '          the relative difference is limited by velocity/z resolution (increase nvgrid, nmu to test).'
            call record_check(check_warn)
        else
            write (check_unit, '(a, a, a, es14.6, a, es14.6, a, f8.3, a, f6.2, a)') ' FAILED:  ', label, ':  '//trim(check_code_label)//' = ', stella_value, &
                '   NEO = ', neo_value, '   (', 100.0 * rel_diff, ' % difference, tolerance ', 100.0 * rel_tol, ' %)'
            if (cancelling) write (check_unit, '(a, es10.3, a, es10.3, a)') '          strongly cancelling: |value| / |value|_no-cancellation = ', &
                ratio, ', but |stella - NEO| / |value|_no-cancellation = ', diff_over_abs, ' is also outside the tolerance.'
            call record_check(check_fail)
        end if

    end subroutine compare_with_neo


    subroutine record_check(status)

        implicit none

        integer, intent(in) :: status

        n_checks = n_checks + 1
        select case (status)
        case (check_pass)
            n_pass = n_pass + 1
        case (check_warn)
            n_warn = n_warn + 1
        case default
            n_fail = n_fail + 1
        end select

    end subroutine record_check


    ! Open <run_name>//ext (proc0 only) and reset the counters. Defaults (NEO run): ext = '.neo_transport_check.dat', title = 'NEO transport check',
    ! code_label = 'stella' (the name printed for the value computed from the data on stella's grids). The SFINCS interface passes its own.
    subroutine open_neo_check_file(ext, title, code_label, description)
        use file_utils, only: open_output_file

        implicit none

        character(len=*), intent(in), optional :: ext, title, code_label, description

        check_ext = neo_check_ext
        check_title = 'NEO transport check'
        check_code_label = 'stella'
        if (present(ext)) check_ext = ext
        if (present(title)) check_title = title
        if (present(code_label)) check_code_label = code_label

        call open_output_file(check_unit, trim(check_ext))
        n_checks = 0; n_pass = 0; n_warn = 0; n_fail = 0; n_notices = 0

        write (check_unit, '(a)') '############################################################################################'
        if (present(description)) then
            write (check_unit, '(a)') ' '//trim(description)
        else
            write (check_unit, '(a)') ' Verification of the NEO data on stella grids against out.neo.transport (and .left, .right)'
        end if
        write (check_unit, '(a, f6.2, a)') ' Relative tolerance: ', 100.0 * neo_check_rel_tol, ' %'
        write (check_unit, '(a)') ' SUCCESS = within tolerance; WARNING = outside tolerance but limited by resolution in a strongly'
        write (check_unit, '(a)') ' cancelling quantity; FAILED = outside tolerance; NOTICE = setup issue (not a diagnostic).'
        write (check_unit, '(a)') '############################################################################################'

    end subroutine open_neo_check_file


    ! Write a setup notice (not a diagnostic) to the check file and count it in the summary.
    subroutine add_neo_check_notice(message)

        implicit none

        character(len=*), intent(in) :: message

        write (check_unit, '(a)') ' NOTICE: '//trim(message)
        n_notices = n_notices + 1

    end subroutine add_neo_check_notice


    ! Write the totals to the file, close it, and print the one-line summary to the screen.
    subroutine close_neo_check_file_and_summarise
        use iso_fortran_env, only: output_unit
        use file_utils, only: close_output_file, run_name

        implicit none

        character(len=512) :: summary
        character(len=32)  :: frac_pass, frac_warn, frac_fail, tol

        write (frac_pass, '(i0, a, i0)') n_pass, '/', n_checks
        write (frac_warn, '(i0, a, i0)') n_warn, '/', n_checks
        write (frac_fail, '(i0, a, i0)') n_fail, '/', n_checks
        write (tol, '(f0.2, a)') 100.0 * neo_check_rel_tol, '%'
        if (tol(1:1) == '.') tol = '0'//trim(tol)

        summary = ' '//trim(check_title)//': '//trim(frac_pass)//' diagnostics pass to within '//trim(tol)//', ' &
                  //trim(frac_warn)//' warnings, '//trim(frac_fail)//' FAILED.'
        if (n_notices > 0) then
            write (tol, '(i0)') n_notices
            summary = trim(summary)//' '//trim(tol)//' setup notice(s).'
        end if

        write (check_unit, '(a)') ''
        write (check_unit, '(a)') '############################################################################################'
        write (check_unit, '(a)') trim(summary)
        write (check_unit, '(a)') '############################################################################################'
        call close_output_file(check_unit)
        check_unit = output_unit

        write (output_unit, '(a)') trim(summary)//' See '//trim(run_name)//trim(check_ext)//' for more details.'

    end subroutine close_neo_check_file_and_summarise


! ================================================================================================================================================================================= !
! ----------------------------------------------- Calculations of out.neo.transport quantities from NEO data on stella grids, in NEO units. --------------------------------------- !
! ================================================================================================================================================================================= !
!
! All on a single field line (ia = 1), using stella's flux-surface average.
!
! ================================================================================================================================================================================= !

    ! Calculate the flux surface average of a zed-dependent quantity.
    real function calculate_fsa(a)
        ! Grids. 
        use grids_z, only: nzgrid

        ! Geometry.
        use geometry, only: dl_over_b

        implicit none

        real, dimension(-nzgrid:), intent(in) :: a
        
        ! Local variables. 
        integer, parameter :: ia = 1

        calculate_fsa = sum(dl_over_b(ia, :) * a) / sum(dl_over_b(ia, :))

    end function calculate_fsa

    ! ================================================================================= !
    ! Calculate the square of the normalised NEO potential: 
    ! < (e Φ_1 / T_norm)^2 >, with neo_phi = e Φ_1 / T_ref.
    real function calculate_phi_sq(neo_phi, conv)
        ! Grids. 
        use grids_z, only: nzgrid

        implicit none

        real, dimension(-nzgrid:), intent(in)      :: neo_phi
        type(stella_to_neo_units_data), intent(in) :: conv

        calculate_phi_sq = calculate_fsa((conv%T * neo_phi)**2)

    end function calculate_phi_sq


    ! ================================================================================= !
    ! Calculate the speices zeroth order parallel flow:  
    ! < u_par,s B > / (v_norm B_unit) = 
    ! < (u_par,s / v_th,s) (B / B_ref) > * stm_s * (v_th,ref / v_norm) * (B_ref / B_unit).
    subroutine calculate_upar_B(neo_h_global, neo_phi, conv, upar_B)
        ! Grids. 
        use grids_z, only: nzgrid
        use grids_species, only: spec, nspec
        
        ! Geometry.
        use geometry, only: bmag

        implicit none

        real, dimension(-nzgrid:, :, :, :), intent(in) :: neo_h_global
        real, dimension(-nzgrid:), intent(in)          :: neo_phi
        type(stella_to_neo_units_data), intent(in)     :: conv
        real, dimension(:), intent(out)                :: upar_B

        ! Local variables. 
        real, dimension(:, :), allocatable :: upar
        real, dimension(:), allocatable    :: phi_copy
        integer :: is
        integer, parameter :: ia = 1

        ! Allocate temporary arrays. 
        allocate(upar(-nzgrid:nzgrid, nspec), phi_copy(-nzgrid:nzgrid))
        
        ! Create a copy of phi. 
        phi_copy = neo_phi

        ! Calculate the first moment.
        call get_neo_moment(neo_h_global, phi_copy, "first", upar)

        ! Calculate the parallel flow in NEO units. 
        do is = 1, nspec
            upar_B(is) = calculate_fsa(upar(:, is) * bmag(ia, :)) * spec(is)%stm * conv%vth * conv%B
        end do

        ! Deallocate temporary arrays. 
        deallocate(upar, phi_copy)

    end subroutine calculate_upar_B


    ! ================================================================================= !
    ! Calculate the first order bootstrap current:
    ! < j_par B > / (e n_norm v_norm B_unit) = 
    ! sum_s Z_s (n_s / n_norm) < u_par,s B > / (v_norm B_unit).
    real function calculate_jpar(upar_B, conv)
        ! Grids.
        use grids_species, only: spec, nspec

        implicit none

        real, dimension(:), intent(in)             :: upar_B
        type(stella_to_neo_units_data), intent(in) :: conv
        
        ! Local variables.
        integer :: is

        calculate_jpar = 0.0

        ! Calculate the bootstrap current in NEO units.
        do is = 1, nspec
            calculate_jpar = calculate_jpar + spec(is)%z * spec(is)%dens * conv%n * upar_B(is)
        end do

    end function calculate_jpar


    ! ================================================================================= !
    ! Calculate the species first-order flow coefficients:
    ! Dimensional: 
    ! K_s = [ < u_par,s B > - (I/ψ') ρ_* (T_s / Z_s) (dlnndr_s - (Z_s / T_s) dΦ_0/dr + dlntdr_s) ] &
    ! / < B² / n_s >  
    ! Dimensionless: 
    ! k_s = - K_s < B² / n_s > / ( dlntdr_s (I/ψ') ρ_* T_s / Z_s )                                    
    subroutine calculate_flow_coefficients(upar_B, neo_equil, conv, conv_local, fprim, tprim, K_big, k_small)
        ! Grids. 
        use grids_z, only: nzgrid
        use grids_species, only: spec, nspec
        
        ! Geometry.
        use geometry, only: bmag, geo_surf, drhodpsip
        
        ! NEO equilibrium data.  
        use NEO_interface, only: neo_equil_data

        implicit none

        real, dimension(:), intent(in)             :: upar_B, fprim, tprim
        type(neo_equil_data), intent(in)           :: neo_equil
        type(stella_to_neo_units_data), intent(in) :: conv, conv_local
        real, dimension(:), intent(out)            :: K_big, k_small

        ! Local variables. 
        real    :: I_div_psip, B2_avg, B2_div_dens, n_s, T_s, diamagnetic_fac
        integer :: is
        integer, parameter :: ia = 1

        ! Calculate necessary geometry factors. 
        I_div_psip = geo_surf%rgeo * drhodpsip
        B2_avg = calculate_fsa((bmag(ia, :) * conv%B)**2)

        ! Calculate the flow coeffecients. 
        do is = 1, nspec
            n_s = spec(is)%dens * conv_local%n
            T_s = spec(is)%temp * conv_local%T
            B2_div_dens = B2_avg / n_s

            diamagnetic_fac = I_div_psip * neo_equil%rho_star(1) * T_s / spec(is)%z

            K_big(is) = (upar_B(is) - diamagnetic_fac * (fprim(is) + tprim(is))) / B2_div_dens

            if (abs(tprim(is)) > epsilon(0.0)) then
                k_small(is) = -K_big(is) * B2_div_dens / (tprim(is) * diamagnetic_fac)
            else
                k_small(is) = 0.0
            end if
        end do

    end subroutine calculate_flow_coefficients


    ! ================================================================================= !
    ! Calculate the species first-order poloidal flow at the outboard midplane:
    subroutine calculate_vpol_th0(K_big, conv, conv_local, vpol_th0)
        ! Grids. 
        use grids_z, only: nzgrid, zed 
        use grids_species, only: spec, nspec

        ! Geometry. 
        use geometry, only: bpol2

        implicit none

        real, dimension(:), intent(in)             :: K_big
        type(stella_to_neo_units_data), intent(in) :: conv, conv_local
        real, dimension(:), intent(out)            :: vpol_th0

        ! Local variables. 
        real    :: Bp_th0
        integer :: is, iz0
        integer, parameter :: ia = 1

        ! Index of the outboard midplane, z = θ = 0.
        iz0 = minloc(abs(zed), dim = 1) - nzgrid - 1

        ! Calculate geometry coeffecients. 
        Bp_th0 = sqrt(max(bpol2(ia, iz0), 0.0)) * conv%B

        ! Calculate species first order toroidal flows at the outboard midplane. 
        do is = 1, nspec
            vpol_th0(is) = K_big(is) * Bp_th0 / (spec(is)%dens * conv_local%n)
        end do

    end subroutine calculate_vpol_th0


    ! ================================================================================= !
    ! Calculate the species first-order toroidal flow at the outboard midplane:
    subroutine calculate_vtor_th0(K_big, neo_equil, conv, conv_local, fprim, tprim, vtor_th0)
        ! Grids. 
        use grids_z, only: nzgrid, zed
        use grids_species, only: spec, nspec
        
        ! Geoemetry. 
        use geometry, only: btor, geo_surf, drhodpsip
        
        ! Neo equilibrium data. 
        use NEO_interface, only: neo_equil_data

        implicit none

        real, dimension(:), intent(in)             :: K_big, fprim, tprim
        type(neo_equil_data), intent(in)           :: neo_equil
        type(stella_to_neo_units_data), intent(in) :: conv, conv_local
        real, dimension(:), intent(out)            :: vtor_th0

        ! Local variables. 
        real    :: Bt_th0, I_div_psip, n_s, T_s
        integer :: is, iz0

        ! Index of the outboard midplane, z = θ = 0.
        iz0 = minloc(abs(zed), dim = 1) - nzgrid - 1

        ! Calculate necessary geometry factors.
        Bt_th0 = btor(iz0) * conv%B
        I_div_psip = geo_surf%rgeo * drhodpsip

        ! Calculate the species first order toroidal flows at the outboard midplane.
        do is = 1, nspec
            n_s = spec(is)%dens * conv_local%n
            T_s = spec(is)%temp * conv_local%T
            vtor_th0(is) = K_big(is) * Bt_th0 / n_s &
            + I_div_psip * neo_equil%rho_star(1) * T_s / spec(is)%z * (fprim(is) + tprim(is)) / Bt_th0
        end do

    end subroutine calculate_vtor_th0


    ! ================================================================================= !
    ! Calculate the species zeroeth-order toroidal flow at the outboard midplane:
    ! v_φ^(0)(θ = 0) / v_norm.
    ! This is zero unless NEO was run with rotation_model = 2.
    real function calculate_vtor0_th0(neo_equil)
        ! Grids. 
        use grids_z, only: nzgrid, zed
        
        ! Geometry. 
        use geometry, only: rmajor
        
        ! NEO equilibrium data. 
        use NEO_interface, only: neo_equil_data

        implicit none

        type(neo_equil_data), intent(in) :: neo_equil
        
        ! Local variables. 
        integer :: iz0

        ! Index of the outboard midplane, z = θ = 0.
        iz0 = minloc(abs(zed), dim = 1) - nzgrid - 1

        ! Calculate the species zeroeth-order toroidal flow at the outboard midplane.
        calculate_vtor0_th0 = neo_equil%angular_frequency(1) * rmajor(iz0)

    end function calculate_vtor0_th0


    ! ================================================================================= !
    ! Calculate the zeroth-order parallel flow of rigid toroidal rotation.
    ! < u_par^(0) B > / (v_norm B_unit): 
    ! u_par^(0) = ω_0 R B_t / B, so < u_par^(0) B > = ω_0 < R B_t >.
    real function calculate_uparB0(neo_equil, conv)
        ! Grids.
        use grids_z, only: nzgrid
        
        ! Geometry. 
        use geometry, only: rmajor, btor
        
        ! NEO equilibrium data. 
        use NEO_interface, only: neo_equil_data

        implicit none

        type(neo_equil_data), intent(in)           :: neo_equil
        type(stella_to_neo_units_data), intent(in) :: conv

        ! Calculate the zeroth-order parallel flow of rigid toroidal rotation.
        calculate_uparB0 = neo_equil%angular_frequency(1) * calculate_fsa(rmajor * btor) * conv%B

    end function calculate_uparB0

    ! ================================================================================= !
    ! A generalised second order species radial flux calculation; used for particle
    ! heat and momentum fluxes. 
    ! flux_s = < ∫ d³v F_0 H_1 E^p v_d,s · ∇r >,    
    ! v_d,s · ∇r = (m_s / (Z_s e B)) (v∥² + v⊥²/2) (b × ∇B · ∇r) / B.
    ! An optional extra weight w(z, v∥) multiplies the integrand (used for the momentum flux).
    subroutine calculate_radial_flux(neo_h_global, neo_equil, conv, conv_local, energy_power, flux, flux_abs, weight)
        ! Grids. 
        use grids_z, only: nzgrid
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac, nvpa, vpa, nmu, mu
        use grids_species, only: spec, nspec
        
        ! Geometry. 
        use geometry, only: bmag, B_times_gradB_dot_gradx, dxdpsi, drhodpsip
        
        ! Calculations.
        use calculations_velocity_integrals, only: integrate_vmu
        
        ! NEO equilibrium data. 
        use NEO_interface, only: neo_equil_data

        implicit none

        real, dimension(-nzgrid:, :, :, :), intent(in) :: neo_h_global
        type(neo_equil_data), intent(in)               :: neo_equil
        ! conv: central ratios (B_ref / B_unit); conv_local: this surface's ratios (local n_s, T_s).
        type(stella_to_neo_units_data), intent(in)     :: conv, conv_local
        integer, intent(in)                            :: energy_power
        real, dimension(:), intent(out)                :: flux, flux_abs
        real, dimension(-nzgrid:, :), intent(in), optional :: weight   ! (z, v∥)


        ! Local variables. 
        real, dimension(:, :), allocatable :: tmp, tmp_abs
        real, dimension(:), allocatable    :: drift_moment, drift_moment_abs
        real    :: D, prefac, energy
        integer :: is, iz, iv, imu
        integer, parameter :: ia = 1

        ! Allocate temporary arrays. 
        allocate(tmp(nvpa, nmu), tmp_abs(nvpa, nmu))
        allocate(drift_moment(-nzgrid:nzgrid), drift_moment_abs(-nzgrid:nzgrid))

        ! Calculate necessary NEO factors. 
        D = dxdpsi / drhodpsip

        ! Calculate species second order flux. 
        do is = 1, nspec
            do iz = -nzgrid, nzgrid
                do imu = 1, nmu
                    do iv = 1, nvpa
                        ! Calculate energy here for convenience. 
                        energy = vpa(iv)**2 + 2.0 * mu(imu) * bmag(ia, iz)

                        ! Calculate appropriate integrand. 
                        tmp(iv, imu) = energy**energy_power * (vpa(iv)**2 + mu(imu) * bmag(ia, iz)) * neo_h_global(iz, iv, imu, is) &
                        * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                        if (present(weight)) tmp(iv, imu) = tmp(iv, imu) * weight(iz, iv)
                        tmp_abs(iv, imu) = abs(tmp(iv, imu))
                    end do
                end do
                ! Integrate over velocity space. 
                call integrate_vmu(tmp, iz, drift_moment(iz))
                call integrate_vmu(tmp_abs, iz, drift_moment_abs(iz))
            end do

            ! (n_s / n_norm) * 2 (T_s / T_norm) / Z_s * ρ_* / (B_ref / B_unit) / D.
            prefac = (spec(is)%dens * conv_local%n) * 2.0 * (spec(is)%temp * conv_local%T) &
            / spec(is)%z * neo_equil%rho_star(1) / conv%B / D

            ! Finally get the fluxes, including the flux surface average. 
            flux(is)     = prefac * calculate_fsa(drift_moment * B_times_gradB_dot_gradx(ia, :))
            flux_abs(is) = abs(prefac) * calculate_fsa(abs(drift_moment_abs * B_times_gradB_dot_gradx(ia, :)))
        end do

        ! Deallocate temporary arrays. 
        deallocate(tmp, tmp_abs, drift_moment, drift_moment_abs)

    end subroutine calculate_radial_flux


    ! ================================================================================= !
    ! Calculate the second order radial particle flux for each species. 
    ! Γ_s / (n_norm v_norm), and the same with no cancellation.
    subroutine calculate_pflux(neo_h_global, neo_equil, conv, conv_local, pflux, pflux_abs)
        ! Grids. 
        use grids_z, only: nzgrid
        
        ! NEO equilibrium data. 
        use NEO_interface, only: neo_equil_data

        implicit none

        real, dimension(-nzgrid:, :, :, :), intent(in) :: neo_h_global
        type(neo_equil_data), intent(in)               :: neo_equil
        type(stella_to_neo_units_data), intent(in)     :: conv, conv_local
        real, dimension(:), intent(out)                :: pflux, pflux_abs

        call calculate_radial_flux(neo_h_global, neo_equil, conv, conv_local, 0, pflux, pflux_abs)

    end subroutine calculate_pflux

    ! ================================================================================= !
    ! Calculate the second order radial energy/heat flux for each species.
    ! Q_s / (n_norm v_norm T_norm), and the same with no cancellation. 
    ! NEO's eflux column is Q_s + ω_0 Π_s.
    subroutine calculate_eflux(neo_h_global, neo_equil, conv, conv_local, eflux, eflux_abs)
        ! Grids. 
        use grids_z, only: nzgrid
        use grids_species, only: spec, nspec
        
        ! NEO equilibriumd data. 
        use NEO_interface, only: neo_equil_data

        implicit none

        real, dimension(-nzgrid:, :, :, :), intent(in) :: neo_h_global
        type(neo_equil_data), intent(in)               :: neo_equil
        type(stella_to_neo_units_data), intent(in)     :: conv, conv_local
        real, dimension(:), intent(out)                :: eflux, eflux_abs

        ! Local variables. 
        integer :: is

        call calculate_radial_flux(neo_h_global, neo_equil, conv, conv_local, 1, eflux, eflux_abs)

        ! m v² / 2 = T_s E, with T_s / T_norm on this surface.
        do is = 1, nspec
            eflux(is)     = eflux(is)     * spec(is)%temp * conv_local%T
            eflux_abs(is) = eflux_abs(is) * spec(is)%temp * conv_local%T
        end do

    end subroutine calculate_eflux


    ! ================================================================================= !
    ! Calculate the second order radial momentum flux for each species.
    ! Π_s / (n_norm T_norm a_norm), and the same with no cancellation. 
    ! For an up-down symmetric equilibrium without rotation: Π_s = 0. 
    ! WARNING: THIS IS NOT A GENERAL CALCULATION, IT IS NOT ACCURATE FOR UP-DOWN
    ! ASYMMETRIC EQUILIBRIA.
    subroutine calculate_mflux(neo_h_global, neo_equil, conv, conv_local, mflux, mflux_abs)
        ! Grids. 
        use grids_z, only: nzgrid
        use grids_velocity, only: vpa, nvpa
        use grids_species, only: spec, nspec
        
        ! Geomtery.
        use geometry, only: rmajor, btor, bmag
        
        ! NEO equilibrium data. 
        use NEO_interface, only: neo_equil_data

        implicit none

        real, dimension(-nzgrid:, :, :, :), intent(in) :: neo_h_global
        type(neo_equil_data), intent(in)               :: neo_equil
        type(stella_to_neo_units_data), intent(in)     :: conv, conv_local
        real, dimension(:), intent(out)                :: mflux, mflux_abs

        ! Local variables. 
        real, dimension(:, :), allocatable :: weight
        real    :: fac
        integer :: is, iz, iv
        integer, parameter :: ia = 1

        ! Allocate temporary arrays. 
        allocate(weight(-nzgrid:nzgrid, nvpa))
       
        do iv = 1, nvpa
            do iz = -nzgrid, nzgrid
                weight(iz, iv) = vpa(iv) * rmajor(iz) * btor(iz) / bmag(ia, iz)
            end do
        end do

        call calculate_radial_flux(neo_h_global, neo_equil, conv, conv_local, 0, mflux, mflux_abs, weight)

        ! m_s v_th,s v_norm / T_norm = sqrt(2 (m_s / m_norm)(T_s / T_norm)), with T_s on this surface.
        do is = 1, nspec
            fac = sqrt(2.0 * (spec(is)%mass * conv_local%m) * (spec(is)%temp * conv_local%T))
            mflux(is)     = mflux(is) * fac
            mflux_abs(is) = mflux_abs(is) * fac
        end do

        ! Deallocate temporary arrays. 
        deallocate(weight)

    end subroutine calculate_mflux


! ================================================================================================================================================================================= !
! ------------------------------------------- Check that a quantity that should vanish by symmetry is consistent with zero in both codes. ----------------------------------------- !
! ================================================================================================================================================================================= !
!
! Both |stella| / value_abs and |NEO| / value_abs must be below zero_tol, where value_abs is the stella calculation with no cancellation. A relative stella/NEO
! difference is meaningless here, since both values are numerical noise.
!
! ================================================================================================================================================================================= !

    subroutine compare_zero_with_neo(label, stella_value, neo_value, value_abs, reason)
        implicit none

        character(len=*), intent(in) :: label, reason
        real, intent(in)             :: stella_value, neo_value, value_abs

        ! Local variables. 
        real :: r_stella, r_neo
        real, parameter :: zero_tol = 1.0e-4

        r_stella = abs(stella_value) / max(value_abs, tiny(1.0))
        r_neo    = abs(neo_value) / max(value_abs, tiny(1.0))

        if (r_stella <= zero_tol .and. r_neo <= zero_tol) then
            write (check_unit, '(a, a, a, es14.6, a, es14.6, a, es10.3, a, es10.3, a)') ' SUCCESS: ', label, ':  '//trim(check_code_label)//' = ', stella_value, &
            '   NEO = ', neo_value, '   (both zero: |stella|, |NEO| / |value|_no-cancellation = ', r_stella, ', ', r_neo, ')'
            call record_check(check_pass)
        else
            write (check_unit, '(a, a, a, es14.6, a, es14.6, a, es10.3, a, es10.3, a)') ' FAILED:  ', label, ':  '//trim(check_code_label)//' = ', stella_value, &
            '   NEO = ', neo_value, '   (expected zero: |stella|, |NEO| / |value|_no-cancellation = ', r_stella, ', ', r_neo, ')'
            write (check_unit, '(a, a)') '          expected zero because: ', reason
            call record_check(check_fail)
        end if

    end subroutine compare_zero_with_neo


! ================================================================================================================================================================================= !
! ------------------------------------------------------- Compute and compare all out.neo.transport quantities on one flux surface. ------------------------------------------------ !
! ================================================================================================================================================================================= !

    subroutine verify_neo_transport_on_surface(surface_label, neo_h_global, neo_phi, transport, equil, conv, conv_local)
        ! Grids. 
        use grids_z, only: nzgrid
        use grids_species, only: nspec, spec
        
        ! Neo equilibrium and transport data. 
        use NEO_interface, only: neo_transport_data, neo_equil_data

        implicit none

        character(len=*), intent(in)                   :: surface_label
        real, dimension(-nzgrid:, :, :, :), intent(in) :: neo_h_global
        real, dimension(-nzgrid:), intent(in)          :: neo_phi
        type(neo_transport_data), intent(in)           :: transport
        ! NEO profiles on this surface (for ω_0).
        type(neo_equil_data), intent(in)               :: equil
        ! conv: central-surface ratios (for quantities normalised to fixed NEO units, e.g. Φ_1).
        ! conv_local: this surface's ratios (for moments of H_1, which are in local n_s, v_th,s units).
        type(stella_to_neo_units_data), intent(in)     :: conv, conv_local

        ! Local variables. 
        real, dimension(nspec) :: upar_B, pflux, pflux_abs, eflux, eflux_abs, mflux, mflux_abs
        real, dimension(nspec) :: fprim_used, tprim_used, K_big, k_small, vpol_th0, vtor_th0
        character(len=128)     :: label
        integer                :: is
        real, parameter        :: rel_tol = neo_check_rel_tol
        ! Absolute tolerance for quantities that vanish in a non-rotating run.
        real, parameter        :: abs_tol = 1.0e-8

        write (check_unit, '(a)') ''
        write (check_unit, '(a)') ' Verifying NEO data on stella grids against out.neo.transport, '//trim(surface_label)//' surface:'

        if (.not. transport%found) then
            write (check_unit, '(a)') ' NOTICE: no transport data for the '//trim(surface_label)//' surface; skipping.'
            n_notices = n_notices + 1
            return
        end if

        ! Equilibrium potential correction. 
        call compare_with_neo('<(e Phi_1/T_norm)^2>', calculate_phi_sq(neo_phi, conv), transport%phi_sq_avg(1), rel_tol)

        ! Zeroth-order toroidal rotation.
        if (allocated(equil%angular_frequency)) then
            if (abs(equil%angular_frequency(1)) > 0.0) then
                write (check_unit, '(a, es14.6, a)') ' NOTICE: NEO was run with zeroth-order rotation, omega_0 = ', equil%angular_frequency(1), &
                ' (v_norm/a); stella assumes no zeroth-order rotation.'
                n_notices = n_notices + 1
            end if
            call compare_with_neo('v_tor^(0)(theta=0)/v_norm', calculate_vtor0_th0(equil), transport%vtor0_th0(1), rel_tol, abs_tol)
            call compare_with_neo('<u_par^(0) B>/(v_norm B_unit)', calculate_uparB0(equil, conv), transport%uparB0(1), rel_tol, abs_tol)

            ! Second-order radial particle flux.
            call calculate_pflux(neo_h_global, equil, conv, conv_local, pflux, pflux_abs)
            do is = 1, nspec
                write (label, '(a, i0)') 'Gamma/(n_norm v_norm), species ', is
                call compare_with_neo(trim(label), pflux(is), transport%pflux(is, 1), rel_tol, value_abs = pflux_abs(is))
            end do

            ! Second-order radial energy flux.
            call calculate_eflux(neo_h_global, equil, conv, conv_local, eflux, eflux_abs)
            do is = 1, nspec
                write (label, '(a, i0)') 'Q/(n_norm v_norm T_norm), species ', is
                call compare_with_neo(trim(label), eflux(is), transport%eflux(is, 1), rel_tol, value_abs = eflux_abs(is))
            end do

            ! Second-order radial momentum flux.
            call calculate_mflux(neo_h_global, equil, conv, conv_local, mflux, mflux_abs)
            do is = 1, nspec
                write (label, '(a, i0)') 'Pi/(n_norm T_norm a_norm), species ', is
                call compare_zero_with_neo(trim(label), mflux(is), transport%mflux(is, 1), mflux_abs(is), &
                'no zeroth-order rotation and an up-down symmetric equilibrium (stella omits NEO''s driftxrot3 term).')
            end do
        else
            write (check_unit, '(a)') ' NOTICE: no NEO profiles for the '//trim(surface_label)//' surface; skipping the zeroth-order flows and radial fluxes.'
            n_notices = n_notices + 1
        end if

        call calculate_upar_B(neo_h_global, neo_phi, conv_local, upar_B)
        call compare_with_neo('<j_par B>/(e n_norm v_norm B_unit)', calculate_jpar(upar_B, conv_local), transport%jpar(1), rel_tol)

        do is = 1, nspec
            write (label, '(a, i0)') '<u_par B>/(v_norm B_unit), species ', is
            call compare_with_neo(trim(label), upar_B(is), transport%uparB(is, 1), rel_tol)
        end do

        ! First-order flow coefficients.
        if (allocated(equil%rho_star)) then
            if (trim(surface_label) == 'central') then
                do is = 1, nspec
                    fprim_used(is) = spec(is)%fprim
                    tprim_used(is) = spec(is)%tprim
                    if (abs(spec(is)%fprim - equil%density_gradient(is, 1)) > 1.0e-4 * max(abs(spec(is)%fprim), 1.0) .or. &
                    abs(spec(is)%tprim - equil%temperature_gradient(is, 1)) > 1.0e-4 * max(abs(spec(is)%tprim), 1.0)) then
                        write (check_unit, '(a, i0, a, 2es13.5, a, 2es13.5)') ' NOTICE: species ', is, &
                        ' gradients differ: stella fprim, tprim = ', spec(is)%fprim, spec(is)%tprim, &
                        '   NEO dlnndr, dlntdr = ', equil%density_gradient(is, 1), equil%temperature_gradient(is, 1)
                        n_notices = n_notices + 1
                    end if
                end do
            else
                ! stella has no profile information off the central surface: use NEO's gradients there.
                fprim_used = equil%density_gradient(:, 1)
                tprim_used = equil%temperature_gradient(:, 1)
            end if

            call calculate_flow_coefficients(upar_B, equil, conv, conv_local, fprim_used, tprim_used, K_big, k_small)
            do is = 1, nspec
                write (label, '(a, i0)') 'k (dimensionless flow coefficient), species ', is
                call compare_with_neo(trim(label), k_small(is), transport%k_small(is, 1), rel_tol)
            end do
            do is = 1, nspec
                write (label, '(a, i0)') 'K/(n_norm v_norm/B_unit) (dimensional flow coefficient), species ', is
                call compare_with_neo(trim(label), K_big(is), transport%K_big(is, 1), rel_tol)
            end do

            ! First-order poloidal flow at the outboard midplane.
            call calculate_vpol_th0(K_big, conv, conv_local, vpol_th0)
            do is = 1, nspec
                write (label, '(a, i0)') 'v_pol(theta=0)/v_norm, species ', is
                call compare_with_neo(trim(label), vpol_th0(is), transport%vpol_th0(is, 1), rel_tol)
            end do

            ! First-order toroidal flow at the outboard midplane.
            call calculate_vtor_th0(K_big, equil, conv, conv_local, fprim_used, tprim_used, vtor_th0)
            do is = 1, nspec
                write (label, '(a, i0)') 'v_tor(theta=0)/v_norm, species ', is
                call compare_with_neo(trim(label), vtor_th0(is), transport%vtor_th0(is, 1), rel_tol)
            end do
        end if

    end subroutine verify_neo_transport_on_surface


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------ End Module. -------------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

end module neoclassical_terms_neo