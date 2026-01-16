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

    public :: initialised_neoclassical_terms_neo

    private

    logical :: include_neoclassical_terms                                           
    integer :: neo_option_switch                                                   ! Should be = 2 for NEO.
    integer :: nradii                                                              ! Can only be 3.
    real    :: drho                                                                ! Typically taken as 0.01.

    integer :: iz, unit
    character(len=128) :: filename

    real, dimension(:, :, :), allocatable :: neo_h, dneo_h_dpsi, dneo_h_dz  
    real, dimension(:, :), allocatable :: neo_phi, dneo_phi_dpsi, dneo_phi_dz
    real, dimension(:, :, :), allocatable :: dneo_h_dvpa, dneo_h_dmu
    
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
        use mp, only: proc0, broadcast
        use iso_fortran_env, only: output_unit
        use parallelisation_layouts, only: vmu_lo
        use grids_z, only: nzgrid
        use grids_velocity, only: nvpa, nmu 
        use grids_species, only: nspec
        use NEO_interface, only: read_basic_neo_files, read_neo_f_and_phi, neo_grid_data, neo_version_data        
        use neoclassical_diagnostics, only: write_dneo_h_dz_diagnostic, write_dneo_phi_dz_diagnostic

        implicit none

        real, dimension(:, :, :, :, :), allocatable :: neo_h_hat_in, neo_h_hat_right_in, neo_h_hat_left_in                ! Holds vectors for reconstructing NEO H_1 on 3 flux surfaces.
        real, dimension(:, :), allocatable :: neo_phi_in, neo_phi_right_in, neo_phi_left_in                               ! Holds NEO ϕ^1_0 on 3 flux surfaces.

        ! Intermediate arrays hold NEO h_hat data evaluated on the stella z grid, but on the NEO velocity grids.
        real, dimension(:, :, :, :, :), allocatable :: neo_h_hat_z_grid, neo_h_hat_right_z_grid, neo_h_hat_left_z_grid 
        real, dimension(:, :), allocatable :: neo_phi_right, neo_phi_left                                                 ! Holds NEO ϕ^1_0 data evaluated on the stella z grid.

        ! Holds NEO H_1 data evaluated on the stella z, v∥​ and μ grids. Since ϕ^1_0 is independent of velocity variables, there are no accompanying arrays for ϕ^1_0 here.
        real, dimension(:, :, :, :, :), allocatable :: neo_h_local, neo_h_local_right, neo_h_local_left

        ! Holds NEO H_1 data evaluated on the stella z, v∥​ and μ grids, compacted into 3 indicdes.
        real, dimension(:, :, :), allocatable :: neo_h_right, neo_h_left                 

        integer :: iz
        integer :: surface_index
        integer :: ierr      

        type(neo_grid_data) :: neo_grid
        type(neo_version_data) :: neo_version

        if (initialised_neoclassical_terms_neo) return ! Initialise only once.
        initialised_neoclassical_terms_neo = .true.

        if (proc0) then
            call read_basic_neo_files(neo_grid, neo_version)
            write(output_unit, '(A)') '! ======================================================================================================== !'
            write(output_unit, '("Reading neo files created on system ",A," at ",A," (commit : ",A,")")') neo_version%system, neo_version%date, neo_version%commit
            write(output_unit, '(A)') '! ======================================================================================================== !'
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
        if (.not. allocated(neo_phi)) allocate(neo_phi(-nzgrid:nzgrid, neo_grid%n_radial))
        if (.not. allocated(dneo_h_dpsi)) allocate(dneo_h_dpsi(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))
        if (.not. allocated(dneo_phi_dpsi)) allocate(dneo_phi_dpsi(-nzgrid:nzgrid, neo_grid%n_radial))
        if (.not. allocated(dneo_h_dz)) allocate(dneo_h_dz(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))
        if (.not. allocated(dneo_phi_dz)) allocate(dneo_phi_dz(-nzgrid:nzgrid, neo_grid%n_radial))
        if (.not. allocated(dneo_h_dvpa)) allocate(dneo_h_dvpa(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))
        if (.not. allocated(dneo_h_dmu)) allocate(dneo_h_dmu(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))

        ! Allocate all temporary arrays needed for initilization. 

        call allocate_temp_arrays

        if (proc0) then
            call read_neo_f_and_phi(neo_h_hat_in, neo_phi_in, neo_grid)                                ! Read in NEO h and ϕ^1_0 data for the central surface. 
            call read_neo_f_and_phi(neo_h_hat_right_in, neo_phi_right_in, neo_grid, suffix = '.right') ! Repeat for the right surface.  
            call read_neo_f_and_phi(neo_h_hat_left_in, neo_phi_left_in, neo_grid, suffix = '.left')    ! Repeat for the left surface.
        end if        
        
        ! Broadcast the read in data.
 
        call broadcast(neo_h_hat_in) ; call broadcast(neo_phi_in)             
        call broadcast(neo_h_hat_right_in) ; call broadcast(neo_phi_right_in)
        call broadcast(neo_h_hat_left_in) ; call broadcast(neo_phi_left_in)
   
        ! Interpolates the NEO h_hat data on to the stella z grid for all three flux surfaces. h_hat data is still on the NEO velocity grids at this stage.

        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_in, neo_grid, 1, neo_h_hat_z_grid, .false.)              ! Calls interpolation on to stella z-grid for the central surface. 
        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_right_in, neo_grid, 1, neo_h_hat_right_z_grid, .false.)  ! Repeat for the right surface.
        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_left_in, neo_grid, 1, neo_h_hat_left_z_grid, .false.)    ! Repeat for the left surface. 

        call get_neo_phi_on_stella_z_grid(neo_phi_in, neo_grid, 1, neo_phi, .false.)                   ! Calls interpolation on to stella z-grid for the central surface.
        call get_neo_phi_on_stella_z_grid(neo_phi_right_in, neo_grid, 1, neo_phi_right, .false.)       ! Repeat for the right surface.
        call get_neo_phi_on_stella_z_grid(neo_phi_left_in, neo_grid, 1, neo_phi_left, .false.)         ! Repeat for the left surface. 

        call get_neo_h_on_stella_grids(neo_h_hat_z_grid, neo_grid, 1, neo_h_local)                ! Now reconstruct H_1 on stellas v∥​ and μ grids for central surface.
        call get_neo_h_on_stella_grids(neo_h_hat_right_z_grid, neo_grid, 1, neo_h_local_right)                                             
        call get_neo_h_on_stella_grids(neo_h_hat_left_z_grid, neo_grid, 1, neo_h_local_left)      

        ! Now compact distribution into 3 indices for use in the GK equation and also for calculating the derivatives.  

        do iz = -nzgrid, nzgrid
            call distribute_vmus_over_procs(neo_h_local(iz, :, :, :, 1), neo_h(iz, :, 1))      
            call distribute_vmus_over_procs(neo_h_local_right(iz, :, :, :, 1), neo_h_right(iz, :, 1))
            call distribute_vmus_over_procs(neo_h_local_left(iz, :, :, :, 1), neo_h_left(iz, :, 1))
        end do        

        ! Now that we have H_1 (not normalised to the Maxwellian here) and ϕ^1_0 for the three flux surfaces, the radial, z, v∥​ and μ derivatives are needed.
        ! Calculate the psi derivative arrays via a simple finite difference method. 

        dneo_h_dpsi = (neo_h_right - neo_h_left) / (2 * drho)
        dneo_phi_dpsi = (neo_phi_right - neo_phi_left) / (2 * drho) 

        ! z derivatives can be obtained via the derivative option of the interpolation routine. 

        call get_dneo_h_dz(neo_h, 1, dneo_h_dz)
        call get_dneo_phi_dz(neo_phi, 1, dneo_phi_dz)

        ! Finally we need the derivatives of the distribution with respect to stellas velocity variables.

        call get_neo_h_velocity_derivs_on_stella_grids(neo_h_hat_z_grid, neo_grid, 1, dneo_h_dvpa, dneo_h_dmu)
   
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

        ! Allocate the ϕ^1_0 arrays (for left and right flux surfaces) on stellas z grid.                                                                                                    
        if (.not. allocated(neo_phi_right)) allocate(neo_phi_right(-nzgrid:nzgrid, neo_grid%n_radial)) 
        if (.not. allocated(neo_phi_left)) allocate(neo_phi_left(-nzgrid:nzgrid, neo_grid%n_radial))    

        ! Allocate the NEO H_1 5D arrays on stellas z and velocity grids.
        if (.not. allocated(neo_h_local)) allocate(neo_h_local(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))   
        if (.not. allocated(neo_h_local_right)) allocate(neo_h_local_right(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))
        if (.not. allocated(neo_h_local_left)) allocate(neo_h_local_left(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))

        ! Allocate the NEO H_1 3D arrays (for left and right flux surfaces) on stellas z and velocity grids.
        if (.not. allocated(neo_h_right)) allocate(neo_h_right(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial)) 
        if (.not. allocated(neo_h_left)) allocate(neo_h_left(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_proc, neo_grid%n_radial))   

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
        if (allocated(neo_phi_right)) deallocate(neo_phi_right) 
        if (allocated(neo_phi_left)) deallocate(neo_phi_left)
        if (allocated(neo_h_local)) deallocate(neo_h_local)       
        if (allocated(neo_h_local_right)) deallocate(neo_h_local_right)
        if (allocated(neo_h_local_left)) deallocate(neo_h_local_left)
        if (allocated(neo_h_right)) deallocate(neo_h_right)   
        if (allocated(neo_h_left)) deallocate(neo_h_left)  

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

      real, dimension(-nzgrid:, :, :), intent(in) :: local
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
        use grids_z, only: nzgrid, zed
        use NEO_interface, only: neo_grid_data
        use periodic_splines, only: periodic_spline, new_periodic_spline, delete_periodic_spline
        use constants, only: twopi

        implicit none

        real, intent(in)  :: neo_h_hat_in(:, :, :, :, :)    
        real, intent(out) :: neo_h_hat_z_grid(-nzgrid:, :, :, :, :)
        type(neo_grid_data), intent(in) :: neo_grid
        integer, intent(in) :: surface_index
        logical, intent(in), optional :: derivative
            
        type(periodic_spline) :: the_spline
        logical :: do_derivative
        integer :: ix, ie, is, iz                  

        if (present(derivative)) then
            do_derivative = derivative
        else
            do_derivative = .false.
        end if

        do ix = 1, neo_grid%n_xi + 1
            do ie = 1, neo_grid%n_energy + 1
                do is = 1, neo_grid%n_species
                    the_spline = new_periodic_spline(neo_grid%theta, neo_h_hat_in(:, ix, ie, is, surface_index), twopi)

                    if (do_derivative) then
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
! ================================================================================================================================================================================= !

    subroutine get_neo_phi_on_stella_z_grid(neo_phi_in, neo_grid, surface_index, neo_phi, derivative)
        use grids_z, only: nzgrid, zed
        use NEO_interface, only: neo_grid_data
        use periodic_splines, only: periodic_spline, new_periodic_spline, delete_periodic_spline
        use constants, only: twopi
        use optionals, only: get_option_with_default
        ! use neoclassical_diagnostics, only: write_neo_phi_on_stella_z_grid_diagnostic
    
        implicit none

        real, intent(in)  :: neo_phi_in(:, :)    
        type(neo_grid_data), intent(in) :: neo_grid
        integer, intent(in) :: surface_index
        real, intent(out) :: neo_phi(-nzgrid:, :)                 
        logical, intent(in), optional :: derivative
        
        logical :: local_derivative
        type(periodic_spline) :: the_spline
        integer :: iz

        the_spline = new_periodic_spline(neo_grid%theta, neo_phi_in(:, surface_index), twopi)

        local_derivative = get_option_with_default(derivative, .false.)

        do iz = -nzgrid, nzgrid
            if (local_derivative) then
                neo_phi(iz, surface_index) = the_spline%derivative(zed(iz))
            else
                neo_phi(iz, surface_index) = the_spline%interpolate(zed(iz))
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
	use grids_z, only: nzgrid, zed
        use grids_velocity, only: vpa, nvpa, mu, nmu
        use geometry, only: bmag                 
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
			xi_in = vpa(iv) / sqrt(vpa(iv)**2 + 2 * bmag(1, iz) * mu(imu)) 
			E_in = (vpa(iv)**2 + 2 * bmag(1, iz) * mu(imu)) / 2

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
        use polynomials, only: get_legendre_array, get_laguerre_array
        use NEO_interface, only: neo_grid_data
    
        implicit none
    
        real, dimension(:, :), intent(in) :: neo_h_hat_z_grid ! A 2-D slice of the 5-dimensional neo_h_hat_z_grid data, at fixed z, species and surface index. 
        real, intent(in) :: xi_in, E_in
        type(neo_grid_data), intent(in) :: neo_grid
    
        integer :: ie, ix

        real :: P(0:neo_grid%n_xi+1)
        real :: L05(0:neo_grid%n_energy+1)
        real :: L15(0:neo_grid%n_energy+1)

        the_sum = 0.0

        call get_legendre_array(xi_in, neo_grid%n_xi+1, P, derivative = .false.)

        ! We need to evaluate the Laguerre function with two different values of k, the first for the first ξ index and the other for all other ξ indices.
        
        call get_laguerre_array(E_in, neo_grid%n_energy+1, 0.5, L05, derivative = .false.)
        call get_laguerre_array(E_in, neo_grid%n_energy+1, 1.5, L15, derivative = .false.)

        ! Now iterate over data and perform the sum. First add in the piece for the first ξ index. 
        do ie = 1, neo_grid%n_energy+1
       	    the_sum = the_sum + P(0) * L05(ie) * neo_h_hat_z_grid(1, ie)

            ! Now add in all the other xi indicies.
            do ix = 2, neo_grid%n_xi+1
                the_sum = the_sum + P(ix) * L15(ie) * neo_h_hat_z_grid(ix, ie)
            end do
        end do
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
        use polynomials, only: get_legendre_array, get_laguerre_array
        use NEO_interface, only: neo_grid_data
    
        implicit none
    
        real, dimension(:, :), intent(in) :: neo_h_hat_z_grid ! A 2-D slice of the 5-dimensional neo_h_hat_z_grid data, at fixed z, species and surface index. 
        real, intent(in) :: xi_in, E_in
        type(neo_grid_data), intent(in) :: neo_grid
    
        integer :: ie, ix

        real :: P(0:neo_grid%n_xi+1)
        real :: dP(0:neo_grid%n_xi+1)
        real :: L05(0:neo_grid%n_energy+1)
        real :: L15(0:neo_grid%n_energy+1)

        the_sum = 0.0

        call get_legendre_array(xi_in, neo_grid%n_xi+1, P, derivative = .true., dP = dP)

        ! We need to evaluate the Laguerre function with two different values of k, the first for the first ξ index and the other for all other ξ indices.
        
        call get_laguerre_array(E_in, neo_grid%n_energy+1, 0.5, L05, derivative = .false.)
        call get_laguerre_array(E_in, neo_grid%n_energy+1, 1.5, L15, derivative = .false.)

        ! Now iterate over data and perform the sum. First add in the piece for the first ξ index. 
        do ie = 1, neo_grid%n_energy+1
       	    the_sum = the_sum + dP(0) * L05(ie) * neo_h_hat_z_grid(1, ie)

            ! Now add in all the other xi indicies.
            do ix = 2, neo_grid%n_xi+1
                the_sum = the_sum + dP(ix) * L15(ie) * neo_h_hat_z_grid(ix, ie)
            end do
        end do
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
        use polynomials, only: get_legendre_array, get_laguerre_array
        use NEO_interface, only: neo_grid_data

        implicit none

        real, dimension(:, :), intent(in) :: neo_h_hat_z_grid ! A 2-D slice of the 5-dimensional neo_h_hat_z_grid data, at fixed z, species and surface index. 
        real, intent(in) :: xi_in, E_in
        type(neo_grid_data), intent(in) :: neo_grid

        integer :: ie, ix

        real :: P(0:neo_grid%n_xi+1)
        real :: L05(0:neo_grid%n_energy+1)
        real :: dL05(0:neo_grid%n_energy+1)
        real :: L15(0:neo_grid%n_energy+1)
        real :: dL15(0:neo_grid%n_energy+1)

        the_sum = 0.0

        call get_legendre_array(xi_in, neo_grid%n_xi+1, P, derivative = .false.)

        ! We need to evaluate the Laguerre function with two different values of k, the first for the first ξ index and the other for all other ξ indices.

        call get_laguerre_array(E_in, neo_grid%n_energy+1, 0.5, L05, derivative = .true., dL=dL05)
        call get_laguerre_array(E_in, neo_grid%n_energy+1, 1.5, L15, derivative = .true., dL=dL15)

        ! Now iterate over data and perform the sum. First add in the piece for the first ξ index. 
        do ie = 1, neo_grid%n_energy+1
            the_sum = the_sum + P(0) * dL05(ie) * neo_h_hat_z_grid(1, ie)

            ! Now add in all the other ξ indicies.
            do ix = 2, neo_grid%n_xi+1
                the_sum = the_sum + P(ix) * dL15(ie) * neo_h_hat_z_grid(ix, ie)
            end do
        end do
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

    subroutine get_neo_h_velocity_derivs_on_stella_grids(neo_h_hat_z_grid, neo_grid, surface_index, dneo_h_dvpa, dneo_h_dmu)
        use grids_z, only: nzgrid, zed
        use grids_velocity, only: vpa, nvpa, mu, nmu
        use geometry, only: bmag
        use NEO_interface, only: neo_grid_data
         
        use, intrinsic :: ieee_arithmetic    ! FOR DEBUGGING. 
 
        implicit none

        real, intent(in)                :: neo_h_hat_z_grid(-nzgrid:, :, :, :, :)
        type(neo_grid_data), intent(in) :: neo_grid
        integer, intent(in)             :: surface_index
        real, intent(out)               :: dneo_h_dvpa(-nzgrid:, :, :)
        real, intent(out)               :: dneo_h_dmu(-nzgrid:, :, :)

        real    :: xi_in, E_in
        integer :: iz, iv, imu, is
        real    :: dneo_h_dxi, dneo_h_dE

        real, allocatable :: dneo_h_dvpa_global(:, :, :, :, :)
        real, allocatable :: dneo_h_dmu_global(:, :, :, :, :)

        ! FOR DEBUGGING: 
    
        real :: dvpa_min, dvpa_max
        real :: dmu_min,  dmu_max
        logical :: bad_data

        bad_data = .false.

        allocate(dneo_h_dvpa_global(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))
        allocate(dneo_h_dmu_global(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))

        do is = 1, neo_grid%n_species
            do iz = -nzgrid, nzgrid
                do iv = 1, nvpa
                    do imu = 1, nmu
                        ! Calculate (ξ_in, E_in) from the (v∥​, μ) stella grid point.
			xi_in = vpa(iv) / sqrt(vpa(iv)**2 + 2 * bmag(1, iz) * mu(imu)) 
			E_in = (vpa(iv)**2 + 2 * bmag(1, iz) * mu(imu)) / 2

                        ! Construct dH_1/dξ|_E at the given (ξ_in, E_in) point. 
                        dneo_h_dxi = get_neo_h_at_xi_deriv_energy(neo_h_hat_z_grid(iz, :, :, is, surface_index), E_in, xi_in, neo_grid)

                        ! Construct dH_1/dE​|_ξ at the given (ξ_in, E_in) point.
                        dneo_h_dE = get_neo_h_at_xi_energy_deriv(neo_h_hat_z_grid(iz, :, :, is, surface_index), E_in, xi_in, neo_grid)

                        ! Construct dH_1/dv∥|_μ at the given (ξ_in, E_in) point. 
                        dneo_h_dvpa_global(iz, iv, imu, is, surface_index) = vpa(iv) * dneo_h_dE + ((1 - xi_in**2)/sqrt(2 * E_in)) * dneo_h_dxi                           
 
                        ! Construct dH_1/dμ|_v∥ at the given (ξ_in, E_in) point.
                        dneo_h_dmu_global(iz, iv, imu, is, surface_index) = bmag(1, iz) * dneo_h_dE - ((bmag(1, iz) * xi_in) / (2 * E_in)) * dneo_h_dxi 
                    end do
                end do
            end do
        end do
 
        ! Collapse derivative arrays into 3 dimensions for use in GKE. 

        do iz = -nzgrid, nzgrid
            call distribute_vmus_over_procs(dneo_h_dvpa_global(iz, :, :, :, 1), dneo_h_dvpa(iz, :, 1))
            call distribute_vmus_over_procs(dneo_h_dmu_global(iz, :, :, :, 1), dneo_h_dmu(iz, :, 1))
        end do 

        ! Deallocate the global arrays. 

        deallocate(dneo_h_dvpa_global)
        deallocate(dneo_h_dmu_global)

        ! FOR DEBUGGING:

        dvpa_min = minval(dneo_h_dvpa)
        dvpa_max = maxval(dneo_h_dvpa)
        dmu_min  = minval(dneo_h_dmu)
        dmu_max  = maxval(dneo_h_dmu)

        if (any(.not. ieee_is_finite(dneo_h_dvpa))) bad_data = .true.
        if (any(.not. ieee_is_finite(dneo_h_dmu)))  bad_data = .true.

        if (bad_data) then
            write(*,*) 'NEO velocity-derivative diagnostic failure'
            write(*,*) 'dvpa min/max = ', dvpa_min, dvpa_max
            write(*,*) 'dmu  min/max = ', dmu_min,  dmu_max
            error stop 'Non-finite values detected in NEO velocity derivatives'
        end if

    end subroutine get_neo_h_velocity_derivs_on_stella_grids


! ================================================================================================================================================================================= !
! -------------------------------- From NEO H_1 on stella z, v∥​ and μ grids, get the derivative of H_1 with respect to z via finite difference methods. --------------------------- !
! ================================================================================================================================================================================= !

    subroutine get_dneo_h_dz(neo_h, surface_index, dneo_h_dz)
        use grids_z, only: nzgrid, zed
        use parallelisation_layouts, only: vmu_lo
 
        implicit none

        real, intent(in)  :: neo_h(-nzgrid:, vmu_lo%llim_proc:, :)
        integer, intent(in) :: surface_index
        real, intent(out) :: dneo_h_dz(-nzgrid:, vmu_lo%llim_proc:, :)         

        real :: dz

        integer :: ivmu, iv, imu, is
        integer :: iz
         
        dz = zed(1) - zed(0)   
        
        ! For the  interior points only. 
        do iz = -nzgrid + 1, nzgrid - 1
            dneo_h_dz(iz, :, surface_index) = (neo_h(iz + 1, :, surface_index) - neo_h(iz - 1, :, surface_index)) / (2.0*dz)
        end do

        ! For the periodic boundary points only.
        dneo_h_dz(-nzgrid, :, surface_index) = (neo_h(-nzgrid + 1, :, surface_index) - neo_h(nzgrid - 1, :, surface_index)) / (2.0 * dz)

        dneo_h_dz(nzgrid, :, surface_index) = dneo_h_dz(-nzgrid, :, surface_index)
    end subroutine get_dneo_h_dz

! ================================================================================================================================================================================= !
! ------------------------------------ From NEO ϕ^1_0 on stella z grid, get the derivative of ϕ^1_0 with respect to z via finite difference methods. ------------------------------ !
! ================================================================================================================================================================================= !

    subroutine get_dneo_phi_dz(neo_phi, surface_index, dneo_phi_dz)
        use grids_z, only: nzgrid, zed
 
        implicit none

        real, intent(in)  :: neo_phi(-nzgrid:, :)
        integer, intent(in) :: surface_index
        real, intent(out) :: dneo_phi_dz(-nzgrid:, :)         

        real :: dz
        integer :: iz
         
        dz = zed(1) - zed(0)                             
        
        ! For the  interior points only. 
        do iz = -nzgrid + 1, nzgrid - 1
            dneo_phi_dz(iz, surface_index) = (neo_phi(iz + 1, surface_index) - neo_phi(iz - 1, surface_index)) / (2.0*dz)
        end do

        ! For the periodic boundary points only.
        dneo_phi_dz(-nzgrid, surface_index) = (neo_phi(-nzgrid + 1, surface_index) - neo_phi(nzgrid - 1, surface_index)) / (2.0 * dz)

        dneo_phi_dz(nzgrid, surface_index) = dneo_phi_dz(-nzgrid, surface_index)
    end subroutine get_dneo_phi_dz 


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------ End Module. -------------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

end module neoclassical_terms_neo
