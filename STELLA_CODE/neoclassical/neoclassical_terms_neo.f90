! ================================================================================================================================================================================= !
! ------------------------ Routines for remapping NEO data on to stella grids and calculating all quantites needed for higher order GK calculations. ------------------------------ !
! ================================================================================================================================================================================= !
! 
! NEO uses pitch angle cosine, ξ = v∥​/v, and normalised energy, E, for velocity coordinates. stella uses v∥​ and μ. A remapping of the NEO data on to the stella grids is required. 
! 
! We must also reconstruct the NEO H_1 from the amplitudes provided by out.neo.f: see https://gacode.io/neo/outputs.html#neo-out-neo-f for details on the Legendre/Laguerre
! representation of H_1 in terms of the amplitudes. 
!
! This representation can also be used to calculate the derivatives of H_1 in v∥​ and μ. The spatial gradient of H_1 (and F_1) may be calculated by a finite differences method. 
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

    public :: h_neo, dh_neo_dpsi, dh_neo_dtheta, dh_neo_denergy, dh_neo_dxi        ! Will represent NEO's distribution and its derivatives in real space and velocity space. 
    public :: phi_neo, dphi_neo_dpsi, dphi_neo_dtheta                              ! Will represents NEO's ϕ^1_0 and its derivatives in real space. 

    public :: initialised_neoclassical_terms_neo

    private

    logical :: include_neoclassical_terms                                           
    integer :: neo_option_switch                                                   ! Should be = 2 for NEO.
    integer :: nradii                                                              ! Can only be 3.
    real    :: drho                                                                ! Typically taken as 0.01.

    integer :: iz, unit
    character(len=128) :: filename

    real(8), dimension(:, :, :), allocatable :: h_neo, dh_neo_dpsi, dh_neo_dtheta 
    real(8), dimension(:, :), allocatable :: phi_neo, dphi_neo_dpsi, dphi_neo_dtheta
    real(8), dimension(:, :, :), allocatable :: dh_neo_denergy 
    real(8), dimension(:, :, :), allocatable :: dh_neo_dxi
    
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
        use grids_z, only: nzgrid
        use grids_velocity, only: nvpa, nmu 
        use grids_species, only: nspec
        use NEO_interface, only: read_basic_neo_files, read_neo_f_and_phi, neo_grid_data, neo_version_data        

        implicit none

        real, dimension(:, :, :, :, :), allocatable :: neo_h_hat_in, neo_h_hat_right_in, neo_h_hat_left_in ! Holds vectors for reconstructing NEO H_1 on 3 flux surfaces.
        real, dimension(:, :), allocatable :: neo_phi_in, neo_phi_right_in, neo_phi_left_in                ! Holds NEO ϕ^1_0 on 3 flux surfaces.

        ! Intermediate arrays hold NEO h_hat data evaluated on the stella z grid, but on the NEO velocity grids.
        real, dimension(:, :, :, :, :), allocatable :: neo_h_hat_z_grid, neo_h_hat_right_z_grid, neo_h_hat_left_z_grid 
        real, dimension(:, :), allocatable :: neo_phi_z_grid, neo_phi_right_z_grid, neo_phi_left_z_grid                ! Holds NEO ϕ^1_0 data evaluated on the stella z grid.
        
        ! Holds NEO H_1 data evaluated on the stella z, v∥​ and μ grids. Since ϕ^1_0 is independent of velocity variables, there are no accompanying arrays for ϕ^1_0 here.
        real, dimension(:, :, :, :, :), allocatable :: neo_h, neo_h_right, neo_h_left

        ! ====================================================================================== !
        ! ------------------------ STILL NEED ARRAYS FOR DERIVATIVES HERE! --------------------- !
        ! ====================================================================================== !

        integer :: surface_index      

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
       
        allocate(neo_h_hat_in(neo_grid%n_theta, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial))       ! Allocate for the h_hat arrays on NEO grids.  
        allocate(neo_h_hat_right_in(neo_grid%n_theta, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial))
        allocate(neo_h_hat_left_in(neo_grid%n_theta, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial))

        allocate(neo_phi_in(neo_grid%n_theta, neo_grid%n_radial)) ! Allocate for the phi arrays on NEO grid. 
        allocate(neo_phi_right_in(neo_grid%n_theta, neo_grid%n_radial))
        allocate(neo_phi_left_in(neo_grid%n_theta, neo_grid%n_radial))

        if (proc0) then
            call read_neo_f_and_phi(neo_h_hat_in, neo_phi_in, neo_grid)                                ! Read in NEO h and ϕ^1_0 data for the central surface. 
            call read_neo_f_and_phi(neo_h_hat_right_in, neo_phi_right_in, neo_grid, suffix = '.right') ! Repeat for the right surface.  
            call read_neo_f_and_phi(neo_h_hat_left_in, neo_phi_left_in, neo_grid, suffix = '.left')    ! Repeat for the left surface.
        end if        
        
        call broadcast(neo_h_hat_in) ; call broadcast(neo_phi_in)             ! Broadcast the central surface data. 
        call broadcast(neo_h_hat_right_in) ; call broadcast(neo_phi_right_in) ! Repeat for the right surface.
        call broadcast(neo_h_hat_left_in) ; call broadcast(neo_phi_left_in)   ! Repeat for the left surface.

        ! Allocates size of the NEO h_hat data on stellas z grid.

        allocate(neo_h_hat_z_grid(-nzgrid:nzgrid, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial))  
        allocate(neo_h_hat_right_z_grid(-nzgrid:nzgrid, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial))  ! Repeat for the right surface.
        allocate(neo_h_hat_left_z_grid(-nzgrid:nzgrid, neo_grid%n_xi+1, neo_grid%n_energy+1, neo_grid%n_species, neo_grid%n_radial)) ! Repeat for the left surface.
   
        ! Interpolates the NEO h_hat data on to the stella z grid for all three flux surfaces. h_hat data is still on the NEO velocity grids at this stage.

        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_in, neo_grid, 1, neo_h_hat_z_grid)                      ! Calls interpolation on to stella z-grid for the central surface. 
        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_right_in, neo_grid, 1, neo_h_hat_right_z_grid, 'right') ! Repeat for the right surface.
        call get_neo_h_hat_on_stella_z_grid(neo_h_hat_left_in, neo_grid, 1, neo_h_hat_left_z_grid, 'left')    ! Repeat for the left surface.

        ! Now repeat this process for the ϕ^1_0 data sets. 

        allocate(neo_phi_z_grid(-nzgrid:nzgrid, neo_grid%n_radial))                                           ! Allocates size of the NEO ϕ^1_0 data on stellas z grid.
        allocate(neo_phi_right_z_grid(-nzgrid:nzgrid, neo_grid%n_radial))                                     ! Repeat for the right surface. 
        allocate(neo_phi_left_z_grid(-nzgrid:nzgrid, neo_grid%n_radial))                                      ! Repeat for the left surface.

        call get_neo_phi_on_stella_z_grid(neo_phi_in, neo_grid, 1, neo_phi_z_grid)                            ! Calls interpolation on to stella z-grid for the central surface.
        call get_neo_phi_on_stella_z_grid(neo_phi_right_in, neo_grid, 1, neo_phi_right_z_grid, 'right')       ! Repeat for the right surface.
        call get_neo_phi_on_stella_z_grid(neo_phi_left_in, neo_grid, 1, neo_phi_left_z_grid, 'left')          ! Repeat for the left surface. 

        ! Now, we need to construct H_1 from the interpolated neo_h_hat data. Allocate the sizes of the datasets on the stella grids. 

        allocate(neo_h(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))
        allocate(neo_h_right(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))
        allocate(neo_h_left(-nzgrid:nzgrid, nvpa, nmu, neo_grid%n_species, neo_grid%n_radial))   

        call get_neo_h_on_stella_grids(neo_h_hat_z_grid, neo_grid, 1, neo_h)                         ! Now reconstruct H_1 on stellas v∥​ and μ grids for central surface.
        call get_neo_h_on_stella_grids(neo_h_hat_right_z_grid, neo_grid, 1, neo_h_right, 'right')    ! Repeat for the right surface.                                         
        call get_neo_h_on_stella_grids(neo_h_hat_left_z_grid, neo_grid, 1, neo_h_left, 'left')       ! Repeat for the left surface.
        
        ! Now that we have H_1 (not normalised to the Maxwellian here) and ϕ^1_0 for the three flux surfaces, the radial, z, v∥​ and μ derivatives are needed. 

    end subroutine init_neoclassical_terms_neo


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------- Finish the neoclassical terms. ---------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine finish_neoclassical_terms_neo
        implicit none

        call deallocate_arrays
    end subroutine finish_neoclassical_terms_neo


! ================================================================================================================================================================================= !
! ------------------------------------------------------------------------------- Utilities. -------------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

! ================================================================================================================================================================================= !
! --------------------------------- Allocates module level arrays for H_1, ϕ^1_0 and their derivatives in real space and velocity space. ------------------------------------------ !
! ================================================================================================================================================================================= !

    subroutine allocate_arrays
        implicit none
    end subroutine allocate_arrays


! ================================================================================================================================================================================= !
! --------------------------------------------------------------------- Deallocates module level arrays. -------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine deallocate_arrays
        implicit none
    end subroutine deallocate_arrays


! ================================================================================================================================================================================= !
! ---------------------------------- Interpolate the NEO h_hat data from the NEO θ grid to the stella z grid for a specified flux surface. ---------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine get_neo_h_hat_on_stella_z_grid(neo_h_hat_in, neo_grid, surface_index, neo_h_hat_z_grid, suffix)
        use grids_z, only: nzgrid, zed
        use NEO_interface, only: neo_grid_data
        use splines, only: linear_interp_periodic
        use constants, only: twopi
    
        implicit none

        real, intent(in)  :: neo_h_hat_in(:, :, :, :, :)    
        real, intent(out) :: neo_h_hat_z_grid(-nzgrid:, :, :, :, :)

        type(neo_grid_data), intent(in) :: neo_grid
        integer, intent(in) :: surface_index
        integer :: ix, ie, is, iz                                      ! ix is for NEO ξ = v∥​/v, ie for NEO E and is for species.
        character(len=*), intent(in), optional :: suffix               ! For saving interpolated data on surfaces. 
        integer :: unit
        character(len=256) :: filename

        do ix = 1, neo_grid%n_xi + 1
            do ie = 1, neo_grid%n_energy + 1
                do is = 1, neo_grid%n_species
                    call linear_interp_periodic(neo_grid%theta, neo_h_hat_in(:, ix, ie, is, surface_index), zed, neo_h_hat_z_grid(:, ix, ie, is, surface_index), twopi) 
		end do
	    end do
	end do
       
        ! ================================================= !
        ! ----------------- Diagnostic. ------------------- !
        ! ================================================= !

        ! unit = 99 

        ! if (present(suffix)) then
            ! write(filename,'("neo_h_hat_z_grid_",A,"_surf_",I0,".dat")') trim(suffix), surface_index
        ! else
            ! write(filename,'("neo_h_hat_z_grid_surf_",I0,".dat")') surface_index
        ! end if

        ! open(unit=unit, file=filename, status='replace', action='write')

        ! write(unit,'(A)') '# Diagnostic output: neo_h_hat_z_grid'
        ! write(unit,'(A,I0)') '# surface_index = ', surface_index
        ! write(unit,'(A)') '# Columns:'
        ! write(unit,'(A)') '#   iz   : z-grid index'
        ! write(unit,'(A)') '#   ix   : NEO xi index'
        ! write(unit,'(A)') '#   ie   : NEO energy index'
        ! write(unit,'(A)') '#   is   : species index'
        ! write(unit,'(A)') '#   z    : stella z coordinate'
        ! write(unit,'(A)') '#   hhat : interpolated neo_h_hat on stella z-grid'
        ! write(unit,'(A)') '#'
        ! write(unit,'(A)') '# iz   il   ie   is        z              hhat'
        ! -----------------------------------

        ! do ix = 1, neo_grid%n_xi + 1
            ! do ie = 1, neo_grid%n_energy + 1
                ! do is = 1, neo_grid%n_species
                    ! do iz = -nzgrid, nzgrid
                        ! write(unit,'(I6,1X,I4,1X,I4,1X,I4,1X,ES16.8,1X,ES16.8)') iz, ix, ie, is, zed(iz), neo_h_hat_z_grid(iz, ix, ie, is, surface_index)
                    ! end do
                ! end do
            ! end do
        ! end do

        ! close(unit)
    end subroutine get_neo_h_hat_on_stella_z_grid


! ================================================================================================================================================================================= !
! --------------------------------- Interpolate the NEO ϕ^1_0 data from the NEO θ grid to the stella z grid for a specified flux surface. ----------------------------------------- !
! ================================================================================================================================================================================= !

    subroutine get_neo_phi_on_stella_z_grid(neo_phi_in, neo_grid, surface_index, neo_phi_z_grid, suffix)
        use grids_z, only: nzgrid, zed
        use NEO_interface, only: neo_grid_data
        use splines, only: linear_interp_periodic
        use constants, only: twopi
    
        implicit none

        real, intent(in)  :: neo_phi_in(:, :)    
        type(neo_grid_data), intent(in) :: neo_grid
        character(len=*), intent(in), optional :: suffix               ! For saving interpolated data on surfaces.
        integer, intent(in) :: surface_index
        real, intent(out) :: neo_phi_z_grid(-nzgrid:, :)        
         
        integer :: unit
        character(len=256) :: filename

        call linear_interp_periodic(neo_grid%theta, neo_phi_in(:, surface_index), zed, neo_phi_z_grid(:,  surface_index), twopi)

        ! ================================================= !
        ! ----------------- Diagnostic. ------------------- !
        ! ================================================= !

        ! unit = 99 

        ! if (present(suffix)) then
            ! write(filename,'("neo_phi_z_grid_",A,"_surf_",I0,".dat")') trim(suffix), surface_index
        ! else
            ! write(filename,'("neo_phi_z_grid",I0,".dat")') surface_index
        ! end if

        ! open(unit=unit, file=filename, status='replace', action='write')

        ! write(unit,'(A)') '# Diagnostic output: neo_phi_z_grid'
        ! write(unit,'(A,I0)') '# surface_index = ', surface_index
        ! write(unit,'(A)') '# Columns:'
        ! write(unit,'(A)') '#   iz   : z-grid index'
        ! write(unit,'(A)') '#   z    : stella z coordinate'
        ! write(unit,'(A)') '#   phi^1_0 : interpolated neo_phi on stella z-grid'
        ! write(unit,'(A)') '#'
        ! write(unit,'(A)') '# iz            z              hhat'
        ! -----------------------------------

        ! do iz = -nzgrid, nzgrid
            ! write(unit,'(I6,ES16.8,1X,ES16.8)') iz, zed(iz), neo_phi_z_grid(iz, surface_index)
        ! end do
        ! close(unit)
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
        character(len=*), intent(in), optional :: suffix               ! For saving interpolated data on surfaces.
        
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

        ! ================================================= !
        ! ----------------- Diagnostic. ------------------- !
        ! ================================================= !

        ! unit = 99

        ! if (present(suffix)) then
            ! write(filename,'("neo_h_stella_grids_",A,"_surf_",I0,".dat")') trim(suffix), surface_index
        ! else
            ! write(filename,'("neo_h_stella_grids_",I0,".dat")') surface_index
        ! end if

        ! open(unit=unit, file=filename, status='replace', action='write')

        ! write(unit,'(A)') '# Diagnostic output: neo_h on stella grids'
        ! write(unit,'(A,I0)') '# surface_index = ', surface_index
        ! write(unit,'(A)') '# Columns:'
        ! write(unit,'(A)') '#   iz   : z-grid index'
        ! write(unit,'(A)') '#   iv   : v_parallel grid index'
        ! write(unit,'(A)') '#   imu  : mu grid index'
        ! write(unit,'(A)') '#   is   : species index'
        ! write(unit,'(A)') '#   z    : stella z coordinate'
        ! write(unit,'(A)') '#   vpa  : parallel velocity'
        ! write(unit,'(A)') '#   mu   : magnetic moment'
        ! write(unit,'(A)') '#   h    : reconstructed neo_h'
        ! write(unit,'(A)') '#'
        ! write(unit,'(A)') '# iz    iv   imu   is        z              vpa              mu               h'

        ! do is = 1, neo_grid%n_species
            ! do iz = -nzgrid, nzgrid
                ! do iv = 1, nvpa
                    ! do imu = 1, nmu
                        ! write(unit,'(I6,1X,I6,1X,I6,1X,I4,1X,ES16.8,1X,ES16.8,1X,ES16.8,1X,ES16.8)') iz, iv, imu, is, zed(iz), vpa(iv), mu(imu), neo_h(iz, iv, imu, is, surface_index)
                    ! end do
                ! end do
            ! end do
        ! end do
        ! close(unit)
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

        call get_legendre_array(xi_in, neo_grid%n_xi+1, P)

        ! We need to evaluate the Laguerre function with two different values of k, the first for the first ξ index and the other for all other ξ indices.
        
        call get_laguerre_array(E_in, neo_grid%n_energy+1, 0.5, L05)
        call get_laguerre_array(E_in, neo_grid%n_energy+1, 1.5, L15)

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
! ------------------------------------------------------------------------------ End Module. -------------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

end module neoclassical_terms_neo
