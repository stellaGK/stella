! ===================================================================================================== !
! ------------------------------------------ NEO INTERFACE -------------------------------------------- !
! ===================================================================================================== !
!
! This module reads in NEO output data associated with the first order neoclassical correction to the 
! equilibrium distributuion, F_1. This is needed as input for second order simulations in stella.   
! 
! NEO output for the distribution correction, H_1, and the electrostatic potential 
! correction, ϕ^1_0, is related to F_1 via
!
! F_1 = H_1 - (e * Z/T) * ϕ^1_0 * F_0
!
! where:
!   - F_1: Total first-order correction to the distribution function.
!   - H_1: NEO's non-adiabatic distribution function correction.
!   - ϕ^1_0: NEO's electrostatic potential correction.
!   - F_0: The leading order Maxwellian distribution.
!   - e, Z, T: Elementary Charge, charge number, and temperature.
!
! stella requires output data from 3 seperate NEO runs (for 3 neighbouring flux surfaces). This is needed
! to calculate the equilibrium gradient drive arising from F_1. Files to be read in are: 
!
! out.neo.f
! out.neo.f.right
! out.neo.f.left
! out.neo.phi
! out.neo.phi.right
! out.neo.phi.left
! out.neo.version
! out.neo.species
! out.neo.grid
! out.neo.version
!
! Files without a "left" or "right" suffix correspond to the central flux surface for which stella
! simulations are to be run.
!
! ===================================================================================================== !

module NEO_interface

    implicit none

    public :: read_basic_neo_files, read_neo_f_and_phi, neo_version_data, neo_equil_data
    public :: neo_grid_data, neo_species_data

    private

! ===================================================================================================== !
! ----------------- Represents the contents of out.neo.version - NEO metadata. ------------------------ ! 
! ===================================================================================================== !

    type neo_version_data
        character(len=:), allocatable :: commit
        character(len=:), allocatable :: system
        character(len=:), allocatable :: date
    end type neo_version_data


! ===================================================================================================== !
! ------------------ Represents the contents of out.neo.grid - NEO grid data. ------------------------- ! 
! ===================================================================================================== !

    type neo_grid_data
        integer :: n_species = -1
        integer :: n_energy = -1
        integer :: n_xi = -1
        integer :: n_theta = -1
        real(8), dimension(:), allocatable :: theta ! Use double precision for floating point numbers. 
        integer :: n_radial = -1
        real(8), dimension(:), allocatable :: radius
    end type neo_grid_data


! ===================================================================================================== !
! ---------------- Represents the contents of out.neo.equil - NEO equilibrium ------------------------- !
! ===================================================================================================== !

    type neo_equil_data
        real(8), dimension(:), allocatable :: radius  
        real(8), dimension(:), allocatable :: radial_electric_field
        real(8), dimension(:), allocatable :: q_safety
        real(8), dimension(:), allocatable :: rho_star
        real(8), dimension(:), allocatable :: major_radius
        real(8), dimension(:), allocatable :: angular_frequency
        real(8), dimension(:), allocatable :: rotation_shear ! 1-dimensional arrays for radially varying quantities. 
     
        real(8), dimension(:, :), allocatable :: density
        real(8), dimension(:, :), allocatable :: temperature
        real(8), dimension(:, :), allocatable :: density_gradient
        real(8), dimension(:, :), allocatable :: temperature_gradient
        real(8), dimension(:, :), allocatable :: collision_frequency ! 2-D arrays for radius and species dependent quantites. 
    end type neo_equil_data


! ===================================================================================================== !
! ----------------- Represents the contents of out.neo.species - NEO species data. -------------------- !
! ===================================================================================================== !

    type neo_species_data
        real(8), dimension(:), allocatable :: mass
        real(8), dimension(:), allocatable :: charge
    end type neo_species_data


! ===================================================================================================== !

contains

! ===================================================================================================== !
! ------------------------------- Utility functions for reading files. -------------------------------- !
! ===================================================================================================== !

! Add these functions to utils? Already exist there? 

! ===================================================================================================== !
! ---------------------------------------- Check file exists. ----------------------------------------- !
! ===================================================================================================== !

    logical function file_exists(filename) result(exists)
        implicit none

        character(len=*), intent(in) :: filename
        inquire(file = trim(filename), exist = exists)
    end function file_exists


! ===================================================================================================== !
! ---------------------------------- Return string from one line. ------------------------------------- !
! ===================================================================================================== !

    function read_line(unit) result(line)
        implicit none

        integer, intent(in) :: unit
        character(len=:), allocatable :: line
        integer, parameter :: chunk_size = 2
        character(len=chunk_size) :: chars
        integer :: err, size_read

        line = ''

        do while (.true.)
            ! Read the next chunk.
            read(unit, '(a)', iostat = err, advance='no', size = size_read) chars

            ! Append read chars to line.
            line = line // chars(:size_read)

            ! Now check if we've reached the end of the record or the end of the file. If so then exit 
            ! the loop.
            if (is_iostat_eor(err) .or. is_iostat_end(err)) exit
        end do
    end function read_line


! ===================================================================================================== !
! -------------------------------------- Reading in NEO Data. ----------------------------------------- !
! ===================================================================================================== !

    subroutine read_basic_neo_files(grid, version, equil, species, basename)
        implicit none

        use optionals, only: get_option_with_default           ! Optionals logic imported from gs2 utils. 
        
        type(neo_grid_data), intent(out) :: grid
        character(len=*), intent(in), optional :: basename
        character(len=:), allocatable :: basename_internal
        type(neo_version_data), intent(out) :: version
        type(neo_equil_data), intent(out), optional :: equil
        type(neo_species_data), intent(out), optional :: species

        basename_internal = get_option_with_default(basename, "out.neo")
        version = get_neo_version_data(basename_internal)
        grid = get_neo_grid_data(basename_internal)

  
    end subroutine read_basic_neo_files


! ===================================================================================================== !
! -- Read out.neo.version - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-version. -- !
! ===================================================================================================== !

    function get_neo_version_data(basename) result(version)
        implicit none

        character(len=*), parameter :: suffix = ".version"
        character(len=*), intent(in) :: basename
        type(neo_version_data) :: version
        character(len=:), allocatable :: filename
        integer :: unit

        filename = basename // suffix

        if ( file_exists(filename) ) then
            open(newunit = unit, file = filename, status="old", action="read")
            version%commit = trim(read_line(unit))
            version%system = trim(read_line(unit))
            version%date = trim(read_line(unit))
            close(unit)
        end if
    end function get_neo_version_data


! ===================================================================================================== !
! ----- Read out.neo.grid - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-grid. ----- !
! ===================================================================================================== !

    function get_neo_grid_data(basename) result(grid)
        implicit none
    
        character(len=*), parameter :: suffix = ".grid"
        character(len=*), intent(in) :: basename
        type(neo_grid_data) :: grid
        character(len=:), allocatable :: filename
        integer :: unit, i

        filename = basename // suffix

        if ( file_exists(filename) ) then
            open(newunit = unit, file = filename, status="old", action="read")

            read(unit, *) grid%n_species
            read(unit, *) grid%n_energy
            read(unit, *) grid%n_xi
            read(unit, *) grid%n_theta

            allocate (grid%theta(grid%n_theta))
                do i = 1, grid%n_theta
                read(unit, *) grid%theta(i)
            end do

            read(unit, *) grid%n_radial

            allocate (grid%radius(grid%n_radial))
            do i = 1, grid%n_radial
                read(unit, *) grid%radius(i)
            end do

            close(unit)
        end if
    end function get_neo_grid_data


! ===================================================================================================== !
! -- Read out.neo.species - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-species. -- !
! ===================================================================================================== !

    function get_neo_species_data(basename, grid) result(species)
        implicit none
    
        character(len=*), parameter :: suffix = ".species"
        character(len=*), intent(in) :: basename
        type(neo_grid_data), intent(in) :: grid ! Needed for grid%n_species
        type(neo_species_data) :: species
        integer :: unit, is
        character(len=:), allocatable :: filename
    
        filename = basename // suffix

        if ( file_exists(filename) ) then
            open(newunit = unit, file = filename, status="old", action="read")
        
            allocate(species%mass(grid%n_species))
            allocate(species%charge(grid%n_species))

            read(unit, *) (species%mass(is), species%charge(is), is = 1, grid%n_species)

            close(unit)
        end if
    end function get_neo_species_data


! ===================================================================================================== !
! ---- Read out.neo.equil - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-equil. ---- !
! ===================================================================================================== !

    function get_neo_equil_data(basename, grid) result(equil)
        implicit none
    
        character(len=*), parameter :: suffix = ".equil"
        character(len=*), intent(in) :: basename
        type(neo_grid_data), intent(in) :: grid
        type(neo_equil_data) :: equil
        integer :: unit, ir, is, err
        character(len=:), allocatable :: filename

        filename = basename // suffix

        if ( file_exists(filename) ) then
            open(newunit = unit, file = filename, status="old", action="read")
            allocate(equil%radius(grid%n_radial))
            allocate(equil%radial_electric_field(grid%n_radial))
            allocate(equil%q_safety(grid%n_radial))
            allocate(equil%rho_star(grid%n_radial))
            allocate(equil%major_radius(grid%n_radial))
            allocate(equil%angular_frequency(grid%n_radial))
            allocate(equil%rotation_shear(grid%n_radial))

            allocate(equil%density(grid%n_species, grid%n_radial))
            allocate(equil%temperature(grid%n_species, grid%n_radial))
            allocate(equil%density_gradient(grid%n_species, grid%n_radial))
            allocate(equil%temperature_gradient(grid%n_species, grid%n_radial))
            allocate(equil%collision_frequency(grid%n_species, grid%n_radial))

            do ir = 1, grid%n_radial
                read(unit, '(e16.8)', advance='no') equil%radius(ir)
                read(unit, '(e16.8)', advance='no') equil%radial_electric_field(ir)
                read(unit, '(e16.8)', advance='no') equil%q_safety(ir)
                read(unit, '(e16.8)', advance='no') equil%rho_star(ir)
                read(unit, '(e16.8)', advance='no') equil%major_radius(ir)
                read(unit, '(e16.8)', advance='no') equil%angular_frequency(ir)
                read(unit, '(e16.8)', advance='no') equil%rotation_shear(ir)

                do is = 1, grid%n_species
                    read(unit, '(e16.8)', advance='no') equil%density(is, ir)
                    read(unit, '(e16.8)', advance='no') equil%temperature(is, ir)
                    read(unit, '(e16.8)', advance='no') equil%density_gradient(is, ir)
                    read(unit, '(e16.8)', advance='no') equil%temperature_gradient(is, ir)
                    read(unit, '(e16.8)', advance='no') equil%collision_frequency(is, ir)
                end do

                ! Advance to the next record. Note we pass iostat here as we don't
                ! care if we hit the end of file.
                read(unit, '()', iostat=err)
            end do

            close(unit)
        end if
    end function get_neo_equil_data


! ===================================================================================================== !
! -------- Read out.neo.f - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-f. -------- !
! -- Note that out.neo.f contains data on the basis functions from which H_1 is to be reconstructed. -- !
! ===================================================================================================== !

    function get_neo_f_data(basename, grid, suffix) result(reshaped_f)
        implicit none
    
        character(len=*), intent(in) :: basename, suffix
        type(neo_grid_data), intent(in) :: grid
        real(8), dimension(:, :, :, :, :), allocatable :: reshaped_f 
        real(8), dimension(:), allocatable :: raw_f                   
        character(len=:), allocatable :: filename
        integer :: unit, i, expected_size
    
        filename = basename // suffix

        expected_size = grid%n_radial * grid%n_species * (grid%n_energy + 1) &
            * (grid%n_xi + 1) * grid%n_theta

        allocate(raw_f(expected_size))

        if (file_exists(filename) ) then
            open(newunit = unit, file = filename, status="old", action="read")

            do i = 1, expected_size
                read(unit, *) raw_f(i) ! Using list-directed read
            end do

            close(unit)
        end if

        reshaped_f = reshape(raw_f, & 
            [grid%n_theta, grid%n_xi+1, grid%n_energy+1, grid%n_species, grid%n_radial])        
    end function get_neo_f_data


! ===================================================================================================== !
! ------ Read out.neo.phi - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-phi. ------ !
! ===================================================================================================== !

    function get_neo_phi_data(basename, grid, suffix) result(phi)
        implicit none
    
        character(len=*), intent(in) :: basename, suffix
        type(neo_grid_data), intent(in) :: grid
        real(8), dimension(:, :), allocatable :: phi
        character(len=:), allocatable :: filename
        integer :: unit, ir
    
        filename = basename // suffix
    
        allocate( phi(grid%n_theta, grid%n_radial))
        
        if (file_exists(filename)) then
            open(newunit = unit, file = filename, status="old", action="read")
            do ir = 1, grid%n_radial
                read(unit, *) phi(:, ir)
            end do
            close(unit)
        end if
    end function get_neo_phi_data


! ===================================================================================================== !
! ----------------------------- Get H_1 and ϕ for a single flux surface. ------------------------------ !
! ===================================================================================================== !

    subroutine read_neo_f_and_phi(neo_f, neo_phi, grid, basename, suffix)
        implicit none

        use iso_fortran_env, only: output_unit
        use optionals, only: get_option_with_default

        real, dimension(:, :, :, :, :), intent(out) :: neo_f
        real, dimension(:, :), intent(out) :: neo_phi
        type(neo_grid_data), intent(in) :: grid
        character(len=*), intent(in), optional :: basename
        character(len=*), intent(in), optional :: suffix
        character(len=:), allocatable :: basename_internal, suffix_internal
        basename_internal = get_option_with_default(basename, "out.neo")
        suffix_internal = get_option_with_default(suffix, "")
    
        neo_f = get_neo_f_data(basename_internal, grid, ".f"//suffix)
        neo_phi = get_neo_phi_data(basename_internal, grid, ".phi"//suffix)
    end subroutine read_neo_f_and_phi

! ===================================================================================================== !
! -------------------------------------- End of NEO interface. ---------------------------------------- !
! ===================================================================================================== !

end module NEO_interface




