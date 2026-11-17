! ================================================================================================================================================================================ !
! --------------------------------------------------------------------------------- NEO INTERFACE -------------------------------------------------------------------------------- !
! ================================================================================================================================================================================ !
!
! This module contains routines for reading in NEO output data, namely the first order non-adiabatic correction to the equilibrium distributuion, H_1, and the first order 
! correction to the equilibrium electrostatic potential, ϕ^1_0. This is needed as input for higher order simulations in stella. The total first order equilibrium distribution
! function, F_1, is given by:
!
! F_1 * F_0 = ( H_1 - Z * ϕ^1_0 ) * F_0
!
! where F_1, H_1 and ϕ^1_0 are dimensionless. 
! 
!   - F_1: Total first-order correction to the equilibrium distribution function.
!   - H_1: NEO's non-adiabatic equilibrium distribution function correction.
!   - ϕ^1_0: NEO's equilibrium electrostatic potential correction.
!   - F_0: The Maxwellian distribution.
!   - e, Z, T: Elementary Charge, species charge number, and species temperature.
!
! stella requires output data from 3 seperate NEO runs for 3 neighbouring flux surfaces. The central surface is the surface of interest for the proceeding gyrokinetic simulation.
! Data from left and right surfaces is  needed to calculate the equilibrium gradient drive arising from F_1. NEO Files to be read in are: 
!
! out.neo.f
! out.neo.f.right
! out.neo.f.left
! out.neo.phi
! out.neo.phi.right
! out.neo.phi.left
! out.neo.version
! out.neo.species
! out.neo.equil
! out.neo.equil.right         < === Profiles on the neighbouring surfaces, for verifying the transport on those surfaces.
! out.neo.equil.left
! out.neo.grid
! out.neo.version
! out.neo.transport          < === Needed for verifying correctness of interpolation & reconstruction.
! out.neo.transport.right
! out.neo.transport.left
!
! Files without a "left" or "right" suffix correspond to the central flux surface for which stella simulations are to be run.
!
! ================================================================================================================================================================================ !

module NEO_interface

    implicit none

    public :: read_basic_neo_files, read_neo_f_and_phi, read_neo_transport, read_neo_equil
    public :: neo_version_data, neo_equil_data, neo_grid_data, neo_species_data, neo_transport_data

    private

! ================================================================================================================================================================================ !
! ------------------------------------------------------------- Represents the contents of out.neo.version - NEO metadata. ------------------------------------------------------- ! 
! ================================================================================================================================================================================ !

    type neo_version_data
        character(len=:), allocatable :: commit
        character(len=:), allocatable :: system
        character(len=:), allocatable :: date
    end type neo_version_data


! ================================================================================================================================================================================ !
! ------------------------------------------------------------- Represents the contents of out.neo.grid - NEO grid data. --------------------------------------------------------- ! 
! ================================================================================================================================================================================ !

    type neo_grid_data
        integer :: n_species = -1
        integer :: n_energy = -1
        integer :: n_xi = -1
        integer :: n_theta = -1
        real, dimension(:), allocatable :: theta 
        integer :: n_radial = -1
        real, dimension(:), allocatable :: radius
    end type neo_grid_data


! ================================================================================================================================================================================ !
! --------------------------------------------------------- Represents the contents of out.neo.equil - NEO equilibrium data. ----------------------------------------------------- !
! ================================================================================================================================================================================ !

    type neo_equil_data
        real, dimension(:), allocatable :: radius  
        real, dimension(:), allocatable :: radial_electric_field
        real, dimension(:), allocatable :: q_safety
        real, dimension(:), allocatable :: rho_star
        real, dimension(:), allocatable :: major_radius
        real, dimension(:), allocatable :: angular_frequency
        real, dimension(:), allocatable :: rotation_shear  
     
        real, dimension(:, :), allocatable :: density
        real, dimension(:, :), allocatable :: temperature
        real, dimension(:, :), allocatable :: density_gradient
        real, dimension(:, :), allocatable :: temperature_gradient
        real, dimension(:, :), allocatable :: collision_frequency 
    end type neo_equil_data


! ================================================================================================================================================================================ !
! --------------------------------------------------------- Represents the contents of out.neo.species - NEO species data. ------------------------------------------------------- !
! ================================================================================================================================================================================ !

    type neo_species_data
        real, dimension(:), allocatable :: mass
        real, dimension(:), allocatable :: charge
    end type neo_species_data

! ================================================================================================================================================================================ !
! --------------------------------------------- Represents the contents of out.neo.transport - NEO transport data (one set per surface). ----------------------------------------- !
! ------------------------------------------------------- See: https://gacode.io/neo/outputs.html#neo-out-neo-transport ---------------------------------------------------------- !
! ================================================================================================================================================================================ !

    type neo_transport_data
        ! .true. only if the file was found and every column was read successfully.
        logical :: found = .false.

        ! \rho = r/a
        real, dimension(:), allocatable :: radius       
        ! (\phi^1_0)^2
        real, dimension(:), allocatable :: phi_sq_avg  
        ! < j_par B > / (e n_norm v_norm B_unit); first order bootstrap current. 
        real, dimension(:), allocatable :: jpar         
        ! v_phi^(0)(theta = 0) / v_norm; zeroth-order toroidal flow at the outboard midplane.  
        real, dimension(:), allocatable :: vtor0_th0    
         ! < u_par^(0) B > / (v_norm B_unit); zeroeth-order parallel flow. 
        real, dimension(:), allocatable :: uparB0        
        ! Gamma_sigma / (n_norm v_norm) ; Second-order radial particle flux.
        real, dimension(:, :), allocatable :: pflux      
        ! Q_sigma / (n_norm v_norm T_norm); Second-order radial energy flux.
        real, dimension(:, :), allocatable :: eflux  
        ! Pi_sigma / (n_norm T_norm a_norm); Second-order radial momentum flux.
        real, dimension(:, :), allocatable :: mflux      
        ! < u_par,sigma B > / (v_norm B_unit); First-order parallel flow.
        real, dimension(:, :), allocatable :: uparB      
        ! k_sigma; First-order dimensionless flow coefficient.        
        real, dimension(:, :), allocatable :: k_small    
        ! K_sigma / (n_norm v_norm / B_unit); First-order dimensional flow coefficient.
        real, dimension(:, :), allocatable :: K_big      
        ! v_theta,sigma(theta = 0) / v_norm; First-order poloidal flow, outboard midplane.
        real, dimension(:, :), allocatable :: vpol_th0   
        ! v_phi,sigma(theta = 0) / v_norm; First-order toroidal flow, outboard midplane.
        real, dimension(:, :), allocatable :: vtor_th0
    end type neo_transport_data


! ================================================================================================================================================================================ !

contains

! ================================================================================================================================================================================ !
! -------------------------------------------------------------------------- Utility functions for reading files. ---------------------------------------------------------------- !
! ================================================================================================================================================================================ !

! ================================================================================================================================================================================ !
! --------------------------------------------------------------------------------- Check file exists. --------------------------------------------------------------------------- !
! ================================================================================================================================================================================ !

    logical function file_exists(filename) result(exists)
        implicit none

        character(len=*), intent(in) :: filename
        inquire(file = trim(filename), exist = exists)
    end function file_exists


! ================================================================================================================================================================================ !
! ----------------------------------------------------------------------------- Return string from one line. --------------------------------------------------------------------- !
! ================================================================================================================================================================================ !

    function read_line(unit) result(line)
        implicit none

        integer, intent(in) :: unit
        character(len=:), allocatable :: line
        integer, parameter :: chunk_size = 2
        character(len=chunk_size) :: chars
        integer :: err, size_read

        line = ''

        do while (.true.)
            read(unit, '(a)', iostat = err, advance='no', size = size_read) chars

            line = line // chars(:size_read)

            if (is_iostat_eor(err) .or. is_iostat_end(err)) exit
        end do
    end function read_line


! ================================================================================================================================================================================ !
! ----------------------------------------------------------------- Number of whitespace-separated fields in a line. ------------------------------------------------------------- !
! ================================================================================================================================================================================ !

    pure integer function count_fields(line) result(n)
        implicit none

        character(len=*), intent(in) :: line
        integer :: i
        logical :: in_field

        n = 0
        in_field = .false.
        do i = 1, len(line)
            if (line(i:i) == ' ' .or. line(i:i) == achar(9)) then
                in_field = .false.
            else if (.not. in_field) then
                in_field = .true.
                n = n + 1
            end if
        end do
    end function count_fields


! ================================================================================================================================================================================ !
! --------------------------------------------------------------------------------- Reading in NEO Data. ------------------------------------------------------------------------- !
! ================================================================================================================================================================================ !

    subroutine read_basic_neo_files(grid, version, equil, species, basename)
        use optionals, only: get_option_with_default       
     
        implicit none
   
        type(neo_grid_data), intent(out) :: grid
        character(len=*), intent(in), optional :: basename
        character(len=:), allocatable :: basename_internal
        type(neo_version_data), intent(out) :: version
        type(neo_equil_data), intent(out), optional :: equil
        type(neo_species_data), intent(out), optional :: species

        basename_internal = get_option_with_default(basename, "out.neo")
        version = get_neo_version_data(basename_internal)
        grid = get_neo_grid_data(basename_internal)
        if (present(species)) species = get_neo_species_data(basename_internal, grid)
        if (present(equil))   equil   = get_neo_equil_data(basename_internal, grid)


  
    end subroutine read_basic_neo_files


! ================================================================================================================================================================================ !
! ------------------------------------------ Read out.neo.version - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-version. ------------------------------------- !
! ================================================================================================================================================================================ !

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


! ================================================================================================================================================================================ !
! ----------------------------------------------- Read out.neo.grid - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-grid. -------------------------------------- !
! ================================================================================================================================================================================ !

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


! ================================================================================================================================================================================ !
! ----------------------------------------- Read out.neo.species - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-species. -------------------------------------- !
! ================================================================================================================================================================================ !

    function get_neo_species_data(basename, grid) result(species)
        implicit none
        character(len=*), parameter :: suffix = ".species"
        character(len=*), intent(in) :: basename
        type(neo_grid_data), intent(in) :: grid
        type(neo_species_data) :: species
        integer :: unit, is
        character(len=:), allocatable :: filename
        filename = basename // suffix

        if ( file_exists(filename) ) then
            open(newunit = unit, file = filename, status="old", action="read")
            allocate(species%mass(grid%n_species))
            allocate(species%charge(grid%n_species))

            do is = 1, grid%n_species
                read(unit, '(e16.8)', advance='no') species%mass(is)
                read(unit, '(e16.8)', advance='no') species%charge(is)
            end do

            close(unit)
        end if
    end function get_neo_species_data


! ================================================================================================================================================================================ !
! ---------------------------------------- Read out.neo.equil - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-equil. ------------------------------------------- !
! ================================================================================================================================================================================ !

    function get_neo_equil_data(basename, grid, extra_suffix) result(equil)
        use iso_fortran_env, only: output_unit
        use optionals, only: get_option_with_default

        implicit none
    
        character(len=*), parameter :: suffix = ".equil"
        character(len=*), intent(in) :: basename
        type(neo_grid_data), intent(in) :: grid
        ! Optional surface suffix, e.g. ".left" -> out.neo.equil.left.
        character(len=*), intent(in), optional :: extra_suffix
        type(neo_equil_data) :: equil
        integer :: unit, ir, is, err
        character(len=:), allocatable :: filename

        filename = basename // suffix // get_option_with_default(extra_suffix, "")

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

                read(unit, '()', iostat=err)
            end do

            close(unit)
        else
            write (output_unit, '(3a)') ' WARNING: ', filename, ' not found.'
        end if
    end function get_neo_equil_data


! ================================================================================================================================================================================ !
! ------------------------------------------------------ Read the NEO equilibrium profiles of one flux surface (suffix '', '.left' or '.right'). --------------------------------- !
! ================================================================================================================================================================================ !

    subroutine read_neo_equil(equil, grid, basename, suffix)
        use optionals, only: get_option_with_default

        implicit none

        type(neo_equil_data), intent(out)      :: equil
        type(neo_grid_data), intent(in)        :: grid
        character(len=*), intent(in), optional :: basename
        character(len=*), intent(in), optional :: suffix
        character(len=:), allocatable :: basename_internal, suffix_internal

        basename_internal = get_option_with_default(basename, "out.neo")
        suffix_internal = get_option_with_default(suffix, "")

        equil = get_neo_equil_data(basename_internal, grid, suffix_internal)
    end subroutine read_neo_equil


! ================================================================================================================================================================================ !
! ------------------------------------------------- Read out.neo.f - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-f. ------------------------------------------ !
! ------------------------------------------ Note that out.neo.f contains data on the basis functions from which H_1 is to be reconstructed. ------------------------------------- !
! ================================================================================================================================================================================ !

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
                read(unit, *) raw_f(i) 
            end do

            close(unit)
        end if

        reshaped_f = reshape(raw_f, & 
            [grid%n_theta, grid%n_xi+1, grid%n_energy+1, grid%n_species, grid%n_radial])        
    end function get_neo_f_data


! ================================================================================================================================================================================ !
! ------------------------------------------------- Read out.neo.phi - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-phi. -------------------------------------- !
! ================================================================================================================================================================================ !

    function get_neo_phi_data(basename, grid, suffix) result(phi)
        implicit none
    
        character(len=*), intent(in) :: basename, suffix
        type(neo_grid_data), intent(in) :: grid
        real, dimension(:, :), allocatable :: phi
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


! ================================================================================================================================================================================ !
! ------------------------------------------------------------------------- Get H_1 and ϕ^1_0 for a single flux surface. ------------------------------------------------------------- !
! ================================================================================================================================================================================ !

    subroutine read_neo_f_and_phi(neo_f, neo_phi, grid, basename, suffix)
        use iso_fortran_env, only: output_unit
        use optionals, only: get_option_with_default

        implicit none

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


! ================================================================================================================================================================================ !
! ---------------------------------------------------- Read the NEO transport coefficients of one flux surface (suffix '', '.left' or '.right'). --------------------------------- !
! ================================================================================================================================================================================ !

    subroutine read_neo_transport(transport, grid, basename, suffix, verbose)
        ! Optionals.
        use optionals, only: get_option_with_default

        implicit none

        type(neo_transport_data), intent(out) :: transport
        type(neo_grid_data), intent(in) :: grid
        character(len=*), intent(in), optional :: basename
        character(len=*), intent(in), optional :: suffix
        
        ! Local variables. 
        ! Print the column check and the values read (pass proc0).
        logical, intent(in), optional :: verbose      
        character(len=:), allocatable :: basename_internal, suffix_internal

        basename_internal = get_option_with_default(basename, "out.neo")
        suffix_internal = get_option_with_default(suffix, "")

        transport = get_neo_transport_data(basename_internal, grid, suffix_internal, &
        get_option_with_default(verbose, .false.))
    end subroutine read_neo_transport


! ================================================================================================================================================================================ !
! ----------------------------------------- Read out.neo.transport - see: http://gafusion.github.io/doc/neo/outputs.html#neo-out-neo-transport. ---------------------------------- !
! ================================================================================================================================================================================ !

    function get_neo_transport_data(basename, grid, suffix, verbose) result(transport)
        ! Fortran. 
        use iso_fortran_env, only: output_unit

        implicit none

        character(len=*), intent(in) :: basename, suffix
        type(neo_grid_data), intent(in) :: grid
        logical, intent(in) :: verbose
        type(neo_transport_data) :: transport

        ! Local variables. 
        character(len=:), allocatable :: filename, first_line
        real, dimension(:), allocatable :: row
        integer :: unit, ir, is, ncol_expected, ncol_found, ios, c0

        filename = basename // ".transport" // suffix
        ncol_expected = 5 + 8 * grid%n_species
        transport%found = .false.

        if (.not. file_exists(filename)) then
            if (verbose) write (output_unit, '(3a)') ' WARNING: ', filename, ' not found: NEO transport comparison skipped.'
            return
        end if

        ! Check the number of columns on the first row before reading anything.
        open(newunit = unit, file = filename, status = "old", action = "read")
        first_line = read_line(unit)
        close(unit)
        ncol_found = count_fields(first_line)
        if (ncol_found /= ncol_expected) then
            if (verbose) write (output_unit, '(3a, i0, a, i0, a, i0, a)') &
            ' WARNING: ', filename, ' has ', ncol_found, &
            ' columns, expected 5 + 8 x ', grid%n_species, &
            ' = ', ncol_expected, ': NEO transport comparison skipped.'
            return
        end if

        ! Allocate transport coeffecient arrays. 
        allocate(transport%radius(grid%n_radial)) 
        allocate(transport%phi_sq_avg(grid%n_radial))
        allocate(transport%jpar(grid%n_radial))
        allocate(transport%vtor0_th0(grid%n_radial))
        allocate(transport%uparB0(grid%n_radial))
        allocate(transport%pflux(grid%n_species, grid%n_radial))
        allocate(transport%eflux(grid%n_species, grid%n_radial))
        allocate(transport%mflux(grid%n_species, grid%n_radial))
        allocate(transport%uparB(grid%n_species, grid%n_radial))
        allocate(transport%k_small(grid%n_species, grid%n_radial)) 
        allocate(transport%K_big(grid%n_species, grid%n_radial))
        allocate(transport%vpol_th0(grid%n_species, grid%n_radial))
        allocate(transport%vtor_th0(grid%n_species, grid%n_radial))
        allocate(row(ncol_expected))

        open(newunit = unit, file = filename, status = "old", action = "read")
        do ir = 1, grid%n_radial
            read(unit, *, iostat = ios) row
            if (ios /= 0) then
                if (verbose) write (output_unit, '(3a, i0, a)') &
                ' WARNING: could not read row ', filename, ' ', ir, &
                ': NEO transport comparison skipped.'
                close(unit)
                return
            end if

            transport%radius(ir)     = row(1)
            transport%phi_sq_avg(ir) = row(2)
            transport%jpar(ir)       = row(3)
            transport%vtor0_th0(ir)  = row(4)
            transport%uparB0(ir)     = row(5)
            do is = 1, grid%n_species
                c0 = 5 + 8 * (is - 1)
                transport%pflux(is, ir)    = row(c0 + 1)
                transport%eflux(is, ir)    = row(c0 + 2)
                transport%mflux(is, ir)    = row(c0 + 3)
                transport%uparB(is, ir)    = row(c0 + 4)
                transport%k_small(is, ir)  = row(c0 + 5)
                transport%K_big(is, ir)    = row(c0 + 6)
                transport%vpol_th0(is, ir) = row(c0 + 7)
                transport%vtor_th0(is, ir) = row(c0 + 8)
            end do
        end do
        close(unit)

        transport%found = .true.

        if (verbose) then
            write (output_unit, '(3a, i0, a)') ' Read ', filename, ' (', ncol_found, ' columns):'
            do ir = 1, grid%n_radial
                write (output_unit, '(a, es14.6)') '############################################################'
                write (output_unit, '(a, es14.6)') '   r/a                           = ', transport%radius(ir)
                write (output_unit, '(a, es14.6)') '   <(e Phi_1/T_norm)^2>          = ', transport%phi_sq_avg(ir)
                write (output_unit, '(a, es14.6)') '   <j_par B>/(e n v_norm B_unit) = ', transport%jpar(ir)
                write (output_unit, '(a, 2es14.6)') '   zeroth-order v_tor, <u_par B> = ', transport%vtor0_th0(ir), transport%uparB0(ir)
                write (output_unit, '(a)') '   species   Gamma          Q              Pi             <u_par B>      k              K' // &
                                           '              v_pol(0)       v_tor(0)'
                do is = 1, grid%n_species
                    write (output_unit, '(i9, 8es15.6)') is, transport%pflux(is, ir), transport%eflux(is, ir), &
                    transport%mflux(is, ir), transport%uparB(is, ir), transport%k_small(is, ir), &
                    transport%K_big(is, ir), transport%vpol_th0(is, ir), transport%vtor_th0(is, ir)
                end do
                write (output_unit, '(a, es14.6)') '############################################################'
            end do
        end if
    end function get_neo_transport_data


! ================================================================================================================================================================================ !
! -------------------------------------------------------------------------------- End of NEO interface. ------------------------------------------------------------------------- ! 
! ================================================================================================================================================================================ !

end module NEO_interface