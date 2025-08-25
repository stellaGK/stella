module namelist_species
    
    use stella_common_types, only: spec_type

    implicit none

    public :: read_namelist_species_options
    public :: read_namelist_species_stella
    public :: species_option_stella, species_option_inputprofs, species_option_euterpe, species_option_multibox
    public :: ion_species, electron_species, slowing_down_species, tracer_species
    public :: read_namelist_euterpe_parameters

    private

    integer, parameter :: species_option_stella = 1
    integer, parameter :: species_option_inputprofs = 2
    integer, parameter :: species_option_euterpe = 3
    integer, parameter :: species_option_multibox = 4

    integer, parameter :: ion_species = 1
    integer, parameter :: electron_species = 2 ! for collision operator
    integer, parameter :: slowing_down_species = 3 ! slowing-down distn
    integer, parameter :: tracer_species = 4 ! for test particle diffusion studies
    
    ! Internal variables
    integer :: in_file
    logical :: dexist


contains
    !****************************************************************************
    !                              SPECIES OPTIONS                              !
    !****************************************************************************
    subroutine read_namelist_species_options(nspec, species_option_switch, &
                                read_profile_variation, write_profile_variation, ecoll_zeff)

        use mp, only: proc0

        implicit none

        integer, intent (out) :: nspec
        logical, intent (out) :: read_profile_variation, write_profile_variation
        integer, intent (out) :: species_option_switch
        logical, intent (out) :: ecoll_zeff
        
        character(20) :: species_option

        if (.not. proc0) return
        call set_default_species_options
        call read_input_file_species_options
        call check_inputs_species_options
    
    contains
        !------------------------ Default input parameters -----------------------
        subroutine set_default_species_options

            implicit none

            nspec = 2
            species_option = 'stella'
            read_profile_variation = .false.
            write_profile_variation = .false.
            ecoll_zeff = .false.

        end subroutine set_default_species_options

        !------------------------ Read input file species options ----------------
        subroutine read_input_file_species_options

            use file_utils, only: input_unit_exist, error_unit
            use text_options, only: text_option, get_option_value

            implicit none

            integer :: ierr

            type(text_option), dimension(4), parameter :: species_opts = (/ &
                                                    text_option('default', species_option_stella), &
                                                    text_option('stella', species_option_stella), &
                                                    text_option('input.profiles', species_option_inputprofs), &
                                                    text_option('euterpe', species_option_euterpe)/)

            namelist /species_options/ nspec, read_profile_variation, write_profile_variation, &
                                    species_option, ecoll_zeff

            in_file = input_unit_exist("species_options", dexist)
            if (dexist) read (unit=in_file, nml=species_options)
            ierr = error_unit()

            call get_option_value(species_option, species_opts, species_option_switch, &
                ierr, "species_option in namelist_species_options")

        end subroutine read_input_file_species_options

        !------------------------ Check inputs species options -------------------
        subroutine check_inputs_species_options

            use file_utils, only: runtype_option_switch, runtype_multibox
            use file_utils, only: error_unit
            use mp, only: mp_abort, job
            use parameters_physics, only: radial_variation

            implicit none

            integer :: ierr

            if (runtype_option_switch == runtype_multibox .and. (job /= 1) .and. radial_variation) then
                !will need to readjust the species parameters in the left/right boxes
                species_option_switch = species_option_multibox
            end if

            if (nspec < 1) then
                ierr = error_unit()
                write (unit=ierr, &
                    fmt="('Invalid nspec in species_options: ', i5)") nspec
                call mp_abort('Invalid nspec in species_options')
            end if

        end subroutine check_inputs_species_options

    end subroutine read_namelist_species_options

    !****************************************************************************
    !                          SPECIES OPTION : STELLA                          !
    !****************************************************************************
    subroutine read_namelist_species_stella (nspec, spec)

        use mp, only: proc0

        implicit none

        integer, intent (in) :: nspec
        type(spec_type), dimension(:), intent (out) :: spec

        real :: z, mass, dens, temp, tprim, fprim, d2ndr2, d2Tdr2, bess_fac
        character(20) :: type
        
        if (.not. proc0) return
        call read_input_file_species_stella

    contains
        !------------------------ Read input file species ------------------------
        subroutine read_input_file_species_stella

            use file_utils, only: input_unit_exist, error_unit
            use file_utils, only: get_indexed_namelist_unit
            use text_options, only: text_option, get_option_value

            implicit none

            integer :: ierr, unit, is
            type(text_option), dimension(9), parameter :: typeopts = (/ &
                                                    text_option('default', ion_species), &
                                                    text_option('ion', ion_species), &
                                                    text_option('electron', electron_species), &
                                                    text_option('e', electron_species), &
                                                    text_option('beam', slowing_down_species), &
                                                    text_option('fast', slowing_down_species), &
                                                    text_option('alpha', slowing_down_species), &
                                                    text_option('slowing-down', slowing_down_species), &
                                                    text_option('trace', tracer_species)/)

            namelist /species_parameters/ z, mass, dens, temp, &
                tprim, fprim, d2ndr2, d2Tdr2, bess_fac, type

            do is = 1, nspec
                call get_indexed_namelist_unit(unit, "species_parameters", is)
                call set_default_species_parameters
                read (unit=unit, nml=species_parameters)
                close (unit=unit)

                spec(is)%z = z
                spec(is)%mass = mass
                spec(is)%dens = dens
                spec(is)%temp = temp
                spec(is)%tprim = tprim
                spec(is)%fprim = fprim
                ! this is (1/n_s)*d^2 n_s / drho^2
                spec(is)%d2ndr2 = d2ndr2
                ! this is (1/T_s)*d^2 T_s / drho^2
                spec(is)%d2Tdr2 = d2Tdr2

                spec(is)%dens_psi0 = dens
                spec(is)%temp_psi0 = temp

                spec(is)%bess_fac = bess_fac

                ! if (present(write_profile_variation)) then
                !     write (filename, "(A,I1)") "specprof_", is
                !     open (1002, file=filename, status='unknown')
                !     write (1002, '(6e13.5)') dens, temp, fprim, tprim, d2ndr2, d2Tdr2
                !     close (1002)
                ! end if
                ! if (present(read_profile_variation)) then
                !     write (filename, "(A,I1)") "specprof_", is
                !     open (1002, file=filename, status='unknown')
                !     read (1002, '(6e13.5)') dens, temp, fprim, tprim, d2ndr2, d2Tdr2
                !     close (1002)

                !     dr = geo_surf%rhoc - geo_surf%rhoc_psi0
                !     spec(is)%dens = dens * (1.0 - dr * fprim)! + 0.5*dr**2*d2ndr2)
                !     spec(is)%temp = temp * (1.0 - dr * tprim)! + 0.5*dr**2*d2Tdr2)
                !     spec(is)%fprim = (fprim - dr * d2ndr2) * (dens / spec(is)%dens)
                !     spec(is)%tprim = (tprim - dr * d2Tdr2) * (temp / spec(is)%temp)
                !     !spec(is)%dens = 1.0
                !     !spec(is)%temp = 1.0
                ! end if
                ierr = error_unit()
                call get_option_value(type, typeopts, spec(is)%type, ierr, "type in species_parameters_x")
            end do

        end subroutine read_input_file_species_stella
        
        !------------------------ Set default species parameters ----------------
        subroutine set_default_species_parameters
        
            implicit none

            z = 1
            mass = 1.0
            dens = 1.0
            temp = 1.0
            tprim = -999.9
            fprim = -999.9
            d2ndr2 = 0.0
            d2Tdr2 = 0.0
            bess_fac = 1.0
            type = "default"

        end subroutine set_default_species_parameters

    end subroutine read_namelist_species_stella

    !****************************************************************************
    !                        SPECIES OPTIONS : euterpe                          !
    !****************************************************************************
    subroutine read_namelist_euterpe_parameters(nradii, data_file)

        use mp, only: proc0

        implicit none

        integer, intent(out) :: nradii
        character(*), intent(out) :: data_file
        
        if (.not. proc0) return
        call set_default_euterpe_parameters
        call read_input_file_euterpe_parameters
    
    contains
        !------------------------ Default input parameters -----------------------
        subroutine set_default_euterpe_parameters

            implicit none

            nradii = 1000
            data_file = 'euterpe.dat'

        end subroutine set_default_euterpe_parameters

        !------------------------ Read input file species options ----------------
        subroutine read_input_file_euterpe_parameters

            use file_utils, only: input_unit_exist, error_unit

            implicit none

            namelist /euterpe_parameters/ nradii, data_file

            in_file = input_unit_exist("euterpe_parameters", dexist)
            if (dexist) read (unit=in_file, nml=euterpe_parameters)

        end subroutine read_input_file_euterpe_parameters

    end subroutine read_namelist_euterpe_parameters

end module namelist_species
