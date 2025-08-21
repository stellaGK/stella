module input_file_neoclassical_input

    implicit none

    public :: read_namelist_neoclassical_input

    public :: neo_option_sfincs

    private

    integer, parameter :: neo_option_sfincs = 1

    ! These variables are used in every single subroutine, so make them global
    integer :: in_file
    logical :: dexist
contains

    !****************************************************************************
    !                             NEOCLASSICAL INPUT                            !
    !****************************************************************************
    subroutine read_namelist_neoclassical_input(include_neoclassical_terms, neo_option_switch, &
                nradii, drho)

        use mp, only: proc0

        implicit none

        integer, intent (out) :: neo_option_switch
        integer, intent (out) :: nradii
        real, intent (out) :: drho
        logical, intent(out) :: include_neoclassical_terms

        character(10) :: neo_option

        if (.not. proc0) return
        call set_default_parameters_neoclassical_input
        call read_input_file_neoclassical_input
        call check_inputs_neoclassical_input

    contains
      
        !------------------------ Default input parameters -----------------------
        subroutine set_default_parameters_neoclassical_input

            implicit none

            ! set to .true. to include neoclassical terms in GK equation
            include_neoclassical_terms = .false.
            ! number of radial points used for radial derivatives
            ! of neoclassical quantities
            nradii = 5.0
            ! spacing in rhoc between points used for radial derivatives
            drho = 0.01
            ! option for obtaining neoclassical distribution function and potential
            neo_option = 'sfincs'

        end subroutine set_default_parameters_neoclassical_input

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_neoclassical_input

            use file_utils, only: input_unit_exist, error_unit
            use text_options, only: text_option, get_option_value

            implicit none

            integer :: ierr
            type(text_option), dimension(2), parameter :: neoopts = (/ &
                                        text_option('default', neo_option_sfincs), &
                                        text_option('sfincs', neo_option_sfincs)/)
      
            namelist /neoclassical_input/ include_neoclassical_terms, neo_option, &
                nradii, drho

            in_file = input_unit_exist("neoclassical_input", dexist)
            if (dexist) read (unit=in_file, nml=neoclassical_input)

            ierr = error_unit()
            call get_option_value &
                (neo_option, neoopts, neo_option_switch, &
                ierr, "neo_option in neoclassical_input")

        end subroutine read_input_file_neoclassical_input

        subroutine check_inputs_neoclassical_input

            implicit none

            if (nradii /= 3 .and. nradii /= 5) then
                write (*, *) 'WARNING: only nradii of 3 or 5 is currently supported in neoclassical_input namelist'
                write (*, *) 'WARNING: forcing nradii=5'
                nradii = 5
            end if

        end subroutine check_inputs_neoclassical_input

    end subroutine read_namelist_neoclassical_input

end module input_file_neoclassical_input