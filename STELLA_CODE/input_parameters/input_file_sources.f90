module input_file_sources

    implicit none

    public :: read_namelist_sources

    public :: source_option_none, source_option_krook, &
             source_option_projection
    private
    
    integer, parameter :: source_option_none = 1, &
                         source_option_krook = 2, &
                         source_option_projection = 3

    ! These variables are used in every single subroutine, so make them global
    integer :: in_file
    logical :: dexist

contains

    !****************************************************************************
    !                                  SOURCES                                  !
    !****************************************************************************
    subroutine read_namelist_sources(source_option_switch, nu_krook, tcorr_source, &
                    ikxmax_source, krook_odd, exclude_boundary_regions, &
                    tcorr_source_qn, exclude_boundary_regions_qn, from_zero, &
                    conserve_momentum, conserve_density)

        use mp, only: proc0

        implicit none

        integer, intent (out) :: source_option_switch
        real, intent (out) :: nu_krook, tcorr_source
        integer, intent (out):: ikxmax_source
        logical, intent (out) :: krook_odd, exclude_boundary_regions
        real, intent (out) :: tcorr_source_qn
        logical, intent (out) :: exclude_boundary_regions_qn
        logical, intent (out) :: from_zero
        logical, intent (out) :: conserve_momentum, conserve_density
        
        character(30) :: source_option

        if (.not. proc0) return
        call set_default_parameters_sources
        call read_input_file_sources
        call check_inputs_sources

    contains
      
        !------------------------ Default input parameters -----------------------
        subroutine set_default_parameters_sources

            use parameters_physics, only: radial_variation
            use kxky_grid_parameters, only: periodic_variation

            implicit none

            source_option = 'none'
            nu_krook = 0.05
            tcorr_source = 0.02
            tcorr_source_qn = 0.0
            ikxmax_source = 1.0
            if (periodic_variation) ikxmax_source = 2 ! kx=0 and kx=1
            krook_odd = .true.! damp only the odd mode that can affect profiles
            exclude_boundary_regions = radial_variation .and. .not. periodic_variation
            ! exclude_boundary_regions_qn = exclude_boundary_regions
            from_zero = .true.
            conserve_momentum = .false.
            conserve_density = .false.

        end subroutine set_default_parameters_sources

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_sources

            use file_utils, only: input_unit_exist, error_unit
            use text_options, only: text_option, get_option_value

            implicit none

            integer :: ierr
            type(text_option), dimension(4), parameter :: sourceopts = &
                                (/text_option('default', source_option_none), &
                                    text_option('none', source_option_none), &
                                    text_option('krook', source_option_krook), &
                                    text_option('projection', source_option_projection)/)

            namelist /sources/ source_option, nu_krook, tcorr_source, &
                            ikxmax_source, krook_odd, exclude_boundary_regions, &
                            tcorr_source_qn, exclude_boundary_regions_qn, from_zero, &
                            conserve_momentum, conserve_density

            in_file = input_unit_exist("sources", dexist)
            if (dexist) read (unit=in_file, nml=sources)

            ierr = error_unit()
            call get_option_value &
                (source_option, sourceopts, source_option_switch, &
                ierr, "source_option in sources")

        end subroutine read_input_file_sources

        subroutine check_inputs_sources

            implicit none

            if (tcorr_source_qn < 0) tcorr_source_qn = tcorr_source

        end subroutine check_inputs_sources

    end subroutine read_namelist_sources

end module input_file_sources