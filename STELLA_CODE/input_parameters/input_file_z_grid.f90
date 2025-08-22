module input_file_z_grid

    implicit none

    public :: read_namelist_z_grid
    public :: read_namelist_z_boundary_condition

    public :: boundary_option_zero, boundary_option_self_periodic
    public :: boundary_option_linked, boundary_option_linked_stellarator

    private

    integer, parameter :: boundary_option_zero = 1, &
                         boundary_option_self_periodic = 2, &
                         boundary_option_linked = 3, &
                         boundary_option_linked_stellarator = 4
    ! These variables are used in every single subroutine, so make them global
    integer :: in_file
    logical :: dexist
contains
    !****************************************************************************
    !                                   Z GRID                                  !
    !****************************************************************************
    subroutine read_namelist_z_grid(nzed, nperiod, ntubes, zed_equal_arc)

        use mp, only: proc0

        implicit none

        integer, intent (out) :: nzed
        integer, intent (out) :: nperiod, ntubes
        logical, intent (out) :: zed_equal_arc

        if (.not. proc0) return
        call set_default_z_grid
        call read_input_file_z_grid
        call check_inputs_z_grid

    contains
        !------------------------ Default input parameters -----------------------
        subroutine set_default_z_grid

            implicit none 

            nzed = 24
            nperiod = 1
            ntubes = 1
            ! if zed_equal_arc = T, then zed is chosen to be arc length
            ! if zed_equal_arc = F, then zed is poloidal (axisymmetric)
            ! or zeta (toroidal) angle
            zed_equal_arc = .false.

        end subroutine set_default_z_grid

        subroutine read_input_file_z_grid

            use file_utils, only: input_unit_exist

            implicit none

            namelist /z_grid/ nzed, nperiod, ntubes, zed_equal_arc

            in_file = input_unit_exist("z_grid", dexist)
            if (dexist) read (unit=in_file, nml=z_grid)

        end subroutine read_input_file_z_grid

        subroutine check_inputs_z_grid

            use parameters_physics, only: full_flux_surface

            implicit none

            ! Make sure <nzed> is an even integer, otherwise the potential explodes
            if (MOD(nzed,2) .eq. 1) nzed = nzed + 1

            ! force use of equal arc grid to ensure gradpar alpha-independent
            ! necessary to obtain efficient numerical solution of parallel streaming
            if (full_flux_surface) zed_equal_arc = .true.
        
        end subroutine check_inputs_z_grid

    end subroutine read_namelist_z_grid

    !****************************************************************************
    !                           Z GRID BOUNDARY OPTIONS                         !
    !****************************************************************************
    subroutine read_namelist_z_boundary_condition (boundary_option_switch, shat_zero, grad_x_grad_y_zero, & 
        dkx_over_dky) ! put rhostar here? 

        use mp, only: proc0

        implicit none

        integer, intent (out) :: boundary_option_switch
        real, intent (out) :: shat_zero, grad_x_grad_y_zero, dkx_over_dky

        character (20) :: boundary_option

        if (.not. proc0) return
        call set_default_z_boundary_condition
        call read_input_file_z_boundary_condition

    contains
        !------------------------ Default input parameters -----------------------
        subroutine set_default_z_boundary_condition

            implicit none
            boundary_option = 'default'
            ! set minimum shat value below which we assume periodic BC
            shat_zero = 1.e-5
            ! set the minimum nabla x . nabla value at the end of the FT which we assume
            ! periodic BC instead of the stellarator symmetric ones
            grad_x_grad_y_zero = 1.e-5
        
            ! set the ratio between dkx and dky, assuming jtwist = 1.
            ! if it is < 0, the code will just use the nfield_periods in the input file
            dkx_over_dky = -1

        end subroutine set_default_z_boundary_condition

        subroutine read_input_file_z_boundary_condition
            use file_utils, only: input_unit_exist, error_unit
            use text_options, only: text_option, get_option_value

            implicit none 

            integer :: ierr 
            type(text_option), dimension(7), parameter :: boundaryopts = &
                            (/text_option('default', boundary_option_zero), &
                                text_option('zero', boundary_option_zero), &
                                text_option('unconnected', boundary_option_zero), &
                                text_option('self-periodic', boundary_option_self_periodic), &
                                text_option('periodic', boundary_option_self_periodic), &
                                text_option('linked', boundary_option_linked), &
                                text_option('stellarator', boundary_option_linked_stellarator)/)

            namelist /z_boundary_condition/ boundary_option, shat_zero, &
                grad_x_grad_y_zero, dkx_over_dky

            in_file = input_unit_exist("z_boundary_condition", dexist)
            if (dexist) read (unit=in_file, nml=z_boundary_condition)

            ierr = error_unit()
            call get_option_value(boundary_option, boundaryopts, boundary_option_switch, &
                ierr, "boundary_option in input_file_z_grid")

        end subroutine read_input_file_z_boundary_condition

    end subroutine read_namelist_z_boundary_condition

end module input_file_z_grid