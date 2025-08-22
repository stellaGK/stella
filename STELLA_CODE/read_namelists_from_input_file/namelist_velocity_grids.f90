module namelist_velocity_grids

    implicit none

    public :: read_namelist_velocity_grids

    private

    ! These variables are used in every single subroutine, so make them global
    integer :: in_file
    logical :: dexist

contains

    !****************************************************************************
    !                               VELOCITY GRIDS                              !
    !****************************************************************************
    subroutine read_namelist_velocity_grids (nvgrid, nmu, vpa_max, vperp_max, &
        equally_spaced_mu_grid, conservative_wgts_vpa)

        use mp, only: proc0

        implicit none 

        integer, intent(out) :: nvgrid
        integer, intent(out) :: nmu
        real, intent(out) :: vpa_max, vperp_max
        logical, intent(out) :: equally_spaced_mu_grid, conservative_wgts_vpa

        call set_default_velocity_grids
        call read_input_file_velocity_grids

    contains
        !------------------------ Default input parameters -----------------------
        subroutine set_default_velocity_grids

            implicit none

            nvgrid = 24
            vpa_max = 3.0
            nmu = 12
            vperp_max = 3.0
            equally_spaced_mu_grid = .false.
            conservative_wgts_vpa = .false.

        end subroutine set_default_velocity_grids

        !------------------------ Read input file velocity grids parameters ----------------
        subroutine read_input_file_velocity_grids

            use file_utils, only: input_unit_exist

            implicit none

            namelist /velocity_grids/ nvgrid, nmu, vpa_max, vperp_max, &
                equally_spaced_mu_grid, conservative_wgts_vpa


            in_file = input_unit_exist("velocity_grids", dexist)
            if (dexist) read (unit=in_file, nml=velocity_grids)

        end subroutine read_input_file_velocity_grids

    end subroutine read_namelist_velocity_grids

end module namelist_velocity_grids