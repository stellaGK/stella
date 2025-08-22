module input_file_diagnostics

    implicit none

    public :: read_namelist_diagnostics
    public :: read_namelist_diagnostics_potential
    public :: read_namelist_diagnostics_omega
    public :: read_namelist_diagnostics_distribution
    public :: read_namelist_diagnostics_fluxes
    public :: read_namelist_diagnostics_moments

    private

    ! These variables are used in every single subroutine, so make them global
    integer :: in_file
    logical :: dexist
contains

    !****************************************************************************
    !                                DIAGNOSTICS                                !
    !****************************************************************************
    subroutine read_namelist_diagnostics (nwrite, navg, nsave, nc_mult, &
                    save_for_restart, write_all, write_all_time_traces, &
                    write_all_spectra_kxkyz, write_all_spectra_kxky, &
                    write_all_velocity_space, write_all_potential, &
                    write_all_omega, write_all_distribution, write_all_fluxes, &
                    write_all_moments)

        use mp, only: proc0

        implicit none

        integer, intent (out) :: nwrite, navg, nsave, nc_mult
        logical, intent (out) :: save_for_restart, write_all
        logical, intent (out) :: write_all_time_traces
        logical, intent (out) :: write_all_spectra_kxkyz
        logical, intent (out) :: write_all_spectra_kxky
        logical, intent (out) :: write_all_velocity_space
        logical, intent (out) :: write_all_potential
        logical, intent (out) :: write_all_omega
        logical, intent (out) :: write_all_distribution
        logical, intent (out) :: write_all_fluxes
        logical, intent (out) :: write_all_moments

        if (.not. proc0) return
        call set_default_parameters_diagnostics
        call read_input_file_diagnostics
        call check_inputs_diagnostics

    contains
      
        !------------------------ Default input parameters -----------------------
        subroutine set_default_parameters_diagnostics

            implicit none

            nwrite = 50.0
            navg = 50.0
            nsave = -1.0
            nc_mult = 1.0

            save_for_restart = .false.
            write_all = .false.
            write_all_time_traces = .true. 
            write_all_spectra_kxkyz = .false.
            write_all_spectra_kxky = .false.

            write_all_velocity_space = .false.
            write_all_potential = .false.
            write_all_omega = .false.
            write_all_distribution = .false.
            write_all_fluxes = .false.
            write_all_moments = .false.

        end subroutine set_default_parameters_diagnostics

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_diagnostics

            use file_utils, only: input_unit_exist, error_unit

            implicit none

            namelist /diagnostics/ nwrite, navg, nsave, nc_mult, &
                            save_for_restart, write_all, &
                            write_all_time_traces, write_all_spectra_kxkyz, &
                            write_all_spectra_kxky, write_all_velocity_space, &
                            write_all_potential, write_all_omega, &
                            write_all_distribution, write_all_fluxes, &
                            write_all_moments

            in_file = input_unit_exist("diagnostics", dexist)
            if (dexist) read (unit=in_file, nml=diagnostics)

        end subroutine read_input_file_diagnostics

        subroutine check_inputs_diagnostics

            implicit none

            if (write_all) then
                write_all_time_traces = .true.
                write_all_spectra_kxkyz = .true.
                write_all_spectra_kxky = .true.
                write_all_velocity_space = .true.
                write_all_potential = .true.
                write_all_omega = .true.
                write_all_distribution = .true.
                write_all_fluxes = .true.
                write_all_moments = .true.
            end if

            if (write_all_time_traces .or. write_all_spectra_kxky) then
                write_all_omega = .true.
            end if

        end subroutine check_inputs_diagnostics

    end subroutine read_namelist_diagnostics

    !****************************************************************************
    !                           DIAGNOSTICS POTENTIAL                           !
    !****************************************************************************
    subroutine read_namelist_diagnostics_potential (write_all_potential, &
                    write_all_time_traces, write_all_spectra_kxkyz, &
                    write_all_spectra_kxky, &
                    write_all_potential_time_traces, write_all_potential_spectra, &
                    write_phi2_vs_time, write_apar2_vs_time, &
                    write_bpar2_vs_time, write_phi_vs_kxkyz, &
                    write_apar_vs_kxkyz, write_bpar_vs_kxkyz, &
                    write_phi2_vs_kxky, write_apar2_vs_kxky, & 
                    write_bpar2_vs_kxky)

        use mp, only: proc0

        implicit none

        logical, intent (in) :: write_all_potential
        logical, intent (in) :: write_all_time_traces, write_all_spectra_kxkyz, &
                                write_all_spectra_kxky

        logical, intent (out) :: write_all_potential_time_traces, write_all_potential_spectra
        logical, intent (out) :: write_phi2_vs_time, write_apar2_vs_time, write_bpar2_vs_time
        logical, intent (out) :: write_phi_vs_kxkyz, write_apar_vs_kxkyz, write_bpar_vs_kxkyz
        logical, intent (out) :: write_phi2_vs_kxky, write_apar2_vs_kxky, write_bpar2_vs_kxky

        if (.not. proc0) return
        call set_default_parameters_diagnostics_potential
        call read_input_file_diagnostics_potential
        call check_inputs_diagnostics_potential

    contains
      
        !------------------------ Default input parameters -----------------------
        subroutine set_default_parameters_diagnostics_potential

            implicit none

            write_all_potential_time_traces = .false.
            write_all_potential_spectra = .false. 

            write_phi2_vs_time = .true.
            write_apar2_vs_time = .true.
            write_bpar2_vs_time = .true.

            write_phi_vs_kxkyz = .false.
            write_apar_vs_kxkyz = .false.
            write_bpar_vs_kxkyz = .false.

            write_phi2_vs_kxky = .false.
            write_apar2_vs_kxky = .false.
            write_bpar2_vs_kxky = .false.

        end subroutine set_default_parameters_diagnostics_potential

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_diagnostics_potential

            use file_utils, only: input_unit_exist, error_unit

            implicit none

            namelist /diagnostics_potential/ write_all_potential_spectra, &
                                write_all_potential_time_traces, &
                                write_phi2_vs_time, write_apar2_vs_time, &
                                write_bpar2_vs_time, write_phi_vs_kxkyz, &
                                write_apar_vs_kxkyz, write_bpar_vs_kxkyz, &
                                write_phi2_vs_kxky, write_apar2_vs_kxky, &
                                write_bpar2_vs_kxky
                                

            in_file = input_unit_exist("diagnostics_potential", dexist)
            if (dexist) read (unit=in_file, nml=diagnostics_potential)

        end subroutine read_input_file_diagnostics_potential

        subroutine check_inputs_diagnostics_potential

            use file_utils, only: error_unit

            implicit none

            if (write_all_time_traces) write_all_potential_time_traces = .true. 

            if (write_all_potential_time_traces) then
                write_phi2_vs_time = .true.
                write_apar2_vs_time = .true.
                write_bpar2_vs_time = .true.
            end if

            if (write_all_spectra_kxkyz .or. write_all_potential_spectra) then
                write_phi_vs_kxkyz = .true.
                write_apar_vs_kxkyz = .true.
                write_bpar_vs_kxkyz = .true.
            end if

            if (write_all_spectra_kxky .or. write_all_potential_spectra) then
                write_phi2_vs_kxky = .true.
                write_apar2_vs_kxky = .true.
                write_bpar2_vs_kxky = .true.
            end if

        end subroutine check_inputs_diagnostics_potential

    end subroutine read_namelist_diagnostics_potential

    !****************************************************************************
    !                           DIAGNOSTICS OMEGA                          !
    !****************************************************************************
    subroutine read_namelist_diagnostics_omega (write_all_omega, &
                    write_omega_vs_kxky, write_omega_avg_vs_kxky)

        use mp, only: proc0

        implicit none

        logical, intent (in) :: write_all_omega
        logical, intent (out) :: write_omega_vs_kxky
        logical, intent (out) :: write_omega_avg_vs_kxky

        if (.not. proc0) return
        call set_default_parameters_diagnostics_omega
        call read_input_file_diagnostics_omega
        call check_inputs_diagnostics_omega

    contains
      
        !------------------------ Default input parameters -----------------------
        subroutine set_default_parameters_diagnostics_omega

            use parameters_physics, only: include_nonlinear

            implicit none

            write_omega_vs_kxky = .not. include_nonlinear
            write_omega_avg_vs_kxky = .not. include_nonlinear

        end subroutine set_default_parameters_diagnostics_omega

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_diagnostics_omega

            use file_utils, only: input_unit_exist, error_unit

            implicit none

            namelist /diagnostics_omega/ write_omega_vs_kxky, write_omega_avg_vs_kxky

            in_file = input_unit_exist("diagnostics_omega", dexist)
            if (dexist) read (unit=in_file, nml=diagnostics_omega)

        end subroutine read_input_file_diagnostics_omega

        subroutine check_inputs_diagnostics_omega

            use file_utils, only: error_unit
            use parameters_physics, only: include_nonlinear

            implicit none

            if (write_all_omega) then
                write_omega_vs_kxky = .true.
                write_omega_avg_vs_kxky = .true.
            end if

            if (include_nonlinear) then
                write_omega_vs_kxky = .false.
                write_omega_avg_vs_kxky = .false.
            end if

        end subroutine check_inputs_diagnostics_omega

    end subroutine read_namelist_diagnostics_omega

    !****************************************************************************
    !                          DIAGNOSTICS DISTRIBUTION                         !
    !****************************************************************************
    subroutine read_namelist_diagnostics_distribution (write_all_distribution, &
                    write_all_spectra_kxkyz, &
                    write_all_velocity_space, &
                    write_g2_vs_vpamus, write_g2_vs_zvpas, &
                    write_g2_vs_zmus, write_g2_vs_kxkyzs, &
                    write_g2_vs_zvpamus, &
                    write_distribution_g, write_distribution_h, &
                    write_distribution_f)

        use mp, only: proc0

        implicit none

        logical, intent (in) :: write_all_distribution
        logical, intent (in) :: write_all_spectra_kxkyz
        logical, intent (in) :: write_all_velocity_space

        logical, intent (out) :: write_g2_vs_vpamus
        logical, intent (out) :: write_g2_vs_zvpas
        logical, intent (out) :: write_g2_vs_zmus
        logical, intent (out) :: write_g2_vs_kxkyzs
        logical, intent (out) :: write_g2_vs_zvpamus
        logical, intent (out) :: write_distribution_g
        logical, intent (out) :: write_distribution_h
        logical, intent (out) :: write_distribution_f

        if (.not. proc0) return
        call set_default_parameters_diagnostics_distribution
        call read_input_file_diagnostics_distribution
        call check_inputs_diagnostics_distribution

    contains
      
        !------------------------ Default input parameters -----------------------
        subroutine set_default_parameters_diagnostics_distribution

            implicit none

            write_g2_vs_vpamus = .false.
            write_g2_vs_zvpas = .false.
            write_g2_vs_zmus = .false.
            write_g2_vs_kxkyzs = .false.
            write_g2_vs_zvpamus = .false.
            write_distribution_g = .true.
            write_distribution_h = .false.
            write_distribution_f = .false.

        end subroutine set_default_parameters_diagnostics_distribution

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_diagnostics_distribution

            use file_utils, only: input_unit_exist, error_unit

            implicit none

            namelist /diagnostics_distribution/ write_g2_vs_vpamus, write_g2_vs_zvpas, &
                                write_g2_vs_zmus, write_g2_vs_kxkyzs, &
                                write_g2_vs_zvpamus, &
                                write_distribution_g, write_distribution_h, &
                                write_distribution_f

            in_file = input_unit_exist("diagnostics_distribution", dexist)
            if (dexist) read (unit=in_file, nml=diagnostics_distribution)

        end subroutine read_input_file_diagnostics_distribution

        subroutine check_inputs_diagnostics_distribution

            use file_utils, only: error_unit

            implicit none

            if (write_all_distribution) then
                write_g2_vs_vpamus = .true.
                write_g2_vs_zvpas = .true.
                write_g2_vs_zmus = .true.
                write_g2_vs_kxkyzs = .true.
                write_g2_vs_zvpamus = .true.
                write_distribution_g = .true.
                write_distribution_h = .true.
                write_distribution_f = .true.
            end if

            if (write_all_spectra_kxkyz) then 
                write_g2_vs_kxkyzs = .true.
                write_distribution_g = .true.
                write_distribution_h = .true.
                write_distribution_f = .true.
            end if 

            if (write_all_velocity_space) then 
                write_g2_vs_vpamus = .true.
                write_distribution_g = .true.
                write_distribution_h = .true.
                write_distribution_f = .true.
            end if


        end subroutine check_inputs_diagnostics_distribution

    end subroutine read_namelist_diagnostics_distribution

    !****************************************************************************
    !                          DIAGNOSTICS FLUXES                         !
    !****************************************************************************
    subroutine read_namelist_diagnostics_fluxes (write_all_fluxes, &
                            write_all_time_traces, write_all_spectra_kxkyz, &
                            write_all_spectra_kxky, &
                            flux_norm, write_fluxes_vs_time, &
                            write_radial_fluxes, write_fluxes_kxkyz, &
                            write_fluxes_kxky)

        use mp, only: proc0

        implicit none

        logical, intent (in) :: write_all_fluxes
        logical, intent (in) :: write_all_time_traces, write_all_spectra_kxkyz, &
                                write_all_spectra_kxky
        logical, intent (out) :: flux_norm
        logical, intent (out) :: write_fluxes_vs_time
        logical, intent (out) :: write_radial_fluxes
        logical, intent (out) :: write_fluxes_kxkyz
        logical, intent (out) :: write_fluxes_kxky

        if (.not. proc0) return
        call set_default_parameters_diagnostics_fluxes
        call read_input_file_diagnostics_fluxes
        call check_inputs_diagnostics_fluxes

    contains
      
        !------------------------ Default input parameters -----------------------
        subroutine set_default_parameters_diagnostics_fluxes

            use parameters_physics, only: radial_variation

            implicit none

            flux_norm = .true. 
            write_fluxes_vs_time = .true.
            write_radial_fluxes = radial_variation
            write_fluxes_kxkyz = .false.
            write_fluxes_kxky = .false.

        end subroutine set_default_parameters_diagnostics_fluxes

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_diagnostics_fluxes

        use file_utils, only: input_unit_exist, error_unit

            implicit none

            namelist /diagnostics_fluxes/ flux_norm, write_fluxes_vs_time, &
                                write_radial_fluxes, write_fluxes_kxkyz, &
                                write_fluxes_kxky

            in_file = input_unit_exist("diagnostics_fluxes", dexist)
            if (dexist) read (unit=in_file, nml=diagnostics_fluxes)

        end subroutine read_input_file_diagnostics_fluxes

        subroutine check_inputs_diagnostics_fluxes

            use file_utils, only: error_unit

            implicit none

            if (write_all_fluxes) then
                write_fluxes_vs_time = .true.
                write_radial_fluxes = .true.
                write_fluxes_kxkyz = .true.
                write_fluxes_kxky = .true.
            end if

            if (write_all_time_traces) then 
                write_fluxes_vs_time = .true.
            end if
            if (write_all_spectra_kxkyz) then 
                write_fluxes_kxkyz = .true.
            end if
            if (write_all_spectra_kxky) then
                write_fluxes_kxky = .true.
            end if

        end subroutine check_inputs_diagnostics_fluxes

    end subroutine read_namelist_diagnostics_fluxes

    !****************************************************************************
    !                            DIAGNOSTICS MOMENTS                            !
    !****************************************************************************
    subroutine read_namelist_diagnostics_moments (write_all_moments, &
                    write_moments, write_radial_moments)

        use mp, only: proc0

        implicit none

        logical, intent (in) :: write_all_moments
        logical, intent (out) :: write_moments
        logical, intent (out) :: write_radial_moments

        if (.not. proc0) return
        call set_default_parameters_diagnostics_moments
        call read_input_file_diagnostics_moments
        call check_inputs_diagnostics_moments

    contains
      
        !------------------------ Default input parameters -----------------------
        subroutine set_default_parameters_diagnostics_moments

            use parameters_physics, only: radial_variation

            implicit none

            write_moments = .false.
            write_radial_moments = radial_variation

        end subroutine set_default_parameters_diagnostics_moments

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_diagnostics_moments

        use file_utils, only: input_unit_exist, error_unit

            implicit none

            namelist /diagnostics_moments/ write_moments, write_radial_moments

            in_file = input_unit_exist("diagnostics_moments", dexist)
            if (dexist) read (unit=in_file, nml=diagnostics_moments)

        end subroutine read_input_file_diagnostics_moments

        subroutine check_inputs_diagnostics_moments

            use file_utils, only: error_unit

            implicit none

            if (write_all_moments) then
                write_moments = .true.
                write_radial_moments = .true.
            end if

        end subroutine check_inputs_diagnostics_moments

        end subroutine read_namelist_diagnostics_moments

end module input_file_diagnostics
