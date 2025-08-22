module input_file_radial_variation

    implicit none

    public :: read_namelist_radial_variation

    public :: krook_option_default, krook_option_flat, &
           krook_option_linear, krook_option_exp, krook_option_exp_rev
    public :: mb_zf_option_default, mb_zf_option_skip_ky0, &
           mb_zf_option_zero_ky0, mb_zf_option_zero_fsa
    public :: lr_debug_option_default, lr_debug_option_L, lr_debug_option_R
           
    private

    integer, parameter:: krook_option_default = 2, &
                        krook_option_flat = 0, &
                        krook_option_linear = 1, &
                        krook_option_exp = 2, &
                        krook_option_exp_rev = 3
    integer, parameter :: mb_zf_option_default = 0, &
                         mb_zf_option_skip_ky0 = 1, &
                         mb_zf_option_zero_ky0 = 2, &
                         mb_zf_option_zero_fsa = 3
    integer, parameter:: lr_debug_option_default = 0, &
                        lr_debug_option_L = 1, &
                        lr_debug_option_R = 2

    integer :: in_file
    logical :: dexist

contains

    !****************************************************************************
    !                               RADIAL VARIATION                            !
    !****************************************************************************
    subroutine read_namelist_radial_variation(ky_solve_real, ky_solve_radial, &
                    include_pressure_variation, include_geometric_variation, &
                    smooth_zf, lr_debug_switch, krook_option_switch, mb_zf_option_switch, &
                    rk_step, nu_krook_mb, mb_debug_step, &
                    krook_exponent, krook_efold, phi_bound, phi_pow, &
                    use_dirichlet_bc, boundary_size, krook_size)

        use mp, only: proc0

        implicit none


        logical, intent (out) :: ky_solve_real
        integer, intent (out) :: ky_solve_radial
        logical, intent (out) :: include_pressure_variation
        logical, intent (out) :: include_geometric_variation
        logical, intent (out) :: smooth_zf
        integer, intent (out) :: lr_debug_switch, krook_option_switch, mb_zf_option_switch
        logical, intent (out) :: rk_step
        real, intent (out) ::  nu_krook_mb
        integer, intent (out) :: mb_debug_step
        real, intent (out) :: krook_exponent, krook_efold
        real, intent (out) :: phi_bound, phi_pow
        logical, intent (out) :: use_dirichlet_bc
        integer, intent (out) :: boundary_size, krook_size        

        character(30) :: zf_option, krook_option, lr_debug_option

        if (.not. proc0) return
        call set_default_parameters_radial_variation
        call read_input_file_radial_variation

    contains
      
    !------------------------ Default input parameters -----------------------
    subroutine set_default_parameters_radial_variation

        implicit none

        ky_solve_radial = 0 
        ky_solve_real =.false. 
        include_pressure_variation = .false. 
        include_geometric_variation = .false. 
        smooth_zf = .false.
        lr_debug_option = "default"
        krook_option = "default"
        zf_option = "default"
        rk_step = .false.
        nu_krook_mb = 0.0
        mb_debug_step = -1.0
        krook_exponent = 0.0
        krook_efold = 3.0
        phi_bound = 0.0
        phi_pow = 0.0
        use_dirichlet_bc = .false.
        boundary_size = 4.0
        krook_size = 0.0
        
    end subroutine set_default_parameters_radial_variation

    !---------------------------- Read input file ----------------------------
      subroutine read_input_file_radial_variation

        use file_utils, only: input_unit_exist, error_unit
        use text_options, only: text_option, get_option_value

        implicit none

        integer :: ierr
        type(text_option), dimension(3), parameter :: lr_db_opts = &
            (/text_option('default', lr_debug_option_default), &
                text_option('L', lr_debug_option_L), &
                text_option('R', lr_debug_option_R)/)
        type(text_option), dimension(5), parameter :: krook_opts = &
            (/text_option('default', krook_option_default), &
                text_option('flat', krook_option_flat), &
                text_option('linear', krook_option_linear), &
                text_option('exp', krook_option_exp), &
                text_option('exp_reverse', krook_option_exp_rev)/)
        type(text_option), dimension(4), parameter :: mb_zf_opts = &
            (/text_option('default', mb_zf_option_default), &
                text_option('skip_ky0', mb_zf_option_skip_ky0), &
                text_option('zero_ky0', mb_zf_option_zero_ky0), &
                text_option('zero_fsa', mb_zf_option_zero_fsa)/)

        namelist /multibox_parameters/ ky_solve_real, ky_solve_radial, &
                    include_pressure_variation, include_geometric_variation, &
                    smooth_zf, lr_debug_switch, krook_option_switch, mb_zf_option_switch, &
                    rk_step, nu_krook_mb, mb_debug_step, &
                    krook_exponent, krook_efold, phi_bound, phi_pow, &
                    use_dirichlet_bc, boundary_size, krook_size

        in_file = input_unit_exist("multibox_parameters", dexist)
        if (dexist) read (unit=in_file, nml=multibox_parameters)  

        ierr = error_unit()
        call get_option_value &
            (krook_option, krook_opts, krook_option_switch, &
             ierr, "krook_option in multibox_parameters")
        call get_option_value &
            (zf_option, mb_zf_opts, mb_zf_option_switch, &
             ierr, "zf_option in multibox_parameters")
        call get_option_value &
            (lr_debug_option, lr_db_opts, lr_debug_switch, &
             ierr, "lr_debug_option in multibox_parameters")

        if (krook_size > boundary_size) krook_size = boundary_size

      end subroutine read_input_file_radial_variation

   end subroutine read_namelist_radial_variation

end module input_file_radial_variation