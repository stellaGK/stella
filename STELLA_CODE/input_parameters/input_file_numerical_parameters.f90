module input_file_numerical_parameters

    implicit none

    public :: read_namelist_time_trace_options
    public :: read_namelist_time_step
    public :: read_namelist_numerical_algorithms
    public :: read_namelist_numerical_upwinding_for_derivatives
    public :: read_namelist_numerical_extra ! To be moved

    public :: delt_option_hand, delt_option_auto
    public :: explicit_algorithm_rk3, explicit_algorithm_rk2, &
                explicit_algorithm_rk4, explicit_algorithm_euler

    private

    integer, parameter :: delt_option_hand = 1, delt_option_auto = 2
    integer, parameter :: explicit_algorithm_rk3 = 1, &
        explicit_algorithm_rk2 = 2, &
        explicit_algorithm_rk4 = 3, &
        explicit_algorithm_euler = 4

    ! Internal variables
    integer :: in_file
    logical :: dexist

contains

    !****************************************************************************
    !                             TIME TRACE OPTIONS                            !
    !****************************************************************************
    subroutine read_namelist_time_trace_options(nstep, tend, autostop, avail_cpu_time)

        use mp, only: proc0

        implicit none

        integer, intent (out) :: nstep 
        real, intent (out) :: tend 
        logical, intent (out) :: autostop
        real, intent (out) :: avail_cpu_time

        if (.not. proc0) return
        call set_default_parameters_time_trace_options
        call read_input_file_time_trace_options

    contains
      
        !------------------------ Default input parameters -----------------------
        subroutine set_default_parameters_time_trace_options

            implicit none

            ! By default
            nstep = -1.0
            tend = -1.0
            autostop = .true.
            avail_cpu_time = 10000000000.0

        end subroutine set_default_parameters_time_trace_options

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_time_trace_options

            use file_utils, only: input_unit_exist
            implicit none

            namelist /time_trace_options/ nstep, tend, autostop, avail_cpu_time
            in_file = input_unit_exist("time_trace_options", dexist)
            if (dexist) read (unit=in_file, nml=time_trace_options)

        end subroutine read_input_file_time_trace_options

    end subroutine read_namelist_time_trace_options

    !****************************************************************************
    !                             TIME STEP OPTIONS                             !
    !****************************************************************************
    subroutine read_namelist_time_step(delt, delt_option_switch, &
                                       delt_max, delt_min, &
                                       cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower)

        use mp, only: proc0

        implicit none

        real, intent (out) :: delt
        integer, intent (out) :: delt_option_switch
        real, intent (out) :: delt_max, delt_min
        real, intent (out) :: cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower

        character(30) :: delt_option

        if (.not. proc0) return
        call set_default_parameters_time_step
        call read_input_file_time_step

    contains
      
        !------------------------ Default input parameters -----------------------
        subroutine set_default_parameters_time_step

            implicit none

            ! By default
            delt = 0.03
            delt_option = "default"
            delt_max = -1.0
            delt_min = 1e-10
            cfl_cushion_upper = 0.5
            cfl_cushion_middle = 0.25
            cfl_cushion_lower = 1e-05

        end subroutine set_default_parameters_time_step

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_time_step

            use file_utils, only: input_unit_exist, error_unit
            use text_options, only: text_option, get_option_value

            implicit none

            integer :: ierr 
            type(text_option), dimension(3), parameter :: deltopts = &
                (/text_option('default', delt_option_auto), &
                text_option('set_by_hand', delt_option_hand), &
                text_option('check_restart', delt_option_auto)/)

            !----------------------------------------------------------------------    
            namelist /time_step/ delt, delt_option, &
                delt_max, delt_min, &
                cfl_cushion_upper, cfl_cushion_middle, cfl_cushion_lower

            in_file = input_unit_exist("time_step", dexist)
            if (dexist) read (unit=in_file, nml=time_step)

            ierr = error_unit()
            call get_option_value &
                (delt_option, deltopts, delt_option_switch, ierr, &
                "delt_option in input_file_numerical_parameters.f90")

        end subroutine read_input_file_time_step

    end subroutine read_namelist_time_step

    !****************************************************************************
    !                           NUMERICAL ALGORITHMS                            !
    !****************************************************************************
    subroutine read_namelist_numerical_algorithms(explicit_algorithm_switch, flip_flop, &
                                                  stream_implicit, stream_iterative_implicit, &
                                                  stream_matrix_inversion, driftkinetic_implicit, &
                                                  mirror_implicit, mirror_semi_lagrange, &
                                                  mirror_linear_interp, drifts_implicit, &
                                                  fully_implicit, fully_explicit, &
                                                  maxwellian_inside_zed_derivative, &
                                                  use_deltaphi_for_response_matrix, & 
                                                  split_parallel_dynamics, &
                                                  maxwellian_normalization)

        use mp, only: proc0

        implicit none
        integer, intent(out) :: explicit_algorithm_switch
        logical, intent(out) :: flip_flop
        logical, intent(out) :: stream_implicit, stream_iterative_implicit, stream_matrix_inversion, driftkinetic_implicit
        logical, intent(out) :: mirror_implicit, mirror_semi_lagrange, mirror_linear_interp
        logical, intent(out) :: drifts_implicit, fully_implicit, fully_explicit
        logical, intent(out) :: maxwellian_inside_zed_derivative, use_deltaphi_for_response_matrix
        logical, intent(out) :: split_parallel_dynamics, maxwellian_normalization

        character(30) :: explicit_algorithm

        if (.not. proc0) return
        call set_default_parameters_numerical_algorithms
        call read_input_file_numerical_algorithms

    contains
      
        !------------------------ Default input parameters -----------------------
        subroutine set_default_parameters_numerical_algorithms

            implicit none

            ! By default
            explicit_algorithm = "rk3"
            flip_flop = .false.  
            stream_implicit = .true.
            stream_iterative_implicit = .false.
            stream_matrix_inversion = .false.
            driftkinetic_implicit = .false.
            mirror_implicit = .true.
            mirror_semi_lagrange = .true.
            mirror_linear_interp = .false.
            drifts_implicit = .false.
            fully_implicit = .false.
            fully_explicit = .false.
            maxwellian_inside_zed_derivative = .false.
            use_deltaphi_for_response_matrix = .false.
            split_parallel_dynamics = .false.
            maxwellian_normalization = .false.

        end subroutine set_default_parameters_numerical_algorithms

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_numerical_algorithms

            use file_utils, only: input_unit_exist, error_unit
            use text_options, only: text_option, get_option_value

            implicit none

            integer :: ierr
            type(text_option), dimension(5), parameter :: explicitopts = &
              (/text_option('default', explicit_algorithm_rk3), &
              text_option('rk3', explicit_algorithm_rk3), &
              text_option('rk2', explicit_algorithm_rk2), &
              text_option('rk4', explicit_algorithm_rk4), &
              text_option('euler', explicit_algorithm_euler)/)
              
            namelist /numerical_algorithms/ explicit_algorithm, flip_flop, &
                                    stream_implicit, stream_iterative_implicit, &
                                    stream_matrix_inversion, driftkinetic_implicit, &
                                    mirror_implicit, mirror_semi_lagrange, &
                                    mirror_linear_interp, drifts_implicit, &
                                    fully_implicit, fully_explicit, &
                                    maxwellian_inside_zed_derivative, &
                                    use_deltaphi_for_response_matrix, & 
                                    split_parallel_dynamics, &
                                    maxwellian_normalization

            in_file = input_unit_exist("numerical_algorithms", dexist)
            if (dexist) read (unit=in_file, nml=numerical_algorithms)

            call get_option_value &
                (explicit_algorithm, explicitopts, explicit_algorithm_switch, &
                ierr, "explicit_algorithm in numerical_parameters")

        end subroutine read_input_file_numerical_algorithms

    end subroutine read_namelist_numerical_algorithms

    !****************************************************************************
    !                NUMERICAL UPWINDING FOR DERIVATIVES OPTIONS                 !
    !****************************************************************************
    subroutine read_namelist_numerical_upwinding_for_derivatives(time_upwind, zed_upwind, vpa_upwind) 

        use mp, only: proc0

        implicit none

        real, intent (out) :: time_upwind, zed_upwind, vpa_upwind

        if (.not. proc0) return
        call set_default_parameters_numerical_upwinding_for_derivatives
        call read_input_file_numerical_upwinding_for_derivatives

    contains
        
        !------------------------ Default input parameters -----------------------
        subroutine set_default_parameters_numerical_upwinding_for_derivatives

            implicit none

            ! By default
            time_upwind = 0.02
            zed_upwind = 0.02
            vpa_upwind = 0.02

        end subroutine set_default_parameters_numerical_upwinding_for_derivatives

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_numerical_upwinding_for_derivatives

            use file_utils, only: input_unit_exist
            implicit none

            namelist /numerical_upwinding_for_derivatives/ time_upwind, zed_upwind, vpa_upwind

            in_file = input_unit_exist("numerical_upwinding_for_derivatives", dexist)
            if (dexist) read (unit=in_file, nml=numerical_upwinding_for_derivatives)

        end subroutine read_input_file_numerical_upwinding_for_derivatives

    end subroutine read_namelist_numerical_upwinding_for_derivatives

    !****************************************************************************
    !                           NUMERICAL EXTRA OPTIONS                         !
    !****************************************************************************

    subroutine read_namelist_numerical_extra(nitt, fphi, rng_seed)

        use mp, only: proc0

        implicit none
        integer, intent (out)  :: nitt
        real, intent (out)  :: fphi
        integer, intent (out)  :: rng_seed

        if (.not. proc0) return
        call set_default_parameters_numerical_extra
        call read_input_file_numerical_extra

    contains
        
        subroutine set_default_parameters_numerical_extra

            implicit none

            ! By default
            nitt = 1
            fphi = 1.0
            rng_seed = -1

        end subroutine set_default_parameters_numerical_extra

        subroutine read_input_file_numerical_extra

            use file_utils, only: input_unit_exist, error_unit
            use text_options, only: text_option, get_option_value

            implicit none

            namelist /numerical_extra/ nitt, fphi, rng_seed

            in_file = input_unit_exist("numerical_extra", dexist)
            if (dexist) read (unit=in_file, nml=numerical_extra)
            
        end subroutine read_input_file_numerical_extra

    end subroutine read_namelist_numerical_extra

end module input_file_numerical_parameters