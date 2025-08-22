module input_file_kxky_grid

    implicit none
    public :: read_namelist_kxky_grid_option
    public :: read_namelist_kxky_grid_box
    public :: read_namelist_kxky_grid_range

    public :: kyspacing_linear, kyspacing_exponential
    public :: grid_option_range, grid_option_box
    private

    integer, parameter :: kyspacing_linear = 1, kyspacing_exponential = 2
    integer, parameter :: grid_option_range = 1, grid_option_box = 2
    ! Internal variables
    integer :: in_file
    logical :: dexist

contains

    !****************************************************************************
    !                           KXKY grid - GRID OPTION                         !
    !****************************************************************************
    subroutine read_namelist_kxky_grid_option(grid_option_switch)

        use mp, only: proc0

        implicit none

        integer, intent (out) :: grid_option_switch
        
        character(20) :: grid_option

        if (.not. proc0) return
        call read_default_kxky_grid_option
        call read_input_file_kxky_grid_option

    contains
        !------------------------ Default input parameters -----------------------
        subroutine read_default_kxky_grid_option
            
            implicit none

            grid_option = 'default'

        end subroutine read_default_kxky_grid_option
        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_kxky_grid_option 
        
            use file_utils, only: input_unit_exist, error_unit
            use text_options, only: text_option, get_option_value

            implicit none

            integer :: ierr
      
            type(text_option), dimension(5), parameter :: grid_options = &
                (/text_option('default', grid_option_range), &
                text_option('range', grid_option_range), &
                text_option('box', grid_option_box), &
                text_option('annulus', grid_option_box), &
                text_option('nonlinear', grid_option_box)/)
        
    
            namelist /kxky_grid_option/ grid_option
      
            in_file = input_unit_exist("kxky_grid_option", dexist)
            if (dexist) read (unit=in_file, nml=kxky_grid_option)
            
            ierr = error_unit()
            call get_option_value(grid_option, grid_options, grid_option_switch, &
                ierr, "grid_option in input_file_kxky_grid", .true.)
      
        end subroutine read_input_file_kxky_grid_option

    end subroutine read_namelist_kxky_grid_option

    !****************************************************************************
    !                           KXKY grid - BOX OPTION                          !
    !****************************************************************************
    subroutine read_namelist_kxky_grid_box (nx, ny, ikx_max, naky_all, naky, nakx, nalpha, &
         x0, y0, jtwist, jtwistfac, phase_shift_angle, &
         centered_in_rho, randomize_phase_shift, periodic_variation, reality)

        use mp, only: proc0
    
        implicit none

        integer, intent (out) :: nx, ny, nalpha
        integer, intent (out) :: jtwist
        integer, intent (out) :: naky, nakx
        integer, intent (out) :: ikx_max, naky_all
        real, intent (out) :: jtwistfac
        real, intent (out) :: phase_shift_angle
        real, intent (out) :: x0, y0
        logical, intent (out) :: centered_in_rho, periodic_variation, randomize_phase_shift
        logical, intent (out) :: reality
    
        if (.not. proc0) return
        call read_default_kxky_grid_box
        call read_input_file_kxky_grid_box

    contains

        !------------------------ Default input parameters -----------------------
        subroutine read_default_kxky_grid_box

            implicit none 

            nx = 1
            ny = 1
            jtwist = -1
            jtwistfac = 1.
            phase_shift_angle = 0.
            x0 = -1.0
            y0 = -1.0
            nalpha = 1
            centered_in_rho = .true.
            randomize_phase_shift = .false.
            periodic_variation = .false.

        end subroutine read_default_kxky_grid_box

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_kxky_grid_box

            use file_utils, only: input_unit_exist, error_unit
            use text_options, only: text_option, get_option_value
            use parameters_physics, only: full_flux_surface
            
            implicit none

            namelist /kxky_grid_box/ nx, ny, jtwist, jtwistfac, x0, y0, &
            centered_in_rho, periodic_variation, &
            randomize_phase_shift, phase_shift_angle

            ! note that jtwist and y0 will possibly be modified
            ! later in init_kxky_grid_box if they make it out
            ! of this subroutine with negative values
            ! it is necessary to wait until later to do this check
            ! because the values to which they may be set will
            ! depend on information from the geometry module,
            ! which itself may rely on ny from here (number of alphas)

            in_file = input_unit_exist("kxky_grid_box", dexist)
            if (dexist) read (unit=in_file, nml=kxky_grid_box)
            
            !> Get the number of de-aliased modes in y and x, using reality to halve the number of ky modes
            naky = (ny - 1) / 3 + 1
            nakx = 2 * ((nx - 1) / 3) + 1

            if (full_flux_surface) nalpha = ny
            
            !> Get the ikx index corresponding to kx_max 
            ikx_max = nakx / 2 + 1

            reality = .true.
            !> Get the total number of ky values, including negative ky; 
            !> this is approximately 2/3 ny because ny includes padding to avoid aliasing
            naky_all = 2 * naky - 1
                        
        end subroutine read_input_file_kxky_grid_box
    
    end subroutine read_namelist_kxky_grid_box

    !****************************************************************************
    !                           KXKY grid - RANGE OPTION                        !
    !****************************************************************************
    subroutine read_namelist_kxky_grid_range (nalpha, naky, nakx, aky_min, aky_max, &
              akx_min, akx_max, theta0_min, theta0_max, &
              kyspacing_option_switch, phase_shift_angle, ikx_max, naky_all)

        use mp, only: mp_abort
    
        implicit none

        integer, intent (out) :: nalpha, naky, nakx
        real, intent (out) :: aky_min, aky_max, akx_min, akx_max
        real, intent (out) :: theta0_min, theta0_max
        integer, intent (out) :: kyspacing_option_switch
        real, intent (out) :: phase_shift_angle
        integer, intent (out) :: ikx_max, naky_all
        
        character(20) :: kyspacing_option

        call read_default_kxky_grid_range
        call read_input_file_kxky_grid_range

    contains
        !------------------------ Default input parameters -----------------------
        subroutine read_default_kxky_grid_range

          implicit none

          kyspacing_option = 'default'

          nalpha = 1
          naky = 1
          nakx = 1
          aky_min = 0.0
          aky_max = 0.0
          !> set these to be nonsense values
          !> so can check later if they've been set
          akx_min = 0.0
          akx_max = -1.0
          theta0_min = 0.0
          theta0_max = -1.0
          phase_shift_angle = 0.
      
        end subroutine read_default_kxky_grid_range

        !---------------------------- Read input file ----------------------------
        subroutine read_input_file_kxky_grid_range

            use file_utils, only: input_unit_exist, error_unit
            use text_options, only: text_option, get_option_value
            use parameters_physics, only: full_flux_surface
      
            implicit none
      
            type(text_option), dimension(3), parameter :: kyspacingopts = &
                 (/text_option('default', kyspacing_linear), &
                 text_option('linear', kyspacing_linear), &
                 text_option('exponential', kyspacing_exponential)/)
      
            integer :: ierr

            namelist /kxky_grid_range/ naky, nakx, &
                 aky_min, aky_max, theta0_min, theta0_max, akx_min, akx_max, kyspacing_option

            ! note that jtwist and y0 will possibly be modified
            ! later in init_kt_grid_range if they make it out
            ! of this subroutine with negative values
            ! it is necessary to wait until later to do this check
            ! because the values to which they may be set will
            ! depend on information from the geometry module,
            ! which itself may rely on ny from here (number of alphas)
            in_file = input_unit_exist("kxky_grid_range", dexist)
            if (dexist) read (unit=in_file, nml=kxky_grid_range)
            
            if (full_flux_surface) then
                write (*, *) '!!! ERROR !!!'
                write (*, *) 'kt_grid "range" option is not supported for full_flux_surface = T. aborting'
                write (*, *) '!!! ERROR !!!'
                call mp_abort('kt_grid "range" option is not supported for full_flux_surface = T. aborting')
            end if
            
            ierr = error_unit()
            call get_option_value(kyspacing_option, kyspacingopts, kyspacing_option_switch, &
                ierr, "kyspacing_option in kxky_grid_range", .true.)

            naky_all = naky
            ikx_max = nakx
            
        end subroutine read_input_file_kxky_grid_range
    
    end subroutine read_namelist_kxky_grid_range

end module input_file_kxky_grid
