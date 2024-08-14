!###############################################################################
!###################### READ PARAMETERS FOR KXKY BOX GRID ######################
!###############################################################################
! Namelist: &parameters_kxky_grids_box
! These flags will allow you to toggle the algorithm choices in stella.
!###############################################################################

module parameters_kxky_grids_box
  
  public :: read_kxky_grids_box
  
    private

    logical :: initialised

  contains
    
    !======================================================================
    !================= READ PARAMETERS FOR KXKY BOX GRID ==================
    !======================================================================
    subroutine read_kxky_grids_box (nx, ny, ikx_max, naky_all, naky, nakx, nalpha, &
         x0, y0, jtwist, jtwistfac, phase_shift_angle, &
         centered_in_rho, randomize_phase_shift, periodic_variation)

        use mp, only: proc0, mp_abort
        use text_options, only: text_option, get_option_value
        use file_utils, only: input_unit, error_unit, input_unit_exist
    
        implicit none

        integer, intent (out) :: nx, ny, nalpha
        integer, intent (out) :: jtwist
        integer, intent (out) :: naky, nakx
        integer, intent (out) :: ikx_max, naky_all
        real, intent (out) :: jtwistfac
        real, intent (out) :: phase_shift_angle
        real, intent (out) :: x0, y0
        logical, intent (out) :: centered_in_rho, periodic_variation, randomize_phase_shift
        
        logical :: error = .false.
    
        if (initialised) return
        
        call read_default_box
        call read_input_file_box

        initialised = .true.
    
    contains

        !**********************************************************************
        !                        SET DEFAULT PARAMETERS                       !
        !**********************************************************************
        ! If not specified in the input file these are the default options that 
        ! will be set for all parameters under the namelist 
        ! &numerical'.
        !**********************************************************************

        subroutine read_default_box

            implicit none 

            nalpha = 1
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

        end subroutine read_default_box


        !**********************************************************************
        !                       READ GRID OPTION FOR KXK                      !
        !**********************************************************************
        ! Read which option to select for the kxky grid layouts
        !**********************************************************************
        subroutine read_input_file_box

            use file_utils, only: input_unit_exist
            use parameters_physics, only: full_flux_surface
            implicit none

            integer :: in_file
            logical :: exist

            namelist /parameters_kxky_grids_box/ nx, ny, jtwist, jtwistfac, x0, y0, &
            centered_in_rho, periodic_variation, &
            randomize_phase_shift, phase_shift_angle

            ! note that jtwist and y0 will possibly be modified
            ! later in init_kt_grids_box if they make it out
            ! of this subroutine with negative values
            ! it is necessary to wait until later to do this check
            ! because the values to which they may be set will
            ! depend on information from the geometry module,
            ! which itself may rely on ny from here (number of alphas)

            in_file = input_unit_exist("parameters_kxky_grids_box", exist)
            if (exist) read (in_file, nml=parameters_kxky_grids_box)

            call check_backwards_compatability_box
            
            !> Get the number of de-aliased modes in y and x, using reality to halve the number of ky modes
            naky = (ny - 1) / 3 + 1
            nakx = 2 * ((nx - 1) / 3) + 1

            !> Get the total number of ky values, including negative ky;
            !> this is approximately 2/3 ny because ny includes padding to avoid aliasing
            naky_all = 2 * naky - 1

            if (full_flux_surface) nalpha = ny
            
            !> Get the ikx index corresponding to kx_max 
            ikx_max = nakx / 2 + 1
            
            !> Get the total number of ky values, including negative ky; 
            !> this is approximately 2/3 ny because ny includes padding to avoid aliasing
            naky_all = 2 * naky - 1
                        
        end subroutine read_input_file_box

        !**********************************************************************
        !                    CHECK BACKWARDS COMPATIBILITY                    !
        !**********************************************************************
        ! Make sure stella either runs or aborts old names for variables or
        ! namelists are used
        !**********************************************************************

        subroutine check_backwards_compatability_box

            use file_utils, only: input_unit, input_unit_exist
            use parameters_numerical, only: print_extra_info_to_terminal
            implicit none
            
            logical :: old_nml_exist
            integer :: in_file

            namelist /kt_grids_box_parameters/ nx, ny, jtwist, jtwistfac, x0, y0, &
            centered_in_rho, periodic_variation, &
            randomize_phase_shift, phase_shift_angle
            
            !> rrm phase_shift angle

            in_file = input_unit_exist("kt_grids_box_parameters", old_nml_exist)
            if (old_nml_exist) then
               read (in_file, nml=kt_grids_box_parameters)
               if (print_extra_info_to_terminal) then 
                  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                  write(*,*) 'Please replace the namelist <kt_grids_box_parameters> in the'
                  write(*,*) 'input file with <parameters_kxky_grids_box> and include all'
                  write(*,*) 'variable names under this new namelist.'
               end if
            ! write(*,*) "Aborting in parameters_kxky_grids_box.f90. & 
            !         The namelist <kt_grids_box_parameters> does not exist. & 
            !         Please replace this with the title <parameters_kxky_grids_box>"
            ! call mp_abort("Aborting in parameters_kxky_grids_box.f90. &
            !         The namelist <kt_grids_box_parameters> does not exist. & 
            !         Please replace this with the title <parameters_kxky_grids_box>")
            end if 
            
        end subroutine check_backwards_compatability_box

    end subroutine read_kxky_grids_box

end module parameters_kxky_grids_box
