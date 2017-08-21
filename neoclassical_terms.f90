module neoclassical_terms

  implicit none

  public :: init_neoclassical_terms
  public :: include_neoclassical_terms
  public :: finish_neoclassical_terms

  logical :: include_neoclassical_terms

  integer :: neo_option_switch
  integer, parameter :: neo_option_sfincs = 1

  logical :: neoinit = .false.

contains

  subroutine init_neoclassical_terms

    use sfincs_interface, only: get_neo_from_sfincs
    
    implicit none
    
    if (neoinit) return
    neoinit = .true.

    call read_parameters
    if (include_neoclassical_terms) then
       select case (neo_option_switch)
       case (neo_option_sfincs)
          call get_neo_from_sfincs
       end select
    end if

  end subroutine init_neoclassical_terms

  subroutine read_parameters

    use mp, only: proc0, broadcast
    use file_utils, only: error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value

    implicit none

    type (text_option), dimension (2), parameter :: neoopts = (/ &
         text_option('default', neo_option_sfincs), &
         text_option('sfincs', neo_option_sfincs) /)
    character (10) :: neo_option

    namelist /neoclassical_input/ include_neoclassical_terms, &
         neo_option

    logical :: exist
    integer :: ierr, in_file

    if (proc0) then
       include_neoclassical_terms = .false.
       neo_option = 'sfincs'

       in_file = input_unit_exist("neoclassical_input", exist)
       if (exist) read (unit=in_file, nml=neoclassical_input)

       ierr = error_unit()
       call get_option_value &
            (neo_option, neoopts, neo_option_switch, &
            ierr, "neo_option in neoclassical_input")
    end if

    call broadcast (include_neoclassical_terms)
    call broadcast (neo_option_switch)

  end subroutine read_parameters

  subroutine finish_neoclassical_terms

    implicit none

    neoinit = .false.

  end subroutine finish_neoclassical_terms

end module neoclassical_terms
