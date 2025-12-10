!> Some utilities for working with optional arguments
module optionals
  use iso_fortran_env, only: real32, real64, real128
  implicit none

  private

  public :: get_option_with_default

  interface get_option_with_default
     module procedure get_option_with_default_real32
     module procedure get_option_with_default_real64
     module procedure get_option_with_default_real128
     module procedure get_option_with_default_complex32
     module procedure get_option_with_default_complex64
     module procedure get_option_with_default_complex128
     module procedure get_option_with_default_integer
     module procedure get_option_with_default_logical
     module procedure get_option_with_default_character
  end interface get_option_with_default

contains

  !> Returns `option` if present or `default` if not.
  elemental function get_option_with_default_real32(option, default) result(option_out)
    implicit none
    real(kind=real32), intent(in), optional :: option
    real(kind=real32), intent(in) :: default
    real(kind=real32) :: option_out

    if( present(option) ) then
       option_out = option
    else
       option_out = default
    end if
  end function get_option_with_default_real32

  !> Returns `option` if present or `default` if not.
  elemental function get_option_with_default_real64(option, default) result(option_out)
    implicit none
    real(kind=real64), intent(in), optional :: option
    real(kind=real64), intent(in) :: default
    real(kind=real64) :: option_out

    if( present(option) ) then
       option_out = option
    else
       option_out = default
    end if
  end function get_option_with_default_real64

  !> Returns `option` if present or `default` if not.
  elemental function get_option_with_default_real128(option, default) result(option_out)
    implicit none
    real(kind=real128), intent(in), optional :: option
    real(kind=real128), intent(in) :: default
    real(kind=real128) :: option_out

    if( present(option) ) then
       option_out = option
    else
       option_out = default
    end if
  end function get_option_with_default_real128

  !> Returns `option` if present or `default` if not.
  elemental function get_option_with_default_complex32(option, default) result(option_out)
    implicit none
    complex(kind=real32), intent(in), optional :: option
    complex(kind=real32), intent(in) :: default
    complex(kind=real32) :: option_out

    if( present(option) ) then
       option_out = option
    else
       option_out = default
    end if
  end function get_option_with_default_complex32

  !> Returns `option` if present or `default` if not.
  elemental function get_option_with_default_complex64(option, default) result(option_out)
    implicit none
    complex(kind=real64), intent(in), optional :: option
    complex(kind=real64), intent(in) :: default
    complex(kind=real64) :: option_out

    if( present(option) ) then
       option_out = option
    else
       option_out = default
    end if
  end function get_option_with_default_complex64

  !> Returns `option` if present or `default` if not.
  elemental function get_option_with_default_complex128(option, default) result(option_out)
    implicit none
    complex(kind=real128), intent(in), optional :: option
    complex(kind=real128), intent(in) :: default
    complex(kind=real128) :: option_out

    if( present(option) ) then
       option_out = option
    else
       option_out = default
    end if
  end function get_option_with_default_complex128

  !> Returns `option` if present or `default` if not.
  elemental function get_option_with_default_integer(option, default) result(option_out)
    implicit none
    integer, intent(in), optional :: option
    integer, intent(in) :: default
    integer :: option_out

    if( present(option) ) then
       option_out = option
    else
       option_out = default
    end if
  end function get_option_with_default_integer

  !> Returns `option` if present or `default` if not.
  elemental function get_option_with_default_logical(option, default) result(option_out)
    implicit none
    logical, intent(in), optional :: option
    logical, intent(in) :: default
    logical :: option_out

    if( present(option) ) then
       option_out = option
    else
       option_out = default
    end if
  end function get_option_with_default_logical

  !> Returns `option` if present or `default` if not.
  !>
  !> @note Unlike the other members of the get_option_with_default
  !> interface, we cannot mark this routine as elemental due to
  !> the allocatable return type.
  pure function get_option_with_default_character(option, default) result(option_out)
    implicit none
    character(len=*), intent(in), optional :: option
    character(len=*), intent(in) :: default
    character(len=:), allocatable :: option_out

    if( present(option) ) then
       option_out = option
    else
       option_out = default
    end if
  end function get_option_with_default_character

end module optionals
