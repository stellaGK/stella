# include "define.inc"

module shmem
  implicit none
  private

  public :: shmem_available
# ifdef SHMEM
  logical, parameter :: shmem_available = .true.
#else
  logical, parameter :: shmem_available = .false.

  public :: shmem_integer_put, shmem_real_put, shmem_logical_put
  public :: shmem_character_put
  public :: shmem_integer_get, shmem_real_get, shmem_logical_get
  public :: shmem_character_get
  public :: shmem_wait

  interface shmem_integer_put
     module procedure shmem_integer_put_array, shmem_integer_put_scalar
  end interface

  interface shmem_real_put
     module procedure shmem_real_put_array, shmem_real_put_scalar
  end interface

  interface shmem_logical_put
     module procedure shmem_logical_put_array, shmem_logical_put_scalar
  end interface

  interface shmem_integer_get
     module procedure shmem_integer_get_array, shmem_integer_get_scalar
  end interface

  interface shmem_real_get
     module procedure shmem_real_get_array, shmem_real_get_scalar
  end interface

  interface shmem_logical_get
     module procedure shmem_logical_get_array, shmem_logical_get_scalar
  end interface

contains

  subroutine shmem_integer_put_array (target, src, len, pe)
    implicit none
    integer, dimension (:), intent (in out) :: target
    integer, dimension (:), intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_integer_put")
  end subroutine shmem_integer_put_array

  subroutine shmem_integer_put_scalar (target, src, len, pe)
    implicit none
    integer, intent (in out) :: target
    integer, intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_integer_put")
  end subroutine shmem_integer_put_scalar

  subroutine shmem_real_put_array (target, src, len, pe)
    implicit none
    real, dimension (:), intent (in out) :: target
    real, dimension (:), intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_real_put")
  end subroutine shmem_real_put_array

  subroutine shmem_real_put_scalar (target, src, len, pe)
    implicit none
    real, intent (in out) :: target
    real, intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_real_put")
  end subroutine shmem_real_put_scalar

  subroutine shmem_complex_put_array (target, src, len, pe)
    implicit none
    complex, dimension (:), intent (in out) :: target
    complex, dimension (:), intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_complex_put")
  end subroutine shmem_complex_put_array

  subroutine shmem_complex_put_scalar (target, src, len, pe)
    implicit none
    complex, intent (in out) :: target
    complex, intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_complex_put")
  end subroutine shmem_complex_put_scalar

  subroutine shmem_logical_put_array (target, src, len, pe)
    implicit none
    logical, dimension (:), intent (in out) :: target
    logical, dimension (:), intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_logical_put")
  end subroutine shmem_logical_put_array

  subroutine shmem_logical_put_scalar (target, src, len, pe)
    implicit none
    logical, intent (in out) :: target
    logical, intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_logical_put")
  end subroutine shmem_logical_put_scalar

  subroutine shmem_character_put (target, src, len, pe)
    implicit none
    character(*), intent (in out) :: target
    character(*), intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_character_put")
  end subroutine shmem_character_put

  subroutine shmem_integer_get_array (target, src, len, pe)
    implicit none
    integer, dimension (:), intent (in out) :: target
    integer, dimension (:), intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_integer_get")
  end subroutine shmem_integer_get_array

  subroutine shmem_integer_get_scalar (target, src, len, pe)
    implicit none
    integer, intent (in out) :: target
    integer, intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_integer_get")
  end subroutine shmem_integer_get_scalar

  subroutine shmem_real_get_array (target, src, len, pe)
    implicit none
    real, dimension (:), intent (in out) :: target
    real, dimension (:), intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_real_get")
  end subroutine shmem_real_get_array

  subroutine shmem_real_get_scalar (target, src, len, pe)
    implicit none
    real, intent (in out) :: target
    real, intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_real_get")
  end subroutine shmem_real_get_scalar

  subroutine shmem_complex_get_array (target, src, len, pe)
    implicit none
    complex, dimension (:), intent (in out) :: target
    complex, dimension (:), intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_complex_get")
  end subroutine shmem_complex_get_array

  subroutine shmem_complex_get_scalar (target, src, len, pe)
    implicit none
    complex, intent (in out) :: target
    complex, intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_complex_get")
  end subroutine shmem_complex_get_scalar

  subroutine shmem_logical_get_array (target, src, len, pe)
    implicit none
    logical, dimension (:), intent (in out) :: target
    logical, dimension (:), intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_logical_get")
  end subroutine shmem_logical_get_array

  subroutine shmem_logical_get_scalar (target, src, len, pe)
    implicit none
    logical, intent (in out) :: target
    logical, intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_logical_get")
  end subroutine shmem_logical_get_scalar

  subroutine shmem_character_get (target, src, len, pe)
    implicit none
    character(*), intent (in out) :: target
    character(*), intent (in) :: src
    integer, intent (in) :: len, pe
    call error ("shmem_character_get")
  end subroutine shmem_character_get

  subroutine shmem_wait (var, cmp)
    implicit none
    integer, intent (in out) :: var
    integer, intent (in) :: cmp
    call error ("shmem_wait")
  end subroutine shmem_wait

  subroutine error (msg)
    implicit none
    character(*), intent (in) :: msg

    print *, "shmem error: "//msg
    call coredump
    stop
  end subroutine error

  subroutine coredump
    real, dimension (:), allocatable :: a
    deallocate (a)
  end subroutine coredump
#endif

end module shmem
