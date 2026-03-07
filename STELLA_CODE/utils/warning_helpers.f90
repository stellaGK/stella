!> Small functions for making semantically-correct code that doesn't
!> raise compiler warnings (or at least, fewer). The main downside of
!> this is requiring a function call, which can destroy vectorisation
!> of loops.
module warning_helpers
  use, intrinsic :: iso_fortran_env, only: real32, real64, real128
  implicit none

  private

  public :: exactly_equal, not_exactly_equal
  public :: is_zero, is_not_zero
  public :: almost_equal
  public :: complex_sq_mod

  !> Are two numbers exactly equal to each other
  interface exactly_equal
    module procedure :: exactly_equal_r32, exactly_equal_r64, exactly_equal_r128
    module procedure :: exactly_equal_c32, exactly_equal_c64, exactly_equal_c128
  end interface exactly_equal

  !> Are two numbers not exactly equal to each other
  interface not_exactly_equal
    module procedure :: not_exactly_equal_r32, not_exactly_equal_r64, not_exactly_equal_r128
    module procedure :: not_exactly_equal_c32, not_exactly_equal_c64, not_exactly_equal_c128
  end interface not_exactly_equal

  !> Is the argument exactly equal to zero
  interface is_zero
    module procedure :: is_zero_r32, is_zero_r64, is_zero_r128
    module procedure :: is_zero_c32, is_zero_c64, is_zero_c128
  end interface is_zero

  !> Is the argument not exactly zero
  interface is_not_zero
    module procedure :: is_not_zero_r32, is_not_zero_r64, is_not_zero_r128
    module procedure :: is_not_zero_c32, is_not_zero_c64, is_not_zero_c128
  end interface is_not_zero

  !> Are two numbers almost equal within some tolerance
  interface almost_equal
    module procedure :: almost_equal_r32, almost_equal_r64, almost_equal_r128
  end interface almost_equal

  !> Is the argument almost zero
  interface almost_zero
    module procedure :: almost_zero_r32, almost_zero_r64, almost_zero_r128
  end interface almost_zero

  !> Returns \(z \bar{z} = |z|^2\), the square modulus of a complex number
  interface complex_sq_mod
    module procedure :: complex_sq_mod_r32, complex_sq_mod_r64, complex_sq_mod_r128
  end interface complex_sq_mod

contains
  logical elemental function exactly_equal_r32(a, b)
    real(real32), intent(in) :: a, b
    exactly_equal_r32 = abs(a - b) <= 0.0
  end function exactly_equal_r32

  logical elemental function exactly_equal_r64(a, b)
    real(real64), intent(in) :: a, b
    exactly_equal_r64 = abs(a - b) <= 0.0
  end function exactly_equal_r64

  logical elemental function exactly_equal_r128(a, b)
    real(real128), intent(in) :: a, b
    exactly_equal_r128 = abs(a - b) <= 0.0
  end function exactly_equal_r128

  logical elemental function exactly_equal_c32(a, b)
    complex(real32), intent(in) :: a, b
    exactly_equal_c32 = exactly_equal(real(a), real(b)) .and. exactly_equal(aimag(a), aimag(b))
  end function exactly_equal_c32

  logical elemental function exactly_equal_c64(a, b)
    complex(real64), intent(in) :: a, b
    exactly_equal_c64 = exactly_equal(real(a), real(b)) .and. exactly_equal(aimag(a), aimag(b))
  end function exactly_equal_c64

  logical elemental function exactly_equal_c128(a, b)
    complex(real128), intent(in) :: a, b
    exactly_equal_c128 = exactly_equal(real(a), real(b)) .and. exactly_equal(aimag(a), aimag(b))
  end function exactly_equal_c128

  logical elemental function not_exactly_equal_r32(a, b)
    real(real32), intent(in) :: a, b
    not_exactly_equal_r32 = .not. exactly_equal(a, b)
  end function not_exactly_equal_r32

  logical elemental function not_exactly_equal_r64(a, b)
    real(real64), intent(in) :: a, b
    not_exactly_equal_r64 = .not. exactly_equal(a, b)
  end function not_exactly_equal_r64

  logical elemental function not_exactly_equal_r128(a, b)
    real(real128), intent(in) :: a, b
    not_exactly_equal_r128 = .not. exactly_equal(a, b)
  end function not_exactly_equal_r128

  logical elemental function not_exactly_equal_c32(a, b)
    complex(real32), intent(in) :: a, b
    not_exactly_equal_c32 = .not. exactly_equal(a, b)
  end function not_exactly_equal_c32

  logical elemental function not_exactly_equal_c64(a, b)
    complex(real64), intent(in) :: a, b
    not_exactly_equal_c64 = .not. exactly_equal(a, b)
  end function not_exactly_equal_c64

  logical elemental function not_exactly_equal_c128(a, b)
    complex(real128), intent(in) :: a, b
    not_exactly_equal_c128 = .not. exactly_equal(a, b)
  end function not_exactly_equal_c128

  logical elemental function is_zero_r32(scalar)
    real(real32), intent(in) :: scalar
    is_zero_r32 = exactly_equal(scalar, 0.0_real32)
  end function is_zero_r32

  logical elemental function is_zero_r64(scalar)
    real(real64), intent(in) :: scalar
    is_zero_r64 = exactly_equal(scalar, 0.0_real64)
  end function is_zero_r64

  logical elemental function is_zero_r128(scalar)
    real(real128), intent(in) :: scalar
    is_zero_r128 = exactly_equal(scalar, 0.0_real128)
  end function is_zero_r128

  logical elemental function is_zero_c32(scalar)
    complex(real32), intent(in) :: scalar
    is_zero_c32 = exactly_equal(scalar, cmplx(0.0, 0.0, kind = real32))
  end function is_zero_c32

  logical elemental function is_zero_c64(scalar)
    complex(real64), intent(in) :: scalar
    is_zero_c64 = exactly_equal(scalar, cmplx(0.0, 0.0, kind = real64))
  end function is_zero_c64

  logical elemental function is_zero_c128(scalar)
    complex(real128), intent(in) :: scalar
    is_zero_c128 = exactly_equal(scalar, cmplx(0.0, 0.0, kind = real128))
  end function is_zero_c128

  logical elemental function is_not_zero_r32(scalar)
    real(real32), intent(in) :: scalar
    is_not_zero_r32 = .not. is_zero(scalar)
  end function is_not_zero_r32

  logical elemental function is_not_zero_r64(scalar)
    real(real64), intent(in) :: scalar
    is_not_zero_r64 = .not. is_zero(scalar)
  end function is_not_zero_r64

  logical elemental function is_not_zero_r128(scalar)
    real(real128), intent(in) :: scalar
    is_not_zero_r128 = .not. is_zero(scalar)
  end function is_not_zero_r128

  logical elemental function is_not_zero_c32(scalar)
    complex(real32), intent(in) :: scalar
    is_not_zero_c32 = .not. is_zero(scalar)
  end function is_not_zero_c32

  logical elemental function is_not_zero_c64(scalar)
    complex(real64), intent(in) :: scalar
    is_not_zero_c64 = .not. is_zero(scalar)
  end function is_not_zero_c64

  logical elemental function is_not_zero_c128(scalar)
    complex(real128), intent(in) :: scalar
    is_not_zero_c128 = .not. is_zero(scalar)
  end function is_not_zero_c128

  !> Are two numbers almost equal
  !>
  !> Uses the same predicate as numpy's `isclose`:
  !>
  !>     abs(a - b) <= (atol + rtol * abs(b))
  logical elemental function almost_equal_r32(a, b, rtol, atol)
    use optionals, only: get_option_with_default
    real(real32), intent(in) :: a, b
    real(real32), intent(in), optional :: rtol
    real(real32), intent(in), optional :: atol

    real(real32) :: rtol_val, atol_val

    rtol_val = get_option_with_default(rtol, 1e-5_real32)
    atol_val = get_option_with_default(atol, 1e-8_real32)

    almost_equal_r32 = (abs(a - b) <= (atol_val + (rtol_val * abs(b))))
  end function almost_equal_r32

  !> Are two numbers almost equal
  !>
  !> Uses the same predicate as numpy's `isclose`:
  !>
  !>     abs(a - b) <= (atol + rtol * abs(b))
  logical elemental function almost_equal_r64(a, b, rtol, atol)
    use optionals, only: get_option_with_default
    real(real64), intent(in) :: a, b
    real(real64), intent(in), optional :: rtol
    real(real64), intent(in), optional :: atol

    real(real64) :: rtol_val, atol_val

    rtol_val = get_option_with_default(rtol, 1e-5_real64)
    atol_val = get_option_with_default(atol, 1e-8_real64)

    almost_equal_r64 = (abs(a - b) <= (atol_val + (rtol_val * abs(b))))
  end function almost_equal_r64

  !> Are two numbers almost equal
  !>
  !> Uses the same predicate as numpy's `isclose`:
  !>
  !>     abs(a - b) <= (atol + rtol * abs(b))
  logical elemental function almost_equal_r128(a, b, rtol, atol)
    use optionals, only: get_option_with_default
    real(real128), intent(in) :: a, b
    real(real128), intent(in), optional :: rtol
    real(real128), intent(in), optional :: atol

    real(real128) :: rtol_val, atol_val

    rtol_val = get_option_with_default(rtol, 1e-5_real128)
    atol_val = get_option_with_default(atol, 1e-8_real128)

    almost_equal_r128 = (abs(a - b) <= (atol_val + (rtol_val * abs(b))))
  end function almost_equal_r128

  !> Is a scalar almost zero
  logical elemental function almost_zero_r32(scalar, rtol, atol)
    use optionals, only: get_option_with_default
    real(real32), intent(in) :: scalar
    real(real32), intent(in), optional :: rtol
    real(real32), intent(in), optional :: atol

    almost_zero_r32 = almost_equal(scalar, 0.0_real32, rtol, &
         atol=get_option_with_default(atol, 0.0_real32))
  end function almost_zero_r32

  logical elemental function almost_zero_r64(scalar, rtol, atol)
    use optionals, only: get_option_with_default
    real(real64), intent(in) :: scalar
    real(real64), intent(in), optional :: rtol
    real(real64), intent(in), optional :: atol

    almost_zero_r64 = almost_equal(scalar, 0.0_real64, rtol, &
         atol=get_option_with_default(atol, 0.0_real64))
  end function almost_zero_r64

  logical elemental function almost_zero_r128(scalar, rtol, atol)
    use optionals, only: get_option_with_default
    real(real128), intent(in) :: scalar
    real(real128), intent(in), optional :: rtol
    real(real128), intent(in), optional :: atol

    almost_zero_r128 = almost_equal(scalar, 0.0_real128, rtol, &
         atol=get_option_with_default(atol, 0.0_real128))
  end function almost_zero_r128

  real(real32) elemental function complex_sq_mod_r32(scalar)
    complex(real32), intent(in) :: scalar
    complex_sq_mod_r32 = real(scalar * conjg(scalar))
  end function complex_sq_mod_r32

  real(real64) elemental function complex_sq_mod_r64(scalar)
    complex(real64), intent(in) :: scalar
    complex_sq_mod_r64 = real(scalar * conjg(scalar))
  end function complex_sq_mod_r64

  real(real128) elemental function complex_sq_mod_r128(scalar)
    complex(real128), intent(in) :: scalar
    complex_sq_mod_r128 = real(scalar * conjg(scalar))
  end function complex_sq_mod_r128

end module warning_helpers
