!> Module for dealing with common polynomials such as
!> Legendre, Laguerre etc.
module polynomials
  implicit none

  private

  public :: evaluate_legendre, evaluate_legendre_derivative
  public :: evaluate_laguerre_function, evaluate_laguerre_function_derivative

contains

  !> Computes the value of the Legendre polynomials
  !> at the passed point for all orders from 0 up to
  !> the size of values minus 1.
  !>
  !> Uses a recurrence relation to compute order > 1
  !>
  !> @todo Add checking - consistency of input and validity of result.
  pure subroutine evaluate_legendre(point, values, monic)
    use optionals, only: get_option_with_default
    implicit none
    real, intent(in) :: point
    real, dimension(0:), intent(out) :: values
    logical, intent(in), optional :: monic
    integer :: order, max_order, n

    max_order = size(values)-1

    values(0) = 1.0
    values(1) = point

    do order = 2, max_order
       ! Define n just for short hand and ease of comparison with
       ! typical expressions of the recurrence relation.
       n = order - 1
       values(order) = ((2*n+1)*point*values(n) - n * values(n-1))/(n+1)
    end do

    if (get_option_with_default(monic, .false.)) values = values / values(0)
  end subroutine evaluate_legendre

  !> Computes the derivative of the Legendre polynomials
  !> at the passed point for all orders from 0 up to
  !> the size of values minus 1.
  !>
  !> Uses a recurrence relation to compute order > 1
  !>
  !> @todo Add checking - consistency of input and validity of result.
  pure subroutine evaluate_legendre_derivative(point, values, monic)
    use warning_helpers, only: not_exactly_equal
    implicit none
    real, intent(in) :: point
    real, dimension(0:), intent(out) :: values
    logical, intent(in), optional :: monic
    real, dimension(:), allocatable :: function_values
    integer :: order, max_order, n

    max_order = size(values)-1
    allocate(function_values(0:max_order))
    call evaluate_legendre(point, function_values, monic)

    values(0) = 0.0
    values(1) = 1.0

    ! Note the recurrence relation is not useful when |point|==1
    ! so handle as a special case
    if (not_exactly_equal(abs(point), 1.0)) then
       do order = 2, max_order
          ! Define n just for short hand and ease of comparison with
          ! typical expressions of the recurrence relation.
          n = order - 1
          values(order) = (point*function_values(order) - function_values(order-1)) * order &
               / (point * point - 1)
       end do
    else
       do order = 2, max_order
          values(order) = order*(order+1) / 2.0

          ! Ensure correct sign for odd orders
          if (mod(order, 2) == 0) values(order) = point*values(order)
       end do
    end if
  end subroutine evaluate_legendre_derivative

  !> Computes the value of the Laguerre function
  !> at the passed point for all orders from 0 up to
  !> the size of values minus 1.
  !>
  !> Uses a recurrence relation to compute order > 1
  !>
  !> @todo Add checking - consistency of input and validity of result.
  pure subroutine evaluate_laguerre_function(point, k, values, monic)
    use optionals, only: get_option_with_default
    implicit none
    real, intent(in) :: point, k
    real, dimension(0:), intent(out) :: values
    logical, intent(in), optional :: monic
    integer :: order, max_order, n

    max_order = size(values)-1

    values(0) = 1.0
    values(1) = 1.0 + k - point

    do order = 2, max_order
       ! Define n just for short hand and ease of comparison with
       ! typical expressions of the recurrence relation.
       n = order
       values(order) = ((2*n-1+k-point)*values(n-1) - (n-1+k) * values(n-2))/n
    end do

    if (get_option_with_default(monic, .false.)) values = values / values(0)
  end subroutine evaluate_laguerre_function


  !> Computes the derivative of the Laguerre function
  !> at the passed point for all orders from 0 up to
  !> the size of values minus 1.
  !>
  !> Uses a recurrence relation to compute order > 1
  !>
  !> @todo Add checking - consistency of input and validity of result.
  pure subroutine evaluate_laguerre_function_derivative(point, k, values, monic)
    implicit none
    real, intent(in) :: point, k
    real, dimension(0:), intent(out) :: values
    logical, intent(in), optional :: monic
    real, dimension(:), allocatable :: function_values
    integer :: order, max_order, n

    max_order = size(values)-1
    allocate(function_values(0:max_order))
    call evaluate_laguerre_function(point, k+1, function_values, monic)

    values(0) = 0.0

    do order = 1, max_order
       ! Define n just for short hand and ease of comparison with
       ! typical expressions of the recurrence relation.
       n = order
       values(order) = -function_values(n-1)
    end do
  end subroutine evaluate_laguerre_function_derivative

end module polynomials
