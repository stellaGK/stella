! ================================================================================================================================================================================= !
! ---------------------------------------- Routines for calculating the values of Legendre and Laguerre polynomials and their derivatives. ---------------------------------------- !
! ================================================================================================================================================================================= !
!
! These are needed for calculations in higher order GK simulations, primarily when constructing NEO's H_1 data on the stella velocity grids.
!
! ================================================================================================================================================================================= !

module polynomials
    implicit none

    public :: get_legendre_array, get_laguerre_array
    ! public :: get_legendre_derivative_array, get_laguerre_derivative_array

    private

contains

! ================================================================================================================================================================================= !
! ------- Calculates all Legendre polynomials at value x, upto a maximum integer, m. Results are returned as the array P with elements P_l(x) where l = 0, 1, ... , m. ------------ !
! ================================================================================================================================================================================= !
!
! Given the base cases P_0(x) = 1 and P_1(x) = x, the three-term recurrence property means we can calculate P_{l+1} as:
!
! (l+1) * P_{l+1} = (2*l + 1) * x * P_l(x) - l * P_{l-1}(x)
!
! If calc_deriv = .true., we compute the derivatives too. The derivatives follow the recurrence relation: 
!
! dP_{l+1}/dx = (l+1) * P_l(x) + x * dP_l(x)/dx 
!
! ================================================================================================================================================================================= !

    pure subroutine get_legendre_array(x, m, P, derivative, dP)
        implicit none

        real, intent(in)              :: x          ! Evaluation point. 
        integer, intent(in)           :: m          ! Maximum polynomial order.
        logical, intent(in), optional :: derivative ! Flag to calculate derivatives.
        real, intent(out), optional   :: dP(0:m)    ! P'(l) = Derivative of order l.
        real, intent(out)             :: P(0:m)     ! P(l) = Legendre polynomial of order l.

        logical :: calc_deriv
        integer :: l

        ! Determine if we should calculate derivatives.
        calc_deriv = .false.
        if (present(derivative)) calc_deriv = derivative

        if (calc_deriv .and. .not. present(dP)) then
            error stop "derivative=.true. but dP not provided"
        end if

        ! Define the two base cases first. 

        P(0) = 1.0
        if (calc_deriv .and. present(dP)) dP(0) = 0.0

        if (m >= 1) then
            P(1) = x
            if (calc_deriv .and. present(dP)) dP(1) = 1.0
        end if

        ! The three-term recurrence property will provide the remaining P(l) that the user requires. 
    
        do l = 1, m - 1
            P(l + 1) = ((2 * l + 1) * x * P(l) - l * P(l - 1)) / (l + 1)

            if (calc_deriv .and. present(dP)) then
                dP(l + 1) = (l + 1) * P(l) + x * dP(l)
            end if
        end do
    end subroutine get_legendre_array


! ================================================================================================================================================================================= !
! - Calculates all associated Laguerre polynomials, of order q, at value y, upto a maximum integer, n. Results are returned as the array L with elements L_k(y) where k = 0, ... n. ! 
! ================================================================================================================================================================================= !
!
! Given the base cases L^q_0(y) = 1 and L^q_1(y) = 1 + q - y, three-term recurrence properties mean we can calculate L^q_{k+1}(y) as:
!
! (k+1) * L^q_{k+1}(y) = (2*k + 1 + q - y) * L^q_k(y) - (k + q) * L^q_{k-1}(y)
!
! ! If calc_deriv = .true., we compute the derivatives too. The derivatives follow the recurrence relation:
!
! dL_{k+1}^q(y)/dy = dL_{k}^q(y)/dy - L_{k}^q(y)
!
! ================================================================================================================================================================================= !

    pure subroutine get_laguerre_array(y, n, q, L, derivative, dL)
        implicit none

        real, intent(in)              :: y          ! Evaluation point.
        integer, intent(in)           :: n          ! Maximum polynomial order
        real, intent(out)             :: L(0:n)     ! L(l) = Associated Laguerre polynomial of order l.
        real, intent(in)              :: q          ! Polynomial order.
        logical, intent(in), optional :: derivative ! Flag to calculate derivatives.
        real, intent(out), optional   :: dL(0:n)    ! L'(l) = Derivative of order l.

        logical :: calc_deriv
        integer :: k

        ! Determine if we should calculate derivatives.
        calc_deriv = .false.
        if (present(derivative)) calc_deriv = derivative

        if (calc_deriv .and. .not. present(dL)) then
            error stop "derivative=.true. but dL not provided"
        end if

        ! Define the two base cases first.

        L(0) = 1.0
        if (calc_deriv .and. present(dL)) dL(0) = 0.0

        if (n >= 1) then
            L(1) = 1 + q - y
            if (calc_deriv .and. present(dL)) dL(1) = -1.0
        end if

        ! The three-term recurrence property will provide the remaining L(k) that the user requires.

        do k = 1, n - 1
            L(k + 1) = ((2 * k + 1 + q - y) * L(k) - (k + q) * L(k - 1)) / (k + 1)

            if (calc_deriv .and. present(dL)) then
                dL(k + 1) = dL(k) - L(k) 
            end if
        end do
    end subroutine get_laguerre_array


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------------- End module. ---------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

end module polynomials
