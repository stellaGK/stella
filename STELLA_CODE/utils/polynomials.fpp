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
! Given the base cases P_0(x) = 1 and P_1(x) = 1, three-term recurrence properties mean we can calculate P_{l+1} as:
!
! (l+1) * P_{l+1} = (2*l + 1) * x * P_l(x) - l * P_{l-1}(x)
!
! ================================================================================================================================================================================= !

    pure subroutine get_legendre_array(x, m, P)
        implicit none

        integer :: l

        real, intent(in)    :: x          ! Evaluation point. 
        integer, intent(in) :: m          ! Maximum polynomial order.
        real, intent(out)   :: P(0:m)     ! P(l) = Legendre polynomial of order l.

        ! Define the two base cases first. 
        
        P(0) = 1.0

        if (m >= 1) then
            P(1) = x
        end if

        ! The three-term recurrence property will provide the remaining P(l) that the user requires. 
    
        do l = 1, m - 1
            P(l + 1) = ((2 * l + 1) * x * P(l) - l * P(l - 1)) / (l + 1)
        end do
    end subroutine get_legendre_array


! ================================================================================================================================================================================= !
! ----- Calculates all Legendre polynomial derivatives w.r.t x, upto a maximum integer, m. Results are returned as the array dP with elements dP_l(x) where l = 0, 1, ... , m. ---- !
! ================================================================================================================================================================================= !
!
!
!
! ================================================================================================================================================================================= !

    ! subroutine get_legendre_derivative_array(x, m, dP)
        ! implicit none
    ! end subroutine get_legendre_derivative_array


! ================================================================================================================================================================================= !
! - Calculates all associated Laguerre polynomials, of order q, at value y, upto a maximum integer, n. Results are returned as the array L with elements L_k(y) where k = 0, ... n. ! 
! ================================================================================================================================================================================= !
!
! Given the base cases L^q_0(y) = 1 and L^q_1(y) = 1 + q - y, three-term recurrence properties mean we can calculate L^q_{k+1}(y) as:
!
! (k+1) * L^q_{k+1}(y) = (2*k + 1 + q - y) * L^q_k(y) - (k + q) * L^q_{k-1}(y)
!
! ================================================================================================================================================================================= !

    pure subroutine get_laguerre_array(y, n, q, L)
        implicit none

        integer :: k 

        real, intent(in)    :: y          ! Evaluation point.
        integer, intent(in) :: n          ! Maximum polynomial order
        real, intent(out)   :: L(0:n)     ! L(l) = Associated Laguerre polynomial of order l.
        real, intent(in)    :: q          ! Polynomial order.

        ! Define the two base cases first.

        L(0) = 1.0

        if (n >= 1) then
            L(1) = 1 + q - y
        end if

        ! The three-term recurrence property will provide the remaining L(k) that the user requires.

        do k = 1, n - 1
            L(k + 1) = ((2 * k + 1 + q - y) * L(k) - (k + q) * L(k - 1)) / (k + 1)
        end do
    end subroutine get_laguerre_array


! ================================================================================================================================================================================= !
! ------------------------------------ Calculates all associated Laguerre polynomial derivatives w.r.t y, of order q, upto a maximum integer, n. ---------------------------------- ! 
! ------------------------------------------------ Results are returned as the array dP with elements dL_k(y) where l = 0, 1, ... , n. -------------------------------------------- !
! ================================================================================================================================================================================= !
!
!
!
! ================================================================================================================================================================================= !

    ! subroutine get_laguerre_derivative_array(y, n, q, dL)
        ! implicit none
    ! end subroutine get_laguerre_derivative_array


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------------- End module. ---------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

end module polynomials
