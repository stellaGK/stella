! ================================================================================================================================================================================= !
! ---------------------------------------- Routines for calculating the values of Legendre and Laguerre polynomials and their derivatives. ---------------------------------------- !
! ================================================================================================================================================================================= !
!
! These are needed for calculations in higher order GK simulations, primarily when constructing NEO's H_1 data on the stella velocity grids.
!
! ================================================================================================================================================================================= !

module polynomials
    implicit none

    public :: get_legendre, get_laguerre
    public :: get_legendre_deriv, get_laguerre_deriv

    private

contains

! ================================================================================================================================================================================= !
! ----------------------------- Calculates all Legendre polynomials at value xi, upto a maximum integer, m - 1. Results are returned as the array P. ------------------------------ !
! ================================================================================================================================================================================= !
!
! 
!
! ================================================================================================================================================================================= !

    pure subroutine get_legendre(xi_in, P)
        implicit none

        real, intent(in)     :: xi_in     ! Evaluation point. 
        real, intent(out)    :: P(0:)     ! P(l) = Legendre polynomial of order l.

        integer :: n, order, max_order
 
        ! Determine the maximum order up to which we calculate the polynomials.
        max_order = size(P) - 1

        ! Define the two base cases first.
        P(0) = 1.0
        P(1) = xi_in 

        ! The three-term recurrence relation  will provide the remaining P(l). 
        do order = 2, max_order
            n = order - 1
            P(order) = ( ( 2*n + 1 ) * xi_in * P(n) - n * P(n - 1) ) / (n + 1)
        end do      
    end subroutine get_legendre


! ================================================================================================================================================================================= !
! ------------------- Calculates all associated Laguerre polynomials, of order q, at value E_in, upto a max integer, size(L) -1. Results are returned as the array L. ------------- ! 
! ================================================================================================================================================================================= !
!
! 
!
! ================================================================================================================================================================================= !

    pure subroutine get_laguerre(E_in, q, L)
        implicit none

        real, intent(in)     :: E_in        ! Evaluation point.
        real, intent(in)     :: q           ! Polynomial order.
        real, intent(out)    :: L(0:)       ! L(l) = Associated Laguerre polynomial of order l

        integer :: n, order, max_order 

        ! Define the two base cases first.

        max_order = size(L) - 1

        L(0) = 1.0
        L(1) = 1.0 + q - E_in

        ! The three-term recurrence property will provide the remaining L(k).
        do order = 2, max_order
            n = order
            L(order) = ( ( 2*n - 1 + q - E_in)*L(n - 1) - (n - 1 + q) * L(n - 2) ) / n
        end do
    end subroutine get_laguerre


! ================================================================================================================================================================================= !
! ---------------------- Calculates all Legendre derivative polynomials at value xi, upto a maximum integer, m - 1. Results are returned as the array dP. ------------------------- !
! ================================================================================================================================================================================= !
!
! 
!
! ================================================================================================================================================================================= !


    pure subroutine get_legendre_deriv(xi_in, dP)
        use warning_helpers, only: not_exactly_equal

        implicit none

        real, intent(in)     :: xi_in     ! Evaluation point. 
        real, intent(out)    :: dP(0:)     ! P(l) = Legendre polynomial of order l.

        real, dimension(:), allocatable :: P
        integer :: n, order, max_order

        ! Determine the maximum order up to which we calculate the polynomials.
        max_order = size(dP) - 1
    
        ! Get the legendre polynomials first.
        allocate(P(0:max_order))
        call get_legendre(xi_in, P)

        ! Define the two base cases first.
        dP(0) = 0.0
        dP(1) = 1.0

        ! The three-term recurrence relation  will provide the remaining dP(l). We must handle the x_in = 1 point with care if it exists.
        if (not_exactly_equal(abs(xi_in), 1.0)) then
            do order = 2, max_order
                n = order - 1
                dP(order) = (xi_in * P(order) - P(n)) * order / (xi_in * xi_in - 1)
            end do
        else
            do order = 2, max_order
                dP(order) = order * ( order + 1 ) / 2.0

                ! Ensure correct sign for odd orders.
                if (mod(order, 2) == 0) dP(order) = xi_in * dP(order)
            end do
        end if

        deallocate(P)
    end subroutine get_legendre_deriv


! ================================================================================================================================================================================= !
! ----------- Calculates all associated Laguerre derivative polynomials, of order q, at value E_in, upto a max integer, size(dL) -1. Results are returned as the array dL. -------- !
! ================================================================================================================================================================================= !
!
! 
!
! ================================================================================================================================================================================= !

    pure subroutine get_laguerre_deriv(E_in, q, dL)
        implicit none

        real, intent(in)     :: E_in        ! Evaluation point.
        real, intent(in)     :: q           ! Polynomial order.
        real, intent(out)    :: dL(0:)      ! dL(l) = Associated Laguerre derivative polynomials.

        integer :: n, order, max_order
        real, dimension(:), allocatable :: L

        ! Determine the maximum order up to which we calculate the polynomials.
        max_order = size(dL) - 1

        ! Get the legendre polynomials first.
        allocate(L(0:max_order))
        call get_laguerre(E_in, q + 1, L)

        ! Define the base case. 
        dL(0) = 0.0
 
        do order = 1, max_order
            n = order - 1
            dL(order) = -L(n)
        end do
  
        deallocate(L)
    end subroutine get_laguerre_deriv


! ================================================================================================================================================================================= !
! ---------------------------------------------------------------------------------- End module. ---------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

end module polynomials
