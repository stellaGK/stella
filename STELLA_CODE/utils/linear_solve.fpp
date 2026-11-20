!###############################################################################
!                      LU decomposition of dense matrices
!###############################################################################
! The LU decomposition, back substitution and inverse of dense matrices are
! computed by LAPACK (?getrf, ?getrs and ?getri). When stella is compiled with
! PLASMA (see plasma_interface), PLASMA's plasma_?getrf, plasma_?getrs and
! plasma_?getri are used instead, except inside OpenMP parallel regions, where
! PLASMA can't be called (see <plasma_usable>).
!
! The factors are stored as by LAPACK: L (unit diagonal) below and U on and
! above the diagonal, and idx(j) is the row that was swapped with row j at
! step j. They are the same for both libraries, so the factors of one can be
! solved by the other.
!
! The LAPACK routines of the precision of stella are called: d/z with double
! precision (the default), s/c otherwise. PLASMA needs double precision.
!
! A singular matrix is still factorised: its zero pivots, i.e., zeros on the
! diagonal of U, are replaced by a small number, so that the solves do not give
! NaNs or Infs. This is the factorisation of a nearby matrix, since LAPACK
! leaves the column of L below a zero pivot zero.
!###############################################################################
module linear_solve

   implicit none

   public :: lu_decomposition
   public :: lu_back_substitution
   public :: lu_inverse

   private

   interface lu_decomposition
      module procedure lu_decomposition_real
      module procedure lu_decomposition_complex
   end interface

   interface lu_back_substitution
      module procedure lu_back_substitution_real
      module procedure lu_back_substitution_real_complex
      module procedure lu_back_substitution_complex
      module procedure lu_back_substitution_matrix_real
      module procedure lu_back_substitution_matrix_complex
   end interface

   interface lu_inverse
      module procedure lu_inverse_real
      module procedure lu_inverse_complex
   end interface

   ! The value of the zero pivots of a singular matrix
   real, parameter :: zero_pivot = 1.0e-20

contains

   !****************************************************************************
   !                              LU decomposition
   !****************************************************************************
   ! On exit, <lu> contains the LU factors of the square matrix <lu>, <idx> the
   ! row swaps and <d> the sign of the permutation (+1 or -1), so that the
   ! determinant is d * product(diagonal of lu).
   !****************************************************************************
   subroutine lu_decomposition_real(lu, idx, d)

#ifdef PLASMA
      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_usable, plasma_dgetrf
#endif

      implicit none

      real, dimension(:, :), intent(in out) :: lu
      integer, dimension(:), intent(out) :: idx
      real, intent(out) :: d

      integer :: n, info
#ifdef PLASMA
      real, dimension(:, :), allocatable :: a
      integer(c_int), dimension(size(lu, 1)) :: ipiv
#endif

      !-------------------------------------------------------------------------

      n = size(lu, 1)

#ifdef PLASMA
      if (plasma_usable()) then
         ! plasma_dgetrf does not return the factors of a singular matrix, which
         ! are then computed by LAPACK
         a = lu
         info = plasma_dgetrf(int(n, c_int), int(n, c_int), lu, int(n, c_int), ipiv)
         if (info == 0) then
            idx(1:n) = ipiv
            d = permutation_sign(idx, n)
            return
         end if
         lu = a
      end if
#endif

      if (kind(lu) == kind(1.0d0)) then
         call dgetrf(n, n, lu, size(lu, 1), idx, info)
      else
         call sgetrf(n, n, lu, size(lu, 1), idx, info)
      end if
      if (info > 0) call replace_zero_pivots_real(lu, n)
      d = permutation_sign(idx, n)

   end subroutine lu_decomposition_real

   subroutine lu_decomposition_complex(lu, idx, d)

#ifdef PLASMA
      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_usable, plasma_zgetrf
#endif

      implicit none

      complex, dimension(:, :), intent(in out) :: lu
      integer, dimension(:), intent(out) :: idx
      real, intent(out) :: d

      integer :: n, info
#ifdef PLASMA
      complex, dimension(:, :), allocatable :: a
      integer(c_int), dimension(size(lu, 1)) :: ipiv
#endif

      !-------------------------------------------------------------------------

      n = size(lu, 1)

#ifdef PLASMA
      if (plasma_usable()) then
         ! plasma_zgetrf does not return the factors of a singular matrix, which
         ! are then computed by LAPACK
         a = lu
         info = plasma_zgetrf(int(n, c_int), int(n, c_int), lu, int(n, c_int), ipiv)
         if (info == 0) then
            idx(1:n) = ipiv
            d = permutation_sign(idx, n)
            return
         end if
         lu = a
      end if
#endif

      if (kind(lu) == kind((1.0d0, 0.0d0))) then
         call zgetrf(n, n, lu, size(lu, 1), idx, info)
      else
         call cgetrf(n, n, lu, size(lu, 1), idx, info)
      end if
      if (info > 0) call replace_zero_pivots_complex(lu, n)
      d = permutation_sign(idx, n)

   end subroutine lu_decomposition_complex

   !****************************************************************************
   ! Replace the zeros on the diagonal of U of a singular matrix, see the top of
   ! this module
   !****************************************************************************
   subroutine replace_zero_pivots_real(lu, n)

      implicit none

      real, dimension(:, :), intent(in out) :: lu
      integer, intent(in) :: n

      integer :: j

      write (*, *) 'singular matrix in lu_decomposition'
      do j = 1, n
         if (lu(j, j) == 0.0) lu(j, j) = zero_pivot
      end do

   end subroutine replace_zero_pivots_real

   subroutine replace_zero_pivots_complex(lu, n)

      implicit none

      complex, dimension(:, :), intent(in out) :: lu
      integer, intent(in) :: n

      integer :: j

      write (*, *) 'singular matrix in lu_decomposition'
      do j = 1, n
         if (lu(j, j) == 0.0) lu(j, j) = zero_pivot
      end do

   end subroutine replace_zero_pivots_complex

   ! The sign of the permutation of the LU decomposition, i.e., -1 for every row swap
   pure real function permutation_sign(idx, n)

      implicit none

      integer, dimension(:), intent(in) :: idx
      integer, intent(in) :: n

      integer :: j

      permutation_sign = 1.0
      do j = 1, n
         if (idx(j) /= j) permutation_sign = -permutation_sign
      end do

   end function permutation_sign

   !****************************************************************************
   !                             Back substitution
   !****************************************************************************
   ! Solves A x = b, given the LU factors <lu> and <idx> of A from
   ! <lu_decomposition>. On exit, <b> contains x. <b> can be a vector or a
   ! matrix with a right-hand side in every column.
   !****************************************************************************
   subroutine lu_back_substitution_real(lu, idx, b)

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      real, dimension(:), intent(in out) :: b

      call getrs_real(lu, idx, size(lu, 1), 1, b, size(b))

   end subroutine lu_back_substitution_real

   subroutine lu_back_substitution_real_complex(lu, idx, b)

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:), intent(in out) :: b

      real, dimension(size(lu, 1), 2) :: br
      integer :: n

      ! Solve for the real and imaginary parts as two right-hand sides
      n = size(lu, 1)
      br(:, 1) = real(b(1:n))
      br(:, 2) = aimag(b(1:n))
      call getrs_real(lu, idx, n, 2, br, n)
      b(1:n) = cmplx(br(:, 1), br(:, 2))

   end subroutine lu_back_substitution_real_complex

   subroutine lu_back_substitution_complex(lu, idx, b)

      implicit none

      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:), intent(in out) :: b

      call getrs_complex(lu, idx, size(lu, 1), 1, b, size(b))

   end subroutine lu_back_substitution_complex

   subroutine lu_back_substitution_matrix_real(lu, idx, b)

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      real, dimension(:, :), intent(in out) :: b

      call getrs_real(lu, idx, size(lu, 1), size(b, 2), b, size(b, 1))

   end subroutine lu_back_substitution_matrix_real

   subroutine lu_back_substitution_matrix_complex(lu, idx, b)

      implicit none

      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:, :), intent(in out) :: b

      call getrs_complex(lu, idx, size(lu, 1), size(b, 2), b, size(b, 1))

   end subroutine lu_back_substitution_matrix_complex

   !****************************************************************************
   ! Solves A x = b for the <nrhs> right-hand sides in the first <n> rows of the
   ! columns of b (leading dimension <ldb>), given the LU factors of A.
   !****************************************************************************
   subroutine getrs_real(lu, idx, n, nrhs, b, ldb)

#ifdef PLASMA
      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_usable, plasma_dgetrs, PlasmaNoTrans
#endif

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      integer, intent(in) :: n, nrhs, ldb
      real, dimension(ldb, nrhs), intent(in out) :: b

      integer :: info
#ifdef PLASMA
      real, dimension(:, :), allocatable :: a
#endif

      !-------------------------------------------------------------------------

#ifdef PLASMA
      if (plasma_usable()) then
         ! plasma_dgetrs does not change the factors, but has no const arguments
         a = lu
         info = plasma_dgetrs(PlasmaNoTrans, int(n, c_int), int(nrhs, c_int), a, int(size(lu, 1), c_int), &
                              int(idx(1:n), c_int), b, int(ldb, c_int))
         return
      end if
#endif

      if (kind(lu) == kind(1.0d0)) then
         call dgetrs('N', n, nrhs, lu, size(lu, 1), idx, b, ldb, info)
      else
         call sgetrs('N', n, nrhs, lu, size(lu, 1), idx, b, ldb, info)
      end if

   end subroutine getrs_real

   subroutine getrs_complex(lu, idx, n, nrhs, b, ldb)

#ifdef PLASMA
      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_usable, plasma_zgetrs, PlasmaNoTrans
#endif

      implicit none

      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      integer, intent(in) :: n, nrhs, ldb
      complex, dimension(ldb, nrhs), intent(in out) :: b

      integer :: info
#ifdef PLASMA
      complex, dimension(:, :), allocatable :: a
#endif

      !-------------------------------------------------------------------------

#ifdef PLASMA
      if (plasma_usable()) then
         ! plasma_zgetrs does not change the factors, but has no const arguments
         a = lu
         info = plasma_zgetrs(PlasmaNoTrans, int(n, c_int), int(nrhs, c_int), a, int(size(lu, 1), c_int), &
                              int(idx(1:n), c_int), b, int(ldb, c_int))
         return
      end if
#endif

      if (kind(lu) == kind((1.0d0, 0.0d0))) then
         call zgetrs('N', n, nrhs, lu, size(lu, 1), idx, b, ldb, info)
      else
         call cgetrs('N', n, nrhs, lu, size(lu, 1), idx, b, ldb, info)
      end if

   end subroutine getrs_complex

   !****************************************************************************
   !                               Matrix inverse
   !****************************************************************************
   ! One shouldn't need to compute the inverse of a matrix if solving the linear
   ! equation A.x = y; a simple back substitution should suffice. Only compute
   ! the inverse of A when absolutely necessary!
   !
   ! Computes the inverse of A from its LU factors <lu> and <idx> (from
   ! <lu_decomposition>).
   !****************************************************************************
   subroutine lu_inverse_real(lu, idx, inverse)

#ifdef PLASMA
      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_usable, plasma_dgetri
#endif

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      real, dimension(:, :), intent(out) :: inverse

      real, dimension(size(lu, 1)) :: work
      integer :: n, info

      !-------------------------------------------------------------------------

      n = size(lu, 1)
      inverse(1:n, 1:n) = lu(1:n, 1:n)

#ifdef PLASMA
      if (plasma_usable()) then
         info = plasma_dgetri(int(n, c_int), inverse, int(size(inverse, 1), c_int), int(idx(1:n), c_int))
         return
      end if
#endif

      if (kind(lu) == kind(1.0d0)) then
         call dgetri(n, inverse, size(inverse, 1), idx, work, n, info)
      else
         call sgetri(n, inverse, size(inverse, 1), idx, work, n, info)
      end if

   end subroutine lu_inverse_real

   subroutine lu_inverse_complex(lu, idx, inverse)

#ifdef PLASMA
      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_usable, plasma_zgetri
#endif

      implicit none

      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:, :), intent(out) :: inverse

      complex, dimension(size(lu, 1)) :: work
      integer :: n, info

      !-------------------------------------------------------------------------

      n = size(lu, 1)
      inverse(1:n, 1:n) = lu(1:n, 1:n)

#ifdef PLASMA
      if (plasma_usable()) then
         info = plasma_zgetri(int(n, c_int), inverse, int(size(inverse, 1), c_int), int(idx(1:n), c_int))
         return
      end if
#endif

      if (kind(lu) == kind((1.0d0, 0.0d0))) then
         call zgetri(n, inverse, size(inverse, 1), idx, work, n, info)
      else
         call cgetri(n, inverse, size(inverse, 1), idx, work, n, info)
      end if

   end subroutine lu_inverse_complex

end module linear_solve
