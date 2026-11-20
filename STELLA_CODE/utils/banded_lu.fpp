!###############################################################################
!                    LU factorisation of complex band matrices
!###############################################################################
! Factorises and solves the square complex band systems of stella, with n rows,
! kl sub-diagonals and ku super-diagonals, stored as in LAPACK: A(i,j) is
! ab(band_row(i, j, kl, ku), j), and the array has band_ldab(kl, ku) rows.
!
! By default, LAPACK's zgbtrf and zgbtrs are used (cgbtrf and cgbtrs in single
! precision). When stella is compiled with
! PLASMA (see plasma_interface), PLASMA's plasma_zgbtrf factorises the matrix
! instead, unless we are inside an OpenMP parallel region, where plasma_zgbtrf
! can't be called (see <plasma_usable>). Then:
!  - PLASMA works on tiles of <band_nb> x <band_nb> elements, and needs more
!    rows in <ab> than LAPACK, hence <band_ldab>.
!  - Within every panel of <band_nb> columns, plasma_zgbtrf also applies the row
!    swaps to the earlier columns of L, as in a dense LU. These swaps are undone
!    afterwards (see <plasma_band_factors_to_lapack>), so that the factors are
!    stored as those of zgbtrf, whichever library computed them. Hence the
!    factors are always solved with zgbtrs, also inside OpenMP parallel regions.
!  - plasma_zgbtrf parallelises every factorisation itself, so the callers that
!    factorise in a loop make their loop serial when <band_lu_factorise_threaded>
!    is true. Otherwise every factorisation in the loop calls zgbtrf.
!###############################################################################
module banded_lu

   implicit none

   private

   public :: band_ldab, band_row
   public :: band_lu_factorise, band_lu_solve
   public :: band_lu_factorise_threaded

#ifdef PLASMA
   ! The tile size of PLASMA for the band routines. The band systems of stella
   ! are small, so the default tile size of PLASMA (256) would waste memory and time
   integer, parameter :: band_nb = 16
#endif

contains

   !****************************************************************************
   !                 Leading dimension of the band storage
   !****************************************************************************
   ! The number of rows of the array <ab> that stores a band matrix with <kl>
   ! sub-diagonals and <ku> super-diagonals and its LU factors.
   !****************************************************************************
   pure integer function band_ldab(kl, ku)

      implicit none

      integer, intent(in) :: kl, ku

      ! The rows needed by zgbtrf
      band_ldab = 2 * kl + ku + 1
#ifdef PLASMA
      ! (tiles in the upper band + tiles in the lower band + 1) * band_nb, see
      ! plasma/compute/zgbtrf.c. With PLASMA this is chosen at runtime, so the
      ! array must be large enough for both
      band_ldab = max(band_ldab, ((kl + ku + band_nb - 1) / band_nb + (kl + band_nb - 1) / band_nb + 1) * band_nb)
#endif

   end function band_ldab

   !****************************************************************************
   ! The row of <ab> that stores the element A(i,j) of the band matrix, which
   ! is stored in column j of <ab>. The first <kl> rows are needed for the fill-in
   ! of the LU factorisation.
   !****************************************************************************
   pure integer function band_row(i, j, kl, ku)

      implicit none

      integer, intent(in) :: i, j, kl, ku

      band_row = kl + ku + 1 + i - j

   end function band_row

   !****************************************************************************
   ! Whether <band_lu_factorise>, called here, parallelises the factorisation
   ! itself (with PLASMA). A loop over factorisations should then be serial.
   !****************************************************************************
   logical function band_lu_factorise_threaded()

      use plasma_interface, only: plasma_usable

      implicit none

      band_lu_factorise_threaded = plasma_usable()

   end function band_lu_factorise_threaded

   !****************************************************************************
   !                       LU factorisation of a band matrix
   !****************************************************************************
   ! On entry, <ab> contains the n x n band matrix A, see <band_row>, and zeros
   ! everywhere else. On exit, it contains the LU factors of A, and <ipiv> the
   ! row swaps, as from zgbtrf. <info> is the status of zgbtrf or plasma_zgbtrf,
   ! i.e., 0 on success and i > 0 if U(i,i) is zero.
   !****************************************************************************
   subroutine band_lu_factorise(n, kl, ku, ab, ipiv, info)

      use mp, only: mp_abort
#ifdef PLASMA
      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_zgbtrf, plasma_usable
#endif

      implicit none

      integer, intent(in) :: n, kl, ku
      complex, dimension(:, :), intent(in out) :: ab
      integer, dimension(:), intent(out) :: ipiv
      integer, intent(out) :: info

#ifdef PLASMA
      integer(c_int), dimension(n) :: ipiv_c
      integer :: nb_old
#endif

      !-------------------------------------------------------------------------

      if (size(ab, 1) < band_ldab(kl, ku)) call mp_abort('band_lu_factorise: ab has too few rows')

#ifdef PLASMA
      if (plasma_usable()) then
         call set_band_nb(nb_old)
         info = int(plasma_zgbtrf(int(n, c_int), int(n, c_int), int(kl, c_int), int(ku, c_int), &
                                  ab, int(size(ab, 1), c_int), ipiv_c))
         call reset_nb(nb_old)
         ipiv(1:n) = ipiv_c
         call plasma_band_factors_to_lapack(n, kl, ku, band_nb, ab, ipiv)
         return
      end if
#endif
      if (kind(ab) == kind((1.0d0, 0.0d0))) then
         call zgbtrf(n, n, kl, ku, ab, size(ab, 1), ipiv, info)
      else
         call cgbtrf(n, n, kl, ku, ab, size(ab, 1), ipiv, info)
      end if

   end subroutine band_lu_factorise

   !****************************************************************************
   !                     Solve with the LU factors of a band matrix
   !****************************************************************************
   ! Solves A x = b, given the factors <ab> and <ipiv> of A from
   ! <band_lu_factorise>. On exit, <b> contains x. This can be called from inside
   ! OpenMP parallel regions.
   !****************************************************************************
   subroutine band_lu_solve(n, kl, ku, ab, ipiv, b)

      implicit none

      integer, intent(in) :: n, kl, ku
      complex, dimension(:, :), intent(in) :: ab
      integer, dimension(:), intent(in) :: ipiv
      complex, dimension(:), intent(in out) :: b

      integer :: info

      !-------------------------------------------------------------------------

      if (kind(ab) == kind((1.0d0, 0.0d0))) then
         call zgbtrs('N', n, kl, ku, 1, ab, size(ab, 1), ipiv, b, size(b), info)
      else
         call cgbtrs('N', n, kl, ku, 1, ab, size(ab, 1), ipiv, b, size(b), info)
      end if

   end subroutine band_lu_solve

#ifdef PLASMA
   !****************************************************************************
   !              Store the band factors of PLASMA as those of LAPACK
   !****************************************************************************
   ! plasma_zgbtrf factorises panels of <nb> columns as a dense LU, so the row
   ! swaps of a panel are also applied to the earlier columns of L in the panel
   ! (see plasma/compute/pzgbtrf.c), and column j of L can have nonzeros down to
   ! row (last row of the panel) + kl. zgbtrf and zgbtrs only apply the swap of
   ! step k to the columns from k onwards, so column j of L holds the multipliers
   ! of step j, in rows j+1 to j+kl. Undoing the swaps of the steps k > j of the
   ! panel on column j, from the last step backwards, gives these multipliers, as
   ! the swaps of step k only exchange rows k and ipiv(k) >= k > j:
   !    P_last ... P_(j+1) l_j(LAPACK) = l_j(PLASMA)
   ! The rest of the factors (U, ipiv) is the same. With nb = 1 nothing changes.
   !****************************************************************************
   subroutine plasma_band_factors_to_lapack(n, kl, ku, nb, ab, ipiv)

      implicit none

      integer, intent(in) :: n, kl, ku, nb
      complex, dimension(:, :), intent(in out) :: ab
      integer, dimension(:), intent(in) :: ipiv

      integer :: kd, j0, j1, j, k, ip, ilast
      complex :: tmp

      !-------------------------------------------------------------------------

      ! The row of <ab> that contains the diagonal, A(i,j) is ab(kd + i - j, j)
      kd = kl + ku + 1

      do j0 = 1, n, nb
         j1 = min(n, j0 + nb - 1)
         ilast = min(n, j1 + kl)
         do j = j0, j1 - 1
            do k = j1, j + 1, -1
               ip = ipiv(k)
               if (ip /= k) then
                  tmp = ab(kd + k - j, j)
                  ab(kd + k - j, j) = ab(kd + ip - j, j)
                  ab(kd + ip - j, j) = tmp
               end if
            end do
            ! Below row j + kl, column j of L is zero now
            if (ilast > j + kl) ab(kd + kl + 1:kd + ilast - j, j) = 0.0
         end do
      end do

   end subroutine plasma_band_factors_to_lapack

   !****************************************************************************
   ! The PLASMA tile size is a global setting, so it is set to <band_nb> for every
   ! call of the band routines and reset afterwards, so that the dense routines
   ! keep the default tile size.
   !****************************************************************************
   subroutine set_band_nb(nb_old)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_get, plasma_set, PlasmaNb

      implicit none

      integer, intent(out) :: nb_old

      integer(c_int) :: info, val

      info = plasma_get(PlasmaNb, val)
      nb_old = val
      info = plasma_set(PlasmaNb, int(band_nb, c_int))

   end subroutine set_band_nb

   subroutine reset_nb(nb_old)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_set, PlasmaNb

      implicit none

      integer, intent(in) :: nb_old

      integer(c_int) :: info

      info = plasma_set(PlasmaNb, int(nb_old, c_int))

   end subroutine reset_nb
#endif

end module banded_lu
