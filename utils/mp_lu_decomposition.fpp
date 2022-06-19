module mp_lu_decomposition

#ifdef ISO_C_BINDING

   implicit none

   public :: lu_decomposition_local
   public :: lu_inverse_local
   public :: lu_matrix_multiply_local
   public :: lu_back_substitution_local

   interface lu_decomposition_local
!    module procedure lu_decomposition_local_real
      module procedure lu_decomposition_local_complex
   end interface

   interface lu_inverse_local
!    module procedure lu_inverse_local_real
      module procedure lu_inverse_local_complex
   end interface

   interface lu_matrix_multiply_local
!    module procedure lu_decomposition_local_real
      module procedure lu_matrix_multiply_local_complex
   end interface

   interface lu_back_substitution_local
!    module procedure lu_decomposition_local_real
      module procedure lu_back_substitution_local_complex
   end interface

contains

   subroutine lu_decomposition_local_complex(mp_comm, root, win, lu, idx, d)

      use mpi

      implicit none

      integer, intent(in) :: win, mp_comm, root
      complex, dimension(:, :), intent(in out) :: lu
      integer, dimension(:), intent(out) :: idx
      real, intent(out) :: d

      real, parameter :: zero = 1.0e-20
      real, dimension(size(lu, 1)) :: vv
      complex, dimension(size(lu, 2)) :: dum

      integer :: i, j, k, n, imax, lo, hi
      integer :: iproc, nproc, ierr
      real :: dmax, tmp

      n = size(lu, 1)

      call mpi_comm_size(mp_comm, nproc, ierr)
      call mpi_comm_rank(mp_comm, iproc, ierr)

      d = 1.0
      !the following is a loop to avoid copying entire matrix
      ! with (cabs(lu))
      do i = 1, n
         vv(i) = maxval(cabs(lu(i, :)))
      end do
      if (any(vv == 0.0)) &
         write (*, *) 'singular matrix in lu_decomposition on process ', iproc
      vv = 1.0 / vv
      do j = 1, n
         !divide up the work using row_limits
         call split_n_tasks(n - j, iproc, nproc, lo, hi, llim=j + 1)

         !pivot if needed
         dmax = -1.0
         do k = j, n
            tmp = vv(k) * abs(lu(k, j))
            if (tmp > dmax) then
               dmax = tmp
               imax = k
            end if
         end do

         if (iproc == root) then
            idx(j) = imax
            if (j /= imax) then
               dum = lu(imax, :)
               lu(imax, :) = lu(j, :)
               lu(j, :) = dum
               vv(imax) = vv(j)
               d = -d
            end if
            if (lu(j, j) == 0.0) lu(j, j) = zero
         else
            if (j /= imax) vv(imax) = vv(j)
         end if

         call mpi_win_fence(0, win, ierr)

         !get the lead multiplier
         do i = lo, hi
            lu(i, j) = lu(i, j) / lu(j, j)
         end do

         call mpi_win_fence(0, win, ierr)

         do k = lo, hi
            do i = j + 1, n
               lu(i, k) = lu(i, k) - lu(i, j) * lu(j, k)
            end do
         end do

         call mpi_win_fence(0, win, ierr)
      end do

   end subroutine lu_decomposition_local_complex

   subroutine lu_inverse_local_complex(mp_comm, win, lu, idx, inverse)

      use linear_solve, only: lu_back_substitution

      implicit none

      integer, intent(in) :: win, mp_comm
      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:, :), intent(out) :: inverse

      integer :: i, n, nproc, iproc, ierr
      integer :: lo, hi

      n = size(lu, 1)

      call mpi_comm_size(mp_comm, nproc, ierr)
      call mpi_comm_rank(mp_comm, iproc, ierr)

      call split_n_tasks(n, iproc, nproc, lo, hi)

      do i = lo, hi
         inverse(:, i) = 0
         inverse(i, i) = 1.0
      end do

      call mpi_win_fence(0, win, ierr)

      do i = lo, hi
         call lu_back_substitution(lu, idx, inverse(:, i))
      end do

      call mpi_win_fence(0, win, ierr)

   end subroutine lu_inverse_local_complex

   subroutine lu_matrix_multiply_local_complex(mp_comm, win, mat, b)

      implicit none

      integer, intent(in) :: win, mp_comm
      complex, dimension(:, :), intent(in) :: mat
      complex, dimension(:), intent(out) :: b
      complex, dimension(size(b)) :: a

      integer :: i, n, nproc, iproc
      integer :: lo, hi, ierr

      n = size(mat, 1)

      call mpi_comm_size(mp_comm, nproc, ierr)
      call mpi_comm_rank(mp_comm, iproc, ierr)

      call split_n_tasks(n, iproc, nproc, lo, hi)

      do i = lo, hi
         a(i) = sum(mat(i, :) * b(:))
      end do

      call mpi_win_fence(0, win, ierr)

      do i = lo, hi
         b(i) = a(i)
      end do

      call mpi_win_fence(0, win, ierr)

   end subroutine lu_matrix_multiply_local_complex

   subroutine lu_back_substitution_local_complex(mp_comm, win, lu, idx, b, &
                                                 time_parallel_streaming)

      use mpi
      use mp, only: mpicmplx
      use job_manage, only: time_message

      implicit none

      integer, intent(in) :: win, mp_comm
      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:), intent(in out) :: b
      real, dimension(:, :), intent(in out) :: time_parallel_streaming

      integer :: i, j, n, ii, ll, lo, hi
      integer :: iproc, nproc, ierr, aproc
      complex :: summ, dot, dot_local

      integer, parameter :: blocksize = 1

      call mpi_comm_size(mp_comm, nproc, ierr)
      call mpi_comm_rank(mp_comm, iproc, ierr)

      n = size(lu, 1)

      !perform pivoting on root node
      ii = 0
      if (iproc == 0) then
         do i = 1, n
            ll = idx(i)
            summ = b(ll)
            b(ll) = b(i)
            b(i) = summ
            if (ii /= 0 .and. summ /= 0.0) then
               ii = i
            end if
         end do
      end if

      call mpi_win_fence(0, win, ierr)

      aproc = nproc

      do i = ii + 1, n
         if (iproc==0) call time_message(.false., time_parallel_streaming(:, 4), ' dot')
         call split_n_tasks(n - i + 1, iproc, nproc, lo, hi, llim = i, &
                            blocksize = blocksize, aproc = aproc)
         do j = lo, hi
            b(j) = b(j) - lu(j, i - 1) * b(i - 1)
         enddo
         if (iproc==0) call time_message(.false., time_parallel_streaming(:, 4), ' dot')
         if (iproc == 0) call time_message(.false., time_parallel_streaming(:, 6), ' fence')
         call mpi_barrier(mp_comm, ierr)
         if (iproc == 0) call time_message(.false., time_parallel_streaming(:, 6), ' fence')
      end do

      call mpi_win_fence(0, win, ierr)

      do i = n, 1, -1
         summ = 1.0 / lu(i, i)
         if (iproc == 0) call time_message(.false., time_parallel_streaming(:, 5), ' dot')
         call split_n_tasks(i - 1, iproc, nproc, lo, hi, llim = 1, &
                            blocksize = blocksize, aproc = aproc)

         do j = lo, hi
            b(j) = b(j) - lu(j, i) * b(i) * summ
         enddo
         if (iproc == 0) call time_message(.false., time_parallel_streaming(:, 5), ' dot')

         if (iproc == 0) call time_message(.false., time_parallel_streaming(:, 7), ' fence')
         if (iproc == 0) b(i) = b(i) * summ
         call mpi_barrier(mp_comm, ierr)
         if (iproc == 0) call time_message(.false., time_parallel_streaming(:, 7), ' fence')

!        if (iproc==0) call time_message(.false., time_parallel_streaming(:, 5), ' dot')
!        call split_n_tasks(n - i, iproc, nproc, lo, hi, llim = i + 1, &
!                           blocksize = blocksize, aproc = aproc)
!        if (lo > hi) then
!           dot_local = 0
!        else
!           dot_local = dot_product(conjg(lu(i, lo:hi)), b(lo:hi))
!        end if
!        if (iproc==0) call time_message(.false., time_parallel_streaming(:, 5), ' dot')

!        if (iproc==0) call time_message(.false., time_parallel_streaming(:, 5), ' reduce')
!        if(aproc.gt.1) then
!           call mpi_reduce(dot_local, dot, 1, mpicmplx, MPI_SUM, 0, mp_comm, ierr)
!        else
!           dot = dot_local
!        endif
!        if (iproc==0) call time_message(.false., time_parallel_streaming(:, 5), ' reduce')

!        if (iproc == 0) b(i) = (b(i) - dot) / lu(i, i)
!        if (iproc==0) call time_message(.false., time_parallel_streaming(:, 7), ' fence')
!        if(aproc.gt.1) call mpi_win_fence(0, win, ierr)
!        if(aproc.gt.1) call mpi_barrier(mp_comm, ierr)
!        if (iproc==0) call time_message(.false., time_parallel_streaming(:, 7), ' fence')
      end do

      call mpi_win_fence(0, win, ierr)

   end subroutine lu_back_substitution_local_complex

   subroutine split_n_tasks(n, iproc, nproc, lo, hi, llim, blocksize, aproc)

      implicit none

      integer, intent(in) :: n, iproc, nproc
      integer, intent(out) :: lo, hi
      integer, optional, intent(in) :: llim
      integer, optional, intent(in) :: blocksize ! minimum amount of cells per node
      integer, optional, intent(out) :: aproc !how many processes have work

      integer :: n_div, n_mod, llim_l, blocksize_l

      llim_l = 1
      if (present(llim)) llim_l = llim

      blocksize_l = 1
      if (present(blocksize)) blocksize_l = blocksize

      n_div = n / nproc
      n_mod = mod(n, nproc)

      if (n_div .lt. blocksize_l) then
         lo = min(iproc * blocksize_l + llim_l, n + llim_l)
         hi = min(lo + blocksize_l - 1, n + llim_l - 1)
         if (present(aproc)) then
            aproc = n / blocksize_l
            if (aproc * blocksize_l .lt. n) aproc = aproc + 1
         endif
      else
         lo = iproc * n_div + min(iproc, n_mod) + llim_l
         hi = lo + n_div - 1
         if (iproc < n_mod) hi = hi + 1
         if (present(aproc)) then
            if(n_div .gt. 0) then
               aproc = nproc
            else
               aproc = n_mod
            endif
         endif
      endif
      


   end subroutine split_n_tasks


#endif

end module mp_lu_decomposition
