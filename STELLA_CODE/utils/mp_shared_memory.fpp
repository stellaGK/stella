# include "define.inc"

!###############################################################################
!                  SHARED-MEMORY DOMAINS AND SHARED WINDOWS
!###############################################################################
!
! Large read-only arrays that are the same on every MPI process, such as the LU
! factors of the response matrix, are stored once per shared-memory domain in an
! MPI shared-memory window, instead of once per process. One process per domain
! (the domain root, rank 0 of <comm_domain>) writes the array, and all processes
! of the domain read it.
!
! The domains are chosen with <sharing>:
!     'numa'  One domain per NUMA domain (default). The NUMA domain of a process
!             is found from the cores it may run on (its affinity mask) and the
!             cores of every NUMA domain in /sys/devices/system/node. A process
!             whose cores lie in more than one NUMA domain (e.g. a process that
!             is not bound to cores) keeps its own copy, so that no process
!             reads the memory of another socket in every time step. This does
!             not depend on how stella is started (mpirun or srun), and the
!             domains may have different numbers of processes.
!     'node'  One domain per node (all processes that can share memory).
!     'none'  Every process keeps its own copy (one domain per process).
!
! When the NUMA domains can not be read (not Linux, no /sys), 'numa' falls back
! to 'none'.
!
!###############################################################################
module mp_shared_memory

#ifndef MPIINC
   use mpi
#endif

   implicit none
   private

   public :: init_shared_memory_domains, finish_shared_memory_domains
   public :: allocate_shared_window, free_shared_window, fence_shared_window
   public :: broadcast_in_domain, broadcast_to_domain_roots
   public :: comm_domain, comm_domain_roots
   public :: domain_rank, domain_size, domain_root
   public :: ndomains, min_domain_size, max_domain_size, nprivate_processes

#ifdef MPIINC
   include 'mpif.h'
#endif

   ! Communicator of the processes of this domain, and of the roots of all domains
   ! (MPI_COMM_NULL on the processes that are not a domain root)
   integer :: comm_domain = MPI_COMM_NULL
   integer :: comm_domain_roots = MPI_COMM_NULL

   ! Rank and size of <comm_domain>; the root (rank 0) writes the shared arrays
   integer :: domain_rank = 0, domain_size = 1
   logical :: domain_root = .true.

   ! Summary over all processes, for the output
   integer :: ndomains = 0, min_domain_size = 0, max_domain_size = 0
   integer :: nprivate_processes = 0

   logical :: initialised = .false.

contains

   !============================================================================
   !                     Split the processes into domains
   !============================================================================
   ! <comm> is the communicator of the processes that build the shared arrays
   ! together (mp_comm), and <sharing> is 'numa', 'node' or 'none'.
   !============================================================================
   subroutine init_shared_memory_domains(comm, sharing)

      implicit none

      integer, intent(in) :: comm
      character(*), intent(in) :: sharing

      integer :: comm_node, node_rank, node_size, rank, color, numa, ierr
      integer, dimension(2) :: sizes
      integer :: is_root, is_private

      ! -------------------------------------------------------------------------

      if (initialised) return
      initialised = .true.

      call mpi_comm_rank(comm, rank, ierr)

      ! All processes that can share memory with this one
      call mpi_comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, comm_node, ierr)
      call mpi_comm_rank(comm_node, node_rank, ierr)
      call mpi_comm_size(comm_node, node_size, ierr)

      ! The domain of this process within its node: <color> 0 for the whole node,
      ! the NUMA domain, or a color of its own for a private copy
      is_private = 0
      select case (trim(sharing))
      case ('node')
         color = 0
      case ('none')
         color = node_rank
         is_private = 1
      case default
         numa = numa_domain_of_process()
         if (numa >= 0) then
            color = numa
         else
            color = 100000 + node_rank
            is_private = 1
         end if
      end select
      call mpi_comm_split(comm_node, color, node_rank, comm_domain, ierr)
      call mpi_comm_free(comm_node, ierr)
      call mpi_comm_rank(comm_domain, domain_rank, ierr)
      call mpi_comm_size(comm_domain, domain_size, ierr)
      domain_root = (domain_rank == 0)

      ! The roots of all domains, which exchange the shared arrays
      color = MPI_UNDEFINED
      if (domain_root) color = 0
      call mpi_comm_split(comm, color, rank, comm_domain_roots, ierr)

      ! Number and sizes of the domains, for the output
      is_root = merge(1, 0, domain_root)
      call mpi_allreduce(is_root, ndomains, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
      call mpi_allreduce(is_private, nprivate_processes, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
      sizes = [-domain_size, domain_size]
      call mpi_allreduce(MPI_IN_PLACE, sizes, 2, MPI_INTEGER, MPI_MAX, comm, ierr)
      min_domain_size = -sizes(1); max_domain_size = sizes(2)

   end subroutine init_shared_memory_domains

   !============================================================================
   !              The NUMA domain of the cores this process may use
   !============================================================================
   ! Returns -1 if the cores lie in more than one NUMA domain, or if the NUMA
   ! domains can not be read.
   !============================================================================
   integer function numa_domain_of_process() result(numa)

#ifdef __linux__
      use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_int64_t
#endif

      implicit none

#ifdef __linux__
      interface
         function sched_getaffinity(pid, cpusetsize, mask) bind(C, name='sched_getaffinity')
            import :: c_int, c_size_t, c_int64_t
            integer(c_int), value :: pid
            integer(c_size_t), value :: cpusetsize
            integer(c_int64_t), dimension(*) :: mask
            integer(c_int) :: sched_getaffinity
         end function sched_getaffinity
      end interface

      ! Affinity mask for up to 64 * 64 = 4096 cores
      integer, parameter :: nwords = 64
      integer(c_int64_t), dimension(nwords) :: mask
      integer :: node, cpu, nmissing, numa_of_cpu
      logical, dimension(0:64 * nwords - 1) :: allowed, in_node
      logical :: exists
      character(len=64) :: filename

      ! -------------------------------------------------------------------------

      numa = -1
      mask = 0
      if (sched_getaffinity(0_c_int, int(8 * nwords, c_size_t), mask) /= 0) return
      do cpu = 0, 64 * nwords - 1
         allowed(cpu) = btest(mask(cpu / 64 + 1), mod(cpu, 64))
      end do
      if (.not. any(allowed)) return

      ! Find the NUMA domain of every allowed core (node numbers may have gaps)
      numa_of_cpu = -1; nmissing = 0
      do node = 0, 1023
         write (filename, '(A, I0, A)') '/sys/devices/system/node/node', node, '/cpulist'
         inquire (file=trim(filename), exist=exists)
         if (.not. exists) then
            nmissing = nmissing + 1
            if (nmissing > 16) exit
            cycle
         end if
         if (.not. read_cpulist(trim(filename), in_node)) return
         if (.not. any(allowed .and. in_node)) cycle
         ! The cores of this process lie in more than one NUMA domain
         if (numa_of_cpu >= 0 .or. any(allowed .and. .not. in_node)) return
         numa_of_cpu = node
      end do
      numa = numa_of_cpu
#else
      numa = -1
#endif

   end function numa_domain_of_process

   !============================================================================
   !        Read a list of cores such as "0-63,128-191" from <filename>
   !============================================================================
   logical function read_cpulist(filename, in_list) result(success)

      implicit none

      character(*), intent(in) :: filename
      logical, dimension(0:), intent(out) :: in_list

      character(len=4096) :: line
      integer :: unit, ios, start, comma, dash, first, last

      ! -------------------------------------------------------------------------

      success = .false.
      in_list = .false.
      open (newunit=unit, file=filename, action='read', status='old', iostat=ios)
      if (ios /= 0) return
      line = ''
      read (unit, '(A)', iostat=ios) line
      close (unit)
      if (ios > 0) return

      start = 1
      do while (start <= len_trim(line))
         comma = index(line(start:), ',')
         if (comma == 0) then
            comma = len_trim(line) + 1
         else
            comma = start + comma - 1
         end if
         dash = index(line(start:comma - 1), '-')
         if (dash == 0) then
            read (line(start:comma - 1), *, iostat=ios) first
            last = first
         else
            dash = start + dash - 1
            read (line(start:dash - 1), *, iostat=ios) first
            if (ios == 0) read (line(dash + 1:comma - 1), *, iostat=ios) last
         end if
         if (ios /= 0) return
         if (first < 0 .or. last > ubound(in_list, 1)) return
         in_list(first:last) = .true.
         start = comma + 1
      end do
      success = .true.

   end function read_cpulist

   !============================================================================
   !            Allocate a shared window of <nbytes> bytes on the root
   !============================================================================
   ! Every process of the domain gets the address of the window in <base>.
   ! Returns MPI_WIN_NULL if the MPI library does not use the unified memory
   ! model, in which case the caller keeps a private copy.
   !============================================================================
   subroutine allocate_shared_window(nbytes, window, base)

      use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr

      implicit none

      integer(kind=MPI_ADDRESS_KIND), intent(in) :: nbytes
      integer, intent(out) :: window
      type(c_ptr), intent(out) :: base

      integer(kind=MPI_ADDRESS_KIND) :: size, memory_model
      integer :: disp_unit, ierr
      logical :: flag

      ! -------------------------------------------------------------------------

      ! Only the root allocates memory, so that it lies on the root's NUMA domain
      size = 0
      if (domain_root) size = max(nbytes, 1_MPI_ADDRESS_KIND)
      disp_unit = 1
      call mpi_win_allocate_shared(size, disp_unit, MPI_INFO_NULL, comm_domain, base, window, ierr)
      if (.not. domain_root) call mpi_win_shared_query(window, 0, size, disp_unit, base, ierr)

      ! Loads and stores of other processes are only visible after a fence
      ! (or a synchronisation) in the unified memory model
      call mpi_win_get_attr(window, MPI_WIN_MODEL, memory_model, flag, ierr)
      if (.not. flag .or. memory_model /= MPI_WIN_UNIFIED) then
         call mpi_win_free(window, ierr)
         window = MPI_WIN_NULL
         base = c_null_ptr
         return
      end if
      call mpi_win_fence(0, window, ierr)

   end subroutine allocate_shared_window

   !============================================================================
   !       Make the writes into the window visible to the whole domain
   !============================================================================
   subroutine fence_shared_window(window)

      implicit none

      integer, intent(in) :: window
      integer :: ierr

      if (window /= MPI_WIN_NULL) call mpi_win_fence(0, window, ierr)

   end subroutine fence_shared_window

   !============================================================================
   !                             Free a window
   !============================================================================
   subroutine free_shared_window(window)

      implicit none

      integer, intent(in out) :: window
      integer :: ierr

      if (window /= MPI_WIN_NULL) call mpi_win_free(window, ierr)
      window = MPI_WIN_NULL

   end subroutine free_shared_window

   !============================================================================
   !          Broadcast from the root of the domain to its processes
   !============================================================================
   subroutine broadcast_in_domain(x)

      implicit none

      integer, dimension(:), intent(in out) :: x
      integer :: ierr

      if (domain_size > 1) call mpi_bcast(x, size(x), MPI_INTEGER, 0, comm_domain, ierr)

   end subroutine broadcast_in_domain

   !============================================================================
   !   Broadcast a complex matrix from the first root to the roots of all domains
   !============================================================================
   subroutine broadcast_to_domain_roots(x, mpicmplx)

      implicit none

      complex, dimension(:, :), intent(in out) :: x
      integer, intent(in) :: mpicmplx
      integer :: ierr

      if (domain_root) call mpi_bcast(x, size(x), mpicmplx, 0, comm_domain_roots, ierr)

   end subroutine broadcast_to_domain_roots

   !============================================================================
   !                          Free the communicators
   !============================================================================
   subroutine finish_shared_memory_domains

      implicit none

      integer :: ierr

      if (.not. initialised) return
      if (comm_domain /= MPI_COMM_NULL) call mpi_comm_free(comm_domain, ierr)
      if (comm_domain_roots /= MPI_COMM_NULL) call mpi_comm_free(comm_domain_roots, ierr)
      comm_domain = MPI_COMM_NULL
      comm_domain_roots = MPI_COMM_NULL
      domain_rank = 0; domain_size = 1; domain_root = .true.
      initialised = .false.

   end subroutine finish_shared_memory_domains

end module mp_shared_memory
