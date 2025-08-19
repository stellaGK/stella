module response_matrix

   use netcdf
   use mpi
#ifdef ISO_C_BINDING
   use, intrinsic :: iso_c_binding, only: c_intptr_t
#endif
   use debug_flags, only: debug => response_matrix_debug

   implicit none

   public :: init_response_matrix, finish_response_matrix
   public :: read_response_matrix
   public :: response_matrix_initialized

   private

   logical :: response_matrix_initialized = .false.
   integer, parameter :: mat_unit = 70
#ifdef ISO_C_BINDING
   integer(c_intptr_t) :: cur_pos
#endif
   character(100) :: message_dgdphi, message_QN, message_lu
   real, dimension(2) :: time_dgdphi
   real, dimension(2) :: time_QN
   real, dimension(2) :: time_lu

contains

   subroutine init_response_matrix

      use linear_solve, only: lu_decomposition
      use arrays_fields, only: response_matrix
      use stella_layouts, only: iv_idx, is_idx
      use parameters_kxky_grids, only: naky
      use mp, only: proc0
      use parameters_numerical, only: mat_gen
#ifdef ISO_C_BINDING
      use arrays_fields, only: response_window
#endif

      implicit none

#ifdef ISO_C_BINDING
      integer :: ierr
#endif
      debug = (debug .and. proc0)

      if (debug) call write_response_matrix_message
      call setup_response_matrix_timings
      call setup_response_matrix_file_io

      if (response_matrix_initialized) return
      response_matrix_initialized = .true.

      if (.not. allocated(response_matrix)) allocate (response_matrix(naky))

#ifdef ISO_C_BINDING
      call setup_shared_memory_window
#endif

      call construct_response_matrix

#ifdef ISO_C_BINDING
      call mpi_win_fence(0, response_window, ierr)
#endif

      if (proc0 .and. mat_gen) then
         close (unit=mat_unit)
      end if

      if (debug) then
         write (*, '(A)') "    ############################################################"
         write (*, '(A)') " "
      end if

   end subroutine init_response_matrix

   subroutine write_response_matrix_message

      write (*, *) " "
      write (*, '(A)') "    ############################################################"
      write (*, '(A)') "                         RESPONSE MATRIX"
      write (*, '(A)') "    ############################################################"

   end subroutine write_response_matrix_message

   subroutine setup_response_matrix_timings

      implicit none

      message_dgdphi = '     calculate dgdphi: '
      message_QN = '     calculate QN:     '
      message_lu = '     calculate LU:     '
      time_dgdphi = 0.0
      time_QN = 0.0
      time_lu = 0.0

   end subroutine setup_response_matrix_timings

   subroutine setup_response_matrix_file_io

      use mp, only: proc0, job
      use parameters_numerical, only: mat_gen
      use system_fortran, only: systemf
      use parameters_kxky_grids, only: naky

      implicit none

      character(len=15) :: job_str
      character(len=100) :: file_name

      ! All matrices handled by processor i_proc and job are stored
      ! on a single file named: response_mat_job.iproc
      if (proc0 .and. mat_gen) then
         call systemf('mkdir -p mat')

         write (job_str, '(I1.1)') job
         file_name = './mat/response_mat_'//trim(job_str)

         open (unit=mat_unit, status='replace', file=file_name, &
               position='rewind', action='write', form='unformatted')
         write (unit=mat_unit) naky
      end if

   end subroutine setup_response_matrix_file_io

#ifdef ISO_C_BINDING
   subroutine setup_shared_memory_window

      use mpi
      use, intrinsic :: iso_c_binding, only: c_intptr_t
      use mp, only: sgproc0, real_size
      use mp, only: create_shared_memory_window
      use arrays_fields, only: response_window
      use fields, only: nfields
      use parameters_kxky_grids, only: naky
      use extended_zgrid, only: neigen, nsegments, nzed_segment
      use extended_zgrid, only: periodic

      implicit none

      integer(kind=MPI_ADDRESS_KIND) :: win_size
      integer :: iky, ie
      integer :: nresponse

      ! Create a single shared memory window for all the response matrices and
      ! permutation arrays.
      ! Creating a window for each matrix/array would lead to performance
      ! degradation on some clusters
      if (response_window == MPI_WIN_NULL) then
         win_size = 0
         if (sgproc0) then
            do iky = 1, naky
               do ie = 1, neigen(iky)
                  if (periodic(iky)) then
                     nresponse = (nsegments(ie, iky) * nzed_segment) * nfields
                  else
                     nresponse = (nsegments(ie, iky) * nzed_segment + 1) * nfields
                  end if
                  win_size = win_size &
                             + int(nresponse, MPI_ADDRESS_KIND) * 4_MPI_ADDRESS_KIND &
                             + int(nresponse**2, MPI_ADDRESS_KIND) * 2 * real_size
               end do
            end do
         end if

         call create_shared_memory_window(win_size, response_window, cur_pos)
      end if

   end subroutine setup_shared_memory_window
#endif

   subroutine construct_response_matrix

      use mp, only: proc0
      use job_manage, only: time_message
      use parameters_numerical, only: mat_gen
      use arrays_fields, only: response_matrix
      use parameters_kxky_grids, only: naky
      use extended_zgrid, only: neigen
#ifdef ISO_C_BINDING
      use arrays_fields, only: response_window
#endif

      implicit none

      integer :: iky, ie
#ifdef ISO_C_BINDING
      integer :: ierr
#endif

      ! for a given ky and set of connected kx values
      ! give a unit impulse to phi at each zed location
      ! in the extended domain and solve for h(zed_extended,(vpa,mu,s))

      do iky = 1, naky

         if (proc0 .and. mat_gen) then
            write (unit=mat_unit) iky, neigen(iky)
         end if

         ! the response matrix for each ky has neigen(ky)
         ! independent sets of connected kx values
         if (.not. associated(response_matrix(iky)%eigen)) &
            allocate (response_matrix(iky)%eigen(neigen(iky)))

         if (debug) call time_message(.false., time_dgdphi, message_dgdphi)

         call calculate_vspace_integrated_response(iky)

         !DSO - This ends parallelization over velocity space.
         !      At this point every processor has int dv dgdphi for a given ky
         !      and so the quasineutrality solve and LU decomposition can be
         !      parallelized locally if need be.
         !      This is preferable to parallelization over ky as the LU
         !      decomposition (and perhaps QN) will be dominated by the
         !      ky with the most connections

         if (debug) then
            call time_message(.true., time_dgdphi, message_dgdphi)
            call time_message(.false., time_QN, message_QN)
         end if

#ifdef ISO_C_BINDING
         call mpi_win_fence(0, response_window, ierr)
#endif

         call apply_field_solve_to_finish_response_matrix(iky)

#ifdef ISO_C_BINDING
         call mpi_win_fence(0, response_window, ierr)
#endif

         if (debug) then
            call time_message(.true., time_QN, message_QN)
            call time_message(.false., time_lu, message_lu)
         end if

         call lu_decompose_response_matrix(iky)

         if (proc0 .and. debug) then
            call time_message(.true., time_lu, message_lu)
         end if

         time_dgdphi = 0
         time_QN = 0
         time_lu = 0

         do ie = 1, neigen(iky)
            if (proc0 .and. mat_gen) then
               write (unit=mat_unit) response_matrix(iky)%eigen(ie)%idx
               write (unit=mat_unit) response_matrix(iky)%eigen(ie)%zloc
            end if
         end do

      end do

   end subroutine construct_response_matrix

   subroutine calculate_vspace_integrated_response(iky)

      use mp, only: proc0
      use parameters_numerical, only: mat_gen
      use physics_parameters, only: include_apar, include_bpar
      use extended_zgrid, only: neigen, ikxmod
      use extended_zgrid, only: nsegments, nzed_segment
      use extended_zgrid, only: periodic
      use extended_zgrid, only: iz_low, iz_up
      use stella_layouts, only: vmu_lo
      use fields, only: nfields

      implicit none

      integer, intent(in) :: iky

      integer :: ie, idx, ikx, iseg
      integer :: iz, izl_offset, izup
      integer :: nz_ext, nresponse, nresponse_per_field
      complex, dimension(:, :), allocatable :: gext
      complex, dimension(:), allocatable :: phi_ext, apar_ext, bpar_ext

      ! loop over the sets of connected kx values
      do ie = 1, neigen(iky)

         ! number of zeds x number of segments on extended zed domain
         nz_ext = nsegments(ie, iky) * nzed_segment + 1

         ! treat zonal mode specially to avoid double counting
         ! as it is periodic
         if (periodic(iky)) then
            nresponse_per_field = nz_ext - 1
         else
            nresponse_per_field = nz_ext
         end if
         nresponse = nresponse_per_field * nfields

         if (proc0 .and. mat_gen) then
            write (unit=mat_unit) ie, nresponse
         end if

         call setup_response_matrix_zloc_idx(iky, ie, nresponse)

         allocate (gext(nz_ext, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (phi_ext(nz_ext))
         allocate (apar_ext(nz_ext))
         allocate (bpar_ext(nz_ext))

         ! idx is the index in the extended zed domain
         ! that we are giving a unit impulse
         idx = 0

         ! loop over segments, starting with 1
         ! first segment is special because it has
         ! one more unique zed value than all others
         ! since domain is [z0-pi:z0+pi], including both endpoints
         ! i.e., one endpoint is shared with the previous segment
         iseg = 1
         ! ikxmod gives the kx corresponding to iseg,ie,iky
         ikx = ikxmod(iseg, ie, iky)
         izl_offset = 0
         ! avoid double-counting of periodic points for zonal mode (and other periodic modes)
         if (periodic(iky)) then
            izup = iz_up(iseg) - 1
         else
            izup = iz_up(iseg)
         end if
         ! no need to obtain response to impulses at negative kx values
         do iz = iz_low(iseg), izup
            idx = idx + 1
            call get_dpdf_dphi_matrix_column(iky, ie, idx, nz_ext, nresponse_per_field, phi_ext, apar_ext, bpar_ext, gext)
            if (include_apar) call get_dpdf_dapar_matrix_column(iky, ie, idx, nz_ext, nresponse_per_field, phi_ext, apar_ext, bpar_ext, gext)
            if (include_bpar) call get_dpdf_dbpar_matrix_column(iky, ie, idx, nz_ext, nresponse_per_field, phi_ext, apar_ext, bpar_ext, gext)
         end do
         ! once we have used one segment, remaining segments
         ! have one fewer unique zed point
         izl_offset = 1
         if (nsegments(ie, iky) > 1) then
            do iseg = 2, nsegments(ie, iky)
               ikx = ikxmod(iseg, ie, iky)
               do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
                  idx = idx + 1
                  call get_dpdf_dphi_matrix_column(iky, ie, idx, nz_ext, nresponse_per_field, phi_ext, apar_ext, bpar_ext, gext)
                  if (include_apar) call get_dpdf_dapar_matrix_column(iky, ie, idx, nz_ext, nresponse_per_field, phi_ext, apar_ext, bpar_ext, gext)
                  if (include_bpar) call get_dpdf_dbpar_matrix_column(iky, ie, idx, nz_ext, nresponse_per_field, phi_ext, apar_ext, bpar_ext, gext)
               end do
               if (izl_offset == 0) izl_offset = 1
            end do
         end if
         deallocate (gext, phi_ext, apar_ext, bpar_ext)
      end do

   end subroutine calculate_vspace_integrated_response

   subroutine setup_response_matrix_zloc_idx(iky, ie, nresponse)

#ifdef ISO_C_BINDING
      use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
      use mp, only: nbytes_real
#endif
      use arrays_fields, only: response_matrix

      implicit none

      integer, intent(in) :: iky, ie, nresponse

#ifdef ISO_C_BINDING
      type(c_ptr) :: cptr

      !exploit MPIs shared memory framework to reduce memory consumption of
      !response matrices

      if (.not. associated(response_matrix(iky)%eigen(ie)%zloc)) then
         cptr = transfer(cur_pos, cptr)
         call c_f_pointer(cptr, response_matrix(iky)%eigen(ie)%zloc, (/nresponse, nresponse/))
         cur_pos = cur_pos + nresponse**2 * 2 * nbytes_real
      end if

      if (.not. associated(response_matrix(iky)%eigen(ie)%idx)) then
         cptr = transfer(cur_pos, cptr)
         call c_f_pointer(cptr, response_matrix(iky)%eigen(ie)%idx, (/nresponse/))
         cur_pos = cur_pos + nresponse * 4
      end if
#else
      ! for each ky and set of connected kx values,
      ! must have a response matrix that is N x N
      ! with N = number of zeds per 2pi segment x number of 2pi segments
      if (.not. associated(response_matrix(iky)%eigen(ie)%zloc)) &
         allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse, nresponse))

      ! response_matrix%idx is needed to keep track of permutations
      ! to the response matrix made during LU decomposition
      ! it will be input to LU back substitution during linear solve
      if (.not. associated(response_matrix(iky)%eigen(ie)%idx)) &
         allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))
#endif

   end subroutine setup_response_matrix_zloc_idx

   subroutine apply_field_solve_to_finish_response_matrix(iky)

#ifdef ISO_C_BINDING
      use mp, only: sgproc0
#endif
      use physics_parameters, only: include_apar, include_bpar
      use extended_zgrid, only: neigen
      use extended_zgrid, only: nsegments, nzed_segment
      use extended_zgrid, only: periodic
      use arrays_fields, only: response_matrix

      implicit none

      integer, intent(in) :: iky

      integer :: ie, idx, offset_apar, offset_bpar
      integer :: nz_ext, nresponse
      character(5) :: dist
      complex, dimension(:), allocatable :: phi_ext, apar_ext, bpar_ext

      ! solve quasineutrality
      ! for local stella, this is a diagonal process, but global stella
      ! may require something more sophisticated
      dist = 'g'

      ! loop over the sets of connected kx values
      do ie = 1, neigen(iky)
#ifdef ISO_C_BINDING
         if (sgproc0) then
#endif
            ! number of zeds x number of segments
            nz_ext = nsegments(ie, iky) * nzed_segment + 1

            ! treat zonal mode specially to avoid double counting
            ! as it is periodic
            if (periodic(iky)) then
               nresponse = nz_ext - 1
            else
               nresponse = nz_ext
            end if

            allocate (phi_ext(nz_ext))
            allocate (apar_ext(nz_ext))
            allocate (bpar_ext(nz_ext))
            !> set up offset_apar and offset_bpar consistently
            !> so that the array slices below are consistent with
            !> the size of the response matrix
            if (include_apar) then
               offset_apar = nresponse
            else
               offset_apar = 0
            end if
            if (include_bpar) then
               offset_bpar = offset_apar + nresponse
            else
               offset_bpar = 0
            end if
            ! obtain the response matrix entries due to unit impulses in phi;
            ! this accounts for terms appearing both in quasineutrality and parallel ampere
            do idx = 1, nresponse
               phi_ext(nz_ext) = 0.0
               phi_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx)
               if (include_apar) then
                  apar_ext(nz_ext) = 0.0
                  apar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx)
               end if
               if (include_bpar) then
                  bpar_ext(nz_ext) = 0.0
                  bpar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx)
               end if
               call get_fields_for_response_matrix(phi_ext, apar_ext, bpar_ext, iky, ie, dist)

               ! next need to create column in response matrix from phi_ext and apar_ext
               ! negative sign because matrix to be inverted in streaming equation
               ! is identity matrix - response matrix
               ! add in contribution from identity matrix
               phi_ext(idx) = phi_ext(idx) - 1.0
               response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx) = -phi_ext(:nresponse)
               if (include_apar) response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx) = -apar_ext(:nresponse)
               if (include_bpar) response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx) = -bpar_ext(:nresponse)
            end do

            if (include_apar) then
               ! obtain the response matrix entries due to unit impulses in apar;
               ! this accounts for terms appearing both in quasineutrality and parallel ampere
               do idx = 1, nresponse
                  phi_ext(nz_ext) = 0.0
                  phi_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx + offset_apar)
                  apar_ext(nz_ext) = 0.0
                  apar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx + offset_apar)
                  if (include_bpar) then
                     bpar_ext(nz_ext) = 0.0
                     bpar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx + offset_apar)
                  end if
                  call get_fields_for_response_matrix(phi_ext, apar_ext, bpar_ext, iky, ie, dist)

                  ! next need to create column in response matrix from phi_ext and apar_ext
                  ! negative sign because matrix to be inverted in streaming equation
                  ! is identity matrix - response matrix
                  ! add in contribution from identity matrix for diagonal entries
                  apar_ext(idx) = apar_ext(idx) - 1.0
                  response_matrix(iky)%eigen(ie)%zloc(:nresponse, offset_apar + idx) = -phi_ext(:nresponse)
                  response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, offset_apar + idx) = -apar_ext(:nresponse)
                  if (include_bpar) response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, offset_apar + idx) = -bpar_ext(:nresponse) 
               end do
            end if
            
            if (include_bpar) then
               ! obtain the response matrix entries due to unit impulses in bpar;
               ! this accounts for terms appearing both in quasineutrality and parallel ampere
               do idx = 1, nresponse
                  phi_ext(nz_ext) = 0.0
                  phi_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx + offset_bpar)
                  if (include_apar) then
                     apar_ext(nz_ext) = 0.0
                     apar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx + offset_bpar)
                  end if
                  bpar_ext(nz_ext) = 0.0
                  bpar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx + offset_bpar)
                  call get_fields_for_response_matrix(phi_ext, apar_ext, bpar_ext, iky, ie, dist)

                  ! next need to create column in response matrix from phi_ext and apar_ext
                  ! negative sign because matrix to be inverted in streaming equation
                  ! is identity matrix - response matrix
                  ! add in contribution from identity matrix for diagonal entries
                  bpar_ext(idx) = bpar_ext(idx) - 1.0
                  response_matrix(iky)%eigen(ie)%zloc(:nresponse, offset_bpar + idx) = -phi_ext(:nresponse)
                  if (include_apar) response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, offset_bpar + idx) = -apar_ext(:nresponse)
                  response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, offset_bpar + idx) = -bpar_ext(:nresponse) 
               end do
            end if

            deallocate (phi_ext, apar_ext, bpar_ext)
#ifdef ISO_C_BINDING
         end if
#endif
      end do

   end subroutine apply_field_solve_to_finish_response_matrix

   subroutine lu_decompose_response_matrix(iky)

#ifdef ISO_C_BINDING
      use mp, only: sgproc0
#endif
      use mp, only: mp_abort
      use arrays_fields, only: response_matrix
      use parameters_numerical, only: lu_option_switch
      use parameters_numerical, only: lu_option_none, lu_option_local, lu_option_global
      use extended_zgrid, only: neigen
      use linear_solve, only: lu_decomposition

      implicit none

      integer, intent(in) :: iky

      integer :: ie
      real :: dum

      ! now we have the full response matrix. Finally, perform its LU decomposition
      select case (lu_option_switch)
      case (lu_option_global)
         call parallel_LU_decomposition_global(iky)
      case (lu_option_local)
#ifdef ISO_C_BINDING
         call parallel_LU_decomposition_local(iky)
#else
         call mp_abort('stella must be built with HAS_ISO_BINDING in order to use local parallel LU decomposition.')
#endif
      case default
         do ie = 1, neigen(iky)
#ifdef ISO_C_BINDING
            if (sgproc0) then
#endif
               ! now that we have the reponse matrix for this ky and set of connected kx values
               !get the LU decomposition so we are ready to solve the linear system
               call lu_decomposition(response_matrix(iky)%eigen(ie)%zloc, &
                                     response_matrix(iky)%eigen(ie)%idx, dum)

#ifdef ISO_C_BINDING
            end if
#endif
         end do
      end select

   end subroutine lu_decompose_response_matrix

   subroutine read_response_matrix

      use arrays_fields, only: response_matrix
      use common_types, only: response_matrix_type
      use parameters_kxky_grids, only: naky
      use extended_zgrid, only: neigen
      use extended_zgrid, only: nsegments
      use extended_zgrid, only: nzed_segment
      use extended_zgrid, only: periodic
      use mp, only: proc0, job, broadcast, mp_abort
      use fields, only: nfields

      implicit none

      integer :: iky, ie, nz_ext
      integer :: iky_dump, neigen_dump, naky_dump, nresponse_dump
      integer :: nresponse, nresponse_per_field
      character(len=15) :: job_str
      character(len=100) :: file_name
      integer :: ie_dump, istat
      logical, parameter :: debug = .false.

!   All matrices handled for the job i_job are read
!   from a single file named: responst_mat.ijob by that
!   jobs root process

      if (proc0) then
         write (job_str, '(I1.1)') job
         file_name = './mat/response_mat.'//trim(job_str)

         open (unit=mat_unit, status='old', file=file_name, &
               action='read', form='unformatted', iostat=istat)
         if (istat /= 0) then
            print *, 'Error opening response_matrix by root processor for job ', job_str
         end if

         read (unit=mat_unit) naky_dump
         if (naky /= naky_dump) call mp_abort('mismatch in naky and naky_dump')
      end if

      if (.not. allocated(response_matrix)) allocate (response_matrix(naky))

      do iky = 1, naky
         if (proc0) then
            read (unit=mat_unit) iky_dump, neigen_dump
            if (iky_dump /= iky .or. neigen_dump /= neigen(iky)) &
               call mp_abort('mismatch in iky_dump/neigen_dump')
         end if

         if (.not. associated(response_matrix(iky)%eigen)) &
            allocate (response_matrix(iky)%eigen(neigen(iky)))

         ! loop over the sets of connected kx values
         do ie = 1, neigen(iky)
            ! number of zeds x number of segments
            nz_ext = nsegments(ie, iky) * nzed_segment + 1

            ! treat zonal mode specially to avoid double counting
            ! as it is periodic
            if (periodic(iky)) then
               nresponse_per_field = nz_ext - 1
            else
               nresponse_per_field = nz_ext
            end if
            nresponse = nresponse_per_field * nfields

            if (proc0) then
               read (unit=mat_unit) ie_dump, nresponse_dump
               if (ie_dump /= ie .or. nresponse /= nresponse_dump) &
                  call mp_abort('mismatch in ie/nresponse_dump')
            end if

            ! for each ky and set of connected kx values,
            ! must have a response matrix that is N x N
            ! with N = number of zeds per 2pi segment x number of 2pi segments
            if (.not. associated(response_matrix(iky)%eigen(ie)%zloc)) &
               allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse, nresponse))

            ! response_matrix%idx is needed to keep track of permutations
            ! to the response matrix made during LU decomposition
            ! it will be input to LU back substitution during linear solve
            if (.not. associated(response_matrix(iky)%eigen(ie)%idx)) &
               allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))
            if (proc0) then
               read (unit=mat_unit) response_matrix(iky)%eigen(ie)%idx
               read (unit=mat_unit) response_matrix(iky)%eigen(ie)%zloc
            end if

            call broadcast(response_matrix(iky)%eigen(ie)%idx)
            call broadcast(response_matrix(iky)%eigen(ie)%zloc)

         end do
      end do

      if (proc0) close (mat_unit)

      if (debug) then
         print *, 'File', file_name, ' successfully read by root proc for job: ', job_str
      end if
   end subroutine read_response_matrix

   subroutine get_dpdf_dphi_matrix_column(iky, ie, idx, nz_ext, nresponse, phi_ext, apar_ext, bpar_ext, pdf_ext)

      use stella_layouts, only: vmu_lo
      use parameters_numerical, only: time_upwind_plus
      use physics_parameters, only: include_apar, include_bpar
      use implicit_solve, only: get_gke_rhs, sweep_g_zext
      use arrays_fields, only: response_matrix
      use extended_zgrid, only: periodic, phase_shift
      use physics_parameters, only: full_flux_surface
#ifdef ISO_C_BINDING
      use mp, only: sgproc0
#endif

      implicit none

      integer, intent(in) :: iky, ie, idx, nz_ext, nresponse
      complex, dimension(:), intent(out) :: phi_ext, apar_ext, bpar_ext
      complex, dimension(:, vmu_lo%llim_proc:), intent(out) :: pdf_ext

      complex, dimension(:), allocatable :: dum
      integer :: ivmu, it
      integer :: offset_apar, offset_bpar

      ! provide a unit impulse to phi^{n+1} (or Delta phi^{n+1}) at the location
      ! in the extended zed domain corresponding to index 'idx'
      ! note that it is sufficient to give a unit real impulse (as opposed to
      ! separately giving real and imaginary impulse) for the following reason:
      ! split homogeneous GKE, L[f] = R[phi], into L[f1] = R[phir] and L[f2] = i*R[phii],
      ! with f = f1 + f2; then phi = df1/dphir * phir + df2/dphii * phii.
      ! however, we see that if phir = phii = 1, L[f1] = R[1] = L[-i*f2],
      ! and thus f2 = i * f1.  This gives phi = df1/dphir * (phir + i * phii) = df1/dphir * phi
      phi_ext = 0.0
      ! how phi^{n+1} enters the GKE depends on whether we are solving for the
      ! non-Boltzmann pdf, h, or the guiding centre pdf, 'g'
      phi_ext(idx) = time_upwind_plus
      
      !> TOGO-GA: check division rather than multiplication -- kept division for now to be consistent with 
      !> parallel_streaming phase shift 
      if (periodic(iky) .and. idx == 1) phi_ext(nz_ext) = phi_ext(1) / phase_shift(iky)

      ! dum is a scratch array that takes the place of the pdf and phi
      ! at the previous time level,
      ! which is set to zero for the response matrix approach
      allocate (dum(nz_ext)); dum = 0.0

      ! set the flux tube index to one
      ! need to check, but think this is okay as the homogeneous equation solved here for the
      ! response matrix construction is the same for all flux tubes in the flux tube train
      it = 1

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc         
         ! calculate the RHS of the GK equation (using dum=0 as the pdf at the previous time level,
         ! and phi_ext as the potential) and store it in pdf_ext
         call get_gke_rhs(ivmu, iky, ie, dum, phi_ext, dum, dum, dum, dum, pdf_ext(:, ivmu))
         ! given the RHS of the GK equation (pdf_ext), solve for the pdf at the
         ! new time level by sweeping in zed on the extended domain;
         ! the rhs is input as 'pdf_ext' and over-written with the updated solution for the pdf
         call sweep_g_zext(iky, ie, it, ivmu, pdf_ext(:, ivmu))
      end do

      deallocate (dum)

      ! we now have the pdf on the extended zed domain at this ky and set of connected kx values
      ! corresponding to a unit impulse in phi at this location
      ! now integrate over velocities to get a square response matrix
      ! (this ends the parallelization over velocity space, so every core should have a
      ! copy of phi_ext)
      call integrate_over_velocity(pdf_ext, phi_ext, apar_ext, bpar_ext, iky, ie)

#ifdef ISO_C_BINDING
      if (sgproc0) then
#endif
         response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx) = phi_ext(:nresponse)
         offset_apar = 0
         if (include_apar) then
            offset_apar = nresponse
            response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx) = apar_ext(:nresponse)
         end if
         if (include_bpar) then
            offset_bpar = offset_apar + nresponse
            response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx) = bpar_ext(:nresponse)
         end if
#ifdef ISO_C_BINDING
      end if
#endif

   end subroutine get_dpdf_dphi_matrix_column

   subroutine get_dpdf_dapar_matrix_column(iky, ie, idx, nz_ext, nresponse, phi_ext, apar_ext, bpar_ext, pdf_ext)

      use stella_layouts, only: vmu_lo
      use parameters_numerical, only: time_upwind_plus
      use physics_parameters, only: include_apar, include_bpar
      use implicit_solve, only: get_gke_rhs, sweep_g_zext
      use arrays_fields, only: response_matrix
      use extended_zgrid, only: periodic
#ifdef ISO_C_BINDING
      use mp, only: sgproc0
#endif

      implicit none

      integer, intent(in) :: iky, ie, idx, nz_ext, nresponse
      complex, dimension(:), intent(out) :: phi_ext, apar_ext, bpar_ext
      complex, dimension(:, vmu_lo%llim_proc:), intent(out) :: pdf_ext

      complex, dimension(:), allocatable :: dum
      integer :: ivmu, it
      integer :: offset_apar, offset_bpar
      
      ! provide a unit impulse to apar^{n+1} (or Delta apar^{n+1}) at the location
      ! in the extended zed domain corresponding to index 'idx'
      ! note that it is sufficient to give a unit real impulse (as opposed to
      ! separately giving real and imaginary impulse) for the following reason:
      ! split homogeneous GKE, L[f] = R[apar], into L[f1] = R[aparr] and L[f2] = i*R[apari],
      ! with f = f1 + f2; then apar = df1/daparr * aparr + df2/dapari * apari.
      ! however, we see that if aparr = apari = 1, L[f1] = R[1] = L[-i*f2],
      ! and thus f2 = i * f1.  This gives apar = df1/daparr * (aparr + i * apari) = df1/daparr * apar
      apar_ext = 0.0
      apar_ext(idx) = 1.0

      if (periodic(iky) .and. idx == 1) apar_ext(nz_ext) = apar_ext(1)

      ! dum is a scratch array that takes the place of the pdf and phi
      ! at the previous time level,
      ! which is set to zero for the response matrix approach
      allocate (dum(nz_ext)); dum = 0.0

      ! set the flux tube index to one
      ! need to check, but think this is okay as the homogeneous equation solved here for the
      ! response matrix construction is the same for all flux tubes in the flux tube train
      it = 1
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! calculate the RHS of the GK equation (using dum=0 as the pdf at the previous time level,
         ! and phi_ext as the potential) and store it in pdf_ext
         call get_gke_rhs(ivmu, iky, ie, dum, dum, apar_ext * time_upwind_plus, apar_ext, dum, dum, pdf_ext(:, ivmu))
         ! given the RHS of the GK equation (pdf_ext), solve for the pdf at the
         ! new time level by sweeping in zed on the extended domain;
         ! the rhs is input as 'pdf_ext' and over-written with the updated solution for the pdf
         call sweep_g_zext(iky, ie, it, ivmu, pdf_ext(:, ivmu))
      end do

      deallocate (dum)

      ! we now have the pdf on the extended zed domain at this ky and set of connected kx values
      ! corresponding to a unit impulse in phi at this location
      ! now integrate over velocities to get a square response matrix
      ! (this ends the parallelization over velocity space, so every core should have a
      ! copy of phi_ext)
      call integrate_over_velocity(pdf_ext, phi_ext, apar_ext, bpar_ext, iky, ie)

#ifdef ISO_C_BINDING
      if (sgproc0) then
#endif
         if (include_apar) then
            offset_apar = nresponse
         else
            offset_apar = 0
         end if
         response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx + nresponse) = phi_ext(:nresponse)
         if (include_apar) then
            response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx + offset_apar) = apar_ext(:nresponse)
         end if
         if (include_bpar) then
            offset_bpar = offset_apar + nresponse
            response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx + offset_apar) = bpar_ext(:nresponse)
         end if
#ifdef ISO_C_BINDING
      end if
#endif

   end subroutine get_dpdf_dapar_matrix_column

   ! modelled on get_dpdf_dphi_matrix_column above
   subroutine get_dpdf_dbpar_matrix_column(iky, ie, idx, nz_ext, nresponse, phi_ext, apar_ext, bpar_ext, pdf_ext)

      use stella_layouts, only: vmu_lo
      use parameters_numerical, only: time_upwind_plus
      use physics_parameters, only: include_apar, include_bpar
      use implicit_solve, only: get_gke_rhs, sweep_g_zext
      use arrays_fields, only: response_matrix
      use extended_zgrid, only: periodic
#ifdef ISO_C_BINDING
      use mp, only: sgproc0
#endif

      implicit none

      integer, intent(in) :: iky, ie, idx, nz_ext, nresponse
      complex, dimension(:), intent(out) :: phi_ext, apar_ext, bpar_ext
      complex, dimension(:, vmu_lo%llim_proc:), intent(out) :: pdf_ext

      complex, dimension(:), allocatable :: dum
      integer :: ivmu, it
      integer :: offset_apar, offset_bpar

      ! provide a unit impulse to phi^{n+1} (or Delta phi^{n+1}) at the location
      ! in the extended zed domain corresponding to index 'idx'
      ! note that it is sufficient to give a unit real impulse (as opposed to
      ! separately giving real and imaginary impulse) for the following reason:
      ! split homogeneous GKE, L[f] = R[phi], into L[f1] = R[phir] and L[f2] = i*R[phii],
      ! with f = f1 + f2; then phi = df1/dphir * phir + df2/dphii * phii.
      ! however, we see that if phir = phii = 1, L[f1] = R[1] = L[-i*f2],
      ! and thus f2 = i * f1.  This gives phi = df1/dphir * (phir + i * phii) = df1/dphir * phi
      bpar_ext = 0.0
      ! how phi^{n+1} enters the GKE depends on whether we are solving for the
      ! non-Boltzmann pdf, h, or the guiding centre pdf, 'g'
      bpar_ext(idx) = time_upwind_plus

      if (periodic(iky) .and. idx == 1) bpar_ext(nz_ext) = bpar_ext(1)

      ! dum is a scratch array that takes the place of the pdf and phi
      ! at the previous time level,
      ! which is set to zero for the response matrix approach
      allocate (dum(nz_ext)); dum = 0.0

      ! set the flux tube index to one
      ! need to check, but think this is okay as the homogeneous equation solved here for the
      ! response matrix construction is the same for all flux tubes in the flux tube train
      it = 1
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! calculate the RHS of the GK equation (using dum=0 as the pdf at the previous time level,
         ! and phi_ext as the potential) and store it in pdf_ext
         call get_gke_rhs(ivmu, iky, ie, dum, dum, dum, dum, dum, bpar_ext, pdf_ext(:, ivmu))
         ! given the RHS of the GK equation (pdf_ext), solve for the pdf at the
         ! new time level by sweeping in zed on the extended domain;
         ! the rhs is input as 'pdf_ext' and over-written with the updated solution for the pdf
         call sweep_g_zext(iky, ie, it, ivmu, pdf_ext(:, ivmu))
      end do

      deallocate (dum)

      ! we now have the pdf on the extended zed domain at this ky and set of connected kx values
      ! corresponding to a unit impulse in phi at this location
      ! now integrate over velocities to get a square response matrix
      ! (this ends the parallelization over velocity space, so every core should have a
      ! copy of phi_ext)
      call integrate_over_velocity(pdf_ext, phi_ext, apar_ext, bpar_ext, iky, ie)

#ifdef ISO_C_BINDING
      if (sgproc0) then
#endif
         offset_apar = 0
         if (include_apar) offset_apar = nresponse
         if (include_bpar) offset_bpar = offset_apar + nresponse
         
         response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx + offset_bpar) = phi_ext(:nresponse)
         if (include_apar) then
            response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx + offset_bpar) = apar_ext(:nresponse)
         end if
         if (include_bpar) then
            response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx + offset_bpar) = bpar_ext(:nresponse)
         end if
#ifdef ISO_C_BINDING
      end if
#endif

   end subroutine get_dpdf_dbpar_matrix_column

   subroutine integrate_over_velocity(g, phi, apar, bpar, iky, ie)

      use stella_layouts, only: vmu_lo
      use physics_parameters, only: include_apar, include_bpar

      implicit none

      complex, dimension(:, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:), intent(out) :: phi, apar, bpar
      integer, intent(in) :: iky, ie

      call integrate_over_velocity_phi(g, phi, iky, ie)
      if (include_apar) call integrate_over_velocity_apar(g, apar, iky, ie)
      if (include_bpar) call integrate_over_velocity_bpar(g, bpar, iky, ie)

   end subroutine integrate_over_velocity

   subroutine integrate_over_velocity_phi(g, phi, iky, ie)

      use stella_layouts, only: vmu_lo
      use species, only: nspec, spec
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: nsegments
      use vpamu_grids, only: integrate_species
      use gyro_averages, only: gyro_average
      use mp, only: sum_allreduce

      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use parameters_numerical, only: driftkinetic_implicit
      use vpamu_grids, only: integrate_species_ffs_rm

      use physics_parameters, only: full_flux_surface

      use gyro_averages, only: j0_B_const

      implicit none

      complex, dimension(:, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:), intent(out) :: phi
      integer, intent(in) :: iky, ie

      integer :: idx, iseg, ikx, iz, ia
      integer :: izl_offset
      real, dimension(nspec) :: wgt
      complex, dimension(:), allocatable :: g0

      integer :: ivmu, imu, iv, is

      ia = 1

      allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      wgt = spec%z * spec%dens_psi0
      phi = 0.

      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      do iz = iz_low(iseg), iz_up(iseg)
         idx = idx + 1
         if (.not. full_flux_surface .and. (.not. driftkinetic_implicit)) then
            call gyro_average(g(idx, :), iky, ikx, iz, g0)
            call integrate_species(g0, iz, wgt, phi(idx), reduce_in=.false.)
         else
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               is = is_idx(vmu_lo, ivmu)
               g0(ivmu) = g(idx, ivmu) * j0_B_const(iky, ikx, iz, ivmu)
            end do
            call integrate_species_ffs_rm(g0, wgt, phi(idx), reduce_in=.false.)
         end if
      end do

      izl_offset = 1
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               if (.not. full_flux_surface .and. (.not. driftkinetic_implicit)) then
                  call gyro_average(g(idx, :), iky, ikx, iz, g0)
                  call integrate_species(g0, iz, wgt, phi(idx), reduce_in=.false.)
               else
                  do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                     iv = iv_idx(vmu_lo, ivmu)
                     imu = imu_idx(vmu_lo, ivmu)
                     is = is_idx(vmu_lo, ivmu)
                     g0(ivmu) = g(idx, ivmu) * j0_B_const(iky, ikx, iz, ivmu)
                  end do
                  call integrate_species_ffs_rm(g0, wgt, phi(idx), reduce_in=.false.)
               end if
            end do
            if (izl_offset == 0) izl_offset = 1
         end do
      end if

      call sum_allreduce(phi)

   end subroutine integrate_over_velocity_phi

   subroutine integrate_over_velocity_bpar(g, bpar, iky, ie)

      use stella_layouts, only: vmu_lo, imu_idx
      use species, only: nspec, spec
      use physics_parameters, only: beta
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: nsegments
      use vpamu_grids, only: integrate_species, mu
      use gyro_averages, only: gyro_average_j1
      use mp, only: sum_allreduce

      implicit none

      complex, dimension(:, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:), intent(out) :: bpar
      integer, intent(in) :: iky, ie

      integer :: idx, iseg, ikx, iz, ia, imu, ivmu
      integer :: izl_offset
      real, dimension(nspec) :: wgt
      complex, dimension(:), allocatable :: g0

      ia = 1

      allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      wgt = -2.0 * beta * spec%temp_psi0 * spec%dens_psi0
      bpar = 0.

      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      do iz = iz_low(iseg), iz_up(iseg)
         idx = idx + 1
         call gyro_average_j1(g(idx, :), iky, ikx, iz, g0)
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            imu = imu_idx(vmu_lo, ivmu)
            g0(ivmu) = g0(ivmu) * mu(imu)
         end do
         call integrate_species(g0, iz, wgt, bpar(idx), reduce_in=.false.)
      end do
      izl_offset = 1
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               call gyro_average_j1(g(idx, :), iky, ikx, iz, g0)
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                  imu = imu_idx(vmu_lo, ivmu)
                  g0(ivmu) = g0(ivmu) * mu(imu)
               end do               
               call integrate_species(g0, iz, wgt, bpar(idx), reduce_in=.false.)
            end do
            if (izl_offset == 0) izl_offset = 1
         end do
      end if

      call sum_allreduce(bpar)

   end subroutine integrate_over_velocity_bpar

   subroutine integrate_over_velocity_apar(g, apar, iky, ie)

      use stella_layouts, only: vmu_lo, iv_idx
      use physics_parameters, only: beta
      use species, only: nspec, spec
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: nsegments
      use vpamu_grids, only: integrate_species
      use vpamu_grids, only: vpa
      use gyro_averages, only: gyro_average
      use mp, only: sum_allreduce

      implicit none

      complex, dimension(:, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:), intent(out) :: apar
      integer, intent(in) :: iky, ie

      integer :: idx, iseg, ikx, iz, ia
      integer :: ivmu, iv
      integer :: izl_offset
      real, dimension(nspec) :: wgt
      complex, dimension(:), allocatable :: g0

      ia = 1

      allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      wgt = spec%z * spec%dens_psi0 * spec%stm_psi0 * beta
      apar = 0.

      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      do iz = iz_low(iseg), iz_up(iseg)
         idx = idx + 1
         call gyro_average(g(idx, :), iky, ikx, iz, g0)
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            g0(ivmu) = g0(ivmu) * vpa(iv)
         end do
         call integrate_species(g0, iz, wgt, apar(idx), reduce_in=.false.)
      end do
      izl_offset = 1
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               call gyro_average(g(idx, :), iky, ikx, iz, g0)
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                  iv = iv_idx(vmu_lo, ivmu)
                  g0(ivmu) = g0(ivmu) * vpa(iv)
               end do
               call integrate_species(g0, iz, wgt, apar(idx), reduce_in=.false.)
            end do
            if (izl_offset == 0) izl_offset = 1
         end do
      end if

      call sum_allreduce(apar)

   end subroutine integrate_over_velocity_apar

   subroutine get_fields_for_response_matrix(phi, apar, bpar, iky, ie, dist)

      use physics_parameters, only: include_apar, include_bpar

      implicit none

      complex, dimension(:), intent(in out) :: phi, apar, bpar
      integer, intent(in) :: iky, ie
      character(*), intent(in) :: dist

      if (include_bpar) then
         call get_phi_and_bpar_for_response_matrix(phi, bpar, iky, ie, dist)
      else
         call get_phi_for_response_matrix(phi, iky, ie, dist)
      end if
      if (include_apar) call get_apar_for_response_matrix(apar, iky, ie, dist)
 
   end subroutine get_fields_for_response_matrix

   subroutine get_phi_for_response_matrix(phi, iky, ie, dist)

      use zgrid, only: nzgrid
      use species, only: spec
      use species, only: has_electron_species
      use geometry, only: dl_over_b
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: nsegments
      use grids_kxky, only: zonal_mode, akx
      use arrays_fields, only: gamtot, gamtot3
      use arrays_fields, only: gamtot_h, gamtot3_h
      use physics_parameters, only: adiabatic_option_switch
      use physics_parameters, only: adiabatic_option_fieldlineavg

      implicit none

      complex, dimension(:), intent(inout) :: phi
      integer, intent(in) :: iky, ie
      character(*), intent(in) :: dist

      integer :: idx, iseg, ikx, iz, ia
      integer :: izl_offset
      complex :: tmp
      real, dimension(:), allocatable :: gamma_fac

      ia = 1

      allocate (gamma_fac(-nzgrid:nzgrid))

      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      if (dist == 'h') then
         gamma_fac = gamtot_h
      else
         gamma_fac = gamtot(iky, ikx, :)
      end if
      if (zonal_mode(iky) .and. abs(akx(ikx)) < epsilon(0.)) then
         phi(:) = 0.0
         return
      end if
      do iz = iz_low(iseg), iz_up(iseg)
         idx = idx + 1
         phi(idx) = phi(idx) / gamma_fac(iz)
      end do
      izl_offset = 1
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            if (dist == 'h') then
               gamma_fac = gamtot_h
            else
               gamma_fac = gamtot(iky, ikx, :)
            end if
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               phi(idx) = phi(idx) / gamma_fac(iz)
            end do
            if (izl_offset == 0) izl_offset = 1
         end do
      end if

      if (.not. has_electron_species(spec) .and. &
          adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         if (zonal_mode(iky)) then
            ! no connections for ky = 0
            iseg = 1
            tmp = sum(dl_over_b(ia, :) * phi)
            if (dist == 'h') then
               phi = phi + tmp * gamtot3_h
            else
               phi = phi + tmp * gamtot3(ikxmod(1, ie, iky), :)
            end if
         end if
      end if

      deallocate (gamma_fac)

   end subroutine get_phi_for_response_matrix

   subroutine get_phi_and_bpar_for_response_matrix(phi, bpar, iky, ie, dist)

      use zgrid, only: nzgrid
      use species, only: spec
      use species, only: has_electron_species
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: nsegments
      use grids_kxky, only: zonal_mode, akx
      use arrays_fields, only: gamtotinv11, gamtotinv13, gamtotinv31, gamtotinv33
      use arrays_fields, only: gamtot_h
      use physics_parameters, only: adiabatic_option_switch
      use physics_parameters, only: adiabatic_option_fieldlineavg
      use mp, only: mp_abort
      
      implicit none

      complex, dimension(:), intent(inout) :: phi, bpar
      integer, intent(in) :: iky, ie
      character(*), intent(in) :: dist

      integer :: idx, iseg, ikx, iz, ia
      integer :: izl_offset
      complex :: antot1, antot3
      real, dimension(:), allocatable :: gammainv11, gammainv13, gammainv31, gammainv33

      ia = 1

      allocate (gammainv11(-nzgrid:nzgrid))
      allocate (gammainv13(-nzgrid:nzgrid))
      allocate (gammainv31(-nzgrid:nzgrid))
      allocate (gammainv33(-nzgrid:nzgrid))

      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      if (dist == 'h') then
         gammainv11 = 1.0/gamtot_h
         gammainv13 = 0.0
         gammainv31 = 0.0
         gammainv33 = 1.0
      else
         gammainv11 = gamtotinv11(iky, ikx, :)
         gammainv13 = gamtotinv13(iky, ikx, :)
         gammainv31 = gamtotinv31(iky, ikx, :)
         gammainv33 = gamtotinv33(iky, ikx, :)
      end if
      if (zonal_mode(iky) .and. abs(akx(ikx)) < epsilon(0.)) then
         phi(:) = 0.0
         bpar(:) = 0.0
         return
      end if
      do iz = iz_low(iseg), iz_up(iseg)
         idx = idx + 1
         antot1 = phi(idx)
         antot3 = bpar(idx)
         phi(idx) = antot1 * gammainv11(iz) + antot3 * gammainv13(iz)
         bpar(idx) = antot1 * gammainv31(iz) + antot3 * gammainv33(iz)
      end do
      izl_offset = 1
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            if (dist == 'h') then
               gammainv11 = 1.0/gamtot_h
               gammainv13 = 0.0
               gammainv31 = 0.0
               gammainv33 = 1.0
            else
               gammainv11 = gamtotinv11(iky, ikx, :)
               gammainv13 = gamtotinv13(iky, ikx, :)
               gammainv31 = gamtotinv31(iky, ikx, :)
               gammainv33 = gamtotinv33(iky, ikx, :)
            end if
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               antot1 = phi(idx)
               antot3 = bpar(idx)
               phi(idx) = antot1 * gammainv11(iz) + antot3 * gammainv13(iz)
               bpar(idx) = antot1 * gammainv31(iz) + antot3 * gammainv33(iz)
            end do
            if (izl_offset == 0) izl_offset = 1
         end do
      end if

      if (.not. has_electron_species(spec) .and. &
          adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         call mp_abort('adiabatic electrons not yet supported for include_bpar = T. aborting.')
      end if

      deallocate (gammainv11, gammainv13, gammainv31, gammainv33)

   end subroutine get_phi_and_bpar_for_response_matrix

   subroutine get_apar_for_response_matrix(apar, iky, ie, dist)

      use zgrid, only: nzgrid
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: nsegments
      use grids_kxky, only: zonal_mode, akx
      use arrays_fields, only: apar_denom
      use arrays_dist_fn, only: kperp2

      implicit none

      complex, dimension(:), intent(in out) :: apar
      integer, intent(in) :: iky, ie
      character(*), intent(in) :: dist

      integer :: idx, iseg, ikx, iz, ia
      integer :: izl_offset
      real, dimension(:), allocatable :: denominator

      ia = 1

      allocate (denominator(-nzgrid:nzgrid))

      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      if (dist == 'g') then
         denominator = kperp2(iky, ikx, ia, :)
      else if (dist == 'gbar') then
         denominator = apar_denom(iky, ikx, :)
      end if
      if (zonal_mode(iky) .and. abs(akx(ikx)) < epsilon(0.)) then
         apar(:) = 0.0
         return
      end if
      do iz = iz_low(iseg), iz_up(iseg)
         idx = idx + 1
         apar(idx) = apar(idx) / denominator(iz)
      end do
      izl_offset = 1
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            if (dist == 'g') then
               denominator = kperp2(iky, ikx, ia, :)
            elseif (dist == 'gbar') then
               denominator = apar_denom(iky, ikx, :)
            end if
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               apar(idx) = apar(idx) / denominator(iz)
            end do
            if (izl_offset == 0) izl_offset = 1
         end do
      end if

      deallocate (denominator)

   end subroutine get_apar_for_response_matrix

   subroutine finish_response_matrix

      use arrays_fields, only: response_matrix
#if !defined ISO_C_BINDING

      implicit none

#else
      use arrays_fields, only: response_window
      use mpi

      implicit none

      integer :: ierr

      if (response_window /= MPI_WIN_NULL) call mpi_win_free(response_window, ierr)
#endif

      if (allocated(response_matrix)) deallocate (response_matrix)
      response_matrix_initialized = .false.

   end subroutine finish_response_matrix

!----------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
! !!!!!!!!!!!!!!PARALLEL LU DECOMPOSITIONS!!!!!!!!!!!!!!!! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
!----------------------------------------------------------!

#ifdef ISO_C_BINDING

   !this subroutine parallelizes the LU decomposition on a single
   !node using MPIs shared memory interface
   !It also splits up jtwist the independent matrices across nodes
   !Ideal speed up: cores_per_node*min(jtwist,ncores)
   subroutine parallel_LU_decomposition_local(iky)

      use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
      use arrays_fields, only: response_matrix
      use mp, only: barrier, broadcast, sum_allreduce
      use mp, only: mp_comm, scope, allprocs, sharedprocs, curr_focus
      use mp, only: scrossdomprocs, sgproc0, mp_abort, real_size
      use mp, only: job, iproc, proc0, nproc, numnodes, inode
      use mp_lu_decomposition, only: lu_decomposition_local
      use job_manage, only: njobs
      use extended_zgrid, only: neigen
      use mpi

      implicit none

      integer, intent(in) :: iky

      integer, dimension(:, :), allocatable :: eig_limits
      integer, dimension(:), allocatable :: job_list
      logical, dimension(:, :), allocatable :: node_jobs

      complex, dimension(:, :), pointer :: lu

      type(c_ptr) :: bptr

      logical :: needs_send = .false.

      integer :: prior_focus, nodes_on_job
      integer :: ijob, j, ie, n, ediv, emod
      integer :: jroot, neig, ierr, win, nroot
      integer(kind=MPI_ADDRESS_KIND) :: win_size
      integer :: disp_unit = 1
      real :: dmax

      prior_focus = curr_focus

      call scope(sharedprocs)

      allocate (node_jobs(0:(numnodes - 1), 0:(njobs - 1))); node_jobs = .false.
      allocate (job_list(0:(nproc - 1))); job_list = 0
      allocate (eig_limits(0:numnodes, 0:(njobs - 1))); eig_limits = 0

      job_list(iproc) = job
      call sum_allreduce(job_list)

      if (proc0) then
         do j = 0, nproc - 1
            node_jobs(inode, job_list(j)) = .true. !create a map of which nodes have which jobs
         end do
      end if

      !make sure all processors have this map
      call scope(allprocs)
      call mpi_allreduce &
         (MPI_IN_PLACE, node_jobs, size(node_jobs), MPI_LOGICAL, MPI_LOR, mp_comm, ierr)
      call scope(sharedprocs)

      do ijob = 0, njobs - 1
         jroot = -1
         do j = 0, nproc - 1
            if (job_list(j) == ijob) then
               jroot = j !the first processor on this job will be the root process
               exit
            end if
         end do

         if (jroot == -1) cycle !no processors on this node are on this job

         if (iproc == jroot) neig = neigen(iky)

         ! broadcast number of matrices
         call broadcast(neig, jroot)

         ! split up neig across nodes that have the current job
         nodes_on_job = count(node_jobs(:, ijob))
         ediv = neig / nodes_on_job
         emod = mod(neig, nodes_on_job)

         eig_limits(0, ijob) = 1
         do j = 1, numnodes
            if (node_jobs(j - 1, ijob)) then
               eig_limits(j, ijob) = eig_limits(j - 1, ijob) + ediv
               if (emod > 0) then
                  eig_limits(j, ijob) = eig_limits(j, ijob) + 1
                  emod = emod - 1
               end if
            else
               eig_limits(j, ijob) = eig_limits(j - 1, ijob)
            end if
         end do

         do ie = eig_limits(inode, ijob), eig_limits(inode + 1, ijob) - 1
            win_size = 0
            if (iproc == jroot) then
               needs_send = .true.
               n = size(response_matrix(iky)%eigen(ie)%idx)
               win_size = int(n * n, MPI_ADDRESS_KIND) * 2 * real_size !complex size
            end if

            !broadcast size of matrix
            call broadcast(n, jroot)

            !allocate the window
            call mpi_win_allocate_shared(win_size, disp_unit, MPI_INFO_NULL, mp_comm, bptr, win, ierr)

            if (iproc /= jroot) then
               !make sure all the procs have the right memory address
               call mpi_win_shared_query(win, jroot, win_size, disp_unit, bptr, ierr)
            end if

            ! bind this c_ptr to our fortran matrix
            call c_f_pointer(bptr, lu, (/n, n/))

            !load the matrix
            if (iproc == jroot) lu = response_matrix(iky)%eigen(ie)%zloc

            !syncronize the processors
            call mpi_win_fence(0, win, ierr)

            ! All the processors have the matrix.
            ! Now perform LU decomposition
            call lu_decomposition_local(mp_comm, jroot, win, lu, &
                                        response_matrix(iky)%eigen(ie)%idx, dmax)

            !copy the decomposed matrix over
            if (iproc == jroot) response_matrix(iky)%eigen(ie)%zloc = lu

            call mpi_win_free(win, ierr)
         end do
      end do

      call scope(scrossdomprocs)

      !copy all the matrices across all nodes
      if (sgproc0) then
         do ie = 1, neigen(iky)
            nroot = 0
            if (needs_send .and. &
                (ie >= eig_limits(inode, job) .and. ie < eig_limits(inode + 1, job))) nroot = iproc
            !first let processors know who is sending the data
            call sum_allreduce(nroot)
            !now send the data
            call broadcast(response_matrix(iky)%eigen(ie)%zloc, nroot)
            call broadcast(response_matrix(iky)%eigen(ie)%idx, nroot)
         end do
      end if

      call scope(prior_focus)

      deallocate (node_jobs, job_list, eig_limits)
   end subroutine parallel_LU_decomposition_local

#endif /* ISO_C_BINDING */

   !this subroutine parallelizes the LU decomposition across
   !all cores. Ideal speed up: ncores
   subroutine parallel_LU_decomposition_global(iky)

      use arrays_fields, only: response_matrix
      use mp, only: barrier, broadcast, sum_allreduce
      use mp, only: mp_comm, scope, allprocs, sharedprocs, curr_focus
      use mp, only: job, iproc, proc0, nproc, mpicmplx
#ifdef ISO_C_BINDING
      use mp, only: sgproc0, scrossdomprocs
#endif
      use job_manage, only: njobs
      use extended_zgrid, only: neigen
      use mpi
      use linear_solve, only: imaxloc

      implicit none

      integer, intent(in) :: iky

      integer, dimension(:), allocatable :: job_roots, eig_roots
      integer, dimension(:), allocatable :: row_limits, eig_limits
      integer, dimension(MPI_STATUS_SIZE) :: status

      real, parameter :: zero = 1.0e-20
      integer, dimension(:), allocatable :: idx
      complex, dimension(:, :), allocatable :: lu
      real, dimension(:), allocatable :: vv
      complex, dimension(:), allocatable :: dum

      integer :: sproc
      logical :: sproc0

      integer :: eig_comm, ceig_comm !c for 'cross'
      integer :: ieig_core, ceig_core, eig_cores
      integer :: ncomm
      integer :: prior_focus
      integer :: ie, ie_hi, r_lo, r_hi
      integer :: ijob, i, j, k, n, n_send, rsize
      integer :: imax, neig, ierr
      integer :: istage, nstage
      integer :: rdiv, rmod
      integer :: ediv, emod

      real :: dmax, tmp

      prior_focus = curr_focus

      sproc = iproc
      sproc0 = proc0

      call scope(allprocs)

      allocate (job_roots(0:njobs - 1)); job_roots = 0

      if (sproc0) job_roots(job) = iproc

      call sum_allreduce(job_roots)

      do ijob = 0, njobs - 1
         if (job == ijob .and. sproc0) then
            neig = neigen(iky)
         end if

         ! broadcast number of matrices for this job
         call broadcast(neig, job_roots(ijob))

         !set up communicator for cores working on a single matrix
         call mpi_comm_split(mp_comm, mod(iproc, neig), iproc, eig_comm, ierr)
         call mpi_comm_size(eig_comm, eig_cores, ierr)
         call mpi_comm_rank(eig_comm, ieig_core, ierr)

         !set up a communicator that crosses the previous one
         call mpi_comm_split(mp_comm, ieig_core, iproc, ceig_comm, ierr)
         call mpi_comm_rank(ceig_comm, ceig_core, ierr)

         call mpi_bcast(ceig_core, 1, MPI_INT, 0, eig_comm, ierr)

         ncomm = min(neig, nproc) !number of communicators

         allocate (eig_roots(0:ncomm - 1)); eig_roots = 0
         allocate (eig_limits(0:ncomm))
         allocate (row_limits(0:eig_cores))

         if (ieig_core == 0) eig_roots(ceig_core) = iproc

         call sum_allreduce(eig_roots)

         ! split up neigen across cores
         ediv = neig / ncomm
         emod = mod(neig, ncomm)

         !how many stages will the LU decomposition take?
         nstage = ediv
         if (emod > 0) nstage = nstage + 1

         !determine which parts of neigen this communicator processes
         eig_limits(0) = 1
         do j = 1, ncomm
            eig_limits(j) = eig_limits(j - 1) + ediv
            if (j <= emod) then
               eig_limits(j) = eig_limits(j) + 1
            end if
         end do

         do istage = 0, nstage - 1
            !transfer the data from job root to root of subcommunicator
            do j = 0, ncomm - 1
               ie = eig_limits(j) + istage
               ie_hi = eig_limits(j + 1) - 1
               if (ie > ie_hi) cycle

               if (iproc == job_roots(ijob) .and. iproc == eig_roots(j)) then !no need for data transfer
                  n = size(response_matrix(iky)%eigen(ie)%idx)
                  allocate (lu(n, n))
                  lu = response_matrix(iky)%eigen(ie)%zloc
               else if (iproc == job_roots(ijob)) then !send data to subroots
                  !send size of matrix
                  n_send = size(response_matrix(iky)%eigen(ie)%idx)
                  call mpi_send(n_send, 1, MPI_INT, eig_roots(j), j, mp_comm, ierr)
                  !send matrix
                  call mpi_send(response_matrix(iky)%eigen(ie)%zloc, &
                                n_send * n_send, mpicmplx, eig_roots(j), nproc + j, mp_comm, ierr)
               else if (iproc == eig_roots(j)) then !subroot gets the data
                  !receive size of matrix
                  call mpi_recv(n, 1, MPI_INT, job_roots(ijob), j, mp_comm, status, ierr)
                  allocate (lu(n, n))
                  !receive matrix
                  call mpi_recv(lu, n * n, mpicmplx, job_roots(ijob), nproc + j, mp_comm, status, ierr)
               end if
            end do

            if (istage >= (eig_limits(ceig_core + 1) - eig_limits(ceig_core))) cycle !nothing for this communicator to do

            !broadcast matrix and its size across the communicator
            call mpi_bcast(n, 1, MPI_INT, 0, eig_comm, ierr)
            if (.not. allocated(lu)) allocate (lu(n, n))
            if (.not. allocated(vv)) allocate (vv(n))

            call mpi_bcast(lu, n * n, mpicmplx, 0, eig_comm, ierr)

            allocate (dum(n))
            allocate (idx(n))

            ! All the processors have the matrix.
            ! Now perform LU decomposition
            vv = maxval(cabs(lu), dim=2)
            if (any(vv == 0.0)) &
               write (*, *) 'singular matrix in lu_decomposition on job ', job, ', process ', iproc
            vv = 1.0 / vv
            do j = 1, n
               !divide up the work using row_limits
               rdiv = (n - j) / eig_cores
               rmod = mod(n - j, eig_cores)
               row_limits(0) = j + 1
               if (rdiv == 0) then
                  row_limits(rmod + 1:) = -1
                  do k = 1, rmod
                     row_limits(k) = row_limits(k - 1) + 1
                  end do
               else
                  do k = 1, eig_cores
                     row_limits(k) = row_limits(k - 1) + rdiv
                     if (k <= rmod) row_limits(k) = row_limits(k) + 1
                  end do
               end if

               !pivot if needed
               dmax = -1.0
               do k = j, n
                  tmp = vv(k) * abs(lu(k, j))
                  if (tmp > dmax) then
                     dmax = tmp
                     imax = k
                  end if
               end do
!         imax = (j-1) + imaxloc(vv(j:n)*cabs(lu(j:n,j)))
               if (j /= imax) then
                  dum = lu(imax, :)
                  lu(imax, :) = lu(j, :)
                  lu(j, :) = dum
                  vv(imax) = vv(j)
               end if
               if (ieig_core == 0) idx(j) = imax

               !get the lead multiplier
               if (lu(j, j) == 0.0) lu(j, j) = zero
               do i = j + 1, n
                  lu(i, j) = lu(i, j) / lu(j, j)
               end do

               r_lo = row_limits(ieig_core)
               r_hi = row_limits(ieig_core + 1) - 1

               do k = r_lo, r_hi
                  do i = j + 1, n
                     lu(i, k) = lu(i, k) - lu(i, j) * lu(j, k)
                  end do
               end do

               do i = 0, eig_cores - 1
                  r_lo = row_limits(i)
                  r_hi = row_limits(i + 1) - 1
                  rsize = (r_hi - r_lo + 1) * (n - j)
                  if (r_lo > r_hi) cycle
                  !call mpi_bcast(lu(j+1:n,r_lo:r_hi),rsize,mpicmplx,i,eig_comm,ierr)
                  do k = r_lo, r_hi
                     call mpi_bcast(lu(j + 1:n, k), n - j, mpicmplx, i, eig_comm, ierr)
                  end do
               end do
            end do
            !LU decomposition ends here

            !copy the decomposed matrix over
            do j = 0, ncomm - 1

               ie = eig_limits(j) + istage
               ie_hi = eig_limits(j + 1) - 1
               if (ie > ie_hi) cycle

               if (iproc == job_roots(ijob) .and. iproc == eig_roots(j)) then !no need for data transfer
                  response_matrix(iky)%eigen(ie)%zloc = lu
                  response_matrix(iky)%eigen(ie)%idx = idx
               else if (iproc == eig_roots(j)) then !subroot sends the data
                  !send indices
                  call mpi_send(idx, n, MPI_INT, job_roots(ijob), j, mp_comm, ierr)
                  !send matrix
                  call mpi_send(lu, n * n, mpicmplx, job_roots(ijob), nproc + j, mp_comm, ierr)
               else if (iproc == job_roots(ijob)) then !receive data from subroot
                  !receive indices
                  call mpi_recv(response_matrix(iky)%eigen(ie)%idx, &
                                n, MPI_INT, eig_roots(j), j, mp_comm, status, ierr)
                  !receive matrix
                  call mpi_recv(response_matrix(iky)%eigen(ie)%zloc, &
                                n * n, mpicmplx, eig_roots(j), nproc + j, mp_comm, status, ierr)
               end if
            end do
            deallocate (vv, lu, idx, dum)
         end do
         deallocate (eig_roots, eig_limits, row_limits)
      end do

#ifdef ISO_C_BINDING
      if (sgproc0) then
         call scope(scrossdomprocs)
         !copy all the matrices across all nodes
         do ie = 1, neigen(iky)
            call broadcast(response_matrix(iky)%eigen(ie)%zloc)
            call broadcast(response_matrix(iky)%eigen(ie)%idx)
         end do
      end if

      call scope(prior_focus)
#else
      call scope(prior_focus)

      !copy all the matrices across all nodes
      do ie = 1, neigen(iky)
         call broadcast(response_matrix(iky)%eigen(ie)%zloc)
         call broadcast(response_matrix(iky)%eigen(ie)%idx)
      end do
#endif
      deallocate (job_roots)
   end subroutine parallel_LU_decomposition_global

end module response_matrix
