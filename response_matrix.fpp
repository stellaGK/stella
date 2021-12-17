module response_matrix

   use netcdf
   use mpi

   implicit none

   public :: init_response_matrix, finish_response_matrix
   public :: read_response_matrix
   public :: response_matrix_initialized

   private

   logical :: response_matrix_initialized = .false.
   integer, parameter :: mat_unit = 70

#if defined MPI && defined ISO_C_BINDING
   integer :: window = MPI_WIN_NULL
#endif

contains

   subroutine init_response_matrix

      use linear_solve, only: lu_decomposition
      use fields_arrays, only: response_matrix
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, is_idx
      use kt_grids, only: naky
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: neigen, ikxmod
      use extended_zgrid, only: nsegments
      use extended_zgrid, only: nzed_segment
      use extended_zgrid, only: periodic
      use job_manage, only: time_message
      use mp, only: proc0, job, mp_abort
      use run_parameters, only: mat_gen, lu_option_switch
      use run_parameters, only: lu_option_none, lu_option_local, lu_option_global
      use run_parameters, only: fphi, fapar, fbpar
      use system_fortran, only: systemf
#if defined MPI && defined ISO_C_BINDING
      use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
      use mp, only: curr_focus, sgproc0, mp_comm, sharedsubprocs, scope, barrier
      use mp, only: real_size, nbytes_real, real_size
      use mpi
#endif

      implicit none

      integer :: iky, ie, iseg, iz
      integer :: ikx
      integer :: nz_ext, nresponse, nresponse_per_field
      integer :: idx
      integer :: izl_offset, izup
#if defined MPI && defined ISO_C_BINDING
      integer :: prior_focus, ierr
      integer :: disp_unit = 1
      integer(kind=MPI_ADDRESS_KIND) :: win_size
      integer*8 :: cur_pos
      type(c_ptr) :: bptr, cptr
#endif
      real :: dum
      complex, dimension(:), allocatable :: fields_ext
      complex, dimension(:, :), allocatable :: gext
      logical :: debug = .false.
      character(100) :: message_dgdphi, message_QN, message_lu
      real, dimension(2) :: time_response_matrix_dgdphi
      real, dimension(2) :: time_response_matrix_QN
      real, dimension(2) :: time_response_matrix_lu

      ! Related to the saving of the the matrices in netcdf format
      character(len=15) :: fmt, job_str
      character(len=100) :: file_name
      integer :: istatus
      istatus = 0

      if (proc0 .and. debug) then
         write (*, *) " "
         write (*, '(A)') "    ############################################################"
         write (*, '(A)') "                         RESPONSE MATRIX"
         write (*, '(A)') "    ############################################################"
      end if
      message_dgdphi = '     calculate dgdphi: '
      message_QN = '     calculate QN:     '
      message_lu = '     calculate LU:     '
      time_response_matrix_dgdphi = 0
      time_response_matrix_QN = 0
      time_response_matrix_lu = 0

!   All matrices handled by processor i_proc and job are stored
!   on a single file named: response_mat_job.iproc
      fmt = '(I5.5)'
      if (proc0 .and. mat_gen) then

         call systemf('mkdir -p mat')

         write (job_str, '(I1.1)') job
         file_name = './mat/response_mat_'//trim(job_str)
         open (unit=mat_unit, status='replace', file=file_name, &
               position='rewind', action='write', form='unformatted')
         write (unit=mat_unit) naky
      end if

      if (response_matrix_initialized) return
      response_matrix_initialized = .true.

      if (.not. allocated(response_matrix)) allocate (response_matrix(naky))

#if defined ISO_C_BINDING && defined MPI

!   Create a single shared memory window for all the response matrices and
!   permutation arrays.
!   Creating a window for each matrix/array would lead to performance
!   degradation on some clusters
      if (window == MPI_WIN_NULL) then
         prior_focus = curr_focus
         call scope(sharedsubprocs)
         win_size = 0
         if (sgproc0) then
            do iky = 1, naky
               do ie = 1, neigen(iky)
                  if (periodic(iky)) then
                     nresponse = nsegments(ie, iky) * nzed_segment
                  else
                     nresponse_per_field = nsegments(ie, iky) * nzed_segment + 1
                  end if

                  ! The response matrix must be at least nfields*nz_ext (or nfields*(nz_ext-1))
                  ! in size. For now, (for simplicity of coding), accept  3 fields; NB leads
                  ! to unnecessary memory cost for <3 fields.
                  ! Should we count the number of fields and use that?
                  nresponse = nresponse_per_field * 3

                  win_size = win_size &
                             + int(nresponse, MPI_ADDRESS_KIND) * 4_MPI_ADDRESS_KIND &
                             + int(nresponse**2, MPI_ADDRESS_KIND) * 2 * real_size
               end do
            end do
         end if
         call mpi_win_allocate_shared(win_size, disp_unit, MPI_INFO_NULL, mp_comm, &
                                      bptr, window, ierr)

         if (.not. sgproc0) then
            !make sure all the procs have the right memory address
            call mpi_win_shared_query(window, 0, win_size, disp_unit, bptr, ierr)
         end if
         call mpi_win_fence(0, window, ierr)

         !the following is a hack that allows us to perform pointer arithmetic in Fortran
         cur_pos = transfer(bptr, cur_pos)

         call scope(prior_focus)
      end if
#endif

      ! for a given ky and set of connected kx values
      ! give a unit impulse to phi at each zed location
      ! in the extended domain and solve for h(zed_extended,(vpa,mu,s))

      do iky = 1, naky

         if (proc0 .and. mat_gen) THEN
            write (unit=mat_unit) iky, neigen(iky)
         end if

         ! the response matrix for each ky has neigen(ky)
         ! independent sets of connected kx values
         if (.not. associated(response_matrix(iky)%eigen)) &
            allocate (response_matrix(iky)%eigen(neigen(iky)))

         if (proc0 .and. debug) then
            call time_message(.false., time_response_matrix_dgdphi, message_dgdphi)
         end if

         ! loop over the sets of connected kx values
         do ie = 1, neigen(iky)

            ! number of zeds x number of segments
            nz_ext = nsegments(ie, iky) * nzed_segment + 1

            ! treat zonal mode specially to avoid double counting
            ! as it is periodic
            if (periodic(iky)) then
               nresponse = nz_ext - 1
            else
               nresponse_per_field = nz_ext
            end if

            ! The response matrix must be at least nfields*nz_ext (or nfields*(nz_ext-1))
            ! in size. For now, (for simplicity of coding), accept  3 fields; NB leads
            ! to unnecessary memory cost for <3 fields.
            ! Should we count the fields and use that?
            nresponse = nresponse_per_field * 3

            if (proc0 .and. mat_gen) then
               write (unit=mat_unit) ie, nresponse
            end if

#if !defined ISO_C_BINDING || !defined MPI
            ! for each ky and set of connected kx values,
            ! must have a response matrix that is N x N
            ! with N = number of zeds per 2pi segment x number of 2pi segments x nfield
            if (.not. associated(response_matrix(iky)%eigen(ie)%zloc)) &
               allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse, nresponse))

            ! response_matrix%idx is needed to keep track of permutations
            ! to the response matrix made during LU decomposition
            ! it will be input to LU back substitution during linear solve
            if (.not. associated(response_matrix(iky)%eigen(ie)%idx)) &
               allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))
#else
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
#endif

            allocate (gext(nz_ext, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            allocate (fields_ext(3 * nz_ext)) ! Multiply by 3 because 3 fields
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
               if (fphi > epsilon(0.0)) then
                  ! Get the gext corresponding to a unit impulse in phi
                  call get_dgdfield_matrix_column(iky, ikx, iz, ie, idx, nz_ext, gext, field="phi")

                  ! we now have g on the extended zed domain at this ky and set of connected kx values
                  ! corresponding to a unit impulse in phi at this location
                  ! now calculate fields_ext=(phi(z+), apar(z+), bpar(z+)) corresponding
                  ! to this gbar, and store in the appropriate response_matrix column.
                  ! (this ends the parallelization over velocity space, so every core should have a
                  !  copy of fields_ext)
                  call get_em_fields_for_response_matrix(gext, iky, ie, nresponse_per_field, fields_ext)

                  ! We now have the fields - use to populate the column in response_matrix,
                  ! which is actually (identity matrix - response matrix)
                  ! NB fields_ext is always 3*nz_ext in size but the response matrix
                  ! is either (3*nz_ext x 3*nz_ext) (non-zonal modes) or
                  ! (3*(nz_ext-1) x 3*(nz_ext-1)) (zonal modes)

                  ! Subtract identity matrix
                  fields_ext(idx) = fields_ext(idx) - 1.0

                  ! negative sign because the matrix to be inverted in streaming equation
                  ! is actually (identity matrix - response matrix)
#if !defined ISO_C_BINDING || !defined MPI
                  response_matrix(iky)%eigen(ie)%zloc(:, idx) = -fields_ext(:nresponse)
#else
                  if (sgproc0) response_matrix(iky)%eigen(ie)%zloc(:, idx) = -fields_ext(:nresponse)
#endif
               end if
               if (fapar > epsilon(0.0)) then
                  ! Get the gext corresponding to a unit impulse in apar
                  call get_dgdfield_matrix_column(iky, ikx, iz, ie, idx, nz_ext, gext, field="apar")
                  ! Calculate ( phi(z+) , apar(z+) , bpar(z+) ) corresponding to
                  ! gext
                  call get_em_fields_for_response_matrix(gext, iky, ie, nresponse_per_field, fields_ext)

                  ! Subtract identity matrix - because this is populating the
                  ! (nresponse_per_field+idx)th column, identity matrix is 1 for (nresponse_per_field+idx),
                  ! zero otherwise
                  fields_ext(nresponse_per_field + idx) = fields_ext(nresponse_per_field + idx) - 1.0
#if !defined ISO_C_BINDING || !defined MPI
                  response_matrix(iky)%eigen(ie)%zloc(:, (nresponse_per_field + idx)) = -fields_ext(:nresponse)
#else
                  if (sgproc0) response_matrix(iky)%eigen(ie)%zloc(:, (nresponse_per_field + idx)) = -fields_ext(:nresponse)
#endif
               end if
               if (fbpar > epsilon(0.0)) then
                  ! Get the gext corresponding to a unit impulse in bpar
                  call get_dgdfield_matrix_column(iky, ikx, iz, ie, idx, nz_ext, gext, field="bpar")
                  ! Calculate ( phi(z+) , apar(z+) , bpar(z+) ) corresponding to
                  ! gext
                  call get_em_fields_for_response_matrix(gext, iky, ie, nresponse_per_field, fields_ext)

                  ! Subtract identity matrix - because this is populating the
                  ! (nresponse_per_field+idx)th column, identity matrix is 1 for (nresponse_per_field+idx),
                  ! zero otherwise
                  fields_ext(2 * nresponse_per_field + idx) = fields_ext(2 * nresponse_per_field + idx) - 1.0

#if !defined ISO_C_BINDING || !defined MPI
                  response_matrix(iky)%eigen(ie)%zloc(:, (2 * nresponse_per_field + idx)) = -fields_ext(:nresponse)
#else
                  if (sgproc0) response_matrix(iky)%eigen(ie)%zloc(:, (2 * nresponse_per_field + idx)) = -fields_ext(:nresponse)
#endif
               end if
            end do
            ! once we have used one segments, remaining segments
            ! have one fewer unique zed point

            izl_offset = 1
            if (nsegments(ie, iky) > 1) then
               do iseg = 2, nsegments(ie, iky)
                  ikx = ikxmod(iseg, ie, iky)
                  do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
                     idx = idx + 1

                     if (fphi > epsilon(0.0)) then
                        ! Get the gext and the fields corresponding to a unit impulse in phi
                        call get_dgdfield_matrix_column(iky, ikx, iz, ie, idx, nz_ext, gext, field="phi")
                        call get_em_fields_for_response_matrix(gext, iky, ie, nresponse_per_field, fields_ext)
                        fields_ext(idx) = fields_ext(idx) - 1.0
#if !defined ISO_C_BINDING || !defined MPI
                        response_matrix(iky)%eigen(ie)%zloc(:, idx) = -fields_ext(:nresponse)
#else
                        if (sgproc0) response_matrix(iky)%eigen(ie)%zloc(:, idx) = -fields_ext(:nresponse)
#endif
                     end if
                     if (fapar > epsilon(0.0)) then
                        ! Get the gext and the fields corresponding to a unit impulse in apar
                        call get_dgdfield_matrix_column(iky, ikx, iz, ie, idx, nz_ext, gext, field="apar")
                        call get_em_fields_for_response_matrix(gext, iky, ie, nresponse_per_field, fields_ext)
                        fields_ext(nresponse_per_field + idx) = fields_ext(nresponse_per_field + idx) - 1.0

#if !defined ISO_C_BINDING || !defined MPI
                        response_matrix(iky)%eigen(ie)%zloc(:, (nresponse_per_field + idx)) = -fields_ext(:nresponse)
#else
                        if (sgproc0) response_matrix(iky)%eigen(ie)%zloc(:, (nresponse_per_field + idx)) = -fields_ext(:nresponse)
#endif
                     end if

                     if (fbpar > epsilon(0.0)) then
                        ! Get the gext and the fields corresponding to a unit impulse in bpar
                        call get_dgdfield_matrix_column(iky, ikx, iz, ie, idx, nz_ext, gext, field="bpar")
                        call get_em_fields_for_response_matrix(gext, iky, ie, nresponse_per_field, fields_ext)
                        fields_ext(2 * nresponse_per_field + idx) = fields_ext(2 * nresponse_per_field + idx) - 1.0

#if !defined ISO_C_BINDING || !defined MPI
                        response_matrix(iky)%eigen(ie)%zloc(:, (2 * nresponse_per_field + idx)) = -fields_ext(:nresponse)
#else
                        if (sgproc0) response_matrix(iky)%eigen(ie)%zloc(:, (2 * nresponse_per_field + idx)) = -fields_ext(:nresponse)
#endif
                     end if

                  end do
                  if (izl_offset == 0) izl_offset = 1
               end do
            end if
            deallocate (gext, fields_ext)
         end do

         !DSO (editted by RJD) - This ends parallelization over velocity space.
         !      At this point every processor has
         !      (response_matrix = I - response matrix) for a given ky
         !      and so the LU decomposition can be parallelized locally if need be.
         !      This is preferable to parallelization over ky as the LU
         !      decomposition will be dominated by the ky with the most connections.
         !      Field equations are already calculated by this point.

         if (proc0 .and. debug) then
            call time_message(.true., time_response_matrix_dgdphi, message_dgdphi)
            call time_message(.false., time_response_matrix_lu, message_lu)
         end if
#ifdef ISO_C_BINDING
         call mpi_win_fence(0, window, ierr)
#endif

         !now we have the full response matrix. Finally, perform its LU decomposition
#ifdef MPI
         select case (lu_option_switch)
         case (lu_option_global)
            call parallel_LU_decomposition_global(iky)
         case (lu_option_local)
#ifdef ISO_C_BINDING
            call parallel_LU_decomposition_local(iky)
#else
            call mp_abort('Stella must be built with HAS_ISO_BINDING in order to use local parallel LU decomposition.')
#endif
         case default
#endif
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
#ifdef MPI
         end select
#endif

         if (proc0 .and. debug) then
            call time_message(.true., time_response_matrix_lu, message_lu)
         end if

         time_response_matrix_dgdphi = 0
         time_response_matrix_QN = 0
         time_response_matrix_lu = 0

         do ie = 1, neigen(iky)
            if (proc0 .and. mat_gen) then
               write (unit=mat_unit) response_matrix(iky)%eigen(ie)%idx
               write (unit=mat_unit) response_matrix(iky)%eigen(ie)%zloc
            end if
         end do

         !if(proc0)  write (*,*) 'job', iky, iproc, response_matrix(iky)%eigen(1)%zloc(5,:)
      end do

#ifdef ISO_C_BINDING
      call mpi_win_fence(0, window, ierr)
#endif

      if (proc0 .and. mat_gen) then
         close (unit=mat_unit)
      end if

      if (proc0 .and. debug) then
         write (*, '(A)') "    ############################################################"
         write (*, '(A)') " "
      end if

   end subroutine init_response_matrix

   subroutine read_response_matrix

      use fields_arrays, only: response_matrix
      use common_types, only: response_matrix_type
      use kt_grids, only: naky
      use extended_zgrid, only: neigen
      use extended_zgrid, only: nsegments
      use extended_zgrid, only: nzed_segment
      use extended_zgrid, only: periodic
      use mp, only: proc0, job, broadcast, mp_abort

      implicit none

      integer :: iky, ie, nz_ext
      integer :: iky_dump, neigen_dump, naky_dump, nresponse_dump
      integer :: nresponse
      character(len=15) :: fmt, job_str
      character(len=100) :: file_name
      integer :: istatus, ie_dump, istat
      logical, parameter :: debug = .false.
      istatus = 0

!   All matrices handled for the job i_job are read
!   from a single file named: responst_mat.ijob by that
!   jobs root process

      if (proc0) then
         fmt = '(I5.5)'
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
               nresponse = nz_ext - 1
            else
               nresponse = nz_ext
            end if

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

   subroutine get_dgdfield_matrix_column(iky, ikx, iz, ie, idx, nz_ext, gext, field)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use extended_zgrid, only: periodic
      use run_parameters, only: maxwellian_inside_zed_derivative
      use parallel_streaming, only: z_tridiagonal_solve
      use parallel_streaming, only: stream_sign
      use run_parameters, only: zed_upwind, time_upwind
      use finite_differences, only: tridag
#if defined ISO_C_BINDING && defined MPI
      use mp, only: mp_abort
#endif

      implicit none

      integer, intent(in) :: iky, ikx, iz, ie, idx, nz_ext
      character(*), intent(in) :: field
      complex, dimension(:, vmu_lo%llim_proc:), intent(in out) :: gext

      integer :: ivmu, iv, is, ia
      complex, dimension(:), allocatable :: a, b, c

      ia = 1

      ! a, b, c represent entries on the LHS of the matrix equation.
      allocate (a(1:nz_ext))
      allocate (b(1:nz_ext))
      allocate (c(1:nz_ext))

      if (.not. maxwellian_inside_zed_derivative) then
         ! get -vpa*b.gradz*Ze/T*F0*d<chi>/dz corresponding to unit impulse in a field;
         ! this could be phi, apar or bpar
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            call get_rhs_homogeneous_equation(iky, ikx, iz, ia, idx, nz_ext, ivmu, gext, field)
            call get_lhs_homogeneous_equation(iky, ikx, ia, nz_ext, ivmu, a, b, c)
            call tridag(1, a, b, c, gext(:, ivmu))
            ! hack for now (duplicates much of the effort from sweep_zed_zonal)
            if (periodic(iky)) then
               call mp_abort("Not set up the zonal flow calculation of ghom in response_matrix")
               call sweep_zed_zonal_response(iv, is, stream_sign(iv), gext(:, ivmu))
            else
               ! invert parallel streaming equation to get g^{n+1} on extended zed grid
               ! (I + (1+alph)/2*dt*vpa)*g_{inh}^{n+1} = RHS = gext
               call z_tridiagonal_solve(iky, ie, ivmu, gext(:, ivmu))
            end if

         end do
      else
         call mp_abort("maxwellian_inside_zed_derivative not currently supported in electromagnetic response matrix calculation")
      end if

   end subroutine get_dgdfield_matrix_column

   subroutine get_lhs_homogeneous_equation(iky, ikx, ia, nz_ext, ivmu, a, b, c)

      use parallel_streaming, only: stream_sign
      use run_parameters, only: drifts_implicit_in_z, stream_implicit
      use run_parameters, only: zed_upwind, time_upwind
      use stella_layouts, only: iv_idx, is_idx
      use stella_layouts, only: vmu_lo

      implicit none

      integer, intent(in) :: iky, ikx, ia, nz_ext, ivmu
      complex, dimension(:), intent(in out) :: a, b, c
      integer :: iv, is
      real, dimension(:), allocatable :: stream_a, stream_b, stream_c
      complex, dimension(:), allocatable :: wdrift_a, wdrift_b, wdrift_c

      ! LHS looks like
      ! dg/dt + vpa.stm.gradpar dg/dz + (ikx.wdriftx_g + iky.wdrifty_g)
      ! Discretised, it looks like
      !

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      if (stream_sign(iv) < 0) then
         a(1) = 0.
         a(2:) = 0.5 * (1 - zed_upwind)
         b(:) = 0.5 * (1 + zed_upwind)
         c(:) = 0.
      else
         a(:) = 0.
         b(:) = 0.5 * (1 + zed_upwind)
         c(1:nz_ext - 1) = 0.5 * (1 - zed_upwind)
         c(nz_ext) = 0.
      end if

      if (stream_implicit) then
         ! Could avoid allocating/deallocating each time by allocating once at the
         ! start and deallocating at the end.
         allocate (stream_a(1:nz_ext))
         allocate (stream_b(1:nz_ext))
         allocate (stream_c(1:nz_ext))
         iv = iv_idx(vmu_lo, ivmu)
         call get_lhs_streaming_term(ia, nz_ext, iv, is, stream_a, stream_b, stream_c)
         a = a + stream_a
         b = b + stream_b
         c = c + stream_c
         deallocate (stream_a)
         deallocate (stream_b)
         deallocate (stream_c)
      end if

      if (drifts_implicit_in_z) then
         allocate (wdrift_a(1:nz_ext))
         allocate (wdrift_b(1:nz_ext))
         allocate (wdrift_c(1:nz_ext))
         call get_lhs_wdrift_term(iky, ikx, ia, nz_ext, ivmu, wdrift_a, wdrift_b, wdrift_c)
         a = a + wdrift_a
         b = b + wdrift_b
         c = c + wdrift_c
         deallocate (wdrift_a)
         deallocate (wdrift_b)
         deallocate (wdrift_c)
      end if

   end subroutine get_lhs_homogeneous_equation

   subroutine get_rhs_homogeneous_equation(iky, ikx, iz, ia, idx, nz_ext, ivmu, gext, field)

      use constants, only: zi
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use kt_grids, only: zonal_mode, nakx, naky!, aky, akx
      use species, only: spec
      use vpamu_grids, only: vpa, mu
      use gyro_averages, only: aj0x, aj1x
      use run_parameters, only: driftkinetic_implicit, drifts_implicit_in_z, stream_implicit
      use run_parameters, only: maxwellian_inside_zed_derivative
      use parallel_streaming, only: stream_sign
      use run_parameters, only: zed_upwind, time_upwind
#if defined ISO_C_BINDING && defined MPI
      use mp, only: mp_abort
#endif

      implicit none

      integer, intent(in) :: iky, ikx, iz, ia, idx, nz_ext, ivmu
      character(*), intent(in) :: field
      complex, dimension(:, vmu_lo%llim_proc:), intent(in out) :: gext

      integer :: iv, imu, is
      real :: stream_fac0, stream_fac1
      complex :: wdrift_fac0, wdrift_fac1, wstar_fac0, wstar_fac1
      real :: gyro_fac

      ! initialize g to zero everywhere along extended zed domain
      gext(:, ivmu) = 0.0
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      ! give unit impulse to phi at this zed location
      ! and compute -vpa*b.gradz*Ze/T*d<phi>/dz*F0 (RHS of streaming part of GKE)

      ! NB:  assuming equal spacing in zed below
      ! here, fac = -dt*(1+alph_t)/2*vpa*Ze/T*F0*(gyro_fac)/dz
      ! b.gradz left out because needs to be centred in zed
      ! gyro_fac is a term corresponding to a gyroaveraged element of chi
      if (driftkinetic_implicit) then
         call mp_abort("driftkinetic_implicit not supported electromagnetically")
         gyro_fac = 1.0
      else
         ! Calculate the factor corresponding to a unit impulse in phi, apar
         ! or bpar
         ! <chi> = J0*phi - 2*vpa*vthermal*J0*apar + 4*mu*(T/Z)*aj1*bpar
         if (field == 'phi') then
            gyro_fac = aj0x(iky, ikx, iz, ivmu)
         else if (field == 'apar') then
            gyro_fac = -2 * vpa(iv) * spec(is)%stm_psi0 * aj0x(iky, ikx, iz, ivmu)
         else if (field == 'bpar') then
            gyro_fac = 4 * mu(imu) * spec(is)%tz_psi0 * aj1x(iky, ikx, iz, ivmu)
         else
            ! If this occurs, we've hit a problem - end gracefully
            call mp_abort("Error in calculation of repsonse matrix")
         end if
      end if

      if (stream_implicit) then
         call get_rhs_streaming_term(iz, ia, ivmu, gyro_fac, stream_fac0, stream_fac1)
      else
         stream_fac0 = 0.
         stream_fac1 = 0.
      end if

      if (drifts_implicit_in_z) then
         call get_rhs_wdrift_term(iky, ikx, iz, ia, ivmu, gyro_fac, wdrift_fac0, wdrift_fac1)
      else
         wdrift_fac0 = 0.
         wdrift_fac1 = 0.
      end if

      ! Could include in the earlier if block, but we might want to separate wdrift
      ! and wstar in the future
      if (drifts_implicit_in_z) then
         call get_rhs_wstar_term(iky, iz, ia, ivmu, gyro_fac, wstar_fac0, wstar_fac1)
      else
         wstar_fac0 = 0.
         wstar_fac1 = 0.
      end if

      ! fac0 is the factor multiplying delphi on the RHS
      ! of the homogeneous GKE at this zed index
      gext(idx, ivmu) = stream_fac0 + wdrift_fac0 + wstar_fac0

      if (stream_sign(iv) < 0) then
         ! XXX_fac1 is the factor multiplying delphi on the RHS
         ! of the homogeneous GKE at the zed index to the right of
         ! this one
         if (idx < nz_ext) gext(idx + 1, ivmu) = stream_fac1 + wdrift_fac1 + wstar_fac1
         ! zonal mode BC is periodic instead of zero, so must
         ! treat specially
         if (zonal_mode(iky)) then
            if (idx == 1) then
               gext(nz_ext, ivmu) = stream_fac0 + wdrift_fac0 + wstar_fac0
            else if (idx == nz_ext - 1) then
               gext(1, ivmu) = stream_fac1 + wdrift_fac1 + wstar_fac1
            end if
         end if
      else
         ! XXX_fac1 is the factor multiplying delphi on the RHS
         ! of the homogeneous GKE at the zed index to the left of
         ! this one
         if (idx > 1) gext(idx - 1, ivmu) = stream_fac1 + wdrift_fac1 + wstar_fac1
         ! zonal mode BC is periodic instead of zero, so must
         ! treat specially
         if (zonal_mode(iky)) then
            if (idx == 1) then
               gext(nz_ext, ivmu) = stream_fac0 + wdrift_fac0 + wstar_fac0
               gext(nz_ext - 1, ivmu) = stream_fac1 + wdrift_fac1 + wstar_fac1
            else if (idx == 2) then
               gext(nz_ext, ivmu) = stream_fac1 + wdrift_fac1 + wstar_fac1
            end if
         end if
      end if

   end subroutine get_rhs_homogeneous_equation

   subroutine get_rhs_streaming_term(iz, ia, ivmu, gyro_fac, stream_fac0, stream_fac1)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use stella_time, only: code_dt
      use zgrid, only: delzed, nzgrid
      use species, only: spec
      use stella_geometry, only: gradpar
      use vpamu_grids, only: vpa
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use parallel_streaming, only: stream_sign
      use run_parameters, only: zed_upwind, time_upwind

      implicit none

      integer, intent(in) :: iz, ia, ivmu
      real, intent(in) :: gyro_fac
      real, intent(in out) :: stream_fac0, stream_fac1

      integer :: iv, imu, is
      real :: stream_fac, gradpar_left, gradpar_right
      real :: maxwell_mu_left, maxwell_mu_right

      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      ! 0.125 to account for two linear interpolations (gradpar and maxwell_mu)
      ! (each contributing a factor of 1/2) and the (1+time_upwind)/2 factor
      stream_fac = -0.125 * (1.+time_upwind) * code_dt * vpa(iv) * spec(is)%stm_psi0 &
                   * gyro_fac * spec(is)%zt / delzed(0) * maxwell_vpa(iv, is) * maxwell_fac(is)

      ! Because our discretisation couples neighboring z gridpoints, we
      ! need to find quantities (gradpar, maxwell_mu) at the gridpoint to the
      ! left and right of this gridpoint. If this gridpoint is the
      ! leftmost/rightmost point of end of the extended z grid, we need tp get the
      ! data point from beyond our extended z-grid. Since (gradpar, maxwell_mu)
      ! are periodic, this is fairly straightforward.

      ! Get quantities on gridpoint to the right.
      if (iz < nzgrid) then
         gradpar_right = gradpar(iz + 1)
         maxwell_mu_right = maxwell_mu(ia, iz + 1, imu, is)
      else
         ! periodic - take from the left
         gradpar_right = gradpar(-nzgrid + 1)
         maxwell_mu_right = maxwell_mu(ia, -nzgrid + 1, imu, is)
      end if

      ! Get quantities on gridpoint to the left.
      if (iz > -nzgrid) then
         gradpar_left = gradpar(iz - 1)
         maxwell_mu_left = maxwell_mu(ia, iz - 1, imu, is)
      else
         ! periodic - take from the right
         gradpar_left = gradpar(nzgrid - 1)
         maxwell_mu_left = maxwell_mu(ia, nzgrid - 1, imu, is)
      end if

      ! In the following, gradpar and maxwell_mu are interpolated separately
      ! to ensure consistency to what is done in parallel_streaming.f90
      if (stream_sign(iv) < 0) then
         stream_fac0 = stream_fac * ((1.+zed_upwind) * gradpar(iz) &
                                     + (1.-zed_upwind) * gradpar_left) &
                       * ((1.+zed_upwind) * maxwell_mu(ia, iz, imu, is) &
                          + (1.-zed_upwind) * maxwell_mu_left)
         stream_fac1 = -stream_fac * ((1.+zed_upwind) * gradpar_right &
                                      + (1.-zed_upwind) * gradpar(iz)) &
                       * ((1.+zed_upwind) * maxwell_mu_right &
                          + (1.-zed_upwind) * maxwell_mu(ia, iz, imu, is))
      else
         stream_fac0 = -stream_fac * ((1.+zed_upwind) * gradpar(iz) &
                                      + (1.-zed_upwind) * gradpar_right) &
                       * ((1.+zed_upwind) * maxwell_mu(ia, iz, imu, is) &
                          + (1.-zed_upwind) * maxwell_mu_right)
         stream_fac1 = stream_fac * ((1.+zed_upwind) * gradpar_left &
                                     + (1.-zed_upwind) * gradpar(iz)) &
                       * ((1.+zed_upwind) * maxwell_mu_left &
                          + (1.-zed_upwind) * maxwell_mu(ia, iz, imu, is))
      end if

   end subroutine get_rhs_streaming_term

   subroutine get_rhs_wdrift_term(iky, ikx, iz, ia, ivmu, gyro_fac, wdrift_fac0, wdrift_fac1)

      use constants, only: zi
      use kt_grids, only: aky, akx
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx
      use zgrid, only: nzgrid
      use parallel_streaming, only: stream_sign
      use run_parameters, only: zed_upwind, time_upwind
      use dist_fn_arrays, only: wdrifty_phi, wdriftx_phi

      implicit none

      integer, intent(in) :: iky, ikx, iz, ia, ivmu
      real, intent(in) :: gyro_fac
      complex, intent(in out) :: wdrift_fac0, wdrift_fac1

      integer :: iv
      real :: wdrift_fac
      complex :: wdrifty_phi_left, wdrifty_phi_right, wdriftx_phi_left, wdriftx_phi_right

      ! 0.125 to account for one linear interpolation (wdriftx,y_phi), a
      ! (1+time_upwind)/2 factor and a (1+/-zed_upwind)/2 factor
      wdrift_fac = 0.125 * (1 + time_upwind) * gyro_fac

      ! Because our discretisation couples neighboring z gridpoints, we
      ! need to find quantities (gradpar, maxwell_mu) at the gridpoint to the
      ! left and right of this gridpoint. If this gridpoint is the
      ! leftmost/rightmost point of end of the extended z grid, we need tp get the
      ! data point from beyond our extended z-grid. Since (gradpar, maxwell_mu)
      ! are periodic, this is fairly straightforward.

      ! RJD: We're inconsistenly applying BCs here - we can't use
      ! periodicity to get wdrifty_phi(iz+1), because wdrifty_phi
      ! not periodic. Instead assume wdrifty_phi(iz+1)=wdrifty_phi(iz);
      ! not correct but hopefully won't cause any problems (expect
      ! distribution function to be small here.)
      ! I belive wdriftx_phi is periodic, but use the same BC
      ! as wdrifty_phi until we can prove this.

      ! Get quantities on gridpoint to the right.
      if (iz < nzgrid) then
         wdrifty_phi_right = wdrifty_phi(ia, iz + 1, ivmu)
         wdriftx_phi_right = wdriftx_phi(ia, iz + 1, ivmu)
      else
         ! Off-grid problem: replicate the rightmost point
         wdrifty_phi_right = wdrifty_phi(ia, nzgrid, ivmu)
         wdriftx_phi_right = wdriftx_phi(ia, nzgrid, ivmu)
      end if

      ! Get quantities on gridpoint to the left.
      if (iz > -nzgrid) then
         wdrifty_phi_left = wdrifty_phi(ia, iz - 1, ivmu)
         wdriftx_phi_left = wdriftx_phi(ia, iz - 1, ivmu)
      else
         ! Off-grid problem: replicate the leftmost point
         wdrifty_phi_left = wdrifty_phi(ia, -nzgrid, ivmu)
         wdriftx_phi_left = wdriftx_phi(ia, -nzgrid, ivmu)
      end if

      iv = iv_idx(vmu_lo, ivmu)
      if (stream_sign(iv) < 0) then
         wdrift_fac0 = wdrift_fac * (1.+zed_upwind) * (zi * aky(iky) &
                                                       * ((1.+zed_upwind) * wdrifty_phi(ia, iz, ivmu) + (1.-zed_upwind) * wdrifty_phi_left) &
                                                       + zi * akx(ikx) * ((1.+zed_upwind) * wdriftx_phi(ia, iz, ivmu) &
                                                                          + (1.-zed_upwind) * wdriftx_phi_left))
         wdrift_fac1 = wdrift_fac * (1.-zed_upwind) * ((zi * aky(iky) &
                                                        * ((1.+zed_upwind) * wdrifty_phi_right + (1.-zed_upwind) * wdrifty_phi(ia, iz, ivmu)) &
                                                        + zi * akx(ikx) * ((1.+zed_upwind) * wdriftx_phi_right &
                                                                           + (1.-zed_upwind) * wdriftx_phi(ia, iz, ivmu))))
      else
         wdrift_fac0 = wdrift_fac * (1.+zed_upwind) * (zi * aky(iky) &
                                                       * ((1.-zed_upwind) * wdrifty_phi_right + (1.+zed_upwind) * wdrifty_phi(ia, iz, ivmu)) &
                                                       + zi * akx(ikx) * ((1.-zed_upwind) * wdrifty_phi_right &
                                                                          + (1.+zed_upwind) * wdriftx_phi(ia, iz, ivmu)))
         wdrift_fac1 = wdrift_fac * (1.-zed_upwind) * (zi * aky(iky) &
                                                       * ((1.-zed_upwind) * wdrifty_phi(ia, iz, ivmu) + (1.+zed_upwind) * wdrifty_phi_left) &
                                                       + zi * akx(ikx) * ((1.-zed_upwind) * wdriftx_phi(ia, iz, ivmu) &
                                                                          + (1.+zed_upwind) * wdriftx_phi_left))
      end if

   end subroutine get_rhs_wdrift_term

   subroutine get_rhs_wstar_term(iky, iz, ia, ivmu, gyro_fac, wstar_fac0, wstar_fac1)

      use constants, only: zi
      use kt_grids, only: aky
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx
      use zgrid, only: nzgrid
      use parallel_streaming, only: stream_sign
      use run_parameters, only: zed_upwind, time_upwind
      use dist_fn_arrays, only: wstar

      implicit none

      integer, intent(in) :: iky, iz, ia, ivmu
      real, intent(in) :: gyro_fac
      complex, intent(in out) :: wstar_fac0, wstar_fac1

      integer :: iv
      real :: wstar_fac
      complex :: wstar_left, wstar_right

      ! 0.125 to account for one linear interpolation (wstar), a
      ! (1+time_upwind)/2 factor and a (1+/-zed_upwind)/2 factor
      wstar_fac = 0.125 * (1 + time_upwind) * gyro_fac

      ! Because our discretisation couples neighboring z gridpoints, we
      ! need to find quantities (gradpar, maxwell_mu) at the gridpoint to the
      ! left and right of this gridpoint. If this gridpoint is the
      ! leftmost/rightmost point of end of the extended z grid, we need tp get the
      ! data point from beyond our extended z-grid. Since (gradpar, maxwell_mu)
      ! are periodic, this is fairly straightforward.
      ! As above - can we show wstar periodic? Until then, set
      ! wstar(iz+1) = wstar(iz)

      ! Get quantities on gridpoint to the right.
      if (iz < nzgrid) then
         wstar_right = wstar(ia, iz + 1, ivmu)
      else
         ! Off-grid problem: replicate the rightmost point
         wstar_right = wstar(ia, nzgrid, ivmu)
      end if

      ! Get quantities on gridpoint to the left.
      if (iz > -nzgrid) then
         wstar_left = wstar(ia, iz - 1, ivmu)
      else
         ! Off-grid problem: replicate the leftmost point
         wstar_left = wstar(ia, -nzgrid, ivmu)
      end if

      iv = iv_idx(vmu_lo, ivmu)
      if (stream_sign(iv) < 0) then
         wstar_fac0 = wstar_fac * (1.+zed_upwind) * (zi * aky(iky) &
                                                     * ((1.+zed_upwind) * wstar(ia, iz, ivmu) + (1.-zed_upwind) * wstar_left))
         wstar_fac1 = wstar_fac * (1.-zed_upwind) * (zi * aky(iky) &
                                                     * ((1.+zed_upwind) * wstar_right + (1.-zed_upwind) * wstar(ia, iz, ivmu)))
      else
         wstar_fac0 = wstar_fac * (1.+zed_upwind) * (zi * aky(iky) &
                                                     * ((1.-zed_upwind) * wstar_right + (1.+zed_upwind) * wstar(ia, iz, ivmu)))
         wstar_fac1 = wstar_fac * (1.-zed_upwind) * (zi * aky(iky) &
                                                     * ((1.-zed_upwind) * wstar(ia, iz, ivmu) + (1.+zed_upwind) * wstar_left))
      end if

   end subroutine get_rhs_wstar_term

   subroutine get_lhs_streaming_term(ia, nz_ext, iv, is, stream_a, stream_b, stream_c)

      use stella_time, only: code_dt
      use zgrid, only: delzed, nzgrid
      use species, only: spec
      use stella_geometry, only: gradpar
      use vpamu_grids, only: vpa
      use parallel_streaming, only: stream_sign
      use run_parameters, only: zed_upwind, time_upwind

      implicit none

      integer, intent(in) :: ia, iv, is, nz_ext
      real, dimension(:), intent(in out) :: stream_a, stream_b, stream_c
      real :: stream_fac
      real, dimension(:), allocatable :: gradpar_left, gradpar_right

      ! Get the elements a, b, c corresponding to the streaming term.
      ! NB gext has dimension(1:nz_ext, vmu), and need (a, b, c) to also have
      ! dimension (1:nz_ext). This provides a slight complication for
      ! our geometrical terms (e.g. gradpar), which have dimension (-nz_grid:nz_grid)
      ! 0.25 factor from (1+ut)/2 and interpolation of gradpar
      stream_fac = 0.25 * (1.+time_upwind) * code_dt * vpa(iv) * spec(is)%stm_psi0 / delzed(0)

      ! NB we could avoid allocating & deallocating every time by calculating
      ! once at the start and deallocating at the end.
      allocate (gradpar_left(nz_ext))
      allocate (gradpar_right(nz_ext))

      ! Get gradpar at the gridpoint to the left i.e. shift everything left by
      ! 1, but also account for the fact that gradpar is (-nz_grid:nz_grid) but
      ! gradpar_left is (1:nz_ext)
      gradpar_left(2:nz_ext) = gradpar(-nzgrid:nzgrid - 1)
      gradpar_left(1) = gradpar(nzgrid) ! Periodicity in gradpar
      ! Repeat for gradpar_right; shift everything right by 1
      gradpar_right(1:nz_ext - 1) = gradpar(-nzgrid + 1:nzgrid)
      gradpar_right(nz_ext) = gradpar(-nzgrid)  ! Periodicity in gradpar
      if (stream_sign(iv) < 0) then
         stream_a(1) = 0.
         stream_a(2:nz_ext) = -stream_fac * ((1 + zed_upwind) * gradpar(-nzgrid + 1:nzgrid) + (1 - zed_upwind) * gradpar_left(2:nz_ext))
         stream_b(:) = stream_fac * ((1 + zed_upwind) * gradpar(-nzgrid:nzgrid) + (1 - zed_upwind) * gradpar_left)
         stream_c(:) = 0.
      else
         stream_a(:) = 0.
         stream_b(:) = -stream_fac * ((1 - zed_upwind) * gradpar_right + (1 + zed_upwind) * gradpar(-nzgrid:nzgrid))
         stream_c(1:nz_ext - 1) = stream_fac * ((1 - zed_upwind) * gradpar_right(1:nz_ext - 1) + (1 + zed_upwind) * gradpar(-nzgrid:nzgrid - 1))
         stream_c(nz_ext) = 0.
      end if

      deallocate (gradpar_left)
      deallocate (gradpar_right)

   end subroutine get_lhs_streaming_term

   subroutine get_lhs_wdrift_term(iky, ikx, ia, nz_ext, ivmu, wdrift_a, wdrift_b, wdrift_c)

      use constants, only: zi
      use kt_grids, only: aky, akx
      use zgrid, only: nzgrid
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx
      use parallel_streaming, only: stream_sign
      use run_parameters, only: zed_upwind, time_upwind
      use dist_fn_arrays, only: wdrifty_g, wdriftx_g

      implicit none

      integer, intent(in) :: iky, ikx, ia, nz_ext, ivmu
      complex, dimension(:), intent(in out) :: wdrift_a, wdrift_b, wdrift_c
      integer :: iv
      real, dimension(:), allocatable  :: wdrifty_g_left, wdrifty_g_right, wdriftx_g_left, wdriftx_g_right
      real :: wdrift_fac

      ! 0.125 factor from (1+ut)/2, (1+/-uz)/2 and interpolations of wdriftx,y
      wdrift_fac = 0.125 * (1.+time_upwind)
      iv = iv_idx(vmu_lo, ivmu)

      allocate (wdrifty_g_left(nz_ext))
      allocate (wdrifty_g_right(nz_ext))
      allocate (wdriftx_g_left(nz_ext))
      allocate (wdriftx_g_right(nz_ext))

      ! Get the left-shifted values of wdrifty, accounting for the
      ! translation from (-nzgrid:nzgrid) to (1:nz_ext)
      wdrifty_g_left(2:nz_ext) = wdrifty_g(ia, -nzgrid:nzgrid - 1, ivmu)
      ! Off-grid problem: replicate the leftmost point
      wdrifty_g_left(1) = wdrifty_g(ia, -nzgrid, ivmu)
      wdrifty_g_right(1:nz_ext - 1) = wdrifty_g(ia, -nzgrid + 1:nzgrid, ivmu)
      ! Off-grid problem: replicate the rightmost point
      wdrifty_g_right(nz_ext) = wdrifty_g(ia, nzgrid, ivmu)

      ! Repeat for wdriftx_g
      wdriftx_g_left(2:nz_ext) = wdriftx_g(ia, -nzgrid:nzgrid - 1, ivmu)
      ! Off-grid problem: replicate the leftmost point
      wdriftx_g_left(1) = wdriftx_g(ia, -nzgrid, ivmu)
      wdriftx_g_right(1:nz_ext - 1) = wdriftx_g(ia, -nzgrid + 1:nzgrid, ivmu)
      ! Off-grid problem: replicate the rightmost point
      wdriftx_g_right(nz_ext) = wdriftx_g(ia, -nzgrid, ivmu)

      if (stream_sign(iv) < 0) then
         wdrift_a(1) = 0.
   wdrift_a(2:) = wdrift_fac * (zi * aky(iky) * ((1 + zed_upwind) * wdrifty_g(ia, -nzgrid + 1:nzgrid, ivmu) + (1 - zed_upwind) * wdrifty_g_left(2:)) &
                               + zi * akx(ikx) * ((1 + zed_upwind) * wdriftx_g(ia, -nzgrid + 1:nzgrid, ivmu) + (1 - zed_upwind) * wdriftx_g_left(2:)))
         wdrift_b(:) = wdrift_fac * (zi * aky(iky) * ((1 + zed_upwind) * wdrifty_g(ia, -nzgrid:nzgrid, ivmu) + (1 - zed_upwind) * wdrifty_g_left) &
                                     + zi * akx(ikx) * ((1 + zed_upwind) * wdriftx_g(ia, -nzgrid:nzgrid, ivmu) + (1 - zed_upwind) * wdriftx_g_left))
         wdrift_c(:) = 0.
      else
         wdrift_a(:) = 0.
        wdrift_b(:) = wdrift_fac * (zi * aky(iky) * ((1 - zed_upwind) * wdrifty_g_right(:) + (1 + zed_upwind) * wdrifty_g(ia, -nzgrid:nzgrid, ivmu)) &
                                   + zi * akx(ikx) * ((1 - zed_upwind) * wdriftx_g_right(:) + (1 + zed_upwind) * wdriftx_g(ia, -nzgrid:nzgrid, ivmu)))
    wdrift_c(1:nz_ext-1) = wdrift_fac*(zi*aky(iky)*((1-zed_upwind)*wdrifty_g_right(1:nz_ext-1) + (1+zed_upwind)*wdrifty_g(ia,-nzgrid:nzgrid-1,ivmu)) &
                    + zi * akx(ikx) * ((1 - zed_upwind) * wdriftx_g_right(1:nz_ext - 1) + (1 + zed_upwind) * wdriftx_g(ia, -nzgrid:nzgrid - 1, ivmu)))
         wdrift_c(nz_ext) = 0.
      end if

      deallocate (wdrifty_g_left)
      deallocate (wdrifty_g_right)
      deallocate (wdriftx_g_left)
      deallocate (wdriftx_g_right)

   end subroutine get_lhs_wdrift_term

! subroutine get_phi_matrix
! end subroutine get_phi_matrix

   subroutine sweep_zed_zonal_response(iv, is, sgn, g)

      use zgrid, only: nzgrid, delzed, nztot
      use run_parameters, only: zed_upwind, time_upwind
      use parallel_streaming, only: stream_c

      implicit none

      integer, intent(in) :: iv, is, sgn
      complex, dimension(:), intent(in out) :: g
      integer :: iz, iz1, iz2
      real :: fac1, fac2
      complex, dimension(:), allocatable :: gcf, gpi

      allocate (gpi(-nzgrid:nzgrid))
      allocate (gcf(-nzgrid:nzgrid))
      ! ky=0 is 2pi periodic (no extended zgrid)
      ! decompose into complementary function + particular integral
      ! zero BC for particular integral
      ! unit BC for complementary function (no source)
      if (sgn < 0) then
         iz1 = -nzgrid; iz2 = nzgrid
      else
         iz1 = nzgrid; iz2 = -nzgrid
      end if
      gpi(iz1) = 0.; gcf(iz1) = 1.
      do iz = iz1 - sgn, iz2, -sgn
         fac1 = 1.0 + zed_upwind + sgn * (1.0 + time_upwind) * stream_c(iz, iv, is) / delzed(0)
         fac2 = 1.0 - zed_upwind - sgn * (1.0 + time_upwind) * stream_c(iz, iv, is) / delzed(0)
         gpi(iz) = (-gpi(iz + sgn) * fac2 + 2.0 * g(iz + nzgrid + 1)) / fac1
         gcf(iz) = -gcf(iz + sgn) * fac2 / fac1
      end do
      ! g = g_PI + (g_PI(pi)/(1-g_CF(pi))) * g_CF
      g = gpi + (spread(gpi(iz2), 1, nztot) / (1.-gcf(iz2))) * gcf
      deallocate (gpi, gcf)

   end subroutine sweep_zed_zonal_response

   subroutine integrate_over_velocity(g, phi, iky, ie)

      use stella_layouts, only: vmu_lo
      use species, only: nspec, spec
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: nsegments
      use kt_grids, only: zonal_mode, akx
      use vpamu_grids, only: integrate_species
      use gyro_averages, only: gyro_average
      use mp, only: sum_allreduce

      implicit none

      complex, dimension(:, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:), intent(out) :: phi
      integer, intent(in) :: iky, ie

      integer :: idx, iseg, ikx, iz, ia
      integer :: izl_offset
      real, dimension(nspec) :: wgt
      complex, dimension(:), allocatable :: g0

      ia = 1

      allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      wgt = spec%z * spec%dens_psi0
      phi = 0.

      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      if (zonal_mode(iky) .and. abs(akx(ikx)) < epsilon(0.)) then
         phi(:) = 0.0
         return
      end if
      do iz = iz_low(iseg), iz_up(iseg)
         idx = idx + 1
         call gyro_average(g(idx, :), iky, ikx, iz, g0)
         call integrate_species(g0, iz, wgt, phi(idx), reduce_in=.false.)
      end do
      izl_offset = 1
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               call gyro_average(g(idx, :), iky, ikx, iz, g0)
               call integrate_species(g0, iz, wgt, phi(idx), reduce_in=.false.)
            end do
            if (izl_offset == 0) izl_offset = 1
         end do
      end if

      call sum_allreduce(phi)

   end subroutine integrate_over_velocity

   subroutine get_fields_for_response_matrix(phi, iky, ie)

      use stella_layouts, only: vmu_lo
      use species, only: spec
      use species, only: has_electron_species
      use stella_geometry, only: dl_over_b
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: nsegments
      use kt_grids, only: zonal_mode, akx
      use dist_fn, only: adiabatic_option_switch
      use dist_fn, only: adiabatic_option_fieldlineavg
      use fields, only: gamtot3
      use fields_arrays, only: gamtot

      implicit none

      complex, dimension(:), intent(inout) :: phi
      integer, intent(in) :: iky, ie

      integer :: idx, iseg, ikx, iz, ia
      integer :: izl_offset
      complex, dimension(:), allocatable :: g0
      complex :: tmp

      ia = 1

      allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      if (zonal_mode(iky) .and. abs(akx(ikx)) < epsilon(0.)) then
         phi(:) = 0.0
         return
      end if
      do iz = iz_low(iseg), iz_up(iseg)
         idx = idx + 1
         phi(idx) = phi(idx) / gamtot(iky, ikx, iz)
      end do
      izl_offset = 1
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               phi(idx) = phi(idx) / gamtot(iky, ikx, iz)
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
            phi = phi + tmp * gamtot3(ikxmod(1, ie, iky), :)
         end if
      end if
      deallocate (g0)

   end subroutine get_fields_for_response_matrix

   ! Given gext == gbar on extended zgrid, calculate (phi, apar, bpar) on the
   ! extended zgrid and return in a one-dimenstional array.
   subroutine get_em_fields_for_response_matrix(gext, iky, ie, nresponse_per_field, fields_ext)

      use stella_layouts, only: vmu_lo
      use fields, only: get_fields_vmulo_single
      use mp, only: mp_abort
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: nsegments
      use kt_grids, only: zonal_mode, akx

      implicit none

      complex, dimension(:, vmu_lo%llim_proc:), intent(in) :: gext
      integer, intent(in) :: iky, ie, nresponse_per_field
      complex, dimension(:), intent(inout) :: fields_ext

      integer :: ikx, iseg, iz, idx, izl_offset
      complex :: phi, apar, bpar

      ! nresponse_per_field is nz_ext for non-zonal modes, (nz_ext-1) for zonal modes (periodic BCs).
      ! This means that for zonal modes, fields_ext is still nz_ext*3 in length,
      ! but the last 3 elements are not calculated here (and are subsequently ignored)

      ! Need to get (phi, apar, bpar) for the extended domain of connected kx
      ! values; this means looping over each segment of the extended domain.
      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      if (zonal_mode(iky) .and. abs(akx(ikx)) < epsilon(0.)) then
         fields_ext(:) = 0.0
         return
      end if

      do iz = iz_low(iseg), iz_up(iseg)
         idx = idx + 1
         call get_fields_vmulo_single(gext(idx, :), iky, ikx, iz, phi, apar, bpar, dist="gbar")
         fields_ext(idx) = phi
         fields_ext(nresponse_per_field + idx) = apar
         fields_ext(2 * nresponse_per_field + idx) = bpar
      end do

      izl_offset = 1
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               call get_fields_vmulo_single(gext(idx, :), iky, ikx, iz, phi, apar, bpar, dist="gbar")
               fields_ext(idx) = phi
               fields_ext(nresponse_per_field + idx) = apar
               fields_ext(2 * nresponse_per_field + idx) = bpar
            end do
            ! RJD: Don't know what the purpose of this line of code is,
            ! but copying integrate_over_velocity
            if (izl_offset == 0) izl_offset = 1
         end do
      end if
   end subroutine get_em_fields_for_response_matrix

   subroutine finish_response_matrix

      use fields_arrays, only: response_matrix
#if !defined ISO_C_BINDING
      implicit none

#else
      use mpi

      implicit none

      integer :: ierr

      if (window /= MPI_WIN_NULL) call mpi_win_free(window, ierr)
#endif

      if (allocated(response_matrix)) deallocate (response_matrix)
      response_matrix_initialized = .false.

   end subroutine finish_response_matrix

#ifdef MPI
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
      use fields_arrays, only: response_matrix
      use mp, only: barrier, broadcast, sum_allreduce
      use mp, only: mp_comm, scope, allprocs, sharedprocs, curr_focus
      use mp, only: scrossdomprocs, sgproc0, mp_abort, real_size
      use mp, only: job, iproc, proc0, nproc, numnodes, inode
      use mp_lu_decomposition, only: lu_decomposition_local
      use job_manage, only: njobs
      use extended_zgrid, only: neigen
      use mpi
      use linear_solve, only: imaxloc

      implicit none

      integer, intent(in) :: iky

      integer, dimension(:, :), allocatable :: eig_limits
      integer, dimension(:), allocatable :: job_list
      integer, dimension(:), allocatable :: row_limits
      logical, dimension(:, :), allocatable :: node_jobs

      real, parameter :: zero = 1.0e-20
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

      allocate (node_jobs(numnodes, njobs)); node_jobs = .false.
      allocate (job_list(nproc)); job_list = 0
      allocate (row_limits(0:nproc))
      allocate (eig_limits(0:numnodes, njobs)); eig_limits = 0

      job_list(iproc + 1) = job
      call sum_allreduce(job_list)

      if (proc0) then
         do j = 1, nproc
            node_jobs(inode + 1, job_list(j) + 1) = .true. !create a map of which nodes have which jobs
         end do
      end if

      !make sure all processors have this map
      call scope(allprocs)
      call mpi_allreduce &
         (MPI_IN_PLACE, node_jobs, size(node_jobs), MPI_LOGICAL, MPI_LOR, mp_comm, ierr)
      call scope(sharedprocs)

      do ijob = 0, njobs - 1
         jroot = -1
         do j = 1, nproc
            if (job_list(j) == ijob) then
               jroot = j - 1 !the first processor on this job will be the root process
               exit
            end if
         end do

         if (jroot == -1) cycle !no processors on this node are on this job

         if (iproc == jroot) neig = neigen(iky)

         ! broadcast number of matrices
         call broadcast(neig, jroot)

         ! split up neig across nodes that have the current job
         nodes_on_job = count(node_jobs(:, ijob + 1))
         ediv = neig / nodes_on_job
         emod = mod(neig, nodes_on_job)

         eig_limits(0, ijob + 1) = 1
         do j = 1, numnodes
            if (node_jobs(j, ijob + 1)) then
               eig_limits(j, ijob + 1) = eig_limits(j - 1, ijob + 1) + ediv
               if (emod > 0) then
                  eig_limits(j, ijob + 1) = eig_limits(j, ijob + 1) + 1
                  emod = emod - 1
               end if
            else
               eig_limits(j, ijob + 1) = eig_limits(j - 1, ijob + 1)
            end if
         end do

         do ie = eig_limits(inode, ijob + 1), eig_limits(inode + 1, ijob + 1) - 1
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
                (ie >= eig_limits(inode, job + 1) .and. ie < eig_limits(inode + 1, job + 1))) nroot = iproc
            !first let processors know who is sending the data
            call sum_allreduce(nroot)
            !now send the data
            call broadcast(response_matrix(iky)%eigen(ie)%zloc, nroot)
            call broadcast(response_matrix(iky)%eigen(ie)%idx, nroot)
         end do
      end if

      call scope(prior_focus)

      deallocate (node_jobs, job_list, row_limits, eig_limits)
   end subroutine parallel_LU_decomposition_local

#endif /* ISO_C_BINDING */

   !this subroutine parallelizes the LU decomposition across
   !all cores. Ideal speed up: ncores
   subroutine parallel_LU_decomposition_global(iky)

      use fields_arrays, only: response_matrix
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

#endif /* MPI */

end module response_matrix
