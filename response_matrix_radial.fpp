module response_matrix_radial

   use mpi

   implicit none

   public :: init_response_matrix, finish_response_matrix
   public :: response_matrix_initialized

   private

   logical :: response_matrix_initialized = .false.

#if defined MPI && defined ISO_C_BINDING
   integer :: window = MPI_WIN_NULL
#endif

contains

   subroutine init_response_matrix

      use linear_solve, only: lu_decomposition
      use fields_arrays, only: response_matrix
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, is_idx
      use kt_grids, only: naky, nakx
      use zgrid, only: ntubes
      use full_xzgrid, only: nelements
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: neigen, ikxmod
      use extended_zgrid, only: nsegments
      use extended_zgrid, only: nzed_segment
      use extended_zgrid, only: periodic
      use job_manage, only: time_message
      use mp, only: proc0, mp_abort
      use run_parameters, only: lu_option_switch
      use run_parameters, only: lu_option_none, lu_option_local, lu_option_global
      use system_fortran, only: systemf
#if defined MPI && defined ISO_C_BINDING
      use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_intptr_t
      use mp, only: curr_focus, sgproc0, mp_comm, sharedsubprocs, scope, barrier
      use mp, only: real_size, nbytes_real
      use mpi
#endif

      implicit none

      integer :: iky, ie, iseg, iz
      integer :: ikx
      integer :: nz_ext, nresponse
      integer :: idx
      integer :: izl_offset, izup
#if defined MPI && defined ISO_C_BINDING
      integer :: prior_focus, ierr
      integer :: disp_unit = 1
      integer(kind=MPI_ADDRESS_KIND) :: win_size
      integer(c_intptr_t) :: cur_pos
      type(c_ptr) :: bptr, cptr
#endif
      real :: dum
      complex, dimension(:), allocatable :: phiext
      complex, dimension(:, :), allocatable :: gext
      logical :: debug = .false.
      character(100) :: message_dgdphi, message_QN, message_lu
      real, dimension(2) :: time_response_matrix_dgdphi
      real, dimension(2) :: time_response_matrix_QN
      real, dimension(2) :: time_response_matrix_lu

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
               nresponse = nelements(iky)
               win_size = win_size &
                          + int(nresponse, MPI_ADDRESS_KIND) * 4_MPI_ADDRESS_KIND &
                          + int(nresponse**2, MPI_ADDRESS_KIND) * 2 * real_size
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

         ! the response matrix for each ky has neigen(ky)
         ! independent sets of connected kx values
         if (.not. associated(response_matrix(iky)%eigen)) &
            allocate (response_matrix(iky)%eigen(1))

         ! treat zonal mode specially to avoid double counting
         ! as it is periodic
         nresponse = nelements(iky)

#if !defined ISO_C_BINDING || !defined MPI
         ! for each ky and set of connected kx values,
         ! must have a response matrix that is N x N
         ! with N = number of zeds per 2pi segment x number of 2pi segments
         if (.not. associated(response_matrix(iky)%eigen(1)%zloc)) &
            allocate (response_matrix(iky)%eigen(1)%zloc(nresponse, nresponse))

         ! response_matrix%idx is needed to keep track of permutations
         ! to the response matrix made during LU decomposition
         ! it will be input to LU back substitution during linear solve
         if (.not. associated(response_matrix(iky)%eigen(1)%idx)) &
            allocate (response_matrix(iky)%eigen(1)%idx(nresponse))
#else
         !exploit MPIs shared memory framework to reduce memory consumption of
         !response matrices

         if (.not. associated(response_matrix(iky)%eigen(1)%zloc)) then
            cptr = transfer(cur_pos, cptr)
            call c_f_pointer(cptr, response_matrix(iky)%eigen(1)%zloc, (/nresponse, nresponse/))
            cur_pos = cur_pos + nresponse**2 * 2 * nbytes_real
         end if

         if (.not. associated(response_matrix(iky)%eigen(1)%idx)) then
            cptr = transfer(cur_pos, cptr)
            call c_f_pointer(cptr, response_matrix(iky)%eigen(1)%idx, (/nresponse/))
            cur_pos = cur_pos + nresponse * 4
         end if
#endif
         response_matrix(iky)%eigen(1)%zloc = 0.0
         response_matrix(iky)%eigen(1)%idx = 0

         allocate (phiext(nresponse))
         ! loop over the sets of connected kx values
         do ie = 1, neigen(iky)

            if (proc0 .and. debug) then
               call time_message(.false., time_response_matrix_dgdphi, message_dgdphi)
            end if

            ! number of zeds x number of segments
            nz_ext = nsegments(ie, iky) * nzed_segment + 1

            allocate (gext(nz_ext, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
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
               call get_dgdphi_matrix_column(iky, ikx, iz, ie, idx, nz_ext, nresponse, phiext, gext)
            end do
            ! once we have used one segments, remaining segments
            ! have one fewer unique zed point
            izl_offset = 1
            if (nsegments(ie, iky) > 1) then
               do iseg = 2, nsegments(ie, iky)
                  ikx = ikxmod(iseg, ie, iky)
                  do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
                     idx = idx + 1
                     call get_dgdphi_matrix_column(iky, ikx, iz, ie, idx, nz_ext, nresponse, phiext, gext)
                  end do
                  if (izl_offset == 0) izl_offset = 1
               end do
            end if
            deallocate (gext)
         end do
         deallocate (phiext)
         !DSO - This ends parallelization over velocity space.
         !      At this point every processor has int dv dgdphi for a given ky
         !      and so the quasineutrality solve and LU decomposition can be
         !      parallelized locally if need be.
         !      This is preferable to parallelization over ky as the LU
         !      decomposition (and perhaps QN) will be dominated by the
         !      ky with the most connections

         if (proc0 .and. debug) then
            call time_message(.true., time_response_matrix_dgdphi, message_dgdphi)
            call time_message(.false., time_response_matrix_QN, message_QN)
         end if

#ifdef ISO_C_BINDING
         call mpi_win_fence(0, window, ierr)
#endif

         ! add quasineutrality operator, i.e. Gamma -  int dV (delta g / delta phi),
         ! where Gamma = int dV (1 - J_0^2 ) (Z^2 e T / n) F_0

#if defined ISO_C_BINDING && defined MPI
         if (sgproc0) then
#endif
            call get_fields_for_response_matrix(iky)
#if defined ISO_C_BINDING && defined MPI
         end if
#endif

#ifdef ISO_C_BINDING
         call mpi_win_fence(0, window, ierr)
#endif

         if (proc0 .and. debug) then
            call time_message(.true., time_response_matrix_QN, message_QN)
            call time_message(.false., time_response_matrix_lu, message_lu)
         end if

         !now we have the full response matrix. Finally, perform its LU decomposition
#ifdef MPI
         select case (lu_option_switch)
         case (lu_option_local)
#ifdef ISO_C_BINDING
            call parallel_LU_decomposition_local(iky)
#else
            call mp_abort('Stella must be built with HAS_ISO_BINDING in order to use local parallel LU decomposition.')
#endif
         case default
#endif
#ifdef ISO_C_BINDING
            if (sgproc0) then
#endif
               ! now that we have the reponse matrix for this ky and set of connected kx values
               !get the LU decomposition so we are ready to solve the linear system
               call lu_decomposition(response_matrix(iky)%eigen(1)%zloc, &
                                     response_matrix(iky)%eigen(1)%idx, dum)

               write (*,*) "CRENCH"

#ifdef ISO_C_BINDING
            end if
#endif
#ifdef MPI
         end select
#endif

         if (proc0 .and. debug) then
            call time_message(.true., time_response_matrix_lu, message_lu)
         end if

         time_response_matrix_dgdphi = 0
         time_response_matrix_QN = 0
         time_response_matrix_lu = 0

      end do

#ifdef ISO_C_BINDING
      call mpi_win_fence(0, window, ierr)
#endif

      if (proc0 .and. debug) then
         write (*, '(A)') "    ############################################################"
         write (*, '(A)') " "
      end if

   end subroutine init_response_matrix

   subroutine get_dgdphi_matrix_column(iky, ikx, iz, ie, idx, nz_ext, nresponse, phiext, gext)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use stella_time, only: code_dt
      use zgrid, only: delzed, nzgrid, ntubes
      use kt_grids, only: zonal_mode
      use extended_zgrid, only: periodic
      use full_xzgrid, only: xz_idx
      use species, only: spec
      use stella_geometry, only: gradpar, dbdzed
      use vpamu_grids, only: vpa, mu
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use fields_arrays, only: response_matrix
      use gyro_averages, only: aj0x
      use run_parameters, only: driftkinetic_implicit
      use run_parameters, only: maxwellian_inside_zed_derivative
      use parallel_streaming, only: stream_tridiagonal_solve
      use parallel_streaming, only: stream_sign
      use run_parameters, only: zed_upwind, time_upwind
#if defined ISO_C_BINDING && defined MPI
      use mp, only: sgproc0
#endif

      implicit none

      integer, intent(in) :: iky, ikx, iz, ie, idx, nz_ext, nresponse
      complex, dimension(:), intent(in out) :: phiext
      complex, dimension(:, vmu_lo%llim_proc:), intent(in out) :: gext

      integer :: ivmu, iv, imu, is, ia, it
      integer :: izp, izm, zm
      real :: mu_dbdzed_p, mu_dbdzed_m
      real :: fac, fac0, fac1, gyro_fac

      if(zonal_mode(iky).and.ikx.eq.1) return

      ia = 1

      if (.not. maxwellian_inside_zed_derivative) then
         ! get -vpa*b.gradz*Ze/T*F0*d<phi>/dz corresponding to unit impulse in phi
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            ! initialize g to zero everywhere along extended zed domain
            gext(:, ivmu) = 0.0
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            ! give unit impulse to phi at this zed location
            ! and compute -vpa*b.gradz*Ze/T*d<phi>/dz*F0 (RHS of streaming part of GKE)

            ! NB:  assuming equal spacing in zed below
            ! here, fac = -dt*(1+alph_t)/2*vpa*Ze/T*F0*J0/dz
            ! b.gradz left out because needs to be centred in zed
            if (driftkinetic_implicit) then
               gyro_fac = 1.0
            else
               gyro_fac = aj0x(iky, ikx, iz, ivmu)
            end if

            ! 0.125 to account for two linear interpolations
            fac = -0.125 * (1.+time_upwind) * code_dt * vpa(iv) * spec(is)%stm_psi0 &
                  * gyro_fac * spec(is)%zt / delzed(0) * maxwell_vpa(iv, is) * maxwell_fac(is)

            ! In the following, gradpar and maxwell_mu are interpolated separately
            ! to ensure consistency to what is done in parallel_streaming.f90

            ! stream_sign < 0 corresponds to positive advection speed
            if (stream_sign(iv) < 0) then
               if (iz > -nzgrid) then
                  ! fac0 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at this zed index
                  fac0 = fac * ((1.+zed_upwind) * gradpar(iz) &
                                + (1.-zed_upwind) * gradpar(iz - 1)) &
                         * ((1.+zed_upwind) * maxwell_mu(ia, iz, imu, is) &
                            + (1.-zed_upwind) * maxwell_mu(ia, iz - 1, imu, is))
                  ! fac1 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at the zed index to the right of
                  ! this one
                  if (iz < nzgrid) then
                     fac1 = fac * ((1.+zed_upwind) * gradpar(iz + 1) &
                                   + (1.-zed_upwind) * gradpar(iz)) &
                            * ((1.+zed_upwind) * maxwell_mu(ia, iz + 1, imu, is) &
                               + (1.-zed_upwind) * maxwell_mu(ia, iz, imu, is))
                  else
                     fac1 = fac * ((1.+zed_upwind) * gradpar(-nzgrid + 1) &
                                   + (1.-zed_upwind) * gradpar(nzgrid)) &
                            * ((1.+zed_upwind) * maxwell_mu(ia, -nzgrid + 1, imu, is) &
                               + (1.-zed_upwind) * maxwell_mu(ia, nzgrid, imu, is))
                  end if
               else
                  ! fac0 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at this zed index
                  fac0 = fac * ((1.+zed_upwind) * gradpar(iz) &
                                + (1.-zed_upwind) * gradpar(nzgrid - 1)) &
                         * ((1.+zed_upwind) * maxwell_mu(ia, iz, imu, is) &
                            + (1.-zed_upwind) * maxwell_mu(ia, nzgrid - 1, imu, is))
                  ! fac1 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at the zed index to the right of
                  ! this one
                  fac1 = fac * ((1.+zed_upwind) * gradpar(iz + 1) &
                                + (1.-zed_upwind) * gradpar(iz)) &
                         * ((1.+zed_upwind) * maxwell_mu(ia, iz + 1, imu, is) &
                            + (1.-zed_upwind) * maxwell_mu(ia, iz, imu, is))
               end if
               gext(idx, ivmu) = fac0
               if (idx < nz_ext) gext(idx + 1, ivmu) = -fac1
               ! zonal mode BC is periodic instead of zero, so must
               ! treat specially
               if (periodic(iky)) then
                  if (idx == 1) then
                     gext(nz_ext, ivmu) = fac0
                  else if (idx == nz_ext - 1) then
                     gext(1, ivmu) = -fac1
                  end if
               end if
            else
               if (iz < nzgrid) then
                  ! fac0 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at this zed index
                  fac0 = fac * ((1.+zed_upwind) * gradpar(iz) &
                                + (1.-zed_upwind) * gradpar(iz + 1)) &
                         * ((1.+zed_upwind) * maxwell_mu(ia, iz, imu, is) &
                            + (1.-zed_upwind) * maxwell_mu(ia, iz + 1, imu, is))
                  ! fac1 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at the zed index to the left of
                  ! this one
                  if (iz > -nzgrid) then
                     fac1 = fac * ((1.+zed_upwind) * gradpar(iz - 1) &
                                   + (1.-zed_upwind) * gradpar(iz)) &
                            * ((1.+zed_upwind) * maxwell_mu(ia, iz - 1, imu, is) &
                               + (1.-zed_upwind) * maxwell_mu(ia, iz, imu, is))
                  else
                     fac1 = fac * ((1.+zed_upwind) * gradpar(nzgrid - 1) &
                                   + (1.-zed_upwind) * gradpar(iz)) &
                            * ((1.+zed_upwind) * maxwell_mu(ia, nzgrid - 1, imu, is) &
                               + (1.-zed_upwind) * maxwell_mu(ia, iz, imu, is))
                  end if
               else
                  ! fac0 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at this zed index
                  fac0 = fac * ((1.+zed_upwind) * gradpar(iz) &
                                + (1.-zed_upwind) * gradpar(-nzgrid + 1)) &
                         * ((1.+zed_upwind) * maxwell_mu(ia, iz, imu, is) &
                            + (1.-zed_upwind) * maxwell_mu(ia, -nzgrid + 1, imu, is))
                  ! fac1 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at the zed index to the left of
                  ! this one
                  fac1 = fac * ((1.+zed_upwind) * gradpar(iz - 1) &
                                + (1.-zed_upwind) * gradpar(iz)) &
                         * ((1.+zed_upwind) * maxwell_mu(ia, iz - 1, imu, is) &
                            + (1.-zed_upwind) * maxwell_mu(ia, iz, imu, is))
               end if
               gext(idx, ivmu) = -fac0
               if (idx > 1) gext(idx - 1, ivmu) = fac1
               ! zonal mode BC is periodic instead of zero, so must
               ! treat specially
               if (periodic(iky)) then
                  if (idx == 1) then
                     gext(nz_ext, ivmu) = -fac0
                     gext(nz_ext - 1, ivmu) = fac1
                  else if (idx == 2) then
                     gext(nz_ext, ivmu) = fac1
                  end if
               end if
            end if

            ! hack for now (duplicates much of the effort from sweep_zed_zonal)
            if (periodic(iky)) then
               call sweep_zed_zonal_response(iv, is, stream_sign(iv), gext(:, ivmu))
            else
               ! invert parallel streaming equation to get g^{n+1} on extended zed grid
               ! (I + (1+alph)/2*dt*vpa)*g_{inh}^{n+1} = RHS = gext
               call stream_tridiagonal_solve(iky, ie, iv, is, gext(:, ivmu))
            end if

         end do
      else
         ! get -vpa*b.gradz*Ze/T*F0*d<phi>/dz corresponding to unit impulse in phi
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            ! initialize g to zero everywhere along extended zed domain
            gext(:, ivmu) = 0.0
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            ! give unit impulse to phi at this zed location
            ! and compute -vpa*b.gradz*Ze/T*d<phi>/dz*F0 (RHS of streaming part of GKE)

            ! NB:  assuming equal spacing in zed below
            ! here, fac = -dt*(1+alph_t)/2*vpa*Ze/T*F0*J0/dz
            ! b.gradz left out because needs to be centred in zed
            if (driftkinetic_implicit) then
               gyro_fac = 1.0
            else
               gyro_fac = aj0x(iky, ikx, iz, ivmu)
            end if

            fac = -0.25 * (1.+time_upwind) * code_dt * vpa(iv) * spec(is)%stm_psi0 &
                  * gyro_fac * spec(is)%zt * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)

            mu_dbdzed_p = 1./delzed(0) + mu(imu) * dbdzed(ia, iz) * (1.+zed_upwind)
            mu_dbdzed_m = 1./delzed(0) + mu(imu) * dbdzed(ia, iz) * (1.-zed_upwind)

            ! stream_sign < 0 corresponds to positive advection speed
            if (stream_sign(iv) < 0) then
               if (iz > -nzgrid) then
                  ! fac0 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at this zed index
                  fac0 = fac * ((1.+zed_upwind) * gradpar(iz) &
                                + (1.-zed_upwind) * gradpar(iz - 1)) * mu_dbdzed_p
                  ! fac1 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at the zed index to the right of
                  ! this one
                  if (iz < nzgrid) then
                     izp = iz + 1
                  else
                     izp = -nzgrid + 1
                  end if
                  fac1 = fac * ((1.+zed_upwind) * gradpar(izp) &
                                + (1.-zed_upwind) * gradpar(iz)) * mu_dbdzed_m
               else
                  ! fac0 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at this zed index
                  fac0 = fac * ((1.+zed_upwind) * gradpar(iz) &
                                + (1.-zed_upwind) * gradpar(nzgrid - 1)) * mu_dbdzed_p
                  ! fac1 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at the zed index to the right of
                  ! this one
                  fac1 = fac * ((1.+zed_upwind) * gradpar(iz + 1) &
                                + (1.-zed_upwind) * gradpar(iz)) * mu_dbdzed_m
               end if
               gext(idx, ivmu) = fac0
               if (idx < nz_ext) gext(idx + 1, ivmu) = -fac1
               ! zonal mode BC is periodic instead of zero, so must
               ! treat specially
               if (periodic(iky)) then
                  if (idx == 1) then
                     gext(nz_ext, ivmu) = fac0
                  else if (idx == nz_ext - 1) then
                     gext(1, ivmu) = -fac1
                  end if
               end if
            else
               if (iz < nzgrid) then
                  ! fac0 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at this zed index
                  fac0 = fac * ((1.+zed_upwind) * gradpar(iz) &
                                + (1.-zed_upwind) * gradpar(iz + 1)) * mu_dbdzed_p
                  ! fac1 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at the zed index to the left of
                  ! this one
                  if (iz > -nzgrid) then
                     izm = iz - 1
                  else
                     izm = nzgrid - 1
                  end if
                  fac1 = fac * ((1.+zed_upwind) * gradpar(izm) &
                                + (1.-zed_upwind) * gradpar(iz)) * mu_dbdzed_m
               else
                  ! fac0 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at this zed index
                  fac0 = fac * ((1.+zed_upwind) * gradpar(iz) &
                                + (1.-zed_upwind) * gradpar(-nzgrid + 1)) * mu_dbdzed_p
                  ! fac1 is the factor multiplying delphi on the RHS
                  ! of the homogeneous GKE at the zed index to the left of
                  ! this one
                  fac1 = fac * ((1.+zed_upwind) * gradpar(iz - 1) &
                                + (1.-zed_upwind) * gradpar(iz)) * mu_dbdzed_m
               end if
               gext(idx, ivmu) = -fac0
               if (idx > 1) gext(idx - 1, ivmu) = fac1
               ! zonal mode BC is periodic instead of zero, so must
               ! treat specially
               if (periodic(iky)) then
                  if (idx == 1) then
                     gext(nz_ext, ivmu) = -fac0
                     gext(nz_ext - 1, ivmu) = fac1
                  else if (idx == 2) then
                     gext(nz_ext, ivmu) = fac1
                  end if
               end if
            end if

            ! hack for now (duplicates much of the effort from sweep_zed_zonal)
            if (periodic(iky)) then
               call sweep_zed_zonal_response(iv, is, stream_sign(iv), gext(:, ivmu))
            else
               ! invert parallel streaming equation to get g^{n+1} on extended zed grid
               ! (I + (1+alph)/2*dt*vpa)*g_{inh}^{n+1} = RHS = gext
               call stream_tridiagonal_solve(iky, ie, iv, is, gext(:, ivmu))
            end if

         end do
      end if

      ! we now have g on the extended zed domain at this ky and set of connected kx values
      ! corresponding to a unit impulse in phi at this location
      ! now integrate over velocities to get a square response matrix
      ! (this ends the parallelization over velocity space, so every core should have a
      !  copy of phiext)
      zm = 0
      if (zonal_mode(iky)) zm = 1
      do it = 1, ntubes
         call integrate_over_velocity_radial(gext, phiext, it, iky, ie)

#if !defined ISO_C_BINDING || !defined MPI
         response_matrix(iky)%eigen(1)%zloc(:, xz_idx(ikx, iz, it, zm) = -phiext(:nresponse)
#else
         if (sgproc0) response_matrix(iky)%eigen(1)%zloc(:, xz_idx(ikx, iz, it, zm)) = -phiext(:nresponse)
#endif
      end do

   end subroutine get_dgdphi_matrix_column

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

   subroutine integrate_over_velocity_radial(g, phi, it_in, iky, ie)

      use dist_fn_arrays, only: kperp2, dkperp2dr
      use stella_layouts, only: vmu_lo, imu_idx, is_idx
      use stella_geometry, only: bmag, dBdrho
      use species, only: nspec, spec
      use kt_grids, only: nakx, rho_d_clamped
      use zgrid, only: nzgrid, ntubes
      use extended_zgrid, only: map_from_extended_zgrid
      use full_xzgrid, only: map_to_full_xzgrid
      use vpamu_grids, only: integrate_species, vperp2
      use gyro_averages, only: aj0x, aj1x
      use mp, only: sum_allreduce
      use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use physics_flags, only: radial_variation

      implicit none

      complex, dimension(:, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:), intent(out) :: phi
      integer, intent(in) :: iky, ie, it_in

      integer :: ikx, iz, ia, it, ivmu, imu, is
      real, dimension(nspec) :: wgt
      complex, dimension(:, :), allocatable :: g0k, g0x
      complex, dimension(:, :, :), allocatable :: g0, g1
      complex, dimension(:, :, :, :), allocatable :: gyro_0, gyro_1

      ia = 1

      allocate (g0k(1, nakx))
      allocate (g0x(1, nakx))

      allocate (gyro_0(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      allocate (gyro_1(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      allocate (g0(nakx, -nzgrid:nzgrid, ntubes))
      allocate (g1(nakx, -nzgrid:nzgrid, ntubes))

      wgt = spec%z * spec%dens_psi0

      gyro_0 = 0.0

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         call map_from_extended_zgrid(it_in, ie, iky, g(:, ivmu), gyro_0(:, :, :, ivmu))
      end do

      gyro_0 = gyro_0 * spread(aj0x(iky, :, :, :), 3, ntubes)
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               call integrate_species(gyro_0(ikx, iz, it, :), iz, wgt, g0(ikx, iz, it), reduce_in=.false.)
            end do
         end do
      end do
      call sum_allreduce(g0)

      g1 = 0.0
      if(radial_variation) then
         gyro_1 = 0.0
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            call map_from_extended_zgrid(it_in, ie, iky, g(:, ivmu), gyro_1(:, :, :, ivmu))
         end do

         gyro_1 = gyro_1 * spread(aj0x(iky, :, :, :), 3, ntubes)

         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  gyro_1(:, iz, it, ivmu) = gyro_1(:, iz, it, ivmu) &
                                         * (-0.5 * aj1x(iky, :, iz, ivmu) / aj0x(iky, :, iz, ivmu) * (spec(is)%smz)**2 &
                                            * (kperp2(iky, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                            * (dkperp2dr(iky, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                            + dBdrho(iz) / bmag(ia, iz))
               end do
            end do
         end do

         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  call integrate_species(gyro_1(ikx, iz, it, :), iz, wgt, g1(ikx, iz, it), reduce_in=.false.)
               end do
            end do
         end do

         call sum_allreduce(g1)

         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0k(1, :) = g1(:, iz, it)
               call transform_kx2x_unpadded(g0k, g0x)
               g0x(1, :) = rho_d_clamped * g0x(1, :)
               call transform_x2kx_unpadded(g0x, g0k)
               g1(:, iz, it) = g0k(1, :)
            end do
         end do
      endif
      g0 = g0 + g1

      call map_to_full_xzgrid(iky, g0, phi)

      deallocate (g0, g1, gyro_0, gyro_1, g0k, g0x)

   end subroutine integrate_over_velocity_radial

   subroutine get_fields_for_response_matrix(iky)

      use fields_arrays, only: response_matrix
      use fields_arrays, only: gamtot, dgamtotdr
      use kt_grids, only: nakx, rho_d_clamped, zonal_mode
      use zgrid, only: nzgrid, ntubes
      use full_xzgrid, only: xz_idx
      use extended_zgrid, only: periodic
      use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use physics_flags, only: radial_variation
!     use species, only: spec
!     use species, only: has_electron_species
!     use stella_geometry, only: dl_over_b
!     use physics_flags, only: adiabatic_option_switch
!     use physics_flags, only: adiabatic_option_fieldlineavg

      implicit none

      integer, intent(in) :: iky

      integer ::ikx, ix, jx, iz, it, ia, zm, pm, idx1, idx2
      complex, dimension(:, :), allocatable :: gamma_mat, g0k, g0x
!     complex :: tmp

      ia = 1
      pm = 0 ; zm = 0
      if (periodic(iky)) pm = 1
      if (zonal_mode(iky)) zm = 1

      allocate (gamma_mat(nakx - zm, nakx - zm))
      allocate (g0k(1, nakx))
      allocate (g0x(1, nakx))

      do iz = -nzgrid, nzgrid - pm
         gamma_mat = 0.0
         do ikx = 1 + zm, nakx
            if(radial_variation) then
               g0k(1, :) = 0.0
               g0k(1, ikx) = dgamtotdr(iky, ikx, iz)

               call transform_kx2x_unpadded(g0k, g0x)
               g0x(1, :) = rho_d_clamped * g0x(1, :)
               call transform_x2kx_unpadded(g0x, g0k)

               !row column
               gamma_mat(:, ikx - zm) = g0k(1, (1 + zm):)
            endif
            gamma_mat(ikx - zm, ikx - zm) = gamma_mat(ikx - zm, ikx - zm) + gamtot(iky, ikx, iz)
         end do
         do it = 1, ntubes
            do ix = 1 + zm, nakx
               do jx = 1 + zm, nakx
                  idx1 = xz_idx(ix, iz, it, zm)
                  idx2 = xz_idx(jx, iz, it, zm)
                  response_matrix(iky)%eigen(1)%zloc(idx1, idx2) = response_matrix(iky)%eigen(1)%zloc(idx1, idx2) &
                                                                   + gamma_mat(ix - zm, jx - zm)
               end do
            end do
         end do
!        call add_submatrix_to_full_matrix(iky, iz, gamma_mat, response_matrix(iky)%eigen(1)%zloc)
      end do

!     if (.not. has_electron_species(spec) .and. &
!         adiabatic_option_switch == adiabatic_option_fieldlineavg) then
!        if (zonal_mode(iky)) then
!           ! no connections for ky = 0
!           iseg = 1
!           tmp = sum(dl_over_b(ia, :) * phi)
!           phi = phi + tmp * gamtot3(ikxmod(1, ie, iky), :)
!        end if
!     end if

      deallocate (gamma_mat, g0k, g0x)

   end subroutine get_fields_for_response_matrix

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
      use mpi

      implicit none

      integer, intent(in) :: iky

      integer, dimension(:, :), allocatable :: eig_limits
      integer, dimension(:), allocatable :: job_list
      integer, dimension(:), allocatable :: row_limits
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

         if (iproc == jroot) neig = 1

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

      !copy all the matrices across all nodes
      if (sgproc0) then
         do ie = 1, 1
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
#endif /* MPI */

end module response_matrix_radial
