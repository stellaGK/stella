!###############################################################################
!############################### RESPONSE MATRIX ###############################
!###############################################################################
! 
! This module computes the response matrix for inverting the parallel streaming
! operator when parallel_streaming is treated implicitly.
! 
!###############################################################################
module response_matrix

   use netcdf
   use mpi
#ifdef ISO_C_BINDING
   use, intrinsic :: iso_c_binding, only: c_intptr_t
#endif

   ! Load debug flags
   use debug_flags, only: debug => response_matrix_debug

   implicit none

   ! Make routines available to other modules
   public :: init_response_matrix, finish_response_matrix
   public :: read_response_matrix
   public :: initialised_response_matrix

   private

#ifdef ISO_C_BINDING
   integer(c_intptr_t) :: cur_pos
#endif
   character(100) :: message_dgdphi, message_QN, message_lu
   real, dimension(2) :: time_dgdphi
   real, dimension(2) :: time_QN
   real, dimension(2) :: time_lu
   
   ! Only initialise once
   logical :: initialised_response_matrix = .false.

contains

!###############################################################################
!########################## INITIALISE RESPONSE MATRIX #########################
!###############################################################################

   !****************************************************************************
   !            Main Routine Called to Initialise the Response Matrix           
   !****************************************************************************
   subroutine init_response_matrix

      use linear_solve, only: lu_decomposition
      use arrays, only: response_matrix
      use parallelisation_layouts, only: iv_idx, is_idx
      use grids_kxky, only: naky
      use mp, only: proc0
      use parallelisation_layouts, only: mat_gen
      use file_units, only: unit_response_matrix
#ifdef ISO_C_BINDING
      use arrays, only: response_window
#endif

      implicit none

      ! Local variables
#ifdef ISO_C_BINDING
      integer :: ierr
#endif

      !-------------------------------------------------------------------------

      ! Debug message -> print to terminal
      if (debug) call write_response_matrix_message (1)
      ! Set up response matrix utils
      call setup_response_matrix_timings
      call setup_response_matrix_file_io

      ! Only initialise once
      if (initialised_response_matrix) return
      initialised_response_matrix = .true.

      ! Allocate response matrix
      if (.not. allocated(response_matrix)) allocate (response_matrix(naky))

#ifdef ISO_C_BINDING
      call setup_shared_memory_window
#endif

      ! Construct the response matrix
      call construct_response_matrix

#ifdef ISO_C_BINDING
      call mpi_win_fence(0, response_window, ierr)
#endif

      ! Write the response matrix to an output file
      if (proc0 .and. mat_gen) then
         close (unit=unit_response_matrix)
      end if

      ! End debug message -> print to terminal
      if (debug) call write_response_matrix_message (2)

   end subroutine init_response_matrix

   !****************************************************************************
   !                         Write message to terminal
   !****************************************************************************
   subroutine write_response_matrix_message (toggle)

      integer, intent (in) :: toggle

      if (toggle == 1) then
         write (*, *) " "
         write (*, '(A)') "    ############################################################"
         write (*, '(A)') "                         RESPONSE MATRIX"
         write (*, '(A)') "    ############################################################"
      elseif (toggle == 2) then
         write (*, '(A)') "    ############################################################"
         write (*, '(A)') " "
      end if

   end subroutine write_response_matrix_message

   !****************************************************************************
   !                 Write timings for dealing with response matrix
   !****************************************************************************
   subroutine setup_response_matrix_timings

      implicit none

      message_dgdphi = '     calculate dgdphi: '
      message_QN = '     calculate QN:     '
      message_lu = '     calculate LU:     '
      time_dgdphi = 0.0
      time_QN = 0.0
      time_lu = 0.0

   end subroutine setup_response_matrix_timings

   !****************************************************************************
   !                     Set up .io file for response matrix                    
   !****************************************************************************
   subroutine setup_response_matrix_file_io

      use mp, only: proc0, job
      use parallelisation_layouts, only: mat_gen
      use system_fortran, only: systemf
      use grids_kxky, only: naky
      use file_units, only: unit_response_matrix

      implicit none

      ! Local variables
      character(len=15) :: job_str
      character(len=100) :: file_name

      !-------------------------------------------------------------------------

      ! All matrices handled by processor i_proc and job are stored
      ! on a single file named: response_mat_job.iproc
      if (proc0 .and. mat_gen) then
         call systemf('mkdir -p mat')
         write (job_str, '(I1.1)') job
         file_name = './mat/response_mat_'//trim(job_str)
         open (unit=unit_response_matrix, status='replace', file=file_name, &
               position='rewind', action='write', form='unformatted')
         write (unit=unit_response_matrix) naky
      end if

   end subroutine setup_response_matrix_file_io

   !****************************************************************************
   !                          Set up shared memory window
   !****************************************************************************
   ! Create a single shared memory window for all the response matrices and
   ! permutation arrays. Creating a window for each matrix/array would lead 
   ! to performance degradation on some clusters
   !****************************************************************************
#ifdef ISO_C_BINDING
   subroutine setup_shared_memory_window

      use mpi
      use, intrinsic :: iso_c_binding, only: c_intptr_t
      use mp, only: sgproc0, real_size
      use mp, only: create_shared_memory_window
      use arrays, only: response_window
      use field_equations, only: nfields
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen, nsegments, nzed_segment
      use grids_extended_zgrid, only: periodic

      implicit none

      ! Local variables
      integer(kind=MPI_ADDRESS_KIND) :: win_size
      integer :: iky, ie
      integer :: nresponse

      !-------------------------------------------------------------------------

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
                  win_size = win_size + int(nresponse, MPI_ADDRESS_KIND) * 4_MPI_ADDRESS_KIND &
                      + int(nresponse**2, MPI_ADDRESS_KIND) * 2 * real_size
               end do
            end do
         end if

         call create_shared_memory_window(win_size, response_window, cur_pos)
      end if

   end subroutine setup_shared_memory_window
#endif


!###############################################################################
!########################## CONSTRUCT RESPONSE MATRIX ##########################
!###############################################################################

   !****************************************************************************
   !               Main routine for constructing response martrix
   !****************************************************************************
   ! Background:
   ! -----------
   ! Different ky modes are independent, so this creates sets of connected kx values:
   !
   ! After one full 2π orbit in zed, the eddy becomes sheared, shifting the mode 
   ! to a higher kx value. To represent this effect, we apply the 
   ! 'twist-and-shift' boundary conditions, which connect these shifted kx values 
   ! along an extended z-grid. For a given ky, the connected kx modes are spaced 
   ! by
   !                 δkx = 2π p ι' * (∂y/∂α) * (∂ψ/∂x) * ky
   ! where ι' is the radial derivative of the rotational transform. 
   !
   ! This construction implies that not all kx modes are mutually connected. 
   ! Instead, for each ky we form distinct "chains" of connected modes, with the 
   ! number of such chains determined by the spacing δkx. 
   !
   ! <neigen> gives the number of distinct chains for a given ky. If we use periodic 
   ! boundary conditions then all modes are connected and <neigen> = nakx
   !
   ! What is done:
   ! ------------- 
   ! In this subroutine we construct the response matrix. This is split into three 
   ! main contributions: 
   !     1) Getting the velocity-space-integrated response. To do this we supply    
   !        a unit impulse at every point in the extended zed domain for a given ky,
   !        and set of connected kx values. With this, we then solve the parallel
   !        streaming equations to find the response of the distribution function
   !        to this unit impulse. This is done for each field that is included 
   !        (phi, apar, bpar). Once this distribution function is obtained
   !        we integrate it with appropriate weights over (vpa, mu) and sum over
   !        species -> This gives the quasineutrality component that acts on g.
   !     2) Apply the field solve to find the matrix that we need to invert. This
   !        is essentially dividing by the correct factor from the field equations. 
   !     3) LU decompose the matrix.
   !****************************************************************************
   subroutine construct_response_matrix

      use mp, only: proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: mat_gen
      use arrays, only: response_matrix
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen
      use file_units, only: unit_response_matrix
#ifdef ISO_C_BINDING
      use arrays, only: response_window
#endif

      implicit none

      ! Local variables
      integer :: iky, ie
#ifdef ISO_C_BINDING
      integer :: ierr
#endif

      !-------------------------------------------------------------------------

      do iky = 1, naky

         if (proc0 .and. mat_gen) then
            write (unit=unit_response_matrix) iky, neigen(iky)
         end if

         ! ---------------------------------------------------------------------
         !                STEP 1) Find distribution function response      
         ! ---------------------------------------------------------------------
         ! For a given ky, we need to associate an %eigen to it - to denote the connected
         ! modes. The response matrix for each ky has neigen(ky).
         if (.not. associated(response_matrix(iky)%eigen)) &
            allocate (response_matrix(iky)%eigen(neigen(iky)))

         ! Write time message
         if (debug) call time_message(.false., time_dgdphi, message_dgdphi)

         call calculate_vspace_integrated_response(iky)

         ! ---------------------------------------------------------------------
         !                   Parallelisation + debug messages
         ! ---------------------------------------------------------------------
         ! This ends parallelisation over velocity space.
         ! At this point every processor has int dv dgdphi for a given ky
         ! and so the quasineutrality solve and LU decomposition can be
         ! parallelised locally if need be.
         ! This is preferable to parallelisation over ky as the LU
         ! decomposition (and perhaps QN) will be dominated by the
         ! ky with the most connections -- (DSO)

         ! Debug time messages
         if (debug) then
            call time_message(.true., time_dgdphi, message_dgdphi)
            call time_message(.false., time_QN, message_QN)
         end if

#ifdef ISO_C_BINDING
         call mpi_win_fence(0, response_window, ierr)
#endif

         ! ---------------------------------------------------------------------
         !                     STEP 2) Apply field solve
         ! ---------------------------------------------------------------------

         call apply_field_solve_to_finish_response_matrix(iky)

         ! ---------------------------------------------------------------------
         !                   Parallelisation + debug messages
         ! ---------------------------------------------------------------------

#ifdef ISO_C_BINDING
         call mpi_win_fence(0, response_window, ierr)
#endif

         ! Debug time messages
         if (debug) then
            call time_message(.true., time_QN, message_QN)
            call time_message(.false., time_lu, message_lu)
         end if

         ! ---------------------------------------------------------------------
         !                   STEP 3) LU decompose the matrix 
         ! ---------------------------------------------------------------------

         call lu_decompose_response_matrix(iky)

         ! ---------------------------------------------------------------------
         !             Finish with correct debug messages and write to files
         ! ---------------------------------------------------------------------
         if (proc0 .and. debug) then
            call time_message(.true., time_lu, message_lu)
         end if

         time_dgdphi = 0
         time_QN = 0
         time_lu = 0

         ! Write response matrix to file 
         do ie = 1, neigen(iky)
            if (proc0 .and. mat_gen) then
               write (unit=unit_response_matrix) response_matrix(iky)%eigen(ie)%idx
               write (unit=unit_response_matrix) response_matrix(iky)%eigen(ie)%zloc
            end if
         end do

      end do

   end subroutine construct_response_matrix

   !============================================================================
   !****************************************************************************
   !                Step 1) Find distribution function response                      
   !****************************************************************************
   !============================================================================
   ! The next routines are all dedicated to step 1).
   !
   !
   ! Within the routines get_dpdf_d * _matrix_column (where * can be phi, apar, 
   ! or bpar) the following steps take place:
   !
   ! 1.A) First initialise unit impulse for each field at each grid point in the 
   !      extended zed domain for a given connected chain.
   ! 1.B) Then solve the inhomogeneous parallel streaming equation for the distribution
   !      function given a field with a unit impulse.
   ! 1.C) Then integrate over velocities to get the appropriate field operator that acts
   !      on the distibution function:
   !           For phi:     2B/sqrt(π) int dvpa int dmu J_0 * g
   !           For apar:    β sum_s Z_s n_s vth 2*B0/sqrt{π} \int d^2v vpar J_0 g
   !           For bpar:    - 2β sum_s n_s T_s 2*B0/sqrt{π} \int d^2v mu J_1/a_s g 
   !============================================================================
   subroutine calculate_vspace_integrated_response(iky)

      use mp, only: proc0
      use parallelisation_layouts, only: mat_gen
      use parameters_physics, only: include_apar, include_bpar
      use grids_extended_zgrid, only: neigen, ikxmod
      use grids_extended_zgrid, only: nsegments, nzed_segment
      use grids_extended_zgrid, only: periodic
      use grids_extended_zgrid, only: iz_low, iz_up
      use parallelisation_layouts, only: vmu_lo
      use field_equations, only: nfields
      use file_units, only: unit_response_matrix

      implicit none

      ! Arguments
      integer, intent(in) :: iky

      ! Local variables
      integer :: ie, idx, ikx, iseg
      integer :: iz, izl_offset, izup
      integer :: nz_ext, nresponse, nresponse_per_field
      complex, dimension(:, :), allocatable :: gext
      complex, dimension(:), allocatable :: phi_ext, apar_ext, bpar_ext

      !-------------------------------------------------------------------------

      ! Loop over the independent chains for a given ky value  <neigen> is an
      ! integer that depends on ky, as different ky modes will have a different
      ! number of eigen chains.
      do ie = 1, neigen(iky)
         !----------------------------------------------------------------------
         !              Set up system to compute response of the pdf
         !----------------------------------------------------------------------
         ! <nz_ext> is the dimension of the extended zed domain, i.e.,
         ! <nz_ext> = (number of zeds) x (number of segments on extended zed domain)
         ! <nzed_segment> = (number of zeds)
         !                -> This is the number of unique zed values in all segments 
         !                   except the first. The first segment has one extra unique 
         !                   zed value (all others have one grid common with the
         !                   previous segment due to periodicity)
         ! <nsegments>    = (number of segments on extended zed domain)
         !                = (Nkx -1)/neigen 
         !                -> <nsegments> is how many 2π segments in z are linked together
         nz_ext = nsegments(ie, iky) * nzed_segment + 1

         ! We need to treat zonal mode specially to avoid double counting the end
         ! points in zed, as it is periodic so these end points are the same.
         ! This is also needed if 'periodic' is chosen for the boundary conditions.
         ! At this point <nresponse> is the number of independednt zed values for
         ! each field. This is the number of unit impulses we need to apply in order
         ! to get the response matrix.
         if (periodic(iky)) then
            nresponse_per_field = nz_ext - 1
         else
            nresponse_per_field = nz_ext
         end if

         ! If electromagnetic, we need to consider the response of phi, apar, and bpar. 
         ! If we have more fields we need to apply more unit impulses -> for each of
         ! phi, apar, and bpar.
         nresponse = nresponse_per_field * nfields

         ! Write <ie>, <nresponse> to file
         if (proc0 .and. mat_gen) then
            write (unit=unit_response_matrix) ie, nresponse
         end if

         ! Allocate response_matrix%eigen%zloc - size of response matrix to invert
         !          response_matrix%idx        - pivot index needed for LU decomposition
         call setup_response_matrix_zloc_idx(iky, ie, nresponse)

         ! Allocate arrays on the extended zed domain
         ! Fields only have 1 index, as we are looping over ky modes, and kx and zed are connected
         ! on via the extended zed domain
         ! The distribution function, gext, has 2 dimensions - the extended zed dimension, and
         ! velocity.
         allocate (phi_ext(nz_ext)); phi_ext = 0.0 
         allocate (apar_ext(nz_ext)); apar_ext = 0.0
         allocate (bpar_ext(nz_ext)); bpar_ext = 0.0
         allocate (gext(nz_ext, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); gext = 0.0

         !----------------------------------------------------------------------
         !              Get Response of the pdf to unit impulses
         !----------------------------------------------------------------------
         ! Steps 1.A), 1.B) and 1.C) occur below
         !----------------------------------------------------------------------

         ! <idx> here is the index in the extended zed domain that we are giving
         ! a unit impulse to.
         idx = 0

         ! Loop over segments, and idex this with <iseg>. 
         ! The first segment is special because it has one more unique
         ! zed value than all the other segments since the domain
         ! is [z0-pi:z0+pi], and we are including both endpoints.
         ! For iseg > 1 one endpoint is shared with the previous segment.
         iseg = 1
         ! ikxmod gives the kx corresponding to (iseg,ie,iky)
         ! i.e. given a ky (iky), which chain are we in (ie), and within that chain
         !      which segment of 2π are we considering -> this is associated with a specific
         !      kx value. 
         ikx = ikxmod(iseg, ie, iky)
         izl_offset = 0

         ! Avoid double-counting of periodic points for zonal mode (and other periodic modes)
         ! Here, define <izup> as the upper zed value within a segment, If the mode is periodic
         ! this will be one less, as it has a repeated point.
         if (periodic(iky)) then
            izup = iz_up(iseg) - 1
         else
            izup = iz_up(iseg)
         end if

         ! Here, we apply a unit impulse at every value of zed in this sements, and find the 
         ! response of gext to this unit impulse. The impulse is provided at the <idx> value.
         ! Note - No need to obtain response to impulses at negative kx values
         do iz = iz_low(iseg), izup
            idx = idx + 1
            call get_dpdf_dphi_matrix_column(iky, ie, idx, nz_ext, nresponse_per_field, phi_ext, apar_ext, bpar_ext, gext)
            if (include_apar) call get_dpdf_dapar_matrix_column(iky, ie, idx, nz_ext, nresponse_per_field, phi_ext, apar_ext, bpar_ext, gext)
            if (include_bpar) call get_dpdf_dbpar_matrix_column(iky, ie, idx, nz_ext, nresponse_per_field, phi_ext, apar_ext, bpar_ext, gext)
         end do

         ! Once we have used one segment, remaining segments have one fewer unique zed point
         ! Note that these routines are identical to the lines above, and do the same thing, but are 
         ! simply for the connected segments.
         izl_offset = 1
         if (nsegments(ie, iky) > 1) then
            ! Loop over remaining segments
            do iseg = 2, nsegments(ie, iky)
               ! Compute the index of kx that is connected in the given eigen chain in this segment 
               ikx = ikxmod(iseg, ie, iky)
               ! Now apply a unit impulse at each zed location. To do this loop over the zed index, but
               ! recall that there is one less zed grid point in these connected segments as they share
               ! a grid point with the previous segment. 
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
         !----------------------------------------------------------------------
      end do

   end subroutine calculate_vspace_integrated_response

   !****************************************************************************
   !                            Allocate zloc and idx 
   !****************************************************************************
   ! Sets up the response matrix storage for a given (iky, ie) eigenmode.  
   !
   ! Allocate the following: 
   ! ----------------------
   ! - response_matrix%eigen%zloc is the dimension of the response matrix for a 
   !   given eigen chain 
   ! - response_matrix%idx is needed to keep track of permutations to the response
   !   matrix made during LU decomposition it will be input to LU back substitution
   !   during linear solve
   !
   ! Memory:
   ! -------
   ! - When ISO_C_BINDING is available, it uses MPI’s shared memory buffer
   !   to avoid redundant allocations across ranks. This is done by mapping
   !   existing memory into Fortran pointers with c_f_pointer.  
   ! - Otherwise, it allocates fresh Fortran arrays.  
   !   cur_pos acts as a memory cursor that steps through a shared memory block, 
   !   ensuring each zloc and idx block points to its own reserved region.
   !****************************************************************************
   subroutine setup_response_matrix_zloc_idx(iky, ie, nresponse)

#ifdef ISO_C_BINDING
      use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
      use mp, only: nbytes_real
#endif
      use arrays, only: response_matrix

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie, nresponse

#ifdef ISO_C_BINDING
      type(c_ptr) :: cptr

      !-------------------------------------------------------------------------
      ! Exploit MPIs shared memory framework to reduce memory consumption 
      ! of the response matrices by mapping existing memory blocks instead 
      ! of allocating separately on each process.
      !-------------------------------------------------------------------------

      ! Associate zloc (a 2D response matrix) if not already connected
      if (.not. associated(response_matrix(iky)%eigen(ie)%zloc)) then
         ! Convert current memory cursor (cur_pos) into a C pointer
         cptr = transfer(cur_pos, cptr)

         ! Map the C pointer into Fortran as a 2D array (nresponse x nresponse)
         call c_f_pointer(cptr, response_matrix(iky)%eigen(ie)%zloc, (/nresponse, nresponse/))
         
         ! Advance cursor: nresponse^2 elements, each complex (2 reals)
         cur_pos = cur_pos + nresponse**2 * 2 * nbytes_real
      end if

      ! Associate idx (pivot indices) if not already connected
      if (.not. associated(response_matrix(iky)%eigen(ie)%idx)) then
         ! Convert current memory cursor to a C pointer
         cptr = transfer(cur_pos, cptr)

         ! Map it into Fortran as a 1D integer array of length nresponse
         call c_f_pointer(cptr, response_matrix(iky)%eigen(ie)%idx, (/nresponse/))

         ! Advance cursor: nresponse integers, each assumed to take 4 bytes
         cur_pos = cur_pos + nresponse * 4
      end if
#else
      !-------------------------------------------------------------------------
      ! If ISO_C_BINDING is not available:
      ! Allocate arrays normally in Fortran
      !-------------------------------------------------------------------------
      ! For each ky and set of connected kx values, so we must have a response
      ! matrix that is N x N , with N = number of zeds per 2pi segment x number 
      ! of 2pi segments

      ! Allocate zloc as an (nresponse x nresponse) matrix if not already allocated
      if (.not. associated(response_matrix(iky)%eigen(ie)%zloc)) &
         allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse, nresponse))

      ! Allocate idx as a length-nresponse vector for LU decomposition pivots
      if (.not. associated(response_matrix(iky)%eigen(ie)%idx)) &
         allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))
#endif

   end subroutine setup_response_matrix_zloc_idx

   !****************************************************************************
   !                       Get dg/dphi response matrix column
   !****************************************************************************
   ! Here, we apply a unit impulse to phi at a given point on the extended zed
   ! domain. This location is indicated via <idx>
   ! 
   ! Variables:
   ! ----------
   ! <iky> = which ky value (these are all independent, but determines the number
   !                         of eigen chains, and the spacings within each eigen chain)
   ! <ie>  = identifies which eigen chain we are looking at. Only modes which are in 
   !         the same eigen chain can communicate with eachother
   ! <idx> = location of unit impulse within the segment 
   ! <nz_ext>   = length of extended zed domain for the given eigen chain 
   ! <phi_ext>  = phi on the extended zed domain
   ! <apar_ext> = apar on the extended zed domain
   ! <bpar_ext> = bpar on the extended zed domain
   ! <pdf_ext>  = distribution function on the extended zed domain
   ! 
   ! What is being done: 
   ! -------------------
   !  - Set phi^{n+1} = 0.0
   !  - Provide a (real) unit impulse to phi^{n+1} (or Delta phi^{n+1}) at the
   !    location in the extended zed domain corresponding to index <idx>.
   !  - If using periodic boundary conditions/or treating the zonal mode, then
   !    we need to make sure phi is also periodic in zed - so apply phase shift 
   !    to opposite mode. This is only important when the unit impulse is applied 
   !    at the boundary, otherwise both end points are automatically zero. 
   !  - With this unit impulse for phi, solve the parallel streaming equation for
   !    the response of the distribution function. 
   !  - Operate on the response of the PDF with the field operator. This is 
   !    done by intgrating over velocities, and multiplying by the correct factors.
   !****************************************************************************
   subroutine get_dpdf_dphi_matrix_column(iky, ie, idx, nz_ext, nresponse, phi_ext, apar_ext, bpar_ext, pdf_ext)
   
#ifdef ISO_C_BINDING
      use mp, only: sgproc0
#endif

      use parallelisation_layouts, only: vmu_lo
      use parameters_numerical, only: time_upwind_plus
      use parameters_physics, only: include_apar, include_bpar
      use gk_implicit_terms, only: get_gke_rhs, sweep_g_zext
      use arrays, only: response_matrix
      use grids_extended_zgrid, only: periodic, phase_shift

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie, idx, nz_ext, nresponse
      complex, dimension(:), intent(out) :: phi_ext, apar_ext, bpar_ext
      complex, dimension(:, vmu_lo%llim_proc:), intent(out) :: pdf_ext

      ! Local variables
      complex, dimension(:), allocatable :: dum
      integer :: ivmu, it
      integer :: offset_apar, offset_bpar

      !-------------------------------------------------------------------------
      !                    1.A) Initialise a unit impulse in phi
      !-------------------------------------------------------------------------
      ! Provide a unit impulse to phi^{n+1} (or Delta phi^{n+1}) at the location
      ! in the extended zed domain corresponding to index 'idx'
      ! note that it is sufficient to give a unit real impulse (as opposed to
      ! separately giving real and imaginary impulse) for the following reason:
      ! split homogeneous GKE, L[f] = R[phi], into L[f1] = R[phir] and L[f2] = i*R[phii],
      ! with f = f1 + f2; then phi = df1/dphir * phir + df2/dphii * phii.
      ! however, we see that if phir = phii = 1, L[f1] = R[1] = L[-i*f2],
      ! and thus f2 = i * f1.  This gives phi = df1/dphir * (phir + i * phii) = df1/dphir * phi
      phi_ext = 0.0
      apar_ext = 0.0
      bpar_ext = 0.0
      
      ! How phi^{n+1} enters the GKE depends on whether we are solving for the
      ! non-Boltzmann pdf, h, or the guiding centre pdf, 'g'
      phi_ext(idx) = time_upwind_plus
      
      ! Need to make sure that if the mode is periodic, then the boundaries match up to
      ! a phase factor. In practice this only matters if the unit impulse is at the 
      ! boundary (i.e. <idx ==1 ) otherwise phi = 0.0 at both boundary points anyway. 
      ! TOGO-GA: check division rather than multiplication -- kept division for now to be consistent with 
      ! parallel_streaming phase shift 
      if (periodic(iky) .and. idx == 1) phi_ext(nz_ext) = phi_ext(1) / phase_shift(iky)

      ! <dum> is a scratch array that takes the place of the pdf and phi
      ! at the previous time level. It also replaces the other fields (apar and bpar).
      ! This is set to zero for the response matrix approach because there is no sources
      ! coming from the previous time step in the response matrix equation.
      allocate (dum(nz_ext)); dum = 0.0

      !-------------------------------------------------------------------------
      !      1.B) Get distribution function response to unit impulse in phi
      !-------------------------------------------------------------------------

      ! Set the flux tube index to one - eed to check, but think this is okay as 
      ! the homogeneous equation solved here for the response matrix construction 
      ! is the same for all flux tubes in the flux tube train
      it = 1

      ! Solve for the response of the distribution function given this unit impulse in <phi>
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! Calculate the RHS of the GK equation (using dum=0 as the pdf at the 
         ! previous time level, and phi_ext as the input potential) and store 
         ! it in pdf_ext
         call get_gke_rhs(ivmu, iky, ie, dum, phi_ext, dum, dum, dum, dum, pdf_ext(:, ivmu))
         ! Given the RHS of the GK equation (pdf_ext), solve for the pdf at the
         ! new time level by sweeping in zed on the extended domain; the RHS is
         ! input as 'pdf_ext' and over-written with the updated solution for the pdf
         call sweep_g_zext(iky, ie, it, ivmu, pdf_ext(:, ivmu))
      end do

      deallocate (dum)

      !-------------------------------------------------------------------------
      !          1.C) Integrate over velocity to get the matrix to invert
      !-------------------------------------------------------------------------
      ! We now have the pdf on the extended zed domain at this ky and set of 
      ! connected kx values corresponding to a unit impulse in phi at this 
      ! location now integrate over velocities to get a square response matrix.
      !
      ! integrate_over_velocity is the operator that acts on the pdf in 
      ! quasineutrality. e.g. if using g for the pdf: 
      !
      !        sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ]
      ! 
      ! (this ends the parallelisation over velocity space, so every core should have a
      ! copy of phi_ext)
      call integrate_over_velocity(pdf_ext, phi_ext, apar_ext, bpar_ext, iky, ie)

      !-------------------------------------------------------------------------
      !                          Fill in the response matrix
      !-------------------------------------------------------------------------
      ! The matrix to invert is just simply the output from the routine above
      ! so just fill in the response matrix with these values. 
      ! Note that apar_ext = 0.0, and bpar_ext = 0.0 below 
#ifdef ISO_C_BINDING
      if (sgproc0) then
#endif
         ! For phi 
         response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx) = phi_ext(:nresponse)

         ! For apar - the response column is shifted by <nresponse>, which is the number
         ! of independent zed points taken up by the <phi> response.
         ! Note that apar_ext = 0.0 here, but we still need to fill in the response
         ! matrix. 
         offset_apar = 0
         if (include_apar) then
            offset_apar = nresponse
            response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx) = apar_ext(:nresponse)
         end if

         ! For bpar - the response column is shifted by another <nresponse>, which is the 
         ! number of independent zed points taken up by the <phi> or <apar> response (depending
         ! on if include_apar = .true.)
         ! Note that bpar_ext = 0.0 here, but we still need to fill in the response
         ! matrix. 
         if (include_bpar) then
            offset_bpar = offset_apar + nresponse
            response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx) = bpar_ext(:nresponse)
         end if
#ifdef ISO_C_BINDING
      end if
#endif
      !-------------------------------------------------------------------------

   end subroutine get_dpdf_dphi_matrix_column

   !****************************************************************************
   !                       Get dg/dapar response matrix column
   !****************************************************************************
   ! Here, we apply a unit impulse to apar at a given point on the extended zed
   ! domain. This location is indicated via <idx>
   ! 
   ! Variables:
   ! ----------
   ! <iky> = which ky value (these are all independent, but determines the number
   !                         of eigen chains, and the spacings within each eigen chain)
   ! <ie>  = identifies which eigen chain we are looking at. Only modes which are in 
   !         the same eigen chain can communicate with eachother
   ! <idx> = location of unit impulse within the segment 
   ! <nz_ext>   = length of extended zed domain for the given eigen chain 
   ! <phi_ext>  = phi on the extended zed domain
   ! <apar_ext> = apar on the extended zed domain
   ! <bpar_ext> = bpar on the extended zed domain
   ! <pdf_ext>  = distribution function on the extended zed domain
   ! 
   ! What is being done: 
   ! -------------------
   !  - Set apar^{n+1} = 0.0
   !  - Provide a (real) unit impulse to apar^{n+1} at the location in the
   !    extended zed domain corresponding to index <idx>.
   !  - If using periodic boundary conditions/or treating the zonal mode, then
   !    we need to make sure apar is also periodic in zed - so apply phase shift 
   !    to opposite mode. This is only important when the unit impulse is applied 
   !    at the boundary, otherwise both end points are automatically zero. 
   !  - With this unit impulse for apar, solve the parallel streaming equation for
   !    the response of the distribution function. 
   !  - Operate on the response of the PDF with the field operator. This is 
   !    done by intgrating over velocities, and multiplying by the correct factors. 
   !****************************************************************************
   subroutine get_dpdf_dapar_matrix_column(iky, ie, idx, nz_ext, nresponse, phi_ext, apar_ext, bpar_ext, pdf_ext)

      use parallelisation_layouts, only: vmu_lo
      use parameters_numerical, only: time_upwind_plus
      use parameters_physics, only: include_apar, include_bpar
      use gk_implicit_terms, only: get_gke_rhs, sweep_g_zext
      use arrays, only: response_matrix
      use grids_extended_zgrid, only: periodic, phase_shift
#ifdef ISO_C_BINDING
      use mp, only: sgproc0
#endif

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie, idx, nz_ext, nresponse
      complex, dimension(:), intent(out) :: phi_ext, apar_ext, bpar_ext
      complex, dimension(:, vmu_lo%llim_proc:), intent(out) :: pdf_ext

      ! Local variables
      complex, dimension(:), allocatable :: dum
      integer :: ivmu, it
      integer :: offset_apar, offset_bpar
      
      !-------------------------------------------------------------------------
      !                   1.A) Initialise a unit impulse in apar
      !-------------------------------------------------------------------------
      ! Provide a unit impulse to apar^{n+1} (or Delta apar^{n+1}) at the location
      ! in the extended zed domain corresponding to index 'idx'
      ! note that it is sufficient to give a unit real impulse (as opposed to
      ! separately giving real and imaginary impulse) for the following reason:
      ! split homogeneous GKE, L[f] = R[apar], into L[f1] = R[aparr] and L[f2] = i*R[apari],
      ! with f = f1 + f2; then apar = df1/daparr * aparr + df2/dapari * apari.
      ! however, we see that if aparr = apari = 1, L[f1] = R[1] = L[-i*f2],
      ! and thus f2 = i * f1.  This gives apar = df1/daparr * (aparr + i * apari) = df1/daparr * apar
      phi_ext = 0.0
      apar_ext = 0.0
      bpar_ext = 0.0
      apar_ext(idx) = 1.0

      ! Need to make sure that if the mode is periodic, then the boundaries match up to
      ! a phase factor. In practice this only matters if the unit impulse is at the 
      ! boundary (i.e. <idx ==1 ) otherwise phi = 0.0 at both boundary points anyway. 
      if (periodic(iky) .and. idx == 1) apar_ext(nz_ext) = apar_ext(1) / phase_shift(iky)

      ! <dum> is a scratch array that takes the place of the pdf and the fields
      ! at the previous time level. It also replaces the other fields (apar and bpar).
      ! This is set to zero for the response matrix approach because there is no sources
      ! coming from the previous time step in the response matrix equation.
      allocate (dum(nz_ext)); dum = 0.0

      !-------------------------------------------------------------------------
      !      1.B) Get distribution function response to unit impulse in apar
      !-------------------------------------------------------------------------

      ! Set the flux tube index to one - eed to check, but think this is okay as 
      ! the homogeneous equation solved here for the response matrix construction 
      ! is the same for all flux tubes in the flux tube train
      it = 1

      ! Solve for the response of the distribution function given this unit impulse in <apar>
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! Calculate the RHS of the GK equation (using dum=0 as the pdf at the 
         ! previous time level, and phi_ext as the input potential) and store 
         ! it in pdf_ext
         call get_gke_rhs(ivmu, iky, ie, dum, dum, apar_ext * time_upwind_plus, apar_ext, dum, dum, pdf_ext(:, ivmu))
         ! Given the RHS of the GK equation (pdf_ext), solve for the pdf at the
         ! new time level by sweeping in zed on the extended domain; the RHS is
         ! input as 'pdf_ext' and over-written with the updated solution for the pdf
         call sweep_g_zext(iky, ie, it, ivmu, pdf_ext(:, ivmu))
      end do

      deallocate (dum)

      !-------------------------------------------------------------------------
      !          1.C) Integrate over velocity to get the matrix to invert
      !-------------------------------------------------------------------------
      ! We now have the pdf on the extended zed domain at this ky and set of 
      ! connected kx values corresponding to a unit impulse in apar at this 
      ! location now integrate over velocities to get a square response matrix.
      !
      ! integrate_over_velocity is the operator that acts on the pdf in 
      ! Parallel Ampere Law. e.g. if using g for the pdf: 
      !
      !          β sum_s Z_s n_s vth 2*B0/sqrt{π} \int d^2v vpar J_0 g
      ! 
      ! (this ends the parallelisation over velocity space, so every core should have a
      ! copy of phi_ext)
      call integrate_over_velocity(pdf_ext, phi_ext, apar_ext, bpar_ext, iky, ie)

      !-------------------------------------------------------------------------
      !                          Fill in the response matrix
      !-------------------------------------------------------------------------
      ! The matrix to invert is just simply the output from the routine above
      ! so just fill in the response matrix with these values. 
      ! Note that phi_ext = 0.0 and bpar_ext = 0.0 below
#ifdef ISO_C_BINDING
      if (sgproc0) then
#endif
         if (include_apar) then
            offset_apar = nresponse
         else
            offset_apar = 0
         end if
         ! For phi
         response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx + nresponse) = phi_ext(:nresponse)
         ! For apar
         if (include_apar) then
            response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx + offset_apar) = apar_ext(:nresponse)
         end if
         ! For bpar
         if (include_bpar) then
            offset_bpar = offset_apar + nresponse
            response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx + offset_apar) = bpar_ext(:nresponse)
         end if
#ifdef ISO_C_BINDING
      end if
#endif

   end subroutine get_dpdf_dapar_matrix_column

   !****************************************************************************
   !                       Get dg/dbpar response matrix column
   !****************************************************************************
   ! Here, we apply a unit impulse to bpar at a given point on the extended zed
   ! domain. This location is indicated via <idx>
   ! 
   ! Variables:
   ! ----------
   ! <iky> = which ky value (these are all independent, but determines the number
   !                         of eigen chains, and the spacings within each eigen chain)
   ! <ie>  = identifies which eigen chain we are looking at. Only modes which are in 
   !         the same eigen chain can communicate with eachother
   ! <idx> = location of unit impulse within the segment 
   ! <nz_ext>   = length of extended zed domain for the given eigen chain 
   ! <phi_ext>  = phi on the extended zed domain
   ! <apar_ext> = apar on the extended zed domain
   ! <bpar_ext> = bpar on the extended zed domain
   ! <pdf_ext>  = distribution function on the extended zed domain
   ! 
   ! What is being done: 
   ! -------------------
   !  - Set bpar^{n+1} = 0.0
   !  - Provide a (real) unit impulse to bpar^{n+1} at the location in the
   !    extended zed domain corresponding to index <idx>.
   !  - If using periodic boundary conditions/or treating the zonal mode, then
   !    we need to make sure bpar is also periodic in zed - so apply phase shift 
   !    to opposite mode. This is only important when the unit impulse is applied 
   !    at the boundary, otherwise both end points are automatically zero. 
   !  - With this unit impulse for bpar, solve the parallel streaming equation for
   !    the response of the distribution function. 
   !  - Operate on the response of the PDF with the field operator. This is 
   !    done by intgrating over velocities, and multiplying by the correct factors. 
   !****************************************************************************
   subroutine get_dpdf_dbpar_matrix_column(iky, ie, idx, nz_ext, nresponse, phi_ext, apar_ext, bpar_ext, pdf_ext)
   
#ifdef ISO_C_BINDING
      use mp, only: sgproc0
#endif

      use parallelisation_layouts, only: vmu_lo
      use parameters_numerical, only: time_upwind_plus
      use parameters_physics, only: include_apar, include_bpar
      use gk_implicit_terms, only: get_gke_rhs, sweep_g_zext
      use arrays, only: response_matrix
      use grids_extended_zgrid, only: periodic, phase_shift

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie, idx, nz_ext, nresponse
      complex, dimension(:), intent(out) :: phi_ext, apar_ext, bpar_ext
      complex, dimension(:, vmu_lo%llim_proc:), intent(out) :: pdf_ext

      ! Local variables
      complex, dimension(:), allocatable :: dum
      integer :: ivmu, it
      integer :: offset_apar, offset_bpar

      !-------------------------------------------------------------------------
      !                  1.A) Initialise a unit impulse in bpar
      !-------------------------------------------------------------------------
      ! Provide a unit impulse to bpar^{n+1} (or Delta bpar^{n+1}) at the location
      ! in the extended zed domain corresponding to index 'idx'
      ! note that it is sufficient to give a unit real impulse (as opposed to
      ! separately giving real and imaginary impulse) for the following reason:
      ! split homogeneous GKE, L[f] = R[bpar], into L[f1] = R[bparr] and L[f2] = i*R[bpari],
      ! with f = f1 + f2; then bpar = df1/dbparr * bparr + df2/dbpari * bpari.
      ! however, we see that if bparr = bpari = 1, L[f1] = R[1] = L[-i*f2],
      ! and thus f2 = i * f1.  This gives bpar = df1/dbparr * (bparr + i * bpari) = df1/dbparr * bpar
      phi_ext = 0.0
      apar_ext = 0.0
      bpar_ext = 0.0
      ! how phi^{n+1} enters the GKE depends on whether we are solving for the
      ! non-Boltzmann pdf, h, or the guiding centre pdf, 'g'
      bpar_ext(idx) = time_upwind_plus

      ! Need to make sure that if the mode is periodic, then the boundaries match up to
      ! a phase factor. In practice this only matters if the unit impulse is at the 
      ! boundary (i.e. <idx ==1 ) otherwise phi = 0.0 at both boundary points anyway. 
      if (periodic(iky) .and. idx == 1) bpar_ext(nz_ext) = bpar_ext(1) / phase_shift(iky)

      ! <dum> is a scratch array that takes the place of the pdf and the fields
      ! at the previous time level. It also replaces the other fields (apar and bpar).
      ! This is set to zero for the response matrix approach because there is no sources
      ! coming from the previous time step in the response matrix equation.
      allocate (dum(nz_ext)); dum = 0.0

      !-------------------------------------------------------------------------
      !      1.B) Get distribution function response to unit impulse in bpar
      !-------------------------------------------------------------------------

      ! Set the flux tube index to one - eed to check, but think this is okay as 
      ! the homogeneous equation solved here for the response matrix construction 
      ! is the same for all flux tubes in the flux tube train
      it = 1

      ! Solve for the response of the distribution function given this unit impulse in <bpar>
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! Calculate the RHS of the GK equation (using dum=0 as the pdf at the 
         ! previous time level, and phi_ext as the input potential) and store 
         ! it in pdf_ext
         call get_gke_rhs(ivmu, iky, ie, dum, dum, dum, dum, dum, bpar_ext, pdf_ext(:, ivmu))
         ! Given the RHS of the GK equation (pdf_ext), solve for the pdf at the
         ! new time level by sweeping in zed on the extended domain; the RHS is
         ! input as 'pdf_ext' and over-written with the updated solution for the pdf
         call sweep_g_zext(iky, ie, it, ivmu, pdf_ext(:, ivmu))
      end do

      deallocate (dum)

      !-------------------------------------------------------------------------
      !         1.C) Integrate over velocity to get the matrix to invert
      !-------------------------------------------------------------------------
      ! We now have the pdf on the extended zed domain at this ky and set of 
      ! connected kx values corresponding to a unit impulse in apar at this 
      ! location now integrate over velocities to get a square response matrix.
      !
      ! integrate_over_velocity is the operator that acts on the pdf in 
      ! Parallel Ampere Law. e.g. if using g for the pdf: 
      !
      !          - 2β sum_s n_s T_s 2*B0/sqrt{π} \int d^2v mu J_1/a_s g 
      ! 
      ! (this ends the parallelisation over velocity space, so every core should have a
      ! copy of phi_ext)
      call integrate_over_velocity(pdf_ext, phi_ext, apar_ext, bpar_ext, iky, ie)

      !-------------------------------------------------------------------------
      !                          Fill in the response matrix
      !-------------------------------------------------------------------------
      ! The matrix to invert is just simply the output from the routine above
      ! so just fill in the response matrix with these values. 
      ! Note that phi_ext = 0.0 and bpar_ext = 0.0 below

#ifdef ISO_C_BINDING
      if (sgproc0) then
#endif
         offset_apar = 0
         if (include_apar) offset_apar = nresponse
         if (include_bpar) offset_bpar = offset_apar + nresponse

         ! For phi         
         response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx + offset_bpar) = phi_ext(:nresponse)
         ! For apar
         if (include_apar) then
            response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx + offset_bpar) = apar_ext(:nresponse)
         end if
         ! For bpar
         if (include_bpar) then
            response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx + offset_bpar) = bpar_ext(:nresponse)
         end if
#ifdef ISO_C_BINDING
      end if
#endif

   end subroutine get_dpdf_dbpar_matrix_column

   !****************************************************************************
   !                          Perform velocity integral
   !****************************************************************************
   ! Perform the appropriate velocity integral for the fields: 
   ! For phi:     2B/sqrt(π) int dvpa int dmu J_0 * g
   ! For apar:    β sum_s Z_s n_s vth 2*B0/sqrt{π} \int d^2v vpar J_0 g
   ! For bpar:    - 2β sum_s n_s T_s 2*B0/sqrt{π} \int d^2v mu J_1/a_s g 
   !****************************************************************************
   subroutine integrate_over_velocity(g, phi, apar, bpar, iky, ie)

      use parallelisation_layouts, only: vmu_lo
      use parameters_physics, only: include_apar, include_bpar

      implicit none

      ! Arguments
      complex, dimension(:, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:), intent(out) :: phi, apar, bpar
      integer, intent(in) :: iky, ie

      !-------------------------------------------------------------------------
      call integrate_over_velocity_phi(g, phi, iky, ie)
      !-------------------------------------------------------------------------
      !                             Electromagnetic 
      !-------------------------------------------------------------------------
      if (include_apar) call integrate_over_velocity_apar(g, apar, iky, ie)
      if (include_bpar) call integrate_over_velocity_bpar(g, bpar, iky, ie)
      !-------------------------------------------------------------------------

   end subroutine integrate_over_velocity

   !****************************************************************************
   !                          Velocity integral for phi
   !****************************************************************************
   !                       2B/sqrt(π) int dvpa int dmu J_0 * g
   !****************************************************************************
   subroutine integrate_over_velocity_phi(g, phi, iky, ie)

      use mp, only: sum_allreduce
      use parallelisation_layouts, only: vmu_lo
      use grids_species, only: nspec, spec
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: nsegments
      use calculations_velocity_integrals, only: integrate_species
      use calculations_gyro_averages, only: gyro_average
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      use parameters_numerical, only: driftkinetic_implicit
      use calculations_velocity_integrals, only: integrate_species_ffs_rm
      use parameters_physics, only: full_flux_surface
      use arrays_gyro_averages, only: j0_B_const

      implicit none

      ! Arguments
      complex, dimension(:, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:), intent(out) :: phi
      integer, intent(in) :: iky, ie

      ! Local variables
      integer :: idx, iseg, ikx, iz, ia
      integer :: izl_offset
      real, dimension(nspec) :: wgt
      complex, dimension(:), allocatable :: g0
      integer :: ivmu, imu, iv, is

      !-------------------------------------------------------------------------
      allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      phi = 0.
      wgt = spec%z * spec%dens_psi0

      ia = 1
      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)

      !-------------------------------------------------------------------------
      !                             Integrals
      !-------------------------------------------------------------------------
      ! Start with the first segment. This needs to be treated seperately becuase
      ! of the additional point in the segment. 
      do iz = iz_low(iseg), iz_up(iseg)
         ! We only need to consider the zed grid point that has a unit impulse
         ! as all other grid points will be zero. 
         idx = idx + 1
         if (.not. full_flux_surface .and. (.not. driftkinetic_implicit)) then
            !-------------------------------------------------------------------
            !                             Flux tube
            !-------------------------------------------------------------------
            ! First get J0 * g
            call gyro_average(g(idx, :), iky, ikx, iz, g0)
            ! Now integrate over vpa, mu and sum over species.
            ! This returns: 2B/sqrt(π) int dvpa int dmu J_0 * g
            call integrate_species(g0, iz, wgt, phi(idx), reduce_in=.false.)
         else
            !-------------------------------------------------------------------
            !                           Full Flux Surface
            !-------------------------------------------------------------------
            ! First multiply by B and gyroaverage, but use the piece that is 
            ! constant in alpha. This is because for FFS the response matrix 
            ! only treats the part that is constant in alpha, and the remainder
            ! is treated using an iterative approach. 
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               is = is_idx(vmu_lo, ivmu)
               g0(ivmu) = g(idx, ivmu) * j0_B_const(iky, ikx, iz, ivmu)
            end do
            ! Integrate over species, using the correct weights in order to 
            ! return the constant-in-alpha component of 
            !              2B/sqrt(π) int dvpa int dmu J_0 * g
            call integrate_species_ffs_rm(g0, wgt, phi(idx), reduce_in=.false.)
            !-------------------------------------------------------------------
         end if
      end do

      ! Now treat all the other segments. This can be looped over, as all the remaining 
      ! segments have the same number of zed points.
      ! Note that this set of routines is identical to above, but is just for the 
      ! remaining segments with iseg =/= 1
      izl_offset = 1
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               if (.not. full_flux_surface .and. (.not. driftkinetic_implicit)) then
                  !-------------------------------------------------------------
                  !                       Flux tube
                  !-------------------------------------------------------------
                  ! First get J0 * g
                  call gyro_average(g(idx, :), iky, ikx, iz, g0)
                  ! Now integrate over vpa, mu and sum over species.
                  ! This returns: 2B/sqrt(π) int dvpa int dmu J_0 * g
                  call integrate_species(g0, iz, wgt, phi(idx), reduce_in=.false.)
               else
                  !-------------------------------------------------------------
                  !                     Full Flux Surface
                  !-------------------------------------------------------------
                  ! First multiply by B and gyroaverage, but use the piece that is 
                  ! constant in alpha. This is because for FFS the response matrix 
                  ! only treats the part that is constant in alpha, and the remainder
                  ! is treated using an iterative approach. 
                  do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                     iv = iv_idx(vmu_lo, ivmu)
                     imu = imu_idx(vmu_lo, ivmu)
                     is = is_idx(vmu_lo, ivmu)
                     g0(ivmu) = g(idx, ivmu) * j0_B_const(iky, ikx, iz, ivmu)
                  end do
                  ! Integrate over species, using the correct weights in order to 
                  ! return the constant-in-alpha component of 
                  !              2B/sqrt(π) int dvpa int dmu J_0 * g
                  call integrate_species_ffs_rm(g0, wgt, phi(idx), reduce_in=.false.)
               end if
               !-------------------------------------------------------------
            end do
            if (izl_offset == 0) izl_offset = 1
         end do
      end if

      !-------------------------------------------------------------------------
      ! Collect the parts of <phi> spread out across processors. 
      call sum_allreduce(phi)
      !-------------------------------------------------------------------------

   end subroutine integrate_over_velocity_phi

   !****************************************************************************
   !                            Velocity integral for apar
   !****************************************************************************
   !             β sum_s Z_s n_s vth 2*B0/sqrt{π} \int d^2v vpar J_0 g          
   !****************************************************************************
   subroutine integrate_over_velocity_apar(g, apar, iky, ie)

      use parallelisation_layouts, only: vmu_lo, iv_idx
      use parameters_physics, only: beta
      use grids_species, only: nspec, spec
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: nsegments
      use calculations_velocity_integrals, only: integrate_species
      use grids_velocity, only: vpa
      use calculations_gyro_averages, only: gyro_average
      use mp, only: sum_allreduce

      implicit none

      ! Arguments
      complex, dimension(:, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:), intent(out) :: apar
      integer, intent(in) :: iky, ie

      ! Local variables
      integer :: idx, iseg, ikx, iz, ia
      integer :: ivmu, iv
      integer :: izl_offset
      real, dimension(nspec) :: wgt
      complex, dimension(:), allocatable :: g0

      !-------------------------------------------------------------------------

      allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      apar = 0.
      wgt = spec%z * spec%dens_psi0 * spec%stm_psi0 * beta

      ia = 1
      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)

      !-------------------------------------------------------------------------
      !                             Integrals
      !-------------------------------------------------------------------------
      ! Start with the first segment. This needs to be treated seperately becuase
      ! of the additional point in the segment. 
      do iz = iz_low(iseg), iz_up(iseg)
         ! We only need to consider the zed grid point that has a unit impulse
         ! as all other grid points will be zero. 
         idx = idx + 1
         ! First get J0 * g
         call gyro_average(g(idx, :), iky, ikx, iz, g0)
         ! Multiply by vpa to give: J0 * g * vpa
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            g0(ivmu) = g0(ivmu) * vpa(iv)
         end do
         ! Now integrate over vpa, mu and sum over species.
         ! This returns: β sum_s Z_s n_s vth 2*B0/sqrt{π} \int d^2v vpar J_0 g
         call integrate_species(g0, iz, wgt, apar(idx), reduce_in=.false.)
      end do

      ! Now treat all the other segments. This can be looped over, as all the remaining 
      ! segments have the same number of zed points.
      ! Note that this set of routines is identical to above, but is just for the 
      ! remaining segments with iseg =/= 1
      izl_offset = 1
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               ! First get J0 * g
               call gyro_average(g(idx, :), iky, ikx, iz, g0)
               ! Multiply by vpa to give: J0 * g * vpa
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                  iv = iv_idx(vmu_lo, ivmu)
                  g0(ivmu) = g0(ivmu) * vpa(iv)
               end do
               ! Now integrate over vpa, mu and sum over species.
               ! This returns: β sum_s Z_s n_s vth 2*B0/sqrt{π} \int d^2v vpar J_0 g
               call integrate_species(g0, iz, wgt, apar(idx), reduce_in=.false.)
            end do
            if (izl_offset == 0) izl_offset = 1
         end do
      end if

      !-------------------------------------------------------------------------
      ! Collect the parts of <apar> spread out across processors. 
      call sum_allreduce(apar)
      !-------------------------------------------------------------------------

   end subroutine integrate_over_velocity_apar

   !****************************************************************************
   !                            Velocity integral for bpar                      
   !****************************************************************************
   !           - 2β sum_s n_s T_s 2*B0/sqrt{π} \int d^2v mu J_1/a_s g          
   !****************************************************************************
   subroutine integrate_over_velocity_bpar(g, bpar, iky, ie)

      use parallelisation_layouts, only: vmu_lo, imu_idx
      use grids_species, only: nspec, spec
      use parameters_physics, only: beta
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: nsegments
      use calculations_velocity_integrals, only: integrate_species
      use grids_velocity, only: mu
      use calculations_gyro_averages, only: gyro_average_j1
      use mp, only: sum_allreduce

      implicit none

      ! Arguments
      complex, dimension(:, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:), intent(out) :: bpar
      integer, intent(in) :: iky, ie

      ! Local variables
      integer :: idx, iseg, ikx, iz, ia, imu, ivmu
      integer :: izl_offset
      real, dimension(nspec) :: wgt
      complex, dimension(:), allocatable :: g0

      !-------------------------------------------------------------------------
      allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      bpar = 0.
      wgt = -2.0 * beta * spec%temp_psi0 * spec%dens_psi0

      ia = 1
      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      !-------------------------------------------------------------------------
      !                             Integrals
      !-------------------------------------------------------------------------
      ! Start with the first segment. This needs to be treated seperately becuase
      ! of the additional point in the segment. 
      do iz = iz_low(iseg), iz_up(iseg)
         ! We only need to consider the zed grid point that has a unit impulse
         ! as all other grid points will be zero. 
         idx = idx + 1
         ! First get J1/a_s * g , where a_s is the argument of the Bessel function.
         ! Note that 'J1' in the code is actually J1/a_s
         call gyro_average_j1(g(idx, :), iky, ikx, iz, g0)
         ! Multiply by mu to get: J1/a_s * g * mu
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            imu = imu_idx(vmu_lo, ivmu)
            g0(ivmu) = g0(ivmu) * mu(imu)
         end do
         ! Now integrate over vpa, mu and sum over species.
         ! This returns: - 2β sum_s n_s T_s 2*B0/sqrt{π} \int d^2v mu J_1/a_s g 
         call integrate_species(g0, iz, wgt, bpar(idx), reduce_in=.false.)
      end do

      ! Now treat all the other segments. This can be looped over, as all the remaining 
      ! segments have the same number of zed points.
      ! Note that this set of routines is identical to above, but is just for the 
      ! remaining segments with iseg =/= 1
      izl_offset = 1
      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            ikx = ikxmod(iseg, ie, iky)
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               ! First get J1/a_s * g , where a_s is the argument of the Bessel function.
               ! Note that 'J1' in the code is actually J1/a_s
               call gyro_average_j1(g(idx, :), iky, ikx, iz, g0)
               ! Multiply by mu to get: J1/a_s * g * mu
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                  imu = imu_idx(vmu_lo, ivmu)
                  g0(ivmu) = g0(ivmu) * mu(imu)
               end do               
               ! Now integrate over vpa, mu and sum over species.
               ! This returns: - 2β sum_s n_s T_s 2*B0/sqrt{π} \int d^2v mu J_1/a_s g 
               call integrate_species(g0, iz, wgt, bpar(idx), reduce_in=.false.)
            end do
            if (izl_offset == 0) izl_offset = 1
         end do
      end if

      !-------------------------------------------------------------------------
      ! Collect the parts of <bpar> spread out across processors. 
      call sum_allreduce(bpar)
      !-------------------------------------------------------------------------

   end subroutine integrate_over_velocity_bpar

   !============================================================================
   !****************************************************************************
   !                         Step 2) Fields response
   !****************************************************************************
   !============================================================================
   ! The next routines are dedicated to getting the fields response
   ! Solve the field equations - for local stella, this is a diagonal process
   !============================================================================

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine apply_field_solve_to_finish_response_matrix(iky)

#ifdef ISO_C_BINDING
      use mp, only: sgproc0
#endif
      use parameters_physics, only: include_apar, include_bpar
      use grids_extended_zgrid, only: neigen
      use grids_extended_zgrid, only: nsegments, nzed_segment
      use grids_extended_zgrid, only: periodic
      use arrays, only: response_matrix

      implicit none

      ! Arguments
      integer, intent(in) :: iky

      ! Local variables
      integer :: ie, idx, offset_apar, offset_bpar
      integer :: nz_ext, nresponse
      character(5) :: dist
      complex, dimension(:), allocatable :: phi_ext, apar_ext, bpar_ext

      !-------------------------------------------------------------------------
      
      dist = 'g'

      ! We are still considering a single ky value. Now we need to loop over the 
      ! eigen chains in the system - i.e. the sets of connected kx values when 
      ! linked on the extended zed domain.
      do ie = 1, neigen(iky)
#ifdef ISO_C_BINDING
         if (sgproc0) then
#endif
            ! Dimension of the extended zed domain:
            ! <nz_ext> = (number of zeds) x (number of segments on extended zed domain)
            nz_ext = nsegments(ie, iky) * nzed_segment + 1

            ! Treat zonal mode, or the case of periodic boundary conditions, 
            ! specially as this is periodic - to avoid double counting the boundary
            ! points
            if (periodic(iky)) then
               nresponse = nz_ext - 1
            else
               nresponse = nz_ext
            end if

            allocate (phi_ext(nz_ext)); phi_ext = 0.0
            allocate (apar_ext(nz_ext)); apar_ext = 0.0
            allocate (bpar_ext(nz_ext)); bpar_ext = 0.0
            
            ! Set up offset_apar and offset_bpar consistently so that the array
            ! slices below are consistent with the size of the response matrix
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
            
            ! Obtain the response matrix entries due to unit impulses in phi;
            ! this accounts for terms appearing both in quasineutrality and parallel ampere

            !-------------------------------------------------------------------
            !             Divide by fields pre-factor in field equation
            !-------------------------------------------------------------------
            ! We have the integrated response of the distribution function due to a unit
            ! impulse in the fields. This is currently stored in <response_matrix%eigen%zloc>.
            ! In order to get the correct form of the response matrix we need to first divide
            ! by the appropriate prefactor in the fields equation. 
            do idx = 1, nresponse
               phi_ext(nz_ext) = 0.0
               phi_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx)
               if (include_apar) then
                  apar_ext(nz_ext) = 0.0
                  !apar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx)
               end if
               if (include_bpar) then
                  bpar_ext(nz_ext) = 0.0
                  !bpar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx)
               end if
               call get_fields_for_response_matrix(phi_ext, apar_ext, bpar_ext, iky, ie, dist)

               !-------------------------------------------------------------------
               !                       Form matrix to be inverted 
               !-------------------------------------------------------------------
               ! Note that in order to 
               ! Next need to create column in response matrix from phi_ext and apar_ext
               ! negative sign because matrix to be inverted in streaming equation
               ! is identity matrix - response matrix
               ! add in contribution from identity matrix
               phi_ext(idx) = phi_ext(idx) - 1.0
               response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx) = -phi_ext(:nresponse)
               if (include_apar) response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx) = 0.0! -apar_ext(:nresponse)
               if (include_bpar) response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx) = 0.0! -bpar_ext(:nresponse)
            end do

            ! Obtain the response matrix entries due to unit impulses in apar;
            ! this accounts for terms appearing both in quasineutrality and parallel ampere
            if (include_apar) then
               do idx = 1, nresponse
                  phi_ext(nz_ext) = 0.0
                  !phi_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx + offset_apar)
                  apar_ext(nz_ext) = 0.0
                  apar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx + offset_apar)
                  if (include_bpar) then
                     bpar_ext(nz_ext) = 0.0
                     !bpar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx + offset_apar)
                  end if
                  call get_fields_for_response_matrix(phi_ext, apar_ext, bpar_ext, iky, ie, dist)

                  ! Next need to create column in response matrix from phi_ext and apar_ext
                  ! negative sign because matrix to be inverted in streaming equation
                  ! is identity matrix - response matrix
                  ! add in contribution from identity matrix for diagonal entries
                  apar_ext(idx) = apar_ext(idx) - 1.0
                  response_matrix(iky)%eigen(ie)%zloc(:nresponse, offset_apar + idx) = 0.0 !-phi_ext(:nresponse)
                  response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, offset_apar + idx) = -apar_ext(:nresponse)
                  if (include_bpar) response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, offset_apar + idx) = 0.0! -bpar_ext(:nresponse) 
               end do
            end if
            
            ! Obtain the response matrix entries due to unit impulses in bpar;
            ! this accounts for terms appearing both in quasineutrality and parallel ampere
            if (include_bpar) then
               do idx = 1, nresponse
                  phi_ext(nz_ext) = 0.0
                  !phi_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx + offset_bpar)
                  if (include_apar) then
                     apar_ext(nz_ext) = 0.0
                     !apar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx + offset_bpar)
                  end if
                  bpar_ext(nz_ext) = 0.0
                  bpar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx + offset_bpar)
                  call get_fields_for_response_matrix(phi_ext, apar_ext, bpar_ext, iky, ie, dist)

                  ! Next need to create column in response matrix from phi_ext and apar_ext
                  ! negative sign because matrix to be inverted in streaming equation
                  ! is identity matrix - response matrix
                  ! add in contribution from identity matrix for diagonal entries
                  bpar_ext(idx) = bpar_ext(idx) - 1.0
                  response_matrix(iky)%eigen(ie)%zloc(:nresponse, offset_bpar + idx) = -phi_ext(:nresponse)
                  if (include_apar) response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, offset_bpar + idx) = 0.0 !-apar_ext(:nresponse)
                  response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, offset_bpar + idx) = -bpar_ext(:nresponse) 
               end do
            end if

            deallocate (phi_ext, apar_ext, bpar_ext)
#ifdef ISO_C_BINDING
         end if
#endif
      end do

   end subroutine apply_field_solve_to_finish_response_matrix

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine get_fields_for_response_matrix(phi, apar, bpar, iky, ie, dist)

      use parameters_physics, only: include_apar, include_bpar

      implicit none

      ! Arguments
      complex, dimension(:), intent(in out) :: phi, apar, bpar
      integer, intent(in) :: iky, ie
      character(*), intent(in) :: dist

      !-------------------------------------------------------------------------

      if (include_bpar) then
         call calculate_phi_and_bpar_for_response_matrix(phi, bpar, iky, ie, dist)
      else
         call calculate_phi_for_response_matrix(phi, iky, ie, dist)
      end if
      if (include_apar) call get_apar_for_response_matrix(apar, iky, ie, dist)
 
   end subroutine get_fields_for_response_matrix

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine calculate_phi_for_response_matrix(phi, iky, ie, dist)

      use grids_z, only: nzgrid
      use grids_species, only: spec
      use grids_species, only: has_electron_species
      use geometry, only: dl_over_b
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: nsegments
      use grids_kxky, only: zonal_mode, akx
      use arrays, only: denominator_fields, denominator_fields_MBR
      use arrays, only: denominator_fields_h, denominator_fields_MBR_h
      use grids_species, only: adiabatic_option_switch
      use grids_species, only: adiabatic_option_fieldlineavg

      implicit none

      ! Arguments
      complex, dimension(:), intent(inout) :: phi
      integer, intent(in) :: iky, ie
      character(*), intent(in) :: dist

      ! Local variables
      integer :: idx, iseg, ikx, iz, ia
      integer :: izl_offset
      complex :: tmp
      real, dimension(:), allocatable :: gamma_fac

      !-------------------------------------------------------------------------

      ia = 1

      allocate (gamma_fac(-nzgrid:nzgrid))

      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      if (dist == 'h') then
         gamma_fac = denominator_fields_h
      else
         gamma_fac = denominator_fields(iky, ikx, :)
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
               gamma_fac = denominator_fields_h
            else
               gamma_fac = denominator_fields(iky, ikx, :)
            end if
            do iz = iz_low(iseg) + izl_offset, iz_up(iseg)
               idx = idx + 1
               phi(idx) = phi(idx) / gamma_fac(iz)
            end do
            if (izl_offset == 0) izl_offset = 1
         end do
      end if

      if (.not. has_electron_species(spec) .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          
         ! No connections for ky = 0
         if (zonal_mode(iky)) then
            iseg = 1
            tmp = sum(dl_over_b(ia, :) * phi)
            if (dist == 'h') then
               phi = phi + tmp * denominator_fields_MBR_h
            else
               phi = phi + tmp * denominator_fields_MBR(ikxmod(1, ie, iky), :)
            end if
         end if
      end if

      deallocate (gamma_fac)

   end subroutine calculate_phi_for_response_matrix

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine calculate_phi_and_bpar_for_response_matrix(phi, bpar, iky, ie, dist)

      use grids_z, only: nzgrid
      use grids_species, only: spec
      use grids_species, only: has_electron_species
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: nsegments
      use grids_kxky, only: zonal_mode, akx
      use arrays, only: denominator_fields_inv11, denominator_fields_inv13, denominator_fields_inv31, denominator_fields_inv33
      use arrays, only: denominator_fields_h
      use grids_species, only: adiabatic_option_switch
      use grids_species, only: adiabatic_option_fieldlineavg
      use mp, only: mp_abort
      
      implicit none

      ! Arguments
      complex, dimension(:), intent(inout) :: phi, bpar
      integer, intent(in) :: iky, ie
      character(*), intent(in) :: dist

      ! Local variables
      integer :: idx, iseg, ikx, iz, ia
      integer :: izl_offset
      complex :: antot1, antot3
      real, dimension(:), allocatable :: gammainv11, gammainv13, gammainv31, gammainv33

      !-------------------------------------------------------------------------

      ia = 1

      allocate (gammainv11(-nzgrid:nzgrid))
      allocate (gammainv13(-nzgrid:nzgrid))
      allocate (gammainv31(-nzgrid:nzgrid))
      allocate (gammainv33(-nzgrid:nzgrid))

      idx = 0; izl_offset = 0
      iseg = 1
      ikx = ikxmod(iseg, ie, iky)
      if (dist == 'h') then
         gammainv11 = 1.0/denominator_fields_h
         gammainv13 = 0.0
         gammainv31 = 0.0
         gammainv33 = 1.0
      else
         gammainv11 = denominator_fields_inv11(iky, ikx, :)
         gammainv13 = denominator_fields_inv13(iky, ikx, :)
         gammainv31 = denominator_fields_inv31(iky, ikx, :)
         gammainv33 = denominator_fields_inv33(iky, ikx, :)
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
               gammainv11 = 1.0/denominator_fields_h
               gammainv13 = 0.0
               gammainv31 = 0.0
               gammainv33 = 1.0
            else
               gammainv11 = denominator_fields_inv11(iky, ikx, :)
               gammainv13 = denominator_fields_inv13(iky, ikx, :)
               gammainv31 = denominator_fields_inv31(iky, ikx, :)
               gammainv33 = denominator_fields_inv33(iky, ikx, :)
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

   end subroutine calculate_phi_and_bpar_for_response_matrix

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine get_apar_for_response_matrix(apar, iky, ie, dist)

      use grids_z, only: nzgrid
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: nsegments
      use grids_kxky, only: zonal_mode, akx
      use arrays, only: apar_denom
      use arrays, only: kperp2

      implicit none

      ! Arguments
      complex, dimension(:), intent(in out) :: apar
      integer, intent(in) :: iky, ie
      character(*), intent(in) :: dist

      ! Local variables
      integer :: idx, iseg, ikx, iz, ia
      integer :: izl_offset
      real, dimension(:), allocatable :: denominator

      !-------------------------------------------------------------------------

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

   !============================================================================
   !****************************************************************************
   !                    Step 3) LU decompose the response matrix
   !****************************************************************************
   !============================================================================
   subroutine lu_decompose_response_matrix(iky)

#ifdef ISO_C_BINDING
      use mp, only: sgproc0
#endif
      use mp, only: mp_abort
      use arrays, only: response_matrix
      use parallelisation_layouts, only: lu_option_switch
      use parallelisation_layouts, only: lu_option_none, lu_option_local, lu_option_global
      use grids_extended_zgrid, only: neigen
      use linear_solve, only: lu_decomposition

      implicit none

      ! Arguments
      integer, intent(in) :: iky

      ! Local variables
      integer :: ie
      real :: dum

      !-------------------------------------------------------------------------

      ! Now we have the full response matrix. Finally, perform its LU decomposition
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
               ! Now that we have the reponse matrix for this ky and set of connected kx values
               ! get the LU decomposition so we are ready to solve the linear system
               call lu_decomposition(response_matrix(iky)%eigen(ie)%zloc, &
                                     response_matrix(iky)%eigen(ie)%idx, dum)

#ifdef ISO_C_BINDING
            end if
#endif
         end do
      end select

   end subroutine lu_decompose_response_matrix



!###############################################################################
!########################### READ RESPONSE MATRIX ##############################
!###############################################################################

   !****************************************************************************
   !                         Read the response matrix                        
   !****************************************************************************
   subroutine read_response_matrix

      use arrays, only: response_matrix
      use common_types, only: response_matrix_type
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen
      use grids_extended_zgrid, only: nsegments
      use grids_extended_zgrid, only: nzed_segment
      use grids_extended_zgrid, only: periodic
      use mp, only: proc0, job, broadcast, mp_abort
      use field_equations, only: nfields
      use file_units, only: unit_response_matrix

      implicit none

      ! Local variables
      integer :: iky, ie, nz_ext
      integer :: iky_dump, neigen_dump, naky_dump, nresponse_dump
      integer :: nresponse, nresponse_per_field
      character(len=15) :: job_str
      character(len=100) :: file_name
      integer :: ie_dump, istat
      logical, parameter :: debug = .false.

      !-------------------------------------------------------------------------

      ! All matrices handled for the job i_job are read
      ! from a single file named: responst_mat.ijob by that
      ! jobs root process
      if (proc0) then
         write (job_str, '(I1.1)') job
         file_name = './mat/response_mat.'//trim(job_str)

         open (unit=unit_response_matrix, status='old', file=file_name, &
               action='read', form='unformatted', iostat=istat)
         if (istat /= 0) then
            print *, 'Error opening response_matrix by root processor for job ', job_str
         end if

         read (unit=unit_response_matrix) naky_dump
         if (naky /= naky_dump) call mp_abort('mismatch in naky and naky_dump')
      end if

      if (.not. allocated(response_matrix)) allocate (response_matrix(naky))

      do iky = 1, naky
         if (proc0) then
            read (unit=unit_response_matrix) iky_dump, neigen_dump
            if (iky_dump /= iky .or. neigen_dump /= neigen(iky)) &
               call mp_abort('mismatch in iky_dump/neigen_dump')
         end if

         if (.not. associated(response_matrix(iky)%eigen)) &
            allocate (response_matrix(iky)%eigen(neigen(iky)))

         ! Loop over the sets of connected kx values
         do ie = 1, neigen(iky)
         
            ! Number of zeds x number of segments
            nz_ext = nsegments(ie, iky) * nzed_segment + 1

            ! Treat zonal mode specially to avoid double counting as it is periodic
            if (periodic(iky)) then
               nresponse_per_field = nz_ext - 1
            else
               nresponse_per_field = nz_ext
            end if
            nresponse = nresponse_per_field * nfields

            if (proc0) then
               read (unit=unit_response_matrix) ie_dump, nresponse_dump
               if (ie_dump /= ie .or. nresponse /= nresponse_dump) &
                  call mp_abort('mismatch in ie/nresponse_dump')
            end if

            ! For each ky and set of connected kx values,
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
               read (unit=unit_response_matrix) response_matrix(iky)%eigen(ie)%idx
               read (unit=unit_response_matrix) response_matrix(iky)%eigen(ie)%zloc
            end if

            call broadcast(response_matrix(iky)%eigen(ie)%idx)
            call broadcast(response_matrix(iky)%eigen(ie)%zloc)

         end do
      end do

      if (proc0) close (unit_response_matrix)

      if (debug) then
         print *, 'File', file_name, ' successfully read by root proc for job: ', job_str
      end if
   end subroutine read_response_matrix


!###############################################################################
!############################# FINISH RESPONSE MATRIX ##########################
!###############################################################################

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine finish_response_matrix

      use arrays, only: response_matrix
#if !defined ISO_C_BINDING

      implicit none

#else
      use arrays, only: response_window
      use mpi

      implicit none

      integer :: ierr

      !-------------------------------------------------------------------------

      if (response_window /= MPI_WIN_NULL) call mpi_win_free(response_window, ierr)
#endif

      if (allocated(response_matrix)) deallocate (response_matrix)
      initialised_response_matrix = .false.

   end subroutine finish_response_matrix

!###############################################################################
!########################## PARALLEL LU DECOMPOSITIONS #########################
!###############################################################################

#ifdef ISO_C_BINDING

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   ! This subroutine parallelises the LU decomposition on a single
   ! node using MPIs shared memory interface
   ! It also splits up jtwist the independent matrices across nodes
   ! Ideal speed up: cores_per_node*min(jtwist,ncores)
   !****************************************************************************
   subroutine parallel_LU_decomposition_local(iky)

      use mpi
      use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
      use arrays, only: response_matrix
      use mp, only: barrier, broadcast, sum_allreduce
      use mp, only: mp_comm, scope, allprocs, sharedprocs, curr_focus
      use mp, only: scrossdomprocs, sgproc0, mp_abort, real_size
      use mp, only: job, iproc, proc0, nproc, numnodes, inode
      use mp_lu_decomposition, only: lu_decomposition_local
      use job_manage, only: njobs
      use grids_extended_zgrid, only: neigen

      implicit none

      ! Arguments
      integer, intent(in) :: iky
      
      ! Local variables
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

      !-------------------------------------------------------------------------

      prior_focus = curr_focus

      call scope(sharedprocs)

      allocate (node_jobs(0:(numnodes - 1), 0:(njobs - 1))); node_jobs = .false.
      allocate (job_list(0:(nproc - 1))); job_list = 0
      allocate (eig_limits(0:numnodes, 0:(njobs - 1))); eig_limits = 0

      job_list(iproc) = job
      call sum_allreduce(job_list)

      ! Create a map of which nodes have which jobs
      if (proc0) then
         do j = 0, nproc - 1
            node_jobs(inode, job_list(j)) = .true. 
         end do
      end if

      ! Make sure all processors have this map
      call scope(allprocs)
      call mpi_allreduce &
         (MPI_IN_PLACE, node_jobs, size(node_jobs), MPI_LOGICAL, MPI_LOR, mp_comm, ierr)
      call scope(sharedprocs)

      do ijob = 0, njobs - 1
         jroot = -1
         do j = 0, nproc - 1
            if (job_list(j) == ijob) then
               jroot = j ! The first processor on this job will be the root process
               exit
            end if
         end do

         if (jroot == -1) cycle ! No processors on this node are on this job

         if (iproc == jroot) neig = neigen(iky)

         ! Broadcast number of matrices
         call broadcast(neig, jroot)

         ! Split up neig across nodes that have the current job
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

            ! Broadcast size of matrix
            call broadcast(n, jroot)

            ! Allocate the window
            call mpi_win_allocate_shared(win_size, disp_unit, MPI_INFO_NULL, mp_comm, bptr, win, ierr)

            ! Make sure all the procs have the right memory address
            if (iproc /= jroot) then
               call mpi_win_shared_query(win, jroot, win_size, disp_unit, bptr, ierr)
            end if

            ! Bind this c_ptr to our fortran matrix
            call c_f_pointer(bptr, lu, (/n, n/))

            ! Load the matrix
            if (iproc == jroot) lu = response_matrix(iky)%eigen(ie)%zloc

            ! Syncronize the processors
            call mpi_win_fence(0, win, ierr)

            ! All the processors have the matrix.
            ! Now perform LU decomposition
            call lu_decomposition_local(mp_comm, jroot, win, lu, &
                                        response_matrix(iky)%eigen(ie)%idx, dmax)

            ! Copy the decomposed matrix over
            if (iproc == jroot) response_matrix(iky)%eigen(ie)%zloc = lu

            call mpi_win_free(win, ierr)
         end do
      end do

      call scope(scrossdomprocs)

      ! Copy all the matrices across all nodes
      if (sgproc0) then
         do ie = 1, neigen(iky)
            nroot = 0
            if (needs_send .and. &
                (ie >= eig_limits(inode, job) .and. ie < eig_limits(inode + 1, job))) nroot = iproc
                
            ! First let processors know who is sending the data
            call sum_allreduce(nroot)
            
            ! Now send the data
            call broadcast(response_matrix(iky)%eigen(ie)%zloc, nroot)
            call broadcast(response_matrix(iky)%eigen(ie)%idx, nroot)
         end do
      end if

      call scope(prior_focus)

      deallocate (node_jobs, job_list, eig_limits)
      
   end subroutine parallel_LU_decomposition_local

#endif /* ISO_C_BINDING */

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   ! This subroutine parallelises the LU decomposition across
   ! all cores. Ideal speed up: ncores
   !****************************************************************************
   subroutine parallel_LU_decomposition_global(iky)

#ifdef ISO_C_BINDING
      use mp, only: sgproc0, scrossdomprocs
#endif

      use arrays, only: response_matrix
      use mp, only: barrier, broadcast, sum_allreduce
      use mp, only: mp_comm, scope, allprocs, sharedprocs, curr_focus
      use mp, only: job, iproc, proc0, nproc, mpicmplx
      use job_manage, only: njobs
      use grids_extended_zgrid, only: neigen
      use mpi
      use linear_solve, only: imaxloc

      implicit none

      ! Arguments
      integer, intent(in) :: iky

      ! Local variables
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

      !-------------------------------------------------------------------------

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

         ! Broadcast number of matrices for this job
         call broadcast(neig, job_roots(ijob))

         ! Set up communicator for cores working on a single matrix
         call mpi_comm_split(mp_comm, mod(iproc, neig), iproc, eig_comm, ierr)
         call mpi_comm_size(eig_comm, eig_cores, ierr)
         call mpi_comm_rank(eig_comm, ieig_core, ierr)

         ! Set up a communicator that crosses the previous one
         call mpi_comm_split(mp_comm, ieig_core, iproc, ceig_comm, ierr)
         call mpi_comm_rank(ceig_comm, ceig_core, ierr)

         call mpi_bcast(ceig_core, 1, MPI_INT, 0, eig_comm, ierr)

         ncomm = min(neig, nproc) ! Number of communicators

         allocate (eig_roots(0:ncomm - 1)); eig_roots = 0
         allocate (eig_limits(0:ncomm))
         allocate (row_limits(0:eig_cores))

         if (ieig_core == 0) eig_roots(ceig_core) = iproc

         call sum_allreduce(eig_roots)

         ! Split up neigen across cores
         ediv = neig / ncomm
         emod = mod(neig, ncomm)

         ! How many stages will the LU decomposition take?
         nstage = ediv
         if (emod > 0) nstage = nstage + 1

         ! Determine which parts of neigen this communicator processes
         eig_limits(0) = 1
         do j = 1, ncomm
            eig_limits(j) = eig_limits(j - 1) + ediv
            if (j <= emod) then
               eig_limits(j) = eig_limits(j) + 1
            end if
         end do

         ! Transfer the data from job root to root of subcommunicator
         do istage = 0, nstage - 1
            do j = 0, ncomm - 1
               ie = eig_limits(j) + istage
               ie_hi = eig_limits(j + 1) - 1
               if (ie > ie_hi) cycle

               if (iproc == job_roots(ijob) .and. iproc == eig_roots(j)) then ! No need for data transfer
                  n = size(response_matrix(iky)%eigen(ie)%idx)
                  allocate (lu(n, n))
                  lu = response_matrix(iky)%eigen(ie)%zloc
               else if (iproc == job_roots(ijob)) then ! Send data to subroots
                  ! Send size of matrix
                  n_send = size(response_matrix(iky)%eigen(ie)%idx)
                  call mpi_send(n_send, 1, MPI_INT, eig_roots(j), j, mp_comm, ierr)
                  ! Send matrix
                  call mpi_send(response_matrix(iky)%eigen(ie)%zloc, &
                                n_send * n_send, mpicmplx, eig_roots(j), nproc + j, mp_comm, ierr)
               else if (iproc == eig_roots(j)) then ! Subroot gets the data
                  ! Receive size of matrix
                  call mpi_recv(n, 1, MPI_INT, job_roots(ijob), j, mp_comm, status, ierr)
                  allocate (lu(n, n))
                  ! Receive matrix
                  call mpi_recv(lu, n * n, mpicmplx, job_roots(ijob), nproc + j, mp_comm, status, ierr)
               end if
            end do

            if (istage >= (eig_limits(ceig_core + 1) - eig_limits(ceig_core))) cycle ! Nothing for this communicator to do

            ! Broadcast matrix and its size across the communicator
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
            
               ! Divide up the work using row_limits
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

               ! Pivot if needed
               dmax = -1.0
               do k = j, n
                  tmp = vv(k) * abs(lu(k, j))
                  if (tmp > dmax) then
                     dmax = tmp
                     imax = k
                  end if
               end do

               if (j /= imax) then
                  dum = lu(imax, :)
                  lu(imax, :) = lu(j, :)
                  lu(j, :) = dum
                  vv(imax) = vv(j)
               end if
               if (ieig_core == 0) idx(j) = imax

               ! Get the lead multiplier
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
            
            !-------------------------------------------------------------------
            ! LU decomposition ends here
            !-------------------------------------------------------------------

            ! Copy the decomposed matrix over
            do j = 0, ncomm - 1

               ie = eig_limits(j) + istage
               ie_hi = eig_limits(j + 1) - 1
               if (ie > ie_hi) cycle

               if (iproc == job_roots(ijob) .and. iproc == eig_roots(j)) then ! No need for data transfer
                  response_matrix(iky)%eigen(ie)%zloc = lu
                  response_matrix(iky)%eigen(ie)%idx = idx
               else if (iproc == eig_roots(j)) then ! Subroot sends the data
                  ! Send indices
                  call mpi_send(idx, n, MPI_INT, job_roots(ijob), j, mp_comm, ierr)
                  ! Send matrix
                  call mpi_send(lu, n * n, mpicmplx, job_roots(ijob), nproc + j, mp_comm, ierr)
               else if (iproc == job_roots(ijob)) then ! Receive data from subroot
                  ! Receive indices
                  call mpi_recv(response_matrix(iky)%eigen(ie)%idx, &
                                n, MPI_INT, eig_roots(j), j, mp_comm, status, ierr)
                  ! Receive matrix
                  call mpi_recv(response_matrix(iky)%eigen(ie)%zloc, &
                                n * n, mpicmplx, eig_roots(j), nproc + j, mp_comm, status, ierr)
               end if
            end do
            deallocate (vv, lu, idx, dum)
         end do
         deallocate (eig_roots, eig_limits, row_limits)
      end do

#ifdef ISO_C_BINDING
      ! Copy all the matrices across all nodes
      if (sgproc0) then
         call scope(scrossdomprocs)
         do ie = 1, neigen(iky)
            call broadcast(response_matrix(iky)%eigen(ie)%zloc)
            call broadcast(response_matrix(iky)%eigen(ie)%idx)
         end do
      end if

      call scope(prior_focus)
#else
      call scope(prior_focus)

      ! Copy all the matrices across all nodes
      do ie = 1, neigen(iky)
         call broadcast(response_matrix(iky)%eigen(ie)%zloc)
         call broadcast(response_matrix(iky)%eigen(ie)%idx)
      end do
#endif
      deallocate (job_roots)
   end subroutine parallel_LU_decomposition_global

end module response_matrix
