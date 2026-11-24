!###############################################################################
!############################### RESPONSE MATRIX ###############################
!###############################################################################
! 
! This module computes the response matrix for inverting the parallel streaming
! operator when parallel_streaming is treated implicitly.
! 
! Method:
! -------
! (Note: to see a certain part of the code, search the step it is listed under, for
! example to see where the homogeneous parallel streaming equation is solved, search 1.B)
! 
! Step 1) Get the matrix that we need to invert in implicit parallel streaming.
!         To do this we take the following steps:
! 
!         1.A) First initialise an unit impulse for each field and for each grid point in
!              the extended zed domain for a given ky and set of connected kx modes.
! 
!         1.B) Next, solve the homogeneous parallel streaming equation for the
!              distribution function, given a field with a unit impulse. This is
!              done for each field (i.e. phi, apar, bpar) that is included.
! 
!         1.C) Next, solve the field equations given the distribution function response.
! 
!              1.C.I)  Integrate over velocities to get the appropriate field operator
!                      that acts on the distibution function:
!                      For phi:     2B/sqrt(π) int dvpa int dmu J_0 * g
!                      For apar:    β sum_s Z_s n_s vth 2*B0/sqrt{π} int d^2v vpar J_0 g
!                      For bpar:  - 2β sum_s n_s T_s 2*B0/sqrt{π} int d^2v mu J_1/a_s g
! 
!              1.C.II) Divide by the appropriate pre-factor in the field equations.
! 
!         1.D) Construct the matrix that we need to invert by including the 
!              identity matrix where appropriate.
! 
! Step 2) LU decompose the matrix.
! 
! Storage:
! --------
! The LU factors are the same on every MPI process, and every process reads them
! in the back substitution of every time step. They are stored once per shared-
! memory domain (one per NUMA domain by default, see <response_matrix_sharing> in
! the namelist <parallelisation> and utils/mp_shared_memory) in one MPI shared-
! memory window for all matrices: the root of the domain builds and LU decomposes
! them (with PLASMA or LAPACK, using its OpenMP threads), and all processes of the
! domain read them. A process that does not share its domain (or if the MPI
! library has no unified shared memory) keeps its own copy, as before. The pivots
! <idx> are small and every process keeps its own copy.
! 
! Files:
! ------
! With <save_response_matrix> the LU factors are written to mat/response_mat.<job>,
! and with <load_response_matrix> they are read from it instead of being built,
! after checking that they belong to this simulation (see the section on the
! response matrix file).
! 
!###############################################################################
module response_matrix

   ! Load libraries
   use netcdf

   ! Load debug flags
   use debug_flags, only: debug => response_matrix_debug

   implicit none

   ! Make routines available to other modules
   public :: init_response_matrix
   public :: finish_response_matrix
   
   ! Only initialise once
   public :: initialised_response_matrix

   private

   ! Local variables
   
   ! Only initialise once
   logical :: initialised_response_matrix = .false.

   ! The shared-memory window that holds the LU factors of all response matrices
   ! of this shared-memory domain, if they are shared
   integer :: response_window
   logical :: shared_response_matrix = .false.
   logical :: allocated_response_matrix = .false.

   ! Whether this process writes and LU decomposes the matrices: the root of the
   ! shared-memory domain if they are shared, else every process
   logical :: builds_response_matrix = .true.

   ! The response matrix is built in blocks of columns of at most this size (bytes)
   integer, parameter :: block_bytes = 32 * 1024 * 1024

   ! The response matrix is only read from its file once, at the start
   logical :: tried_loading_response_matrix = .false.

   ! First record of the response matrix file, to recognise files of other versions
   character(len=*), parameter :: response_matrix_file_format = 'stella response matrix file, format 2 (November 2026)'

   ! The LU factors read from the file must give back the unit vector of the
   ! columns built for this simulation to this accuracy (see check_loaded_response_matrix)
   real, parameter :: load_tolerance = 1.e-6

contains

!###############################################################################
!########################## INITIALISE RESPONSE MATRIX #########################
!###############################################################################

   !============================================================================
   !            Main routine called to initialise the response matrix           
   !============================================================================
   subroutine init_response_matrix

      use debug_flags, only: print_extra_info_to_terminal
      use parallelisation_layouts, only: verbose
      use parallelisation_layouts, only: save_response_matrix, load_response_matrix

      implicit none

      logical :: loaded

      ! ------------------------------------------------------------------------
      !              Set up memory and timings for response matrix
      ! ------------------------------------------------------------------------

      ! Print info to command prompt
      if (verbose .and. print_extra_info_to_terminal) then
         write (*, '(A)') '############################################################'
         write (*, '(A)') '                       RESPONSE MATRIX                      '
         write (*, '(A)') '############################################################'
         write (*, '(A)') ''
      end if

      ! Only initialise once
      if (initialised_response_matrix) return
      initialised_response_matrix = .true.

      ! Allocate the response matrices, shared by the processes of each
      ! shared-memory domain if possible
      call allocate_response_matrix

      ! ------------------------------------------------------------------------
      !            Read the response matrix from a file, if requested
      ! ------------------------------------------------------------------------
      ! This is only done at the start of the simulation: when the time step
      ! changes, the response matrix is built again. If the file does not belong
      ! to this simulation, a warning is printed and the matrix is built instead.
      loaded = .false.
      if (load_response_matrix .and. .not. tried_loading_response_matrix) then
         tried_loading_response_matrix = .true.
         loaded = read_response_matrix_file()
      end if

      ! ------------------------------------------------------------------------
      !                      Construct the response matrix
      ! ------------------------------------------------------------------------

      ! This is the main routine for computing the response matrix, where all
      ! the calculations are done. The file is written again whenever the
      ! matrix is built, so that it holds the matrix of the last time step.
      if (.not. loaded) then
         call construct_response_matrix
         if (save_response_matrix) call write_response_matrix_file
      end if

   end subroutine init_response_matrix

!===============================================================================
!===============================================================================
!========================== CONSTRUCT RESPONSE MATRIX ==========================
!===============================================================================
!===============================================================================

   !============================================================================
   !               Main routine for constructing response martrix
   !============================================================================
   ! 
   ! Background:
   ! -----------
   ! Different ky modes are independent, so this creates sets of connected kx values:
   ! 
   ! After one full 2π orbit in zed, the eddy becomes sheared, shifting the mode
   ! to a higher kx value. To represent this effect, we apply the twist-and-shift
   ! boundary conditions, which connect these shifted kx values along an extended
   ! z-grid. For a given ky, the connected kx modes are spaced by,
   !                 δkx = 2π p ∂ι/∂r * (∂y/∂α) * (∂ψ/∂x) * ky
   ! where ∂ι/∂r is the radial derivative of the rotational transform.
   ! 
   ! This construction implies that not all kx modes are mutually connected. It
   ! is clear that for larger ky-values, the step size δkx becomes larger.
   ! Therefore, for each ky, we form distinct chains of connected modes, with the
   ! number of such chains determined by the spacing δkx.
   ! 
   ! The variable <neigen> gives the number of distinct chains of "eigenmodes"
   ! for a given ky, if the twist-and-shift boundary conditions are employed.
   ! If, instead, we use periodic boundary conditions, then all modes are 
   ! unconnected, and we have <neigen> = nakx for each ky-value.
   ! 
   ! What is done:
   ! ------------- 
   ! In this subroutine we construct the response matrix. This is split into two 
   ! main contributions: 
   !     1) Getting the velocity-space-integrated response. To do this we supply
   !        a unit impulse at every point in the extended zed domain for a given ky,
   !        and for a given set of connected kx values. With this, we then solve the
   !        parallel streaming equations to find the response of the distribution
   !        function to this unit impulse. This is done for each field that is
   !        included (phi, apar, bpar). Once this distribution function is obtained
   !        we integrate it with appropriate weights over (vpa, mu) and sum over
   !        species. This gives the quasineutrality component that acts on g.
   !        Then, we apply the field solve to find the matrix that we need to invert.
   !        This is essentially dividing by the correct factor from the field equations.
   !     2) LU decompose the matrix.
   !============================================================================
   subroutine construct_response_matrix

      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen

      implicit none

      ! Local variables
      integer :: iky, ie

      ! ------------------------------------------------------------------------
      
      ! Loop over all ky modes, as these are all independent of one another.
      do iky = 1, naky

         ! ---------------------------------------------------------------------
         !              STEP 1) Find the response matrix to invert              
         ! ---------------------------------------------------------------------
         
         ! The unit impulses (columns of the response matrix) of each (iky, ie) are
         ! independent and are computed in parallel over OpenMP threads inside
         ! calculate_response_matrix_to_invert. The (iky, ie) themselves could
         ! also be computed independently, as they are not needed together until
         ! the LU decomposition.

         ! Loop over the independent chains of modes for a given ky value. Note
         ! that <neigen> is an integer that depends on ky, as different ky modes
         ! will have a different number of eigenmode chains. For twist-and-shift,
         ! the larger ky is, the larger δkx is, and the more independent chains
         ! there are, while all kx-modes are typically connected for the smallest ky.
         ! Hence neigen(ky_min) is typically 1 while neigen(ky_max) is typically nakx.
         do ie = 1, neigen(iky)
            call calculate_response_matrix_to_invert(iky, ie)
         end do 
         
         ! (DSO) -- Here, the parallelisation over velocity space ends. At this
         ! point every processor has the response matrix for a given ky, and the
         ! LU decomposition could be parallelised locally if need be.
         ! This is preferable to parallelisation over ky as the LU
         ! decomposition (and perhaps QN) will be dominated by the
         ! ky with the most connections if twist-and-shift is used.

         ! ---------------------------------------------------------------------
         !                   STEP 2) LU decompose the matrix 
         ! ---------------------------------------------------------------------

         ! LU decompose the matrix
         call lu_decompose_response_matrix(iky)

      end do

   end subroutine construct_response_matrix

   !============================================================================
   !                Step 1) Find distribution function response
   !============================================================================
   ! This subroutine computes the matrix we need to invert. The exact routines
   ! called will depend on which fields are included in the simulation
   ! (phi, apar, bpar).
   !
   ! Each column <icol> of the matrix to invert is the response to a unit impulse
   ! in one field (phi, apar or bpar) at one point <idx> of the extended zed domain:
   !        icol = offset_field + idx,    idx = 1, ..., nresponse_per_field
   ! with offset_field = 0 for phi, offset_apar for apar and offset_bpar for bpar.
   ! All columns are independent of each other, apart from the sum over the
   ! velocity space points of the different MPI ranks in the field solve.
   !
   ! Parallelisation:
   ! ----------------
   !     a) For each column (in parallel over OpenMP threads, each thread with
   !        its own work arrays), apply the unit impulse, solve the homogeneous
   !        gyrokinetic equation for the local velocity points (1.A and 1.B) and
   !        integrate the response of the pdf over the local velocity points (1.C.I).
   !        The (rank-local) velocity integrals are stored directly in the column
   !        of <zloc>, whose rows are ordered in the same way (phi, apar, bpar),
   !        so no extra buffer is needed and the pdf itself is never stored.
   !     b) Sum <zloc> over the MPI ranks with ONE sum_allreduce per (iky, ie)
   !        (outside the OpenMP parallel region), instead of one per unit impulse.
   !     c) For each column (in parallel), divide by the prefactors of the field
   !        equations (1.C.II) and turn it into the column of the matrix to
   !        invert (1.D), in place.
   ! The velocity integrals are summed in the same order as in the original
   ! one-impulse-at-a-time implementation, so the matrix only changes by the
   ! round-off of the (larger) MPI reduction.
   !============================================================================
   subroutine calculate_response_matrix_to_invert(iky, ie)

      use timers, only: region_start, region_end
      use arrays, only: response_matrix

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie

      ! Local variables
      integer :: nresponse, icol, icol_first, ncol_block, ncols
      complex, dimension(:, :), allocatable :: block

      !-------------------------------------------------------------------------

      ! Start the timer
      call region_start('response_matrix')

      ! The columns are computed in blocks of <ncol_block> columns, so that every
      ! process only needs a private copy of one block, and the root of the
      ! shared-memory domain copies the finished columns into <zloc>
      nresponse = get_nresponse(iky, ie)
      ncol_block = max(1, min(nresponse, block_bytes / max(1, nresponse * storage_size(block) / 8)))
      allocate (block(nresponse, ncol_block))

      do icol_first = 1, nresponse, ncol_block
         ncols = min(ncol_block, nresponse - icol_first + 1)
         call compute_response_matrix_columns(iky, ie, [(icol, icol = icol_first, icol_first + ncols - 1)], block(:, :ncols))
         if (builds_response_matrix) response_matrix(iky)%eigen(ie)%zloc(:, icol_first:icol_first + ncols - 1) = block(:, :ncols)
      end do

      deallocate (block)

      ! Stop the timer
      call region_end('response_matrix')

   end subroutine calculate_response_matrix_to_invert

   !============================================================================
   !        Columns <icols> of the matrix to invert of chain <ie> of <iky>
   !============================================================================
   ! Computes the columns <icols> of the matrix to invert into <columns>, on
   ! every MPI process (see calculate_response_matrix_to_invert for the steps).
   !============================================================================
   subroutine compute_response_matrix_columns(iky, ie, icols, columns)

      use mp, only: sum_allreduce, mp_abort
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: beta
      use grids_extended_zgrid, only: nzext
      use grids_extended_zgrid, only: periodic
      use grids_extended_zgrid, only: map_to_iz_ikx_from_izext
      use grids_species, only: nspec, spec
      use grids_species, only: has_electron_species
      use grids_species, only: adiabatic_option_switch
      use grids_species, only: adiabatic_option_fieldlineavg
      use field_equations, only: nfields

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie
      integer, dimension(:), intent(in) :: icols
      complex, dimension(:, :), intent(out) :: columns

      ! Local variables
      integer :: nz_ext, nresponse_per_field
      integer :: offset_apar, offset_bpar
      integer :: k, ncols
      integer, dimension(:), allocatable :: iz_from_izext, ikx_from_izext
      real, dimension(:, :), allocatable :: wgts_fields

      !-------------------------------------------------------------------------

      ! ------------------------------------------------------------------------
      !              Set up system to compute response of the pdf
      ! ------------------------------------------------------------------------
      ! <nz_ext> is the number of points of this chain on the extended zed domain:
      ! for the twist-and-shift chains (number of zeds) x (number of 2pi segments) + 1,
      ! since every segment shares its first point with the previous one
      ! (see <build_chain_points_from_segments> in grids_extended_zgrid).
      ! ------------------------------------------------------------------------

      ! Calculate the number of z-points on the extended z-domain
      nz_ext = nzext(ie, iky)

      ! We need to treat the zonal mode specially to avoid double counting the end
      ! points in zed, as it is periodic so these end points are the same.
      ! This is also needed if 'periodic' is chosen for the boundary conditions.
      ! At this point <nresponse_per_field> is the number of independent zed values
      ! for each field. This is the number of unit impulses we need to apply in order
      ! to get the response matrix, for each of the fields phi, apar and bpar.
      if (periodic(iky)) then
         nresponse_per_field = nz_ext - 1
      else
         nresponse_per_field = nz_ext
      end if

      ! Offsets of the apar and bpar blocks in the matrix to invert (rows and columns)
      offset_apar = 0; offset_bpar = 0
      if (include_apar) offset_apar = nresponse_per_field
      if (include_bpar) offset_bpar = offset_apar + nresponse_per_field

      ! This combination is not supported, abort here (outside of the OpenMP
      ! parallel region) rather than inside the field solve of every column
      if (include_bpar .and. .not. has_electron_species(spec) .and. &
          adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         call mp_abort('adiabatic electrons not yet supported for include_bpar = T. aborting.')
      end if

      ! Map from the extended zed domain to (iz, ikx); shared by all threads
      allocate (iz_from_izext(nz_ext), ikx_from_izext(nz_ext))
      call map_to_iz_ikx_from_izext(iky, ie, iz_from_izext, ikx_from_izext)

      ! Species weights of the velocity integrals for phi, apar and bpar
      allocate (wgts_fields(nspec, 3))
      wgts_fields(:, 1) = spec%z * spec%dens_psi0
      wgts_fields(:, 2) = spec%z * spec%dens_psi0 * spec%stm_psi0 * beta
      wgts_fields(:, 3) = -2.0 * beta * spec%temp_psi0 * spec%dens_psi0

      ! ------------------------------------------------------------------------
      !    Get the matrix we need to invert depending on the fields included
      ! ------------------------------------------------------------------------
      ncols = size(icols)

      ! a) Response of the pdf to the unit impulses, integrated over the
      !    velocity points of this rank (1.A, 1.B and 1.C.I)
      !$omp parallel do default(none) schedule(dynamic, 1) &
      !$omp private(k) &
      !$omp shared(ncols, icols, iky, ie, nz_ext, nresponse_per_field, nfields, iz_from_izext, ikx_from_izext, wgts_fields, columns)
      do k = 1, ncols
         call get_velocity_integrated_pdf_response(iky, ie, icols(k), nz_ext, nresponse_per_field, nfields, &
            iz_from_izext, ikx_from_izext, wgts_fields, columns(:, k))
      end do
      !$omp end parallel do

      ! b) Sum over the velocity points on all MPI ranks. This ends the
      !    parallelisation over velocity space, so every rank has the
      !    velocity integrals of all columns.
      call sum_allreduce(columns)

      ! c) Divide by the field equations prefactors and turn the columns into
      !    the columns of the matrix to invert (1.C.II and 1.D)
      !$omp parallel do default(none) schedule(static) &
      !$omp private(k) &
      !$omp shared(ncols, icols, iky, ie, nz_ext, nresponse_per_field, offset_apar, offset_bpar, columns)
      do k = 1, ncols
         call fill_response_matrix_column(icols(k), nz_ext, nresponse_per_field, &
            offset_apar, offset_bpar, columns(:, k), iky, ie)
      end do
      !$omp end parallel do

      ! Deallocate temporary arrays
      deallocate (wgts_fields, iz_from_izext, ikx_from_izext)

   end subroutine compute_response_matrix_columns

   !============================================================================
   !          Velocity-integrated response of the pdf to a unit impulse
   !============================================================================
   ! Here, we apply a unit impulse to phi, apar or bpar at a given point on the
   ! extended zed domain, solve for the response of the distribution function
   ! for each velocity point of this MPI rank, and integrate it over the local
   ! velocity points. This routine is called inside an OpenMP parallel region,
   ! with one column <icol> per thread, so all work arrays are local.
   !
   ! Variables:
   ! ----------
   ! <iky> = which ky value (these are all independent, but determines the number
   !                         of eigen chains, and the spacings within each eigen chain)
   ! <ie>  = identifies which eigen chain we are looking at. Only modes which are in
   !         the same eigen chain can communicate with eachother
   ! <icol> = column of the response matrix, which determines the field <ifield>
   !          (1 = phi, then apar and bpar if included) and the location <idx>
   !          of the unit impulse on the extended zed domain
   ! <nz_ext>   = length of extended zed domain for the given eigen chain
   ! <nfields>  = number of fields (phi, apar, bpar) included
   ! <moments>  = velocity integrals of the pdf response for each field
   !              (local to this MPI rank, i.e. not yet summed over ranks),
   !              laid out like a column of the response matrix
   !
   ! What is being done:
   ! -------------------
   ! 1.A) First initialise (real) unit impulse phi^{n+1} (or apar^{n+1}, bpar^{n+1})
   !      at the grid point <idx> in the extended zed domain for a given connected chain.
   !      If using periodic boundary conditions/or treating the zonal mode, then
   !      we need to make sure the field is also periodic in zed - so apply phase shift
   !      to opposite mode. This is only important when the unit impulse is applied
   !      at the boundary, otherwise both end points are automatically zero.
   ! 1.B) Then solve the homogeneous parallel streaming equation for the response of
   !      the distribution function given a field with a unit impulse.
   ! 1.C.I) Then integrate over velocities to get the appropriate field operator that acts
   !      on the distibution function:
   !           For phi:     2B/sqrt(π) int dvpa int dmu J_0 * g
   !           For apar:    β sum_s Z_s n_s vth 2*B0/sqrt{π} int d^2v vpar J_0 g
   !           For bpar:    - 2β sum_s n_s T_s 2*B0/sqrt{π} int d^2v mu J_1/a_s g
   !      The integrals are accumulated directly for each velocity point, so that
   !      the pdf response is only needed for one velocity point at a time.
   !============================================================================
   subroutine get_velocity_integrated_pdf_response(iky, ie, icol, nz_ext, nresponse_per_field, nfields, &
         iz_from_izext, ikx_from_izext, wgts_fields, moments)

      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      use parameters_numerical, only: time_upwind_plus
      use parameters_numerical, only: driftkinetic_implicit
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: full_flux_surface
      use gk_implicit_terms, only: get_gke_rhs, sweep_g_zext
      use grids_extended_zgrid, only: periodic, phase_shift
      use grids_velocity, only: vpa, mu
      use grids_velocity, only: wgts_mu, wgts_vpa, wgts_mu_bare
      use arrays_gyro_averages, only: aj0x, aj1x, j0_B_const, aj_index

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie, icol, nz_ext, nresponse_per_field, nfields
      integer, dimension(:), intent(in) :: iz_from_izext, ikx_from_izext
      real, dimension(:, :), intent(in) :: wgts_fields
      complex, dimension(nresponse_per_field, nfields), intent(out) :: moments

      ! Local variables
      integer :: ifield, idx, izext, iz, ikx, ia, it
      integer :: ivmu, iv, imu, is
      integer :: islot_apar, islot_bpar
      logical :: impulse_phi, impulse_apar, impulse_bpar
      logical :: gyro_average_phi
      complex :: g0
      complex, dimension(:), allocatable :: field_impulse, field_impulse_upwinded
      complex, dimension(:), allocatable :: dum, pdf

      !-------------------------------------------------------------------------

      ! Which field (1 = phi, 2 = apar if included else bpar, 3 = bpar)
      ! and where on the extended zed domain gets the unit impulse
      ifield = (icol - 1) / nresponse_per_field + 1
      idx = icol - (ifield - 1) * nresponse_per_field
      impulse_phi = (ifield == 1)
      impulse_apar = include_apar .and. ifield == 2
      impulse_bpar = include_bpar .and. .not. (impulse_phi .or. impulse_apar)

      ! Position of the apar and bpar integrals in <moments>
      islot_apar = 2
      islot_bpar = 2
      if (include_apar) islot_bpar = 3

      ia = 1
      it = 1
      gyro_average_phi = .not. full_flux_surface .and. (.not. driftkinetic_implicit)

      ! ------------------------------------------------------------------------
      !            1.A) Initialise the unit impulse for the field
      ! ------------------------------------------------------------------------
      ! Apply the unit impulse at <idx>. For phi and bpar this is the time-centered
      ! field, i.e. time_upwind_plus * field^{n+1}; for apar we need both apar^{n+1}
      ! (for gbar -> g) and the time-centered time_upwind_plus * apar^{n+1}.
      ! If periodic, make sure the field is periodic in zed via the phase shift.
      allocate (field_impulse(nz_ext)); field_impulse = 0.0
      allocate (dum(nz_ext)); dum = 0.0
      allocate (pdf(nz_ext))
      if (impulse_apar) then
         field_impulse(idx) = 1.0
      else
         field_impulse(idx) = time_upwind_plus
      end if
      if (periodic(iky) .and. idx == 1) field_impulse(nz_ext) = field_impulse(1) / phase_shift(iky)
      if (impulse_apar) then
         allocate (field_impulse_upwinded(nz_ext))
         field_impulse_upwinded = field_impulse * time_upwind_plus
      end if

      moments = 0.0

      ! Loop over the velocity points of this MPI rank
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc

         ! ---------------------------------------------------------------------
         !        1.B) Solve the homogeneous GKE for the pdf response
         ! ---------------------------------------------------------------------
         ! Calculate the RHS of the homogeneous GK equation (i.e. the pdf is zero
         ! and only the unit impulse in one field is included) and store it in <pdf>
         if (impulse_phi) then
            call get_gke_rhs(ivmu, iky, ie, dum, field_impulse, dum, dum, dum, dum, pdf)
         else if (impulse_apar) then
            call get_gke_rhs(ivmu, iky, ie, dum, dum, field_impulse_upwinded, field_impulse, dum, dum, pdf)
         else
            call get_gke_rhs(ivmu, iky, ie, dum, dum, dum, dum, dum, field_impulse, pdf)
         end if

         ! Given the RHS of the GK equation, solve for the pdf at the new time
         ! level by sweeping in zed on the extended domain. The rhs is input
         ! as <pdf> and over-written with the updated solution for the pdf.
         call sweep_g_zext(iky, ie, it, ivmu, pdf)

         ! ---------------------------------------------------------------------
         !        1.C.I) Integrate the pdf over the local velocity points
         ! ---------------------------------------------------------------------
         ! These are the same operations, in the same order, as those of
         ! gyro_average() and integrate_species() for each zed point.
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         ! For phi:     2B/sqrt(π) int dvpa int dmu J_0 * g
         if (gyro_average_phi) then
            do izext = 1, nresponse_per_field
               iz = iz_from_izext(izext)
               ikx = ikx_from_izext(izext)
               g0 = aj0x(iky, ikx, iz, aj_index(ivmu)) * pdf(izext)
               moments(izext, 1) = moments(izext, 1) + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g0 * wgts_fields(is, 1)
            end do
         else
            do izext = 1, nresponse_per_field
               iz = iz_from_izext(izext)
               ikx = ikx_from_izext(izext)
               g0 = pdf(izext) * j0_B_const(iky, ikx, iz, ivmu)
               moments(izext, 1) = moments(izext, 1) + 2.0 * wgts_mu_bare(imu) * wgts_vpa(iv) * g0 * wgts_fields(is, 1)
            end do
         end if

         ! For apar:    β sum_s Z_s n_s vth 2*B0/sqrt{π} \int d^2v vpar J_0 g
         if (include_apar) then
            do izext = 1, nresponse_per_field
               iz = iz_from_izext(izext)
               ikx = ikx_from_izext(izext)
               g0 = aj0x(iky, ikx, iz, aj_index(ivmu)) * pdf(izext)
               g0 = g0 * vpa(iv)
               moments(izext, islot_apar) = moments(izext, islot_apar) &
                  + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g0 * wgts_fields(is, 2)
            end do
         end if

         ! For bpar:    - 2β sum_s n_s T_s 2*B0/sqrt{π} \int d^2v mu J_1/a_s g
         if (include_bpar) then
            do izext = 1, nresponse_per_field
               iz = iz_from_izext(izext)
               ikx = ikx_from_izext(izext)
               g0 = aj1x(iky, ikx, iz, aj_index(ivmu)) * pdf(izext)
               g0 = g0 * mu(imu)
               moments(izext, islot_bpar) = moments(izext, islot_bpar) &
                  + wgts_mu(ia, iz, imu) * wgts_vpa(iv) * g0 * wgts_fields(is, 3)
            end do
         end if

      end do

      deallocate (field_impulse, dum, pdf)
      if (allocated(field_impulse_upwinded)) deallocate (field_impulse_upwinded)

   end subroutine get_velocity_integrated_pdf_response

   !============================================================================
   !              Fill one column of the response matrix to invert
   !============================================================================
   ! <column> (column <icol> of the matrix) holds the velocity integrals of the
   ! pdf response to its unit impulse (summed over all MPI ranks). Divide by the
   ! appropriate prefactors of the field equations (1.C.II) and overwrite it with
   ! the column of the matrix to invert (1.D). Called inside an OpenMP parallel
   ! region, with each thread working on different columns.
   !============================================================================
   subroutine fill_response_matrix_column(icol, nz_ext, nresponse_per_field, &
         offset_apar, offset_bpar, column, iky, ie)

      use parameters_physics, only: include_apar, include_bpar

      implicit none

      ! Arguments
      integer, intent(in) :: icol, nz_ext, nresponse_per_field
      integer, intent(in) :: offset_apar, offset_bpar
      complex, dimension(:), intent(in out) :: column
      integer, intent(in) :: iky, ie

      ! Local variables
      integer :: ifield, idx, nresponse
      complex, dimension(:), allocatable :: phi_ext, apar_ext, bpar_ext

      !-------------------------------------------------------------------------

      ! Field and location of the unit impulse of this column
      ifield = (icol - 1) / nresponse_per_field + 1
      idx = icol - (ifield - 1) * nresponse_per_field
      nresponse = nresponse_per_field

      ! Velocity integrals of the pdf response for each field
      allocate (phi_ext(nz_ext), apar_ext(nz_ext), bpar_ext(nz_ext))
      phi_ext = 0.0; apar_ext = 0.0; bpar_ext = 0.0
      phi_ext(:nresponse) = column(:nresponse)
      if (include_apar) apar_ext(:nresponse) = column(offset_apar + 1:offset_apar + nresponse)
      if (include_bpar) bpar_ext(:nresponse) = column(offset_bpar + 1:offset_bpar + nresponse)

      ! ------------------------------------------------------------------------
      !            1.C.II) Divide by the appropriate field prefactors
      ! ------------------------------------------------------------------------
      call solve_field_equations_using_pdf_response(phi_ext, apar_ext, bpar_ext, iky, ie, nz_ext)

      ! ------------------------------------------------------------------------
      !                      1.D) Compute the matrix to invert
      ! ------------------------------------------------------------------------
      ! Next need to create column in response matrix from phi_ext, apar_ext and bpar_ext
      ! The negative sign occurs because the matrix acts on the RHS of the streaming equation.
      ! For the location where the unit impulse is applied, the matrix to invert is:
      !                       (identity matrix - response matrix)
      ! For all other locations the matrix to invert is:
      !                             (- response_matrix)
      ! So add in contribution from identity matrix for the <idx> location for the field we are solving for:
      if (ifield == 1) then
         phi_ext(idx) = phi_ext(idx) - 1.0
      else if (include_apar .and. ifield == 2) then
         apar_ext(idx) = apar_ext(idx) - 1.0
      else
         bpar_ext(idx) = bpar_ext(idx) - 1.0
      end if

      ! But everywhere else, simply add a negative sign:
      column(:nresponse) = -phi_ext(:nresponse)
      if (include_apar) column(offset_apar + 1:nresponse + offset_apar) = -apar_ext(:nresponse)
      if (include_bpar) column(offset_bpar + 1:nresponse + offset_bpar) = -bpar_ext(:nresponse)

      deallocate (phi_ext, apar_ext, bpar_ext)

   end subroutine fill_response_matrix_column

   !============================================================================
   !      Solve the field equations with the distribution function response
   !============================================================================
   ! We have the response of the distribution function from a unit impulse in the
   ! fields, integrated over velocity space (1.C.I, see
   ! get_velocity_integrated_pdf_response). These velocity integrals come in as
   ! <phi>, <apar> and <bpar>:
   !         For phi:     2B/sqrt(π) int dvpa int dmu J_0 * g
   !         For apar:    β sum_s Z_s n_s vth 2*B0/sqrt{π} \int d^2v vpar J_0 g
   !         For bpar:    - 2β sum_s n_s T_s 2*B0/sqrt{π} \int d^2v mu J_1/a_s g
   ! 1.C.II) Divide by the appropriate pre-factor in the field equations
   !============================================================================
   subroutine solve_field_equations_using_pdf_response(phi, apar, bpar, iky, ie, nz_ext)

      use parameters_physics, only: include_apar, include_bpar

      implicit none

      ! Arguments
      complex, dimension(:), intent(in out) :: phi, apar, bpar
      integer, intent(in) :: iky, ie, nz_ext

      ! ------------------------------------------------------------------------
      !                  1.C.II) Divide by the appropriate prefactor
      ! ------------------------------------------------------------------------
      ! Call the appropriate subroutines depending on which fields are being simulated.
      !
      ! In these routine we divide by the appropriate pre-factors that appear in the
      ! field equations in front of the fields. The exact factor will depend on
      ! which fields are being simulated (as phi and bpar are coupled), and also
      ! on the species options (e.g. adiabatic, modified Boltzmann etc.)
      !
      ! These routines cannot just be called from the field_equations modules
      ! as we need to get the correct factor given the location on the extended
      ! zed grid, so much of these routines are about mapping the prefactors onto
      ! the extended zed domain.
      ! ------------------------------------------------------------------------
      if (include_bpar) then
         call calculate_phi_and_bpar_for_response_matrix
      else
         call calculate_phi_for_response_matrix
      end if
      if (include_apar) call get_apar_for_response_matrix
      ! ------------------------------------------------------------------------

   contains

!===============================================================================
!                    1.C.II) Divide by the appropriate prefactor
!===============================================================================

      !*************************************************************************
      !              Divide by the appropriate denominator for phi
      !*************************************************************************
      ! This is the case used if bpar is not included in the simulation, as the 
      ! field equations for phi and bpar are coupled. 
      !*************************************************************************
      subroutine calculate_phi_for_response_matrix

         use grids_z, only: nzgrid
         use grids_species, only: spec
         use grids_species, only: has_electron_species
         use geometry, only: dl_over_b
         use grids_extended_zgrid, only: izext_iz, izext_ikx
         use grids_extended_zgrid, only: periodic, phase_shift
         use grids_kxky, only: zonal_mode, akx
         use arrays, only: denominator_fields, denominator_fields_MBR
         use grids_species, only: adiabatic_option_switch
         use grids_species, only: adiabatic_option_fieldlineavg

         implicit none

         ! Local variables
         integer :: idx, ikx, iz, ia
         complex :: tmp

         ! ---------------------------------------------------------------------

         ia = 1
         idx = 0
         
         
         ! ---------------------------------------------------------------------
         !                             iky = ikx = 0 mode
         ! ---------------------------------------------------------------------
         ! Stella does not evolve the iky = ikx = 0 mode. Need to identify this
         ! mode and make sure it is set to zero.
         ! ---------------------------------------------------------------------

         ! <ikx> of the first point of the chain
         ikx = izext_ikx(1, ie, iky)
         if (zonal_mode(iky) .and. abs(akx(ikx)) < epsilon(0.)) then
            phi(:) = 0.0
            return
         end if

         ! ---------------------------------------------------------------------
         !                      Divide by correct field factor
         ! ---------------------------------------------------------------------
         ! Loop over all connected segments in a chain. 
         do idx = 1, nz_ext - merge(1, 0, periodic(iky))
            ! (iz, ikx) location of this point of the chain on the extended z-grid
            iz = izext_iz(idx, ie, iky)
            ikx = izext_ikx(idx, ie, iky)
            phi(idx) = phi(idx) / denominator_fields(iky, ikx, iz)
         end do

         ! Treat the periodic point correct by dividing by the phase shift
         if (periodic(iky)) phi(nz_ext) = phi (1) / phase_shift(iky)



         ! ---------------------------------------------------------------------
         !                           Adiabatic electrons
         ! ---------------------------------------------------------------------
         
         ! If using adiabatic electrons, or the modified Boltzmann response
         ! then we need to modify the denominator appropriately. 
         if (.not. has_electron_species(spec) .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         
            ! No connections for ky = 0
            if (zonal_mode(iky)) then
               tmp = sum(dl_over_b(ia, :) * phi)
               phi = phi + tmp * denominator_fields_MBR(izext_ikx(1, ie, iky), :)
            end if
            
         end if


      end subroutine calculate_phi_for_response_matrix

      !*************************************************************************
      !           Divide by the appropriate denominator for phi + bpar
      !*************************************************************************
      ! This is the case used if bpar is included in the simulation, as the 
      ! field equations for phi and bpar are coupled. 
      !*************************************************************************
      subroutine calculate_phi_and_bpar_for_response_matrix

         use mp, only: mp_abort
         use grids_z, only: nzgrid
         use grids_kxky, only: zonal_mode, akx
         use grids_species, only: has_electron_species
         use grids_extended_zgrid, only: izext_iz, izext_ikx
         use grids_extended_zgrid, only: periodic, phase_shift
         use grids_species, only: spec
         use grids_species, only: adiabatic_option_switch
         use grids_species, only: adiabatic_option_fieldlineavg
         use arrays, only: denominator_fields_inv11
         use arrays, only: denominator_fields_inv13
         use arrays, only: denominator_fields_inv31
         use arrays, only: denominator_fields_inv33
         
         implicit none

         ! Local variables
         integer :: idx, ikx, iz, ia
         complex :: antot1, antot3

         ! ---------------------------------------------------------------------

         ia = 1


         idx = 0

         ! ---------------------------------------------------------------------
         !                             iky = ikx = 0 mode
         ! ---------------------------------------------------------------------
         ! Stella does not evolve the iky = ikx = 0 mode. Need to identify this
         ! mode and make sure it is set to zero.
         ! ---------------------------------------------------------------------

         ! <ikx> of the first point of the chain
         ikx = izext_ikx(1, ie, iky)
         if (zonal_mode(iky) .and. abs(akx(ikx)) < epsilon(0.)) then
            phi(:) = 0.0
            bpar(:) = 0.0
            return
         end if

         ! ---------------------------------------------------------------------
         !                      Divide by correct field factor
         ! ---------------------------------------------------------------------
         ! Loop over all connected segments in a chain.
         do idx = 1, nz_ext - merge(1, 0, periodic(iky))
            ! (iz, ikx) location of this point of the chain on the extended z-grid
            iz = izext_iz(idx, ie, iky)
            ikx = izext_ikx(idx, ie, iky)
            antot1 = phi(idx)
            antot3 = bpar(idx)
            phi(idx) = antot1 * denominator_fields_inv11(iky, ikx, iz) + antot3 * denominator_fields_inv13(iky, ikx, iz)
            bpar(idx) = antot1 * denominator_fields_inv31(iky, ikx, iz) + antot3 * denominator_fields_inv33(iky, ikx, iz)
         end do

         ! Treat the periodic point correct by dividing by the phase shift
         if (periodic(iky)) phi(nz_ext) = phi(1) / phase_shift(iky)
         if (periodic(iky)) bpar(nz_ext) = bpar(1) / phase_shift(iky)
         
         

         ! ---------------------------------------------------------------------
         !                           Adiabatic electrons
         ! ---------------------------------------------------------------------
         ! If using adiabatic electrons, or the modified Boltzmann response
         ! then we need to modify the denominator appropriately. 
         if (.not. has_electron_species(spec) .and. &
            adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            call mp_abort('adiabatic electrons not yet supported for include_bpar = T. aborting.')
         end if


      end subroutine calculate_phi_and_bpar_for_response_matrix

      !*************************************************************************
      !              Divide by the appropriate denominator for apar
      !*************************************************************************
      ! This is the routine used if apar is included. Note that it is separated 
      ! from the phi and bpar equations, as the field equation for apar is decoupled.
      !*************************************************************************
      subroutine get_apar_for_response_matrix

         use grids_z, only: nzgrid
         use grids_extended_zgrid, only: izext_iz, izext_ikx
         use grids_extended_zgrid, only: periodic, phase_shift
         use grids_kxky, only: zonal_mode, akx
         use arrays, only: kperp2

         implicit none

         ! Local variables
         integer :: idx, ikx, iz, ia

         ! ---------------------------------------------------------------------

         ia = 1
         idx = 0

         
         ! ---------------------------------------------------------------------
         !                             iky = ikx = 0 mode
         ! ---------------------------------------------------------------------
         ! Stella does not evolve the iky = ikx = 0 mode. Need to identify this
         ! mode and make sure it is set to zero.
         ! ---------------------------------------------------------------------

         ! <ikx> of the first point of the chain
         ikx = izext_ikx(1, ie, iky)
         if (zonal_mode(iky) .and. abs(akx(ikx)) < epsilon(0.)) then
            apar(:) = 0.0
            return
         end if

         ! ---------------------------------------------------------------------
         !                      Divide by correct field factor
         ! ---------------------------------------------------------------------
         ! Loop over all connected segments in a chain. 
         do idx = 1, nz_ext - merge(1, 0, periodic(iky))
            ! (iz, ikx) location of this point of the chain on the extended z-grid
            iz = izext_iz(idx, ie, iky)
            ikx = izext_ikx(idx, ie, iky)
            apar(idx) = apar(idx) / kperp2(iky, ikx, ia, iz)
         end do
         
         ! Treat the periodic point correct by dividing by the phase shift
         if (periodic(iky)) apar(nz_ext) = apar(1) / phase_shift(iky)

         


      end subroutine get_apar_for_response_matrix

   end subroutine solve_field_equations_using_pdf_response

   !============================================================================
   !                 Size of the response matrix of (iky, ie)
   !============================================================================
   ! For each ky and set of connected kx values, the matrix is N x N, with
   ! N = (number of zeds per 2pi segment x number of 2pi segments) x nfields.
   ! The zonal mode (and periodic boundary conditions) do not count the end
   ! points in zed twice, as they are the same.
   !============================================================================
   integer function get_nresponse(iky, ie) result(nresponse)

      use grids_extended_zgrid, only: nzext, periodic
      use field_equations, only: nfields

      implicit none

      integer, intent(in) :: iky, ie

      if (periodic(iky)) then
         nresponse = (nzext(ie, iky) - 1) * nfields
      else
         nresponse = nzext(ie, iky) * nfields
      end if

   end function get_nresponse

   !============================================================================
   !                     Allocate all response matrices
   !============================================================================
   ! - response_matrix%eigen%zloc is the (nresponse x nresponse) matrix of the
   !   chain <ie> of <iky>, which is overwritten by its LU factors. If the
   !   processes of a shared-memory domain share it, it points into one shared
   !   window for all matrices, allocated by the root of the domain (a window
   !   for every matrix would be slow on some clusters).
   ! - response_matrix%idx keeps track of the permutations made during the LU
   !   decomposition; it is input to the LU back substitution. Every process
   !   keeps its own copy.
   !============================================================================
   subroutine allocate_response_matrix

      use, intrinsic :: iso_c_binding, only: c_ptr, c_intptr_t, c_f_pointer
      use mpi, only: MPI_ADDRESS_KIND, MPI_WIN_NULL
      use mp, only: proc0, mp_comm
      use mp_shared_memory, only: init_shared_memory_domains, allocate_shared_window
      use mp_shared_memory, only: domain_size, ndomains, min_domain_size, max_domain_size
      use mp_shared_memory, only: nprivate_processes, domain_root
      use parallelisation_layouts, only: response_matrix_sharing
      use debug_flags, only: print_extra_info_to_terminal
      use parallelisation_layouts, only: verbose
      use arrays, only: response_matrix
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen

      implicit none

      ! Alignment of every matrix in the window (bytes)
      integer(kind=MPI_ADDRESS_KIND), parameter :: alignment = 64

      integer :: iky, ie, nresponse
      integer(kind=MPI_ADDRESS_KIND) :: nbytes, offset
      integer(c_intptr_t) :: address
      type(c_ptr) :: base
      complex :: z
      character(len=12) :: size_mb

      ! -------------------------------------------------------------------------

      ! The matrices are allocated once, and reused when the time step changes
      if (allocated_response_matrix) return
      allocated_response_matrix = .true.

      if (.not. allocated(response_matrix)) allocate (response_matrix(naky))
      do iky = 1, naky
         if (.not. associated(response_matrix(iky)%eigen)) allocate (response_matrix(iky)%eigen(neigen(iky)))
      end do

      ! Split the processes into shared-memory domains
      call init_shared_memory_domains(mp_comm, response_matrix_sharing)

      ! Size of all matrices, each aligned to <alignment> bytes
      nbytes = 0
      do iky = 1, naky
         do ie = 1, neigen(iky)
            nresponse = get_nresponse(iky, ie)
            nbytes = nbytes + aligned(int(nresponse, MPI_ADDRESS_KIND)**2 * storage_size(z) / 8)
         end do
      end do

      ! One shared window for all matrices of the domain
      response_window = MPI_WIN_NULL
      shared_response_matrix = .false.
      if (domain_size > 1) then
         call allocate_shared_window(nbytes, response_window, base)
         shared_response_matrix = (response_window /= MPI_WIN_NULL)
      end if
      builds_response_matrix = domain_root .or. .not. shared_response_matrix

      ! Print the storage to the command prompt
      if (proc0 .and. verbose .and. print_extra_info_to_terminal) then
         write (size_mb, '(F12.1)') real(nbytes) / 1024.**2
         if (shared_response_matrix) then
            write (*, '(A, A, A, I0, A, I0, A, I0, A)') 'The response matrices (', trim(adjustl(size_mb)), &
               ' MiB) are shared by the processes of each of the ', ndomains, ' shared-memory domains (', &
               min_domain_size, ' to ', max_domain_size, ' processes).'
            if (nprivate_processes > 0) write (*, '(I0, A)') nprivate_processes, &
               ' processes keep their own copy, since their cores lie in more than one NUMA domain.'
         else
            write (*, '(A, A, A)') 'Every process keeps its own copy of the response matrices (', &
               trim(adjustl(size_mb)), ' MiB).'
         end if
         write (*, *) ''
      end if

      ! Point every matrix into the window, or allocate it
      offset = 0
      do iky = 1, naky
         do ie = 1, neigen(iky)
            nresponse = get_nresponse(iky, ie)
            if (shared_response_matrix) then
               address = transfer(base, address) + int(offset, c_intptr_t)
               call c_f_pointer(transfer(address, base), response_matrix(iky)%eigen(ie)%zloc, [nresponse, nresponse])
               offset = offset + aligned(int(nresponse, MPI_ADDRESS_KIND)**2 * storage_size(z) / 8)
            else
               allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse, nresponse))
            end if
            allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))
         end do
      end do

   contains

      integer(kind=MPI_ADDRESS_KIND) function aligned(n)
         integer(kind=MPI_ADDRESS_KIND), intent(in) :: n
         aligned = ((n + alignment - 1) / alignment) * alignment
      end function aligned

   end subroutine allocate_response_matrix
 
!===============================================================================
!                    Step 2) LU decompose the response matrix
!===============================================================================
! This section LU decomposes the response matrix for each eigen chain.
! It prepares the matrix for efficient inversion during the implicit solve.
! The decomposition uses <lu_decomposition> (PLASMA or LAPACK when available,
! see utils/linear_solve). The root of every shared-memory domain factorises
! the copy of its domain.
!===============================================================================

   !============================================================================
   !                       Main LU decomposition routine 
   !============================================================================
   subroutine lu_decompose_response_matrix(iky)

      use timers, only: region_start, region_end
      use arrays, only: response_matrix
      use grids_extended_zgrid, only: neigen
      use linear_solve, only: lu_decomposition
      use mp_shared_memory, only: fence_shared_window, broadcast_in_domain

      implicit none

      ! Arguments
      integer, intent(in) :: iky

      ! Local variables
      integer :: ie
      real :: dum

      ! ------------------------------------------------------------------------
      
      ! Start the timer
      call region_start('lu_decomposition')

      ! Now that we have the reponse matrix for this ky and set of connected kx values
      ! get the LU decomposition so we are ready to solve the linear system. Only
      ! the root of the shared-memory domain has written the matrix; it decomposes
      ! it (with its OpenMP threads) and gives the pivots to the other processes.
      do ie = 1, neigen(iky)
         if (builds_response_matrix) call lu_decomposition(response_matrix(iky)%eigen(ie)%zloc, response_matrix(iky)%eigen(ie)%idx, dum)
      end do
      if (shared_response_matrix) then
         call fence_shared_window(response_window)
         do ie = 1, neigen(iky)
            call broadcast_in_domain(response_matrix(iky)%eigen(ie)%idx)
         end do
      end if
      
      ! Stop the timer
      call region_end('lu_decomposition')

   end subroutine lu_decompose_response_matrix

!###############################################################################
!################### WRITE AND READ THE RESPONSE MATRIX FILE ###################
!###############################################################################
! With <save_response_matrix> = .true., the LU factors and pivots of all
! response matrices are written to the Fortran unformatted file
! mat/response_mat.<job> (by proc0, every time the matrix is built), and with
! <load_response_matrix> = .true. they are read from it at the start of the
! simulation instead of being built. The file holds:
!     <response_matrix_file_format>
!     naky, nfields, code_dt
!     for every ky:          iky, neigen(iky)
!        for every chain:    ie, nresponse
!                            idx (pivots)
!                            zloc (LU factors)
! The matrix depends on every input of the simulation (grids, geometry, species,
! physics, time step and its centring), while only the sizes and the time step
! are stored. A matrix that is read is therefore only used if
!   1) the file has this format, the same naky, nfields and time step, and
!      the same number and sizes of the chains, and
!   2) it solves the response matrix of this simulation: for every chain, the
!      first and the middle column of every field are built (a unit impulse
!      each, which costs a small part of building the whole matrix), and the
!      loaded LU factors must give back the unit vector to <load_tolerance>.
! Otherwise a warning says why, and the matrix is built instead.
!###############################################################################

   !============================================================================
   !                     Name of the response matrix file
   !============================================================================
   function response_matrix_file_name() result(file_name)

      use mp, only: job

      implicit none

      character(len=100) :: file_name

      write (file_name, '(A, I0)') './mat/response_mat.', job

   end function response_matrix_file_name

   !============================================================================
   !                   Write the response matrices to the file
   !============================================================================
   subroutine write_response_matrix_file

      use mp, only: proc0
      use system_fortran, only: systemf
      use arrays, only: response_matrix
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen
      use grids_time, only: code_dt
      use field_equations, only: nfields
      use file_units, only: unit_response_matrix

      implicit none

      integer :: iky, ie

      ! ------------------------------------------------------------------------

      ! proc0 is the root of its shared-memory domain, so it holds the factors
      if (.not. proc0) return

      call systemf('mkdir -p mat')
      open (unit=unit_response_matrix, file=trim(response_matrix_file_name()), status='replace', &
            action='write', form='unformatted')
      write (unit_response_matrix) response_matrix_file_format
      write (unit_response_matrix) naky, nfields, code_dt
      do iky = 1, naky
         write (unit_response_matrix) iky, neigen(iky)
         do ie = 1, neigen(iky)
            write (unit_response_matrix) ie, get_nresponse(iky, ie)
            write (unit_response_matrix) response_matrix(iky)%eigen(ie)%idx
            write (unit_response_matrix) response_matrix(iky)%eigen(ie)%zloc
         end do
      end do
      close (unit_response_matrix)

   end subroutine write_response_matrix_file

   !============================================================================
   !       Read the response matrices from the file, and check them
   !============================================================================
   ! Returns .true. if the matrices were read and belong to this simulation.
   ! proc0 first reads the whole file into its own copy (the shared window of its
   ! domain, or its private copy), then the factors are sent to the other copies.
   !============================================================================
   logical function read_response_matrix_file() result(loaded)

      use mp, only: proc0, broadcast, mpicmplx
      use mp_shared_memory, only: broadcast_to_domain_roots, fence_shared_window
      use arrays, only: response_matrix
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen
      use grids_time, only: code_dt
      use field_equations, only: nfields
      use file_units, only: unit_response_matrix
      use timers, only: region_start, region_end

      implicit none

      ! Local variables
      character(len=len(response_matrix_file_format)) :: format_in_file
      integer :: naky_in_file, nfields_in_file, iky_in_file, neigen_in_file
      integer :: ie_in_file, nresponse_in_file
      real :: code_dt_in_file
      integer :: iky, ie, ios
      character(len=300) :: problem
      character(len=100) :: file_name

      ! ------------------------------------------------------------------------

      call region_start('response_matrix')
      file_name = response_matrix_file_name()
      problem = ''

      ! 1) proc0 reads the file and checks its format and sizes
      if (proc0) then
         open (unit=unit_response_matrix, file=trim(file_name), status='old', action='read', &
               form='unformatted', iostat=ios)
         if (ios /= 0) then
            problem = 'it does not exist'
         else
            read (unit_response_matrix, iostat=ios) format_in_file
            if (ios /= 0 .or. format_in_file /= response_matrix_file_format) then
               problem = 'it was not written in this format (it may have been written by another version of stella)'
            end if
            if (problem == '') then
               read (unit_response_matrix, iostat=ios) naky_in_file, nfields_in_file, code_dt_in_file
               if (ios /= 0) then
                  problem = 'it could not be read'
               else if (naky_in_file /= naky) then
                  write (problem, '(A, I0, A, I0)') 'it has naky = ', naky_in_file, ' instead of ', naky
               else if (nfields_in_file /= nfields) then
                  write (problem, '(A, I0, A, I0)') 'it has ', nfields_in_file, ' fields instead of ', nfields
               else if (abs(code_dt_in_file - code_dt) > 1.e-10 * abs(code_dt)) then
                  write (problem, '(A, ES12.5, A, ES12.5)') 'it was made with the time step code_dt = ', &
                     code_dt_in_file, ' instead of ', code_dt
               end if
            end if
            do iky = 1, naky
               if (problem /= '') exit
               read (unit_response_matrix, iostat=ios) iky_in_file, neigen_in_file
               if (ios /= 0 .or. iky_in_file /= iky .or. neigen_in_file /= neigen(iky)) then
                  write (problem, '(A, I0)') 'it has another number of chains of connected kx modes for iky = ', iky
                  exit
               end if
               do ie = 1, neigen(iky)
                  read (unit_response_matrix, iostat=ios) ie_in_file, nresponse_in_file
                  if (ios /= 0 .or. ie_in_file /= ie .or. nresponse_in_file /= get_nresponse(iky, ie)) then
                     write (problem, '(A, I0, A, I0)') 'it has another size of the response matrix for iky = ', &
                        iky, ' and chain ', ie
                     exit
                  end if
                  read (unit_response_matrix, iostat=ios) response_matrix(iky)%eigen(ie)%idx
                  if (ios == 0) read (unit_response_matrix, iostat=ios) response_matrix(iky)%eigen(ie)%zloc
                  if (ios /= 0) then
                     problem = 'it could not be read'
                     exit
                  end if
               end do
            end do
            close (unit_response_matrix)
         end if
      end if

      ! 2) Every process needs the pivots, and every copy of the factors (the root
      !    of every shared-memory domain, or every process) is filled
      call broadcast(problem)
      if (problem == '') then
         do iky = 1, naky
            do ie = 1, neigen(iky)
               call broadcast(response_matrix(iky)%eigen(ie)%idx)
               if (shared_response_matrix) then
                  call broadcast_to_domain_roots(response_matrix(iky)%eigen(ie)%zloc, mpicmplx)
               else
                  call broadcast(response_matrix(iky)%eigen(ie)%zloc)
               end if
            end do
         end do
         if (shared_response_matrix) call fence_shared_window(response_window)

         ! 3) The factors must solve the response matrix of this simulation
         call check_loaded_response_matrix(problem)
      end if

      loaded = (problem == '')
      if (proc0) then
         if (loaded) then
            write (*, '(A)') 'The response matrix has been read from '//trim(file_name)//'.'
         else
            write (*, '(A)') 'WARNING: the response matrix is not read from '//trim(file_name)//', since '// &
               trim(problem)//'. It is built instead.'
         end if
      end if
      call region_end('response_matrix')

   end function read_response_matrix_file

   !============================================================================
   !      Check that the loaded LU factors solve this simulation's matrix
   !============================================================================
   ! For every chain, the first and the middle column of every field are built
   ! for this simulation and solved with the loaded LU factors, which must give
   ! back the unit vector of that column. <problem> is set if they do not.
   !============================================================================
   subroutine check_loaded_response_matrix(problem)

      use mp, only: max_allreduce
      use arrays, only: response_matrix
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen
      use field_equations, only: nfields
      use linear_solve, only: lu_back_substitution

      implicit none

      character(len=*), intent(in out) :: problem

      integer :: iky, ie, ifield, k, nresponse, nresponse_per_field
      integer, dimension(2 * 3) :: icols
      complex, dimension(:, :), allocatable :: columns
      real :: error

      ! ------------------------------------------------------------------------

      error = 0.
      do iky = 1, naky
         do ie = 1, neigen(iky)
            nresponse = get_nresponse(iky, ie)
            nresponse_per_field = nresponse / nfields
            do ifield = 1, nfields
               icols(2 * ifield - 1) = (ifield - 1) * nresponse_per_field + 1
               icols(2 * ifield) = (ifield - 1) * nresponse_per_field + (nresponse_per_field + 1) / 2
            end do
            allocate (columns(nresponse, 2 * nfields))
            call compute_response_matrix_columns(iky, ie, icols(:2 * nfields), columns)
            do k = 1, 2 * nfields
               call lu_back_substitution(response_matrix(iky)%eigen(ie)%zloc, response_matrix(iky)%eigen(ie)%idx, columns(:, k))
               columns(icols(k), k) = columns(icols(k), k) - 1.0
               error = max(error, maxval(abs(columns(:, k))))
            end do
            deallocate (columns)
         end do
      end do
      call max_allreduce(error)

      if (error > load_tolerance) then
         write (problem, '(A, ES9.2, A)') 'its LU factors do not solve the response matrix of this simulation (error ', &
            error, '), so the input file was changed after the file was written'
      end if

   end subroutine check_loaded_response_matrix


!###############################################################################
!############################# FINISH RESPONSE MATRIX ##########################
!###############################################################################

   subroutine finish_response_matrix

      use arrays, only: response_matrix
      use mp_shared_memory, only: free_shared_window, finish_shared_memory_domains

      implicit none

      integer :: iky, ie

      if (allocated(response_matrix)) then
         do iky = 1, size(response_matrix)
            if (.not. associated(response_matrix(iky)%eigen)) cycle
            do ie = 1, size(response_matrix(iky)%eigen)
               if (shared_response_matrix) then
                  nullify (response_matrix(iky)%eigen(ie)%zloc)
               else if (associated(response_matrix(iky)%eigen(ie)%zloc)) then
                  deallocate (response_matrix(iky)%eigen(ie)%zloc)
               end if
               if (associated(response_matrix(iky)%eigen(ie)%idx)) deallocate (response_matrix(iky)%eigen(ie)%idx)
            end do
            deallocate (response_matrix(iky)%eigen)
         end do
         deallocate (response_matrix)
      end if
      if (shared_response_matrix) call free_shared_window(response_window)
      if (allocated_response_matrix) call finish_shared_memory_domains
      shared_response_matrix = .false.
      allocated_response_matrix = .false.
      builds_response_matrix = .true.
      initialised_response_matrix = .false.

   end subroutine finish_response_matrix

end module response_matrix
