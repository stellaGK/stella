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
   public :: read_response_matrix
   
   ! Only initialise once
   public :: initialised_response_matrix

   private

   ! Local variables
   
   ! Only initialise once
   logical :: initialised_response_matrix = .false.

contains

!###############################################################################
!########################## INITIALISE RESPONSE MATRIX #########################
!###############################################################################

   !============================================================================
   !            Main routine called to initialise the response matrix           
   !============================================================================
   subroutine init_response_matrix

      use mp, only: proc0
      use debug_flags, only: print_extra_info_to_terminal
      use parallelisation_layouts, only: verbose
      use linear_solve, only: lu_decomposition
      use parallelisation_layouts, only: iv_idx, is_idx
      use parallelisation_layouts, only: mat_gen
      use file_units, only: unit_response_matrix
      use arrays, only: response_matrix
      use grids_kxky, only: naky

      implicit none

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
      
      ! Set up response matrix utils
      call setup_response_matrix_file_io

      ! Only initialise once
      if (initialised_response_matrix) return
      initialised_response_matrix = .true.

      ! Allocate response matrix
      if (.not. allocated(response_matrix)) allocate (response_matrix(naky))
      
      ! ------------------------------------------------------------------------
      !                      Construct the response matrix                      
      ! ------------------------------------------------------------------------
      
      ! This is the main routine for computing the response matrix, where all
      ! the calculations are done.
      call construct_response_matrix

      ! ------------------------------------------------------------------------
      !                         Close the initialisation                        
      ! ------------------------------------------------------------------------
      
      ! Write the response matrix to an output file if <mat_gen> = True
      if (proc0 .and. mat_gen) then
         close (unit=unit_response_matrix)
      end if
      
   end subroutine init_response_matrix
   
   !============================================================================
   !                     Set up .io file for response matrix                    
   !============================================================================
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

      ! ------------------------------------------------------------------------

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

      use mp, only: proc0
      use parallelisation_layouts, only: mat_gen
      use arrays, only: response_matrix
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen
      use file_units, only: unit_response_matrix

      implicit none

      ! Local variables
      integer :: iky, ie

      ! ------------------------------------------------------------------------
      
      ! Loop over all ky modes, as these are all independent of one another.
      do iky = 1, naky

         ! Write the ky-value to the output file
         if (proc0 .and. mat_gen) then
            write (unit=unit_response_matrix) iky, neigen(iky)
         end if

         ! ---------------------------------------------------------------------
         !              STEP 1) Find the response matrix to invert              
         ! ---------------------------------------------------------------------
         
         ! For a given ky, we need to associate an %eigen to it - to denote the
         ! connected modes. The response matrix for each ky has neigen(ky).
         if (.not. associated(response_matrix(iky)%eigen)) allocate (response_matrix(iky)%eigen(neigen(iky)))

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

         ! ---------------------------------------------------------------------
         !                  Write the response matrix to a file                 
         ! ---------------------------------------------------------------------

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

      use mp, only: proc0, sum_allreduce, mp_abort
      use timers, only: region_start, region_end
      use parallelisation_layouts, only: mat_gen
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
      use file_units, only: unit_response_matrix
      use arrays, only: response_matrix

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie

      ! Local variables
      integer :: nz_ext, nresponse, nresponse_per_field
      integer :: offset_apar, offset_bpar
      integer :: icol
      integer, dimension(:), allocatable :: iz_from_izext, ikx_from_izext
      real, dimension(:, :), allocatable :: wgts_fields

      !-------------------------------------------------------------------------

      ! Start the timer
      call region_start('response_matrix')

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
      ! At this point <nresponse> is the number of independednt zed values for
      ! each field. This is the number of unit impulses we need to apply in order
      ! to get the response matrix.
      if (periodic(iky)) then
         nresponse_per_field = nz_ext - 1
      else
         nresponse_per_field = nz_ext
      end if

      ! If electromagnetic, we need to consider the response of phi, apar, and bpar.
      ! If we have more fields we need to apply more unit impulses -> one for each
      ! of the fields phi, apar, and bpar.
      nresponse = nresponse_per_field * nfields

      ! Write <ie> and <nresponse> to the output file
      if (proc0 .and. mat_gen) then
         write (unit=unit_response_matrix) ie, nresponse
      end if

      ! Allocate response_matrix%eigen%zloc = size of response matrix to invert
      ! Allocate response_matrix%idx = pivot index needed for LU decomposition
      call setup_response_matrix_zloc_idx(iky, ie, nresponse)

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

      ! a) Response of the pdf to the unit impulses, integrated over the
      !    velocity points of this rank (1.A, 1.B and 1.C.I)
      !$omp parallel do default(none) schedule(dynamic, 1) &
      !$omp private(icol) &
      !$omp shared(nresponse, iky, ie, nz_ext, nresponse_per_field, nfields, iz_from_izext, ikx_from_izext, wgts_fields, response_matrix)
      do icol = 1, nresponse
         call get_velocity_integrated_pdf_response(iky, ie, icol, nz_ext, nresponse_per_field, nfields, &
            iz_from_izext, ikx_from_izext, wgts_fields, response_matrix(iky)%eigen(ie)%zloc(:, icol))
      end do
      !$omp end parallel do

      ! b) Sum over the velocity points on all MPI ranks. This ends the
      !    parallelisation over velocity space, so every rank has the
      !    velocity integrals of all columns.
      call sum_allreduce(response_matrix(iky)%eigen(ie)%zloc)

      ! c) Divide by the field equations prefactors and turn the columns into
      !    the columns of the matrix to invert (1.C.II and 1.D)
      !$omp parallel do default(none) schedule(static) &
      !$omp private(icol) &
      !$omp shared(nresponse, iky, ie, nz_ext, nresponse_per_field, offset_apar, offset_bpar)
      do icol = 1, nresponse
         call fill_response_matrix_column(iky, ie, icol, nz_ext, nresponse_per_field, offset_apar, offset_bpar)
      end do
      !$omp end parallel do

      ! Deallocate temporary arrays
      deallocate (wgts_fields, iz_from_izext, ikx_from_izext)

      ! Stop the timer
      call region_end('response_matrix')

   end subroutine calculate_response_matrix_to_invert

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
   ! Column <icol> of <zloc> holds the velocity integrals of the pdf response
   ! to its unit impulse (summed over all MPI ranks). Divide by the appropriate
   ! prefactors of the field equations (1.C.II) and overwrite it with the column
   ! of the matrix to invert (1.D). Called inside an OpenMP parallel region, with
   ! each thread working on different columns of <zloc>.
   !============================================================================
   subroutine fill_response_matrix_column(iky, ie, icol, nz_ext, nresponse_per_field, &
         offset_apar, offset_bpar)

      use parameters_physics, only: include_apar, include_bpar
      use arrays, only: response_matrix

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie, icol, nz_ext, nresponse_per_field
      integer, intent(in) :: offset_apar, offset_bpar

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
      phi_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(:nresponse, icol)
      if (include_apar) apar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:offset_apar + nresponse, icol)
      if (include_bpar) bpar_ext(:nresponse) = response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:offset_bpar + nresponse, icol)

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
      response_matrix(iky)%eigen(ie)%zloc(:nresponse, icol) = -phi_ext(:nresponse)
      if (include_apar) response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, icol) = -apar_ext(:nresponse)
      if (include_bpar) response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, icol) = -bpar_ext(:nresponse)

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
   !                            Allocate zloc and idx                           
   !============================================================================
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
   ! Every MPI process stores its own copy of the response matrix.
   !============================================================================
   subroutine setup_response_matrix_zloc_idx(iky, ie, nresponse)

      use arrays, only: response_matrix

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie, nresponse

      ! For each ky and set of connected kx values, so we must have a response matrix
      ! that is N x N , with N = number of zeds per 2pi segment x number of 2pi segments

      ! Allocate zloc as an (nresponse x nresponse) matrix if not already allocated
      if (.not. associated(response_matrix(iky)%eigen(ie)%zloc)) &
         allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse, nresponse))

      ! Allocate idx as a length-nresponse vector for LU decomposition pivots
      if (.not. associated(response_matrix(iky)%eigen(ie)%idx)) &
         allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))

   end subroutine setup_response_matrix_zloc_idx
 
!===============================================================================
!                    Step 2) LU decompose the response matrix
!===============================================================================
! This section LU decomposes the response matrix for each eigen chain.
! It prepares the matrix for efficient inversion during the implicit solve.
! The decomposition uses <lu_decomposition> (PLASMA or LAPACK when available,
! see utils/linear_solve). Every MPI process factorises its own copy.
!===============================================================================

   !============================================================================
   !                       Main LU decomposition routine 
   !============================================================================
   subroutine lu_decompose_response_matrix(iky)

      use timers, only: region_start, region_end
      use arrays, only: response_matrix
      use grids_extended_zgrid, only: neigen
      use linear_solve, only: lu_decomposition

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
      ! get the LU decomposition so we are ready to solve the linear system
      do ie = 1, neigen(iky)
         call lu_decomposition(response_matrix(iky)%eigen(ie)%zloc, response_matrix(iky)%eigen(ie)%idx, dum)
      end do
      
      ! Stop the timer
      call region_end('lu_decomposition')

   end subroutine lu_decompose_response_matrix

!###############################################################################
!########################### READ RESPONSE MATRIX ##############################
!###############################################################################

   !============================================================================
   !                           Read the response matrix                         
   !============================================================================
   subroutine read_response_matrix

      use arrays, only: response_matrix
      use common_types, only: response_matrix_type
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen
      use grids_extended_zgrid, only: nzext
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

      ! ------------------------------------------------------------------------

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
            nz_ext = nzext(ie, iky)

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

   subroutine finish_response_matrix

      use arrays, only: response_matrix

      implicit none

      if (allocated(response_matrix)) deallocate (response_matrix)
      initialised_response_matrix = .false.

   end subroutine finish_response_matrix

end module response_matrix
