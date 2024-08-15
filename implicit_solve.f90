module implicit_solve

  use debug_flags, only: debug => implicit_solve_debug
  
   implicit none

   public :: time_implicit_advance
   public :: sweep_zed_zonal
   public :: advance_implicit_terms
   public :: get_gke_rhs
   public :: sweep_g_zext

   private

   real, dimension(2, 3) :: time_implicit_advance = 0.

contains

   subroutine advance_implicit_terms(g, phi, apar, bpar)

     use mp, only: proc0
     
      use job_manage, only: time_message
      use stella_layouts, only: vmu_lo
      use parameters_physics, only: include_apar, include_bpar
      use zgrid, only: nzgrid, ntubes
      use parameters_kxky_grids, only: naky, nakx
      use arrays_dist_fn, only: g1, g2
      use parameters_numerical, only: stream_matrix_inversion
      use parameters_numerical, only: use_deltaphi_for_response_matrix
      use parameters_numerical, only: tupwnd_p => time_upwind_plus
      use parameters_numerical, only: tupwnd_m => time_upwind_minus
      use parameters_numerical, only: fphi
      use fields, only: advance_fields, fields_updated
      use extended_zgrid, only: map_to_extended_zgrid, map_from_extended_zgrid
      use extended_zgrid, only: nsegments, nzed_segment

      use parameters_numerical, only: driftkinetic_implicit
      use gyro_averages, only: gyro_average
      use stella_layouts, only: iv_idx, imu_idx, is_idx

      use fields, only: get_fields_source
      use parameters_numerical, only: nitt

      use ffs_solve, only: get_source_ffs_itteration, get_drifts_ffs_itteration
      use species, only: has_electron_species
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar

      integer :: nz_ext
      complex, dimension(:, :, :, :), allocatable :: phi_old, apar_old, bpar_old
      complex, dimension(:, :, :, :), allocatable :: phi_source, apar_source, bpar_source
      character(5) :: dist_choice

      complex, dimension(:, :, :, :, :), allocatable :: phi_source_ffs
      complex, dimension(:, :, :, :), allocatable :: fields_source_ffs
      complex, dimension(:, :, :, :, :), allocatable :: drifts_source_ffs

      integer :: itt
      logical :: modify

      debug = debug .and. proc0
      if(debug) write (*,*) 'No debug messages for implicit_solve.f90 yet'
      
      if(driftkinetic_implicit) then 
         if (.not. allocated(phi_source_ffs)) allocate (phi_source_ffs(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         phi_source_ffs = 0.0
         if (.not. allocated(fields_source_ffs)) allocate (fields_source_ffs(naky, nakx, -nzgrid:nzgrid, ntubes))
         fields_source_ffs = 0.0
         if (.not. allocated(drifts_source_ffs)) allocate (drifts_source_ffs(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         drifts_source_ffs = 0.0
      end if
         
      if (proc0) call time_message(.false., time_implicit_advance(:, 1), ' Implicit time advance')

      allocate (phi_source(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (apar_source(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (bpar_source(naky, nakx, -nzgrid:nzgrid, ntubes))

      ! save the incoming pdf and fields, as they will be needed later
      allocate (phi_old(naky, nakx, -nzgrid:nzgrid, ntubes))
      phi_old = phi
      allocate (apar_old(naky, nakx, -nzgrid:nzgrid, ntubes))
      apar_old = apar
      allocate (bpar_old(naky, nakx, -nzgrid:nzgrid, ntubes))
      bpar_old = bpar

      !> dist_choice indicates whether the non-Boltzmann part of the pdf (h) is evolved
      !> in parallel streaming or if the guiding centre distribution (g = <f>) is evolved
      dist_choice = 'g'
      
      !> store g_n from pervious time step 
      g1 = g
      g2 = g

      call advance_fields(g2, phi, apar, bpar, dist=trim(dist_choice))
      phi_old = phi

      !> if using delphi formulation for response matrix, then phi = phi^n replaces
      !> phi^{n+1} in the inhomogeneous GKE; else set phi_{n+1} to zero in inhomogeneous equation
      !> NB: for FFS phi_source = 0.0 for inhomogeneous step - this is always set below
      !> for fluxtube this ordering doesn't matter
      if (use_deltaphi_for_response_matrix) then
         phi_source = phi
         if (include_bpar) bpar_source = bpar
      else
         phi_source = tupwnd_m * phi
         if (include_bpar) bpar_source = tupwnd_m * bpar
      end if
         
      !!> until fixed
      itt = 1
      do while (itt <= nitt)
         !> save the incoming pdf and phi, as they will be needed later
         !> this will become g^{n+1, i} -- the g from the previous iteration         
         if (driftkinetic_implicit .and. (itt .NE. 1)) then
            call advance_fields(g2, phi, apar, bpar, dist=trim(dist_choice))
            phi_old = phi
         end if 

         if (include_apar) then
            ! when solving for the 'inhomogeneous' piece of the pdf,
            ! use part of apar weighted towards previous time level
            apar_source = tupwnd_m * apar
            ! set apar=0, as in update_pdf it is used as the contribution from
            ! apar^{n+1}, which should not be part of the 'inhomogeneous' GKE eqn
            apar = 0.0
         end if

         !> if using delphi formulation for response matrix, then phi = phi^n replaces
         !> phi^{n+1} in the inhomogeneous GKE; else set phi_{n+1} to zero in inhomogeneous equation
         ! solve for the 'inhomogeneous' piece of the pdf
         if (driftkinetic_implicit) then
            call get_source_ffs_itteration (phi_old, g2, phi_source_ffs)
!!!!!!!            call get_drifts_ffs_itteration (phi_old, g2, drifts_source_ffs)
!!            phi_source_ffs = phi_source_ffs + drifts_source_ffs
            phi_source = 0.0
            !> set the g on the RHS to be g from the previous time step  
            !> FFS will have a RHS source term
            !> modify being passes in will make sure that this source is included
            modify = .true.
            call update_pdf(modify)
         else
            call update_pdf
         end if

         !> We now have g_{inh}^{n+1, i+1} stored in g
         !> calculate associated fields (phi_{inh}^{n+1, i+1})
         fields_updated = .false.
         if (driftkinetic_implicit) then
            !> For FFS we want to re-solve for bar{phi}
            !> NB the 'g' here is g_inh^{n+1, i+1}
            call advance_fields(g, phi, apar, bpar, dist=trim(dist_choice), implicit_solve=.true.) 
            !> g2 = g^{n+1, i}
            !> phi_old = phi^{n+1, i} 
            call get_fields_source(g2, phi_old, fields_source_ffs) 
            phi = phi + fields_source_ffs
         else
            call advance_fields(g, phi, apar, bpar, dist=trim(dist_choice)) 
         end if

         !> solve response_matrix*(phi^{n+1}-phi^{n*}) = phi_{inh}^{n+1}-phi^{n*}
         !> phi = phi_{inh}^{n+1}-phi^{n*} is input and overwritten by phi = phi^{n+1}-phi^{n*}
         if (use_deltaphi_for_response_matrix) phi = phi - phi_old
         if (use_deltaphi_for_response_matrix .and. include_bpar) bpar = bpar - bpar_old
         if (proc0) call time_message(.false., time_implicit_advance(:, 3), ' (back substitution)')
         !> for Drift kinetic implicit this is full phi^{n+1, i+1}
         call invert_parstream_response(phi, apar, bpar)
         if (proc0) call time_message(.false., time_implicit_advance(:, 3), ' (back substitution)')
         
         !> If using deltaphi formulation, must account for fact that phi = phi^{n+1}-phi^{n*}, but
         !> tupwnd_p should multiply phi^{n+1}
         if (use_deltaphi_for_response_matrix) phi = phi + phi_old
         if (use_deltaphi_for_response_matrix .and. include_bpar) bpar = bpar + bpar_old

         if (driftkinetic_implicit) then
            phi_source = phi
         else
	    ! get time-centered phi
            phi_source = tupwnd_m * phi_old + tupwnd_p * phi
            ! get time-centered apar
            if (include_apar) apar_source = tupwnd_m * apar_old + tupwnd_p * apar
            ! get time-centered bpar
            if (include_bpar) bpar_source = tupwnd_m * bpar_old + tupwnd_p * bpar
         end if

         ! solve for the final, updated pdf now that we have phi^{n+1}.
         if (driftkinetic_implicit) then
            !> Pass in modify to include RHS source term
            !> solving for full g = g_{inh} + g_{hom} 
            modify = .true.
            call update_pdf(modify)
            g2 = g 
         else
            call update_pdf
         end if
         
         !! change
         !!error = 0.0 
         itt = itt + 1
      end do

      deallocate (phi_old, apar_old, bpar_old)
      deallocate (phi_source, apar_source, bpar_source)
      if(driftkinetic_implicit) deallocate (fields_source_ffs) 
      if(driftkinetic_implicit) deallocate (phi_source_ffs, drifts_source_ffs) 

      if (proc0) call time_message(.false., time_implicit_advance(:, 1), ' Stream advance')

   contains

      subroutine update_pdf(mod)

         use extended_zgrid, only: neigen

         integer :: ie, it, iky, ivmu
         integer :: ulim
         complex, dimension(:), allocatable :: pdf1, pdf2
         complex, dimension(:), allocatable :: phiext, bparext
         complex, dimension(:), allocatable :: aparext, aparext_new, aparext_old
         complex, dimension(:), allocatable :: phiffsext

         logical, optional, intent(in) :: mod

         ! start the timer for the pdf update
         if (proc0) call time_message(.false., time_implicit_advance(:, 2), ' (bidiagonal solve)')

         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            ! solve for the pdf, given the sources for phi and the pdf on the RHS of the GK equation
            ! we do this on set of connected zed segments at a time
            do iky = 1, naky
               do it = 1, ntubes
                  do ie = 1, neigen(iky)
                     ! nz_ext is the number of grid points in the extended zed domain
                     nz_ext = nsegments(ie, iky) * nzed_segment + 1
                     ! pdf1 and pdf2 will be scratch arrays needed to compute the pdf itself,
                     ! as well as contributions to the GK equation
                     allocate (pdf1(nz_ext), pdf2(nz_ext))
                     ! phiext should contain the appropriate contribution to the time-centred phi;
                     ! for the 'inhomogeneous' GKE, it should have time_upwind_minus * phi^{n};
                     ! for the 'homogeneous' GKE, it should have time_upwind_plus * phi^{n+1};
                     ! and for the full GKE, it should be the sum of these two
                     allocate (phiext(nz_ext))
                     ! bpar is treated like phi above MRH
                     allocate (bparext(nz_ext))
                     ! if advancing apar, aparext should contain the appropriate contribution
                     ! to the time-centred apar;
                     ! for the 'inhomogeneous' GKE, it should have time_upwind_minus * apar^{n};
                     ! for the 'homogeneous' GKE, it should have time_upwind_plus * apar^{n+1};
                     ! and for the full GKE, it should be the sum of these two
                     allocate (aparext(nz_ext)); aparext = 0.0
                     ! if advancing apar, aparext_new should be zero if advancing the 'inhomogeneous'
                     ! piece of g or apar^{n+1} otherwise
                     allocate (aparext_new(nz_ext)); aparext_new = 0.0
                     ! if advancing apar, aparext_old will contain the apar originally passed into
                     ! the implicit time advance; needed to convert from g^{n} to gbar^{n}
                     allocate (aparext_old(nz_ext)); aparext_old = 0.0
                     ! map the incoming pdf 'g1' onto the extended zed domain and call it 'pdf1'
                     call map_to_extended_zgrid(it, ie, iky, g1(iky, :, :, :, ivmu), pdf1, ulim)
                     ! map the incoming potential 'phi_source' onto the extended zed domain and call it 'phiext'
                     call map_to_extended_zgrid(it, ie, iky, phi_source(iky, :, :, :), phiext, ulim)
                     ! map incoming parallel magnetic vector potetial 'apar_sosurce' onto
                     ! extended zed domain and call 'aparext'
                     if (include_apar) then
                        call map_to_extended_zgrid(it, ie, iky, apar_source(iky, :, :, :), aparext, ulim)
                        call map_to_extended_zgrid(it, ie, iky, apar(iky, :, :, :), aparext_new, ulim)
                        call map_to_extended_zgrid(it, ie, iky, apar_old(iky, :, :, :), aparext_old, ulim)
                     end if
                     ! map incoming bpar "bpar_source" onto the extended zed domain and call it "bparext"
                     if (include_bpar) call map_to_extended_zgrid(it, ie, iky, bpar_source(iky, :, :, :), bparext, ulim) 
                     if (present(mod)) then
                        !> For implicit FFS - Need to map the incoming source term on the RHS
                        !> (i.e. the piece that is treated explicitly) 
                        allocate (phiffsext(nz_ext))
                        call map_to_extended_zgrid(it, ie, iky, phi_source_ffs(iky, :, :, :, ivmu), phiffsext, ulim)
                        call get_gke_rhs(ivmu, iky, ie, pdf1, phiext, aparext, aparext_new, aparext_old, bparext, pdf2, phiffsext)
                     else
                        ! calculate the RHS of the GK equation (using pdf1 and phi_source as the
                        ! pdf and potential, respectively) and store it in pdf2
                        call get_gke_rhs(ivmu, iky, ie, pdf1, phiext, aparext, aparext_new, aparext_old, bparext, pdf2)
                     end if
                     ! given the RHS of the GK equation (pdf2), solve for the pdf at the
                     ! new time level by sweeping in zed on the extended domain;
                     ! the rhs is input as 'pdf2' and over-written with the updated solution for the pdf
                     call sweep_g_zext(iky, ie, it, ivmu, pdf2)
                     ! map the pdf 'pdf2' from the extended zed domain
                     ! to the standard zed domain; the mapped pdf is called 'g'
                     call map_from_extended_zgrid(it, ie, iky, pdf2, g(iky, :, :, :, ivmu))
                     deallocate (pdf1, pdf2, phiext, aparext, aparext_new, aparext_old, bparext)
                     if (present(mod)) deallocate(phiffsext)
                  end do
               end do
            end do
         end do

         ! stop the timer for the pdf update
         if (proc0) call time_message(.false., time_implicit_advance(:, 2), ' (bidiagonal solve)')

      end subroutine update_pdf

   end subroutine advance_implicit_terms

   !> get_gke_rhs calculates the RHS of the GK equation.
   !> as the response matrix approach requires separate solution of the 'inhomogeneous' GKE,
   !> the homogeneous GKE (to obtain the response matrix itself),
   !> and the full GKE, which RHS is obtained depends on the input values
   !> for 'pdf', 'phi', 'apar', 'aparnew' and 'aparold'
   subroutine get_gke_rhs(ivmu, iky, ie, pdf, phi, apar, aparnew, aparold, bpar, rhs, phi_ffs)

      implicit none

      integer, intent(in) :: ivmu, iky, ie
      complex, dimension(:), intent(in) :: pdf
      complex, dimension(:), intent(in) :: phi, apar, aparnew, aparold, bpar
      complex, dimension(:), intent(out) :: rhs

      integer :: nz_ext
      complex, dimension(:), allocatable :: rhs_fields

      complex, dimension(:), optional, intent(in) :: phi_ffs

      ! obtain the RHS of the GK eqn for given fields

      ! nz_ext is the number of grid points in this extended zed domain
      nz_ext = size(pdf)
      ! rhs_fields will be the contribution to the GKE RHS from the given fields
      allocate (rhs_fields(nz_ext))

      ! NB: rhs is used as a scratch array in get_contributions_from_fields
      ! so be careful not to move get_contributions_from_pdf before it, or rhs will be over-written
      call get_contributions_from_fields(phi, apar, aparnew, bpar, ivmu, iky, ie, rhs, rhs_fields)

      if (present(phi_ffs)) then
         call get_contributions_from_pdf(pdf, aparold, ivmu, iky, ie, rhs, phi_ffs)
      else
         call get_contributions_from_pdf(pdf, aparold, ivmu, iky, ie, rhs) 
      end if

      ! construct RHS of GK eqn
      rhs = rhs + rhs_fields

      deallocate (rhs_fields)

   end subroutine get_gke_rhs

   !> get_contributions_from_fields takes as input the appropriately averaged
   !> electrostatic potential phi and magnetic vector potential components apar
   !> and returns in rhs the sum of the source terms
   !> involving phi and apar that appear on the RHS of the GK equation when g is the pdf
   subroutine get_contributions_from_fields(phi, apar, aparnew, bpar, ivmu, iky, ie, scratch, rhs)

      use parameters_physics, only: include_apar, include_bpar
      use extended_zgrid, only: map_to_iz_ikx_from_izext

      implicit none

      complex, dimension(:), intent(in) :: phi, apar, aparnew, bpar
      integer, intent(in) :: ivmu, iky, ie
      complex, dimension(:), intent(out) :: scratch, rhs

      integer, dimension(:), allocatable :: iz_from_izext, ikx_from_izext
      integer :: nz_ext

      ! nz_ext is the number of grid points in the extended zed domain
      nz_ext = size(phi)

      ! determine the mapping from the extended domain zed index (izext) to the
      ! zed and kx domain indices (iz, ikx)
      allocate (iz_from_izext(nz_ext))
      allocate (ikx_from_izext(nz_ext))
      call map_to_iz_ikx_from_izext(iky, ie, iz_from_izext, ikx_from_izext)

      ! calculate the contributions to the RHS of the GKE due to source terms proportional to phi
      call get_contributions_from_phi(phi, ivmu, iky, iz_from_izext, ikx_from_izext, scratch, rhs)
      ! if advancing apar, get its contribution to the RHS of the GKE and add to phi contribution
      if (include_apar) then
         call get_contributions_from_apar(apar, aparnew, ivmu, iky, iz_from_izext, ikx_from_izext, scratch, rhs)
      end if
      ! if advancing bpar, get its contribution to the RHS of the GKE and add to phi contribution
      if (include_bpar) then 
         call get_contributions_from_bpar(bpar, ivmu, iky, iz_from_izext, ikx_from_izext, scratch, rhs)
      end if

      deallocate (iz_from_izext, ikx_from_izext)

   end subroutine get_contributions_from_fields

   !> get_contributions_from_phi takes as input the appropriately averaged
   !> electrostatic potential phi and returns in rhs the sum of the source terms
   !> involving phi that appear on the RHS of the GK equation when g is the pdf
   subroutine get_contributions_from_phi(phi, ivmu, iky, iz_from_izext, ikx_from_izext, scratch, rhs)

      use stella_time, only: code_dt
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac, maxwell_mu_avg
      use vpamu_grids, only: vpa
      use parameters_numerical, only: driftkinetic_implicit, maxwellian_normalization
      use parameters_numerical, only: maxwellian_inside_zed_derivative
      use parameters_numerical, only: drifts_implicit
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use neoclassical_terms, only: include_neoclassical_terms
      use neoclassical_terms, only: dfneo_dvpa

      implicit none

      complex, dimension(:), intent(in) :: phi
      integer, intent(in) :: ivmu, iky
      integer, dimension(:), intent(in) :: iz_from_izext, ikx_from_izext
      complex, dimension(:), intent(out) :: scratch, rhs

      real, dimension(:), allocatable :: z_scratch
      integer :: ia, iz, iv, imu, is
      integer :: nz_ext

      ia = 1
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      ! nz_ext is the number of grid points in the extended zed domain
      nz_ext = size(phi)

      ! allocate a 1d array in zed for use as a scratch array
      allocate (z_scratch(-nzgrid:nzgrid))

      ! set scratc to be phi or <phi> depending on whether parallel streaming is
      ! implicit or only implicit in the kperp = 0 (drift kinetic) piece
      if (driftkinetic_implicit) then
         scratch = phi
      else
         !> get <phi>
         call gyro_average_zext(iky, ivmu, ikx_from_izext, iz_from_izext, phi, scratch)
      end if

      call add_streaming_contribution_phi
      if (drifts_implicit) call add_drifts_contribution_phi

      deallocate (z_scratch)

   contains

      subroutine add_streaming_contribution_phi

         use extended_zgrid, only: fill_zext_ghost_zones
         use parallel_streaming, only: get_zed_derivative_extended_domain
         use parallel_streaming, only: center_zed
         use parallel_streaming, only: gradpar_c, stream_sign

         use parameters_numerical, only: driftkinetic_implicit

         integer :: izext
         complex :: scratch_left, scratch_right

         ! fill ghost zones beyond ends of extended zed domain for <phi>
         ! and store values in scratch_left and scratch_right
         call fill_zext_ghost_zones(iky, scratch, scratch_left, scratch_right)

         ! obtain the zed derivative of <phi> (stored in scratch) and store in rhs
         call get_zed_derivative_extended_domain(iv, scratch, scratch_left, scratch_right, rhs)

         ! center Maxwellian factor in mu
         ! and store in dummy variable z_scratch
         if(driftkinetic_implicit) then
            z_scratch = maxwell_mu_avg(ia, :, imu, is)
            call center_zed(iv, z_scratch, -nzgrid)
         else
            if (.not. maxwellian_normalization) then
               ! center Maxwellian factor in mu
               ! and store in dummy variable z_scratch
               z_scratch = maxwell_mu(ia, :, imu, is)
               call center_zed(iv, z_scratch, -nzgrid)
            else
               z_scratch = 1.0
            end if
         end if
         ! multiply by Maxwellian factor
         do izext = 1, nz_ext
            rhs(izext) = rhs(izext) * z_scratch(iz_from_izext(izext))
         end do         

         ! NB: could do this once at beginning of simulation to speed things up
         ! this is vpa*Z/T*exp(-vpa^2)
         z_scratch = vpa(iv) * spec(is)%zt
         if (.not. maxwellian_normalization) z_scratch = z_scratch * maxwell_vpa(iv, is) * maxwell_fac(is)
         ! if including neoclassical correction to equilibrium distribution function
         ! then must also account for -vpa*dF_neo/dvpa*Z/T
         ! CHECK TO ENSURE THAT DFNEO_DVPA EXCLUDES EXP(-MU*B/T) FACTOR !!
         if (include_neoclassical_terms) then
            do iz = -nzgrid, nzgrid
               z_scratch(iz) = z_scratch(iz) - 0.5 * dfneo_dvpa(ia, iz, ivmu) * spec(is)%zt
            end do
            call center_zed(iv, z_scratch, -nzgrid)
         end if

         if (stream_sign(iv) > 0) then
            z_scratch = z_scratch * gradpar_c(:, -1) * code_dt * spec(is)%stm_psi0
         else
            z_scratch = z_scratch * gradpar_c(:, 1) * code_dt * spec(is)%stm_psi0
         end if

         do izext = 1, nz_ext
            rhs(izext) = -z_scratch(iz_from_izext(izext)) * rhs(izext)
         end do

      end subroutine add_streaming_contribution_phi

      subroutine add_drifts_contribution_phi

         use constants, only: zi
         use grids_kxky, only: aky, akx
         use arrays_dist_fn, only: wstar, wdriftx_phi, wdrifty_phi
         use parallel_streaming, only: center_zed
         use extended_zgrid, only: periodic

         integer :: izext, iz, ikx

         ! 'scratch' starts out as the gyro-average of phi, evaluated at zed grid points
         do izext = 1, nz_ext
            ikx = ikx_from_izext(izext)
            iz = iz_from_izext(izext)
            scratch(izext) = zi * scratch(izext) * (akx(ikx) * wdriftx_phi(ia, iz, ivmu) &
                                                    + aky(iky) * (wdrifty_phi(ia, iz, ivmu) + wstar(ia, iz, ivmu)))
         end do
         call center_zed(iv, scratch, 1, periodic(iky))

         rhs = rhs + scratch

      end subroutine add_drifts_contribution_phi

   end subroutine get_contributions_from_phi

   !> get_contributions_from_bpar takes as input the appropriately averaged
   !> electrostatic potential bpar and returns in rhs the sum of the source terms
   !> involving bpar that appear on the RHS of the GK equation when g is the pdf
   subroutine get_contributions_from_bpar(bpar, ivmu, iky, iz_from_izext, ikx_from_izext, scratch, rhs)

      use stella_time, only: code_dt
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: vpa, mu
      use parameters_numerical, only: driftkinetic_implicit, maxwellian_normalization
      use parameters_numerical, only: maxwellian_inside_zed_derivative
      use parameters_numerical, only: drifts_implicit
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use neoclassical_terms, only: include_neoclassical_terms
      use neoclassical_terms, only: dfneo_dvpa
      use extended_zgrid, only: map_to_iz_ikx_from_izext

      implicit none

      complex, dimension(:), intent(in) :: bpar
      integer, intent(in) :: ivmu, iky
      integer, dimension(:), intent(in) :: iz_from_izext, ikx_from_izext
      complex, dimension(:), intent(out) :: scratch, rhs
      complex, dimension(:), allocatable :: scratch2
      
      real, dimension(:), allocatable :: z_scratch
      integer :: ia, iz, iv, imu, is
      integer :: nz_ext

      ia = 1
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      ! nz_ext is the number of grid points in the extended zed domain
      nz_ext = size(bpar)
      
      ! allocate a 1d array in zed for use as a scratch array
      allocate (z_scratch(-nzgrid:nzgrid))
      ! allocate a 1d array to replace the rhs array as a scratch array
      allocate (scratch2(nz_ext))

      ! set scratch to be bpar or <bpar> depending on whether parallel streaming is
      ! implicit or only implicit in the kperp = 0 (drift kinetic) piece
      if (driftkinetic_implicit) then
         scratch = bpar
      else
         call gyro_average_j1_zext(iky, ivmu, ikx_from_izext, iz_from_izext, bpar, scratch)
      end if
      
      call add_streaming_contribution_bpar
      if (drifts_implicit) call add_drifts_contribution_bpar

      deallocate (z_scratch, scratch2)

   contains

      subroutine add_streaming_contribution_bpar

         use extended_zgrid, only: fill_zext_ghost_zones
         use parallel_streaming, only: get_zed_derivative_extended_domain
         use parallel_streaming, only: center_zed
         use parallel_streaming, only: gradpar_c, stream_sign

         integer :: izext
         complex :: scratch_left, scratch_right

         ! fill ghost zones beyond ends of extended zed domain for <bpar>
         ! and store values in scratch_left and scratch_right
         call fill_zext_ghost_zones(iky, scratch, scratch_left, scratch_right)

         ! obtain the zed derivative of <bpar> (stored in scratch) and store in scratch2
         call get_zed_derivative_extended_domain(iv, scratch, scratch_left, scratch_right, scratch2)

         if (.not. maxwellian_normalization) then
            ! center Maxwellian factor in mu
            ! and store in dummy variable z_scratch
            z_scratch = maxwell_mu(ia, :, imu, is)
            call center_zed(iv, z_scratch, -nzgrid)
            ! multiply by Maxwellian factor
            do izext = 1, nz_ext
               scratch2(izext) = scratch2(izext) * z_scratch(iz_from_izext(izext))
            end do
         end if

         ! NB: could do this once at beginning of simulation to speed things up
         ! this is vpa*Z/T*exp(-vpa^2)
         z_scratch = vpa(iv) * 4. * mu(imu)
         if (.not. maxwellian_normalization) z_scratch = z_scratch * maxwell_vpa(iv, is) * maxwell_fac(is)
         ! if including neoclassical correction to equilibrium distribution function
         ! then must also account for -vpa*dF_neo/dvpa*4*mu
         ! CHECK TO ENSURE THAT DFNEO_DVPA EXCLUDES EXP(-MU*B/T) FACTOR !!
         if (include_neoclassical_terms) then
            do iz = -nzgrid, nzgrid
               z_scratch(iz) = z_scratch(iz) - 0.5 * dfneo_dvpa(ia, iz, ivmu) * 4. * mu(imu)
            end do
            call center_zed(iv, z_scratch, -nzgrid)
         end if

         if (stream_sign(iv) > 0) then
            z_scratch = z_scratch * gradpar_c(:, -1) * code_dt * spec(is)%stm_psi0
         else
            z_scratch = z_scratch * gradpar_c(:, 1) * code_dt * spec(is)%stm_psi0
         end if

         do izext = 1, nz_ext
            scratch2(izext) = -z_scratch(iz_from_izext(izext)) * scratch2(izext)
         end do
         ! add scratch2 to rhs
         rhs = rhs + scratch2
      end subroutine add_streaming_contribution_bpar

      subroutine add_drifts_contribution_bpar

         use constants, only: zi
         use grids_kxky, only: aky, akx
         use arrays_dist_fn, only: wstar, wdriftx_bpar, wdrifty_bpar
         use parallel_streaming, only: center_zed
         use extended_zgrid, only: periodic

         integer :: izext, iz, ikx
         real :: constant_factor
         ! 'scratch' starts out as the gyro-average of bpar, evaluated at zed grid points
         constant_factor = 4. * mu(imu) * spec(is)%tz
         do izext = 1, nz_ext
            ikx = ikx_from_izext(izext)
            iz = iz_from_izext(izext)
            ! the bpar part of Zs <chi>/Ts = 4 mu J1 bpar / bs, and wdrifty_bpar and wdriftx_bpar contain the 4 mu factor
            ! the 4 mu Ts/Zs factor must be included explicitly in the wstar term here 
            scratch(izext) = zi * scratch(izext) * (akx(ikx) * wdriftx_bpar(ia, iz, ivmu) &
                                                    + aky(iky) * (wdrifty_bpar(ia, iz, ivmu) + constant_factor * wstar(ia, iz, ivmu)))
         end do
         call center_zed(iv, scratch, 1, periodic(iky))

         rhs = rhs + scratch

      end subroutine add_drifts_contribution_bpar

   end subroutine get_contributions_from_bpar

   !> get_contributions_from_apar takes as input the appropriately averaged
   !> parallel component of the vector potential, apar, and returns in rhs the sum of the source terms
   !> involving apar that appear on the RHS of the GK equation when g is the pdf
   subroutine get_contributions_from_apar(apar, aparnew, ivmu, iky, iz_from_izext, ikx_from_izext, scratch, rhs)

      use parameters_numerical, only: driftkinetic_implicit, drifts_implicit
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

      implicit none

      complex, dimension(:), intent(in) :: apar, aparnew
      integer, intent(in) :: ivmu, iky
      integer, dimension(:), intent(in) :: iz_from_izext, ikx_from_izext
      complex, dimension(:), intent(out) :: scratch
      complex, dimension(:), intent(in out) :: rhs

      complex, dimension(:), allocatable :: scratch2
      integer :: ia, iv, imu, is
      integer :: nz_ext

      ia = 1
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      ! nz_ext is number of grid points in extended zed domain
      nz_ext = size(scratch)

      allocate (scratch2(nz_ext))

      ! set scratch to be apar or <apar> depending on whether parallel streaming is
      ! implicit or only implicit in the kperp = 0 (drift kinetic) piece
      if (driftkinetic_implicit) then
         scratch = apar
         scratch2 = aparnew
      else
         call gyro_average_zext(iky, ivmu, ikx_from_izext, iz_from_izext, apar, scratch)
         call gyro_average_zext(iky, ivmu, ikx_from_izext, iz_from_izext, aparnew, scratch2)
      end if

      call add_gbar_to_g_contribution_apar(scratch2, iky, ia, iv, imu, is, nz_ext, iz_from_izext, rhs)
      if (drifts_implicit) call add_drifts_contribution_apar(scratch, iky, ia, ivmu, iv, is, nz_ext, iz_from_izext, rhs)

      deallocate (scratch2)

   end subroutine get_contributions_from_apar

   !> adds the contributions to the GKE RHS that comes from switching from
   !> gbar^{n+1} = g^{n+1} + (Ze/T)*(vpa/c)*<Apar^{n+1}>*F0 to g^{n+1} = <f^{n+1}>
   !> in the time derivative;
   ! as it involves apar^{n+1}, it should not be present in the solution of the
   ! 'inhomogeneous' GKE; this should have been accounted for by passing in
   ! aparnew=0 so that scratch2 will be zero below
   subroutine add_gbar_to_g_contribution_apar(scratch2, iky, ia, iv, imu, is, nz_ext, iz_from_izext, rhs)

      use parameters_numerical, only: maxwellian_normalization
      use vpamu_grids, only: vpa, maxwell_vpa, maxwell_mu, maxwell_fac
      use parallel_streaming, only: center_zed
      use species, only: spec
      use extended_zgrid, only: periodic

      implicit none

      complex, dimension(:), intent(in out) :: scratch2
      integer, intent(in) :: iky, ia, iv, imu, is, nz_ext
      integer, dimension(:), intent(in) :: iz_from_izext
      complex, dimension(:), intent(in out) :: rhs

      integer :: izext, iz
      real :: constant_factor

      ! avoid repeated multiplication in below izext loop
      constant_factor = -2.0 * spec(is)%zt_psi0 * spec(is)%stm_psi0 * vpa(iv)

      ! incoming 'scratch2' is <apar^{n+1}>
      do izext = 1, nz_ext
         iz = iz_from_izext(izext)
         scratch2(izext) = constant_factor * scratch2(izext)
      end do

      ! if the pdf is not normalized by a Maxwellian then the source term contains a Maxwellian factor
      if (.not. maxwellian_normalization) then
         do izext = 1, nz_ext
            iz = iz_from_izext(izext)
            scratch2(izext) = scratch2(izext) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
         end do
      end if
      call center_zed(iv, scratch2, 1, periodic(iky))
      rhs = rhs + scratch2

   end subroutine add_gbar_to_g_contribution_apar

   subroutine add_drifts_contribution_apar(scratch, iky, ia, ivmu, iv, is, nz_ext, iz_from_izext, rhs)

      use constants, only: zi
      use species, only: spec
      use grids_kxky, only: aky
      use arrays_dist_fn, only: wstar
      use parallel_streaming, only: center_zed
      use extended_zgrid, only: periodic
      use vpamu_grids, only: vpa

      implicit none

      complex, dimension(:), intent(in out) :: scratch, rhs
      integer, intent(in) :: iky, ia, ivmu, iv, is, nz_ext
      integer, dimension(:), intent(in) :: iz_from_izext

      integer :: izext, iz
      complex :: constant_factor

      constant_factor = -2.0 * zi * spec(is)%stm_psi0 * vpa(iv) * aky(iky)

      do izext = 1, nz_ext
         iz = iz_from_izext(izext)
         scratch(izext) = constant_factor * scratch(izext) * wstar(ia, iz, ivmu)
      end do
      call center_zed(iv, scratch, 1, periodic(iky))
      rhs = rhs + scratch

   end subroutine add_drifts_contribution_apar

   subroutine gbar_to_g_zext(pdf, apar, facapar, iky, ivmu, ikx_from_izext, iz_from_izext)

      use species, only: spec
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use parameters_numerical, only: maxwellian_normalization
      use vpamu_grids, only: vpa, maxwell_vpa, maxwell_mu, maxwell_fac

      implicit none

      complex, dimension(:), intent(in out) :: pdf
      integer, intent(in) :: ivmu, iky
      integer, dimension(:), intent(in) :: ikx_from_izext, iz_from_izext
      complex, dimension(:), intent(in) :: apar
      real, intent(in) :: facapar

      integer :: iv, imu, is
      integer :: izext, iz, ia
      integer :: nz_ext

      complex, dimension(:), allocatable :: field, gyro_field

      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      nz_ext = size(apar)

      allocate (field(nz_ext))
      allocate (gyro_field(nz_ext))

      ia = 1

      field = 2.0 * spec(is)%zt * spec(is)%stm_psi0 * vpa(iv) * facapar * apar
      if (.not. maxwellian_normalization) then
         do izext = 1, nz_ext
            iz = iz_from_izext(izext)
            field(izext) = field(izext) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
         end do
      end if
      call gyro_average_zext(iky, ivmu, ikx_from_izext, iz_from_izext, field, gyro_field)
      pdf = pdf - gyro_field

      deallocate (field, gyro_field)

   end subroutine gbar_to_g_zext

   subroutine gyro_average_zext(iky, ivmu, ikx_from_izext, iz_from_izext, fld, gyro_fld)

      use gyro_averages, only: gyro_average

      implicit none

      integer, intent(in) :: iky, ivmu
      integer, dimension(:), intent(in) :: ikx_from_izext, iz_from_izext
      complex, dimension(:), intent(in) :: fld
      complex, dimension(:), intent(out) :: gyro_fld

      integer :: izext, nz_ext

      nz_ext = size(fld)

      do izext = 1, nz_ext
         call gyro_average(fld(izext), iky, ikx_from_izext(izext), iz_from_izext(izext), ivmu, gyro_fld(izext))
      end do

   end subroutine gyro_average_zext

   subroutine gyro_average_j1_zext(iky, ivmu, ikx_from_izext, iz_from_izext, fld, gyro_fld)

      use gyro_averages, only: gyro_average_j1

      implicit none

      integer, intent(in) :: iky, ivmu
      integer, dimension(:), intent(in) :: ikx_from_izext, iz_from_izext
      complex, dimension(:), intent(in) :: fld
      complex, dimension(:), intent(out) :: gyro_fld

      integer :: izext, nz_ext

      nz_ext = size(fld)

      do izext = 1, nz_ext
         call gyro_average_j1(fld(izext), iky, ikx_from_izext(izext), iz_from_izext(izext), ivmu, gyro_fld(izext))
      end do

   end subroutine gyro_average_j1_zext

   !> get_contributions_from_pdf takes as an argument the evolved pdf
   !> (either guiding centre distribution g=<f> or maxwellian-normlized, non-Boltzmann distribution h/F0=f/F0+(Ze*phi/T))
   !> and the scratch array rhs, and returns the source terms that depend on the pdf in rhs
   subroutine get_contributions_from_pdf(pdf, apar, ivmu, iky, ie, rhs, source_ffs)

      use constants, only: zi
      use stella_time, only: code_dt
      use parameters_physics, only: include_apar
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use grids_kxky, only: aky, akx
      use vpamu_grids, only: vpa
      use stella_layouts, only: vmu_lo, iv_idx, is_idx
      use parameters_numerical, only: time_upwind_minus
      use parameters_numerical, only: drifts_implicit
      use parallel_streaming, only: get_zed_derivative_extended_domain, center_zed
      use parallel_streaming, only: gradpar_c, stream_sign
      use arrays_dist_fn, only: wdriftx_g, wdrifty_g
      use extended_zgrid, only: fill_zext_ghost_zones
      use extended_zgrid, only: map_to_iz_ikx_from_izext
      use extended_zgrid, only: periodic

      implicit none

      complex, dimension(:), intent(in) :: pdf, apar
      integer, intent(in) :: ivmu, iky, ie
      complex, dimension(:), intent(out) :: rhs

      real, dimension(:), allocatable :: gradpar_fac
      complex, dimension(:), allocatable :: dpdf_dz
      real :: constant_factor
      integer :: iv, is, iz
      integer :: ia, ikx

      integer, dimension(:), allocatable :: iz_from_izext, ikx_from_izext
      integer :: nz_ext, izext
      complex :: pdf_left, pdf_right

      complex, dimension(:), optional, intent(in) :: source_ffs

      ia = 1
      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      ! nz_ext is the number of grid points in the extended zed domain
      nz_ext = size(pdf)

      ! determine the mapping from the extended domain zed index (izext) to the
      ! zed and kx domain indices (iz, ikx)
      allocate (iz_from_izext(nz_ext))
      allocate (ikx_from_izext(nz_ext))
      call map_to_iz_ikx_from_izext(iky, ie, iz_from_izext, ikx_from_izext)

      ! fill ghost zones beyond ends of extended zed domain for the pdf
      ! and store values in scratch_left and scratch_right
      call fill_zext_ghost_zones(iky, pdf, pdf_left, pdf_right)

      ! obtain the zed derivative of the pdf and store in dpdf_dz
      allocate (dpdf_dz(nz_ext))
      call get_zed_derivative_extended_domain(iv, pdf, pdf_left, pdf_right, dpdf_dz)

      ! compute the z-independent factor appearing in front of the d(pdf)/dz term on the RHS of the Gk equation
      constant_factor = -code_dt * spec(is)%stm_psi0 * vpa(iv) * time_upwind_minus
      ! use the correctly centred (b . grad z) pre-factor for this sign of vpa
      allocate (gradpar_fac(-nzgrid:nzgrid))
      if (stream_sign(iv) > 0) then
         gradpar_fac = gradpar_c(:, -1) * constant_factor
      else
         gradpar_fac = gradpar_c(:, 1) * constant_factor
      end if

      rhs = pdf
      ! if advancing apar, need to use gbar rather than g=<f> for part of source on RHS of GKE,
      ! so convert from g to gbar
      if (include_apar) then
         call gbar_to_g_zext(rhs, apar, -1.0, iky, ivmu, ikx_from_izext, iz_from_izext)
      end if

      if (drifts_implicit) then
         do izext = 1, nz_ext
            ikx = ikx_from_izext(izext)
            iz = iz_from_izext(izext)
            rhs(izext) = rhs(izext) + pdf(izext) * zi * time_upwind_minus &
                         * (wdriftx_g(ia, iz, ivmu) * akx(ikx) + wdrifty_g(ia, iz, ivmu) * aky(iky))
!            rhs(izext) = rhs(izext) * (1.0 + zi * time_upwind_minus &
!                                       * (wdriftx_g(ia, iz, ivmu) * akx(ikx) + wdrifty_g(ia, iz, ivmu) * aky(iky)))
         end do
      end if

      ! cell-center the terms involving the pdf
      call center_zed(iv, rhs, 1, periodic(iky))

      ! construct the source term on the RHS of the GK equation coming from
      ! the pdf evaluated at the previous time level
      if(present(source_ffs)) then 
         do izext = 1, nz_ext
            rhs(izext) = rhs(izext) + gradpar_fac(iz_from_izext(izext)) * dpdf_dz(izext) + source_ffs(izext)
         end do
      else
         do izext = 1, nz_ext
            rhs(izext) = rhs(izext) + gradpar_fac(iz_from_izext(izext)) * dpdf_dz(izext)
         end do
      end if

      deallocate (dpdf_dz)
      deallocate (gradpar_fac)
      deallocate (iz_from_izext, ikx_from_izext)

   end subroutine get_contributions_from_pdf

   subroutine sweep_g_zext(iky, ie, it, ivmu, pdf)

      use constants, only: zi
      use zgrid, only: nzgrid, ntubes, delzed
      use parameters_numerical, only: drifts_implicit
      use parameters_numerical, only: zed_upwind_plus, zed_upwind_minus
      use parameters_numerical, only: time_upwind_plus
      use parameters_kxky_grids, only: nakx
      use grids_kxky, only: akx, aky
      use arrays_dist_fn, only: wdriftx_g, wdrifty_g
      use extended_zgrid, only: map_to_extended_zgrid
      use extended_zgrid, only: periodic, phase_shift
      use parallel_streaming, only: stream_sign, stream_c
      use parallel_streaming, only: center_zed
      use stella_layouts, only: vmu_lo, iv_idx, is_idx

      implicit none

      integer, intent(in) :: iky, ie, it, ivmu
      complex, dimension(:), intent(in out) :: pdf

      complex, dimension(:), allocatable :: wdrift_ext, pdf_cf
      complex, dimension(:, :), allocatable :: wdrift
      complex :: wd_factor, fac1, phase_factor
      real :: zupwnd_p, zupwnd_m, tupwnd_p
      real :: stream_term
      integer :: iz, ikx, ia
      integer :: iv, is
      integer :: ulim, sgn, iz1, iz2

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      sgn = stream_sign(iv)
      ! avoid repeated calculation of constants
      zupwnd_p = 2.0 * zed_upwind_plus
      zupwnd_m = 2.0 * zed_upwind_minus
      tupwnd_p = 2.0 * time_upwind_plus

      ! if treating magentic drifts implicitly in time,
      ! get the drift frequency on the extended zed grid
      if (drifts_implicit) then
         allocate (wdrift(nakx, -nzgrid:nzgrid))
         ia = 1
         ! sum up the kx and ky contributions to the magnetic drift frequency
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               wdrift(ikx, iz) = -zi * (wdriftx_g(ia, iz, ivmu) * akx(ikx) + wdrifty_g(ia, iz, ivmu) * aky(iky))
            end do
         end do
         ! obtain the drift frequency on the extended zed domain
         allocate (wdrift_ext(size(pdf)))
         call map_to_extended_zgrid(it, ie, iky, spread(wdrift, 3, ntubes), wdrift_ext, ulim)
         ! NB: need to check if passing periodic(iky) is the right thing to do here
         call center_zed(iv, wdrift_ext, 1, periodic(iky))
      else
         ulim = size(pdf)
      end if
      ! determine the starting and ending indices for sweep over the extended zed grid.
      ! as we are using a zero-incoming BC, these indices depend on the sign of the advection velocity
      ! note that sgn < 0 actually corresponds to positive advection velocity
      if (sgn < 0) then
         iz1 = 1; iz2 = ulim
      else
         iz1 = ulim; iz2 = 1
      end if
      ! the case of periodic BC must be treated separately from the zero-incoming-BC case
      if (periodic(iky)) then
         ! to enforce periodicity, decompose the pdf into a particular integral
         ! and complementary function.

         ! calculate the particular integral, with zero BC, and store in pdf
         iz = sgn * nzgrid
         pdf(iz1) = 0.0
         call get_updated_pdf(iz, iv, is, sgn, iz1, iz2, wdrift_ext, pdf)
         ! calculate the complementary function, with unit BC, and store in pdf_cf
         allocate (pdf_cf(ulim))
         iz = sgn * nzgrid
         pdf_cf = 0.0; pdf_cf(iz1) = 1.0
         call get_updated_pdf(iz, iv, is, sgn, iz1, iz2, wdrift_ext, pdf_cf)
         ! construct pdf = pdf_PI + (pdf_PI(zend)/(1-pdf_CF(zend))) * pdf_CF
         phase_factor = phase_shift(iky)**(-sgn)
         pdf = pdf + (phase_factor * pdf(iz2) / (1.0 - phase_factor * pdf_cf(iz2))) * pdf_cf
         deallocate (pdf_cf)
      else
         ! specially treat the most upwind grid point
         iz = sgn * nzgrid
         wd_factor = 1.0
         if (drifts_implicit) wd_factor = 1.0 + 0.5 * tupwnd_p * wdrift_ext(iz1)
         stream_term = tupwnd_p * stream_c(iz, iv, is) / delzed(0)
         fac1 = zupwnd_p * wd_factor + sgn * stream_term
         pdf(iz1) = pdf(iz1) * 2.0 / fac1
         ! now that we have the pdf at the most upwind point, sweep over the
         ! rest of the extended zed domain to obtain the pdf(z)
         call get_updated_pdf(iz, iv, is, sgn, iz1, iz2, wdrift_ext, pdf)
      end if

      if (drifts_implicit) deallocate (wdrift, wdrift_ext)

   end subroutine sweep_g_zext

   subroutine get_updated_pdf(iz, iv, is, sgn, iz1, iz2, wdrift_ext, pdf)

      use zgrid, only: nzgrid, delzed
      use parameters_numerical, only: drifts_implicit
      use parameters_numerical, only: zed_upwind_plus, zed_upwind_minus
      use parameters_numerical, only: time_upwind_plus
      use parallel_streaming, only: stream_c

      implicit none

      integer, intent(in out) :: iz
      integer, intent(in) :: iv, is, sgn, iz1, iz2
      complex, dimension(:), intent(in) :: wdrift_ext
      complex, dimension(:), intent(in out) :: pdf

      integer :: izext
      real :: stream_term
      real :: tupwnd_p, zupwnd_p, zupwnd_m
      complex :: wd_factor, fac1, fac2

      tupwnd_p = 2.0 * time_upwind_plus
      zupwnd_p = 2.0 * zed_upwind_plus
      zupwnd_m = 2.0 * zed_upwind_minus

      ! wd_factor will be modified from below unity to account for magnetic drifts
      ! if the drifts are treated implicitly
      wd_factor = 1.0

      do izext = iz1 - sgn, iz2, -sgn
         if (iz == -sgn * nzgrid) then
            iz = sgn * nzgrid - sgn
         else
            iz = iz - sgn
         end if
         if (drifts_implicit) wd_factor = 1.0 + 0.5 * tupwnd_p * wdrift_ext(izext)
         stream_term = tupwnd_p * stream_c(iz, iv, is) / delzed(0)
         fac1 = zupwnd_p * wd_factor + sgn * stream_term
         fac2 = zupwnd_m * wd_factor - sgn * stream_term
         pdf(izext) = (-pdf(izext + sgn) * fac2 + 2.0 * pdf(izext)) / fac1
      end do

   end subroutine get_updated_pdf

   subroutine sweep_zed_zonal(iky, iv, is, sgn, g, llim)

      use zgrid, only: nzgrid, delzed
      use extended_zgrid, only: phase_shift
      use parameters_numerical, only: zed_upwind, time_upwind
      use parallel_streaming, only: stream_c

      implicit none

      integer, intent(in) :: iky, iv, is, sgn, llim
      complex, dimension(llim:), intent(in out) :: g

      integer :: iz, izext, iz1, iz2, npts, ulim
      real :: fac1, fac2
      complex :: pf
      complex, dimension(:), allocatable :: gcf, gpi

      npts = size(g)
      ulim = llim + npts - 1

      allocate (gpi(llim:ulim))
      allocate (gcf(llim:ulim))
      ! ky=0 is 2pi periodic (no extended zgrid)
      ! decompose into complementary function + particular integral
      ! zero BC for particular integral
      ! unit BC for complementary function (no source)
      if (sgn < 0) then
         iz1 = llim; iz2 = ulim
      else
         iz1 = ulim; iz2 = llim
      end if
      pf = phase_shift(iky)**(-sgn)
      gpi(iz1) = 0.; gcf(iz1) = 1.
      iz = sgn * nzgrid
      do izext = iz1 - sgn, iz2, -sgn
         if (iz == -sgn * nzgrid) then
            iz = sgn * nzgrid - sgn
         else
            iz = iz - sgn
         end if
         fac1 = 1.0 + zed_upwind + sgn * (1.0 + time_upwind) * stream_c(iz, iv, is) / delzed(0)
         fac2 = 1.0 - zed_upwind - sgn * (1.0 + time_upwind) * stream_c(iz, iv, is) / delzed(0)
         gpi(izext) = (-gpi(izext + sgn) * fac2 + 2.0 * g(izext)) / fac1
         gcf(izext) = -gcf(izext + sgn) * fac2 / fac1
      end do
      ! g = g_PI + (g_PI(pi)/(1-g_CF(pi))) * g_CF
      g = gpi + (pf * gpi(iz2) / (1.0 - pf * gcf(iz2))) * gcf
      deallocate (gpi, gcf)

   end subroutine sweep_zed_zonal

   !> use the LU-decomposed response matrix and the contributions from the
   !> 'inhomogeneous' fields (phi, apar) to solve for (phi^{n+1}, apar^{n+1})
   subroutine invert_parstream_response(phi, apar, bpar)

      use linear_solve, only: lu_back_substitution
      use parameters_physics, only: include_apar, include_bpar
      use zgrid, only: nzgrid, ntubes
      use extended_zgrid, only: neigen
      use extended_zgrid, only: nsegments
      use extended_zgrid, only: nzed_segment
      use extended_zgrid, only: map_to_extended_zgrid
      use extended_zgrid, only: map_from_extended_zgrid
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: periodic, phase_shift
      use parameters_kxky_grids, only: naky
      use fields, only: nfields
      use arrays_fields, only: response_matrix

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar

      integer :: iky, ie, it, ulim
      integer :: ikx
      integer :: offset_apar, offset_bpar
      integer :: nresponse_per_field, nresponse, nzp
      complex, dimension(:), allocatable :: fld_ext, phi_ext, apar_ext, bpar_ext
      complex, dimension(:), allocatable :: fld

      ! put the fields onto the extended zed grid and use LU back substitution
      do iky = 1, naky
         ! avoid double counting of periodic endpoints for zonal (and any other periodic) modes
         if (periodic(iky)) then
            nzp = 2 * nzgrid
            allocate (fld(nzp * nfields))
            ! set offset integers for array slices involving apar and bpar
            if (include_apar) then
               offset_apar = nzp
            else
               offset_apar = 0
            end if
            if (include_bpar) offset_bpar = offset_apar + nzp
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  ikx = ikxmod(1, ie, iky)
                  ! construct the field vector, consisting of phi (and apar if evolving)
                  fld(:nzp) = phi(iky, ikx, -nzgrid:nzgrid - 1, it)
                  if (include_apar) fld(offset_apar + 1:nzp + offset_apar) = apar(iky, ikx, -nzgrid:nzgrid - 1, it)
                  if (include_bpar) fld(offset_bpar + 1:nzp + offset_bpar) = bpar(iky, ikx, -nzgrid:nzgrid - 1, it)
                  ! use LU back substitution to solve the linear response matrix system
                  call lu_back_substitution(response_matrix(iky)%eigen(ie)%zloc, &
                                            response_matrix(iky)%eigen(ie)%idx, fld)
                  ! unpack phi (and apar if evolving) from the field vector;
                  ! also apply phase shift at periodic point
                  phi(iky, ikx, -nzgrid:nzgrid - 1, it) = fld(:nzp)
                  phi(iky, ikx, nzgrid, it) = phi(iky, ikx, -nzgrid, it) / phase_shift(iky)
                  if (include_apar) then
                     apar(iky, ikx, -nzgrid:nzgrid - 1, it) = fld(offset_apar + 1:nzp + offset_apar)
                     apar(iky, ikx, nzgrid, it) = apar(iky, ikx, -nzgrid, it) / phase_shift(iky)
                  end if
                  if (include_bpar) then
                     bpar(iky, ikx, -nzgrid:nzgrid - 1, it) = fld(offset_bpar + 1:nzp + offset_bpar)
                     bpar(iky, ikx, nzgrid, it) = bpar(iky, ikx, -nzgrid, it) / phase_shift(iky)
                  end if
               end do
            end do
            deallocate (fld)
         else
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  ! solve response_matrix*(phi^{n+1}, apar^{n+1}) = (phi_{inh}^{n+1}, apar_{inh}^{n+1})
                  nresponse_per_field = nsegments(ie, iky) * nzed_segment + 1
                  nresponse = nresponse_per_field * nfields
                  ! fld_ext will contain phi or (phi, apar), depending on whether apar is evolved
                  allocate (fld_ext(nresponse))
                  ! set offset integers for array slices involving apar and bpar
                  if (include_apar) then
                     offset_apar = nresponse_per_field
                  else
                     offset_apar = 0
                  end if
                  if (include_bpar) offset_bpar = offset_apar + nresponse_per_field
                  
                  ! phi_ext contains phi on the extended zed domain
                  allocate (phi_ext(nresponse_per_field))
                  call map_to_extended_zgrid(it, ie, iky, phi(iky, :, :, :), phi_ext, ulim)
                  ! include the phi_ext contribution to fld_ext
                  fld_ext(:nresponse_per_field) = phi_ext

                  ! if apar evolved, obtain apar on the extended zed domain (apar_ext)
                  ! and include its contribution to fld_ext
                  if (include_apar) then
                     allocate (apar_ext(nresponse_per_field))
                     call map_to_extended_zgrid(it, ie, iky, apar(iky, :, :, :), apar_ext, ulim)
                     fld_ext(offset_apar + 1:nresponse_per_field + offset_apar) = apar_ext
                  end if
                  if (include_bpar) then
                     allocate (bpar_ext(nresponse_per_field))
                     call map_to_extended_zgrid(it, ie, iky, bpar(iky, :, :, :), bpar_ext, ulim)
                     fld_ext(offset_bpar + 1:nresponse_per_field + offset_bpar) = bpar_ext
                  end if

                  ! use LU back substitution to solve linear response matrix system
                  call lu_back_substitution(response_matrix(iky)%eigen(ie)%zloc, &
                                            response_matrix(iky)%eigen(ie)%idx, fld_ext)

                  ! get phi_ext contribution from fld_ext and map back to normal (z, kx) grid
                  phi_ext = fld_ext(:nresponse_per_field)
                  call map_from_extended_zgrid(it, ie, iky, phi_ext, phi(iky, :, :, :))

                  ! if advancing apar, get apar_ext from fld_ext and map back to (z, kx) grid
                  if (include_apar) then
                     apar_ext = fld_ext(offset_apar + 1:nresponse_per_field + offset_apar)
                     call map_from_extended_zgrid(it, ie, iky, apar_ext, apar(iky, :, :, :))
                  end if
                  if (include_bpar) then
                     bpar_ext = fld_ext(offset_bpar + 1:nresponse_per_field + offset_bpar)
                     call map_from_extended_zgrid(it, ie, iky, bpar_ext, bpar(iky, :, :, :))
                  end if

                  deallocate (fld_ext)
                  deallocate (phi_ext)
                  if (allocated(apar_ext)) deallocate (apar_ext)
                  if (allocated(bpar_ext)) deallocate (bpar_ext)
               end do
            end do
         end if
      end do

   end subroutine invert_parstream_response

end module implicit_solve
