! ================================================================================================================================================================================= !
! ---------------------------------- Routines for mapping SFINCS data on to stella grids and calculating all quantites needed for HO GK simulations. ------------------------------ !
! ================================================================================================================================================================================= !
!
! TO DO (OR): DOCUMENTATION TO BE ADDED ...
!
! ================================================================================================================================================================================= !

module neoclassical_terms_sfincs

   ! Load debug flags.
   use debug_flags, only: debug => neoclassical_terms_debug

   ! Neoclassical flag and arrays.
   use arrays_neoclassical, only: include_neoclassical_terms
   use arrays_neoclassical, only: dneo_h_dz, dneo_phi_dz 
   use arrays_neoclassical, only: dneo_h_drho, dneo_phi_drho
   use arrays_neoclassical, only: dneo_h_dvpa, dneo_h_dvpa_global, neo_vpa_fac, neo_vpa_fac_global
   use arrays_neoclassical, only: dneo_h_dmu, dneo_h_dmu_global, neo_mu_fac, neo_mu_fac_global 
   use arrays_neoclassical, only: dneo_h_dz, dneo_phi_dz, neo_zed_fac
   use arrays_neoclassical, only: neo_h, neo_phi
   use arrays_neoclassical, only: dneo_phi_dalpha, neo_alpha_fac
   
   implicit none

   ! Make routines available to other modules.
   public :: init_neoclassical_terms_sfincs
   ! True when SFINCS is the selected drift kinetic provider.
   public :: sfincs_is_enabled

   private

   integer :: nradii
   real :: drho
   integer :: neo_option_switch
   integer, parameter :: neo_option_sfincs = 1

   ! Only initialise once
   logical :: initialised_neoclassical_terms = .false.

contains

! ================================================================================================================================================================================= !
! -------------------------------------------------- Logical function telling stella when to include SFINCS corrections and routines. --------------------------------------------- !
! ================================================================================================================================================================================= !

      logical function sfincs_is_enabled()
         ! Neoclassical. 
         use neoclassical_terms_neo, only: neo_is_enabled

         sfincs_is_enabled = include_neoclassical_terms .and. .not. neo_is_enabled()
      end function sfincs_is_enabled

! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

   subroutine init_neoclassical_terms_sfincs
      ! Parallelisation.
      use parallelisation_layouts, only: vmu_lo

      !MP.
      use mp, only: proc0, broadcast

      ! Grids. 
      use grids_z, only: nzgrid
      use grids_kxky, only: nalpha
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: nspec

      ! For getting SFINCS data. 
      use sfincs_interface, only: get_neo_from_sfincs
      use namelist_neoclassical_input, only: read_namelist_neoclassical_input

      ! Neoclassical diagnostics & moments.
      use neoclassical_diagnostics, only: write_neo_phi_on_stella_z_grid_diagnostic
      use neoclassical_diagnostics, only: write_distribution_moment_diagnostic
      use neoclassical_diagnostics, only: write_distribution_moment_diagnostic
      use neoclassical_terms_neo, only: get_neo_moment, get_neo_h_velocity_derivative_moment

      implicit none

      ! Local variables. 
      real, dimension(:, :, :, :, :, :), allocatable :: f_neoclassical
      real, dimension(:, :, :), allocatable :: phi_neoclassical
      real, dimension(:, :, :, :, :), allocatable :: dfneo_dalpha_local

      ! DIAGNOSTICS. 
      real, dimension(:, :), allocatable :: neo_dens
      real, dimension(:, :), allocatable :: neo_u_par
      real, dimension(:, :), allocatable :: neo_dens_vpa_deriv
      real, dimension(:, :), allocatable :: neo_dens_mu_deriv

      real, dimension(nspec) :: fsa_u_B

      integer :: iz, ialpha

      ! Only initialise once.
      if (initialised_neoclassical_terms) return
      initialised_neoclassical_terms = .true.

      ! Read the input parameters in the input file.
      call read_namelist_neoclassical_input(include_neoclassical_terms, neo_option_switch, nradii, drho)

      ! Broadcast the input parameters to all processors.
      call broacast_arrays

      if (include_neoclassical_terms) then
         ! Allocate all temporary arrays needed for initilization. 
         call allocate_temp_arrays
         
         select case (neo_option_switch)
         case (neo_option_sfincs)
            ! Get H_1 and \Phi^1_0 on all three surfaces. 
            call get_neo_from_sfincs(nradii, drho, f_neoclassical, phi_neoclassical, dfneo_dalpha_local, dneo_phi_dalpha)
         end select

         if (debug) write (6, *) 'neoclassical_terms::init_neoclassical_terms::get_h_and_phi_central_surface'
         call get_h_and_phi_on_central_surface(f_neoclassical(:, :, :, :, :, 0), phi_neoclassical(:, :, 0), &
         neo_h, neo_phi)
         if (debug) write (6, *) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_dzed'
         call get_neo_zed_factor(f_neoclassical(:, :, :, :, :, 0), phi_neoclassical(:, :, 0), &
         dneo_h_dz, dneo_phi_dz, neo_zed_fac)
         if (debug) write (6, *) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_dvpa'
         call get_neo_vpa_factors(f_neoclassical(:, :, :, :, :, 0), phi_neoclassical(:, :, 0), &
            dneo_h_dvpa, dneo_h_dvpa_global, neo_vpa_fac, neo_vpa_fac_global)
         if (debug) write (6, *) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_dmu'
         call get_neo_mu_factors(f_neoclassical(:, :, :, :, :, 0), phi_neoclassical(:, :, 0), &
            dneo_h_dmu, dneo_h_dmu_global, neo_mu_fac, neo_mu_fac_global)
         if (debug) write (6, *) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_drho'
         call get_dfneo_drho(f_neoclassical, dneo_h_drho)
         if (debug) write (6, *) 'neoclassical_terms::init_neoclassical_terms::get_dphineo_drho'
         call get_dphineo_drho(phi_neoclassical, dneo_phi_drho)
         ! TODO (OR): neo_alpha_fac currently holds the raw dF_1/dalpha from SFINCS, not the full alpha factor appearing in the GKE.
         do iz = -nzgrid, nzgrid
            do ialpha = 1, nalpha
               call distribute_vmus_over_procs(dfneo_dalpha_local(ialpha, iz, :, :, :), neo_alpha_fac(ialpha, iz, :))
            end do
         end do
         
         if (debug) write (6, *) 'neoclassical_terms::init_neoclassical_terms::write_neoclassical'
         call write_neoclassical(f_neoclassical, phi_neoclassical)

         ! VERIFICATION. Recompute the out.neo.transport quantities from the SFINCS data on stella grids, convert to NEO units and compare with
         ! NEO on each surface. The full check goes to <run_name>.sfincs_transport_check.dat; a one-line summary is printed. Skipped if the NEO
         ! output files are not in the run directory.
         if (debug) write (6, *) 'neoclassical_terms::init_neoclassical_terms::verify_sfincs_against_neo'
         if (proc0) call verify_sfincs_against_neo(f_neoclassical, phi_neoclassical)

         ! NEO COMPARISON DIAGNOSTICS.
         ! SFINCS PHI ON CENTRAL, LEFT AND RIGHT SURFACES. 
         if (proc0) then
            call write_neo_phi_on_stella_z_grid_diagnostic(phi_neoclassical(1, :, 0), "sfincs_phi_on_stella_z_grid_central_surface")
            call write_neo_phi_on_stella_z_grid_diagnostic(phi_neoclassical(1, :, -1), "sfincs_phi_on_stella_z_grid_left_surface")
            call write_neo_phi_on_stella_z_grid_diagnostic(phi_neoclassical(1, :, 1), "sfincs_phi_on_stella_z_grid_right_surface")
         end if

         ! DIAGNOSTIC. 
         ! SFINCS PHI ON CENTRAL SURFACE ZED DERIVATIVE. 
         if (proc0) then
            call write_neo_phi_on_stella_z_grid_diagnostic(dneo_phi_dz(1, :), "sfincs_dneo_phi_dz_on_stella_z_grid_central_surface")
         end if

         ! DIAGNOSTIC. 
         ! SFINCS PHI ON CENTRAL SURFACE RHO DERIVATIVE.
         if (proc0) then
            call write_neo_phi_on_stella_z_grid_diagnostic(dneo_phi_drho(1, :), "sfincs_dneo_phi_drho_on_stella_z_grid_central_surface")
         end if

         ! DIAGNOSTIC.
         ! Get the zeroth order moment of the neoclassical distribution. 
         call get_neo_moment(f_neoclassical(1, :, :, :, :, 0), &
            phi_neoclassical(1, :, 0), "zeroeth", neo_dens)
         if (proc0) then
            call write_distribution_moment_diagnostic(neo_dens, &
               "sfincs_zeroeth_moment_on_stella_z_grid_central_surface")
         end if

         ! DIAGNOSTIC.
         ! Get the first order moment of the neoclassical distribution. 
         call get_neo_moment(f_neoclassical(1, :, :, :, :, 0), &
            phi_neoclassical(1, :, 0), "first", neo_u_par)
         if (proc0) then
            call write_distribution_moment_diagnostic(neo_u_par, &
               "sfincs_first_moment_on_stella_z_grid_central_surface")
         end if 

         ! DIAGNOSTIC. 
         call get_neo_h_velocity_derivative_moment(dneo_h_dvpa_global(1, :, :, :, :), f_neoclassical(1, :, :, :, :, 0), phi_neoclassical(1, :, 0), "vpa", neo_dens_vpa_deriv)
         call get_neo_h_velocity_derivative_moment(dneo_h_dmu_global(1, :, :, :, :), f_neoclassical(1, :, :, :, :, 0), phi_neoclassical(1, :, 0), "mu", neo_dens_mu_deriv)

         ! DIAGNOSTIC.
         if (proc0) then
            call write_distribution_moment_diagnostic(neo_dens_vpa_deriv, "sfincs_distribution_vpa_deriv_first_moment_on_stella_z_grid_central_surface") 
            call write_distribution_moment_diagnostic(neo_dens_mu_deriv, "sfincs_distribution_mu_deriv_first_moment_on_stella_z_grid_central_surface")
         end if

         ! Deallocate all temporary arrays.
         call deallocate_temp_arrays
      end if

   contains

      subroutine broacast_arrays

         use mp, only: broadcast

         implicit none

         call broadcast(include_neoclassical_terms)
         call broadcast(neo_option_switch)
         call broadcast(nradii)
         call broadcast(drho)

      end subroutine broacast_arrays


      ! Allocates temporary arrays for H_1, ϕ^1_0 and their derivatives. 
      subroutine allocate_temp_arrays
         implicit none

         if (.not. allocated(f_neoclassical)) &
         allocate(f_neoclassical(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec, -nradii / 2:nradii / 2))
         if (.not. allocated(phi_neoclassical)) & 
         allocate(phi_neoclassical(nalpha, -nzgrid:nzgrid, -nradii / 2:nradii / 2))
         if (.not. allocated(dfneo_dalpha_local)) &
         allocate(dfneo_dalpha_local(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec))

         ! DIAGNOSTICS. 
         if (.not. allocated(neo_dens)) allocate(neo_dens(-nzgrid:nzgrid, nspec))
         if (.not. allocated(neo_u_par)) allocate(neo_u_par(-nzgrid:nzgrid, nspec))
         if (.not. allocated(neo_dens_vpa_deriv)) allocate(neo_dens_vpa_deriv(-nzgrid:nzgrid, nspec))
         if (.not. allocated(neo_dens_mu_deriv)) allocate(neo_dens_mu_deriv(-nzgrid:nzgrid, nspec))

      end subroutine allocate_temp_arrays


      ! Deallocate temporary arrays.
      subroutine deallocate_temp_arrays
         implicit none

         if (allocated(f_neoclassical)) deallocate(f_neoclassical)
         if (allocated(phi_neoclassical)) deallocate(phi_neoclassical)
         if (allocated(dfneo_dalpha_local)) deallocate(dfneo_dalpha_local)
    
         ! DIAGNOSTICS. 
         if (allocated(neo_dens)) deallocate(neo_dens)
         if (allocated(neo_u_par)) deallocate(neo_u_par)
         if (allocated(neo_dens_vpa_deriv)) deallocate(neo_dens_vpa_deriv)
         if (allocated(neo_dens_mu_deriv)) deallocate(neo_dens_mu_deriv)

       end subroutine deallocate_temp_arrays

   end subroutine init_neoclassical_terms_sfincs

! ================================================================================================================================================================================= !
! -------------------------------------------- Compare transport quantities from the SFINCS data on stella grids with NEO's out.neo.transport, in NEO units. ---------------------- !
! ================================================================================================================================================================================= !
!
! The same calculations as the NEO verification (neoclassical_terms_neo: verify_neo_transport_on_surface) are applied to the SFINCS H_1 and Φ_1 on stella's grids, using the NEO
! profiles (out.neo.equil, out.neo.species) to convert to NEO units, and compared with out.neo.transport. This needs a NEO run with the same geometry, profiles, collisionality
! and radial spacing drho in the run directory: out.neo.grid, out.neo.species, out.neo.equil(.left/.right) and out.neo.transport(.left/.right).
!
! Normalisation on the left/right surfaces. NEO's H_1 there is normalised to that surface's own Maxwellian (local n_s and v_th,s), whereas the SFINCS H_1 is normalised to the
! central n_s and v_th,s (sfincs_vspace_to_stella_vspace). Both are functions of v / v_th,s on their own surface evaluated on stella's central (v∥, μ) grid. Multiplying the
! SFINCS H_1 by (n_s,c / n_s,l) (T_s,l / T_s,c)^(3/2), with SFINCS's own surface profiles n_s,l = n_s,c (1 - Δρ fprim_s) and T_s,l = T_s,c (1 - Δρ tprim_s)
! (pass_inputoptions_to_sfincs), puts it in NEO's local normalisation, so the NEO conversions (conv_local from that surface's out.neo.equil) apply unchanged.
! Φ_1 is e Φ_1 / T_ref on every surface in both codes, so it needs no rescaling. The central surface needs no rescaling either.
!
! ================================================================================================================================================================================= !

   subroutine verify_sfincs_against_neo(fneo, phineo)
      ! Fortran.
      use iso_fortran_env, only: output_unit

      ! Grids.
      use grids_z, only: nzgrid
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: nspec, spec

      ! Physics parameters and geometry.
      use parameters_physics, only: rhostar
      use geometry, only: geo_surf

      ! NEO data.
      use NEO_interface, only: read_basic_neo_files, read_neo_equil, read_neo_transport
      use NEO_interface, only: neo_grid_data, neo_version_data, neo_equil_data, neo_species_data, neo_transport_data

      ! Transport verification routines (shared with the NEO interface).
      use neoclassical_terms_neo, only: stella_to_neo_units_data, get_stella_to_neo_units, verify_neo_transport_on_surface
      use neoclassical_terms_neo, only: open_neo_check_file, close_neo_check_file_and_summarise, add_neo_check_notice

      ! SFINCS settings.
      use sfincs_interface, only: get_sfincs_radial_electric_field

      implicit none

      real, dimension(:, -nzgrid:, :, :, :, -nradii / 2:), intent(in) :: fneo
      real, dimension(:, -nzgrid:, -nradii / 2:), intent(in)          :: phineo

      ! Local variables.
      type(neo_grid_data)                             :: neo_grid
      type(neo_version_data)                          :: neo_version
      type(neo_species_data)                          :: neo_species
      type(neo_equil_data), dimension(-1:1)           :: neo_equil
      type(neo_transport_data), dimension(-1:1)       :: neo_transport
      type(stella_to_neo_units_data), dimension(-1:1) :: conv
      real, dimension(:, :, :, :), allocatable        :: h_local

      character(len=7), dimension(-1:1), parameter :: surface_label = [character(len=7) :: 'left', 'central', 'right']
      character(len=6), dimension(-1:1), parameter :: suffix = [character(len=6) :: '.left', '', '.right']

      character(len=256) :: message
      logical :: found, er_calculated
      real    :: dphihat_drn, delrho, rho_star_ratio_sfincs, rho_star_ratio_neo, dr_neo
      integer :: irad, irad_lo, irad_hi, is
      integer, parameter :: ia = 1
      real, parameter :: tol = 1.0e-4

      ! Without NEO output in the run directory there is nothing to compare with.
      inquire (file='out.neo.grid', exist=found)
      if (.not. found) then
         write (output_unit, '(a)') ' SFINCS vs NEO transport check skipped: no NEO output (out.neo.grid) in the run directory.'
         return
      end if

      ! NEO only has the central surface and one surface either side; SFINCS may have more (nradii = 5), or only the central one (nradii = 1).
      irad_lo = max(-1, -nradii / 2)
      irad_hi = min(1, nradii / 2)

      ! Read the NEO files (quietly: missing files are reported as notices in the check file).
      call read_basic_neo_files(neo_grid, neo_version, neo_equil(0), neo_species)
      do irad = irad_lo, irad_hi
         if (irad /= 0) call read_neo_equil(neo_equil(irad), neo_grid, suffix = trim(suffix(irad)))
         call read_neo_transport(neo_transport(irad), neo_grid, suffix = trim(suffix(irad)), verbose = .false.)
      end do

      call open_neo_check_file(ext = '.sfincs_transport_check.dat', title = 'SFINCS vs NEO transport check', code_label = 'SFINCS', &
         description = 'Transport quantities from the SFINCS data on stella grids (NEO units) against out.neo.transport (and .left, .right)')

      ! ------------------------------------------------------------------------------------------------------------------------------------ !
      ! Setup checks: things that make SFINCS and NEO differ for physical reasons, or that make the comparison inconsistent.
      ! ------------------------------------------------------------------------------------------------------------------------------------ !

      ! NEO (and stella) assume no equilibrium radial electric field.
      call get_sfincs_radial_electric_field(er_calculated, dphihat_drn)
      if (er_calculated) then
         write (message, '(a, es13.5, a)') 'SFINCS solved for the ambipolar radial electric field (dPhiHat/drN = ', dphihat_drn, &
            '); NEO and stella assume dPhi_0/dr = 0, so flows and fluxes will differ.'
         call add_neo_check_notice(message)
      else if (abs(dphihat_drn) > 0.0) then
         write (message, '(a, es13.5, a)') 'SFINCS was run with dPhiHat/drN = ', dphihat_drn, &
            '; NEO and stella assume dPhi_0/dr = 0, so flows and fluxes will differ.'
         call add_neo_check_notice(message)
      end if

      ! ρ_* on the central surface: SFINCS uses stella's rhostar (Delta), NEO the value in out.neo.equil.
      if (allocated(neo_equil(0)%rho_star)) then
         if (abs(rhostar - neo_equil(0)%rho_star(1)) > tol * abs(neo_equil(0)%rho_star(1))) then
            write (message, '(a, es13.5, a, es13.5)') 'rho_star differs on the central surface: stella/SFINCS = ', rhostar, &
               ', NEO = ', neo_equil(0)%rho_star(1)
            call add_neo_check_notice(message)
         end if
      end if

      do irad = irad_lo, irad_hi
         if (irad == 0) cycle
         if (.not. (allocated(neo_equil(irad)%radius) .and. allocated(neo_equil(0)%radius))) cycle
         delrho = irad * drho

         ! Radial spacing: the left/right surfaces must be the same in both codes.
         dr_neo = neo_equil(irad)%radius(1) - neo_equil(0)%radius(1)
         if (abs(dr_neo - delrho) > 1.0e-3 * abs(delrho)) then
            write (message, '(a, a, a, es13.5, a, es13.5)') 'the ', trim(surface_label(irad)), ' surface is at a different radius: SFINCS delta rho = ', &
               delrho, ', NEO delta r/a = ', dr_neo
            call add_neo_check_notice(message)
         end if

         ! ρ_* off the central surface: SFINCS rescales Delta by (r/q)_local / (r/q)_central (pass_inputoptions_to_sfincs); NEO uses out.neo.equil.
         if (allocated(neo_equil(irad)%rho_star) .and. allocated(neo_equil(0)%rho_star)) then
            rho_star_ratio_sfincs = ((geo_surf%rhoc + delrho) / (geo_surf%qinp * (1.0 + delrho * geo_surf%shat / geo_surf%rhoc))) &
                                    / (geo_surf%rhoc / geo_surf%qinp)
            rho_star_ratio_neo = neo_equil(irad)%rho_star(1) / neo_equil(0)%rho_star(1)
            if (abs(rho_star_ratio_sfincs - rho_star_ratio_neo) > tol * abs(rho_star_ratio_neo)) then
               write (message, '(a, a, a, es13.5, a, es13.5)') 'rho_star(', trim(surface_label(irad)), ') / rho_star(central): SFINCS = ', &
                  rho_star_ratio_sfincs, ', NEO = ', rho_star_ratio_neo
               call add_neo_check_notice(message)
            end if
         end if
      end do

      ! ------------------------------------------------------------------------------------------------------------------------------------ !
      ! stella -> NEO unit ratios on each surface (central ratios, and each surface's local ratios from its out.neo.equil).
      ! ------------------------------------------------------------------------------------------------------------------------------------ !
      do irad = irad_lo, irad_hi
         call get_stella_to_neo_units(conv(irad), neo_equil(irad), neo_species, neo_grid%n_species, trim(surface_label(irad)))
      end do

      ! ------------------------------------------------------------------------------------------------------------------------------------ !
      ! Compare on each surface.
      ! ------------------------------------------------------------------------------------------------------------------------------------ !
      allocate (h_local(-nzgrid:nzgrid, nvpa, nmu, nspec))

      do irad = irad_lo, irad_hi
         delrho = irad * drho
         ! Put the SFINCS H_1 in NEO's local normalisation (factor 1 on the central surface).
         do is = 1, nspec
            h_local(:, :, :, is) = fneo(ia, :, :, :, is, irad) &
                                   * (1.0 - delrho * spec(is)%tprim)**1.5 / (1.0 - delrho * spec(is)%fprim)
         end do

         call verify_neo_transport_on_surface(trim(surface_label(irad)), h_local, phineo(ia, :, irad), neo_transport(irad), &
                                              neo_equil(irad), conv(0), conv(irad))
      end do

      deallocate (h_local)

      call close_neo_check_file_and_summarise

   end subroutine verify_sfincs_against_neo

   !****************************************************************************
   !                Distributed local velocity data over processors. 
   !****************************************************************************
   
   subroutine distribute_vmus_over_procs(local, distributed)

      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx

      implicit none

      real, dimension(:, :, :), intent(in) :: local
      real, dimension(vmu_lo%llim_proc:), intent(out) :: distributed

      integer :: ivmu, iv, imu, is

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         distributed(ivmu) = local(iv, imu, is)
      end do

   end subroutine distribute_vmus_over_procs


   !****************************************************************************
   ! Isolate the non-adiabatic distribution and phi on the central surface;
   ! needed in gk_neo_drive.f90.  
   !****************************************************************************
   subroutine get_h_and_phi_on_central_surface(fneo, phineo, neo_h, neo_phi)
      ! Parallelisation.
      use parallelisation_layouts, only: vmu_lo
     
      ! Grids.
      use grids_z, only: nzgrid
      use grids_kxky, only: nalpha

      implicit none

      real, dimension(:, -nzgrid:, :, :, :), intent(in) :: fneo
      real, dimension(:, -nzgrid:), intent(in) :: phineo
      real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(out) :: neo_h
      real, dimension(:, -nzgrid:), intent(out) :: neo_phi

      ! Local variables.
      integer :: ia, iz, imu, is

      ! Store phi on the central surface in the module level array. 
      neo_phi = phineo

      ! Distribute fneo over processors in the module level array.
      do iz = -nzgrid, nzgrid
         do ia = 1, nalpha
            call distribute_vmus_over_procs(fneo(ia, iz, :, :, :), neo_h(ia, iz, :))
         end do
      end do

   end subroutine get_h_and_phi_on_central_surface


   !****************************************************************************
   ! Get the different parallel velocity derivative arrays that are used 
   ! in the HO corrections.  
   !****************************************************************************
   subroutine get_neo_vpa_factors(fneo, phineo, dneo_h_dvpa, dneo_h_dvpa_global, neo_vpa_fac, neo_vpa_fac_global)
      ! Parallelisation.
      use parallelisation_layouts, only: vmu_lo

      ! Grids. 
      use grids_z, only: nzgrid
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: dvpa, vpa
      use grids_species, only: nspec, spec
      use grids_kxky, only: nalpha

      ! Calculations.
      use calculations_finite_differences, only: fd5pt
         
      implicit none

      real, dimension(:, -nzgrid:, :, :, :), intent(in) :: fneo
      real, dimension(:, -nzgrid:), intent(in) :: phineo
      real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(out) :: dneo_h_dvpa, neo_vpa_fac
      real, dimension(:, -nzgrid:, :, :, :), intent(out) :: dneo_h_dvpa_global, neo_vpa_fac_global

      ! Local variables.
      integer :: ia, iz, imu, is
      real, dimension(:), allocatable :: tmp1, tmp2

      ! Allocate temporary arrays. 
      allocate(tmp1(nvpa), tmp2(nvpa))
      
      ! First calculate neo_vpa_fac_global.
      ! Calculate the vpa derivative of the non-adiabatic distribution first, then subtract the 
      ! maxwellian vpa derivative factor. 
      ! This result is neo_vpa_fac_global on the central surface.
      do is = 1, nspec
         do imu = 1, nmu
            do iz = -nzgrid, nzgrid
               do ia = 1, nalpha
                  ! Hack to avoid dealing with negative indices in fd5pt.
                  tmp1 = fneo(ia, iz, :, imu, is)
                  call fd5pt(tmp1, tmp2, dvpa)
                  dneo_h_dvpa_global(ia, iz, :, imu, is) = tmp2

                  neo_vpa_fac_global(ia, iz, :, imu, is) = dneo_h_dvpa_global(ia, iz, :, imu, is) &
                  - 2.0 * vpa(:) * fneo(ia, iz, :, imu, is) + 2.0 * vpa(:) * spec(is)%zt * phineo(ia, iz)  
               end do
            end do
         end do
      end do

      ! Distribute neo_vpa_fac_global over processors to get neo_vpa_fac.
      do iz = -nzgrid, nzgrid
         do ia = 1, nalpha
            call distribute_vmus_over_procs(neo_vpa_fac_global(ia, iz, :, :, :), neo_vpa_fac(ia, iz, :))
            call distribute_vmus_over_procs(dneo_h_dvpa_global(ia, iz, :, :, :), dneo_h_dvpa(ia, iz, :))
         end do
      end do

      ! Deallocate temporary arrays.
      deallocate (tmp1, tmp2)

   end subroutine get_neo_vpa_factors


   !****************************************************************************
   ! Get the different parallel velocity derivative arrays that are used 
   ! in the HO corrections.  
   !****************************************************************************
   subroutine get_neo_mu_factors(fneo, phineo, dneo_h_dmu, dneo_h_dmu_global, neo_mu_fac, neo_mu_fac_global)
      ! Parallelisation.
      use parallelisation_layouts, only: vmu_lo

      ! Grids. 
      use grids_z, only: nzgrid
      use grids_velocity, only: vpa, nvpa 
      use grids_velocity, only: nmu, dmu
      use grids_species, only: nspec, spec
      use grids_kxky, only: nalpha

      ! Geometry.
      use geometry, only: bmag 

      ! Calculations. 
      use calculations_finite_differences, only: fd3pt
         
      implicit none

      real, dimension(:, -nzgrid:, :, :, :), intent(in) :: fneo
      real, dimension(:, -nzgrid:), intent(in) :: phineo
      real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(out) :: dneo_h_dmu, neo_mu_fac
      real, dimension(:, -nzgrid:, :, :, :), intent(out) :: dneo_h_dmu_global, neo_mu_fac_global

      ! Local variables.
      integer :: ia, iz, iv, is
      real, dimension(:), allocatable :: tmp1, tmp2

      ! Allocate temporary arrays. 
      allocate(tmp1(nmu), tmp2(nmu))
      
      do is = 1, nspec
         do iv = 1, nvpa
            do iz = -nzgrid, nzgrid
               do ia = 1, nalpha
                  ! Hack to avoid dealing with negative indices in fd5pt.
                  tmp1 = fneo(ia, iz, iv, :, is)
                  call fd3pt(tmp1, tmp2, dmu)
                  dneo_h_dmu_global(ia, iz, iv, :, is) = tmp2

                  neo_mu_fac_global(ia, iz, iv, :, is) = dneo_h_dmu_global(ia, iz, iv, :, is) &
                  - 2.0 * bmag(ia, iz) * fneo(ia, iz, iv, :, is) &
                  + 2.0 * bmag(ia, iz) * spec(is)%zt * phineo(ia, iz)  
               end do
            end do
         end do
      end do

      ! Distribute neo_mu_fac_global over processors to get neo_mu_fac. 
      do iz = -nzgrid, nzgrid
         do ia = 1, nalpha
            call distribute_vmus_over_procs(neo_mu_fac_global(ia, iz, :, :, :), neo_mu_fac(ia, iz, :))
            call distribute_vmus_over_procs(dneo_h_dmu_global(ia, iz, :, :, :), dneo_h_dmu(ia, iz, :))
         end do
      end do

      ! Deallocate temporary arrays.
      deallocate (tmp1, tmp2)

   end subroutine get_neo_mu_factors


   !****************************************************************************
   ! Get the different zed derivative arrays that are used 
   ! in the HO corrections.  
   !****************************************************************************
   subroutine get_neo_zed_factor(fneo, phineo, dfneo_dz, dphineo_dz, neo_zed_fac)
      ! Parallelisation. 
      use parallelisation_layouts, only: vmu_lo

      ! Grids.
      use grids_z, only: nztot, nzgrid, delzed
      use grids_velocity, only: vpa, mu, nvpa, nmu
      use grids_species, only: spec, nspec
      use grids_kxky, only: nalpha

      ! Geometry. 
      use geometry, only: bmag, dbdzed 

      ! Calculations. 
      use calculations_finite_differences, only: fd5pt

      implicit none    

      real, dimension(:, -nzgrid:, :, :, :), intent(in)            :: fneo
      real, dimension(:, -nzgrid:), intent(in)                     :: phineo
      real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(out) :: dfneo_dz
      real, dimension(:, -nzgrid:), intent(out)                    :: dphineo_dz
      real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(out) :: neo_zed_fac

      ! Local variables. 
      integer :: iv, imu, is, iz, ia
      real, dimension(:), allocatable :: tmp1, tmp2
      real, dimension(:, :, :, :, :), allocatable :: dfneo_dz_global, neo_zed_fac_global

      ! Allocate temporary arrays. 
      allocate(tmp1(nztot), tmp2(nztot))
      allocate(neo_zed_fac_global(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec))
      allocate(dfneo_dz_global(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec))

      ! First calculate the zed derivative of the non-adiabatic distribution. 
      do is = 1, nspec
         do imu = 1, nmu
            do iv = 1, nvpa
               do ia = 1, nalpha
                  ! Hack to avoid dealing with negative indices in fd5pt.
                  tmp1 = fneo(ia, :, iv, imu, is)
                  call fd5pt(tmp1, tmp2, delzed(0))
                  ! Store the result in neo_zed_fac_global.
                  dfneo_dz_global(ia, :, iv, imu, is) = tmp2
               end do
            end do
         end do
      end do

      ! Now calculate he zed derivative of the neoclassical potential. 
      ! Hack to avoid dealing with negative indices in fd5pt.
      do ia = 1, nalpha
        tmp1 = phineo(ia, :)
         call fd5pt(tmp1, tmp2, delzed(0))
         dphineo_dz(ia, :) = tmp2
      end do 

      ! Now get the zed factor appearing in the GKE. 
      do is = 1, nspec
         do iv = 1, nvpa
            do imu = 1, nmu
               do iz = -nzgrid, nzgrid
                  do ia = 1, nalpha 
                     neo_zed_fac_global(ia, iz, iv, imu, is) =  &
                     dfneo_dz_global(ia, iz, iv, imu, is) &
                     - spec(is)%zt * dphineo_dz(ia, iz) - 2.0 * mu(imu) * dbdzed(ia, iz) &
                     * ( fneo(ia, iz, iv, imu, is) - spec(is)%zt * phineo(ia, iz) )
                  end do
               end do
            end do
         end do
      end do

      ! Distribute end result over processors. 
      do iz = -nzgrid, nzgrid
         do ia = 1, nalpha
            call distribute_vmus_over_procs(dfneo_dz_global(ia, iz, :, :, :), dfneo_dz(ia, iz, :))
            call distribute_vmus_over_procs(neo_zed_fac_global(ia, iz, :, :, :), neo_zed_fac(ia, iz, :))
         end do
      end do

      ! Deallocate temporary arrays.
      deallocate(dfneo_dz_global)
      deallocate(neo_zed_fac_global)
      deallocate(tmp1, tmp2)

    end subroutine get_neo_zed_factor


   !****************************************************************************
   !                                      Title
   !****************************************************************************
   !****************************************************************************
   ! Radial derivative of the SFINCS H_1 on the central surface.
   !****************************************************************************
   subroutine get_dfneo_drho(fneo, dfneo)

      use calculations_finite_differences, only: fd3pt, fd5pt
      use grids_z, only: nzgrid
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: nspec, spec                                   
      use parallelisation_layouts, only: vmu_lo
      use grids_kxky, only: nalpha

      implicit none

      real, dimension(:, -nzgrid:, :, :, :, -nradii/2:), intent(in) :: fneo
      real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(out) :: dfneo

      integer :: ia, iz, iv, imu, is, irad                                   
      real :: delrho                                                         
      real, dimension(:), allocatable :: tmp1, tmp2
      real, dimension(:, :), allocatable :: to_local_normalisation           
      real, dimension(:, :, :, :, :), allocatable :: dfneo_local

      ! Allocate local arrays. 
      allocate(tmp1(nradii), tmp2(nradii))
      allocate(to_local_normalisation(nspec, nradii))                       
      allocate(dfneo_local(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec))

      ! Factor converting the SFINCS H_1 on each surface 
      ! to that surface's own Maxwellian normalisation, as is done in NEO.
      do irad = -nradii / 2, nradii / 2
         delrho = irad * drho
         do is = 1, nspec
            to_local_normalisation(is, irad + nradii / 2 + 1) = (1.0 - delrho * spec(is)%tprim)**1.5 &
            / (1.0 - delrho * spec(is)%fprim)
         end do
      end do

      do is = 1, nspec
         do imu = 1, nmu
            do iv = 1, nvpa
               do iz = -nzgrid, nzgrid
                  do ia = 1, nalpha
                     ! hack to avoid dealing with negative indices in fd5pt
                     tmp1 = fneo(ia, iz, iv, imu, is, :) * to_local_normalisation(is, :)  
                     if (nradii == 5) then
                        call fd5pt(tmp1, tmp2, drho)
                     else
                        call fd3pt(tmp1, tmp2, drho)
                     end if
                     dfneo_local(ia, iz, iv, imu, is) = tmp2(nradii / 2 + 1)
                  end do
               end do
            end do
         end do
      end do

      do iz = -nzgrid, nzgrid
         do ia = 1, nalpha
            call distribute_vmus_over_procs(dfneo_local(ia, iz, :, :, :), dfneo(ia, iz, :))
         end do
      end do

      ! Deallocate local arrays.
      deallocate(dfneo_local)
      deallocate(to_local_normalisation)                                    
      deallocate(tmp1, tmp2)

   end subroutine get_dfneo_drho


   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine get_dphineo_drho(phineo, dphineo)

      use calculations_finite_differences, only: fd3pt, fd5pt
      use grids_z, only: nzgrid
      use grids_kxky, only: nalpha

      implicit none

      real, dimension(:, -nzgrid:, -nradii/2:), intent(in) :: phineo
      real, dimension(:, -nzgrid:), intent(out) :: dphineo

      integer :: iz, ia
      real, dimension(:), allocatable :: tmp1, tmp2

      allocate (tmp1(nradii), tmp2(nradii))

      do iz = -nzgrid, nzgrid
         do ia = 1, nalpha
            ! hack to avoid dealing with negative indices in fd5pt
            tmp1 = phineo(ia, iz, :)
            if (nradii == 5) then
               call fd5pt(tmp1, tmp2, drho)
            else
               call fd3pt(tmp1, tmp2, drho)
            end if
            dphineo(ia, iz) = tmp2(nradii / 2 + 1)
         end do
      end do

      deallocate (tmp1, tmp2)

   end subroutine get_dphineo_drho

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine write_neoclassical(fnc, phinc)

      use mp, only: proc0
      use mp, only: send, receive
      use file_utils, only: open_output_file, close_output_file
      use grids_z, only: nzgrid, zed
      use grids_velocity, only: vpa, mu
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
      use parallelisation_layouts, only: idx_local, proc_id
      use grids_kxky, only: nalpha

      implicit none

      real, dimension(:, -nzgrid:, :, :, :, -nradii/2:), intent(in) :: fnc
      real, dimension(:, -nzgrid:, -nradii/2:), intent(in) :: phinc

      integer :: neo_unit
      integer :: irad, iz, ivmu, iv, imu, is, ia
      real, dimension(:, :), allocatable :: dfdv_local, dfdr_local, dfdz_local

      allocate (dfdv_local(nalpha, -nzgrid:nzgrid))
      allocate (dfdr_local(nalpha, -nzgrid:nzgrid))
      allocate (dfdz_local(nalpha, -nzgrid:nzgrid))

      if (proc0) then
         call open_output_file(neo_unit, '.neoclassical')
         write (neo_unit, '(3a8,10a13)') '#1.rad', '2.spec', '3.alpha', '4.zed', '5.mu', &
            '6.vpa', '7.f_neo', '8.dfdvpa_neo', '9.dfdrho_neo', '10.dfdzed_neo', '11.phi_neo', &
            '12.dphidrho', '13.dphidzed'
      end if
      do irad = -nradii / 2, nradii / 2
         do ivmu = vmu_lo%llim_world, vmu_lo%ulim_world
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            if (idx_local(vmu_lo, iv, imu, is)) then
               if (proc0) then
                  dfdv_local = dneo_h_dvpa(:, :, ivmu)
                  dfdr_local = dneo_h_drho(:, :, ivmu)
                  dfdz_local = dneo_h_dz(:, :, ivmu)
               else
                  call send(dneo_h_dvpa(:, :, ivmu), 0)
                  call send(dneo_h_drho(:, :, ivmu), 0)
                  call send(dneo_h_dz(:, :, ivmu), 0)
               end if
            else if (proc0) then
               call receive(dfdv_local, proc_id(vmu_lo, ivmu))
               call receive(dfdr_local, proc_id(vmu_lo, ivmu))
               call receive(dfdz_local, proc_id(vmu_lo, ivmu))
            end if
            if (proc0) then
               do iz = -nzgrid, nzgrid
                  do ia = 1, nalpha
                     write (neo_unit, '(3i8,10e13.5)') irad, is, ia, zed(iz), mu(imu), vpa(iv), &
                        fnc(ia, iz, iv, imu, is, irad), &
                        dfdv_local(ia, iz), &
                        dfdr_local(ia, iz), &
                        dfdz_local(ia, iz), &
                        phinc(ia, iz, irad), dneo_phi_drho(ia, iz), dneo_phi_dz(ia, iz)
                  end do
               end do
            end if
         end do
         if (proc0) write (neo_unit, *)
      end do
      if (proc0) call close_output_file(neo_unit)

      deallocate (dfdv_local, dfdr_local, dfdz_local)

   end subroutine write_neoclassical


end module neoclassical_terms_sfincs