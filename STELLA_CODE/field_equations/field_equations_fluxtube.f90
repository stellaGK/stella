!###############################################################################
!################### ADVANCE FIELDS IN FLUXTUBE SIMULATION #####################
!###############################################################################
! 
! Evolve the fields in time using the field equations in a fluxtube.
! In this module, the quasi-neutrality condition is solved, and if electromagnetic
! effects or radial variation are included, then the appropriate routines are called.
!
! For electrostatic fluxtube simulations this is: 
!     In terms of the distribution function g: 
!     sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g_s + (Z_s n_s/T_s) * (Gamma0 - 1) * phi ] = 0

!     If h is being used as the distribution function, then QN is: 
!     sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * h_s - (Z_s n_s/T_s) * F_0 * phi ] = 0
! (the elelctromagnetic version of QN is in field_equations_electromagnetic.fpp)
! 
! Here we used the guiding-center disitribution function g_s and the non-adiabatic part h_s
!     g_s = <delta f_s>_R = h_s - Z_s n_s/T_s * <phi>_R * F_0
! 
! The arguments of the Bessel functions are
!     J_0(a_k) = J_0(k_perp * rho_s)
!     Gamma0(b_k) = Gamma0(k_perp**2 * rho_s**2 / 2)
! 
! This equation can be rewritten in order to obtain the electrostatic potential:
!     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] / [ sum_s (Z_s² n_s/T_s) (1 - Gamma0) ]
!     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * h ] / [ sum_s (Z_s² n_s/T_s) ]
! 
! The denominators are constants and are calculated when initialising stella
!     denominator_fields[iky,ikz,iz] = sum_s (Z_s² n_s/T_s) (1 - Gamma0)
!     denominator_fields_h = sum_s (Z_s² n_s/T_s)
! 
! The integral over velocity space and species is calculated in stella as
!     integrate_species( . ) = sum_s (2B/sqrt(pi)) int dvpa int dmu ( . )
! 
! To summarize, the fields (= electrostatic potential) can be calculated as
!     phi = integrate_species( J_0 * g_s ) / denominator_fields
!     phi = integrate_species( J_0 * h_s ) / denominator_fields_h
! 
!###############################################################################
module field_equations_fluxtube

   ! Load debug flags
   use debug_flags, only: debug => fields_fluxtube_debug
   
   implicit none

   ! Make routines available to other modules
   public :: advance_fields_fluxtube
   public :: advance_fields_fluxtube_using_field_equations
   public :: init_field_equations_fluxtube

   private

   ! Advance fields for g(kx,ky,z,ivpamus) and g(vpa,mu,ikxkyzs), depending on the layout
   interface advance_fields_fluxtube_using_field_equations
      module procedure advance_fields_fluxtube_using_field_equations_vmulo
      module procedure advance_fields_fluxtube_using_field_equations_kxkyzlo
   end interface

contains

!###############################################################################
!################### ADVANCE FIELDS IN FLUXTUBE SIMULATION #####################
!###############################################################################

   !****************************************************************************
   !                    ADVANCE FIELDS IN FLUXTUBE SIMULATION              
   !****************************************************************************
   ! The fields (electrostatic potential and electromagnetic fields) are evolved
   ! in time though the field equations when simulating a fluxtube domain.
   !
   ! Electrostatic quasi-neutrality equation is the 
   ! minimum field equagtion that must be solved. This can be extended depending on 
   ! the physics that is being included. 
   ! 
   ! Note that Apar and Bpar are only advanced when using electromagnetic stella, 
   ! so these are in field_equations_electromagnetic.fpp
   !****************************************************************************
   subroutine advance_fields_fluxtube(g, phi, apar, bpar, dist, skip_fsa)

      ! Parallelisation
      use mp, only: proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo
      use redistribute, only: scatter
      use initialise_redistribute, only: kxkyz2vmu
      use parallelisation_layouts, only: fields_kxkyz
      use timers, only: time_field_solve

      ! Arrays
      use arrays_distribution_function, only: gvmu

      ! Grids
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(inout) :: phi, apar, bpar
      character(*), intent(in) :: dist
      logical, optional, intent(in) :: skip_fsa
      
      !-------------------------------------------------------------------------
      
      ! Note that fields_kxkyz = F is the default, so this is the default routine used
      if (.not. fields_kxkyz) then

         ! This will call advance_fields_fluxtube_using_field_equations_vmulo
         if (debug) write (*, *) 'field_equations_fluxtube::advance_fields_fluxtube::vmulo'
         call advance_fields_fluxtube_using_field_equations(g, phi, apar, bpar, dist, skip_fsa)

      ! This is the less common option where fields_kxkyz = T  
      ! Note that to use this option it has to be specified by the user
      else if (fields_kxkyz) then 
      
         ! First gather (vpa,mu) onto processor for v-space operations
         ! and parallelise over (kx,ky,z). This changes g -> gvmu
         if (debug) write (*, *) 'field_equations_fluxtube::advance_fields_fluxtube::scatter'
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' field_equations_redist')
         call scatter(kxkyz2vmu, g, gvmu)
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' field_equations_redist')
         
         ! Given gvmu with vpa and mu local, calculate the corresponding fields
         ! This will call advance_fields_fluxtube_using_field_equations_kxkyzlo
         if (debug) write (*, *) 'field_equations_fluxtube::advance_fields_fluxtube::kxkyzlo'
         call advance_fields_fluxtube_using_field_equations(gvmu, phi, apar, bpar, dist, skip_fsa)
         
      end if

   end subroutine advance_fields_fluxtube

   !****************************************************************************
   !  ADVANCE FIELDS USING THE QUASINEUTRALITY EQUATION AND G(KX,KY,Z,IVPAMUS)  
   !****************************************************************************
   ! If <g> is parallelised over (vpa,mu,s) then this subroutine is called.
   ! This is the more common version used compared with parallelising over
   ! (kx,ky,z) and is the default for stella.
   ! 
   ! Here we calculate the electrostatic potential based on the quasi-neutrality equation:
   !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] / [ sum_s (Z_s²ns/Ts) (1 - Gamma0) ]
   !     phi = integrate_species( J_0 * g ) / denominator_fields
   ! 
   ! TODO-GA: remove apar from this and make it only needed for EM stella
   !****************************************************************************
   subroutine advance_fields_fluxtube_using_field_equations_vmulo(g, phi, apar, bpar, dist, skip_fsa)

      ! Parallelisation
      use mp, only: mp_abort, proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx
      use timers, only: time_field_solve
      
      ! Arrays
      use arrays_distribution_function, only: g_scratch
      
      ! Parameters
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: radial_variation
      use parameters_physics, only: fphi
      
      ! Grids
      use grids_z, only: nzgrid
      use grids_species, only: spec
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      use calculations_velocity_integrals, only: integrate_species
      
      ! Routines from other field modules
      use field_equations_electromagnetic, only: advance_fields_electromagnetic
      use field_equations_radialvariation, only: add_radial_correction_int_species
      use field_equations_radialvariation, only: calculate_phi_for_radial_variation
      
      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar, bpar
      character(*), intent(in) :: dist 
      logical, optional, intent(in) :: skip_fsa
      
      ! Local variables
      logical :: skip_fsa_local
      
      !-------------------------------------------------------------------------
      
      ! Debug message
      if (debug) write (*, *) 'field_equations_fluxtube::vmulo'
      
      ! Used for the Dougherty collision operator
      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ! Initialise the electrostatic potential phi
      phi = 0.
      
      ! Note that this advances phi for the electrostatic, fluxtube case.
      ! If electromagnetic effects are included then phi will be advanced below
      if (fphi > epsilon(0.0) .and. .not. include_bpar) then
      
         ! Debug message
         if (debug) write (*, *) 'field_equations_fluxtube::vmulo::electrostatic'
         ! Start timer
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         
         ! First gyro-average the distribution function g at each phase space location
         ! and store this as g_scratch = <g>_R = J_0 g in k-space
         call gyro_average(g, g_scratch)
         
         ! If we are allowing for radial variation then we must modify <g>.
         ! This is done in the field_equations_radialvariation module, but 
         ! is not needed for standard stella (as these are false by default).
         if (radial_variation) call add_radial_correction_int_species(g_scratch)
         
         ! Calculate phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] = integrate_species( J_0 * g )
         if (debug) write (*, *) 'field_equations_fluxtube::vmulo::integrate_species_phi'
         call integrate_species(g_scratch, spec%z * spec%dens_psi0, phi)

         ! Stop timer
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         
         ! Calculate phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] / [ sum_s (Z_s² n_s/T_s) (1 - Gamma0) ]
         ! by dividing with denominator_fields[iky,ikz,iz] = sum_s (Z_s² n_s/T_s) (1 - Gamma0) in the calculate_phi() routine
         if (.not. radial_variation) then
            if (debug) write (*, *) 'field_equations_fluxtube::vmulo::calculate_phi'
            call calculate_phi(phi, dist, skip_fsa_local)
            
         ! The denominator_fields[iky,ikz,iz] factor is modified by radial varation effects
         else if (radial_variation) then
            if (debug) write (*, *) 'field_equations_fluxtube::vmulo::calculate_phi_for_radial_variation'
            call calculate_phi_for_radial_variation(phi, dist , skip_fsa)
            
         end if
      end if

      ! If we have electromagnetic effects we need to calculate the fields A_parallel and B_parallel
      if (include_apar .or. include_bpar) then
         if (debug) write (*, *) 'field_equations_fluxtube::vmulo::electromagnetic'
         apar = 0.; bpar = 0.
         call advance_fields_electromagnetic(g, phi, apar, bpar, dist)
      end if

   end subroutine advance_fields_fluxtube_using_field_equations_vmulo

   !****************************************************************************
   !  ADVANCE FIELDS USING THE QUASINEUTRALITY EQUATION AND G(MU,VPA,IKXKYZS)  
   !****************************************************************************
   ! If <g> is parallelised over (kx,ky,z,s) then this subroutine is called.
   ! 
   ! Here we calculate the electrostatic potential based on the quasi-neutrality equation:
   !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] / [ sum_s (Z_s²n_s/T_s) (1 - Gamma0) ]
   !     phi = integrate_species( J_0 * g ) / denominator_fields
   ! 
   ! TODO-GA: remove apar from this and make it only needed for EM stella
   !****************************************************************************
   subroutine advance_fields_fluxtube_using_field_equations_kxkyzlo(g, phi, apar, bpar, dist, skip_fsa)

      ! Parallelisation
      use mp, only: proc0
      use mp, only: sum_allreduce, mp_abort
      use job_manage, only: time_message
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use timers, only: time_field_solve
      
      ! Parameters
      use parameters_physics, only: fphi
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: radial_variation
      
      ! Grids
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: spec
      use grids_z, only: nzgrid
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      use calculations_velocity_integrals, only: integrate_vmu
      
      ! Routines from other field modules
      use field_equations_electromagnetic, only: advance_fields_electromagnetic
      use field_equations_radialvariation, only: calculate_phi_for_radial_variation

      implicit none

      ! Arguments
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(inout) :: phi, apar, bpar
      logical, optional, intent(in) :: skip_fsa
      character(*), intent(in) :: dist
      
      ! Local variables
      complex :: tmp
      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      logical :: skip_fsa_local
      
      !-------------------------------------------------------------------------
      
      ! Debug message
      if (debug) write (*, *) 'field_equations_fluxtube::kxkyzlo'
      
      ! Used for the Dougherty collision operator
      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ! Assume we only have one field line
      ia = 1
      
      ! Initialise the electrostatic potential phi
      phi = 0.

      ! Note that this advances phi for the electrostatic, fluxtube case.
      ! If electromagnetic effects are included then phi will be advanced below
      if (fphi > epsilon(0.0) .and. .not. include_bpar) then

         ! Debug message
         if (debug) write (*, *) 'field_equations_fluxtube::kxkyzlo::electrostatic'
         ! Start timer
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         
         ! Allocate temporary arrays
         allocate (g0(nvpa, nmu))
         
         ! Iterate over the (kx,ky,z,mu,vpa,s) points
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            
            ! First gyro-average the distribution function g at each phase space location
            ! and store this as g0 = <g>_R = J_0 g
            call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
            
            ! Calculate phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ]
            wgt = spec(is)%z * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp
            
         end do
         
         ! Deallocate temporary array
         deallocate (g0)
         
         ! Sum the values on all processors and send them to <proc0>
         call sum_allreduce(phi)

         ! Calculate phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] / [ sum_s (Z_s² n_s/T_s) (1 - Gamma0) ]
         ! by dividing with denominator_fields[iky,ikz,iz] = sum_s (Z_s² n_s/T_s) (1 - Gamma0) in the calculate_phi() routine
         if (.not. radial_variation) then 
            if (debug) write (*, *) 'field_equations_fluxtube::kxkyzlo::calculate_phi'
            call calculate_phi(phi, dist, skip_fsa_local)

         ! The denominator_fields[iky,ikz,iz] factor is modified by radial varation effects
         else if (radial_variation) then
            if (debug) write (*, *) 'field_equations_fluxtube::kxkyzlo::calculate_phi_for_radial_variation'
            call calculate_phi_for_radial_variation (phi, dist , skip_fsa)
         end if
         
      end if

      ! If we have electromagnetic effects we need to calculate the fields A_parallel and B_parallel
      if (include_apar .or. include_bpar) then 
         if (debug) write (*, *) 'field_equations_fluxtube::kxkyzlo::electromagnetic'
         bpar = 0.; apar = 0.
         call advance_fields_electromagnetic(g, phi, apar, bpar, dist)
      end if

   end subroutine advance_fields_fluxtube_using_field_equations_kxkyzlo

   !****************************************************************************
   !******************************* CALCULATE PHI ******************************
   !****************************************************************************
   ! Here we calculate the electrostatic potential based on the quasi-neutrality equation:
   !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] / [ sum_s (Z_s² n_s/T_s) (1 - Gamma0) ]
   !     phi = integrate_species( J_0 * g ) / denominator_fields
   ! 
   ! Note that the 'phi' variable passed in is:
   !    integrate_species( J_0 * g ) = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ]
   ! 
   ! This routine divides by the appropriate <denominator_fields> factor depending on if we
   ! have kinetic or adiabatic electrons, and also on whether we are using 'g'
   ! or 'h' as our distribution function that we are evolving.
   !****************************************************************************
   subroutine calculate_phi(phi, dist, skip_fsa)

      ! Parallelisation
      use mp, only: proc0, mp_abort
      use job_manage, only: time_message
      use multibox, only: mb_calculate_phi
      use timers, only: time_field_solve
      
      ! Arrays
      use arrays, only: denominator_fields
      use arrays, only: denominator_fields_MBR
      use arrays, only: denominator_fields_h
      use arrays, only: denominator_fields_MBR_h
      
      ! Grids
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: zonal_mode
      use grids_kxky, only: nakx, naky
      use grids_species, only: spec
      
      ! Adiabatic electrons
      use grids_species, only: has_electron_species
      use grids_species, only: adiabatic_option_switch
      use grids_species, only: adiabatic_option_fieldlineavg
      
      ! Geometry
      use geometry, only: dl_over_b

      implicit none

      ! Arguments
      character(*), intent(in) :: dist
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
      logical, optional, intent(in) :: skip_fsa
      
      ! Local variables
      real, dimension(:, :, :, :), allocatable :: denominator_fields_t
      integer :: ia, it, ikx
      complex :: tmp
      logical :: skip_fsa_local
      logical :: has_elec, adia_elec

      !-------------------------------------------------------------------------
      
      ! Debug message
      if (debug) write (*, *) 'field_equations_fluxtube::calculate_phi'
      
      ! Used for the Dougherty collision operator
      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ! Assume we only have one field line
      ia = 1
      
      ! Check if we need to add a Boltzmann response for an adiabatic species
      has_elec = has_electron_species(spec)
      adia_elec = .not. has_elec .and. (adiabatic_option_switch == adiabatic_option_fieldlineavg)

      ! Start timer
      if (proc0) call time_message(.false., time_field_solve(:, 4), ' calculate_phi')
      
      ! ------------------------------------------------------------------------------------------
      !                         Using h as the distribution function
      ! ------------------------------------------------------------------------------------------
      ! If we are using the non-adiabatic part h of the distribution function then
      ! sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * h - (Z_s n_s/T_s) phi ] = 0
      ! phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * h ] / [ sum_s (Z_s² n_s/T_s) ]
      ! denominator_fields_h = sum_s (Z_s² n_s/T_s)
      if (dist == 'h') then
         if (debug) write(*, *) 'field_equations_quasineutrality::fluxtube::calculate_phi::dist==h'
         phi = phi / denominator_fields_h
         
      ! ------------------------------------------------------------------------------------------
      !                    Using g or gbar as the distribution function
      ! ------------------------------------------------------------------------------------------
      ! If we are using the guiding-center distribution function g then
      ! sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g + (Z_s n_s/T_s) (Gamma0 - 1) phi ] = 0
      ! phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] / [ sum_s (Z_s² n_s/T_s) (1 - Gamma0) ]
      ! denominator_fields[iky,ikz,iz] = sum_s (Z_s² n_s/T_s) (1 - Gamma0)
      !
      ! To avoid any issues with division by zero we set phi = 0.0 if the denominator is too small. Only thing
      ! that makes sense is to set phi = 0.0 if the prefactor for phi in QN is also zero. 
      else if (dist == 'g' .or. dist == 'gbar') then
         if (debug) write(*, *) 'field_equations_quasineutrality::fluxtube::calculate_phi::dist==gbar'
         allocate (denominator_fields_t(naky, nakx, -nzgrid:nzgrid, ntubes))
         denominator_fields_t = spread(denominator_fields, 4, ntubes)
         where (denominator_fields_t < epsilon(0.0))
            phi = 0.0
         elsewhere
            phi = phi / denominator_fields_t
         end where
         deallocate (denominator_fields_t)
         
      ! ------------------------------------------------------------------------------------------
      !                             Other - for safety just abort
      ! ------------------------------------------------------------------------------------------
      ! Abort if <dist> is not recognized.
      else
         if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
         call mp_abort('unknown dist option in get_fields. aborting')
         return
      end if

      ! ------------------------------------------------------------------------------------------
      !                              Set kx = ky = 0.0 mode to zero
      ! ------------------------------------------------------------------------------------------
      ! The kx = ky = 0.0 mode is not evolved by stella so make sure this term is set to zero.
      if (debug) write(*, *) 'field_equations_quasineutrality::fluxtube::calculate_phi::set kxky=0.0 to zero'
      if (any(denominator_fields(1, 1, :) < epsilon(0.))) phi(1, 1, :, :) = 0.0
      if (proc0) call time_message(.false., time_field_solve(:, 4), ' calculate_phi')

      ! ------------------------------------------------------------------------------------------
      !                                   Adiabatic electrons
      ! ------------------------------------------------------------------------------------------
      ! Handle adiabatic electrons only if needed. If using kinetic electrons, then stella will not
      ! use these routines.
      ! If using adiabatic electrons then we need to modify phi. We do not explicitly have the
      ! electron distribution function, and use g_e = e <phi>_R * F_0e / T_e + const.
      ! where the const. is chosen to be either: 
      !        - zero if using adiabatic
      !                 sum_(s not e) Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] 
      !                                                        = n_e e^2 / T_e * phi 
      !        - such that the flux-surface average of g_e is zero if using Modified Boltzmann Response (MBR).
      !          In order for this to be satisfied we find 
      !                          \delta n_e = e n_e (phi - <phi>_R) / T_e
      !          where <phi>_FSA is the flux surface average of phi.
      !                 sum_(s not e) Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] 
      !                                                        = n_e e^2 / T_e (phi - <phi>_FSA)
      !
      ! More information about how these equations are solved, and how the denominators are modified
      ! can be found below in the 'initialise' subroutine - scroll down to Adiabatic electrons.
      ! ------------------------------------------------------------------------------------------
      ! For g: 
      !     denominator_fields = sum_s (Z_s² n_s/T_s) (1 - Gamma0)
      !     denominator_fields_MBR = 1 / (T_e/n_e - <1/denominator_field>_FSA ) - 1
      ! For h:
      !     denominator_fields_h = sum_s (Z_s² n_s/T_s) + (n_e/T_e)
      !     denominator_fields_MBR_h = (T_e/n_e - <1/ denominator_fields_h>_FSA ) - 1
      !  => denominator_fields_MBR_h = n_e/T_e / sum_(s not e) (Z_s² n_s/T_s)
      !
      ! ------------------------------------------------------------------------------------------

      ! Start timer 
      if (proc0) call time_message(.false., time_field_solve(:, 5), 'calculate_phi_adia_elec')

      if (adia_elec .and. zonal_mode(1) .and. .not. skip_fsa_local) then
         ! ---------------------------------------------------------------------------------------
         !                 Adiabatic electrons - Using h as the distribution function
         ! ---------------------------------------------------------------------------------------
         if (dist == 'h') then
            ! Debug message
            if (debug) write(*, *) 'field_equations_quasineutrality::fluxtube::calculate_phi::adiabatic_electrons::dist==h'

            ! Add the Boltzmann response for adiabatic electrons
            do it = 1, ntubes
               do ikx = 1, nakx
                  tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                  phi(1, ikx, :, it) = phi(1, ikx, :, it) + tmp * denominator_fields_MBR_h
               end do
            end do

         ! ---------------------------------------------------------------------------------------
         !              Adiabatic electrons - Using g or gbar as the distribution function
         ! ---------------------------------------------------------------------------------------
         else if (dist == 'g' .or. dist == 'gbar') then
            if (debug) write(*, *) 'field_equations_quasineutrality::fluxtube::calculate_phi::adiabatic_electrons::dist==gbar'
            do ikx = 1, nakx
               do it = 1, ntubes
                  tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                  phi(1, ikx, :, it) = phi(1, ikx, :, it) + tmp * denominator_fields_MBR(ikx, :)
               end do
            end do

         ! ---------------------------------------------------------------------------------------
         !                      Adiabatic electrons - Abort if unknown dist fn
         ! ---------------------------------------------------------------------------------------
         else
            if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
            call mp_abort('unknown dist option in get_fields. aborting')
         end if
      end if
      
      ! Stop timer
      if (proc0) call time_message(.false., time_field_solve(:, 5), 'calculate_phi_adia_elec')

   end subroutine calculate_phi

!###############################################################################
!################################# INITIALISE ##################################
!###############################################################################

   !****************************************************************************
   !************************** INITIALISE DENOMINATORS *************************
   !****************************************************************************
   ! The electrostatic potential phi is calculated based on the QN condition
   !     sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g + (Z_s n_s/T_s) (Gamma0 - 1) phi ] = 0
   !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] / [ sum_s (Z_s² n_s/T_s) (1 - Gamma0) ]
   !
   !     denominator_fields[iky,ikz,iz] = sum_s (Z_s² n_s/T_s) (1 - Gamma0)
   !
   ! The denominators needed to calculate <phi> are initialised in this routine.
   ! 
   ! If adiabatic electrons are used then this factor is modified and we use 
   ! <denominator_fields_MBR> which includes the Modified Boltmann Response.
   ! 
   ! Use the following calculations:
   !      (1 - Gamma0(b_k) = (2B/sqrt(pi)) int dvpa int dmu (1 - J_0(a_k)²) exp(v²)
   !      integrate_vmu( . ) = (2B/sqrt(pi)) int dvpa int dmu ( . )
   !      b_k = k²_perp*rho²_s/2
   !      a_k = k_perp*rho_s
   !
   ! ---------------------------------------------------------------------------
   !                                   DEFINITIONS
   ! ---------------------------------------------------------------------------
   ! If kinetic electrons are used : 
   ! For g:
   !     denominator_fields = sum_s (Z_s² n_s/T_s) (1 - Gamma0)
   ! For h: 
   !     denominator_fields_h = sum_s (Z_s² n_s/T_s)
   !
   ! If adiabatic electrons are used : 
   ! For g:
   !     denominator_fields = sum_(s not e) (Z_s² n_s/T_s) (1 - Gamma0) + (n_e/T_e)
   ! For h:
   !     denominator_fields_h = sum_s (Z_s² n_s/T_s) + (n_e/T_e)
   ! 
   ! If a Modofied Boltzmann Response (MBR) is used :
   ! For g: 
   !     denominator_fields_MBR = 1 / (T_e/n_e - <1/denominator_field>_FSA ) - 1
   ! For h:
   !     denominator_fields_MBR_h = (T_e/n_e - <1/ denominator_fields_h>_FSA ) - 1
   !  => denominator_fields_MBR_h = n_e/T_e / sum_(s not e) (Z_s² n_s/T_s)
   ! 
   !****************************************************************************
   subroutine init_field_equations_fluxtube

      ! Parallelisation
      use mp, only: sum_allreduce
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, onlY: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      
      ! Arrays
      use arrays, only: denominator_fields, denominator_fields_MBR
      use arrays, only: denominator_fields_h, denominator_fields_MBR_h, efac, efacp
      use arrays_gyro_averages, only: aj0v

      ! Parameters
      use parameters_numerical, only: maxwellian_normalization
      use parameters_physics, only: fphi
      
      ! Grids
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: spec
      use grids_kxky, only: zonal_mode, akx
      use grids_kxky, only: nakx
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac

      ! Adiabatic electrons
      use grids_species, only: has_electron_species
      use grids_species, only: ion_species
      use grids_species, only: tite, nine
      use grids_species, only: adiabatic_option_switch
      use grids_species, only: adiabatic_option_fieldlineavg

      ! Calculations
      use calculations_velocity_integrals, only: integrate_vmu
      
      ! Geometry
      use geometry, only: dl_over_b

      implicit none

      ! Local variables
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      real :: tmp, wgt
      real, dimension(:, :), allocatable :: g0
      
      !-------------------------------------------------------------------------

      ! Assume we only have a single field line
      ia = 1

      ! Calculate the denominators needed for electrostatic simulations
      if (fphi > epsilon(0.0)) then
      
         ! Allocate temporary arrays
         if (debug) write(*, *) 'field_equations_fluxtube::init_field_equations_fluxtube::init_denominator_fields'
         allocate (g0(nvpa, nmu))


         !----------------------------------------------------------------------
         !--------------- Guiding-center distribution function g ---------------
         !----------------------------------------------------------------------
         ! If we use the guiding-center distribution function g then
         !     sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g + (Z_s n_s/T_s) (Gamma0 - 1) phi ] = 0
         !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] / [ sum_s (Z_s² n_s/T_s) (1 - Gamma0) ]
         !
         !     denominator_fields[iky,ikz,iz] = sum_s (Z_s² n_s/T_s) (1 - Gamma0)
         !
         !     (1 - Gamma0(b_k) = (2B/sqrt(pi)) int dvpa int dmu (1 - J_0(a_k)²) exp(v²)
         !----------------------------------------------------------------------
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            
            ! <denominator_fields> does not depend on flux tube index, so only compute for one flux tube index
            it = it_idx(kxkyz_lo, ikxkyz)
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! Calculate (1 - J_0(a_k)²) exp(v²) for each (kx,ky,z)
            g0 = spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa)

            ! Multiply by the Maxwellian if needed
            if (.not. maxwellian_normalization) then
               g0 = g0 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
            end if

            ! Calculate denominator_fields[iky,ikz,iz] = sum_s (Z_s² n_s/T_s) (1 - Gamma0)
            ! with (1 - Gamma0(b_k) = (2B/sqrt(pi)) int dvpa int dmu (1 - J_0(a_k)²) exp(v²)
            wgt = spec(is)%z * spec(is)%z * spec(is)%dens_psi0 / spec(is)%temp
            call integrate_vmu(g0, iz, tmp)
            denominator_fields(iky, ikx, iz) = denominator_fields(iky, ikx, iz) + tmp * wgt
            
         end do
         
         ! Sum the values on all processors and send them to <proc0>
         call sum_allreduce(denominator_fields)
         
         ! Avoid divide by zero when kx=ky=0; We do not evolve this mode, so the value is irrelevant
         if (zonal_mode(1) .and. akx(1) < epsilon(0.) .and. has_electron_species(spec)) then
            denominator_fields(1, 1, :) = 0.0
         end if


         !----------------------------------------------------------------------
         !           Non-adiabatic part h of the distribution function 
         !----------------------------------------------------------------------
         ! If we are using the non-adiabatic part h of the distribution function then
         !     sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * h - (Z_s n_s/T_s) phi ] = 0
         !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * h ] / [ sum_s (Z_s² n_s/T_s) ]
         !
         !     denominator_fields_h = sum_s (Z_s² n_s/T_s)
         !----------------------------------------------------------------------
         denominator_fields_h = sum(spec%z * spec%z * spec%dens / spec%temp)
         

         !**********************************************************************
         !                             Adiabatic electrons 
         !**********************************************************************
         ! If using adiabatic electrons then we need to modify phi by adding a 
         ! Boltzmann response. We do not explicitly have the electron distribution
         ! function, and use 
         !                 g_e = e <phi>_R * F_0e / T_e + const. 
         ! where the const. is chosen: 
         !     - to be zero if adabatic (no field line average) is used.
         !     - is set such that the flux-surface average of g_e is zero if a 
         !       modified Boltzmann response is used.
         !
         ! In order for this to be satisfied we find 
         !                       \delta n_e = e n_e (phi - <phi>_R) / T_e
         ! where <phi>_FSA is the flux surface average of phi. For a flux tube simulation
         ! this reduces to a field-line average.
         !
         !
         ! ----------------------------------------------------------------------
         !           Use g as distribution function - no field line average 
         ! ----------------------------------------------------------------------
         ! sum_(s not e) Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g + (Gamma0 - 1) (Z_s/T_s) phi ] 
         !                     = n_e/T_e phi
         !
         ! Re-arrange for phi:
         !    phi = sum_(s not e) Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] / [sum_(s not e) (Z_s² n_s/T_s) (1 - Gamma0) + n_e/T_e]
         !
         ! Re-define the denominator:
         !     denominator_fields = sum_(s not e) (Z_s² n_s/T_s) (1 - Gamma0) + n_e/T_e
         ! 
         ! ----------------------------------------------------------------------
         !           Use h as distribution function - no field line average 
         ! ----------------------------------------------------------------------
         ! sum_(s not e) Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * h - (Z_s/T_s) phi ] 
         !                     = n_e/T_e phi
         !
         ! Re-arrange for phi:
         !    phi = sum_(s not e) Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] / [sum_(s not e) (Z_s² n_s/T_s) + n_e/T_e]
         !
         ! Re-define the denominator:
         !     denominator_fields_h = sum_(s not e) (Z_s² n_s/T_s) + n_e/T_e
         !
         ! ----------------------------------------------------------------------
         !                   Use g as distribution function - MBR 
         ! ----------------------------------------------------------------------
         ! If we substitue delta n_e into QN, then the QN equation for g becomes: 
         !
         ! sum_(s not e) Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g + (Gamma0 - 1) (Z_s/T_s) phi ] 
         !                     = n_e/T_e (phi - <phi>_FSA)
         !
         ! Here the FSA corresponds to the ky = 0.0 mode, so this only modifies the zonal component of phi. 
         ! < . >_FSA =  [\int dy \int dl/B (.)] / [\int dy \int dl/B]
         ! this term has an assumed Kronecker delta in ky, so that it only affects the ky = 0.0 mode.
         !
         ! Firstly collect the contribution that affects all ky:
         !     sum_(s not e) Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g = 
         !        [sum_(s not e) (1 - Gamma0) * (Z_s² n_s/T_s) + n_e/T_e ] * phi - n_e/T_e <phi>_FSA 
         ! 
         ! Re-define the denominator:
         !     denominator_fields = sum_(s not e) (Z_s² n_s/T_s) (1 - Gamma0) + n_e/T_e
         !
         ! and divide everything by denominator_fields:
         !     phi - n_e/ T_e * <phi>_FSA / denominator_fields = 
         !              sum_(s not e) Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g] / denominator_fields
         !
         ! Now, we need to deal with the <phi>_FSA term. for this we FSA the whole QN equation:
         !     (1 - <n_e/ T_e / denominator_field>_FSA ) * <phi>_FSA  = 
         !              <sum_(s not e) Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g] / denominator_fields>_FSA
         ! 
         ! Define 
         !     denominator_fields_MBR = 1 / (T_e/n_e - <1/denominator_field>_FSA ) - 1
         ! 
         ! such that
         !     <phi>_FSA = [T_e/n_e*<sum_(s not e) Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g]/denominator_fields>_FSA]*(denominator_fields_MBR + 1)
         ! which is the ky = 0.0 component of phi.
         !
         ! Now we have everything we need to compute phi:
         ! For ky not 0.0:
         !     phi = sum_i Z_i n_i [ (2B/sqrt(pi)) int dvpa int dmu J_0 * gi ] / denominator_fields
         ! For ky = 0.0:
         !     phi = sum_i Z_i n_i [ (2B/sqrt(pi)) int dvpa int dmu J_0 * gi ] / denominator_fields
         !           - n_e / T_e * <sum_(s not e) Z_s n_s * <phi>_FSA * denominator_fields_MBR
         ! =>  phi = phi + <phi>_FSA * denominator_fields_MBR
         !
         ! ----------------------------------------------------------------------
         !                   Use h as distribution function - MBR 
         ! ----------------------------------------------------------------------
         ! If we substitue delta n_e into QN, then the QN equation for h becomes: 
         !
         ! sum_(s not e) Z_s n_s [(2B/sqrt(pi)) int dvpa int dmu J_0 * h_s - (Z_s /T_s) phi ] 
         !                                = n_e/T_e (phi - <phi>_FSA)
         ! 
         ! Here the FSA corresponds to the ky = 0.0 mode, so this only modifies the zonal component of phi.
         ! Rearrange: 
         !    sum_(s not e) Z_s n_s [(2B/sqrt(pi)) int dvpa int dmu J_0 * h_s] / [sum_(s not e) (Z_s² n_s/T_s) + (n_e/T_e)]
         !                               = phi - n_e/T_e * <phi>_FSA / [sum_(s not e) (Z_s² n_s/T_s) + (n_e/T_e)]
         !
         ! Re-define the denominator to include the adiabatic electron contribution:
         !     denominator_fields_h = sum_(s not e) (Z_s² n_s/T_s) + (n_e/T_e)
         ! 
         !     phi - n_e/T_e * <phi>_FSA / denominator_fields_h = 
         !              sum_(s not e) Z_s n_s [(2B/sqrt(pi)) int dvpa int dmu J_0 * h_s] / denominator_fields_h
         !
         ! Now, we need to deal with the <phi>_FSA term. for this we FSA the whole QN equation:
         !     (1 - <n_e/T_e / denominator_fields_h>_FSA ) * <phi>_FSA  =
         !              <sum_(s not e) Z_s n_s [(2B/sqrt(pi)) int dvpa int dmu J_0 * h_s] / denominator_fields_h>_FSA  
         ! 
         ! For the ky = 0.0 mode we also define: 
         !     denominator_fields_MBR_h = (T_e/n_e - <1/ denominator_fields_h>_FSA ) - 1
         !  => denominator_fields_MBR_h = n_e/T_e / sum_(s not e) (Z_s² n_s/T_s)
         !
         ! Then for ky not 0.0:
         !     phi = sum_(s not e) Z_s n_s [(2B/sqrt(pi)) int dvpa int dmu J_0 * h_s] / denominator_fields_h
         ! For ky = 0.0:
         !     phi = sum_(s not e) Z_s n_s [(2B/sqrt(pi)) int dvpa int dmu J_0 * h_s] / denominator_fields_h 
         !          - [n_e/T_e * <sum_(s not e) Z_s n_s [(2B/sqrt(pi)) int dvpa int dmu J_0 * h_s] / denominator_fields_h>_FSA] * denominator_fields_MBR_h
         ! =>  phi = phi + <phi>_FSA * denominator_fields_MBR_h
         ! ----------------------------------------------------------------------
         !     efac = n_e/T_e = (Ti/Te)/Ti / [ (ni/ne)/ni ] = <tite> / <nine> * ni / Ti
         !**********************************************************************

         if (.not. has_electron_species(spec)) then
            if (debug) write(*, *) 'field_equations_quasineutrality::fluxtube::init::init_denominator_fields_MBR'
            !---------------------------------------------------------------------
            !                 Adiabatic electrons - no field line average
            !---------------------------------------------------------------------
            ! denominator_fields = sum_(s not e) (Z_s² n_s/T_s) (1 - Gamma0) + n_e/T_e
            ! denominator_fields_h = sum_(s not e) (Z_s² n_s/T_s) + (n_e/T_e)
            !---------------------------------------------------------------------
            ! Calculate <efac> = n_e/T_e = <tite> / <nine> * ni / Ti
            efac = tite / nine * (spec(ion_species)%dens / spec(ion_species)%temp)
            efacp = efac * (spec(ion_species)%tprim - spec(ion_species)%fprim)
            
            ! Add the contribution of adiabatic electrons to <denominator_fields>
            ! Calculate denominator_fields = sum_(s not e) (Z_s² n_s/T_s) (1 - Gamma0) + (n_e/T_e)
            ! This is the factor that multiplies phi in the final expression. 
            denominator_fields = denominator_fields + efac
            
            ! Add the contribution of adiabatic electrons to <denominator_fields_h>
            ! Calculate denominator_fields_h = sum_(s not e) (Z_s² n_s/T_s) + (n_e/T_e)
            denominator_fields_h = denominator_fields_h + efac
            
            !---------------------------------------------------------------------
            !                 Modified Boltzmann Response (MBR)
            !---------------------------------------------------------------------
            ! For the Modified Boltzmann Response (MBR) we need to consider the flux-surface average
            ! contribution. 
            ! Use delta n_e = (n_e/T_e) ( phi - <phi>_FSA)
            ! 
            ! For h:
            !     denominator_fields_MBR_h = n_e/T_e / sum_(s not e) (Z_s² n_s/T_s)
            ! For g: 
            !     denominator_fields_MBR = 1/ (T_e/n_e - <1/denominator_field>_FSA )
            ! =>  denominator_fields_MBR = 1/ (T_e/n_e - <1/sum_(s not e) (Z_s² n_s/T_s) (1 - Gamma0) + (n_e/T_e)>_FSA)
            !---------------------------------------------------------------------
            if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
               if (zonal_mode(1)) then
                  ! denominator_fields_MBR_h = n_e/T_e / sum_(s not e) (Z_s² n_s/T_s)
                  denominator_fields_MBR_h = efac / (sum(spec%zt * spec%z * spec%dens))
                  do ikx = 1, nakx
                     ! tmp = T_e / n_e - int (dl/B)/(sum_(s not e) (Z_s² n_s/T_s) * (1- Gamma0))
                     tmp = 1./efac - sum(dl_over_b(ia, :) / denominator_fields(1, ikx, :))
                     ! denominator_fields_MBR = 1/ (T_e/n_e - <1/denominator_field>_FSA )
                     denominator_fields_MBR(ikx, :) = 1./(denominator_fields(1, ikx, :) * tmp)
                  end do
                  
                  ! Avoid dividing by zero for kx=ky=0 mode, which we do not need anyway
                  if (akx(1) < epsilon(0.)) then
                     denominator_fields_MBR(1, :) = 0.0
                  end if
               end if
            end if

         end if

         ! Deallocate temporary arrays
         if (allocated(g0)) deallocate (g0)

      end if

   end subroutine init_field_equations_fluxtube

end module field_equations_fluxtube
