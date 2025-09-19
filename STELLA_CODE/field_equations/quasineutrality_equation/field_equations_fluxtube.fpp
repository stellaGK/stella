!###############################################################################
!############## ADVANCE FIELDS USING THE QUASINEUTRALITY EQUATION ##############
!###############################################################################
! 
! Evolve the fields in time using the quasi-neutrality condition
!     sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g + (Zs/Ts) (Gamma0 - 1) phi ] = 0
!     sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * h - (Zs/Ts) phi ] = 0
! 
! Here we used the guiding-center disitribution function g and the non-adiabatic part h
!     g = <delta f>_theta = h_s - Zs/Ts <phi>_theta F0
! 
! The arguments of the Bessel functions are
!     J0(a_k) = J0(k_perp * rho_s)
!     Gamma0(b_k) = Gamma0(k_perp**2 * rho_s**2 / 2)
! 
! This equation can be rewritten in order to obtain the electrostatic potential:
!     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ] / [ sum_s (Zs²ns/Ts) (1 - Gamma0) ]
!     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * h ] / [ sum_s (Zs²ns/Ts) ]
! 
! The denominators are constants and are calculated when initialising stella
!     denominator_QN[iky,ikz,iz] = sum_s (Zs²ns/Ts) (1 - Gamma0)
!     denominator_QN_h = sum_s (Zs²ns/Ts)
! 
! The integral over velocity space and species is calculated in stella as
!     integrate_species( . ) = sum_s (2B/sqrt(pi)) int dvpa int dmu ( . )
! 
! To summarize, the fields (= electrostatic potential) can be calculated as
!     phi = integrate_species( J0 * g ) / denominator_QN
!     phi = integrate_species( J0 * h ) / denominator_QN_h
! 
!###############################################################################
module field_equations_fluxtube

   ! Load debug flags
   use debug_flags, only: debug => fields_fluxtube_debug
   
   implicit none

   ! Make routines available to other modules
   public :: advance_fields_using_field_equations_fluxtube
   public :: advance_fields_using_field_equations_quasineutrality
   public :: init_field_equations_fluxtube

   private

   ! Advance fields for g(kx,ky,z,ivpamus) and g(vpa,mu,ikxkyzs)
   interface advance_fields_using_field_equations_quasineutrality
      module procedure advance_fields_using_field_equations_quasineutrality_vmlo
      module procedure advance_fields_using_field_equations_quasineutrality_kxkyzlo
   end interface

contains

!###############################################################################
!############## ADVANCE FIELDS USING THE QUASINEUTRALITY EQUATION ##############
!###############################################################################

   !****************************************************************************
   !             ADVANCE FIELDS USING THE QUASINEUTRALITY EQUATION              
   !****************************************************************************
   ! The fields (electrostatic potential and electromagnetic fields) are evolved
   ! in time though the quasi-neutrality equation.
   ! 
   ! Note that Apar and Bpar are only advanced when using electromagnetic stella, 
   ! so these are in field_equations_electromagnetic.fpp
   !****************************************************************************
   subroutine advance_fields_using_field_equations_fluxtube(g, phi, apar, bpar, dist, skip_fsa)

      ! Parallelisation
      use mp, only: proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo
      use redistribute, only: scatter
      use calculations_redistribute, only: kxkyz2vmu
      
      ! Arrays
      use arrays_distribution_function, only: gvmu
      use arrays, only: time_field_solve
      
      ! Parameters
      use parallelisation_layouts, only: fields_kxkyz
      
      ! Grids
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(inout) :: phi, apar, bpar
      character(*), intent(in) :: dist
      logical, optional, intent(in) :: skip_fsa
      
      !-------------------------------------------------------------------------
      
      ! Note that fields_kxkyz = F is the default
      if (.not. fields_kxkyz) then
     
         ! This will call advance_fields_using_field_equations_quasineutrality_vmlo
         if (debug) write (*, *) 'field_equations_quasineutrality::fluxtube::advance_fields::vmlo'
         call advance_fields_using_field_equations_quasineutrality(g, phi, apar, bpar, dist, skip_fsa)
         
      ! Note that to use this option it has to be specified by the user
      ! First gather (vpa,mu) onto processor for v-space operations
      ! v-space operations are field solve, dg/dvpa, and collisions
      else if (fields_kxkyz) then 
      
         ! Scatter and time it
         if (debug) write (*, *) 'field_equations_quasineutrality::fluxtube::advance_fields::scatter'
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' field_equations_quasineutrality_redist')
         call scatter(kxkyz2vmu, g, gvmu)
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' field_equations_quasineutrality_redist')
         
         ! Given gvmu with vpa and mu local, calculate the corresponding fields
         ! This will call advance_fields_using_field_equations_quasineutrality_kxkyzlo
         if (debug) write (*, *) 'field_equations_quasineutrality::fluxtube::advance_fields::kxkyzlo'
         call advance_fields_using_field_equations_quasineutrality(gvmu, phi, apar, bpar, dist, skip_fsa)
         
      end if

   end subroutine advance_fields_using_field_equations_fluxtube

   !****************************************************************************
   !  ADVANCE FIELDS USING THE QUASINEUTRALITY EQUATION AND G(KX,KY,Z,IVPAMUS)  
   !****************************************************************************
   ! If <g> is parallelised over (vpa,mu,s) then this subroutine is called.
   ! This is the more common version used compared with parallelising over
   ! (kx,ky,z) and is the default for stella.
   ! 
   ! Here we calculate the electrostatic potential based on the quasi-neutrality equation:
   !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ] / [ sum_s (Zs²ns/Ts) (1 - Gamma0) ]
   !     phi = integrate_species( J0 * g ) / denominator_QN
   ! 
   ! TODO-GA: remove apar from this and make it only needed for EM stella
   !****************************************************************************
   subroutine advance_fields_using_field_equations_quasineutrality_vmlo(g, phi, apar, bpar, dist, skip_fsa)

      ! Parallelisation
      use mp, only: mp_abort, proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx
      
      ! Arrays
      use arrays_distribution_function, only: g_scratch
      use arrays, only: time_field_solve
      
      ! Parameters
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: radial_variation
      use parameters_physics, only: fphi
      
      ! Grids
      use grids_z, only: nzgrid
      use grids_species, only: spec
      use calculations_velocity_integrals, only: integrate_species
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      
      ! Routines from other field modules
      use field_equations_electromagnetic, only: advance_fields_using_QN_electromagnetic
      use field_equations_quasineutrality_radial_variation, only: add_radial_correction_int_species
      use field_equations_quasineutrality_radial_variation, only: calculate_phi_for_radial_variation
      
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
      if (debug) write (*, *) 'field_equations_quasineutrality::fluxtube::advance_fields::vmlo'
      
      ! Used for the Dougherty collision operator
      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ! Initialise the electrostatic potential phi
      phi = 0.
      
      ! Note that this advances phi for the electrostatic, fluxtube case.
      ! If electromagnetic effects are included then phi will be advanced below
      if (fphi > epsilon(0.0) .and. .not. include_bpar) then
      
         ! Start timer
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         
         ! First gyro-average the distribution function g at each phase space location
         ! and store this as g_scratch = <g>_theta = J_0 g in k-space
         call gyro_average(g, g_scratch)
         
         ! If we are allowing for radial variation then we must modify <g>.
         ! This is done in the field_equations_radialvariation module, but 
         ! is not needed for standard stella (as these are false by default).
         if (radial_variation) call add_radial_correction_int_species(g_scratch)
         
         ! Calculate phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ] = integrate_species( J0 * g )
         if (debug) write (*, *) 'field_equations_quasineutrality::fluxtube::vmulo::integrate_species_phi'
         call integrate_species(g_scratch, spec%z * spec%dens_psi0, phi)
         
         ! Stop timer
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         
         ! Calculate phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ] / [ sum_s (Zs²ns/Ts) (1 - Gamma0) ]
         ! by dividing with denominator_QN[iky,ikz,iz] = sum_s (Zs²ns/Ts) (1 - Gamma0) in the calculate_phi() routine
         if (.not. radial_variation) then
            if (debug) write (*, *) 'field_equations_quasineutrality::fluxtube::vmulo::calculate_phi'
            call calculate_phi(phi, dist, skip_fsa_local)
            
         ! The denominator_QN[iky,ikz,iz] factor is modified by radial varation effects
         else if (radial_variation) then
            if (debug) write (*, *) 'field_equations_quasineutrality::fluxtube::vmulo::calculate_phi_for_radial_variation'
            call calculate_phi_for_radial_variation(phi, dist , skip_fsa)
            
         end if
         
      end if

      ! If we have electromagnetic effects we need to calculate the fields A_parallel and B_parallel
      if (include_apar .or. include_bpar) then
         if (debug) write (*, *) 'field_equations_quasineutrality::fluxtube::vmulo::advance_fields::electromagnetic'
         apar = 0.; bpar = 0.
         call advance_fields_using_QN_electromagnetic(g, phi, apar, bpar, dist)
      end if

   end subroutine advance_fields_using_field_equations_quasineutrality_vmlo

   !****************************************************************************
   !  ADVANCE FIELDS USING THE QUASINEUTRALITY EQUATION AND G(MU,VPA,IKXKYZS)  
   !****************************************************************************
   ! If <g> is parallelised over (kx,ky,z,s) then this subroutine is called.
   ! 
   ! Here we calculate the electrostatic potential based on the quasi-neutrality equation:
   !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ] / [ sum_s (Zs²ns/Ts) (1 - Gamma0) ]
   !     phi = integrate_species( J0 * g ) / denominator_QN
   ! 
   ! TODO-GA: remove apar from this and make it only needed for EM stella
   !****************************************************************************
   subroutine advance_fields_using_field_equations_quasineutrality_kxkyzlo(g, phi, apar, bpar, dist, skip_fsa)

      ! Parallelisation
      use mp, only: proc0
      use mp, only: sum_allreduce, mp_abort
      use job_manage, only: time_message
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      
      ! Arrays
      use arrays, only: time_field_solve
      
      ! Parameters
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: fphi
      use parameters_physics, only: radial_variation
      
      ! Grids
      use grids_velocity, only: nvpa, nmu
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_species, only: spec
      use grids_z, only: nzgrid
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      
      ! Routines from other field modules
      use field_equations_electromagnetic, only: advance_fields_using_QN_electromagnetic
      use field_equations_quasineutrality_radial_variation, only: calculate_phi_for_radial_variation

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
      if (debug) write (*, *) 'field_equations_quasineutrality::fluxtube::kxkyzlo'
      
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
      
         ! Start timer
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         if (debug) write (*, *) 'fields::get_field_equations_fluxtube_kxkyzlo_int_dv_g'
         
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
            ! and store this as g0 = <g>_theta = J_0 g
            call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
            
            ! Calculate phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ]
            wgt = spec(is)%z * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp
            
         end do
         
         ! Deallocate temporary arrays
         deallocate (g0)
         
         ! Sum the values on all processors and send them to <proc0>
         call sum_allreduce(phi)

         ! Calculate phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ] / [ sum_s (Zs²ns/Ts) (1 - Gamma0) ]
         ! by dividing with denominator_QN[iky,ikz,iz] = sum_s (Zs²ns/Ts) (1 - Gamma0) in the calculate_phi() routine
         if (.not. radial_variation) then 
            if (debug) write (*, *) 'field_equations_quasineutrality::fluxtube::kxkyzlo::calculate_phi'
            call calculate_phi(phi, dist, skip_fsa_local)

         ! The denominator_QN[iky,ikz,iz] factor is modified by radial varation effects
         else if (radial_variation) then
            if (debug) write (*, *) 'field_equations_quasineutrality::fluxtube::kxkyzlo::calculate_phi_for_radial_variation'
            call calculate_phi_for_radial_variation (phi, dist , skip_fsa)
         end if
         
      end if

      ! If we have electromagnetic effects we need to calculate the fields A_parallel and B_parallel
      if (include_apar .or. include_bpar) then 
         if (debug) write (*, *) 'field_equations_quasineutrality::fluxtube::kxkyzlo::advance_fields_electromagnetic'
         bpar = 0.; apar = 0.
         call advance_fields_using_QN_electromagnetic(g, phi, apar, bpar, dist)
      end if

   end subroutine advance_fields_using_field_equations_quasineutrality_kxkyzlo

   !****************************************************************************
   !******************************* CALCULATE PHI ******************************
   !****************************************************************************
   ! Here we calculate the electrostatic potential based on the quasi-neutrality equation:
   !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ] / [ sum_s (Zs²ns/Ts) (1 - Gamma0) ]
   !     phi = integrate_species( J0 * g ) / denominator_QN
   ! 
   ! Note that the 'phi' variable passed in is:
   !    integrate_species( J0 * g ) = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ]
   ! 
   ! This routine divides by the appropriate <denominator_QN> factor depending on if we
   ! have kinetic or adiabatic electrons, and also on whether we are using 'g'
   ! or 'h' as our distribution function that we are evolving.
   !****************************************************************************
   subroutine calculate_phi(phi, dist, skip_fsa)

      ! Parallelisation
      use mp, only: proc0, mp_abort
      use job_manage, only: time_message
      use multibox, only: mb_calculate_phi
      
      ! Arrays
      use arrays, only: denominator_QN
      use arrays, only: denominator_QN_MBR
      use arrays, only: denominator_QN_h
      use arrays, only: denominator_QN_MBR_h
      use arrays, only: time_field_solve
      
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
      real, dimension(:, :, :, :), allocatable :: denominator_QN_t
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
      
      ! If we are using the non-adiabatic part h of the distribution function then
      ! sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * h - (Zs/Ts) phi ] = 0
      ! phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * h ] / [ sum_s (Zs²ns/Ts) ]
      ! denominator_QN_h = sum_s (Zs²ns/Ts)
      if (dist == 'h') then
         if (debug) write(*, *) 'field_equations_quasineutrality::fluxtube::calculate_phi::dist==h'
         phi = phi / denominator_QN_h
         
      ! If we are using the guiding-center distribution function g then
      ! sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g + (Zs/Ts) (Gamma0 - 1) phi ] = 0
      ! phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ] / [ sum_s (Zs²ns/Ts) (1 - Gamma0) ]
      ! denominator_QN[iky,ikz,iz] = sum_s (Zs²ns/Ts) (1 - Gamma0)
      else if (dist == 'g' .or. dist == 'gbar') then
         if (debug) write(*, *) 'field_equations_quasineutrality::fluxtube::calculate_phi::dist==gbar'
         allocate (denominator_QN_t(naky, nakx, -nzgrid:nzgrid, ntubes))
         denominator_QN_t = spread(denominator_QN, 4, ntubes)
         where (denominator_QN_t < epsilon(0.0))
            phi = 0.0
         elsewhere
            phi = phi / denominator_QN_t
         end where
         deallocate (denominator_QN_t)
         
      ! Abort if <dist> is not recognized.
      else
         if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
         call mp_abort('unknown dist option in get_fields. aborting')
         return
      end if

      ! The kx = ky = 0.0 mode is not evolved by stella so make sure this term is set to zero.
      if (debug) write(*, *) 'field_equations_quasineutrality::fluxtube::calculate_phi::set kxky=0.0 to zero'
      if (any(denominator_QN(1, 1, :) < epsilon(0.))) phi(1, 1, :, :) = 0.0
      if (proc0) call time_message(.false., time_field_solve(:, 4), ' calculate_phi')

      ! Handle adiabatic electrons only if needed.
      ! TODO - write notes!
      if (proc0) call time_message(.false., time_field_solve(:, 5), 'calculate_phi_adia_elec')
      if (adia_elec .and. zonal_mode(1) .and. .not. skip_fsa_local) then
         if (dist == 'h') then
            if (debug) write(*, *) 'field_equations_quasineutrality::fluxtube::calculate_phi::adiabatic_electrons::dist==h'
            do it = 1, ntubes
               do ikx = 1, nakx
                  tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                  phi(1, ikx, :, it) = phi(1, ikx, :, it) + tmp * denominator_QN_MBR_h
               end do
            end do
         else if (dist == 'g' .or. dist == 'gbar') then
            if (debug) write(*, *) 'field_equations_quasineutrality::fluxtube::calculate_phi::adiabatic_electrons::dist==gbar'
            do ikx = 1, nakx
               do it = 1, ntubes
                  tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                  phi(1, ikx, :, it) = phi(1, ikx, :, it) + tmp * denominator_QN_MBR(ikx, :)
               end do
            end do
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
   ! The electrostatic potential phi is calculated based on the quasi-neutrality condition
   !     sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g + (Zs/Ts) (Gamma0 - 1) phi ] = 0
   !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ] / [ sum_s (Zs²ns/Ts) (1 - Gamma0) ]
   !     denominator_QN[iky,ikz,iz] = sum_s (Zs²ns/Ts) (1 - Gamma0)
   ! The denominators needed to calculate <phi> are initialised in this routine.
   ! 
   ! If adiabatic electrons are used then this factor is modified and we use 
   ! <denominator_QN_MBR> which includes the Modified Boltmann Response.
   ! 
   ! Use the following calculations:
   !      (1 - Gamma0(b_k) = (2B/sqrt(pi)) int dvpa int dmu (1 - J0(a_k)²) exp(v²)
   !      integrate_vmu( . ) = (2B/sqrt(pi)) int dvpa int dmu ( . )
   !      b_k = k²_perp*rho²_s/2
   !      a_k = k_perp*rho_s
   !****************************************************************************
   subroutine init_field_equations_fluxtube

      ! Parallelisation
      use mp, only: sum_allreduce
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, onlY: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      
      ! Arrays
      use arrays, only: denominator_QN, denominator_QN_MBR
      use arrays, only: denominator_QN_h, denominator_QN_MBR_h, efac, efacp
      
      ! Parameters
      use parameters_physics, only: fphi
      
      ! Grids
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: spec
      use grids_kxky, only: zonal_mode, akx
      use grids_kxky, only: nakx
      
      ! Calculations
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use parameters_numerical, only: maxwellian_normalization
      use arrays_gyro_averages, only: aj0v
      
      ! Adiabatic electrons
      use grids_species, only: has_electron_species
      use grids_species, only: ion_species
      use grids_species, only: tite, nine
      use grids_species, only: adiabatic_option_switch
      use grids_species, only: adiabatic_option_fieldlineavg
      
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
         if (debug) write(*, *) 'field_equations_fluxtube::init_field_equations_fluxtube::init_denominator_QN'
         allocate (g0(nvpa, nmu))


         !----------------------------------------------------------------------
         !--------------- Guiding-center distribution function g ---------------
         !----------------------------------------------------------------------
         ! If we use the guiding-center distribution function g then
         !     sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g + (Zs/Ts) (Gamma0 - 1) phi ] = 0
         !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ] / [ sum_s (Zs²ns/Ts) (1 - Gamma0) ]
         !     denominator_QN[iky,ikz,iz] = sum_s (Zs²ns/Ts) (1 - Gamma0)
         !     (1 - Gamma0(b_k) = (2B/sqrt(pi)) int dvpa int dmu (1 - J0(a_k)²) exp(v²)
         !----------------------------------------------------------------------
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            
            ! <denominator_QN> does not depend on flux tube index, so only compute for one flux tube index
            it = it_idx(kxkyz_lo, ikxkyz)
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! Calculate (1 - J0(a_k)²) exp(v²) for each (kx,ky,z)
            ! TODO-GA: remove maxwellian_normalisation flag
            g0 = spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa)
            if (.not. maxwellian_normalization) then
               g0 = g0 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
            end if

            ! Calculate denominator_QN[iky,ikz,iz] = sum_s (Zs²ns/Ts) (1 - Gamma0)
            ! with (1 - Gamma0(b_k) = (2B/sqrt(pi)) int dvpa int dmu (1 - J0(a_k)²) exp(v²)
            wgt = spec(is)%z * spec(is)%z * spec(is)%dens_psi0 / spec(is)%temp
            call integrate_vmu(g0, iz, tmp)
            denominator_QN(iky, ikx, iz) = denominator_QN(iky, ikx, iz) + tmp * wgt
            
         end do
         
         ! Sum the values on all processors and send them to <proc0>
         call sum_allreduce(denominator_QN)
         
         ! Avoid divide by zero when kx=ky=0; We do not evolve this mode, so the value is irrelevant
         if (zonal_mode(1) .and. akx(1) < epsilon(0.) .and. has_electron_species(spec)) then
            denominator_QN(1, 1, :) = 0.0
         end if


         !----------------------------------------------------------------------
         !---------- Non-adiabatic part h of the distribution function ---------
         !----------------------------------------------------------------------
         ! If we are using the non-adiabatic part h of the distribution function then
         !     sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * h - (Zs/Ts) phi ] = 0
         !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * h ] / [ sum_s (Zs²ns/Ts) ]
         !     denominator_QN_h = sum_s (Zs²ns/Ts)
         !----------------------------------------------------------------------
         denominator_QN_h = sum(spec%z * spec%z * spec%dens / spec%temp)


         !----------------------------------------------------------------------
         !-------- Adiabatic electrons use a Modified Boltmann Response --------
         !----------------------------------------------------------------------
         ! If the electrons are treated adiabatic then we assume that the non-adiabatic part, h, is negligible.
         !     ge = he - (Ze/Te) <phi>_theta F0e = - (Ze/Te) <phi>_theta F0e
         ! As a result, the perturbation of the electron density in the quasi-neutrality condition is given by,
         !     delta n_e = - (Ze/Te) n_e phi = (n_e/T_e) phi
         ! If we are using the non-adiabatic part h of the distribution function then
         !     sum_i Z_i n_i [ (2B/sqrt(pi)) int dvpa int dmu J0 * hi - (Zi/Ti) phi ] - (n_e/T_e) phi = 0
         !     phi = sum_i Z_i n_i [ (2B/sqrt(pi)) int dvpa int dmu J0 * hi ] / [ sum_i (Zi²ni/Ti) + (n_e/T_e) ]
         !     denominator_QN_h = sum_i (Zi²ni/Ti) + (n_e/T_e)
         !     efac = n_e/T_e = (Ti/Te)/Ti / [ (ni/ne)/ni ] = <tite> / <nine> * ni / Ti
         !----------------------------------------------------------------------
         if (.not. has_electron_species(spec)) then
            if (debug) write(*, *) 'field_equations_quasineutrality::fluxtube::init::init_denominator_QN_MBR'
            
            ! Calculate <efac> = n_e/T_e = <tite> / <nine> * ni / Ti
            efac = tite / nine * (spec(ion_species)%dens / spec(ion_species)%temp)
            efacp = efac * (spec(ion_species)%tprim - spec(ion_species)%fprim)
            
            ! Add the contribution of adiabatic electrons to <denominator_QN>
            ! Calculate denominator_QN = sum_i (Zi²ni/Ti) (1 - Gamma0) + (n_e/T_e)
            denominator_QN = denominator_QN + efac
            
            ! Add the contribution of adiabatic electrons to <denominator_QN_h>
            ! Calculate denominator_QN_h = sum_i (Zi²ni/Ti) + (n_e/T_e)
            denominator_QN_h = denominator_QN_h + efac
            
            ! For the Modified Boltzmann Response (MBR) we replace the perturbed electron density by
            !     delta n_e = (n_e/T_e) ( phi - <phi>_FSA) = (n_e/T_e) ( phi - int dell/B phi)
            ! With < . >_FSA the flux-surface-average, which reduces to a field-line average in the flux-tube approximation.
            ! If we are using the non-adiabatic part h of the distribution function then ... TODO - write notes!
            !     phi = integrate_species(J0 * hi) / denominator_QN_h + (n_e/T_e) / sum_i (Zi²ni/Ti) * int dell/B phi
            !     denominator_QN_MBR_h = (n_e/T_e) / sum_i (Zi²ni/Ti)
            !     denominator_QN_MBR = ??? ... TODO - write notes!
            if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
               if (zonal_mode(1)) then
                  denominator_QN_MBR_h = efac / (sum(spec%zt * spec%z * spec%dens))
                  do ikx = 1, nakx
                     tmp = 1./efac - sum(dl_over_b(ia, :) / denominator_QN(1, ikx, :))
                     denominator_QN_MBR(ikx, :) = 1./(denominator_QN(1, ikx, :) * tmp)
                  end do
                  
                  ! Avoid dividing by zero for kx=ky=0 mode, which we do not need anyway
                  if (akx(1) < epsilon(0.)) then
                     denominator_QN_MBR(1, :) = 0.0
                  end if
               end if
            end if
            
         end if

         ! Deallocate temporary arrays
         if (allocated(g0)) deallocate (g0)

      end if

   end subroutine init_field_equations_fluxtube

end module field_equations_fluxtube
