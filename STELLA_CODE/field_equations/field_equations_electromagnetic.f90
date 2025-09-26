!###############################################################################
!###################### ADVANCE ELECTROMAGNETIC FIELDS #########################
!###############################################################################
! 
! Module for advancing and initialising the fields when electromagnetic 
! effects are included, i.e., when evolving the apar and bpar fields.
! 
! Here, the addition to QN is solved for, when bpar is included.
! Parallel and perpendicular Ampere's Law are also solved.
!###############################################################################
module field_equations_electromagnetic

   ! Load debug flags
   use debug_flags, only: debug => fields_debug

   implicit none

   ! Make routines available to other modules
   public :: init_field_equations_electromagnetic
   public :: allocate_field_equations_electromagnetic
   public :: finish_field_equations_electromagnetic
   public :: advance_fields_electromagnetic
   public :: advance_apar
   
   private

   interface advance_fields_electromagnetic
      module procedure advance_fields_electromagnetic_kxkyzlo
      module procedure advance_fields_electromagnetic_vmulo
   end interface

   real, dimension(:, :, :), allocatable :: denominator_fields13, denominator_fields31, denominator_fields33
contains

!###############################################################################
!####################### ADVANCE ELECTROMAGNETIC FIELDS ########################
!###############################################################################

   !****************************************************************************
   !                       ADVANCE ELECTROMAGNETIC FIELDS    
   !****************************************************************************
   ! The fields (electrostatic potential and electromagnetic fields) are evolved
   ! in time though the quasineutrality equation and Ampere's law.
   ! 
   ! If we are parallelising over (vpa,mu) then this subroutine is called
   ! This is the more common version used compared with parallelising over 
   ! (kx,ky,z) and is the default for stella.
   ! 
   ! This advances the fields when Electromagnetic effects are included, so 
   ! we advance <phi>, <B_parallel>, and <A_parallel>.
   !
   ! Currently, stella only supports electromagnetic effects when simulating
   ! a flux-tube. This module is therefore limited to flux-tube electromagnetic
   ! effects. 
   !****************************************************************************
   subroutine advance_fields_electromagnetic_vmulo(g, phi, apar, bpar, dist)

      ! Parallelisation
      use mp, only: proc0, mp_abort
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx
      
      ! Arrays
      use arrays_distribution_function, only: g_scratch
      use arrays, only: time_field_solve
      
      ! Parameters
      use parameters_physics, only: beta
      use parameters_physics, only: fphi
      use parameters_physics, only: radial_variation
      use parameters_physics, only: include_apar, include_bpar
      
      ! Grids
      use grids_species, only: spec
      use grids_velocity, only: vpa, mu
      use grids_z, only: nzgrid
      
      ! Calculations
      use calculations_velocity_integrals, only: integrate_species
      use calculations_gyro_averages, only: gyro_average
      use calculations_gyro_averages, only: gyro_average_j1
      
      ! Routines from other quasineutrality modules
      use field_equations_radialvariation, only: add_radial_correction_int_species

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar, bpar
      character(*), intent(in) :: dist

      ! Local variables
      integer :: imu, iv, ivmu 
      
      !-------------------------------------------------------------------------
      
      ! Calculate the phi and bpar fields
      if (fphi > epsilon(0.0) .and. include_bpar) then
      
         ! Start timer
         if (debug) write (*, *) 'field_equations_electromagnetic::advance_fields::vmulo::include_bpar'
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g int_dv_g_vperp2')
         
         ! Gyroaverage the distribution function g at each phase space location
         call gyro_average(g, g_scratch)
         
         ! <g> requires modification if radial profile variation is included
         if (radial_variation) call add_radial_correction_int_species(g_scratch)
         
         ! Integrate <g> over velocity space and sum over species
         ! store result in phi, which will be further modified below to account for polarization term
         if (debug) write (*, *) 'field_equations_electromagnetic::advance_fields::vmulo::integrate_species_phi'
         call integrate_species(g_scratch, spec%z * spec%dens_psi0, phi)
         
         ! Gyroaverage the distribution function g at each phase space location
         call gyro_average_j1(g, g_scratch)
         
         ! Multiply by mu factor from vperp2
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            imu = imu_idx(vmu_lo, ivmu)
            g_scratch(:, :, :, :, ivmu) = g_scratch(:, :, :, :, ivmu) * mu(imu)
         end do
         
         ! <g> requires modification if radial profile variation is included; not supported for bpar MRH
         ! Integrate <g> over velocity space and sum over species
         ! store result in bpar, which will be further modified below to account for polarization term
         if (debug) write (*, *) 'field_equations_electromagnetic::advance_fields::vmulo::integrate_species_bpar'
         call integrate_species(g_scratch, -2.0 * beta * spec%temp_psi0 * spec%dens_psi0, bpar)
         
         ! End timer
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g int_dv_g_vperp2')
         
         ! Get phi and bpar
         call calculate_phi_and_bpar(phi, bpar, dist)
         
      end if

      ! Initialise the apar field
      apar = 0.
      
      ! Evolve the apar field
      if (include_apar) then
      
         ! Start timer
         if (debug) write (*, *) 'field_equations_electromagnetic::advance_fields::vmulo::include_apar'
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         
         ! If fphi > 0, then g_scratch = <g> already calculated above
         call gyro_average(g, g_scratch)
         
         ! For parallel Amperes Law, need to calculate parallel current rather than density,
         ! so multiply <g> by vpa before integrating. First get the vpa index, then multiply with vpa
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            g_scratch(:, :, :, :, ivmu) = g_scratch(:, :, :, :, ivmu) * vpa(iv)
         end do
         
         ! Integrate vpa*<g> over velocity space and sum over species
         ! store result in apar, which will be further modified below to account for apar pre-factor
         if (debug) write (*, *) 'field_equations_electromagnetic::advance_fields::vmulo::integrate_species_apar'
         call integrate_species(g_scratch, spec%z * spec%dens_psi0 * spec%stm_psi0 * beta, apar)
         
         ! End timer
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         
         ! Divide the apar obtained above by the appropriate Apar pre-factor;
         ! this is just kperp2 if g = <f> is used or apar_denom = (kperp2 + ...)
         ! if gbar = g + <vpa*apar/c> * Ze/T * F_0 is used
         call get_apar(apar, dist)
         
      end if

   end subroutine advance_fields_electromagnetic_vmulo

   !****************************************************************************
   !**************************** GET FIELDS KXKYZLO ****************************
   !****************************************************************************
   subroutine advance_fields_electromagnetic_kxkyzlo(g, phi, apar, bpar, dist)

      ! Parallelisation
      use mp, only: proc0
      use mp, only: sum_allreduce, mp_abort
      use job_manage, only: time_message
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      
      ! Arrays
      use arrays, only: kperp2 
      use arrays, only: apar_denom, time_field_solve
      
      ! Parameters
      use parameters_physics, only: beta
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: fphi
     
      ! Grids
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: vpa, mu
      use grids_species, only: spec
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      use calculations_velocity_integrals, only: integrate_vmu

      implicit none

      ! Arguments
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(inout) :: phi, apar, bpar
      character(*), intent(in) :: dist
      
      ! Local variables
      complex :: tmp
      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      
      !-------------------------------------------------------------------------
      
      ! Assume we only have one field line
      ia = 1

      ! Calculate the phi and bpar fields
      if (fphi > epsilon(0.0) .and. include_bpar) then
      
         ! Start timer
         if (debug) write (*, *) 'field_equations_electromagnetic::advance_fields::kxkyzlo::include_bpar'
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g int_dv_g_vperp2')
         
         ! Allocate temporary arrays
         allocate (g0(nvpa, nmu))
          
         ! Iterate over the (kx,ky,z,mu,vpa,s) points
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            
            ! Integrate g to get sum_s Z_s n_s J0 g and store in phi
            call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
            wgt = spec(is)%z * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp
            
            ! Integrate g to get - 2 beta sum_s n_s T_s J1 mu g and store in bpar
            call gyro_average_j1(spread(mu, 1, nvpa) * g(:, :, ikxkyz), ikxkyz, g0)
            wgt = -2.0 * beta* spec(is)%z * spec(is)%dens_psi0 * spec(is)%temp_psi0
            call integrate_vmu(g0, iz, tmp)
            bpar(iky, ikx, iz, it) = bpar(iky, ikx, iz, it) + wgt * tmp
            
         end do
         
         ! Deallocate temporary arrays
         deallocate (g0)
         
         ! Sum the values on all processors and send them to <proc0>
         call sum_allreduce(phi)
         call sum_allreduce(bpar)
         
         ! End timer
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g int_dv_g_vperp2')

         ! Get phi and bpar
         call calculate_phi_and_bpar(phi, bpar, dist)
         
      end if

      ! Initialise the apar field
      apar = 0.
      
      ! Evolve the apar field
      if (include_apar) then
      
         ! Debug 
         if (debug) write (*, *) 'field_equations_electromagnetic::advance_fields::kxkyzlo::include_apar'
         
         ! Allocate temporary arrays
         allocate (g0(nvpa, nmu))
         
         ! Iterate over the (kx,ky,z,mu,vpa,s) points
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            
            call gyro_average(spread(vpa, 2, nmu) * g(:, :, ikxkyz), ikxkyz, g0)
            wgt = 2.0 * beta * spec(is)%z * spec(is)%dens * spec(is)%stm
            call integrate_vmu(g0, iz, tmp)
            apar(iky, ikx, iz, it) = apar(iky, ikx, iz, it) + tmp * wgt
            
         end do
         
         ! Sum the values on all processors and send them to <proc0>
         call sum_allreduce(apar)
         
         if (dist == 'h') then
            apar = apar / spread(kperp2(:, :, ia, :), 4, ntubes)
         else if (dist == 'gbar') then
            apar = apar / spread(apar_denom, 4, ntubes)
         else if (dist == 'gstar') then
            write (*, *) 'APAR NOT SETUP FOR GSTAR YET. aborting.'
            call mp_abort('APAR NOT SETUP FOR GSTAR YET. aborting.')
         else
            if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
            call mp_abort('unknown dist option in get_fields. aborting')
         end if
         
         ! Deallocate temporary arrays
         deallocate (g0)
         
      end if

   end subroutine advance_fields_electromagnetic_kxkyzlo

   !****************************************************************************
   !**************************** GET APAR AND BPAR *****************************
   !****************************************************************************
   subroutine calculate_phi_and_bpar(phi, bpar, dist)

      ! Parallelisation
      use mp, only: proc0, mp_abort
      use job_manage, only: time_message
      
      ! Arrays
      use arrays, only: denominator_fields_inv11, denominator_fields_inv13, denominator_fields_inv33, denominator_fields_inv31
      use arrays, only: denominator_fields_h, time_field_solve

      ! Grids
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: nakx, naky

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, bpar
      character(*), intent(in) :: dist
      
      ! Local variables
      integer :: ia, it, ikx, iky, iz
      complex :: antot1, antot3
      
      !-------------------------------------------------------------------------
      
      ! Start timer
      if (debug) write (*, *) 'field_equations_electromagnetic::calculate_phi_and_bpar'
      if (proc0) call time_message(.false., time_field_solve(:, 4), ' calculate_phi_and_bpar')

      ! Assume we only have one field line
      ia = 1
      
      if (dist == 'gbar' .or. dist == 'g') then
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  do iky = 1, naky
                     antot1 = phi(iky,ikx,iz,it)
                     antot3 = bpar(iky,ikx,iz,it)
                     phi(iky,ikx,iz,it) = denominator_fields_inv11(iky,ikx,iz)*antot1 + denominator_fields_inv13(iky,ikx,iz)*antot3
                     bpar(iky,ikx,iz,it) = denominator_fields_inv31(iky,ikx,iz)*antot1 + denominator_fields_inv33(iky,ikx,iz)*antot3
                  end do
               end do
            end do
         end do
         
      else if (dist == 'h') then
         ! divide sum ( Zs int J0 h d^3 v) by sum(Zs^2 ns / Ts)
         phi = phi / denominator_fields_h
         ! do nothing for bpar because
         ! bpar = - 2 * beta * sum(Ts ns int (J1/bs) mu h d^3 v)
         ! which is already stored in bpar when dist = 'h'.
         
      else
         if (proc0) write (*, *) 'Unknown dist option in get_fields. Aborting.'
         call mp_abort('Unknown dist option in get_fields. Aborting.')
         return
         
      end if

   end subroutine calculate_phi_and_bpar

   !****************************************************************************
   !                                    GET APAR
   !****************************************************************************
   ! Get_apar solves pre-factor * Apar = beta_ref * sum_s Z_s n_s vth_s int d3v vpa * J0 * pdf
   ! for apar, with pdf being either g or gbar (specified by dist input).
   ! the input apar is the RHS of the above equation and is overwritten by the true apar
   ! the pre-factor depends on whether g or gbar is used (kperp2 in former case, with additional
   ! term appearing in latter case)
   !****************************************************************************
   subroutine get_apar(apar, dist)

      ! Parallelisation
      use mp, only: proc0, mp_abort
      
      ! Arrays
      use arrays, only: kperp2
      use arrays, only: apar_denom
      
      ! Grids
      use grids_z, only: nzgrid, ntubes
      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: apar
      character(*), intent(in) :: dist

      ! Local variables
      integer :: ia
      
      !-------------------------------------------------------------------------
      
      ! Assume we only have one field line
      ia = 1
      
      if (dist == 'g') then
         where (spread(kperp2(:, :, ia, :), 4, ntubes) > epsilon(0.0))
            apar = apar / spread(kperp2(:, :, ia, :), 4, ntubes)
         elsewhere
            apar = 0.0
         end where
      else if (dist == 'gbar') then
         apar = apar / spread(apar_denom, 4, ntubes)
      else
         if (proc0) write (*, *) 'unknown dist option in get_apar. aborting'
         call mp_abort('unkown dist option in get_apar. aborting')
      end if

   end subroutine get_apar

   !****************************************************************************
   !                                 Advance Apar
   !****************************************************************************
   ! Advance the apar field when electromagnetic effects are included.
   ! This has its own subroutine because Apar is solved using parallel Ampere's Law
   ! which does not involve phi or bpar.
   !****************************************************************************
   subroutine advance_apar(g, dist, apar)

      ! Parallelisation
      use mp, only: mp_abort, sum_allreduce
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      
      ! Parameters
      use parameters_physics, only: include_apar
      use parameters_physics, only: beta
      
      ! Grids
      use grids_species, only: spec
      use grids_z, only: nzgrid
      use grids_velocity, only: nvpa, nmu, vpa
      use calculations_velocity_integrals, only: integrate_vmu
      
      ! Calculations
      use calculations_gyro_averages, only: gyro_average

      implicit none

      ! Arguments
      character(*), intent(in) :: dist
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: apar

      ! Local variables
      integer :: ikxkyz, iky, ikx, iz, it, is
      real :: wgt
      complex :: tmp
      complex, dimension(:, :), allocatable :: scratch
      
      !-------------------------------------------------------------------------
      
      ! Initialise the apar field
      apar = 0.
      
      if (include_apar) then
         allocate (scratch(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            call gyro_average(spread(vpa, 2, nmu) * g(:, :, ikxkyz), ikxkyz, scratch)
            wgt = beta * spec(is)%z * spec(is)%dens_psi0 * spec(is)%stm_psi0
            call integrate_vmu(scratch, iz, tmp)
            apar(iky, ikx, iz, it) = apar(iky, ikx, iz, it) + tmp * wgt
         end do
         ! apar for different species may be spread over processors at this point, so
         ! broadcast to all procs and sum over species
         call sum_allreduce(apar)
         ! divide by the appropriate apar pre-factor to get apar
         call get_apar(apar, dist)
         deallocate (scratch)
      end if

   end subroutine advance_apar

!###############################################################################
!############################ INITALISE & FINALISE #############################
!###############################################################################

   !****************************************************************************
   !*************************** INITALISE THE FIELDS ***************************
   !****************************************************************************
   ! Fill arrays needed for the electromagnetic calculations
   !
   !                                    Definitions
   ! ---------------------------------------------------------------------------
   ! (formally gamtot - computed in field_equations_fluxtube):
   ! denominator_fields = sum_s (Z_s^2 n_s / T_s) int d^3v J0^2 F_0s
   ! 
   ! ---------------------------------------------------------------------------
   !           These are temporary arrays needed to get the full matricies
   ! ---------------------------------------------------------------------------
   ! (formally gamtot13):
   ! denominator_fields13 = - 4 * beta * sum_s (Z_s n_s) int d^3v mu (J0 J1/gamma) F_0s
   !
   ! apar_denom = k_perp^2 + 2 beta sum_s (Z_s^2 n_s / T_s) int d^3v (v_perp^2/2) J0^2 F_0s
   !            = k_perp^2 + sum_s (Z_s^2 n_s / T_s) Gamma_1s
   !
   ! (formally gamtot31):
   ! denominator_fields31 = - 4 * beta * sum_s (Z_s n_s) int d^3v (v_perp^2) (J0 J1/gamma) F_0s
   !
   ! (formally gamtot33):
   ! denominator_fields33 = 1.0 + 8 * beta * sum_s (n_s T_s) int d^3v (mu^2) (J1/gamma)^2 F_0s
   !
   ! ---------------------------------------------------------------------------
   ! These are the elements of the matrix inverse needed to solve for phi and bpar
   ! ---------------------------------------------------------------------------
   ! denominator_fields_inv11 = 1/((denominator_fields13*denominator_fields31)/denominator_fields33)
   ! 
   ! denominator_fields_inv13 = - denominator_fields13 /(denominator_fields*denominator_fields33 - denominator_fields13*denominator_fields31)
   !
   ! denominator_fields_inv33 = denominator_fields/(denominator_fields*denominator_fields33 - denominator_fields13*denominator_fields31)
   !
   ! denominator_fields_inv31 = - denominator_fields31/(denominator_fields*denominator_fields33 - denominator_fields13*denominator_fields31)
   !****************************************************************************
   subroutine init_field_equations_electromagnetic (nfields)

      ! Parallelisation
      use mp, only: sum_allreduce
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, onlY: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      
      ! Arrays
      use arrays, only: kperp2
      use arrays, only: denominator_fields
      use arrays, only: denominator_fields_inv11, denominator_fields_inv13, denominator_fields_inv31, denominator_fields_inv33
      use arrays, only: apar_denom
      
      ! Parameters
      use parameters_physics, only: include_apar, include_bpar
      use grids_kxky, only : nakx, naky
      use parameters_physics, only: beta
      use parameters_physics, only: fphi
      
      ! Grids
      use grids_species, only: spec
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: vpa, mu
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_z, only: nzgrid
      
      ! Calculations
      use arrays_gyro_averages, only: aj0v, aj1v

      implicit none

      integer, intent (inout) :: nfields
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      real :: tmp, wgt, denom_tmp
      real, dimension(:, :), allocatable :: g0

      !-------------------------------------------------------------------------
      
      ! Tally the number of fields to be solved -> phi, apar, bpar
      if (include_apar) nfields = nfields + 1
      if (include_bpar) nfields = nfields + 1

      ! Allocate arrays needed for electromagnetic calculations
      call allocate_field_equations_electromagnetic

      if (.not. (include_apar .or. include_bpar)) return
      
      ! Assume we only have one field line -> not set up for full flux annulus yet
      ia = 1
      
      ! Paralel Ampere's Law apar_denom
      !
      if (include_apar) then
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            ! apar_denom does not depend on flux tube index,
            ! so only compute for one flux tube index
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            g0 = spread(maxwell_vpa(:, is) * vpa**2, 2, nmu) * maxwell_fac(is) &
                 * spread(maxwell_mu(ia, iz, :, is) * aj0v(:, ikxkyz)**2, 1, nvpa)
            wgt = 2.0 * beta * spec(is)%z * spec(is)%z * spec(is)%dens / spec(is)%mass
            call integrate_vmu(g0, iz, tmp)
            apar_denom(iky, ikx, iz) = apar_denom(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(apar_denom)
         apar_denom = apar_denom + kperp2(:, :, ia, :)

         deallocate (g0)
      end if 

      ! 
      if (include_bpar) then
         ! denominator_fields33
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            ! denominator_fields33 does not depend on flux tube index,
            ! so only compute for one flux tube index
            ! denominator_fields33 = 1 + 8 * beta * sum_s (n*T* integrate_vmu(mu*mu*exp(-v^2) *(J1/gamma)*(J1/gamma)))
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            g0 = spread((mu(:) * mu(:) * aj1v(:, ikxkyz) * aj1v(:, ikxkyz)), 1, nvpa) &
                 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
            wgt = 8.0 * spec(is)%temp * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            denominator_fields33(iky, ikx, iz) = denominator_fields33(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(denominator_fields33)

         denominator_fields33 = 1.0 + beta * denominator_fields33

         !denominator_fields13
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            ! denominator_fields13 does not depend on flux tube index,
            ! so only compute for one flux tube index
            ! denominator_fields13 = -4 * sum_s (Z*n* integrate_vmu(mu*exp(-v^2) * J0 *J1/gamma))
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            g0 = spread((mu(:) * aj0v(:, ikxkyz) * aj1v(:, ikxkyz)), 1, nvpa) &
                 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
            wgt = -4.0 * spec(is)%z * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            denominator_fields13(iky, ikx, iz) = denominator_fields13(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(denominator_fields13)
         denominator_fields31 = -0.5 * beta * denominator_fields13 
         deallocate (g0)
      end if

            
      ! Compute coefficients for even part of field solve (phi, bpar)
      if (fphi > epsilon(0.0) .and. include_bpar) then
         do iz = -nzgrid,nzgrid 
            do ikx = 1, nakx
               do iky = 1, naky 
                  ! denominator_fields_inv11
                  denom_tmp = denominator_fields(iky,ikx,iz) - ((denominator_fields13(iky,ikx,iz)*denominator_fields31(iky,ikx,iz))/denominator_fields33(iky,ikx,iz))
                  if (denom_tmp < epsilon(0.0)) then
                     denominator_fields_inv11(iky,ikx,iz) = 0.0
                  else
                     denominator_fields_inv11(iky,ikx,iz) = 1.0/denom_tmp
                  end if
                  ! denominator_fields_inv13, denominator_fields_inv31, denominator_fields_inv33
                  denom_tmp = denominator_fields(iky,ikx,iz)*denominator_fields33(iky,ikx,iz) - denominator_fields13(iky,ikx,iz)*denominator_fields31(iky,ikx,iz)
                  if (denom_tmp < epsilon(0.0)) then
                     denominator_fields_inv13(iky,ikx,iz) = 0.0
                     denominator_fields_inv31(iky,ikx,iz) = 0.0
                     denominator_fields_inv33(iky,ikx,iz) = 0.0
                  else
                     denominator_fields_inv13(iky,ikx,iz) = -denominator_fields13(iky,ikx,iz)/denom_tmp
                     denominator_fields_inv33(iky,ikx,iz) = denominator_fields(iky,ikx,iz)/denom_tmp
                     denominator_fields_inv31(iky,ikx,iz) = -denominator_fields31(iky,ikx,iz)/denom_tmp
                  end if
               end do
            end do
         end do
      end if

      ! Deallocate any arrays that were only needed for the initialisation
      call deallocate_field_equations_electromagnetic_temporary
      
   end subroutine init_field_equations_electromagnetic

   !****************************************************************************
   !***************************** ALLOCATE ARRAYS ******************************
   !****************************************************************************
   ! Allocate arrays needed for solving electromagnetic fields
   ! This includes Apar and Bpar
   ! If include_apar or include_bpar is false, allocate a dummy array of size (1,1,1,1)
   ! This is because these arrays are passed as arguments in many places, so need to 
   ! be allocated, but we don't want to allocate large arrays if they are not needed
   !****************************************************************************
   subroutine allocate_field_equations_electromagnetic

      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, nakx
      use parameters_physics, only: include_apar, include_bpar
      use arrays_fields, only: apar, apar_old
      use arrays_fields, only: bpar, bpar_old
      use arrays, only: denominator_fields_inv11, denominator_fields_inv13, denominator_fields_inv31, denominator_fields_inv33
      use arrays, only: apar_denom

      implicit none

      !----------------------------------------------------------------------
      if (include_apar) then
         if (.not. allocated(apar)) then; allocate (apar(naky, nakx, -nzgrid:nzgrid, ntubes)); apar = 0. ; end if
         if (.not. allocated(apar_old)) then; allocate (apar_old(naky, nakx, -nzgrid:nzgrid, ntubes)); apar_old = 0. ; end if
         if (.not. allocated(apar_denom)) then; allocate (apar_denom(naky, nakx, -nzgrid:nzgrid)); apar_denom = 0. ; end if
      else
         if (.not. allocated(apar)) then; allocate (apar(1, 1, 1, 1)); apar = 0. ; end if
         if (.not. allocated(apar_old)) then; allocate (apar_old(1, 1, 1, 1)); apar_old = 0. ; end if
         if (.not. allocated(apar_denom)) then; allocate (apar_denom(1, 1, 1)); apar_denom = 0. ; end if
      end if

      if (include_bpar) then
         if (.not. allocated(bpar)) then; allocate (bpar(naky, nakx, -nzgrid:nzgrid, ntubes)); bpar = 0. ; end if
         if (.not. allocated(bpar_old)) then; allocate (bpar_old(naky, nakx, -nzgrid:nzgrid, ntubes)); bpar_old = 0. ; end if
         if (.not. allocated(denominator_fields33)) then; allocate (denominator_fields33(naky, nakx, -nzgrid:nzgrid)); denominator_fields33 = 0. ; end if
         if (.not. allocated(denominator_fields13)) then; allocate (denominator_fields13(naky, nakx, -nzgrid:nzgrid)); denominator_fields13 = 0. ; end if
         if (.not. allocated(denominator_fields31)) then; allocate (denominator_fields31(naky, nakx, -nzgrid:nzgrid)); denominator_fields31 = 0. ; end if
         if (.not. allocated(denominator_fields_inv11)) then; allocate (denominator_fields_inv11(naky, nakx, -nzgrid:nzgrid)); denominator_fields_inv11 = 0. ; end if
         if (.not. allocated(denominator_fields_inv31)) then; allocate (denominator_fields_inv31(naky, nakx, -nzgrid:nzgrid)); denominator_fields_inv31 = 0. ; end if
         if (.not. allocated(denominator_fields_inv13)) then; allocate (denominator_fields_inv13(naky, nakx, -nzgrid:nzgrid)); denominator_fields_inv13 = 0. ; end if
         if (.not. allocated(denominator_fields_inv33)) then; allocate (denominator_fields_inv33(naky, nakx, -nzgrid:nzgrid)); denominator_fields_inv33 = 0. ; end if
      else
         if (.not. allocated(bpar)) then; allocate (bpar(1, 1, 1, 1)); bpar = 0. ; end if
         if (.not. allocated(bpar_old)) then; allocate (bpar_old(1, 1, 1, 1)); bpar_old = 0. ; end if
         if (.not. allocated(denominator_fields33)) then; allocate (denominator_fields33(1, 1, 1)); denominator_fields33 = 0. ; end if
         if (.not. allocated(denominator_fields13)) then; allocate (denominator_fields13(1, 1, 1)); denominator_fields13 = 0. ; end if
         if (.not. allocated(denominator_fields31)) then; allocate (denominator_fields31(1, 1, 1)); denominator_fields31 = 0. ; end if
         if (.not. allocated(denominator_fields_inv11)) then; allocate (denominator_fields_inv11(1, 1, 1)); denominator_fields_inv11 = 0. ; end if
         if (.not. allocated(denominator_fields_inv31)) then; allocate (denominator_fields_inv31(1, 1, 1)); denominator_fields_inv31 = 0. ; end if
         if (.not. allocated(denominator_fields_inv13)) then; allocate (denominator_fields_inv13(1, 1, 1)); denominator_fields_inv13 = 0. ; end if
         if (.not. allocated(denominator_fields_inv33)) then; allocate (denominator_fields_inv33(1, 1, 1)); denominator_fields_inv33 = 0. ; end if
      end if

   end subroutine allocate_field_equations_electromagnetic

   !****************************************************************************
   !******************** FINISH THE ELECTROMAGNETIC FIELDS *********************
   !****************************************************************************
   subroutine finish_field_equations_electromagnetic

      use arrays_fields, only: apar
      use arrays_fields, only: apar_old, bpar_old
      use arrays_fields, only: bpar
      use arrays, only: denominator_fields_inv11, denominator_fields_inv13, denominator_fields_inv31, denominator_fields_inv33
      use arrays, only: apar_denom
      
      implicit none

      !----------------------------------------------------------------------

      ! Deallocate electromagnetic arrays
      if (allocated(apar)) deallocate (apar)
      if (allocated(bpar)) deallocate (bpar)
      if (allocated(apar_old)) deallocate(apar_old)
      if (allocated(bpar_old)) deallocate(bpar_old)
      if (allocated(apar_denom)) deallocate (apar_denom)
      if (allocated(denominator_fields_inv11)) deallocate(denominator_fields_inv11)
      if (allocated(denominator_fields_inv31)) deallocate(denominator_fields_inv31)
      if (allocated(denominator_fields_inv13)) deallocate(denominator_fields_inv13)
      if (allocated(denominator_fields_inv33)) deallocate(denominator_fields_inv33)
      
   end subroutine finish_field_equations_electromagnetic

   subroutine deallocate_field_equations_electromagnetic_temporary

      implicit none 

      if (allocated(denominator_fields33)) deallocate (denominator_fields33)
      if (allocated(denominator_fields13)) deallocate (denominator_fields13)
      if (allocated(denominator_fields31)) deallocate (denominator_fields31)

   end subroutine deallocate_field_equations_electromagnetic_temporary

end module field_equations_electromagnetic
