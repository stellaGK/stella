! ================================================================================================================================================================================== !
! ------------------------------------- The following routines are responsible for calculating and advancing the field equations in the HO theory ---------------------------------- ! 
! ================================================================================================================================================================================== !

module field_equations_fluxtube_neoclassical   
    implicit none

    ! Make routines available to other modules.
    public :: init_neo_electrostatic_fields
    public :: init_neo_electromagnetic_fields
    public :: finish_neo_electrostatic_fields
    public :: finish_neo_electromagnetic_fields

    public :: advance_fields_fluxtube_using_neo_field_equations

   ! Advance fields for g(kx,ky,z,ivpamus) and g(vpa,mu,ikxkyzs), depending on the layout.
   interface advance_fields_fluxtube_using_neo_field_equations
      module procedure advance_fields_fluxtube_using_neo_field_equations_vmulo
      module procedure advance_fields_fluxtube_using_neo_field_equations_kxkyzlo
   end interface

   private

contains


! ================================================================================================================================================================================ !
! ------------------------------------------------------ Advances the field equations based on the g(kx,ky,z,ivpamus) layout. ---------------------------------------------------- !
! ================================================================================================================================================================================ !

   subroutine advance_fields_fluxtube_using_neo_field_equations_vmulo(g, phi, apar, bpar, dist, skip_fsa)
      ! Parallelisation.
      use mp, only: mp_abort, proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx
      use timers, only: time_field_solve
      
      ! Arrays.
      use arrays_distribution_function, only: phi_gyro
      
      ! Parameters.
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: fphi
      
      ! Grids.
      use grids_z, only: nzgrid
      use grids_species, only: spec
      
      ! Calculations.
      use calculations_gyro_averages, only: gyro_average
      use calculations_velocity_integrals, only: integrate_species
      
      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar, bpar
      character(*), intent(in) :: dist 
      logical, optional, intent(in) :: skip_fsa
      
      ! Local variables.
      logical :: skip_fsa_local
      
      ! =============================================================================== !
      
      ! Used for the Dougherty collision operator.
      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ! Initialise the electrostatic potential phi.
      phi = 0.
      
      if (fphi > epsilon(0.0)) then
          ! Start timer.
          if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         
          ! First gyro-average the distribution function g at each phase space location and store this as phi_gyro = <g>_R = J_0 g in k-space.
          call gyro_average(g, phi_gyro)
        
          ! Calculate phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ] = integrate_species( J_0 * g ).
          call integrate_species(phi_gyro, spec%z * spec%dens_psi0, phi)

          ! Stop timer.
          if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')         
      end if
         
      ! ========================================================================================================================================================= !
      ! ---------- Now that we have the moments of gneo, we need to get the correct fields. This will depend on which fields we are running with. --------------- ! 
      ! ========================================================================================================================================================= !

      ! If we are running electrostatic, we only need to get phi and there is not coupling between fields. Calculate phi by dividing with denominator_fields_neo[iky,ikz,iz].  
      if (fphi > epsilon(0.0) .and. .not. include_apar .and. .not. include_bpar) then     
          call calculate_neo_phi(phi, dist, skip_fsa_local)
      end if  

      ! TO D0 - add routine which couples phi and apar when bpar is not present. 

   end subroutine advance_fields_fluxtube_using_neo_field_equations_vmulo


! ================================================================================================================================================================================ !
! ------------------------------------------------------ Advances the field equations based on the g(kx,ky,z,ikxkyzs) layout. ---------------------------------------------------- !
! ================================================================================================================================================================================ !

   subroutine advance_fields_fluxtube_using_neo_field_equations_kxkyzlo(g, phi, apar, bpar, dist, skip_fsa)
      ! Parallelisation.
      use mp, only: proc0
      use mp, only: sum_allreduce, mp_abort
      use job_manage, only: time_message
      use parallelisation_layouts, only: kxkyz_lo
      use parallelisation_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use timers, only: time_field_solve
      
      ! Parameters.
      use parameters_physics, only: fphi
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: radial_variation
      
      ! Grids.
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: spec
      use grids_z, only: nzgrid
      
      ! Calculations.
      use calculations_gyro_averages, only: gyro_average
      use calculations_velocity_integrals, only: integrate_vmu

      implicit none

      ! Arguments.
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(inout) :: phi, apar, bpar
      logical, optional, intent(in) :: skip_fsa
      character(*), intent(in) :: dist
      
      ! Local variables.
      complex :: tmp
      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      logical :: skip_fsa_local
      
      ! =============================================================================== !
      
      ! Used for the Dougherty collision operator.
      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ! Assume we only have one field line.
      ia = 1
      
      ! Initialise the electrostatic potential phi.
      phi = 0.

      if (fphi > epsilon(0.0)) then
         ! Start timer.
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         
         ! Allocate temporary arrays.
         allocate (g0(nvpa, nmu))
         
         ! Iterate over the (kx,ky,z,mu,vpa,s) points.
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            
            ! First gyro-average the distribution function g at each phase space location and store this as g0 = <g>_R = J_0 g.
            call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
            
            ! Calculate phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ].
            wgt = spec(is)%z * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp           
         end do
         
         ! Deallocate temporary array.
         deallocate (g0)
         
         ! Sum the values on all processors and send them to <proc0>.
         call sum_allreduce(phi)         
      end if

      ! ========================================================================================================================================================= !
      ! ---------- Now that we have the moments of gneo, we need to get the correct fields. This will depend on which fields we are running with. --------------- ! 
      ! ========================================================================================================================================================= !

      ! If we are running electrostatic, we only need to get phi and there is no coupling between fields. Calculate phi by dividing with denominator_fields_neo[iky,ikz,iz].  
      if (fphi > epsilon(0.0) .and. .not. include_apar .and. .not. include_bpar) then
          call calculate_neo_phi(phi, dist, skip_fsa_local)
      end if

      ! TO D0 - add routine which couples phi and apar when bpar is not present.

   end subroutine advance_fields_fluxtube_using_neo_field_equations_kxkyzlo


! ================================================================================================================================================================================ !
! -------------------------------------------------------------------------- Calculate NEO phi. ---------------------------------------------------------------------------------- !
! ================================================================================================================================================================================ !

   subroutine calculate_neo_phi(phi, dist, skip_fsa)
      ! Parallelisation.
      use mp, only: proc0, mp_abort
      use job_manage, only: time_message
      use multibox, only: mb_calculate_phi
      use timers, only: time_field_solve
      
      ! Arrays.
      use arrays, only: denominator_fields
      use arrays, only: denominator_fields_MBR
      use arrays, only: denominator_fields_h
      use arrays, only: denominator_fields_MBR_h
      
      ! Grids.
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: zonal_mode
      use grids_kxky, only: nakx, naky
      use grids_species, only: spec
      
      ! Adiabatic electrons.
      use grids_species, only: has_electron_species
      use grids_species, only: adiabatic_option_switch
      use grids_species, only: adiabatic_option_fieldlineavg
      
      ! Geometry.
      use geometry, only: dl_over_b

      ! HO corrections. 
      use arrays, only: denominator_fields_neo

      implicit none

      ! Arguments.
      character(*), intent(in) :: dist
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
      logical, optional, intent(in) :: skip_fsa
      
      ! Local variables.
      real, dimension(:, :, :, :), allocatable :: denominator_fields_t_neo
      integer :: ia, it, ikx
      complex :: tmp
      logical :: skip_fsa_local
      logical :: has_elec, adia_elec

      ! ============================================================== !
      
      ! Used for the Dougherty collision operator.
      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ! Assume we only have one field line.
      ia = 1
      
      ! Check if we need to add a Boltzmann response for an adiabatic species.
      has_elec = has_electron_species(spec)
      adia_elec = .not. has_elec .and. (adiabatic_option_switch == adiabatic_option_fieldlineavg)

      ! Start timer.
      if (proc0) call time_message(.false., time_field_solve(:, 4), ' calculate_phi')

      ! ================================================================================================================================================= !
      ! --------------------------------------------------- Using g_neo as the distribution function. --------------------------------------------------- !
      ! ================================================================================================================================================= !
      !                                                                                                                                                   !
      ! If we are using the g_neo distribution, then:                                                                                                     ! 
      !                                                                                                                                                   ! 
      !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ]                                                                                ! 
      !     / [ sum_s (Z_s² n_s/T_s) (2B/sqrt(pi)) int dvpa int dmu exp(-v²) (1 - J_0^2) (1 + F_1 - dH_1/dμ|_v∥ 0.5/B_0)]                                 ! 
      !                                                                                                                                                   ! 
      ! denominator_fields_neo[iky,ikz,iz] = [ sum_s (Z_s² n_s/T_s) (2B/sqrt(pi)) int dvpa int dmu exp(-v²) (1 - J_0^2) (1 + F_1 - dH_1/dμ|_v∥ 0.5/B_0)]  ! 
      !                                                                                                                                                   !
      ! To avoid any issues with division by zero we set phi = 0.0 if the denominator is too small. Only thing that makes sense is to set phi = 0.0 if    !
      ! the prefactor for phi in QN is also zero.                                                                                                         !
      !                                                                                                                                                   !
      ! ================================================================================================================================================= !     

      if (dist == 'gneo') then  
          allocate (denominator_fields_t_neo(naky, nakx, -nzgrid:nzgrid, ntubes))
          denominator_fields_t_neo = spread(denominator_fields_neo, 4, ntubes)
          where (denominator_fields_t_neo < epsilon(0.0))
              phi = 0.0
          elsewhere
              phi = phi / denominator_fields_t_neo
          end where
          deallocate (denominator_fields_t_neo)
      else
         ! Abort if <dist> is not recognized.
         if (proc0) write (*, *) 'unknown dist option in calculate_neo_phi. aborting'
         call mp_abort('unknown dist option in calculate_neo_phi. aborting')
         return
      end if

      ! The kx = ky = 0.0 mode is not evolved by stella so make sure this term is set to zero.
      if (any(denominator_fields(1, 1, :) < epsilon(0.))) phi(1, 1, :, :) = 0.0
      if (proc0) call time_message(.false., time_field_solve(:, 4), ' calculate_phi')

      ! TO DO - Add the correct HO MBR response when using the gneo distribution. 
      if (adia_elec .and. zonal_mode(1) .and. .not. skip_fsa_local) then
          if (dist == 'gneo') then
              do ikx = 1, nakx
                  do it = 1, ntubes
                      tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                      phi(1, ikx, :, it) = phi(1, ikx, :, it) + tmp * denominator_fields_MBR(ikx, :)
                  end do
              end do
          else
              if (proc0) write (*, *) 'unknown dist option in calculate_neo_phi. aborting'
              call mp_abort('unknown dist option in calcuate_neo_phi. aborting')
          end if
      end if
   end subroutine calculate_neo_phi


! ================================================================================================================================================================================ !
! --------------------------------------------  Provides the extended denominator for phi when running electrostatic HO simulations ---------------------------------------------- !
! ================================================================================================================================================================================ !

   subroutine init_neo_electrostatic_fields
      ! Parallelisation.
      use mp, only: sum_allreduce
      use parallelisation_layouts, only: kxkyz_lo, iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      
      ! Arrays.
      use arrays, only: denominator_fields_neo, efac
      use arrays_gyro_averages, only: aj0v

      ! Parameters.
      use parameters_physics, only: fphi
      
      ! Grids.
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: spec
      use grids_kxky, only: zonal_mode, akx
      use grids_kxky, only: nakx
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use grids_z, only: nzgrid

      ! Adiabatic electrons.
      use grids_species, only: has_electron_species
      use grids_species, only: ion_species
      use grids_species, only: tite, nine
      use grids_species, only: adiabatic_option_switch
      use grids_species, only: adiabatic_option_fieldlineavg

      ! Calculations.
      use calculations_velocity_integrals, only: integrate_vmu
      
      ! Geometry.
      use geometry, only: dl_over_b, bmag

      ! NEO data. 
      use neoclassical_terms_neo, only: neo_h_global, neo_phi, dneo_h_dmu_global, neo_dens 

      implicit none

      ! Local variables.
      integer :: ikxkyz, iz, it, ikx, iky, ia, is
      real :: tmp, wgt
      real, dimension(:, :), allocatable :: g0

      ! Assume we only have a single field line.
      ia = 1

      ! Allocate denominator_fields_neo. 
      call allocate_neo_electrostatic_fields

      ! Calculate the denominators needed for electrostatic neoclassical simulations.
      if (fphi > epsilon(0.0)) then
      
         ! Allocate temporary arrays.
         allocate (g0(nvpa, nmu))

         ! ======================================================================================================================================================== ! 
         ! When we are evolving NEO's higher order corrections, we use the distribution function g_neo. The corresponding electrostatic QN condition for phi is:    ! 
         !                                                                                                                                                          ! 
         !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g_neo ]                                                                                       ! 
         !     / [ sum_s (Z_s² n_s/T_s) (2B/sqrt(pi)) int dvpa int dmu exp(-v²) (1 - J_0^2) (1 + F_1 - dH_1/dμ|_v∥ 0.5/B_0)]                                        ! 
         !                                                                                                                                                          ! 
         !     denominator_fields_neo[iky,ikz,iz] = [ sum_s (Z_s² n_s/T_s) (2B/sqrt(pi)) int dvpa int dmu exp(-v²) (1 - J_0^2) (1 + F_1 - dH_1/dμ|_v∥ 0.5/B_0)]     ! 
         !                                                                                                                                                          ! 
         ! ======================================================================================================================================================== !
         
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            ! <denominator_fields_neo> does not depend on flux tube index, so only compute for one flux tube index.
            it = it_idx(kxkyz_lo, ikxkyz)
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! Calculate (1 - J_0(a_k)²) exp(v²) for each (kx,ky,z).
            g0 = spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa) * spread(maxwell_vpa(:, is), 2, nmu) &
               * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)

            ! Multiply this by (1 + F_1 - dH_1/dμ|_v∥ 0.5/B_0) for each (kx,ky,z). 
            g0 = g0 * ( 1 + neo_h_global(iz, :, :, is, 1) - spec(is)%z * neo_phi(iz) - ( 0.5/bmag(ia, iz) ) * dneo_h_dmu_global(iz, :, :, is, 1) )
            
            ! Calculate denominator_fields_neo[iky,ikz,iz].
            wgt = spec(is)%z * spec(is)%z * spec(is)%dens_psi0 / spec(is)%temp
            call integrate_vmu(g0, iz, tmp)
            denominator_fields_neo(iky, ikx, iz) = denominator_fields_neo(iky, ikx, iz) + tmp * wgt
         end do
         
         ! Sum the values on all processors and send them to <proc0>.
         call sum_allreduce(denominator_fields_neo)
         
         ! Avoid divide by zero when kx=ky=0; We do not evolve this mode, so the value is irrelevant.
         if (zonal_mode(1) .and. akx(1) < epsilon(0.) .and. has_electron_species(spec)) then
            denominator_fields_neo(1, 1, :) = 0.0
         end if

         ! ======================================================================================================================================================== ! 
         !                                                                                                                                                          ! 
         ! When using adiabatic electrons, denominator_fields_neo acquires a factor associated with the adiabatic response:                                         !
         !                                                                                                                                                          ! 
         ! denominator_fields_neo[iky,ikz,iz] = denominator_fields_neo[iky,ikz,iz] +2B/sqrt(pi) int dvpa int dmu exp(-v²) (F_{1,e} - Z_e * ϕ^1_0)                   !
         !                                                                                                                                                          !
         ! This factor is the lowest order moment of the electron F_1 distribution and is computed in neoclassical_terms_neo.f90 as "neo_dens".                     !
         ! Here is it just added to the denominator.                                                                                                                !
         !                                                                                                                                                          !
         ! ======================================================================================================================================================== !

         if (.not. has_electron_species(spec)) then
            do iz = -nzgrid, nzgrid
                denominator_fields_neo(:, :, iz) = denominator_fields_neo(:, :, iz) + efac * (1 + neo_dens(iz, 2))
            end do 
         end if

         ! Deallocate temporary arrays.
         if (allocated(g0)) deallocate (g0)
      end if
   end subroutine init_neo_electrostatic_fields



! ================================================================================================================================================================================ !
! --------------------------------------------- Provides phi and apar when running electromagnetic HO simulations, in the abscence of bpar. -------------------------------------- !
! ================================================================================================================================================================================ !

    subroutine get_phi_and_apar_neoclassical_corrections(phi, apar, dist)
        ! Parallelisation.
        use mp, only: proc0, mp_abort
        use job_manage, only: time_message
        use timers, only: time_field_solve

        ! Arrays.
        use arrays, only: denominator_fields_neo_inv11, denominator_fields_neo_inv12
        use arrays, only: denominator_fields_neo_inv21, denominator_fields_neo_inv22 

        ! Grids.
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: nakx, naky

        implicit none

        ! Arguments.
        complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar
        character(*), intent(in) :: dist

        ! Local variables.
        integer :: ia, it, ikx, iky, iz
        complex :: antot1, antot2

        ! ======================================================================================================================================================== ! 
        ! When we are including the neoclassical higher order corrections, we use the distribution function g_neo. Due to the presence of F_1, all fluctuating     !
        ! fields now couple to on another in the field equations. In the abscence of bpar, this reduces to a 2 x 2 matrix problem where the solutions provide      !
        ! phi and apar. Symbolically, the fields are given by:                                                                                                     !
        !                                                                                                                                                          !
        ! phi = denominator_fields_neo_inv11*antot1 + denominator_fields_neo_inv12*antot2                                                                          !
        ! apar = denominator_fields_neo_inv21*antot2 + denominator_fields_neo_inv22*antot2                                                                         !                               
        !                                                                                                                                                          ! 
        ! ======================================================================================================================================================== !

        ! Assume we only have one field line
        ia = 1
      
        if (dist == 'gneo') then
            do it = 1, ntubes
                do iz = -nzgrid, nzgrid
                    do ikx = 1, nakx
                        do iky = 1, naky
                            antot1 = phi(iky,ikx,iz,it)
                            antot2 = apar(iky,ikx,iz,it)

                            phi(iky,ikx,iz,it) = denominator_fields_neo_inv11(iky,ikx,iz)*antot1 + denominator_fields_neo_inv12(iky,ikx,iz)*antot2
                            apar(iky,ikx,iz,it) = - denominator_fields_neo_inv21(iky,ikx,iz)*antot1 + denominator_fields_neo_inv22(iky,ikx,iz)*antot2
                        end do
                    end do
                end do
            end do
        else
            if (proc0) write (*, *) 'Unknown dist option in get_fields. Aborting.'
            call mp_abort('Unknown dist option in get_fields. Aborting.')
            return
        end if

    end subroutine get_phi_and_apar_neoclassical_corrections


! =================================================================================================================================================================================== !

    subroutine init_neo_electromagnetic_fields
        ! Parallelisation.
        use mp, only: sum_allreduce
        use parallelisation_layouts, only: kxkyz_lo
        use parallelisation_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx

        ! Arrays.
        use arrays, only: kperp2
        use arrays, only: denominator_fields_neo
        use arrays, only: denominator_fields_neo_inv22, denominator_fields_neo_inv12, denominator_fields_neo_inv11, denominator_fields_neo_inv21

        ! Parameters.
        use parameters_physics, only: include_apar, include_bpar, beta
        use parameters_physics, only: fphi
      
        ! Grids.
        use grids_kxky, only : nakx, naky
        use grids_z, only: nzgrid
        use grids_species, only: spec
        use grids_velocity, only: vpa, mu, nvpa, nmu
        use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
        ! Calculations.
        use calculations_velocity_integrals, only: integrate_vmu
        use arrays_gyro_averages, only: aj0v, aj1v

        implicit none

        ! Local variables.
        integer :: ikxkyz, iz, it, ikx, iky, is, ia
        real :: tmp, wgt, denom_tmp
        real, dimension(:, :), allocatable :: g0
        real, dimension(:, :, :), allocatable :: denominator_fields_neo_11, denominator_fields_neo_12
        real, dimension(:, :, :), allocatable :: denominator_fields_neo_21, denominator_fields_neo_22 
        
        if (.not. (include_apar .or. include_bpar)) return

        ! Allocate the temporary arrays needed to get the inverse matrix elements. 
        allocate (denominator_fields_neo_11(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_11 = 0.0
        allocate (denominator_fields_neo_12(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_12 = 0.0
        allocate (denominator_fields_neo_21(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_21 = 0.0
        allocate (denominator_fields_neo_22(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_22 = 0.0

        ! Allocate the arrays needed for higher order electromagnetic field solve calculations.
        call allocate_neo_electromagnetic_fields

        ia = 1

        ! TO DO - first compute the temporarary arrays: denominator_fields_neo_11, denominator_fields_neo_12, denominator_fields_neo_21 and denominator_fields_neo_22

        ! Now compute the inverse matrix elements. These are the factors that are actually needed in the field solve.
        if (fphi > epsilon(0.0) .and. include_apar) then
            do iz = -nzgrid,nzgrid 
                do ikx = 1, nakx
                    do iky = 1, naky 
                        denom_tmp = denominator_fields_neo_11(iky,ikx,iz)*denominator_fields_neo_22(iky,ikx,iz) &
                        - denominator_fields_neo_12(iky,ikx,iz)*denominator_fields_neo_21(iky,ikx,iz)
                        if (denom_tmp < epsilon(0.0)) then
                            denominator_fields_neo_inv11(iky,ikx,iz) = 0.0
                            denominator_fields_neo_inv12(iky,ikx,iz) = 0.0
                            denominator_fields_neo_inv21(iky,ikx,iz) = 0.0
                            denominator_fields_neo_inv22(iky,ikx,iz) = 0.0
                        else
                            denominator_fields_neo_inv11(iky,ikx,iz) = denominator_fields_neo_22(iky,ikx,iz)/denom_tmp  
                            denominator_fields_neo_inv12(iky,ikx,iz) = -denominator_fields_neo_12(iky,ikx,iz)/denom_tmp
                            denominator_fields_neo_inv21(iky,ikx,iz) = -denominator_fields_neo_21(iky,ikx,iz)/denom_tmp
                            denominator_fields_neo_inv22(iky,ikx,iz) = denominator_fields_neo_11(iky,ikx,iz)/denom_tmp
                        end if
                    end do
                end do
            end do
        end if
       
        ! Deallocate the temporary arrays now that we have the inverse matrix elements.
        deallocate(denominator_fields_neo_11)
        deallocate(denominator_fields_neo_12)
        deallocate(denominator_fields_neo_21)
        deallocate(denominator_fields_neo_22)
    end subroutine init_neo_electromagnetic_fields 


! =================================================================================================================================================================================== !
! --------------------------------------------------- Allocate the arrays needed for higher order electrostatic field solve calculations. ------------------------------------------- !
! =================================================================================================================================================================================== !

    subroutine allocate_neo_electrostatic_fields
        use grids_z, only: nzgrid
        use grids_kxky, only: naky, nakx        

        ! Arrays. 
        use arrays, only: denominator_fields_neo

        implicit none

        if (.not. allocated(denominator_fields_neo)) then; allocate (denominator_fields_neo(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo = 0.; end if 
    end subroutine allocate_neo_electrostatic_fields

! =================================================================================================================================================================================== !
! ------------------------------------------------- Allocate the arrays needed for higher order electromagnetic field solve calculations. ------------------------------------------- !
! =================================================================================================================================================================================== !

    subroutine allocate_neo_electromagnetic_fields
        ! Grids.
        use grids_z, only: nzgrid
        use grids_kxky, only: naky, nakx
      
        ! Parameters. 
        use parameters_physics, only: include_apar, include_bpar

        ! Arrays. 
        use arrays_fields, only: apar, apar_old
        use arrays, only: denominator_fields_neo_inv11, denominator_fields_neo_inv12, denominator_fields_neo_inv13
        use arrays, only: denominator_fields_neo_inv21, denominator_fields_neo_inv22, denominator_fields_neo_inv23
        use arrays, only: denominator_fields_neo_inv31, denominator_fields_neo_inv32, denominator_fields_neo_inv33
  
        implicit none

        if (include_apar) then
            if (.not. allocated(denominator_fields_neo_inv11)) then; allocate (denominator_fields_neo_inv11(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv11 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv12)) then; allocate (denominator_fields_neo_inv12(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv12 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv21)) then; allocate (denominator_fields_neo_inv21(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv21 = 0. ; end if 
            if (.not. allocated(denominator_fields_neo_inv22)) then; allocate (denominator_fields_neo_inv22(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv22 = 0. ; end if 
        end if        
   
        if (include_bpar) then
            if (.not. allocated(denominator_fields_neo_inv13)) then; allocate (denominator_fields_neo_inv13(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv13 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv31)) then; allocate (denominator_fields_neo_inv31(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv31 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv23)) then; allocate (denominator_fields_neo_inv23(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv23 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv32)) then; allocate (denominator_fields_neo_inv32(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv32 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv33)) then; allocate (denominator_fields_neo_inv33(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv33 = 0. ; end if
        end if
    end subroutine allocate_neo_electromagnetic_fields


! =================================================================================================================================================================================== !
! ------------------------------------------------- Deallocate the arrays needed for higher order electrostatic field solve calculations. ------------------------------------------- !
! =================================================================================================================================================================================== !

    subroutine finish_neo_electrostatic_fields
        ! Arrays.
        use arrays, only: denominator_fields_neo
        
        implicit none

        if (allocated(denominator_fields_neo)) deallocate(denominator_fields_neo)
    end subroutine finish_neo_electrostatic_fields


! =================================================================================================================================================================================== !
! ----------------------------------------------- Deallocate the arrays needed for higher order electromagnetic field solve calculations. ------------------------------------------- !
! =================================================================================================================================================================================== !

    subroutine finish_neo_electromagnetic_fields
        ! Arrays.
        use arrays, only: denominator_fields_neo_inv11, denominator_fields_neo_inv12, denominator_fields_neo_inv13
        use arrays, only: denominator_fields_neo_inv21, denominator_fields_neo_inv22, denominator_fields_neo_inv23
        use arrays, only: denominator_fields_neo_inv31, denominator_fields_neo_inv32, denominator_fields_neo_inv33

        implicit none

        if (allocated(denominator_fields_neo_inv11)) deallocate(denominator_fields_neo_inv11)
        if (allocated(denominator_fields_neo_inv12)) deallocate(denominator_fields_neo_inv12)
        if (allocated(denominator_fields_neo_inv13)) deallocate(denominator_fields_neo_inv13)
        if (allocated(denominator_fields_neo_inv21)) deallocate(denominator_fields_neo_inv21)
        if (allocated(denominator_fields_neo_inv22)) deallocate(denominator_fields_neo_inv22)
        if (allocated(denominator_fields_neo_inv23)) deallocate(denominator_fields_neo_inv23)
        if (allocated(denominator_fields_neo_inv31)) deallocate(denominator_fields_neo_inv31)
        if (allocated(denominator_fields_neo_inv32)) deallocate(denominator_fields_neo_inv32)
        if (allocated(denominator_fields_neo_inv33)) deallocate(denominator_fields_neo_inv33)
    end subroutine finish_neo_electromagnetic_fields

! =================================================================================================================================================================================== !

end module field_equations_fluxtube_neoclassical
