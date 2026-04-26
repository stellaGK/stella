! ================================================================================================================================================================================ !
! 
! The following routines calculate the corrections needed for solving the field equations in the HO theory. 
!
! ================================================================================================================================================================================ !

module field_equations_fluxtube_neoclassical   
   implicit none

   ! Make routines available to other modules.
   public :: get_phi_neoclassical_correction
   public :: get_phi_and_apar_neoclassical_corrections
   public :: get_phi_apar_and_bpar_neoclassical_corrections
   
   public :: init_field_equations_electromagnetic_neo
   public :: allocate_field_equations_electromagnetic_neo
   public :: finish_field_equations_electromagnetic_neo

   private

contains

! ================================================================================================================================================================================ !
! --------------------------------------------  Provides the extended denominator for phi when running electrostatic HO simulations ---------------------------------------------- !
! ================================================================================================================================================================================ !

   subroutine get_phi_neoclassical_correction
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

         ! ======================================================================================================================================================== ! 
         !                                                                                                                                                          ! 
         ! When using adiabatic electrons, it is usually appropriate to include a field line averaged term ... TO DO !                                              ! 
         !                                                                                                                                                          ! 
         ! ======================================================================================================================================================== !
          
            if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
               ! if (zonal_mode(1)) then
                  ! do ikx = 1, nakx
                     ! tmp = T_e / n_e - int (dl/B)/(sum_(s not e) (Z_s² n_s/T_s) * (1- Gamma0))
                     ! tmp = 1./efac - sum(dl_over_b(ia, :) / denominator_fields(1, ikx, :))
                     ! denominator_fields_MBR = 1/ (T_e/n_e - <1/denominator_field>_FSA )
                     ! denominator_fields_MBR(ikx, :) = 1./(denominator_fields(1, ikx, :) * tmp)
                  ! end do
                  
                  ! Avoid dividing by zero for kx=ky=0 mode, which we do not need anyway.
                  ! if (akx(1) < epsilon(0.)) then
                     ! denominator_fields_MBR(1, :) = 0.0
                  ! end if
               ! end if
            end if
         end if

         ! Deallocate temporary arrays.
         if (allocated(g0)) deallocate (g0)
      end if
   end subroutine get_phi_neoclassical_correction


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
                            apar(iky,ikx,iz,it) = denominator_fields_neo_inv21(iky,ikx,iz)*antot1 + denominator_fields_neo_inv22(iky,ikx,iz)*antot2
                        end do
                    end do
                end do
            end do
        else
            if (proc0) write (*, *) 'Unknown dist option in get_phi_and_apar_neoclassical_corrections. Aborting.'
            call mp_abort('Unknown dist option in get_phi_and_apar_neoclassical_corrections. Aborting.')
            return
        end if

    end subroutine get_phi_and_apar_neoclassical_corrections


! ================================================================================================================================================================================ !
! --------------------------------------------------- Provides phi, apar and bpar when running fully electromagnetic HO simulations. --------------------------------------------- !
! ================================================================================================================================================================================ !

    subroutine get_phi_apar_and_bpar_neoclassical_corrections(phi, apar, bpar, dist)
        use grids_z, only: nzgrid  

        implicit none

        ! Arguments.
        complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
        character(*), intent(in) :: dist

    end subroutine get_phi_apar_and_bpar_neoclassical_corrections

! =================================================================================================================================================================================== !

    subroutine init_field_equations_electromagnetic_neo
        ! Parallelisation.
        use mp, only: sum_allreduce
        use parallelisation_layouts, only: kxkyz_lo
        use parallelisation_layouts, onlY: iz_idx, it_idx, ikx_idx, iky_idx, is_idx

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
      
        ! Geometry.
        use geometry, only: bmag
    
        ! Calculations.
        use calculations_velocity_integrals, only: integrate_vmu
        use arrays_gyro_averages, only: aj0v, aj1v

        ! Neoclassical data.
        use neoclassical_terms_neo, only: neo_h_global, neo_phi, dneo_h_dmu_global, dneo_h_dvpa_global

        implicit none

        ! Local variables.
        integer :: ikxkyz, iz, it, ikx, iky, is, ia
        real :: tmp, wgt, denom_tmp
        real, dimension(:, :), allocatable :: g0
        real, dimension(:, :, :), allocatable :: denominator_fields_neo_12
        real, dimension(:, :, :), allocatable :: denominator_fields_neo_21, denominator_fields_neo_22 

        ! Allocate the arrays needed for higher order electromagnetic field solve calculations.
        call allocate_field_equations_electromagnetic_neo

        if (.not. (include_apar .or. include_bpar)) return

        ! Allocate the temporary arrays needed to get the inverse matrix elements. 
        allocate (denominator_fields_neo_12(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_12 = 0.0
        allocate (denominator_fields_neo_21(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_21 = 0.0
        allocate (denominator_fields_neo_22(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_22 = 0.0

        ia = 1

        ! First compute the temporarary arrays: denominator_fields_neo_11, denominator_fields_neo_12, denominator_fields_neo_21 and denominator_fields_neo_22.
        ! Note that denominator_fields_neo_11 has already been calculated in the previous subroutine as denominator_fields_neo. 

        ! ======================================================================================================================================================== ! 
        ! Begin by calculating denominator_fields_neo_12 = DOCUMENTATION: TO DO  
        !
        ! ======================================================================================================================================================== !

        ! Begin with denominator_fields_neo_12.
        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! Calculate exp(v²) * vpa(iv) for each (kx,ky,z).
            g0 = spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is) * spread(vpa(:), 2, nmu)

            ! Multiply this by (1 - J_0(a_k)²) * (dH_1/dμ|_v∥ / 2B_0 - F_1) + (F_1 - dH_1/dv∥|_μ / 2v∥) for each (kx,ky,z). 
            g0 = g0 * ( spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa) * ( 0.5 * dneo_h_dmu_global(iz, :, :, is, 1) / bmag(ia, iz) &
            - neo_h_global(iz, :, :, is, 1) + spec(is)%z * neo_phi(iz) ) & 
            + neo_h_global(iz, :, :, is, 1) - spec(is)%z * neo_phi(iz) - 0.5 * dneo_h_dmu_global(iz, :, :, is, 1) / spread(vpa(:), 2, nmu) )

            ! Calculate denominator_fields_neo_12[iky,ikz,iz].
            wgt = 2.0 * spec(is)%z * spec(is)%z * spec(is)%dens_psi0 * spec(is)%stm / spec(is)%temp
            call integrate_vmu(g0, iz, tmp)
            denominator_fields_neo_12(iky, ikx, iz) = denominator_fields_neo_12(iky, ikx, iz) + tmp * wgt
        end do

        ! Sum the values on all processors and send them to <proc0>.
        call sum_allreduce(denominator_fields_neo_12)

        ! ======================================================================================================================================================== ! 
        ! Now we move on to denominator_fields_neo_21 = DOCUMENTATION: TO DO
        !
        ! ======================================================================================================================================================== !
 
        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! Calculate (1 - J_0(a_k)²) * exp(v²) * vpa(iv) for each (kx,ky,z).
            g0 = spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa) * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is) &
            * spread(vpa(:), 2, nmu)

            ! Multiply this by (F_1 - dH_1/dμ|_v∥ / 2B_0) for each (kx,ky,z). 
            g0 = g0 * ( neo_h_global(iz, :, :, is, 1) - spec(is)%z * neo_phi(iz) - 0.5 * dneo_h_dmu_global(iz, :, :, is, 1) / bmag(ia, iz) )

            ! Calculate denominator_fields_neo_21[iky,ikz,iz].
            wgt = beta * spec(is)%z * spec(is)%z * spec(is)%dens_psi0 * spec(is)%stm / spec(is)%temp
            call integrate_vmu(g0, iz, tmp)
            denominator_fields_neo_21(iky, ikx, iz) = denominator_fields_neo_21(iky, ikx, iz) + tmp * wgt
        end do

        ! Sum the values on all processors and send them to <proc0>.
        call sum_allreduce(denominator_fields_neo_21)

        ! ======================================================================================================================================================== ! 
        ! Finally calculate denominator_fields_neo_22 = DOCUMENTATION: TO DO 
        !
        ! ======================================================================================================================================================== !

        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! Calculate exp(v²) * vpa(iv)² for each (kx,ky,z).
            g0 = spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is) * spread(vpa(:), 2, nmu) * spread(vpa(:), 2, nmu)

            ! Multiply this by (1 - J_0(a_k)²) * (dH_1/dμ|_v∥ / 2B_0 - F_1) + (F_1 - dH_1/dv∥|_μ / 2v∥) for each (kx,ky,z). 
            g0 = g0 * ( spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa) * ( 0.5 * dneo_h_dmu_global(iz, :, :, is, 1) / bmag(ia, iz) &
            - neo_h_global(iz, :, :, is, 1) + spec(is)%z * neo_phi(iz) ) &
            + neo_h_global(iz, :, :, is, 1) - spec(is)%z * neo_phi(iz) - 0.5 * dneo_h_dmu_global(iz, :, :, is, 1) / spread(vpa(:), 2, nmu) )

            ! Calculate denominator_fields_neo_22[iky,ikz,iz].
            wgt = 2.0 * beta * spec(is)%z * spec(is)%z * spec(is)%dens_psi0 * spec(is)%stm * spec(is)%stm / spec(is)%temp
            call integrate_vmu(g0, iz, tmp)
            denominator_fields_neo_22(iky, ikx, iz) = denominator_fields_neo_22(iky, ikx, iz) + tmp * wgt

            ! Add the kperp2 factor.
            denominator_fields_neo_22(iky, ikx, iz) = denominator_fields_neo_22(iky, ikx, iz) + kperp2(iky, ikx, ia, iz)
        end do

        ! Sum the values on all processors and send them to <proc0>.
        call sum_allreduce(denominator_fields_neo_22)

        ! Now compute the inverse matrix elements. These are the factors that are actually needed in the field solve.
        if (fphi > epsilon(0.0) .and. include_apar) then
            do iz = -nzgrid,nzgrid 
                do ikx = 1, nakx
                    do iky = 1, naky 
                        denom_tmp = denominator_fields_neo(iky,ikx,iz)*denominator_fields_neo_22(iky,ikx,iz) &
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
                            denominator_fields_neo_inv22(iky,ikx,iz) = denominator_fields_neo(iky,ikx,iz)/denom_tmp
                        end if
                    end do
                end do
            end do
        end if
       
        ! Deallocate the temporary arrays now that we have the inverse matrix elements.
        deallocate(denominator_fields_neo_12)
        deallocate(denominator_fields_neo_21)
        deallocate(denominator_fields_neo_22)        
    end subroutine init_field_equations_electromagnetic_neo


! =================================================================================================================================================================================== !
! ------------------------------------------------- Allocate the arrays needed for higher order electromagnetic field solve calculations. ------------------------------------------- !
! =================================================================================================================================================================================== !

    subroutine allocate_field_equations_electromagnetic_neo
        ! Grids.
        use grids_z, only: nzgrid, ntubes
        use grids_kxky, only: naky, nakx
      
        ! Parameters. 
        use parameters_physics, only: include_apar

        ! Arrays. 
        use arrays_fields, only: apar, apar_old
        use arrays, only: apar_denom
        use arrays, only: denominator_fields_neo_inv11, denominator_fields_neo_inv12
        use arrays, only: denominator_fields_neo_inv21, denominator_fields_neo_inv22

        implicit none

        ! If include_apar or include_bpar is false, allocate a dummy array of size (1,1,1,1).
        ! This is because these arrays are passed as arguments in many places, so need to be allocated, but we don't want to allocate large arrays if they are not needed.

        if (include_apar) then
            if (.not. allocated(apar)) then; allocate (apar(naky, nakx, -nzgrid:nzgrid, ntubes)); apar = 0. ; end if
            if (.not. allocated(apar_old)) then; allocate (apar_old(naky, nakx, -nzgrid:nzgrid, ntubes)); apar_old = 0. ; end if
            if (.not. allocated(apar_denom)) then; allocate (apar_denom(naky, nakx, -nzgrid:nzgrid)); apar_denom = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv11)) then; allocate (denominator_fields_neo_inv11(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv11 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv12)) then; allocate (denominator_fields_neo_inv12(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv12 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv21)) then; allocate (denominator_fields_neo_inv21(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv21 = 0. ; end if 
            if (.not. allocated(denominator_fields_neo_inv22)) then; allocate (denominator_fields_neo_inv22(naky, nakx, -nzgrid:nzgrid)); denominator_fields_neo_inv22 = 0. ; end if 
            ! Don't think there is any physical motivation to run simulations with finite beta and bpar but no apar included.
            ! Therefore do not need to consider this scenario in setting up the field solver? 
        else
            if (.not. allocated(apar)) then; allocate (apar(1, 1, 1, 1)); apar = 0. ; end if
            if (.not. allocated(apar_old)) then; allocate (apar_old(1, 1, 1, 1)); apar_old = 0. ; end if
            if (.not. allocated(apar_denom)) then; allocate (apar_denom(1, 1, 1)); apar_denom = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv11)) then; allocate (denominator_fields_neo_inv11(1, 1, 1)); denominator_fields_neo_inv11 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv12)) then; allocate (denominator_fields_neo_inv12(1, 1, 1)); denominator_fields_neo_inv12 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv21)) then; allocate (denominator_fields_neo_inv21(1, 1, 1)); denominator_fields_neo_inv21 = 0. ; end if
            if (.not. allocated(denominator_fields_neo_inv22)) then; allocate (denominator_fields_neo_inv22(1, 1, 1)); denominator_fields_neo_inv22 = 0. ; end if
        end if
              
   
    end subroutine allocate_field_equations_electromagnetic_neo


! =================================================================================================================================================================================== !
! ----------------------------------------------- Deallocate the arrays needed for higher order electromagnetic field solve calculations. ------------------------------------------- !
! =================================================================================================================================================================================== !

    subroutine finish_field_equations_electromagnetic_neo
        implicit none
    end subroutine finish_field_equations_electromagnetic_neo

! =================================================================================================================================================================================== !

end module field_equations_fluxtube_neoclassical
