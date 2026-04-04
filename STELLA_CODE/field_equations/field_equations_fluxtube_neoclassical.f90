! ================================================================================================================================================================================ !
! 
! The following routines calculate the corrections needed for solving the field equations in the HO theory. 
!
! ================================================================================================================================================================================ !

module field_equations_fluxtube_neoclassical
   ! Load debug flags.
   ! use debug_flags, only: debug => fields_fluxtube_neoclassical_debug
   
   implicit none

   ! Make routines available to other modules.
   public :: get_phi_neoclassical_correction
   ! public :: get_phi_and_apar_corrections
   ! public :: get_phi_apar_and_bpar_corrections

   private

contains

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
         !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ]                                                                                       ! 
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
         !  When using adiabatic electrons, denominator_fields_neo acquires a factor associated with the adiabatic response:                                        !
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

! =================================================================================================================================================================================== !

   subroutine get_phi_and_apar_neoclassical_corrections
   end subroutine get_phi_and_apar_neoclassical_corrections

! =================================================================================================================================================================================== !

   subroutine get_phi_apar_and_bpar_neoclassical_corrections
   end subroutine get_phi_apar_and_bpar_neoclassical_corrections

! =================================================================================================================================================================================== !

end module field_equations_fluxtube_neoclassical
