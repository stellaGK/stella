!###############################################################################
!                                                                               
!###############################################################################
! 
! Module for ...
! 
!###############################################################################
module quasineutrality_equation_collisions

   ! Load debug flags
   use debug_flags, only: debug => fields_debug

   implicit none

   ! Make routines available to other modules
   public :: get_fields_by_spec
   public :: get_fields_by_spec_idx

   private

contains

!###############################################################################
!########################## ADVANCE FIELDS BY SPECIES ##########################
!###############################################################################
! Note that these advance fields routines are only needed when advancing the
! collision operators
!###############################################################################

   !============================================================================
   !========================= ADVANCE FIELDS BY SPEC ===========================
   !============================================================================
   ! This is used in collisions_dougherty.f90
   !============================================================================
   subroutine get_fields_by_spec(g, fld, skip_fsa)

      use mp, only: sum_allreduce
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use calculations_gyro_averages, only: gyro_average
      use parameters_physics, only: fphi
      use geometry, only: dl_over_b
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: nvpa, nmu
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_kxky, only: nakx
      use grids_kxky, only: zonal_mode
      use grids_species, only: spec, nspec, has_electron_species
      use grids_species, only: adiabatic_option_switch
      use grids_species, only: adiabatic_option_fieldlineavg
      use arrays, only: denominator_QN_MBR_h, denominator_QN_h

      implicit none

      ! Arguments
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld
      logical, optional, intent(in) :: skip_fsa

      ! Local variables
      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      logical :: skip_fsa_local
      complex, dimension(nspec) :: tmp

      !----------------------------------------------------------------------

      if (debug) write (*, *) 'dist_fn::advance_stella::get_fields_by_spec'
      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ! Assume we only have one field line
      ia = 1

      ! Initialise field
      fld = 0.
      
      if (fphi > epsilon(0.0)) then
      
         ! Allocate temporary arrays
         allocate (g0(nvpa, nmu))
         
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            wgt = spec(is)%z * spec(is)%dens_psi0
            call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
            g0 = g0 * wgt
            call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
         end do
         call sum_allreduce(fld)

         fld = fld / denominator_QN_h

         if (.not. has_electron_species(spec) .and. (.not. skip_fsa_local) .and. &
             adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            if (zonal_mode(1)) then
               do ikx = 1, nakx
                  do it = 1, ntubes
                     do is = 1, nspec
                        tmp(is) = sum(dl_over_b(ia, :) * fld(1, ikx, :, it, is))
                        fld(1, ikx, :, it, is) = fld(1, ikx, :, it, is) + tmp(is) * denominator_QN_MBR_h
                     end do
                  end do
               end do
            end if
         end if

         ! Deallocate temporary arrays
         deallocate (g0)
         
      end if

   end subroutine get_fields_by_spec

   !============================================================================
   !========================= ADVANCE FIELDS BY SPEC ===========================
   !============================================================================
   ! This is used in collisions_fokkerplanck.f90
   ! Note that is looks identical to the routine above - we don't know why they
   ! are separated
   ! 
   ! apply phi_isa[ ] to all species indices contained in g
   ! ie get phi_isa[g_is1], phi_isa[g_is2], phi_isa[g_is3] ...
   !============================================================================
   subroutine get_fields_by_spec_idx(isa, g, fld)

      use mp, only: sum_allreduce
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use calculations_gyro_averages, only: gyro_average
      use parameters_physics, only: fphi
      use geometry, only: dl_over_b, bmag
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: vperp2, nvpa, nmu
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_kxky, only: nakx
      use grids_kxky, only: zonal_mode
      use grids_species, only: spec, nspec, has_electron_species
      use grids_species, only: adiabatic_option_switch
      use grids_species, only: adiabatic_option_fieldlineavg
      use arrays, only: kperp2
      use spfunc, only: j0
      use arrays, only: denominator_QN_h, denominator_QN_MBR_h
      
      implicit none

      ! Arguments
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld
      integer, intent(in) :: isa

      ! Local variables
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia, imu
      complex, dimension(nspec) :: tmp
      real :: wgt
      real :: arg
      
      !----------------------------------------------------------------------

      ! Assume we only have one field line
      ia = 1

      ! Initialise field
      fld = 0.
      
      if (fphi > epsilon(0.0)) then
      
         ! Allocate temporary arrays
         allocate (g0(nvpa, nmu))
         
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            wgt = spec(isa)%z * spec(isa)%dens
            do imu = 1, nmu
               ! AVB: changed this for use of j0, check
               arg = spec(isa)%bess_fac * spec(isa)%smz_psi0 * sqrt(vperp2(ia, iz, imu) * kperp2(iky, ikx, ia, iz)) / bmag(ia, iz)
               g0(:, imu) = g(:, imu, ikxkyz) * j0(arg) ! AVB: gyroaverage
            end do
            g0 = g0 * wgt
            call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
         end do
         
         ! Sum the values on all processors and send them to <proc0>
         call sum_allreduce(fld)

         fld = fld / denominator_QN_h

         ! Add a Boltzmann response for adiabatic species
         if (.not. has_electron_species(spec) .and. &
             adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            if (zonal_mode(1)) then
               do ikx = 1, nakx
                  do it = 1, ntubes
                     do is = 1, nspec
                        tmp(is) = sum(dl_over_b(ia, :) * fld(1, ikx, :, it, is))
                        fld(1, ikx, :, it, is) = fld(1, ikx, :, it, is) + tmp(is) * denominator_QN_MBR_h
                     end do
                  end do
               end do
            end if
         end if

         ! Deallocate temporary arrays
         deallocate (g0)
      end if

   end subroutine get_fields_by_spec_idx

end module quasineutrality_equation_collisions
