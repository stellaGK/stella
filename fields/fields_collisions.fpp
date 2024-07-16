module fields_collisions 

    public :: get_fields_by_spec
    public :: get_fields_by_spec_idx

    private

    contains
!###############################################################################
!########################## ADVANCE FIELDS BY SPECIES ##########################
!###############################################################################
!> Note that these advance fields routines are only needed when advancing the  
!> collision operators
!###############################################################################

    !============================================================================
    !========================= ADVANCE FIELDS BY SPEC ===========================
    !============================================================================
    !> This is used in coll_dougherty.f90 
    !============================================================================
    subroutine get_fields_by_spec(g, fld, skip_fsa)

        use mp, only: sum_allreduce
        use stella_layouts, only: kxkyz_lo
        use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
        use gyro_averages, only: gyro_average
        use run_parameters, only: fphi
        use stella_geometry, only: dl_over_b
        use zgrid, only: nzgrid, ntubes
        use vpamu_grids, only: nvpa, nmu
        use vpamu_grids, only: integrate_vmu
        use kt_grids, only: nakx
        use kt_grids, only: zonal_mode
        use species, only: spec, nspec, has_electron_species
        use physics_flags, only: adiabatic_option_switch
        use physics_flags, only: adiabatic_option_fieldlineavg
  
        use fields_arrays, only: gamtot3_h, gamtot_h
  
        implicit none
  
        complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld
        logical, optional, intent(in) :: skip_fsa
  
        real :: wgt
        complex, dimension(:, :), allocatable :: g0
        integer :: ikxkyz, iz, it, ikx, iky, is, ia
        logical :: skip_fsa_local
        complex, dimension(nspec) :: tmp
  
        skip_fsa_local = .false.
        if (present(skip_fsa)) skip_fsa_local = skip_fsa
  
  !      if (debug) write (*, *) 'dist_fn::advance_stella::get_fields_by_spec'
  
        ia = 1
  
        fld = 0.
        if (fphi > epsilon(0.0)) then
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
  
           fld = fld / gamtot_h
  
           if (.not. has_electron_species(spec) .and. (.not. skip_fsa_local) .and. &
               adiabatic_option_switch == adiabatic_option_fieldlineavg) then
              if (zonal_mode(1)) then
                 do ikx = 1, nakx
                    do it = 1, ntubes
                       do is = 1, nspec
                          tmp(is) = sum(dl_over_b(ia, :) * fld(1, ikx, :, it, is))
                          fld(1, ikx, :, it, is) = fld(1, ikx, :, it, is) + tmp(is) * gamtot3_h
                       end do
                    end do
                 end do
              end if
           end if
  
           deallocate (g0)
        end if
  
    end subroutine get_fields_by_spec
  

    !============================================================================
    !========================= ADVANCE FIELDS BY SPEC ===========================
    !============================================================================
    !> This is used in coll_fokkerplanck.f90
    !> Note that is looks identical to the routine above - we don't know why they
    !> are separated 
    !============================================================================
     subroutine get_fields_by_spec_idx(isa, g, fld)
  
        ! apply phi_isa[ ] to all species indices contained in g
        ! ie get phi_isa[g_is1], phi_isa[g_is2], phi_isa[g_is3] ...
  
        use mp, only: sum_allreduce
        use stella_layouts, only: kxkyz_lo
        use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
        use gyro_averages, only: gyro_average
        use run_parameters, only: fphi
        use stella_geometry, only: dl_over_b, bmag
        use zgrid, only: nzgrid, ntubes
        use vpamu_grids, only: vperp2, nvpa, nmu
        use vpamu_grids, only: integrate_vmu
        use kt_grids, only: nakx
        use kt_grids, only: zonal_mode
        use species, only: spec, nspec, has_electron_species
        use physics_flags, only: adiabatic_option_switch
        use physics_flags, only: adiabatic_option_fieldlineavg
        use dist_fn_arrays, only: kperp2
        use spfunc, only: j0

        use fields_arrays, only: gamtot_h, gamtot3_h
        implicit none
  
        complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
        complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld
        integer, intent(in) :: isa
  
        complex, dimension(:, :), allocatable :: g0
        integer :: ikxkyz, iz, it, ikx, iky, is, ia, imu
        complex, dimension(nspec) :: tmp
        real :: wgt
        real :: arg
  
        ia = 1
  
        fld = 0.
        if (fphi > epsilon(0.0)) then
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
           call sum_allreduce(fld)
  
           fld = fld / gamtot_h
  
           if (.not. has_electron_species(spec) .and. &
               adiabatic_option_switch == adiabatic_option_fieldlineavg) then
              if (zonal_mode(1)) then
                 do ikx = 1, nakx
                    do it = 1, ntubes
                       do is = 1, nspec
                          tmp(is) = sum(dl_over_b(ia, :) * fld(1, ikx, :, it, is))
                          fld(1, ikx, :, it, is) = fld(1, ikx, :, it, is) + tmp(is) * gamtot3_h
                       end do
                    end do
                 end do
              end if
           end if
  
           deallocate (g0)
        end if
  
     end subroutine get_fields_by_spec_idx

end module fields_collisions