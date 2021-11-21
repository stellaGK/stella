module zf_diagnostics

  implicit none

  public :: init_zf_diagnostics, finish_zf_diagnostics
  public :: zf_diag_data, zf_staging, ncalc
  public :: calculate_zf_stress
  public :: clear_zf_staging_array
  public :: zf_prl_str,         &
            zf_prl_str_rad,     &
            zf_prl_str_rad_phi, &
            zf_mirror,          &    
            zf_mirror_rad,      &
            zf_mgn_drft,        &
            zf_mgn_drft_rad,    &
            zf_exb_nl,          &
            zf_exb_nl_rad,      &
            zf_source,          &    
            zf_comm


  private

  integer, parameter :: ncalc  = 12
  integer, parameter :: zf_prl_str         =  1, &
                        zf_prl_str_rad     =  2, &
                        zf_prl_str_rad_phi =  3, &
                        zf_mirror          =  4, & !only depends on endpoints in vpa space
                        zf_mirror_rad      =  5, &
                        zf_mgn_drft        =  6, &
                        zf_mgn_drft_rad    =  7, &
                        zf_exb_nl          =  8, &
                        zf_exb_nl_rad      =  9, &
                        zf_source          = 10, &
                        zf_comm            = 11

  real, dimension(:,:), allocatable :: zf_diag_data

  complex, dimension (:,:,:,:), allocatable :: zf_staging
  ! nx, iz, it, ivmu
  complex, dimension (:,:,:), allocatable :: dphi_ky0
  ! nx, iz, it
  complex, dimension (:), allocatable :: phi_fsa, dphi_fsa
  ! nx

  integer :: nwrite_local

contains

  subroutine init_zf_diagnostics (nwrite)

    use kt_grids, only: nakx
    use zgrid, only: nzgrid, ntubes
    use stella_layouts, only: vmu_lo

    implicit none

    integer, intent (in) :: nwrite

    nwrite_local = nwrite
    if (.not.allocated(zf_diag_data)) then
      allocate (zf_diag_data(ncalc,nakx))
    endif
    zf_diag_data = 0.

    if (.not.allocated(zf_staging)) then
      allocate (zf_staging(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      call clear_zf_staging_array
    endif

    if (.not.allocated(dphi_ky0)) allocate (dphi_ky0(nakx,-nzgrid:nzgrid,ntubes))
    if (.not.allocated(dphi_fsa)) allocate (dphi_fsa(nakx))
    if (.not.allocated(phi_fsa)) allocate (phi_fsa(nakx))

  end subroutine init_zf_diagnostics

  subroutine finish_zf_diagnostics

    implicit none

    if (allocated(zf_staging)) deallocate (zf_staging)
    if (allocated(zf_diag_data)) deallocate (zf_diag_data)
    if (allocated(dphi_ky0)) deallocate (dphi_ky0)
    if (allocated(dphi_fsa)) deallocate (dphi_fsa)
    if (allocated(phi_fsa)) deallocate (phi_fsa)

  end subroutine finish_zf_diagnostics

  subroutine clear_zf_staging_array 

    implicit none
    zf_staging = 0.

  end subroutine clear_zf_staging_array 


  subroutine calculate_zf_stress (zf_index,radial) 

    use kt_grids, only: nakx, rho_d_clamped
    use zgrid, only: nzgrid, ntubes
    use stella_time, only: istep
    use stella_layouts, only: vmu_lo
    use stella_geometry, only: dl_over_b, d_dl_over_b_drho
    use physics_flags, only: radial_variation
    use fields_arrays, only: phi
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

    implicit none

    integer, intent (in) :: zf_index
    logical, optional, intent (in) :: radial

    complex, dimension (:,:), allocatable :: g0k, g0x
    integer :: ivmu, ia, iz, it
    logical radial_

    ia = 1

    radial_ = .false.
    if (present(radial)) radial_ = radial

    if (mod(istep,nwrite_local).eq.0) then
      allocate (g0k(1,nakx))
      allocate (g0x(1,nakx))
      
      if (radial_) then
        do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
          do it = 1, ntubes
            do iz = -nzgrid, nzgrid
              g0k(1,:) = zf_staging(:,iz,it,ivmu)
              call transform_kx2x_unpadded (g0k,g0x)
              g0x(1,:) = g0x(1,:)*rho_d_clamped
              call transform_x2kx_unpadded (g0x,g0k)
              zf_staging(:,iz,it,ivmu) = g0k(1,:)
            enddo
          enddo
        enddo
      endif

      call get_fields_vmulo (zf_staging, dphi_ky0)

      !flux surface average
      phi_fsa = 0.; dphi_fsa = 0.
      do it = 1, ntubes
        do iz = -nzgrid, nzgrid
          if(radial_variation) then
            !phi
            g0k(1,:) = phi(1,:,iz,it)
            call transform_kx2x_unpadded (g0k,g0x)
            g0x(1,:) = g0x(1,:)*(dl_over_b(ia,iz) + d_dl_over_b_drho(ia,iz)*rho_d_clamped)
            call transform_x2kx_unpadded (g0x,g0k)
            phi_fsa = phi_fsa + g0k(1,:)

            !dphi
            g0k(1,:) = dphi_ky0(:,iz,it)
            call transform_kx2x_unpadded (g0k,g0x)
            g0x(1,:) = g0x(1,:)*(dl_over_b(ia,iz) + d_dl_over_b_drho(ia,iz)*rho_d_clamped)
            call transform_x2kx_unpadded (g0x,g0k)
            dphi_fsa = dphi_fsa + g0k(1,:)
          else
            phi_fsa  =  phi_fsa + dl_over_b(ia,iz)*phi(1,:,iz,it)
            dphi_fsa = dphi_fsa + dl_over_b(ia,iz)*dphi_ky0(:,iz,it)
          endif
        enddo
      enddo

      !multiply by conj(phi)
      zf_diag_data(zf_index,:) = real(conjg(phi_fsa(:))*dphi_fsa(:))

    endif



    call clear_zf_staging_array

  end subroutine calculate_zf_stress

  subroutine get_fields_vmulo (g, phi)

    use mp, only: mp_abort, sum_allreduce
    use job_manage, only: time_message
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: imu_idx, is_idx
    use gyro_averages, only: gyro_average, aj0x, aj1x
    use run_parameters, only: fphi, fapar
    use stella_geometry, only: dBdrho, bmag
    use physics_flags, only: radial_variation
    use dist_fn_arrays, only: kperp2, dkperp2dr
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: integrate_species, vperp2
    use kt_grids, only: nakx, naky, rho_d_clamped
    use run_parameters, only: ky_solve_radial
    use species, only: spec
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

    implicit none
    
    complex, dimension (:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,-nzgrid:,:), intent (out) :: phi

    integer :: ivmu, iz, it, ia, imu, is, ikx
    complex, dimension (:,:), allocatable :: gyro_g
    complex, dimension (:,:), allocatable :: g0k, g0x

    ia = 1
    phi = 0.
    if (fphi > epsilon(0.0)) then
       allocate (g0k(1,nakx))
       allocate (g0x(1,nakx))
       allocate (gyro_g(nakx,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
       do it = 1, ntubes
         do iz = -nzgrid, nzgrid
           do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
             is = is_idx(vmu_lo,ivmu)
             imu = imu_idx(vmu_lo,ivmu)
             gyro_g(:,ivmu) = aj0x(1,:,iz,ivmu)*g(:,iz,it,ivmu)
             g0k = 0.0
             if(radial_variation.and.ky_solve_radial.gt.0) then
               g0k(1,:) = gyro_g(:,ivmu) &
                     * (-0.5*aj1x(1,:,iz,ivmu)/aj0x(1,:,iz,ivmu)*(spec(is)%smz)**2 &
                     * (kperp2(1,:,ia,iz)*vperp2(ia,iz,imu)/bmag(ia,iz)**2) &
                     * (dkperp2dr(1,:,ia,iz) - dBdrho(iz)/bmag(ia,iz)) &
                     + dBdrho(iz)/bmag(ia,iz))

               call transform_kx2x_unpadded (g0k,g0x)
               g0x(1,:) = g0x(1,:)*rho_d_clamped
               call transform_x2kx_unpadded (g0x,g0k)
             endif
             gyro_g(:,ivmu) = gyro_g(:,ivmu) + g0k(1,:)

           end do
           do ikx = 1, nakx
             call integrate_species (gyro_g(ikx,:), iz, spec%z*spec%dens_psi0, phi(ikx,iz,it),reduce_in=.false.)
           enddo
         end do
       end do
       deallocate (g0k, g0x, gyro_g)

       call sum_allreduce(phi)

       call get_phi(phi)

    end if
    
  end subroutine get_fields_vmulo

  subroutine get_phi (phi)

    use mp, only: mp_abort, job
#if defined MPI && ISO_C_BINDING
    use mp, only: sgproc0, comm_sgroup
#endif
    use job_manage, only: time_message
    use physics_flags, only: full_flux_surface, radial_variation
    use run_parameters, only: ky_solve_radial, ky_solve_real
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: swap_kxky_ordered, nakx, naky, rho_d_clamped, zonal_mode
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
    use stella_geometry, only: dl_over_b, d_dl_over_b_drho
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg
    use species, only: spec, has_electron_species
    use multibox, only: copy_size
    use fields_arrays, only: gamtot, gamtot3, phi_solve, phizf_solve, phi_ext
    use fields_arrays, only: phi_proj
    use file_utils, only: runtype_option_switch, runtype_multibox
    use fields_arrays, only: exclude_boundary_regions_qn, exp_fac_qn, tcorr_source_qn
#if defined MPI && ISO_C_BINDING
    use fields_arrays, only: qn_window
    use mp_lu_decomposition, only: lu_matrix_multiply_local
#endif
    use linear_solve, only: lu_back_substitution

    implicit none

    complex, dimension (:,-nzgrid:,:), intent (in out) :: phi
    integer :: ia, it, iz, ikx, zmi, zm
    integer :: inmat
    complex, dimension (:,:), allocatable :: g0k, g1k, g0x, g0a
    complex :: tmp
    logical :: has_elec, adia_elec
#if defined MPI && ISO_C_BINDING
    integer :: ierr
#endif

    ia = 1; zm = 0
    has_elec  = has_electron_species(spec)
    adia_elec = .not.has_elec  &
                .and.adiabatic_option_switch.eq.adiabatic_option_fieldlineavg

    if (full_flux_surface) return !not supported yet

    if ((radial_variation.and.ky_solve_radial.gt.0                & 
         .and.runtype_option_switch.ne.runtype_multibox)          &
                               .or.                               &!DSO -> sorry for this if statement
        (radial_variation.and.ky_solve_radial.gt.0.and.job.eq.1   &
         .and.runtype_option_switch.eq.runtype_multibox           &
         .and..not.ky_solve_real)) then
      allocate (g0k(1,nakx))
      allocate (g0x(1,nakx))
      allocate (g0a(1,nakx))

      do it = 1, ntubes
        do iz = -nzgrid, nzgrid
          if(ky_solve_radial.eq.0) then
              phi(:,iz,it) = phi(:,iz,it)/gamtot(1,:,iz)
          elseif (.not.(adia_elec.and.zonal_mode(1))) then
            zmi=zm !zero mode may or may not be included in matrix
            call lu_back_substitution(phi_solve(1,iz)%zloc, &
                                      phi_solve(1,iz)%idx, phi((1+zmi):,iz,it))
            if(zmi.gt.0) phi(zmi,iz,it) = 0.0
          endif
        enddo
      enddo
        
      if(ky_solve_radial.eq.0.and.any(gamtot(1,1,:).lt.epsilon(0.))) &
        phi(1,:,:) = 0.0

      deallocate (g0k,g0x,g0a)
    else if (radial_variation.and.ky_solve_radial.gt.0.and.job.eq.1 &
             .and.runtype_option_switch.eq.runtype_multibox) then
      !call mb_get_phi(phi,has_elec,adia_elec)
    else
      phi = phi/spread(gamtot(1,:,:),3,ntubes)
      if(any(gamtot(1,1,:).lt.epsilon(0.))) phi(1,:,:) = 0.0
    end if

    if(any(gamtot(1,1,:).lt.epsilon(0.))) phi(1,:,:) = 0.0

    if (adia_elec.and.zonal_mode(1)) then
      if(radial_variation.and.ky_solve_radial.gt.0.and.job.eq.1 &
          .and.runtype_option_switch.eq.runtype_multibox.and.ky_solve_real) then
          !this is already taken care of in mb_get_phi
      elseif((radial_variation.and.ky_solve_radial.gt.0               &
              .and.runtype_option_switch.ne.runtype_multibox)         &
                                     .or.                             &
             (radial_variation.and.ky_solve_radial.gt.0.and.job.eq.1  &
              .and.runtype_option_switch.eq.runtype_multibox          &
              .and..not.ky_solve_real))  then
        allocate (g0k(1,nakx))
        allocate (g1k(1,nakx))
        allocate (g0x(1,nakx))

        do it = 1, ntubes
          ! calculate <<g>_psi>_T
          g1k = 0.0
          do iz = -nzgrid, nzgrid-1
            g0k(1,:) = phi(:,iz,it)
            call transform_kx2x_unpadded (g0k,g0x)
            g0x(1,:) = (dl_over_b(ia,iz) + d_dl_over_b_drho(ia,iz)*rho_d_clamped)*g0x(1,:)
            if (exclude_boundary_regions_qn) then
              g0x(1,:) = sum(g0x(1,(copy_size+1):(nakx-copy_size))) &
                          / (nakx - 2*copy_size)
              g0x(1,1:copy_size) = 0.0
              g0x(1,(nakx-copy_size+1):) = 0.0
            else
              g0x(1,:) = sum(g0x(1,:))/nakx
            endif

            call transform_x2kx_unpadded(g0x,g0k)

            g1k = g1k + g0k
          enddo

          if (tcorr_source_qn.lt.epsilon(0.0)) then
            do iz = -nzgrid, nzgrid-1
              phi(:,iz,it) = phi(:,iz,it) - g1k(1,:)
            enddo
          else
            do iz = -nzgrid, nzgrid-1
              phi(:,iz,it) = phi(:,iz,it) &
                               - (1.-exp_fac_qn)*g1k(1,:) - exp_fac_qn*phi_proj(:,1,it)
            enddo
          endif

#if defined MPI && ISO_C_BINDING
          if (sgproc0) then
#endif
            do iz = -nzgrid, nzgrid-1
              do ikx = 1, nakx
                inmat = ikx + nakx*(iz+nzgrid)
                phi_ext(inmat) = phi(ikx,iz,it)
              enddo
            enddo
#if defined MPI && ISO_C_BINDING
          endif
          call mpi_win_fence (0, qn_window, ierr)
#endif
        
#if defined MPI && ISO_C_BINDING
          call lu_matrix_multiply_local (comm_sgroup, 0, qn_window, phizf_solve%zloc, phi_ext)
          call mpi_win_fence (0, qn_window, ierr)
#else
          call lu_back_substitution (phizf_solve%zloc,phizf_solve%idx, phi_ext)
#endif

          do iz = -nzgrid, nzgrid-1
            do ikx = 1, nakx
              inmat = ikx + nakx*(iz+nzgrid)
              phi(ikx,iz,it) = phi_ext(inmat)
            enddo
          enddo

          !enforce periodicity
          phi(:,nzgrid,it) = phi(:,-nzgrid,it)

        enddo
        deallocate(g0k,g1k,g0x)
      else
        if(radial_variation) then
          do ikx = 1, nakx
            do it = 1, ntubes
              tmp = sum(dl_over_b(ia,:)*phi(ikx,:,it))
              phi(ikx,:,it) = phi(ikx,:,it) + tmp*gamtot3(ikx,:)
            end do
          end do
        endif
      end if
    end if
    
  end subroutine get_phi

end module zf_diagnostics
