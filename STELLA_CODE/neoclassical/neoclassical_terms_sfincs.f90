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
   use arrays_neoclassical, only: dneo_h_dz, dneo_h_dvpa, dneo_h_drho
   use arrays_neoclassical, only: dneo_phi_dz, dneo_phi_drho, dneo_phi_dalpha
   use arrays_neoclassical, only: neo_alpha_fac

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
         sfincs_is_enabled = include_neoclassical_terms .and. neo_option_switch == 1
      end function sfincs_is_enabled

! ================================================================================================================================================================================= !
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- !
! ================================================================================================================================================================================= !

   subroutine init_neoclassical_terms_sfincs

      use grids_z, only: nzgrid
      use grids_kxky, only: nalpha
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: nspec
      use parallelisation_layouts, only: vmu_lo
      use sfincs_interface, only: get_neo_from_sfincs

      use namelist_neoclassical_input, only: read_namelist_neoclassical_input

      ! Neoclassical arrays (allocation is a no-op if already done).
      use arrays_neoclassical, only: allocate_neoclassical_arrays

      implicit none

      real, dimension(:, :, :, :, :, :), allocatable :: f_neoclassical
      real, dimension(:, :, :), allocatable :: phi_neoclassical
      real, dimension(:, :, :, :, :), allocatable :: dfneo_dalpha_local

      integer :: iz, ialpha

      ! Only initialise once.
      if (initialised_neoclassical_terms) return
      initialised_neoclassical_terms = .true.

      ! Read the input parameters in the input file.
      call read_namelist_neoclassical_input(include_neoclassical_terms, neo_option_switch, nradii, drho)

      ! Broadcast the input parameters to all processors.
      call broacast_arrays

      if (include_neoclassical_terms) then
         ! The neoclassical arrays used in the GKE are owned by arrays_neoclassical. They should already be allocated
         ! during initialisation; this call only guarantees it.
         call allocate_neoclassical_arrays

         allocate (f_neoclassical(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec, -nradii / 2:nradii / 2))
         allocate (phi_neoclassical(nalpha, -nzgrid:nzgrid, -nradii / 2:nradii / 2))
         allocate (dfneo_dalpha_local(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec))
         select case (neo_option_switch)
         case (neo_option_sfincs)
            call get_neo_from_sfincs(nradii, drho, f_neoclassical, phi_neoclassical, dfneo_dalpha_local, dneo_phi_dalpha)
         end select
         if (debug) write (6, *) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_dzed'
         call get_dfneo_dzed(f_neoclassical(:, :, :, :, :, 0), dneo_h_dz)
         if (debug) write (6, *) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_dvpa'
         call get_dfneo_dvpa(f_neoclassical(:, :, :, :, :, 0), dneo_h_dvpa)
         if (debug) write (6, *) 'neoclassical_terms::init_neoclassical_terms::get_dfneo_drho'
         call get_dfneo_drho(f_neoclassical, dneo_h_drho)
         if (debug) write (6, *) 'neoclassical_terms::init_neoclassical_terms::get_dphineo_dzed'
         call get_dphineo_dzed(phi_neoclassical(:, :, 0), dneo_phi_dz)
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
         deallocate (f_neoclassical, phi_neoclassical, dfneo_dalpha_local)
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

   end subroutine init_neoclassical_terms_sfincs

   !****************************************************************************
   !                                      Title
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
   !                                      Title
   !****************************************************************************
   subroutine get_dfneo_dvpa(fneo, dfneo)

      use calculations_finite_differences, only: fd5pt
      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: dvpa
      use grids_species, only: nspec
      use grids_kxky, only: nalpha

      implicit none

      real, dimension(:, -nzgrid:, :, :, :), intent(in) :: fneo
      real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(out) :: dfneo

      integer :: ia, iz, imu, is
      real, dimension(:), allocatable :: tmp1, tmp2
      real, dimension(:, :, :, :, :), allocatable :: dfneo_local

      allocate (tmp1(nvpa), tmp2(nvpa))
      allocate (dfneo_local(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec))

      do is = 1, nspec
         do imu = 1, nmu
            do iz = -nzgrid, nzgrid
               do ia = 1, nalpha
                  ! hack to avoid dealing with negative indices in fd5pt
                  tmp1 = fneo(ia, iz, :, imu, is)
                  call fd5pt(tmp1, tmp2, dvpa)
                  dfneo_local(ia, iz, :, imu, is) = tmp2
               end do
            end do
         end do
      end do

      do iz = -nzgrid, nzgrid
         do ia = 1, nalpha
            call distribute_vmus_over_procs(dfneo_local(ia, iz, :, :, :), dfneo(ia, iz, :))
         end do
      end do

      deallocate (dfneo_local)
      deallocate (tmp1, tmp2)

   end subroutine get_dfneo_dvpa

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine get_dfneo_dzed(fneo, dfneo)

      use calculations_finite_differences, only: fd5pt
      use grids_z, only: nztot, nzgrid, delzed
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: nspec
      use parallelisation_layouts, only: vmu_lo
      use grids_kxky, only: nalpha

      implicit none

      real, dimension(:, -nzgrid:, :, :, :), intent(in) :: fneo
      real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(out) :: dfneo

      integer :: iv, imu, is, iz, ia
      real, dimension(:), allocatable :: tmp1, tmp2
      real, dimension(:), allocatable :: dfneo_local(:, :, :, :, :)

      allocate (tmp1(nztot), tmp2(nztot))
      allocate (dfneo_local(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec))

      do is = 1, nspec
         do imu = 1, nmu
            do iv = 1, nvpa
               do ia = 1, nalpha
                  ! hack to avoid dealing with negative indices in fd5pt
                  tmp1 = fneo(ia, :, iv, imu, is)
                  call fd5pt(tmp1, tmp2, delzed(0))
                  dfneo_local(ia, :, iv, imu, is) = tmp2
               end do
            end do
         end do
      end do

      do iz = -nzgrid, nzgrid
         do ia = 1, nalpha
            call distribute_vmus_over_procs(dfneo_local(ia, iz, :, :, :), dfneo(ia, iz, :))
         end do
      end do

      deallocate (dfneo_local)
      deallocate (tmp1, tmp2)

   end subroutine get_dfneo_dzed

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine get_dfneo_drho(fneo, dfneo)

      use calculations_finite_differences, only: fd3pt, fd5pt
      use grids_z, only: nzgrid
      use grids_velocity, only: nvpa, nmu
      use grids_species, only: nspec
      use parallelisation_layouts, only: vmu_lo
      use grids_kxky, only: nalpha

      implicit none

      real, dimension(:, -nzgrid:, :, :, :, -nradii/2:), intent(in) :: fneo
      real, dimension(:, -nzgrid:, vmu_lo%llim_proc:), intent(out) :: dfneo

      integer :: ia, iz, iv, imu, is
      real, dimension(:), allocatable :: tmp1, tmp2
      real, dimension(:, :, :, :, :), allocatable :: dfneo_local

      allocate (tmp1(nradii), tmp2(nradii))
      allocate (dfneo_local(nalpha, -nzgrid:nzgrid, nvpa, nmu, nspec))

      do is = 1, nspec
         do imu = 1, nmu
            do iv = 1, nvpa
               do iz = -nzgrid, nzgrid
                  do ia = 1, nalpha
                     ! hack to avoid dealing with negative indices in fd5pt
                     tmp1 = fneo(ia, iz, iv, imu, is, :)
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

      deallocate (dfneo_local)
      deallocate (tmp1, tmp2)

   end subroutine get_dfneo_drho

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine get_dphineo_dzed(phineo, dphineo)

      use calculations_finite_differences, only: fd5pt
      use grids_z, only: nztot, nzgrid, delzed
      use grids_kxky, only: nalpha

      implicit none

      real, dimension(:, -nzgrid:), intent(in) :: phineo
      real, dimension(:, -nzgrid:), intent(out) :: dphineo

      integer :: ia
      real, dimension(:), allocatable :: tmp1, tmp2

      !-------------------------------------------------------------------------

      allocate (tmp1(nztot), tmp2(nztot))

      do ia = 1, nalpha
         ! hack to avoid dealing with negative indices in fd5pt
         tmp1 = phineo(ia, :)
         call fd5pt(tmp1, tmp2, delzed(0))
         dphineo(ia, :) = tmp2
      end do

      deallocate (tmp1, tmp2)

   end subroutine get_dphineo_dzed

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