module fields_em

   public :: allocate_fields_em
   public :: finish_fields_em

   private

   real, dimension(:, :, :), allocatable ::  apar_denom

   !> TODO-GA: add debug flag

contains

!###############################################################################
!###################### ADVANCE ELECTROMAGNETIC FIELDS #########################
!###############################################################################

   !============================================================================
   !=========================== ADVANCE APAR KXKYZLO ===========================
   !============================================================================
   subroutine get_fields_apar_kxkyzlo(g, apar, dist)

      use mp, only: proc0
      use mp, only: sum_allreduce, mp_abort
      use job_manage, only: time_message
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use dist_fn_arrays, only: kperp2
      use gyro_averages, only: gyro_average
      use run_parameters, only: fapar
      use physics_parameters, only: beta
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: vpa
      use vpamu_grids, only: integrate_vmu
      use species, only: spec

      implicit none
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: apar
      character(*), intent(in) :: dist
      complex :: tmp

      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia

      ia = 1

      apar = 0.
      if (fapar > epsilon(0.0)) then
         allocate (g0(nvpa, nmu))
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
         deallocate (g0)
      end if

   end subroutine get_fields_apar_kxkyzlo

   !============================================================================
   !============================ ADVANCE APAR VMULO ============================
   !============================================================================
   subroutine get_fields_apar_vmulo(g, apar, dist)
      use mp, only: mp_abort
      use run_parameters, only: fapar
      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: apar
      character(*), intent(in) :: dist

      apar = 0.
      if (fapar > epsilon(0.0)) then
         ! FLAG -- NEW LAYOUT NOT YET SUPPORTED !!
         call mp_abort('APAR NOT YET SUPPORTED FOR NEW FIELD SOLVE. ABORTING.')
         !        allocate (g0(-nvgrid:nvgrid,nmu))
         !        do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         !           iz = iz_idx(kxkyz_lo,ikxkyz)
         !           ikx = ikx_idx(kxkyz_lo,ikxkyz)
         !           iky = iky_idx(kxkyz_lo,ikxkyz)
         !           is = is_idx(kxkyz_lo,ikxkyz)
         !           g0 = spread(aj0v(:,ikxkyz),1,nvpa)*spread(vpa,2,nmu)*g(:,:,ikxkyz)
         !           wgt = 2.0*beta*spec(is)%z*spec(is)%dens*spec(is)%stm
         !           call integrate_vmu (g0, iz, tmp)
         !           apar(iky,ikx,iz) = apar(iky,ikx,iz) + tmp*wgt
         !        end do
         !        call sum_allreduce (apar)
         !        if (dist == 'h') then
         !           apar = apar/kperp2
         !        else if (dist == 'gbar') then
         !           apar = apar/apar_denom
         !        else if (dist == 'gstar') then
         !           write (*,*) 'APAR NOT SETUP FOR GSTAR YET. aborting.'
         !           call mp_abort('APAR NOT SETUP FOR GSTAR YET. aborting.')
         !        else
         !           if (proc0) write (*,*) 'unknown dist option in get_fields. aborting'
         !           call mp_abort ('unknown dist option in get_fields. aborting')
         !        end if
         !        deallocate (g0)
      end if

   end subroutine get_fields_apar_vmulo


!###############################################################################
!############################ INITALISE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !=========================== INITALISE THE FIELDS ===========================
   !============================================================================ 
   !> Fill arrays needed for the electromagnetic calculations
   !============================================================================ 
   subroutine init_fields_em

      use mp, only: sum_allreduce
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, onlY: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use dist_fn_arrays, only: kperp2
      use gyro_averages, only: aj0v
      use run_parameters, only: fapar
      use physics_parameters, only: beta
      use species, only: spec
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: vpa
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: integrate_vmu
      implicit none

      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      real :: tmp, wgt
      real, dimension(:, :), allocatable :: g0

      if (fapar > epsilon(0.)) then
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

   end subroutine init_fields_em

   !============================================================================
   !============================= ALLOCATE ARRAYS ==============================
   !============================================================================ 
   !> Allocate arrays needed for solving electromagnetic fields 
   !> This includes Apar and Bpar
   !============================================================================ 
   subroutine allocate_fields_em

      use zgrid, only: nzgrid, ntubes
      use stella_layouts, only: vmu_lo
      use kt_grids, only: naky, nakx

      implicit none

      !> TODO-GA: add allocation of apar here
      if (.not. allocated(apar_denom)) then
         allocate (apar_denom(naky, nakx, -nzgrid:nzgrid)); apar_denom = 0.
      end if

   end subroutine allocate_fields_em

   !============================================================================
   !==================== FINISH THE ELECTROMAGNETIC FIELDS =====================
   !============================================================================
   subroutine finish_fields_em

      implicit none

      !>TODO-GA:
      !if (allocated(apar)) deallocate (apar)
      if (allocated(apar_denom)) deallocate (apar_denom)

   end subroutine finish_fields_em
end module fields_em
