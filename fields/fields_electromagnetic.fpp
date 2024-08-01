module fields_em

   use debug_flags, only: debug => fields_electromagnetic_debug

   public :: allocate_fields_em
   public :: finish_fields_em

   public :: get_fields_em_kxkyzlo
   public :: get_fields_em_vmulo

   public :: advance_apar
   
   private

   real, dimension(:, :, :), allocatable ::  apar_denom

   !> TODO-GA: add debug flag

contains

!###############################################################################
!###################### ADVANCE ELECTROMAGNETIC FIELDS #########################
!###############################################################################

   !============================================================================
   !============================ GET FIELDS KXKYZLO ============================
   !============================================================================
   subroutine get_fields_em_kxkyzlo(g, apar, dist)

      use mp, only: proc0
      use mp, only: sum_allreduce, mp_abort
      use job_manage, only: time_message
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use dist_fn_arrays, only: kperp2
      use gyro_averages, only: gyro_average
      use physics_parameters, only: beta
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: vpa
      use vpamu_grids, only: integrate_vmu
      use species, only: spec

      use parameters_physics, only: include_apar, include_bpar

      implicit none
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: apar
      character(*), intent(in) :: dist
      complex :: tmp

      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia

      ia = 1

      if (fphi > epsilon(0.0) .and. include_bpar) then
         if (debug) write (*, *) 'fields_electromagnetic::get_fields_em_kxkyzlo::include_bpar'
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g int_dv_g_vperp2')
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            !> integrate g to get sum_s Z_s n_s J0 g and store in phi
            call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
            wgt = spec(is)%z * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp
            !> integrate g to get - 2 beta sum_s n_s T_s J1 mu g and store in bpar
            call gyro_average_j1(spread(mu, 1, nvpa) * g(:, :, ikxkyz), ikxkyz, g0)
            wgt = -2.0 * beta* spec(is)%z * spec(is)%dens_psi0 * spec(is)%temp_psi0
            call integrate_vmu(g0, iz, tmp)
            bpar(iky, ikx, iz, it) = bpar(iky, ikx, iz, it) + wgt * tmp
         end do
         deallocate (g0)
         call sum_allreduce(phi)
         call sum_allreduce(bpar)
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g int_dv_g_vperp2')

         call get_phi_and_bpar(phi, bpar, dist)
      end if

      apar = 0.
      if (include_apar) then
         if (debug) write (*, *) 'fields_electromagnetic::get_fields_em_kxkyzlo::include_apar'
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

   end subroutine get_fields_em_kxkyzlo

   !============================================================================
   !============================= GET FIELDS VMULO =============================
   !============================================================================
   subroutine get_fields_em_vmulo(g, apar, dist)
      use mp, only: mp_abort
      use run_parameters, only: fapar
      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: apar
      character(*), intent(in) :: dist

      if (fphi > epsilon(0.0) .and. include_bpar) then
         if (debug) write (*, *) 'fields_electromagnetic::get_fields_em_vmulo::include_bpar'
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g int_dv_g_vperp2')
         ! gyroaverage the distribution function g at each phase space location
         call gyro_average(g, g_scratch)
         ! <g> requires modification if radial profile variation is included
         if (radial_variation) call add_radial_correction_int_species(g_scratch)
         ! integrate <g> over velocity space and sum over species
         !> store result in phi, which will be further modified below to account for polarization term
         if (debug) write (*, *) 'dist_fn::advance_stella::get_fields_vmulo::integrate_species_phi'
         call integrate_species(g_scratch, spec%z * spec%dens_psi0, phi)
         ! gyroaverage the distribution function g at each phase space location
         call gyro_average_j1(g, g_scratch)
         ! multiply by mu factor from vperp2
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            imu = imu_idx(vmu_lo, ivmu)
            g_scratch(:, :, :, :, ivmu) = g_scratch(:, :, :, :, ivmu) * mu(imu)
         end do
         ! <g> requires modification if radial profile variation is included
         ! not supported for bpar MRH
         ! integrate <g> over velocity space and sum over species
         !> store result in bpar, which will be further modified below to account for polarization term
         if (debug) write (*, *) 'fields_electromagnetic::get_fields_em_vmulo::integrate_species_bpar'
         call integrate_species(g_scratch, -2.0 * beta * spec%temp_psi0 * spec%dens_psi0, bpar)
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g int_dv_g_vperp2')
         call get_phi_and_bpar(phi, bpar, dist)
      end if

      apar = 0.
      if (include_apar) then
         if (debug) write (*, *) 'fields_electromagnetic::get_fields_em_vmulo::include_apar'
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         ! if fphi > 0, then g_scratch = <g> already calculated above
         !if (fphi < epsilon(0.0)) call gyro_average(g, g_scratch)
         ! MRH remove optimisation for ease of including bpar
         call gyro_average(g, g_scratch)
         ! for parallel Amperes Law, need to calculate parallel current rather than density,
         ! so multiply <g> by vpa before integrating
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            ! get the vpa index
            iv = iv_idx(vmu_lo, ivmu)
            ! multiply by vpa
            g_scratch(:, :, :, :, ivmu) = g_scratch(:, :, :, :, ivmu) * vpa(iv)
         end do
         ! integrate vpa*<g> over velocity space and sum over species
         !> store result in apar, which will be further modified below to account for apar pre-factor
         if (debug) write (*, *) 'fields_electromagnetic::get_fields_em_vmulo::integrate_species_apar'
         call integrate_species(g_scratch, spec%z * spec%dens_psi0 * spec%stm_psi0 * beta, apar)
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         ! divide the apar obtained above by the appropriate Apar pre-factor;
         ! this is just kperp2 if g = <f> is used or apar_denom = (kperp2 + ...)
         ! if gbar = g + <vpa*apar/c> * Ze/T * F_0 is used
         call get_apar(apar, dist)
      end if

   end subroutine get_fields_em_vmulo

   !============================================================================
   !============================ GET APAR AND BPAR= ============================
   !============================================================================
   subroutine get_phi_and_bpar(phi, bpar, dist)

      use mp, only: proc0, mp_abort
      use job_manage, only: time_message
      use zgrid, only: nzgrid, ntubes 
      use arrays_kxky, only: nakx, naky
      use arrays_fields, only: gamtotinv11, gamtotinv13, gamtotinv33, gamtotinv31

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, bpar
      integer :: ia, it, ikx, iky, iz
      complex :: antot1, antot3
      
      character(*), intent(in) :: dist

      if (debug) write (*, *) 'fields_electromagnetic::get_phi_and_bpar'

      ia = 1
      if (proc0) call time_message(.false., time_field_solve(:, 4), ' get_phi_and_bpar')
      if (dist == 'gbar' .or. dist == 'g') then
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  do iky = 1, naky
                     antot1 = phi(iky,ikx,iz,it)
                     antot3 = bpar(iky,ikx,iz,it)
                     phi(iky,ikx,iz,it) = gamtotinv11(iky,ikx,iz)*antot1 + gamtotinv13(iky,ikx,iz)*antot3
                     bpar(iky,ikx,iz,it) = gamtotinv31(iky,ikx,iz)*antot1 + gamtotinv33(iky,ikx,iz)*antot3
                  end do
               end do
            end do
         end do
      else if (dist == 'h') then
         !> divide sum ( Zs int J0 h d^3 v) by sum(Zs^2 ns / Ts) 
         phi = phi / gamtot_h
         !> do nothing for bpar because
         !> bpar = - 2 * beta * sum(Ts ns int (J1/bs) mu h d^3 v)
         !> which is already stored in bpar when dist = 'h'.         
      else
         if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
         call mp_abort('unknown dist option in get_fields. aborting')
         return
      end if

   end subroutine get_phi_and_bpar

   !> Get_apar solves pre-factor * Apar = beta_ref * sum_s Z_s n_s vth_s int d3v vpa * J0 * pdf
   !> for apar, with pdf being either g or gbar (specified by dist input).
   !> the input apar is the RHS of the above equation and is overwritten by the true apar
   !> the pre-factor depends on whether g or gbar is used (kperp2 in former case, with additional
   !> term appearing in latter case)
   subroutine get_apar(apar, dist)

      use mp, only: proc0, mp_abort
      use zgrid, only: nzgrid, ntubes
      use arrays_dist_fn, only: kperp2

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: apar
      character(*), intent(in) :: dist

      integer :: ia

      ! this subroutine only considers flux tubes, so set ia = 1
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

   subroutine advance_apar(g, dist, apar)

      use mp, only: mp_abort, sum_allreduce
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use parameters_physics, only: include_apar
      use parameters_physics, only: beta
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: nvpa, nmu, vpa
      use vpamu_grids, only: integrate_vmu
      use gyro_averages, only: gyro_average

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      character(*), intent(in) :: dist
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: apar

      integer :: ikxkyz, iky, ikx, iz, it, is
      real :: wgt
      complex :: tmp
      complex, dimension(:, :), allocatable :: scratch

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

      use parameters_physics, only: include_apar, include_bpar

      use arrays_fields, only: gamtot, dgamtotdr, gamtot3
      use arrays_fields, only: gamtot13, gamtot31, gamtot33
      use arrays_fields, only: gamtotinv11, gamtotinv13, gamtotinv31, gamtotinv33

      implicit none

      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      real :: tmp, wgt
      real, dimension(:, :), allocatable :: g0

      if (include_apar) nfields = nfields + 1
      if (include_bpar) nfields = nfields + 1

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

      if (include_bpar) then
         ! gamtot33
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            ! gamtot33 does not depend on flux tube index,
            ! so only compute for one flux tube index
            ! gamtot33 = 1 + 8 * beta * sum_s (n*T* integrate_vmu(mu*mu*exp(-v^2) *(J1/gamma)*(J1/gamma)))
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            g0 = spread((mu(:) * mu(:) * aj1v(:, ikxkyz) * aj1v(:, ikxkyz)), 1, nvpa) &
                 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
            wgt = 8.0 * spec(is)%temp * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            gamtot33(iky, ikx, iz) = gamtot33(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(gamtot33)

         gamtot33 = 1.0 + beta * gamtot33

                  !gamtot13
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            ! gamtot13 does not depend on flux tube index,
            ! so only compute for one flux tube index
            ! gamtot13 = -4 * sum_s (Z*n* integrate_vmu(mu*exp(-v^2) * J0 *J1/gamma))
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            g0 = spread((mu(:) * aj0v(:, ikxkyz) * aj1v(:, ikxkyz)), 1, nvpa) &
                 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
            wgt = -4.0 * spec(is)%z * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            gamtot13(iky, ikx, iz) = gamtot13(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(gamtot13)
         gamtot31 = -0.5 * beta * gamtot13 
         deallocate (g0)
      end if

            
      if (fphi > epsilon(0.0) .and. include_bpar) then
         !> compute coefficients for even part of field solve (phi, bpar)
         do iz =-nzgrid,nzgrid 
            do ikx = 1, nakx
               do iky = 1, naky 
                  !> gamtotinv11
                  denom_tmp = gamtot(iky,ikx,iz) - ((gamtot13(iky,ikx,iz)*gamtot31(iky,ikx,iz))/gamtot33(iky,ikx,iz))
                  if (denom_tmp < epsilon(0.0)) then
                     gamtotinv11(iky,ikx,iz) = 0.0
                  else
                     gamtotinv11(iky,ikx,iz) = 1.0/denom_tmp
                  end if
                  !> gamtotinv13, gamtotinv31, gamtotinv33
                  denom_tmp = gamtot(iky,ikx,iz)*gamtot33(iky,ikx,iz) - gamtot13(iky,ikx,iz)*gamtot31(iky,ikx,iz)
                  if (denom_tmp < epsilon(0.0)) then
                     gamtotinv13(iky,ikx,iz) = 0.0
                     gamtotinv31(iky,ikx,iz) = 0.0
                     gamtotinv33(iky,ikx,iz) = 0.0
                  else
                     gamtotinv13(iky,ikx,iz) = -gamtot13(iky,ikx,iz)/denom_tmp
                     gamtotinv33(iky,ikx,iz) = gamtot(iky,ikx,iz)/denom_tmp
                     gamtotinv31(iky,ikx,iz) = -gamtot31(iky,ikx,iz)/denom_tmp
                  end if
               end do
            end do
         end do
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

      use parameters_physics, only: include_apar, include_bpar
      use arrays_fields, only: apar, apar_old
      use arrays_fields, only: bpar, bpar_old

      use arrays_fields, only: gamtot, dgamtotdr, gamtot3
      use arrays_fields, only: gamtot13, gamtot31, gamtot33
      use arrays_fields, only: gamtotinv11, gamtotinv13, gamtotinv31, gamtotinv33

      implicit none

      if (include_apar) then 
         if (.not. allocated(apar_denom)) allocate (apar_denom(naky, nakx, -nzgrid:nzgrid)); apar_denom = 0.
      end if 

      if (include_bpar) then 
         allocate (gamtot33(naky, nakx, -nzgrid:nzgrid)); gamtot33 = 0.
         allocate (gamtot13(naky, nakx, -nzgrid:nzgrid)); gamtot13 = 0.
         allocate (gamtot31(naky, nakx, -nzgrid:nzgrid)); gamtot31 = 0.
         allocate (gamtotinv11(naky, nakx, -nzgrid:nzgrid)); gamtotinv11 = 0.  
         allocate (gamtotinv31(naky, nakx, -nzgrid:nzgrid)); gamtotinv31 = 0.    
         allocate (gamtotinv13(naky, nakx, -nzgrid:nzgrid)); gamtotinv13 = 0.
         allocate (gamtotinv33(naky, nakx, -nzgrid:nzgrid)); gamtotinv33 = 0.
      end if 

   end subroutine allocate_fields_em

   !============================================================================
   !==================== FINISH THE ELECTROMAGNETIC FIELDS =====================
   !============================================================================
   subroutine finish_fields_em

      use arrays_fields, only: gamtot, dgamtotdr, gamtot3
      use arrays_fields, only: gamtot13, gamtot31, gamtot33
      use arrays_fields, only: gamtotinv11, gamtotinv13, gamtotinv31, gamtotinv33

      implicit none

      !>TODO-GA:
      !if (allocated(apar)) deallocate (apar)
      if (allocated(apar_denom)) deallocate (apar_denom)
      if (allocated(gamtot33)) deallocate (gamtot33)
      if (allocated(gamtot13)) deallocate (gamtot13)
      if (allocated(gamtot31)) deallocate (gamtot31)
      if (allocated(gamtotinv11)) deallocate (gamtotinv11)
      if (allocated(gamtotinv31)) deallocate (gamtotinv31)
      if (allocated(gamtotinv13)) deallocate (gamtotinv13)
      if (allocated(gamtotinv33)) deallocate (gamtotinv33)

   end subroutine finish_fields_em

end module fields_em
