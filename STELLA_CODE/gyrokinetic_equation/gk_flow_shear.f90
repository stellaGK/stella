!###############################################################################
!###############################################################################
!###############################################################################
! 
! The local version of stella currently implements equilibrium ﬂow shear using 
! the discrete wavenumber-shift method formulated by Hammett et al.
!     - See [2006 - Hammett]
!     - See section 3.5. Equilibrium ﬂow shear in [2022 - St-Onge]
! 
! Default parameters related to the shear flow
! 
!   flow_shear
!     prp_shear_enabled = .false.
!     hammett_flow_shear = .true.
!     g_exb = 0.0
!     g_exbfac = 1.0
!     omprimfac = 1.0
!   /
! 
! Note that while be default <hammett_flow_shear> = True, the advance_parallel_
! flow_shear() routine will only be called if |<g_exb>| > 0. Therefore, parallel 
! flow shear is not included in the simulations by default.
! 
! This module also implements perpendicular and parallel flow shear, related to 
! a mean toroidal flow R*Omega_zeta, implemented for radial variation.
! 
! For notes on the mathematics see [2022 - St-Onge - A novel approach to radially 
! global gyrokinetic simulation using the flux-tube code stella]. E.g., equation 28c:
!     omega_{zeta,k,s} = -k_y * sqrt(m_s/T_s) * q*a/r * I/B * exp(-v²) * gamma_E
! 
! The shear rate is defined as
!     gamma_E = (r/q) (dOmega_zeta(dr) (a/v_{th,ref})
! 
! The advance_parallel_flow_shear() routine add the following term to the gyrokinetic
! equation, as can be seen the 6th term in equation (27) in [2022 - St-Onge]
!     - i * omega_{zeta,k,s} * J0 * phi * code_dt
! 
! To calculate this term we have calculated:
!     <prl_shear> = - sqrt(m_s/T_s) * q*a/r * I/B * exp(-v²) * gamma_E * code_dt (?)
!                 = omega_{zeta,k,s} * code_dt / k_y
! 
!###############################################################################
module gk_flow_shear

   implicit none

   ! Make routines available to other modules
   public :: initialised_flow_shear
   public :: init_flow_shear, finish_flow_shear
   public :: read_parameters_flow_shear
   public :: prl_shear, prl_shear_p, prp_shear
   public :: advance_parallel_flow_shear, advance_perp_flow_shear
   public :: v_edge, v_shift
   public :: shift_times

   ! Input parameters
   public :: prp_shear_enabled               ! Enables perpendicular flow shear
   public :: hammett_flow_shear
   public :: g_exb, g_exbfac, omprimfac 

   private

   complex, dimension(:, :), allocatable :: upwind_advect
   real, dimension(:), allocatable :: prp_shear, shift_times
   integer :: shift_sign, shift_start
   real :: v_edge, v_shift = 0.
   
   ! Input parameters
   logical :: prp_shear_enabled
   logical :: hammett_flow_shear 
   real :: g_exb, g_exbfac, omprimfac
   
   ! Parallel flow shear
   real, dimension(:, :, :), allocatable :: prl_shear, prl_shear_p

   ! Only initialise once
   logical :: initialised_flow_shear = .false.

contains


   !****************************************************************************
   !                              Read input file                               
   !****************************************************************************
   subroutine read_parameters_flow_shear

      ! Parallelisation
      use mp, only: broadcast
      
      ! Read namelists from input file
      use namelist_flow_shear, only: read_namelist_flow_shear

      implicit none

      !-------------------------------------------------------------------------
      
      ! Read the "flow_shear" namelist in the input file
      call read_namelist_flow_shear(prp_shear_enabled, hammett_flow_shear, g_exb, g_exbfac, omprimfac)
      
      ! Broadcast the input parameters
      call broadcast(prp_shear_enabled)
      call broadcast(hammett_flow_shear)
      call broadcast(g_exb)
      call broadcast(g_exbfac)
      call broadcast(omprimfac)

   end subroutine read_parameters_flow_shear

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine init_flow_shear
   
      ! Parallelisation
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use mp, only: job, send, receive, crossdomprocs, subprocs, scope
      use file_utils, only: runtype_option_switch, runtype_multibox
      use job_manage, only: njobs
      
      ! Grids
      use stella_time, only: code_dt
      use grids_species, only: spec
      use grids_z, only: nzgrid
      use grids_kxky, only: x, x_d, akx, aky, zonal_mode, box
      use grids_kxky, only: nalpha, nx, nakx, naky, ikx_max
      use grids_velocity, only: vperp2, vpa, mu
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      
      ! Geometry
      use geometry, only: q_as_x, geo_surf, bmag, btor, rmajor, dBdrho, dIdrho
      use geometry, only: dydalpha, drhodpsi
      
      ! Calculations and flags
      use constants, only: zi, pi
      use parameters_numerical, only: maxwellian_normalization
      use parameters_physics, only: radial_variation
      
      ! Flow shear arrays
      use arrays_store_useful, only: shift_state

      implicit none

      ! Local variables
      integer :: is, imu, iv, ivmu, iz, ia
      real, dimension(:, :), allocatable :: energy

      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_flow_shear) return
      initialised_flow_shear = .true.

      if (abs(g_exb * g_exbfac) > epsilon(0.)) prp_shear_enabled = .true.
      if (runtype_option_switch == runtype_multibox .and. job == 1) then
         hammett_flow_shear = .false.
      end if

      if (runtype_option_switch == runtype_multibox) then
         call scope(crossdomprocs)
         if (job == 1) then
            call send(g_exbfac * g_exb * x(1), 0, 120)
            call send(g_exbfac * g_exb * x(nx), njobs - 1, 121)
            v_shift = 0.0
         elseif (job == 0) then
            call receive(v_edge, 1, 120)
            v_shift = v_edge - g_exbfac * g_exb * x(1)
         elseif (job == njobs - 1) then
            call receive(v_edge, 1, 121)
            v_shift = v_edge - g_exbfac * g_exb * x(nx)
         end if
         call scope(subprocs)
      end if

      ! Assume we only have a single field line
      ia = 1

      !---------------------------- Allocate arrays ----------------------------
      
      ! Allocate temporary array
      if (radial_variation) allocate (energy(nalpha, -nzgrid:nzgrid))

      ! Allocate the parallel flow shear
      if (.not. allocated(prl_shear)) then
         allocate (prl_shear(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         prl_shear = 0.0
      end if

      ! Allocate the parallel flow shear for radial variation / multi box runs
      if (radial_variation) then
         if (.not. allocated(prl_shear_p)) then
            allocate (prl_shear_p(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         end if
      end if 
      
      !--------------------- Calculate parallel flow shear ---------------------

      ! Iterate over the (z, mu,vpa,s) points
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         is = is_idx(vmu_lo, ivmu)
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         
         ! Calculate parallel flow shear
         !     omega_{zeta,k,s} = -k_y * sqrt(m_s/T_s) * q*a/r * I/B * exp(-v²) * gamma_E
         !     gamma_E = (r/q) (dOmega_zeta(dr) (a/v_{th,ref})
         ! TODO - make formula match code
         do iz = -nzgrid, nzgrid
            prl_shear(ia, iz, ivmu) = -omprimfac * g_exb * code_dt * vpa(iv) * spec(is)%stm_psi0 &
              * dydalpha * drhodpsi * (geo_surf%qinp_psi0 / geo_surf%rhoc_psi0) &
              * (btor(iz) * rmajor(iz) / bmag(ia, iz)) * (spec(is)%mass / spec(is)%temp)
         end do
         
         ! Add the Mawellian exp(v²)
         if (.not. maxwellian_normalization) then
            do iz = -nzgrid, nzgrid
               prl_shear(ia, iz, ivmu) = prl_shear(ia, iz, ivmu) &
                   * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
            end do
         end if
         
         ! Add correction due to radial variation
         if (radial_variation) then
            energy = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
            prl_shear_p(:, :, ivmu) = prl_shear(:, :, ivmu) * (dIdrho / spread(rmajor * btor, 1, nalpha) &
                  - spread(dBdrho, 1, nalpha) / bmag &
                  - spec(is)%fprim - spec(is)%tprim * (energy - 2.5) &
                  - 2.*mu(imu) * spread(dBdrho, 1, nalpha))
         end if
      end do

      ! The definition of parallel flow shear depends on the definition of psi (or x)
      if (q_as_x) prl_shear = prl_shear / geo_surf%shat_psi0

      ! Deallocate the temporary arrays
      if (radial_variation) deallocate (energy)

      !------------------- Calculate perpendicular flow shear ------------------

      ! Allocate module arrays
      if (.not. allocated(shift_times)) allocate (shift_times(naky))
      if (.not. allocated(upwind_advect)) allocate (upwind_advect(naky, nakx))
      
      ! Allocate arrays stored in arrays_store_useful.f90
      if (.not. allocated(shift_state)) then
         allocate (shift_state(naky))
         shift_state = 0.
      end if

      if (nakx > 1 .and. abs(g_exb * g_exbfac) > 0) then
         shift_times = abs(akx(2) / (aky * g_exb * g_exbfac))
      end if
      
      if (zonal_mode(1)) shift_times(1) = huge(0.)

      if (g_exb * g_exbfac > 0.) then
         shift_sign = -1
         shift_start = ikx_max
      else
         shift_sign = 1
         shift_start = ikx_max + 1
      end if

      if (box) upwind_advect = exp(-zi * g_exbfac * g_exb * code_dt * spread(aky, 2, nakx) * spread(x_d, 1, naky))

   end subroutine init_flow_shear

   !****************************************************************************
   !                        Advance parallel flow shear                        
   !****************************************************************************
   subroutine advance_parallel_flow_shear(gout)

      ! Parallelisation
      use mp, only: proc0, mp_abort
      use stella_layouts, only: vmu_lo
      
      ! Fields
      use arrays_store_fields, only: phi, apar, bpar
      
      ! Physics flags
      use parameters_physics, only: full_flux_surface
      
      ! Calculations
      use calculations_kxky_derivatives, only: get_dchidy
      
      ! Grids
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: nakx, naky

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      ! Local variables
      complex, dimension(:, :), allocatable :: g0k
      integer :: ivmu, iz, it, ia

      !-------------------------------------------------------------------------

      ! Assume we only have a single field line
      ia = 1

      ! Allocate temporary arrays
      allocate (g0k(naky, nakx))

      ! Abort for full-flux-surface simulations
      if (full_flux_surface) then
         if (proc0) write (*, *) '!!!WARNING: flow shear not currently supported for full_flux_surface=T!!!'
         call mp_abort("flow shear not currently supported for full_flux_surface=T.")
      end if
      
      ! Iterate over the (z,mu,vpa,tube,s) grid
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
            
               ! Take the derivate w.r.t. y in Fourier space
               ! <g0k> = i ky J_0 ϕ_k = d<chi>_theta/dy = get_dchidy(iz, ivmu, phi, apar, bpar, g0)
               call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k)

               ! Add, - omega_{zeta,k,s} * J0 * phi to the right-hand-side of the gyrokinetic equation
               ! TODO - check the sign of this term?
               gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + prl_shear(ia, iz, ivmu) * g0k
               
            end do
         end do
      end do

      ! Deallocate temporary arrays
      deallocate (g0k)

   end subroutine advance_parallel_flow_shear

   !****************************************************************************
   !                      Advance perpendicular flow shear                      
   !****************************************************************************
   subroutine advance_perp_flow_shear(g)

      ! Parallelisation
      use stella_layouts, only: vmu_lo
      
      ! Calculations
      use constants, only: zi
      use calculations_transforms, only: transform_kx2x_unpadded
      use calculations_transforms, only: transform_x2kx_unpadded
      
      ! Grids
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: aky, zonal_mode
      use grids_kxky, only: nakx, naky, ikx_max
      use stella_time, only: code_dt
      
      ! Radial variation
      use file_utils, only: runtype_option_switch
      use file_utils, only: runtype_multibox
      
      ! Flow shear arrays
      use arrays_store_useful, only: shift_state

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      
      complex, dimension(:, :), allocatable :: g0k, g0x
      real :: shift_fac
      integer :: ivmu, iz, it, iky

      !-------------------------------------------------------------------------

      ! Add the perpendicular flow shear if it is enabled
      if (.not. prp_shear_enabled) return

      ! Allocate temporary arrays
      allocate (g0k(naky, nakx))
      allocate (g0x(naky, nakx))

      ! TODO (DSO) - This assumes the timestep is small enough so that a shift is never more than a single cell
      if (hammett_flow_shear) then
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do iky = 1, naky
                  
                     if (zonal_mode(iky)) cycle
                     
                     if (shift_state(iky) > 0.5 * shift_times(iky)) then
                        
                        ! Shift everything left by one
                        if (shift_sign < 0) then
                           g(iky, (ikx_max + 1):(nakx - 1), iz, it, ivmu) = g(iky, ikx_max + 2:, iz, it, ivmu)
                           g(iky, nakx, iz, it, ivmu) = g(iky, 1, iz, it, ivmu)
                           g(iky, :ikx_max - 1, iz, it, ivmu) = g(iky, 2:ikx_max, iz, it, ivmu)
                           
                        ! Shift everything right by one
                        else
                           g(iky, 2:ikx_max, iz, it, ivmu) = g(iky, 1:(ikx_max - 1), iz, it, ivmu)
                           g(iky, 1, iz, it, ivmu) = g(iky, nakx, iz, it, ivmu)
                           g(iky, ikx_max + 2:, iz, it, ivmu) = g(iky, (ikx_max + 1):(nakx - 1), iz, it, ivmu)
                        end if
                        
                        g(iky, shift_start, iz, it, ivmu) = 0.
                        
                     end if
                  end do
               end do
            end do
         end do
         
      else
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid

                  g0k = g(:, :, iz, it, ivmu)

                  call transform_kx2x_unpadded(g0k, g0x)
                  g0x = upwind_advect * g0x

                  call transform_x2kx_unpadded(g0x, g0k)

                  do iky = 1, naky
                     if (zonal_mode(iky)) cycle
                     if (shift_state(iky) > shift_times(iky)) then
                        g0k(iky, shift_start) = 0.0
                     end if
                  end do

                  g(:, :, iz, it, ivmu) = g0k

               end do
            end do
         end do
      end if

      shift_fac = 1.0
      if (hammett_flow_shear) shift_fac = 0.5

      do iky = 1, naky
         if (zonal_mode(iky)) cycle
         if (shift_state(iky) > shift_fac * shift_times(iky)) then
            shift_state(iky) = shift_state(iky) - shift_times(iky)
         end if
      end do

      shift_state = shift_state + code_dt
      if (zonal_mode(1)) shift_state(1) = 0.

      if (runtype_option_switch == runtype_multibox) then
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g(:, :, iz, it, ivmu) = g(:, :, iz, it, ivmu) * exp(-code_dt * zi * spread(aky, 2, nakx) * v_shift)
               end do
            end do
         end do
      end if

      deallocate (g0k, g0x)

   end subroutine advance_perp_flow_shear

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine finish_flow_shear
   
      use arrays_store_useful, only: shift_state

      implicit none

      if (allocated(prl_shear)) deallocate (prl_shear)
      if (allocated(prl_shear_p)) deallocate (prl_shear_p)
      if (allocated(shift_times)) deallocate (shift_times)
      if (allocated(shift_state)) deallocate (shift_state)
      if (allocated(upwind_advect)) deallocate (upwind_advect)

      initialised_flow_shear = .false.

   end subroutine finish_flow_shear

end module gk_flow_shear
