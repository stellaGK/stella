module flow_shear

  implicit none

  public :: flow_shear_initialized
  public :: init_flow_shear, finish_flow_shear
  public :: prl_shear, prl_shear_p, prp_shear
  public :: advance_flow_shear_explicit, advance_flow_shear_implicit
  public :: shift_state, shift_times
  public :: prp_shear_enabled, hammett_flow_shear
! public :: add_flow_shear_radial_variation

  private

  logical :: flow_shear_initialized = .false.
  logical :: prp_shear_enabled = .true.
  logical :: hammett_flow_shear = .false.
  logical :: prp_shear_implicit = .true.

  complex, dimension (:), allocatable :: shift_in, shift_out
  complex, dimension (:,:), allocatable :: upwind_advect
  real, dimension (:,:), allocatable ::    upwind_diss
  real, dimension (:,:,:), allocatable :: prl_shear, prl_shear_p
  real, dimension (:), allocatable :: prp_shear, shift_times, shift_state

  integer :: shift_sign, shift_start

  real :: v_edge, v_shift = 0.

! integer :: shear_option_switch
! integer, parameter:: shear_option_triangle  = 0, &
!                      shear_option_sine      = 1, &
!                      shear_option_flat      = 2 

contains

  subroutine init_flow_shear
    use stella_layouts, only: vmu_lo
    use stella_layouts, only: iv_idx, imu_idx, is_idx
    use stella_time, only: code_dt
    use species, only: spec
    use constants, only: zi, pi
    use zgrid, only: nzgrid
    use kt_grids, only: x, x_d, x0, nalpha, nx, nakx, naky, akx, aky, ikx_max, zonal_mode
    use stella_geometry, only: q_as_x, geo_surf, bmag, btor, rmajor, dBdrho, dIdrho
    use stella_geometry, only: dydalpha, drhodpsi
    use physics_parameters, only: g_exb, g_exbfac, omprimfac
    use vpamu_grids, only: vperp2, vpa, mu
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use physics_flags, only: radial_variation
    use file_utils, only: runtype_option_switch, runtype_multibox
    use job_manage, only: njobs
    use mp, only: job, send, receive, crossdomprocs, subprocs, scope

    implicit none

    integer :: is, imu, iv, ivmu, iz, ia, shift
    real, dimension (:,:), allocatable :: energy

!   type (text_option), dimension (4), parameter :: shear_opts = &
!     (/ text_option('default',  shear_option_triangle), &
!        text_option('triangle', shear_option_triangle) , &
!        text_option('sine',     shear_option_sine), &
!        text_option('flat',     shear_option_flat)/)
!   character(30) :: zf_option, krook_option, shear_option, LR_debug_option

!   namelist /multibox_parameters/ boundary_size, krook_size, shear_option,& 
!                                  smooth_ZFs, zf_option, LR_debug_option, &
!                                  krook_option, RK_step, nu_krook_mb, &
!                                  mb_debug_step, dealiased_shear
!     call get_option_value & 
!       (shear_option, shear_opts, shear_option_switch, & 
!        ierr, "shear_option in multibox_parameters")


    !perpendicular shear is currently done in multibox

    if (flow_shear_initialized) return
    flow_shear_initialized = .true.

    if(abs(g_exb*g_exbfac) < epsilon(0.)) prp_shear_enabled = .false.

    if (.not.allocated(shift_in))  allocate(shift_in(nakx))
    if (.not.allocated(shift_out)) allocate(shift_out(nakx))

    if(runtype_option_switch .eq. runtype_multibox .and. job.ne.1) then 
      prp_shear_implicit = .false.
    endif

    shift_in = exp(pi*zi*x0*akx)

    if(runtype_option_switch .eq. runtype_multibox) then 

      call scope(crossdomprocs)
      if(job == 1) then
        call send(g_exbfac*g_exb*x(1) ,0,120)
        call send(g_exbfac*g_exb*x(nx),njobs-1,121)
      elseif (job == 0) then
        shift=nx/6
        shift_in = exp(2.0*pi*zi*x0*akx/3.0)
        call receive(v_edge, 1, 120)
        v_shift=v_edge-g_exbfac*g_exb*x(shift)
      elseif (job == njobs-1) then
        shift=nx/6
        shift_in = exp(-2.0*pi*zi*x0*akx/3.0)
        call receive(v_edge, 1, 121)
        v_shift=v_edge-g_exbfac*g_exb*x(nx-shift)
      endif
      call scope(subprocs)
    endif
    
    if(runtype_option_switch .eq. runtype_multibox) then 
      if(job == 0) then
      elseif(job==njobs-1) then
      endif
    endif
    shift_out = 1./shift_in

    ia=1

    !parallel flow shear

    allocate (energy(nalpha,-nzgrid:nzgrid))

    if (.not.allocated(prl_shear)) &
      allocate (prl_shear(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

    if (radial_variation.and..not.allocated(prl_shear_p)) &
      allocate (prl_shear_p(nalpha,-nzgrid:nzgrid,vmu_lo%llim_proc:vmu_lo%ulim_alloc))


    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      is = is_idx(vmu_lo,ivmu)
      iv = iv_idx(vmu_lo,ivmu)
      imu = imu_idx(vmu_lo,ivmu)
      do iz = -nzgrid,nzgrid
        prl_shear(ia,iz,ivmu) = -omprimfac*g_exb*code_dt*vpa(iv)*spec(is)%stm_psi0 &
               * dydalpha*drhodpsi  &
               *(geo_surf%qinp_psi0/geo_surf%rhoc_psi0)  & 
               *(btor(iz)*rmajor(iz)/bmag(ia,iz))*(spec(is)%mass/spec(is)%temp) &
               * maxwell_vpa(iv,is)*maxwell_mu(ia,iz,imu,is)*maxwell_fac(is)
      enddo
      if(radial_variation) then
       energy = (vpa(iv)**2 + vperp2(:,:,imu))*(spec(is)%temp_psi0/spec(is)%temp)
       prl_shear_p(:,:,ivmu) = prl_shear(:,:,ivmu)*(dIdrho/spread(rmajor*btor,1,nalpha) &
                               - spread(dBdrho,1,nalpha)/bmag &
                               - spec(is)%fprim - spec(is)%tprim*(energy-2.5)  &
                               - 2.*mu(imu)*spread(dBdrho,1,nalpha))
       endif
    enddo

    if(q_as_x) prl_shear = prl_shear/geo_surf%shat_psi0

    deallocate(energy)

    !perpendicular flow shear

    if (.not.allocated(shift_times)) allocate (shift_times(naky)); shift_times = 0.
    if (.not.allocated(shift_state)) allocate (shift_state(naky)); shift_state = 0.
!   if (.not.allocated(upwind_advect)) allocate (upwind_advect(naky,nx)); upwind_advect = 0.
    if (.not.allocated(upwind_advect)) allocate (upwind_advect(naky,nakx)); upwind_advect = 0.
    if (.not.allocated(upwind_diss))   allocate (upwind_diss(naky,nx));    upwind_diss  = 0.

    shift_times = abs(akx(2)/(aky*g_exb*g_exbfac))
    if(zonal_mode(1)) shift_times(1) = 0.

    if(g_exb*g_exbfac > 0.) then
      shift_sign = -1
      shift_start = ikx_max
    else
      shift_sign=1
      shift_start = ikx_max + 1
    endif

!   upwind_advect = exp(-zi*g_exbfac*g_exb*code_dt*spread(aky,2,nx)*spread(x,1,naky))
    upwind_advect = exp(-zi*g_exbfac*g_exb*code_dt*spread(aky,2,nakx)*spread(x_d,1,naky))

!   upwind_diss = exp(-abs(g_exbfac*g_exb*code_dt*spread(aky,2,nx)) * &
!                    spread(1. - cos(akx(2)*x),1,naky)/akx(2))
!   upwind_diss = exp(-1.5*abs(g_exbfac*g_exb*code_dt*spread(aky,2,nx)) * &
!                    spread(3. - 4.*cos(akx(2)*x)+cos(akx(3)*x),1,naky)/(6.*akx(2)))
!   upwind_diss = exp(-nx*abs(g_exbfac*g_exb*code_dt*spread(aky,2,nx)) * &
!                    spread(70. - 112.*cos(akx(2)*x) + 56.*cos(akx(3)*x)           &
!                               -  16.*cos(akx(4)*x) +  2.*cos(akx(5)*x),1,naky)   &
!                     /(280.*akx(2)))


  end subroutine init_flow_shear

  subroutine advance_flow_shear_explicit (g, gout)

    use stella_layouts, only: vmu_lo
    use constants, only: zi, pi
    use physics_flags, only: nonlinear
    use physics_parameters, only: g_exb
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky, nx, akx,aky, ikx_max, zonal_mode
    use stella_time, only: code_dt
    use fields, only: get_dchidy
    use fields_arrays, only: phi, apar
    use finite_differences, only: fifth_order_upwind
    use finite_differences, only: third_order_upwind
    use finite_differences, only: first_order_upwind
    use finite_differences, only: second_order_centered
    use finite_differences, only: four_point_triangle

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gout
    

    complex, dimension (:,:), allocatable :: g0k, g0x
    complex, dimension (:), allocatable :: g_shift, dg

    integer :: ivmu, iz, it, ia, iky, sshear

    ia = 1

    sshear = 1
    if(g_exb < 0) sshear = -1

    allocate (g0k(naky,nakx))
    allocate (g0x(naky,nx))
    allocate (dg(nakx),g_shift(nakx))


    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      do it = 1, ntubes
        do iz = -nzgrid, nzgrid
          call get_dchidy (iz, ivmu, phi(:,:,iz,it), apar(:,:,iz,it), g0k)
            
          !parallel flow shear
          gout(:,:,iz,it,ivmu) = gout(:,:,iz,it,ivmu) + prl_shear(ia,iz,ivmu)*g0k

          !perpendicular flow shear

          !call get_dgdy (g(:,:,iz,it,ivmu), g0k)
          if(.not.prp_shear_implicit.and.prp_shear_enabled) then 
            g0k = spread(shift_in,1,naky)*g(:,:,iz,it,ivmu)
            do iky=1, naky

              if(zonal_mode(iky)) cycle

              !shift the x grid so that it's properly centered
              g_shift(1:(ikx_max-1)) = g0k(iky,(ikx_max+1):)
              g_shift(ikx_max:)      = g0k(iky,1:ikx_max)
!
!             call first_order_upwind (1,g_shift,akx(2),sshear,dg)
!             call third_order_upwind (1,g_shift,akx(2),sshear,dg)
!             call second_order_centered (1,g_shift,akx(2),dg)
!             call four_point_triangle (1,g_shift,akx(2),dg)

              call fifth_order_upwind (1,g_shift,akx(2),sshear,dg)

              !shift the x grid so that kx=0 is at the first index
              g_shift(1:ikx_max)    = dg(ikx_max:)
              g_shift((ikx_max+1):) = dg(1:(ikx_max-1))
 
              !positive shear leads to advection in the _negative_ kx direction
              gout(iky,:,iz,it,ivmu) = gout(iky,:,iz,it,ivmu) + code_dt*aky(iky)*g_exb*shift_out*(g_shift) &
                                                              - code_dt*zi*aky(iky)*v_shift*g(iky,:,iz,it,ivmu)
            enddo
          endif
        enddo
      enddo
    enddo

    deallocate (g0k, g0x,dg, g_shift)


  end subroutine advance_flow_shear_explicit

  subroutine advance_flow_shear_implicit (g)

    use stella_layouts, only: vmu_lo
    use constants, only: zi
    use physics_parameters, only: g_exb, g_exbfac, g_exb_efold
!   use stella_transforms, only: transform_kx2x_solo, transform_x2kx_solo
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky, nx, ikx_max, zonal_mode
    use stella_time, only: code_dt

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
    

    complex, dimension (:,:), allocatable :: g0k, g0x

    integer :: ivmu, iz, it, iky

    if(.not.prp_shear_implicit .or..not. prp_shear_enabled) return

    allocate (g0k(naky,nakx))
    allocate (g0x(naky,nakx))
!   allocate (g0x(naky,nx))

    if(hammett_flow_shear) then
      !TODO (DSO) - This assumes the timestep is small enough so that a shift is never
      !             more than a single cell
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            do iky=1,naky
              if(zonal_mode(iky)) cycle
              if(shift_state(iky) > 0.5*shift_times(iky)) then
                if(shift_sign < 0) then
                  !shift everything left by one
                  g(iky,(ikx_max+1):(nakx-1),iz,it,ivmu) = g(iky,ikx_max+2:,iz,it,ivmu)
                  g(iky,nakx,iz,it,ivmu) = g(iky,1,iz,it,ivmu)
                  g(iky,:ikx_max-1,iz,it,ivmu) = g(iky,2:ikx_max,iz,it,ivmu)
                else
                  !shift everything right by one
                  g(iky,2:ikx_max,iz,it,ivmu) = g(iky,1:(ikx_max-1),iz,it,ivmu)
                  g(iky,1,iz,it,ivmu) = g(iky,nakx,iz,it,ivmu)
                  g(iky,ikx_max+2:,iz,it,ivmu) = g(iky,(ikx_max+1):(nakx-1),iz,it,ivmu)
                endif
                g(iky,shift_start,iz,it,ivmu) = 0.
              endif
            enddo
          enddo
        enddo
      enddo
    else
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid

            g0k = g(:,:,iz,it,ivmu)

            call transform_kx2x_unpadded (g0k, g0x)
            g0x = upwind_advect*g0x
!           g0x = upwind_diss  *g0x
         
            call transform_x2kx_unpadded (g0x, g0k)

            do iky = 1, naky
              if(zonal_mode(iky)) cycle
              if(shift_state(iky) > shift_times(iky)) then
                g0k(iky,shift_start)   = 0.0
              endif
            enddo
            if(g_exb < 0) then
!              g0k(:,ikx_max-1) = exp( 0.5*g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max-1)
!              g0k(:,ikx_max)   = exp(     g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max)
!              g0k(:,ikx_max+1) = exp( 0.5*g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max+1)
!             g0k(:,ikx_max+1) = exp(     g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max+1)
!             g0k(:,ikx_max+2) = exp( 0.5*g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max+2)
            else 
!              g0k(:,ikx_max)   = exp(-0.5*g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max)
!              g0k(:,ikx_max+1) = exp(    -g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max+1)
!              g0k(:,ikx_max+2) = exp(-0.5*g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max+2)
!             g0k(:,ikx_max-1) = exp(-0.5*g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max-1)
!             g0k(:,ikx_max)   = exp(    -g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max)
            endif

            g(:,:,iz,it,ivmu) = g0k

          enddo
        enddo
      enddo
    endif

    do iky=1,naky
      if(zonal_mode(iky)) cycle
      if(shift_state(iky) > 0.5*shift_times(iky)) then
        shift_state(iky) = shift_state(iky) - shift_times(iky)
      endif
    enddo

    shift_state = shift_state + code_dt
    if(zonal_mode(1)) shift_state(1) = 0.

    deallocate(g0k,g0x)

  end subroutine advance_flow_shear_implicit

  subroutine finish_flow_shear

    implicit none

    if (allocated(prl_shear)) deallocate (prl_shear)
    if (allocated(prl_shear_p)) deallocate (prl_shear_p)
    if (allocated(shift_times)) deallocate (shift_times)
    if (allocated(shift_state)) deallocate (shift_state)
    if (allocated(shift_in)) deallocate (shift_in)
    if (allocated(shift_out)) deallocate (shift_out)
    if (allocated(upwind_advect)) deallocate (upwind_advect)
    if (allocated(upwind_diss))   deallocate (upwind_diss)

    flow_shear_initialized = .false.

  end subroutine finish_flow_shear

end module flow_shear
