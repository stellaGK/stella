module flow_shear

  implicit none

  public :: flow_shear_initialized
  public :: shift_state
  public :: init_flow_shear, finish_flow_shear
  public :: prl_shear, prl_shear_p, prp_shear
  public :: advance_flow_shear_explicit, advance_flow_shear_implicit
! public :: add_flow_shear_radial_variation

  private

  logical :: flow_shear_initialized = .false.
  logical :: hammett_flow_shear = .true.

  real, dimension (:,:,:), allocatable :: prl_shear, prl_shear_p
  real, dimension (:), allocatable :: prp_shear, shift_times, shift_state

  integer :: shift_sign, shift_start, shift_end

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
    use zgrid, only: nzgrid
    use kt_grids, only: nalpha, naky, akx, aky, ikx_max
    use stella_geometry, only: q_as_x, geo_surf, bmag, btor, rmajor, dBdrho, dIdrho
    use physics_parameters, only: g_exb, g_exbfac, omprimfac
    use vpamu_grids, only: vperp2, vpa, mu
    use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
    use physics_flags, only: radial_variation

    implicit none

    integer :: is, imu, iv, ivmu, iz, ia
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

    shift_times = abs(akx(2)/(aky*g_exb*g_exbfac))

    if(g_exb*g_exbfac > 0.) then
      shift_sign = -1
      shift_start = ikx_max
      shift_end   = ikx_max + 1
    else
      shift_sign=1
      shift_start = ikx_max + 1
      shift_end   = ikx_max
    endif


  end subroutine init_flow_shear

  subroutine advance_flow_shear_explicit (g, gout)

    use stella_layouts, only: vmu_lo
    use constants, only: zi, pi
    use physics_flags, only: nonlinear
    use physics_parameters, only: g_exb
    use stella_transforms, only: transform_kx2x_solo, transform_x2kx_solo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky, nx, akx,aky, ikx_max, x, x_d, x0
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
    complex, dimension (:), allocatable :: g_shift, dg, shift_in, shift_out

    integer :: ivmu, iz, it, ia, iky, sshear, bb

    real dp

    ia = 1
    bb = 3

    sshear = 1
    if(g_exb < 0) sshear = -1

    allocate (g0k(naky,nakx))
    allocate (g0x(naky,nx))
    allocate (dg(nakx),g_shift(nakx))
    allocate(shift_in(nakx))
    allocate(shift_out(nakx))

    shift_in = exp(pi*zi*x0*akx)
    shift_out = 1./shift_in

    if (.not.allocated(prp_shear)) then 
      allocate(prp_shear(nx))

      prp_shear = x
  
      dp = (prp_shear(bb+1)-prp_shear(nx-bb))/(2.*bb+1)
      do iky = 1, bb
        prp_shear(iky)      =  dp*(iky-0.5)
        prp_shear(nx-iky+1) = -dp*(iky-0.5)

      enddo
      do iky = 1, nx
        write (*,*) prp_shear(iky)
      enddo

      g0x = spread(prp_shear,1,naky)
      call transform_x2kx_solo(g0x,g0k)
      call transform_kx2x_solo(g0k,g0x)

      prp_shear = real(g0x(1,:))
    endif

    do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      do it = 1, ntubes
        do iz = -nzgrid, nzgrid
          call get_dchidy (iz, ivmu, phi(:,:,iz,it), apar(:,:,iz,it), g0k)
            
          !parallel flow shear
          gout(:,:,iz,it,ivmu) = gout(:,:,iz,it,ivmu) + prl_shear(ia,iz,ivmu)*g0k

          !perpendicular flow shear
          !if we use advance_ExB_nonlinearity, we can piggy-back on the 
          !FFTs there and save on some computations
!         if (.not.nonlinear) then
!           call get_dgdy (g(:,:,iz,it,ivmu), g0k)

!           !inverse and forward transforms
!           call transform_kx2x_solo (g0k, g0x)
!           !g0x = -spread(shear,1,naky)*g0x
!           g0x = -g_exb*spread(prp_shear,1,naky)*g0x
!           call transform_x2kx_solo (g0x, g0k)

!           gout(:,:,iz,it,ivmu) = gout(:,:,iz,it,ivmu) + code_dt*g0k
!         endif

          !call get_dgdy (g(:,:,iz,it,ivmu), g0k)
!         g0k = spread(shift_in,1,naky)*g(:,:,iz,it,ivmu)
          do iky=1, naky
            !call third_order_upwind (1,g0k(iky,:),akx(2),sshear,dg)

!           !shift the x grid so that it's properly centered
!           g_shift(1:(ikx_max-1)) = g0k(iky,(ikx_max+1):)
!           g_shift(ikx_max:)      = g0k(iky,1:ikx_max)
!
!           call first_order_upwind (1,g_shift,akx(2),sshear,dg)
!           call third_order_upwind (1,g_shift,akx(2),sshear,dg)
!           call fifth_order_upwind (1,g_shift,akx(2),sshear,dg)
!           call second_order_centered (1,g_shift,akx(2),dg)
!           call four_point_triangle (1,g_shift,akx(2),dg)
!
!           !shift the x grid so that kx=0 is at the first index
!           g_shift(1:ikx_max)    = dg(ikx_max:)
!           g_shift((ikx_max+1):) = dg(1:(ikx_max-1))
!           
!
!           !positive shear leads to advection in the _negative_ kx direction
!           gout(iky,:,iz,it,ivmu) = gout(iky,:,iz,it,ivmu) + code_dt*aky(iky)*g_exb*shift_out*(g_shift)
          enddo
        enddo
      enddo
    enddo

    deallocate (g0k, g0x,dg, g_shift)
    deallocate (shift_in, shift_out)


  end subroutine advance_flow_shear_explicit

  subroutine advance_flow_shear_implicit (g)

    use stella_layouts, only: vmu_lo
    use constants, only: zi
    use physics_parameters, only: g_exb, g_exbfac, g_exb_efold
    use stella_transforms, only: transform_kx2x_solo, transform_x2kx_solo
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky, nx, akx, aky, x, ikx_max
    use stella_time, only: code_dt

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: g
    

    complex, dimension (:,:), allocatable :: g0k, g0x

    integer :: ivmu, iz, it, iky

    allocate (g0k(naky,nakx))
    allocate (g0x(naky,nx))

    if(hammett_flow_shear) then
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            do iky=1,naky
              if(shift_state(iky) > shift_times(iky)) then
                if(shift_sign < 0) then
                  !shift everything left by one
                  g(iky,(ikx_max+1):(nakx-1),iz,it,ivmu) = g(iky,ikx_max+2:,iz,it,ivmu)
                  g(iky,nakx,iz,it,ivmu) = g(iky,1,iz,it,ivmu)
                  g(iky,:ikx_max-1,iz,it,ivmu) = g(iky,2:ikx_max,iz,it,ivmu)
                else
                  !shift everything right by one
                  g(iky,2:ikx_max,iz,it,ivmu) = g(iky,1:(ikx_max-1),iz,it,ivmu)
                  g(iky,1,iz,it,ivmu) = g(iky,nakx,iz,it,ivmu)
                  g(iky,ikx_max+1:,iz,it,ivmu) = g(iky,ikx_max:(nakx-1),iz,it,ivmu)
                endif
                g(iky,shift_start,iz,it,ivmu) = 0.
              endif
            enddo
          enddo
        enddo
      enddo

      do iky=1,naky
        if(shift_state(iky) > shift_times(iky)) &
           shift_state(iky) = shift_state(iky) - shift_times(iky)
      enddo

      shift_state = shift_state + code_dt
    else
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid

            g0k = g(:,:,iz,it,ivmu)

            call transform_kx2x_solo (g0k, g0x)
            g0x = exp(-zi*g_exbfac*g_exb*code_dt*spread(aky,2,nx)*spread(x,1,naky))*g0x
         
            call transform_x2kx_solo (g0x, g0k)

            if(g_exb < 0) then
              g0k(:,ikx_max-1) = exp( 0.5*g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max-1)
              g0k(:,ikx_max)   = exp(     g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max)
!             g0k(:,ikx_max+1) = exp(     g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max+1)
!             g0k(:,ikx_max+2) = exp( 0.5*g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max+2)
            else 
!             g0k(:,ikx_max-1) = exp(-0.5*g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max-1)
!             g0k(:,ikx_max)   = exp(    -g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max)
              g0k(:,ikx_max+1) = exp(    -g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max+1)
              g0k(:,ikx_max+2) = exp(-0.5*g_exb_efold*g_exb*code_dt*aky/akx(2))*g0k(:,ikx_max+2)
            endif

            g(:,:,iz,it,ivmu) = g0k

          enddo
        enddo
      enddo
    endif

    deallocate(g0k,g0x)

  end subroutine advance_flow_shear_implicit

  subroutine finish_flow_shear

    implicit none

    if (allocated(prl_shear)) deallocate (prl_shear)
    if (allocated(prl_shear_p)) deallocate (prl_shear_p)

    flow_shear_initialized = .false.

  end subroutine finish_flow_shear

end module flow_shear
