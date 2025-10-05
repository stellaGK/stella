!###############################################################################
!                               HYPER DISSIPATION                               
!###############################################################################
! 
! This module adds hyper dissipation to the gyrokinetic equation.
! 
! A small amount of hyper-viscosity is typically employed in simulations to 
! avoid spectral pile-up. The form for hyper-viscosity is
! 
!   dg/dt = - <d_hyper> * k_perp**4 / k_perp_max**4 * g
! 
!--------------------------------- Questions -----------------------------------
! 
! TODO - If <use_physical_ksqr> = .false. which is used for FFS, the form
! of kperp2 assumes a simple s-alpha model, is this correct for stellarators?
! 
! TODO - Write documentation for advance_hyper_vpa() and advance_hyper_zed().
! 
!--------------------------------- Input file ----------------------------------
! 
!&dissipation_and_collisions_options
!   hyper_dissipation = .false
!/
!&hyper_dissipation
!   d_hyper = 0.05
!   d_zed = 0.05
!   d_vpa = 0.05
!   hyp_zed = .false.
!   hyp_vpa = .false.
!   use_physical_ksqr = .true.
!   scale_to_outboard = .false.
!/
! 
!###############################################################################
module dissipation_hyper

   implicit none

   ! Initialise the hyper dissipation in dissipation_and_collisions.f90
   public :: read_parameters_hyper
   public :: init_hyper
   
   ! Routines for the implicit gyrokinetic equation
   public :: advance_hyper_dissipation
   
   ! Routines for the explicit gyrokinetic equation
   public :: advance_hyper_vpa
   public :: advance_hyper_zed
   
   ! Variables specified in the input file
   public :: D_hyper, D_zed, D_vpa
   public :: hyp_vpa, hyp_zed

   private

   ! Variables specified in the input file
   logical :: use_physical_ksqr, scale_to_outboard
   real :: D_hyper, D_zed, D_vpa
   logical :: hyp_vpa, hyp_zed
   
   ! Variables calculated at initialisation
   real :: k2max  ! Maximum of the perpendicular wanumber squared, i.e., max(kperp**2)
   real :: tfac   ! Factor to define the ballooning angle theta0

contains

!###############################################################################
!################################ READ PARAMETERS ##############################
!###############################################################################

   !****************************************************************************
   !            Read the hyper dissipation namelist in the input file
   !****************************************************************************
   subroutine read_parameters_hyper

      use mp, only: broadcast
      use namelist_dissipation, only: read_namelist_hyper_dissipation

      implicit none

      !-------------------------------------------------------------------------

      ! Read hyper dissipation namelist from input file
      call read_namelist_hyper_dissipation (D_hyper, D_zed, D_vpa, &
         hyp_zed, hyp_vpa, use_physical_ksqr, scale_to_outboard)

      ! Broadcast the parameters to all processors
      call broadcast(use_physical_ksqr)
      call broadcast(scale_to_outboard)
      call broadcast(D_hyper)
      call broadcast(D_zed)
      call broadcast(D_vpa)
      call broadcast(hyp_vpa)
      call broadcast(hyp_zed)

   end subroutine read_parameters_hyper

!###############################################################################
!######################### INITIALISE HYPER DISSIPATION ########################
!###############################################################################

   !****************************************************************************
   !                        Initialise hyper dissipation                        
   !****************************************************************************
   subroutine init_hyper

      ! Grids
      use grids_kxky, only: ikx_max, nakx, naky
      use grids_kxky, only: aky, akx, theta0
      use grids_z, only: nzgrid, zed
      
      ! Geometry
      use geometry, only: geo_surf
      use geometry, only: q_as_x
      use arrays, only: kperp2
      
      ! Warning for VMEC and <scale_to_outboard> = True
      use geometry, only: geo_option_switch, geo_option_vmec
      use vmec_geometry, only: zeta_center, alpha0

      implicit none

      ! Local variables
      integer :: iky, ikx, iz, ia
      real :: temp

      !-------------------------------------------------------------------------

      ! Assume we only have a single field line
      ia = 1
      
      ! Add a warning for the <scale_to_outboard> flag since it assumes alpha=0 and zeta_center=0
      if (scale_to_outboard) then
         if (geo_option_switch==geo_option_vmec) then
            if (zeta_center/=0.0 .or. alpha0/=0.0) then
               write (*,*) 'Warning: <scale_to_outboard> = True, but z=0 does not correspond to the outboard midplance for the chosen field line.'
            end if
         end if
      end if

      ! Use kperp**2(kx,ky,alpha,z) to determine max(kperp**2)
      if (use_physical_ksqr) then
      
         ! Get max(kperp**2) at z=0 which typically corresponds to the outboard midplane
         ! However if alpha!=0  or zeta_center!=0, then z=0 does not correspond to the outboard midplane
         if (scale_to_outboard) then
            k2max = maxval(kperp2(:, :, ia, 0))
            
         ! Get max(kperp**2) along the entire field line
         else
            k2max = maxval(kperp2)
         end if
         
      ! Avoid spatially dependent kperp (through the geometric coefficients).
      ! Nonetheless, kperp is still allowed to vary along zed with global magnetic shear.
      ! This is useful for full_flux_annulus and radial_variation runs.
      ! TODO-GA: Is the s-alpha model correct for FFS? Or is scale_to_outboard set to True always?
      else
      
         ! Define the ballooning angle as theta0 = kx / (ky*shat)
         ! and we define tfac = (1 / theta0**2) * (kx**2 / ky**2) = shat**2
         tfac = geo_surf%shat**2

         ! If x = q we define the ballooning angle as theta0 = kx / ky
         ! and we define tfac = (1 / theta0**2) * (kx**2 / ky**2) = 1
         if (q_as_x) tfac = 1.0

         ! Get max(kperp**2) at the outboard midplane
         ! Here we basically assume that |∇x|² = |∇y|² = 1 and ∇x . ∇y = 0
         if (scale_to_outboard) then
            k2max = akx(ikx_max)**2 + aky(naky)**2
            
         ! Get max(kperp**2) along the field line
         ! Here we use a simple s-alpha model which gives:
         ! k_perp**2 = ky**2 ( 1 + hat{s}**2 (theta - theta0)**2 )
         else
            k2max = -1.0
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  do iky = 1, naky
                     temp = aky(iky)**2 * (1.0 + tfac * (zed(iz) - theta0(iky, ikx))**2)
                     if (temp > k2max) k2max = temp
                  end do
               end do
            end do
         end if
      end if
      
      ! If we failed to set max(kperp**2), set it to 1.0
      if (k2max < epsilon(0.0)) k2max = 1.0

   end subroutine init_hyper
   
!###############################################################################
!########################### ADVANCE HYPER DISSIPATION #########################
!###############################################################################

   !****************************************************************************
   !             Add hyper dissipation to the gyrokinetic equation              
   !****************************************************************************
   subroutine advance_hyper_dissipation(g)

      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      
      ! Grids
      use grids_kxky, only: aky, akx, naky
      use grids_kxky, only: theta0, zonal_mode
      use grids_time, only: code_dt
      use grids_z, only: nzgrid, ntubes, zed
      
      ! Geometry
      use arrays, only: kperp2

      implicit none

      ! The distribution function g(kx,ky,z,tube,imuvpas)
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g

      ! Local variables
      integer :: ia, ivmu, iz, it, iky

      !-------------------------------------------------------------------------

      ! Assume we only have a single field line
      ia = 1

     ! Add hyper-dissipation of the form dg/dt = -D * (k/kmax)^4 * g
     if (use_physical_ksqr) then
     
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            g(:, :, :, :, ivmu) = g(:, :, :, :, ivmu) / (1. + code_dt * (spread(kperp2(:, :, ia, :), 4, ntubes) / k2max)**2 * D_hyper)
         end do
         
      ! Avoid spatially dependent kperp if <use_physical_ksqr> = .false.
      ! Use a simple s-alpha model: k_perp**2 = ky**2 ( 1 + hat{s}**2 (theta - theta0)**2 )
      ! Add in hyper-dissipation of the form dg/dt = -D * (k/kmax)^4 * g
      else
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do iky = 1, naky
                     if (zonal_mode(iky)) then
                        g(iky, :, iz, it, ivmu) = g(iky, :, iz, it, ivmu) / (1.+code_dt * (akx(:)**2 / k2max)**2 * D_hyper)
                     else
                        g(iky, :, iz, it, ivmu) = g(iky, :, iz, it, ivmu) / &
                            (1. + code_dt * (aky(iky)**2 * (1.0 + tfac * (zed(iz) - theta0(iky, :))**2) / k2max)**2 * D_hyper)
                     end if
                  end do
               end do
            end do
         end do 
      end if

   end subroutine advance_hyper_dissipation
 
   !****************************************************************************
   !                                      Title
   !****************************************************************************
   ! Computes the fourth derivative of g in vpa and returns this in dgdvpa
   ! multiplied by the vpa diffusion coefficient
   !****************************************************************************
   subroutine advance_hyper_vpa(g, dgdvpa)
   
      use grids_time, only: code_dt
      use grids_z, only: nzgrid
      use parallelisation_layouts, only: vmu_lo, kxkyz_lo
      use redistribute, only: gather, scatter
      use initialise_redistribute, only: kxkyz2vmu
      use grids_velocity, only: nmu, nvpa, dvpa

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dgdvpa

      ! Local variables
      complex, dimension(:, :, :), allocatable :: g0v, g1v

      !-------------------------------------------------------------------------

      allocate (g0v(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      allocate (g1v(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

      call scatter(kxkyz2vmu, g, g0v)
      call get_dgdvpa_fourth_order(g0v, g1v)
      call gather(kxkyz2vmu, g1v, dgdvpa)
      dgdvpa = -code_dt * D_vpa * dvpa**4 / 16 * dgdvpa
      
      deallocate (g0v)
      deallocate (g1v)

   end subroutine advance_hyper_vpa

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine advance_hyper_zed(g, dgdz)
   ! computes the fourth derivative of g in z and returns this in
   ! dgdz multiplied by the z hyper diffusion coefficient
      use grids_time, only: code_dt
      use grids_z, only: nzgrid, delzed
      use parallelisation_layouts, only: vmu_lo

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dgdz

      !-------------------------------------------------------------------------

      call get_dgdz_fourth_order(g, dgdz)
      dgdz = -code_dt * D_zed * delzed(0)**4 / 16 * dgdz

   end subroutine advance_hyper_zed
   
!###############################################################################
!################################## CALCULATIONS ###############################
!###############################################################################

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine get_dgdvpa_fourth_order(g, gout)

      use calculations_finite_differences, only: fourth_derivate_second_centered_vpa
      use parallelisation_layouts, only: kxkyz_lo, iz_idx, is_idx
      use grids_velocity, only: nvpa, nmu, dvpa

      implicit none

      ! Arguments
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(inout) :: gout

      ! Local variables
      integer :: ikxkyz, imu, iz, is
      complex, dimension(:), allocatable :: tmp

      !-------------------------------------------------------------------------

      allocate (tmp(nvpa))
      
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            call fourth_derivate_second_centered_vpa(1, g(:, imu, ikxkyz), dvpa, tmp)
            gout(:, imu, ikxkyz) = tmp
         end do
      end do

      deallocate (tmp)
      
   end subroutine get_dgdvpa_fourth_order

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine get_dgdz_fourth_order(g, dgdz)

      use calculations_finite_differences, only: fourth_derivative_second_centered_zed
      use parallelisation_layouts, only: vmu_lo
      use parallelisation_layouts, only: iv_idx
      use grids_z, only: nzgrid, delzed, ntubes
      use grids_extended_zgrid, only: neigen, nsegments
      use grids_extended_zgrid, only: iz_low, iz_up
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: fill_zed_ghost_zones
      use grids_extended_zgrid, only: periodic
      use grids_kxky, only: naky
      use parallelisation_layouts, only: iv_idx, imu_idx, is_idx

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(inout) :: dgdz

      ! Local variables
      integer :: iseg, ie, iky, iv, it, ivmu, imu
      complex, dimension(2) :: gleft, gright

      !-------------------------------------------------------------------------
      
      ! FLAG -- assuming delta zed is equally spaced below!
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         do iky = 1, naky
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  do iseg = 1, nsegments(ie, iky)
                  
                     ! Fill in ghost zones at boundaries in g(z)
                     call fill_zed_ghost_zones(it, iseg, ie, iky, g(:, :, :, :, ivmu), gleft, gright)
                     
                     ! Get dg/dz
                     call fourth_derivative_second_centered_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
                        g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it, ivmu), &
                        delzed(0), gleft, gright, periodic(iky), &
                        dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it, ivmu))
                        
                  end do
               end do
            end do
         end do
      end do
   end subroutine get_dgdz_fourth_order

end module dissipation_hyper
