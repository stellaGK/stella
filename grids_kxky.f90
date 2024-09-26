module grids_kxky

   implicit none

   public :: init_grids_kxky, finish_grids_kxky

   public :: aky, akx
   public :: aky_all, aky_all_ordered
   public :: theta0, zed0
   public :: zonal_mode

   public :: x, x_d, y
   public :: dy, dx
   public :: boundary_size, copy_size, krook_size 
   
   public :: rho, rho_d, rho_clamped, rho_d_clamped
   public :: g0x
   public :: box
   
   private 

   real, dimension(:), allocatable :: aky, akx
   real, dimension(:), allocatable :: aky_all, aky_all_ordered
   real, dimension(:, :), allocatable :: theta0, zed0
   real, dimension(:), allocatable :: x, x_d, y
   logical, dimension(:), allocatable :: zonal_mode

   !> For radial variation
   real, dimension(:), allocatable :: rho, rho_d, rho_clamped, rho_d_clamped
   complex, dimension(:, :), allocatable :: g0x

   !> Not allocated here - need to work out where to put these 
   integer :: boundary_size, copy_size, krook_size

   !> Required for flux calculations
   real :: dx, dy

   !> Internal Calculations
   real :: dkx, dky, dx_d
   logical :: box
   
   logical :: initialised

contains

   !======================================================================
   !======================= INITIALISE KXKY GRIDS ========================
   !======================================================================
   subroutine init_grids_kxky

       use mp, only: mp_abort
       use common_types, only: flux_surface_type
       use zgrid, only: init_zgrid
       use parameters_kxky_grids, only: gridopt_switch, gridopt_range, gridopt_box
       use parameters_kxky_grids, only: naky
       
       implicit none
 
       if (initialised) return

 
       call init_zgrid
       
       select case (gridopt_switch)
       case (gridopt_range)
          call init_grids_kxky_range
       case (gridopt_box)
          call init_grids_kxky_box
       end select

       !> determine if iky corresponds to zonal mode
       if (.not. allocated(zonal_mode)) allocate (zonal_mode(naky))
       zonal_mode = .false.
       if (abs(aky(1)) < epsilon(0.)) zonal_mode(1) = .true.
 
       initialised = .true.

   contains

       !**********************************************************************
       !                       INITIALISE RANGE KXKY GRID                    !
       !**********************************************************************
       !> DESCRIPTION
       !**********************************************************************

       subroutine init_grids_kxky_range

           use common_types, only: flux_surface_type
           use geometry, only: geo_surf, q_as_x
           use zgrid, only: shat_zero

           use parameters_kxky_grids, only: naky, nakx, aky_min, aky_max, &
                akx_min, akx_max, theta0_min, theta0_max, &
                kyspacing_option_switch, kyspacing_linear, kyspacing_exponential, &
                ikx_max, naky_all
           
           implicit none
   
           integer :: i, j
           real :: dkx, dky, dtheta0, tfac
           real :: zero
   
           box = .false.
           call allocate_arrays
           
           ! NB: we are assuming here that all ky are positive
           ! when running in range mode
           dky = 0.0
           
           if (naky > 1) then
           select case (kyspacing_option_switch)
           case (kyspacing_linear)
               dky = (aky_max - aky_min) / real(naky - 1)
               aky = (/(aky_min + dky * real(i), i=0, naky - 1)/)
           case (kyspacing_exponential)
               dky = (log(aky_max) - log(aky_min)) / real(naky - 1)
               aky = (/(exp(log(aky_min) + dky * real(i)), i=0, naky - 1)/)
           end select
           else
           aky = (/(aky_min, i=0, naky - 1)/)
           end if
   
           ! set default kx and theta0 to 0
           akx = 0.0; theta0 = 0.0
   
           if (q_as_x) then
           tfac = 1.0
           else
           tfac = geo_surf%shat
           end if
   
           zero = 100.*epsilon(0.)
   
           ! if theta0_min and theta0_max have been specified,
           ! use them to determine akx_min and akx_max
           if (theta0_max > theta0_min - zero) then
              if (geo_surf%shat > epsilon(0.)) then
                 akx_min = theta0_min * tfac * aky(1)
                 akx_max = theta0_max * tfac * aky(1)
              else
                 akx_min = theta0_max * tfac * aky(1)
                 akx_max = theta0_min * tfac * aky(1)
              end if
           end if
   
           ! shat_zero is minimum shat value below which periodic BC is enforced
           if (abs(geo_surf%shat) > shat_zero) then  ! ie assumes boundary_option .eq. 'linked'
           ! if kx_min and akx_max specified in input
              ! instead of theta0_min and theta0_max,
              ! use them to get theta0_min and theta0_max
              if (theta0_min > theta0_max + zero .and. abs(aky(1)) > zero) then
                 theta0_min = akx_min / (tfac * aky(1))
                 theta0_max = akx_max / (tfac * aky(1))
                 dtheta0 = 0.0
                 if (nakx > 1) dtheta0 = (theta0_max - theta0_min) / real(nakx - 1)
   
                 do j = 1, naky
                    theta0(j, :) &
                         = (/(theta0_min + dtheta0 * real(i), i=0, nakx - 1)/)
                 end do
                 akx= theta0(1, :) * tfac * aky(1)
              else if (akx_max > akx_min - zero .or. nakx == 1) then
                 dkx = 0.0
                 if (nakx > 1) dkx = (akx_max - akx_min) / real(nakx - 1)
                 akx = (/(akx_min + dkx * real(i), i=0, nakx - 1)/)
                 
                 dtheta0 = 0.0
                 if (nakx > 1) dtheta0 = (theta0_max - theta0_min) / real(nakx - 1)
                 
                 if (geo_surf%shat > epsilon(0.)) then
                    do j = 1, naky
                       theta0(j, :) &
                            = (/(theta0_min + dtheta0 * real(i), i=0, nakx - 1)/)
                    end do
                 else
                    do j = 1, naky
                       theta0(j, :) &
                            = (/(theta0_min + dtheta0 * real(i), i=nakx - 1, 0, -1)/)
                    end do
                 end if
              else
                 call mp_abort('ky=0 is inconsistent with akx_min different from akx_max. aborting.')
              end if
   
           else
           ! here assume boundary_option .eq. 'periodic'
           ! used for periodic finite kx ballooning space runs with shat=0
              dkx = 0.0
              if (nakx > 1) dkx = (akx_max - akx_min) / real(nakx - 1)
              akx = (/(akx_min + dkx * real(i), i=0, nakx - 1)/)
           end if
   
           zed0 = theta0 * geo_surf%zed0_fac

           ikx_max = nakx
           naky_all = naky
           
       end subroutine init_grids_kxky_range

    
       !**********************************************************************
       !                       INITIALISE BOX KXKY GRID                    !
       !**********************************************************************
       !> DESCRIPTION
       !**********************************************************************
       subroutine init_grids_kxky_box
           
           use mp, only: mp_abort, proc0, broadcast
           use common_types, only: flux_surface_type
           use constants, only: pi, zi
           use geometry, only: geo_surf, twist_and_shift_geo_fac, dydalpha
           use geometry, only: q_as_x, get_x_to_rho, dxdpsi, drhodpsi
           use geometry, only: geo_option_switch, geo_option_vmec
           use parameters_physics, only: rhostar
           use parameters_physics, only: full_flux_surface, radial_variation
           use file_utils, only: runtype_option_switch, runtype_multibox
           use zgrid, only: nperiod
           use zgrid, only: boundary_option_switch, boundary_option_linked
           use zgrid, only: boundary_option_linked_stellarator
           use ran, only: ranf

           use write_radial_grid, only: dump_radial_grid

           use parameters_kxky_grids, only: nx, ny, ikx_max, naky_all, naky, nakx, &
                x0, y0, jtwist, jtwistfac, phase_shift_angle, ikx_twist_shift, &
                centered_in_rho, randomize_phase_shift, periodic_variation
           
           implicit none
     
           integer :: ikx, iky
           integer :: ikyneg
           real :: x_shift, dqdrho, pfac, norm
     
           box = .true.
           call allocate_arrays
           
           !> set jtwist and y0 for cases where they have not been specified
           !> and for which it makes sense to set them automatically
           if (jtwist < 1) then
              jtwist = max(1, int(abs(twist_and_shift_geo_fac) + 0.5))
              jtwist = max(1, int(jtwistfac * jtwist + 0.5))
           end if
           !> signed version of jtwist, with sign determined by, e.g., magnetic shear
           ikx_twist_shift = -jtwist * int(sign(1.0, twist_and_shift_geo_fac))
     
           if (y0 < 0.) then
              if (full_flux_surface) then
                 !> if simulating a flux annulus, then
                 !> y0 determined by the physical
                 !> extent of the device
                 if (rhostar > 0.) then
                    !y0 = 1./(rhostar*geo_surf%rhotor)
                    y0 = geo_surf%rhotor / rhostar
                 else
                    call mp_abort('must set rhostar if simulating a full flux surface. aborting.')
                 end if
              else
                 !> if simulating a flux tube
                 !> makes no sense to have y0 < 0.0
                 !> so abort
                 call mp_abort('y0 negative only makes sense when simulating a flux annulus.  aborting.')
              end if
           end if

           !> get the grid spacing in ky and then in kx using twist-and-shift BC
           dky = 1./y0
           
           ! kx = ky * twist_shift_geo_fac / jtwist for every linked boundary condition
           ! except for the periodic ones
           select case (boundary_option_switch)
           case (boundary_option_linked)
              dkx = (2 * nperiod - 1) * dky * abs(twist_and_shift_geo_fac) / real(jtwist)
           case (boundary_option_linked_stellarator)
              dkx = dky * abs(twist_and_shift_geo_fac) / real(jtwist)
           case default
              if (x0 < epsilon(0.0)) then
                 dkx = dky
              else
                 dkx = 1./x0
              end if
           end select
           
           x0 = 1./dkx

           !> ky goes from zero to aky_max
           do iky = 1, naky
              aky(iky) = real(iky - 1) * dky
           end do
     
           !> aky_all contains all ky values (positive and negative),
           !> stored in the same order as akx (0 -> aky_max, -aky_max -> -dky)
           !> first set arrays equal for ky >= 0
           aky_all(:naky) = aky
           !> aky_all_ordered contains all ky values, stored from
           !> most negative to most positive (-aky_max -> aky_max)
           aky_all_ordered(naky:naky_all) = aky
           !> next fill in ky < 0
           do iky = naky + 1, naky_all
              !> this is the ky index corresponding to +ky in original array
              ikyneg = naky_all - iky + 2
              aky_all(iky) = -aky(ikyneg)
           end do
           aky_all_ordered(:naky - 1) = aky_all(naky + 1:)
     
           !> kx goes from zero to kx_max down to zero...
           do ikx = 1, ikx_max
              akx(ikx) = real(ikx - 1) * dkx
           end do
           !> and then from -kx_max to -|akx_min|
           do ikx = ikx_max + 1, nakx
              akx(ikx) = real(ikx - nakx - 1) * dkx
           end do
     
           !> set theta0=0 for ky=0
           theta0(1, :) = 0.0
           if (q_as_x) then
              do ikx = 1, nakx
                 !> theta0 = kx/ky
                 theta0(2:, ikx) = akx(ikx) / aky(2:)
              end do
           else
              do ikx = 1, nakx
                 !> theta0 = kx/ky/shat
                 theta0(2:, ikx) = akx(ikx) / (aky(2:) * geo_surf%shat)
              end do
           end if

           norm = 1.
           if (naky > 1) norm = aky(2)
           if (rhostar > 0.) then
              if (geo_option_switch == geo_option_vmec) then
                 phase_shift_angle = -2.*pi * (2 * nperiod - 1) * dydalpha / (rhostar * geo_surf%qinp_psi0)
              else
                 phase_shift_angle = -2.*pi * (2 * nperiod - 1) * geo_surf%qinp_psi0 * dydalpha / rhostar
              end if
           else if (randomize_phase_shift) then
              if (proc0) phase_shift_angle = 2.*pi * ranf() / norm
              call broadcast(phase_shift_angle)
           else
              phase_shift_angle = 2.*pi * phase_shift_angle / norm
           end if
     
           !> MAB: a lot of the radial variation coding below should probably be tidied away
           !> into one or more separate subroutines
    
           dx = (2 * pi * x0) / nx
           dy = (2 * pi * y0) / ny
     
           x_shift = pi * x0
           pfac = 1.0
           if (periodic_variation) pfac = 0.5
           if (centered_in_rho) then
              if (q_as_x) then
                 dqdrho = geo_surf%shat * geo_surf%qinp / geo_surf%rhoc
                 x_shift = pi * x0 * (1.0 - 0.5 * pfac * rhostar * pi * x0 * geo_surf%d2qdr2 / (dqdrho**2 * dxdpsi))
              else
                 x_shift = pi * x0 * (1.0 - 0.5 * pfac * rhostar * pi * x0 * geo_surf%d2psidr2 * drhodpsi**2 / dxdpsi)
              end if
           end if
     
           do ikx = 1, nx
              if (radial_variation .or. runtype_option_switch == runtype_multibox) then
                 if (periodic_variation) then
                    if (ikx <= nx / 2) then
                       x(ikx) = (ikx - 1) * dx - 0.5 * x_shift
                    else
                       x(ikx) = x(nx - ikx + 1)
                    end if
                 else
                    x(ikx) = (ikx - 0.5) * dx - x_shift
                 end if
              else
                 x(ikx) = (ikx - 1) * dx
              end if
           end do
     
           dx_d = (2 * pi * x0) / nakx
           do ikx = 1, nakx
              if (radial_variation .or. runtype_option_switch == runtype_multibox) then
                 if (periodic_variation) then
                    if (ikx <= (nakx + 1) / 2) then
                       x_d(ikx) = (ikx - 1) * dx_d - 0.5 * x_shift
                    else
                       x_d(ikx) = x_d(nakx - ikx + 1)
                    end if
                 else
                    x_d(ikx) = (ikx - 0.5) * dx_d - x_shift
                 end if
              else
                 x_d(ikx) = (ikx - 1) * dx_d
              end if
           end do
     
           call get_x_to_rho(1, x, rho)
           call get_x_to_rho(1, x_d, rho_d)

           rho_clamped = rho
           rho_d_clamped = rho_d
           
           zed0 = theta0 * geo_surf%zed0_fac
     
           if (radial_variation) call dump_radial_grid (x, rho, nx)
     
           if (radial_variation .and. (any((rho + geo_surf%rhoc) < 0.0) &
                                       .or. any((rho + geo_surf%rhoc) > 1.0))) then
              call mp_abort('rho(x) is beyond range [0,1]. Try changing rhostar or q/psi profiles')
           end if
     
           do iky = 1, ny
              y(iky) = (iky - 1) * dy
           end do

           call broadcast_parameters
           
       end subroutine init_grids_kxky_box
       
       !**********************************************************************
       !                     ALLOCATE ARRAYS FOR KXKY GRIDS                  !
       !**********************************************************************
       !> DESCRIPTION
       !**********************************************************************

       subroutine allocate_arrays
         
         use parameters_kxky_grids, only: nakx, naky, naky_all, nx , ny 
         
         implicit none
         
         if(.not. allocated(akx)) allocate (akx(nakx))
         if(.not. allocated(aky)) allocate (aky(naky))
         if(.not. allocated(aky_all)) allocate(aky_all(naky_all))
         if(.not. allocated(aky_all_ordered)) allocate (aky_all_ordered(naky_all))
         if(.not. allocated(theta0)) allocate (theta0(naky, nakx))
         if(.not. allocated(zed0)) allocate (zed0(naky, nakx))
         
         if (box) then
            if (.not. allocated(x_d)) allocate (x_d(nakx))
            if (.not. allocated(rho)) allocate (rho(nx))
            if (.not. allocated(rho_d)) allocate (rho_d(nakx))
            if (.not. allocated(rho_clamped)) allocate (rho_clamped(nx))
            if (.not. allocated(rho_d_clamped)) allocate (rho_d_clamped(nakx))
            
            if (.not. allocated(x)) allocate (x(nx))
            if (.not. allocated(y)) allocate (y(ny))
            
         end if
         
       end subroutine allocate_arrays
       
       
     end subroutine init_grids_kxky
          
     subroutine broadcast_parameters
       
       use mp, only: broadcast

       implicit none

       ! call broadcast(akx)
       ! call broadcast(aky)
       ! call broadcast(aky_all)
       ! call broadcast(aky_all_ordered)
       ! call broadcast(theta0)
       ! call broadcast(zed0)

!       if(box) then
       call broadcast(x_d)
       call broadcast(rho)
       call broadcast(rho_d)
       call broadcast(rho_clamped)
       call broadcast(rho_d_clamped)
!          call broadcast(x)
!          call broadcast(y)
!       end if
             
     end subroutine broadcast_parameters
   
   !======================================================================
   !======================= INITIALISE KXKY GRIDS ========================
   !======================================================================

   subroutine finish_grids_kxky

       use parameters_kxky_grids, only: reality 
       implicit none
 
       if (allocated(aky)) deallocate (aky)
       if (allocated(aky_all)) deallocate (aky_all)
       if (allocated(aky_all_ordered)) deallocate (aky_all_ordered)
       if (allocated(akx)) deallocate (akx)
       if (allocated(theta0)) deallocate (theta0)
       if (allocated(zed0)) deallocate (zed0)
 
       if (allocated(x)) deallocate (x)
       if (allocated(y)) deallocate (y)
 
       if (allocated(x_d)) deallocate (x_d)
       if (allocated(rho)) deallocate (rho)
       if (allocated(rho_d)) deallocate (rho_d)
       if (allocated(rho_clamped)) deallocate (rho_clamped)
       if (allocated(rho_d_clamped)) deallocate (rho_d_clamped)
 
       if (allocated(g0x)) deallocate (g0x)

       reality = .false. 
       initialised = .false.
     end subroutine finish_grids_kxky

end module grids_kxky
