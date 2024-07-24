module stella_geometry

   use common_types, only: flux_surface_type

   implicit none

   public :: init_geometry, finish_init_geometry, finish_geometry
   public :: communicate_geo_multibox
   public :: grho, grho_norm, grad_x
   public :: bmag, dbdzed, btor, bmag_psi0
   public :: gradpar, gradpar_eqarc, b_dot_grad_z, zed_eqarc
   public :: cvdrift, cvdrift0
   public :: gbdrift, gbdrift0
   public :: dcvdriftdrho, dcvdrift0drho
   public :: dgbdriftdrho, dgbdrift0drho
   public :: gds2, gds21, gds22, gds23, gds24, gds25, gds26
   public :: dgds2dr, dgds21dr, dgds22dr
   public :: exb_nonlin_fac, exb_nonlin_fac_p
   public :: jacob, djacdrho
   public :: drhodpsi, drhodpsi_psi0
   public :: dl_over_b, d_dl_over_b_drho
   public :: dBdrho, d2Bdrdth, dgradpardrho, dIdrho
   public :: geo_surf
   public :: Rmajor
   public :: alpha
   public :: theta_vmec
   public :: zeta
   public :: zed_scalefac
   public :: dxdXcoord, dydalpha
   public :: sign_torflux
   public :: aref, bref
   public :: twist_and_shift_geo_fac
   public :: q_as_x, get_x_to_rho, gfac
   public :: dVolume
   public :: grad_x_grad_y_end
   public :: x_displacement_fac
   public :: geo_option_switch, geo_option_vmec

   public :: gradzeta_gradx, gradzeta_grady
   public :: gradpar_zeta

   private

   type(flux_surface_type) :: geo_surf
   
   ! TODO-GA: import this from run_parameters as "use run_parameters, only: print_extra_info_to_terminal"
   logical :: print_extra_info_to_terminal = .false.

   real :: grad_x_grad_y_end
   real :: aref, bref
   real :: dxdXcoord, dydalpha
   real :: dqdrho
   real :: dIdrho
   real :: grho_norm
   real :: drhodpsi, drhodpsi_psi0, shat, qinp
   real :: exb_nonlin_fac, exb_nonlin_fac_p
   real :: gradpar_eqarc
   real :: zed_scalefac
   real :: twist_and_shift_geo_fac, gfac
   real, dimension(:), allocatable :: zed_eqarc
   real, dimension(:), allocatable :: gradpar
   real, dimension(:, :), allocatable :: b_dot_grad_z
   real, dimension(:, :), allocatable :: bmag, bmag_psi0, dbdzed
   real, dimension(:, :), allocatable :: twist_and_shift_geo_fac_full
   real, dimension(:, :), allocatable :: cvdrift, cvdrift0
   real, dimension(:, :), allocatable :: gbdrift, gbdrift0
   real, dimension(:, :), allocatable :: dcvdriftdrho, dcvdrift0drho
   real, dimension(:, :), allocatable :: dgbdriftdrho, dgbdrift0drho
   real, dimension(:, :), allocatable :: gds2, gds21, gds22, gds23, gds24, gds25, gds26
   real, dimension(:, :), allocatable :: dgds2dr, dgds21dr
   real, dimension(:, :), allocatable :: dgds22dr
   real, dimension(:, :), allocatable :: theta_vmec
   real, dimension(:, :), allocatable :: jacob, djacdrho, grho, grad_x
   real, dimension(:, :), allocatable :: dl_over_b, d_dl_over_b_drho
   real, dimension(:, :, :), allocatable :: dVolume
   real, dimension(:, :), allocatable :: x_displacement_fac
   real, dimension(:), allocatable :: dBdrho, d2Bdrdth, dgradpardrho
   real, dimension(:), allocatable :: btor, Rmajor
   real, dimension(:), allocatable :: alpha
   real, dimension(:, :), allocatable :: zeta

   integer :: sign_torflux
   integer :: geo_option_switch
   integer, parameter :: geo_option_local = 1
   integer, parameter :: geo_option_inputprof = 2
   integer, parameter :: geo_option_vmec = 3
   integer, parameter :: geo_option_multibox = 4

   logical :: overwrite_geometry
   logical :: overwrite_bmag, overwrite_gradpar
   logical :: overwrite_gds2, overwrite_gds21, overwrite_gds22
   logical :: overwrite_gds23, overwrite_gds24
   logical :: overwrite_gbdrift, overwrite_cvdrift, overwrite_gbdrift0
   logical :: q_as_x
   character(100) :: geo_file

   logical :: vmec_chosen = .false.
   logical :: geoinit = .false.
   logical :: set_bmag_const

   real, dimension(:, :), allocatable :: gradzeta_gradx, gradzeta_grady
   real, dimension(:, :), allocatable :: gradpar_zeta
contains

   subroutine init_geometry(nalpha, naky)

      use constants, only: pi
      use mp, only: proc0
      use millerlocal, only: read_local_parameters, get_local_geo
      use millerlocal, only: communicate_parameters_multibox
      use vmec_geo, only: read_vmec_parameters, get_vmec_geo
      use vmec_to_stella_geometry_interface_mod, only: desired_zmin
      use inputprofiles_interface, only: read_inputprof_geo
      use zgrid, only: nzed, nzgrid
      use zgrid, only: zed, delzed
      use zgrid, only: shat_zero, zed_equal_arc
      use zgrid, only: grad_x_grad_y_zero
      use zgrid, only: boundary_option_switch, boundary_option_self_periodic
      use zgrid, only: boundary_option_linked, boundary_option_linked_stellarator
      use zgrid, only: dkx_over_dky
      use file_utils, only: get_unused_unit
      use physics_flags, only: include_geometric_variation, const_alpha_geo

      implicit none

      logical, parameter :: debug = .false.

      integer, intent(in) :: nalpha, naky
      logical :: stellarator_symmetric_BC

      real :: dpsidrho, dpsidrho_psi0
      real :: new_zeta_min
      integer :: iy, ia, iz
      integer :: dxdXcoord_sign, dydalpha_sign
      real :: field_period_ratio, bmag_z0
      real, dimension(:, :), allocatable :: grad_alpha_grad_alpha
      real, dimension(:, :), allocatable :: grad_alpha_grad_psi
      real, dimension(:, :), allocatable :: grad_psi_grad_psi
      real, dimension(:, :), allocatable :: gbdrift_alpha, cvdrift_alpha
      real, dimension(:, :), allocatable :: gbdrift0_psi, cvdrift0_psi

      if (geoinit) return
      geoinit = .true.

      ! B = grad alpha x grad psi
      ! for tokamak calculations, alpha = zeta - q * theta
      ! and psi = psi_poloidal
      ! for stellarator calculations, alpha = theta - iota * zeta
      ! and psi = -psi_toroidal

      ! default is no re-scaling of zed
      zed_scalefac = 1.0

      if (proc0) then
         call read_parameters
         select case (geo_option_switch)
         case (geo_option_local)
            ! read in Miller local parameters
            call read_local_parameters(nzed, nzgrid, geo_surf)
            ! allocate geometry arrays
            call allocate_arrays(nalpha, nzgrid)
            ! use Miller local parameters to get
            ! geometric coefficients needed by stella
            call get_local_geo(nzed, nzgrid, zed, zed_equal_arc, &
                               dpsidrho, dpsidrho_psi0, dIdrho, grho(1, :), &
                               bmag(1, :), bmag_psi0(1, :), &
                               gds2(1, :), gds21(1, :), gds22(1, :), &
                               gds23(1, :), gds24(1, :), gradpar, &
                               gbdrift0(1, :), gbdrift(1, :), cvdrift0(1, :), cvdrift(1, :), &
                               dBdrho, d2Bdrdth, dgradpardrho, btor, rmajor, &
                               dcvdrift0drho(1, :), dcvdriftdrho(1, :), &
                               dgbdrift0drho(1, :), dgbdriftdrho(1, :), &
                               dgds2dr(1, :), dgds21dr(1, :), &
                               dgds22dr(1, :), djacdrho(1, :))
            !> b_dot_grad_z is the alpha-dependent b . grad z,
            !> and gradpar is the constant-in-alpha part of it.
            !> for axisymmetric systems, b_dot_grad_z is independent of alpha.
            b_dot_grad_z(1, :) = gradpar
            ! note that psi here is the enclosed poloidal flux divided by 2pi
            drhodpsi = 1./dpsidrho
            drhodpsi_psi0 = 1./dpsidrho_psi0
            ! dxdXcoord = a*Bref*dx/dpsi = sign(dx/dpsi) * a*q/r
            dxdXcoord_sign = 1
            sign_torflux = -1
            if (q_as_x) then
               dxdXcoord = dxdXcoord_sign * dpsidrho
            else
               dxdXcoord = dxdXcoord_sign * geo_surf%qinp_psi0 / geo_surf%rhoc_psi0
            end if
            ! dydalpha = (dy/dalpha) / a = sign(dydalpha) * (dpsi/dr) / (a*Bref)
            dydalpha_sign = 1
            dydalpha = dydalpha_sign * dpsidrho
            !> | grad x |  = dx/drho * | grad rho | = dx/dpsi * dpsi/drho * | grad rho |
            !> = q/rho * dpsidrho * grho
            grad_x = grho * dxdXcoord * dpsidrho_psi0
            ! abs(twist_and_shift_geo_fac) is dkx/dky * jtwist
            ! minus its sign gives the direction of the shift in kx
            ! to be used for twist-and-shift BC
            if (q_as_x) then
               twist_and_shift_geo_fac = 2.0 * pi
            else
               twist_and_shift_geo_fac = 2.0 * pi * geo_surf%shat_psi0
            end if
            ! aref and bref should not be needed, so set to 1
            aref = 1.0; bref = 1.0
            zeta(1, :) = zed * geo_surf%qinp

            !!GA - Ensuring all arrays are filled with the alpha = 1 information
            bmag = spread(bmag(1, :), 1, nalpha)
            bmag_psi0 = spread(bmag_psi0(1, :), 1, nalpha)
            gds2 = spread(gds2(1, :), 1, nalpha)
            gds21 = spread(gds21(1, :), 1, nalpha)
            gds22 = spread(gds22(1, :), 1, nalpha)
            gds23 = spread(gds23(1, :), 1, nalpha)
            gds24 = spread(gds24(1, :), 1, nalpha)
            gbdrift0 = spread(gbdrift0(1, :), 1, nalpha)
            gbdrift = spread(gbdrift(1, :), 1, nalpha)
            cvdrift0 = spread(cvdrift0(1, :), 1, nalpha)
            cvdrift = spread(cvdrift(1, :), 1, nalpha)
            dcvdrift0drho = spread(dcvdrift0drho(1, :), 1, nalpha)
            dcvdriftdrho = spread(dcvdriftdrho(1, :), 1, nalpha)
            dgbdrift0drho = spread(dgbdrift0drho(1, :), 1, nalpha)
            dgbdriftdrho = spread(dgbdriftdrho(1, :), 1, nalpha)
            dgds2dr = spread(dgds2dr(1, :), 1, nalpha)
            dgds21dr = spread(dgds21dr(1, :), 1, nalpha)
            dgds22dr = spread(dgds22dr(1, :), 1, nalpha)
            djacdrho = spread(djacdrho(1, :), 1, nalpha)
            b_dot_grad_z = spread(b_dot_grad_z(1, :), 1, nalpha)
            zeta = spread(zeta(1, :), 1, nalpha)

            !!GA - for correct momentum diagnostics
            !>R^2 * grad zeta . grad y / B^2
            gradzeta_grady = geo_surf%rhoc / (geo_surf%qinp * bmag**2)
            gradzeta_gradx = 0.0
            !>R^2 * b . grad zeta = R^2 * b. grad theta * dzeta/dtheta
            gradpar_zeta = geo_surf%rmaj * spread(btor, 1, nalpha) / bmag

         case (geo_option_multibox)
            ! read in Miller local parameters
            call read_local_parameters(nzed, nzgrid, geo_surf)

            ! allocate geometry arrays
            call allocate_arrays(nalpha, nzgrid)

            call communicate_parameters_multibox(surf=geo_surf)

            ! use Miller local parameters to get
            ! geometric coefficients needed by stella
            call get_local_geo(nzed, nzgrid, zed, zed_equal_arc, &
                               dpsidrho, dpsidrho_psi0, dIdrho, grho(1, :), &
                               bmag(1, :), bmag_psi0(1, :), &
                               gds2(1, :), gds21(1, :), gds22(1, :), &
                               gds23(1, :), gds24(1, :), gradpar, &
                               gbdrift0(1, :), gbdrift(1, :), cvdrift0(1, :), cvdrift(1, :), &
                               dBdrho, d2Bdrdth, dgradpardrho, btor, rmajor, &
                               dcvdrift0drho(1, :), dcvdriftdrho(1, :), &
                               dgbdrift0drho(1, :), dgbdriftdrho(1, :), &
                               dgds2dr(1, :), dgds21dr(1, :), &
                               dgds22dr(1, :), djacdrho(1, :))
            !> b_dot_grad_z is the alpha-dependent b . grad z,
            !> and gradpar is the constant-in-alpha part of it.
            !> for axisymmetric systems, b_dot_grad_z is independent of alpha.
            b_dot_grad_z(1, :) = gradpar
            ! note that psi here is the enclosed poloidal flux divided by 2pi
            drhodpsi = 1./dpsidrho
            drhodpsi_psi0 = 1./dpsidrho_psi0
            ! dxdXcoord = a*Bref*dx/dpsi = sign(dx/dpsi) * a*q/r
            dxdXcoord_sign = 1
            sign_torflux = -1
            if (q_as_x) then
               dxdXcoord = dxdXcoord_sign / drhodpsi_psi0
            else
               dxdXcoord = dxdXcoord_sign * geo_surf%qinp_psi0 / geo_surf%rhoc_psi0
            end if
            ! dydalpha = (dy/dalpha) / a = sign(dydalpha) * (dpsi/dr) / (a*Bref)
            dydalpha_sign = 1
            dydalpha = dydalpha_sign / drhodpsi_psi0
            !> | grad x |  = dx/drho * | grad rho | = dx/dpsi * dpsi/drho * | grad rho |
            !> = q/rho * dpsidrho * grho
            grad_x = grho * dxdXcoord * dpsidrho_psi0
            ! abs(twist_and_shift_geo_fac) is dkx/dky * jtwist
            ! minus its sign gives the direction of the shift in kx
            ! to be used for twist-and-shift BC
            if (q_as_x) then
               twist_and_shift_geo_fac = 2.0 * pi
            else
               twist_and_shift_geo_fac = 2.0 * pi * geo_surf%shat_psi0
            end if
            ! aref and bref should not be needed, so set to 1
            aref = 1.0; bref = 1.0
            zeta(1, :) = zed * geo_surf%qinp

            !!GA - Ensuring all arrays are filled with the alpha = 1 information
            bmag = spread(bmag(1, :), 1, nalpha)
            bmag_psi0 = spread(bmag_psi0(1, :), 1, nalpha)
            gds2 = spread(gds2(1, :), 1, nalpha)
            gds21 = spread(gds21(1, :), 1, nalpha)
            gds22 = spread(gds22(1, :), 1, nalpha)
            gds23 = spread(gds23(1, :), 1, nalpha)
            gds24 = spread(gds24(1, :), 1, nalpha)
            gbdrift0 = spread(gbdrift0(1, :), 1, nalpha)
            gbdrift = spread(gbdrift(1, :), 1, nalpha)
            cvdrift0 = spread(cvdrift0(1, :), 1, nalpha)
            cvdrift = spread(cvdrift(1, :), 1, nalpha)
            dcvdrift0drho = spread(dcvdrift0drho(1, :), 1, nalpha)
            dcvdriftdrho = spread(dcvdriftdrho(1, :), 1, nalpha)
            dgbdrift0drho = spread(dgbdrift0drho(1, :), 1, nalpha)
            dgbdriftdrho = spread(dgbdriftdrho(1, :), 1, nalpha)
            dgds2dr = spread(dgds2dr(1, :), 1, nalpha)
            dgds21dr = spread(dgds21dr(1, :), 1, nalpha)
            dgds22dr = spread(dgds22dr(1, :), 1, nalpha)
            djacdrho = spread(djacdrho(1, :), 1, nalpha)
            b_dot_grad_z = spread(b_dot_grad_z(1, :), 1, nalpha)
            zeta = spread(zeta(1, :), 1, nalpha)

            !!GA - for correct momentum diagnostics
            !>R^2 * grad zeta . grad y / B^2
            gradzeta_grady = geo_surf%rhoc / (geo_surf%qinp * bmag**2)
            gradzeta_gradx = 0.0
            !>R^2 * b . grad zeta = R^2 * b. grad theta * dzeta/dtheta
            gradpar_zeta = geo_surf%rmaj * spread(btor, 1, nalpha) / bmag

         case (geo_option_inputprof)
            ! first read in some local parameters
            ! only thing needed really is rhoc
            call read_local_parameters(nzed, nzgrid, geo_surf)
            ! allocate geometry arrays
            call allocate_arrays(nalpha, nzgrid)
            ! now overwrite local parameters
            ! with those from input.profiles file
            ! use rhoc from input as surface
            call read_inputprof_geo(geo_surf)
            call get_local_geo(nzed, nzgrid, zed, zed_equal_arc, &
                               dpsidrho, dpsidrho_psi0, dIdrho, grho(1, :), &
                               bmag(1, :), bmag_psi0(1, :), &
                               gds2(1, :), gds21(1, :), gds22(1, :), &
                               gds23(1, :), gds24(1, :), gradpar, &
                               gbdrift0(1, :), gbdrift(1, :), cvdrift0(1, :), cvdrift(1, :), &
                               dBdrho, d2Bdrdth, dgradpardrho, btor, rmajor, &
                               dcvdrift0drho(1, :), dcvdriftdrho(1, :), &
                               dgbdrift0drho(1, :), dgbdriftdrho(1, :), &
                               dgds2dr(1, :), dgds21dr(1, :), &
                               dgds22dr(1, :), djacdrho(1, :))
            !> b_dot_grad_z is the alpha-dependent b . grad z,
            !> and gradpar is the constant-in-alpha part of it.
            !> for axisymmetric systems, b_dot_grad_z is independent of alpha.
            b_dot_grad_z(1, :) = gradpar
            ! psi here is enclosed poloidal flux divided by 2pi
            drhodpsi = 1./dpsidrho
            drhodpsi_psi0 = 1./dpsidrho_psi0
            ! dxdXcoord = a*Bref*dx/dpsi = sign(dx/dpsi) * a*q/r
            dxdXcoord_sign = 1
            dxdXcoord = dxdXcoord_sign * geo_surf%qinp / geo_surf%rhoc
            sign_torflux = -1
            ! dydalpha = (dy/dalpha) / a = sign(dydalpha) * (dpsi/dr) / (a*Bref)
            dydalpha_sign = 1
            dydalpha = dydalpha_sign * dpsidrho
            !> | grad x |  = dx/drho * | grad rho | = dx/dpsi * dpsi/drho * | grad rho |
            !> = q/rho * dpsidrho * grho
            grad_x = grho * dxdXcoord * dpsidrho_psi0
            ! abs(twist_and_shift_geo_fac) is dkx/dky * jtwist
            ! minus its sign gives the direction of the shift in kx
            ! to be used for twist-and-shift BC
            twist_and_shift_geo_fac = 2.0 * pi * geo_surf%shat
            ! aref and bref should not be needed so set to 1
            aref = 1.0; bref = 1.0

            zeta(1, :) = zed * geo_surf%qinp

            !!GA - Ensuring all arrays are filled with the alpha = 1 information
            bmag = spread(bmag(1, :), 1, nalpha)
            bmag_psi0 = spread(bmag_psi0(1, :), 1, nalpha)
            gds2 = spread(gds2(1, :), 1, nalpha)
            gds21 = spread(gds21(1, :), 1, nalpha)
            gds22 = spread(gds22(1, :), 1, nalpha)
            gds23 = spread(gds23(1, :), 1, nalpha)
            gds24 = spread(gds24(1, :), 1, nalpha)
            gbdrift0 = spread(gbdrift0(1, :), 1, nalpha)
            gbdrift = spread(gbdrift(1, :), 1, nalpha)
            cvdrift0 = spread(cvdrift0(1, :), 1, nalpha)
            cvdrift = spread(cvdrift(1, :), 1, nalpha)
            dcvdrift0drho = spread(dcvdrift0drho(1, :), 1, nalpha)
            dcvdriftdrho = spread(dcvdriftdrho(1, :), 1, nalpha)
            dgbdrift0drho = spread(dgbdrift0drho(1, :), 1, nalpha)
            dgbdriftdrho = spread(dgbdriftdrho(1, :), 1, nalpha)
            dgds2dr = spread(dgds2dr(1, :), 1, nalpha)
            dgds21dr = spread(dgds21dr(1, :), 1, nalpha)
            dgds22dr = spread(dgds22dr(1, :), 1, nalpha)
            djacdrho = spread(djacdrho(1, :), 1, nalpha)
            b_dot_grad_z = spread(b_dot_grad_z(1, :), 1, nalpha)
            zeta = spread(zeta(1, :), 1, nalpha)

            !!GA - for correct momentum diagnostics
            !>R^2 * grad zeta . grad y / B^2
            gradzeta_grady = geo_surf%rhoc / (geo_surf%qinp * bmag**2)
            gradzeta_gradx = 0.0
            !>R^2 * b . grad zeta = R^2 * b. grad theta * dzeta/dtheta
            gradpar_zeta = geo_surf%rmaj * spread(btor, 1, nalpha) / bmag

         case (geo_option_vmec)
            vmec_chosen = .true.
            !> read in input parameters for vmec
            !> nalpha may be specified via input file
            if (debug) write (*, *) 'init_geometry::read_vmec_parameters'
            call read_vmec_parameters
            !> allocate geometry arrays
            if (debug) write (*, *) 'init_geometry::allocate_arrays'
            call allocate_arrays(nalpha, nzgrid)
            if (debug) write (*, *) 'init_geometry::allocate_temporary_arrays'
            call allocate_temporary_arrays(nalpha, nzgrid)
            !> get geometry coefficients from vmec
            if (debug) write (*, *) 'init_geometry::get_vmec_geo'
            !> abs(twist_and_shift_geo_fac) is dkx/dky * jtwist
            !> minus its sign gives the direction of the shift in kx
            !> to be used for twist-and-shift BC
            allocate (twist_and_shift_geo_fac_full(nalpha, -nzgrid:nzgrid))
            twist_and_shift_geo_fac_full = 0
            !> The code will start computing for the initial zgrid
            stellarator_symmetric_BC = .false.
            call get_vmec_geo(new_zeta_min, stellarator_symmetric_BC, &
                              nzgrid, nalpha, naky, geo_surf, grho, bmag, gradpar, &
                              b_dot_grad_z, grad_alpha_grad_alpha, &
                              grad_alpha_grad_psi, grad_psi_grad_psi, &
                              gds23, gds24, gds25, gds26, gbdrift_alpha, gbdrift0_psi, &
                              cvdrift_alpha, cvdrift0_psi, sign_torflux, &
                              theta_vmec, zed_scalefac, aref, bref, alpha, zeta, &
                              field_period_ratio, x_displacement_fac, &
                              gradzeta_grady, gradzeta_gradx, gradpar_zeta)

            if (print_extra_info_to_terminal) then 	
               write (*, '(A)') "############################################################"
               write (*, '(A)') "                     BOUNDARY CONDITIONS"
               write (*, '(A)') "############################################################"
               write (*, *)
            end if

            if ((boundary_option_switch == boundary_option_linked_stellarator) .and. (dkx_over_dky > 0)) then
               ! If stellarator_symmetric_BC and dkx_over_dky>0 the code will compute the new
               ! zeta to work with
               stellarator_symmetric_BC = .true.
               ! twist_and_shift_geo_fac_full for the whole zgrid
               twist_and_shift_geo_fac_full = -2.*(geo_surf%rhotor * geo_surf%rhotor) * (grad_alpha_grad_psi) / (grad_psi_grad_psi)
               ! New min(zeta) to construct the new_number_of_periods_stella
               call desired_zmin(nalpha, nzgrid, zeta, twist_and_shift_geo_fac_full, dkx_over_dky, new_zeta_min)
               ! Final grid of zeta using new_zeta_min
               call get_vmec_geo(new_zeta_min, stellarator_symmetric_BC, &
                                 nzgrid, nalpha, naky, geo_surf, grho, bmag, gradpar, &
                                 b_dot_grad_z, grad_alpha_grad_alpha, &
                                 grad_alpha_grad_psi, grad_psi_grad_psi, &
                                 gds23, gds24, gds25, gds26, gbdrift_alpha, gbdrift0_psi, &
                                 cvdrift_alpha, cvdrift0_psi, sign_torflux, &
                                 theta_vmec, zed_scalefac, aref, bref, alpha, zeta, &
                                 field_period_ratio, x_displacement_fac, &
                                 gradzeta_grady, gradzeta_gradx, gradpar_zeta)
               ! Restart the variable twist_and_shift_geo_fac_full
               twist_and_shift_geo_fac_full = 0
            end if

            !> Define <dxdXcoord> = (ρref/a)(dx̃/dψ̃) and use dx/dψ = 1/(a*ρ0*Bref)
            !> dxdXcoord = (ρref/a)(dx̃/dψ̃) = (ρref/a)(d(x/ρref)/d(ψ/(a^2 Bref)) = a Bref (dx/dψ) = 1/ρ0
            dxdXcoord = 1 / geo_surf%rhotor

            !> Define <dydalpha> = (ρref/a)(dỹ/dα̃) and use dy/dα = a*ρ0
            !> dydalpha = (ρref/a)(dỹ/dα̃) = (ρref/a)(d(y/ρref)/dα) = (1/a) (dy/dα) = ρ0
            dydalpha = geo_surf%rhotor

            !> Define <drhodpsi> = dρ0/dψ̃ and use ρ0 = sqrt(psi_t/psi_{t,LCFS}) and ψ = sgn(psi_t)*psi_t
            !> Use dρ0/dpsi_t = d(sqrt(psi_t/psi_{t,LCFS}))/dpsi_t = sgn(psi_t)/(a^2 ρ0 Bref) with Bref = 2 |psi_{t,LCFS}| / a^2
            !> drhodpsi = dρ0/dψ̃ = dρ0/d(ψ/(a^2 Bref)) = a^2 Bref sgn(psi_t) dρ0/dpsi_t = 1/ρ0
            drhodpsi = 1 / geo_surf%rhotor
            drhodpsi_psi0 = drhodpsi

            bmag_psi0 = bmag
            grad_x_grad_y_end = grad_alpha_grad_psi(1, nzgrid) * (aref * aref * bref)
            select case (boundary_option_switch)
            case (boundary_option_linked_stellarator)
               !to be used for stellarator symmetric twist-and-shift BC
               !twist_and_shift_geo_fac = -nabla x. nabla y /|nabla x|^2
               if (print_extra_info_to_terminal) write (*, *) 'Stellarator symmetric twist and shift BC selected'
               twist_and_shift_geo_fac_full = -2.*(geo_surf%rhotor * geo_surf%rhotor) * (grad_alpha_grad_psi) / (grad_psi_grad_psi)
               if ((abs(grad_x_grad_y_end) <= grad_x_grad_y_zero) .and. (print_extra_info_to_terminal)) write (*, *) 'Using periodic boundary conditions as grad_x_grad_y_end < grad_x_grad_y_zero'
               if (print_extra_info_to_terminal) write (*, *)
            case default
               ! to be used for std twist-and-shift BC
               ! twist_and_shift_geo_fac =
               ! -2.*pi*geo_surf%shat*geo_surf%qinp*drhodpsi*dydalpha/(dxdpsi*geo_surf%rhotor)
               if (print_extra_info_to_terminal) write (*, *) 'Standard twist and shift BC selected'
               twist_and_shift_geo_fac_full = -2.*pi * geo_surf%shat * drhodpsi * dydalpha / (geo_surf%qinp * dxdXcoord * geo_surf%rhotor) &
                                              * field_period_ratio
               if ((abs(geo_surf%shat) <= shat_zero) .and. (print_extra_info_to_terminal)) write (*, *) 'Using periodic boundary conditions as shat < shat_zero'
               if (print_extra_info_to_terminal) write (*, *)
            end select
            twist_and_shift_geo_fac = twist_and_shift_geo_fac_full(1, nzgrid)
            deallocate (twist_and_shift_geo_fac_full)

            !> grad_x = | grad x |
            grad_x = sqrt(abs(grad_psi_grad_psi * dxdXcoord**2))

            !> gds2 = |grad y|^2 = |grad alpha|^2 * (dy/dalpha)^2
            !> note that rhotor = sqrt(psi/psi_LCFS)
            gds2 = grad_alpha_grad_alpha * dydalpha**2

            !> Define <gds21> = hat{s} ∇x . ∇y
            !> Use (dx/dψ)*(dy/dα) = 1/(a ρ0 Bref) * (a ρ0) = 1/Bref
            !> Use ∇x . ∇y = (dx/dψ)(dy/dα) ∇ψ . ∇α = (1/Bref) ∇ψ . ∇α = <grad_alpha_grad_psi>
            gds21 = grad_alpha_grad_psi * geo_surf%shat

            !> gds22 = shat^2 * |grad x|^2 = shat^2 * |grad psi_t|^2 * (dx/dpsi_t)^2
            gds22 = (geo_surf%shat * grad_x)**2
            !gds22 = geo_surf%shat**2 * grad_psi_grad_psi * dxdXcoord**2

            !> gbdrift_alpha and cvdrift_alpha contain
            !> the grad-B and curvature drifts projected onto
            !> the grad alpha direction
            !> need the projections on grad y
            gbdrift = gbdrift_alpha * dydalpha
            cvdrift = cvdrift_alpha * dydalpha

            !> gbdrift0_psi and cvdrift0_psi contain
            !> the grad-B and curvature drifts projected onto
            !> the grad psi direction
            !> need the projections on grad x
            gbdrift0 = gbdrift0_psi * dxdXcoord
            cvdrift0 = cvdrift0_psi * dxdXcoord

            call deallocate_temporary_arrays

            !> can test FFS implementation by setting all geometric coefficients
            !> to their values at a given alpha; i.e., make the system axisymmetric
            if (const_alpha_geo) call set_ffs_geo_coefs_constant(nalpha)
         end select

         if (overwrite_geometry) call overwrite_selected_geometric_coefficients(nalpha)

         ! exb_nonlin_fac is equivalent to kxfac/2 in gs2
         if (q_as_x) then
            dqdrho = geo_surf%shat * geo_surf%qinp / geo_surf%rhoc
            exb_nonlin_fac = -0.5 * sign_torflux * dxdXcoord * dydalpha * drhodpsi * dqdrho
            !the following will get multiplied by exb_nonlin_fac in advance_exb_nonlinearity
            exb_nonlin_fac_p = geo_surf%d2qdr2 / dqdrho - geo_surf%d2psidr2 * drhodpsi
         else
            exb_nonlin_fac = -0.5 * sign_torflux * dxdXcoord * dydalpha
            exb_nonlin_fac_p = 0.0
         end if
      end if

      if (.not. proc0) call allocate_arrays(nalpha, nzgrid)

      if (debug .and. proc0) write (*, *) 'init_geometry::broadcast_arrays'
      call broadcast_arrays

      gfac = 1.0
      if (.not. include_geometric_variation) gfac = 0.0

      ! should reduce to 2*pi*shat in axisymmetric case
      ! but not in non-axisymmetric case
!    twist_and_shift_geo_fac = geo_surf%shat*(gds21(1,-nzgrid)/gds22(1,-nzgrid)-gds21(1,nzgrid)/gds22(1,nzgrid))

      ! FLAG DSO - the followiing assumes a linear relation from q to rho, but this will
      !            not be correct if d2qdrho != 0
      dqdrho = geo_surf%shat * geo_surf%qinp / geo_surf%rhoc

      ! this old definition of jacob should have been fine, as it was only ever used in both the numerator
      ! and denominator of averages -- and do any constant factors cancelled out
      !jacob = 1.0/abs(drhodpsi*spread(gradpar,1,nalpha)*bmag)

      !> jacob is the Jacobian from Cartesian coordinates to (y,x,z) coordinates
      !> is ((grad y x grad x) . grad z)^(-1) = Lref*(dalpha/dy)*(dpsi/dx)/(Lref*Bref)*(B/Bref . grad z)^(-1)
      !> Lref*(dalpha/dy) = 1/dydalpha; (dpsi/dx)/(Lref*Bref) = 1 / dxdXcoord ; (B/Bref . grad z) = gradpar*bmag
      jacob = -sign_torflux / (dydalpha * dxdXcoord * b_dot_grad_z * bmag)
      !    jacob = 1.0/(dydalpha*dxdXcoord*spread(gradpar,1,nalpha)*bmag)

      ! this is dl/B
      dl_over_b = spread(delzed, 1, nalpha) * jacob

      ! the next line is to avoid double counting the end points for ky = 0 modes (which leads to destabilization
      ! of the zonal modes for certain input parameters)
      ! FLAG DSO - while this is correct for ky = 0 modes and sufficient for output, if dl_over_b is applied to
      ! non-zero ky modes, a more sophisticated approach will be required that takes into account the sign of
      ! v_parallel
      dl_over_b(:, nzgrid) = 0.

      ! this is the correction to flux-surface-averaging for adiabatic electrons
      d_dl_over_b_drho = spread(delzed, 1, nalpha) * djacdrho
      d_dl_over_b_drho(:, nzgrid) = 0
      d_dl_over_b_drho = d_dl_over_b_drho - dl_over_b &
                         * spread(sum(d_dl_over_b_drho, dim=2) / sum(dl_over_b, dim=2), 2, 2 * nzgrid + 1)
      d_dl_over_b_drho = gfac * d_dl_over_b_drho / spread(sum(dl_over_b, dim=2), 2, 2 * nzgrid + 1)

      ! normalize dl/B by int dl/B
      dl_over_b = dl_over_b / spread(sum(dl_over_b, dim=2), 2, 2 * nzgrid + 1)

      ! this is what we use to normalize the fluxes
      grho_norm = sum(dl_over_b(1, :) * grho(1, :))

      ! would probably be better to compute this in the various
      ! geometry subroutine (Miller, vmec, etc.), as there
      ! B is likely calculated on a finer z-grid
      do iy = 1, nalpha
         call get_dzed(nzgrid, delzed, bmag(iy, :), dbdzed(iy, :))
      end do

      select case (boundary_option_switch)
      case (boundary_option_linked)
         ! if magnetic shear almost zero, override parallel
         ! boundary condition so that it is periodic if using the standard
         ! twist and shift bc, in which kx_shift is proportional to shat
         if (abs(geo_surf%shat) <= shat_zero) &
            boundary_option_switch = boundary_option_self_periodic
      case (boundary_option_linked_stellarator)
         ! if magnetic nabla x. nabla y is almost zero, override parallel
         ! boundary condition so that it is periodic if using the stellarator
         ! symmetric twist and shift bc, in which kx_shift is proportional to nabla
         ! x. nabla y
         if (abs(grad_x_grad_y_end) <= grad_x_grad_y_zero) &
            boundary_option_switch = boundary_option_self_periodic
      case default
         ! low shear is fine for periodic/zero parallel BCs, so do nothing
      end select

      ! theta_eqarc is parallel coordinate such that
      ! b . grad theta_eqarc = constant
      ! and theta_eqarc = theta at +/- pi
      ! b . grad theta_eqarc = b . grad theta dtheta_eqarc/dtheta
      ! --> dtheta_eqarc/dtheta = b . grad theta_eqarc / b . grad theta
      ! --> 2*pi = b . grad theta_eqarc * int_0^{2pi} dtheta 1/(b.grad theta)
      ! this gives b . grad theta_eqarc, from which we get
      ! theta_eqarc = theta_min + int_{0}^{theta} dtheta' b . grad theta_eqarc / b . grad theta'
      call get_gradpar_eqarc(gradpar, zed, delzed, gradpar_eqarc)
      call get_zed_eqarc(gradpar, delzed, zed, gradpar_eqarc, zed_eqarc)

      if (proc0) call write_geometric_coefficients(nalpha)

      ! AVB: temporary, set bmag = constant in z for Spitzer problem
      if (set_bmag_const) then
         bmag_z0 = bmag(1, 0)
         print *, ''
         print *, '! SETTING BMAG = CONSTANT IN Z'
         print *, ''
         do ia = 1, nalpha
            do iz = -nzgrid, nzgrid
               bmag(ia, iz) = bmag_z0
            end do
         end do
      end if

   contains

      subroutine allocate_temporary_arrays(nalpha, nzgrid)

         implicit none

         integer, intent(in) :: nalpha, nzgrid

         allocate (grad_alpha_grad_alpha(nalpha, -nzgrid:nzgrid))
         allocate (grad_alpha_grad_psi(nalpha, -nzgrid:nzgrid))
         allocate (grad_psi_grad_psi(nalpha, -nzgrid:nzgrid))
         allocate (gbdrift_alpha(nalpha, -nzgrid:nzgrid))
         allocate (cvdrift_alpha(nalpha, -nzgrid:nzgrid))
         allocate (gbdrift0_psi(nalpha, -nzgrid:nzgrid))
         allocate (cvdrift0_psi(nalpha, -nzgrid:nzgrid))

      end subroutine allocate_temporary_arrays

      subroutine deallocate_temporary_arrays

         implicit none

         deallocate (grad_alpha_grad_alpha)
         deallocate (grad_alpha_grad_psi)
         deallocate (grad_psi_grad_psi)
         deallocate (gbdrift_alpha)
         deallocate (cvdrift_alpha)
         deallocate (gbdrift0_psi)
         deallocate (cvdrift0_psi)

      end subroutine deallocate_temporary_arrays

      subroutine overwrite_selected_geometric_coefficients(nalpha)

         use file_utils, only: get_unused_unit
         use zgrid, only: nzgrid

         implicit none

         integer, intent(in) :: nalpha
         integer :: geofile_unit
         character(100) :: dum_char
         real :: dum_real

         integer :: ia, iz
         real :: bmag_file, gradpar_file
         real :: gds2_file, gds21_file, gds22_file, gds23_file, gds24_file
         real :: gbdrift_file, cvdrift_file, gbdrift0_file

         call get_unused_unit(geofile_unit)
         open (geofile_unit, file=trim(geo_file), status='old', action='read')

         read (geofile_unit, fmt=*) dum_char
         read (geofile_unit, fmt=*) dum_char
         read (geofile_unit, fmt=*) dum_char

         ! overwrite bmag, gradpar, gds2, gds21, gds22, gds23, gds24, gbdrift, cvdrift, gbdrift0, and cvdrift0
         ! with values from file
         do ia = 1, nalpha
            do iz = -nzgrid, nzgrid
               read (geofile_unit, fmt='(13e12.4)') dum_real, dum_real, dum_real, bmag_file, gradpar_file, &
                  gds2_file, gds21_file, gds22_file, gds23_file, &
                  gds24_file, gbdrift_file, cvdrift_file, gbdrift0_file
               if (overwrite_bmag) bmag(ia, iz) = bmag_file
               if (overwrite_gradpar) gradpar(iz) = gradpar_file
               if (overwrite_gradpar) b_dot_grad_z(1, iz) = gradpar_file ! assuming we are only reading in for a single alpha. Usually, gradpar is the average of all b_dot_grad_z values.
               if (overwrite_gds2) gds2(ia, iz) = gds2_file
               if (overwrite_gds21) gds21(ia, iz) = gds21_file
               if (overwrite_gds22) gds22(ia, iz) = gds22_file
               if (overwrite_gds23) gds23(ia, iz) = gds23_file
               if (overwrite_gds24) gds24(ia, iz) = gds24_file
               if (overwrite_gbdrift) gbdrift(ia, iz) = gbdrift_file
               if (overwrite_cvdrift) cvdrift(ia, iz) = cvdrift_file
               if (overwrite_gbdrift0) gbdrift0(ia, iz) = gbdrift0_file
            end do
         end do
         cvdrift0 = gbdrift0

         close (geofile_unit)

      end subroutine overwrite_selected_geometric_coefficients

      subroutine set_ffs_geo_coefs_constant(nalpha)

         implicit none

         integer, intent(in) :: nalpha

         call set_coef_constant(gbdrift0, nalpha)
         call set_coef_constant(cvdrift0, nalpha)
         call set_coef_constant(gbdrift, nalpha)
         call set_coef_constant(cvdrift, nalpha)
         call set_coef_constant(grad_x, nalpha)
         call set_coef_constant(grho, nalpha)
         call set_coef_constant(bmag, nalpha)
         call set_coef_constant(bmag_psi0, nalpha)
         call set_coef_constant(gds2, nalpha)
         call set_coef_constant(gds21, nalpha)
         call set_coef_constant(gds22, nalpha)
         call set_coef_constant(gds23, nalpha)
         call set_coef_constant(gds24, nalpha)
         call set_coef_constant(gds25, nalpha)
         call set_coef_constant(gds26, nalpha)
         call set_coef_constant(theta_vmec, nalpha)
         call set_coef_constant(x_displacement_fac, nalpha)
         call set_coef_constant(zeta, nalpha)
         call set_coef_constant(b_dot_grad_z, nalpha)
         !!GA
         call set_coef_constant(gradzeta_grady, nalpha)
         call set_coef_constant(gradzeta_gradx, nalpha)
         call set_coef_constant(gradpar_zeta, nalpha)
         ! following coefficients calculated later using above coefficients
         !      call set_coef_constant (dbdzed, nalpha)
         !      call set_coef_constant (jacob, nalpha)
         !      call set_coef_constant (dl_over_b, nalpha)

      end subroutine set_ffs_geo_coefs_constant

      subroutine set_coef_constant(coef, nalpha)

         implicit none

         real, dimension(:, -nzgrid:), intent(in out) :: coef
         integer, intent(in) :: nalpha

         coef = spread(coef(1, :), 1, nalpha)

      end subroutine set_coef_constant

   end subroutine init_geometry

   subroutine allocate_arrays(nalpha, nzgrid)

      implicit none

      integer, intent(in) :: nalpha, nzgrid

      if (.not. allocated(bmag)) allocate (bmag(nalpha, -nzgrid:nzgrid)); bmag = 0.0 
      if (.not. allocated(bmag_psi0)) allocate (bmag_psi0(nalpha, -nzgrid:nzgrid)); bmag_psi0 = 0.0
      if (.not. allocated(gds2)) allocate (gds2(nalpha, -nzgrid:nzgrid)); gds2 = 0.0
      if (.not. allocated(gds21)) allocate (gds21(nalpha, -nzgrid:nzgrid)); gds21 = 0.0
      if (.not. allocated(gds22)) allocate (gds22(nalpha, -nzgrid:nzgrid)); gds22 = 0.0
      if (.not. allocated(gds23)) allocate (gds23(nalpha, -nzgrid:nzgrid)); gds23 = 0.0
      if (.not. allocated(gds24)) allocate (gds24(nalpha, -nzgrid:nzgrid)); gds24 = 0.0
      if (.not. allocated(gds25)) allocate (gds25(nalpha, -nzgrid:nzgrid)); gds25 = 0.0
      if (.not. allocated(gds26)) allocate (gds26(nalpha, -nzgrid:nzgrid)); gds26 = 0.0
      if (.not. allocated(dgds2dr)) allocate (dgds2dr(nalpha, -nzgrid:nzgrid)); dgds2dr = 0.0
      if (.not. allocated(dgds21dr)) allocate (dgds21dr(nalpha, -nzgrid:nzgrid)); dgds21dr = 0.0
      if (.not. allocated(dgds22dr)) allocate (dgds22dr(nalpha, -nzgrid:nzgrid)); dgds22dr = 0.0
      if (.not. allocated(gbdrift)) allocate (gbdrift(nalpha, -nzgrid:nzgrid)); gbdrift = 0.0
      if (.not. allocated(gbdrift0)) allocate (gbdrift0(nalpha, -nzgrid:nzgrid)); gbdrift0 = 0.0
      if (.not. allocated(cvdrift)) allocate (cvdrift(nalpha, -nzgrid:nzgrid)); cvdrift = 0.0
      if (.not. allocated(cvdrift0)) allocate (cvdrift0(nalpha, -nzgrid:nzgrid)); cvdrift0 = 0.0
      if (.not. allocated(dgbdriftdrho)) allocate (dgbdriftdrho(nalpha, -nzgrid:nzgrid)); dgbdriftdrho = 0.0
      if (.not. allocated(dcvdriftdrho)) allocate (dcvdriftdrho(nalpha, -nzgrid:nzgrid)); dcvdriftdrho = 0.0
      if (.not. allocated(dgbdrift0drho)) allocate (dgbdrift0drho(nalpha, -nzgrid:nzgrid)); dgbdrift0drho = 0.0
      if (.not. allocated(dcvdrift0drho)) allocate (dcvdrift0drho(nalpha, -nzgrid:nzgrid)); dcvdrift0drho = 0.0
      if (.not. allocated(dbdzed)) allocate (dbdzed(nalpha, -nzgrid:nzgrid)); dbdzed = 0.0
      if (.not. allocated(theta_vmec)) allocate (theta_vmec(nalpha, -nzgrid:nzgrid)); theta_vmec = 0.0
      if (.not. allocated(jacob)) allocate (jacob(nalpha, -nzgrid:nzgrid)); jacob = 0.0
      if (.not. allocated(djacdrho)) allocate (djacdrho(nalpha, -nzgrid:nzgrid)); djacdrho = 0.0
      if (.not. allocated(grho)) allocate (grho(nalpha, -nzgrid:nzgrid)); grho = 0.0
      if (.not. allocated(grad_x)) allocate (grad_x(nalpha, -nzgrid:nzgrid)); grad_x = 0.0
      if (.not. allocated(dl_over_b)) allocate (dl_over_b(nalpha, -nzgrid:nzgrid)); dl_over_b = 0.0
      if (.not. allocated(d_dl_over_b_drho)) allocate (d_dl_over_b_drho(nalpha, -nzgrid:nzgrid)); d_dl_over_b_drho = 0.0
      if (.not. allocated(b_dot_grad_z)) allocate (b_dot_grad_z(nalpha, -nzgrid:nzgrid)); b_dot_grad_z = 0.0

      if (.not. allocated(gradpar)) allocate (gradpar(-nzgrid:nzgrid)); gradpar = 0.0
      if (.not. allocated(zed_eqarc)) allocate (zed_eqarc(-nzgrid:nzgrid)); zed_eqarc = 0.0
      if (.not. allocated(btor)) allocate (btor(-nzgrid:nzgrid)); btor = 0.0
      if (.not. allocated(rmajor)) allocate (rmajor(-nzgrid:nzgrid)); rmajor = 0.0
      if (.not. allocated(dBdrho)) allocate (dBdrho(-nzgrid:nzgrid)); dBdrho = 0.0
      if (.not. allocated(d2Bdrdth)) allocate (d2Bdrdth(-nzgrid:nzgrid)); d2Bdrdth = 0.0
      if (.not. allocated(dgradpardrho)) allocate (dgradpardrho(-nzgrid:nzgrid)); dgradpardrho = 0.0

      if (.not. allocated(alpha)) allocate (alpha(nalpha)); alpha = 0.0
      if (.not. allocated(zeta)) allocate (zeta(nalpha, -nzgrid:nzgrid)); zeta = 0.0

      if (.not. allocated(x_displacement_fac)) allocate (x_displacement_fac(nalpha, -nzgrid:nzgrid)); x_displacement_fac = 0.0

      !!GA - Needed for the momentum flux diagnostic for non-axisymmetric devices
      if (.not. allocated(gradzeta_grady)) allocate (gradzeta_grady(nalpha, -nzgrid:nzgrid)); gradzeta_grady = 0.0
      if (.not. allocated(gradzeta_gradx)) allocate (gradzeta_gradx(nalpha, -nzgrid:nzgrid)); gradzeta_gradx = 0.0
      if (.not. allocated(gradpar_zeta)) allocate (gradpar_zeta(nalpha, -nzgrid:nzgrid)); gradpar_zeta = 0.0

   end subroutine allocate_arrays

   subroutine read_parameters

      use text_options
      use file_utils, only: error_unit, input_unit_exist
      use file_utils, only: runtype_option_Switch, runtype_multibox
      use mp, only: job
      use physics_flags, only: radial_variation

      implicit none

      integer :: in_file, ierr
      logical :: exist

      character(20) :: geo_option
      type(text_option), dimension(5), parameter :: geoopts = (/ &
                                                    text_option('default', geo_option_local), &
                                                    text_option('miller', geo_option_local), &
                                                    text_option('local', geo_option_local), &
                                                    text_option('input.profiles', geo_option_inputprof), &
                                                    text_option('vmec', geo_option_vmec)/)

      namelist /geo_knobs/ geo_option, geo_file, overwrite_bmag, overwrite_gradpar, &
         overwrite_gds2, overwrite_gds21, overwrite_gds22, overwrite_gds23, overwrite_gds24, &
         overwrite_gbdrift, overwrite_cvdrift, overwrite_gbdrift0, q_as_x, set_bmag_const

      geo_option = 'local'
      overwrite_bmag = .false.
      overwrite_gradpar = .false.
      overwrite_gds2 = .false.
      overwrite_gds21 = .false.
      overwrite_gds22 = .false.
      overwrite_gds23 = .false.
      overwrite_gds24 = .false.
      overwrite_gbdrift = .false.
      overwrite_cvdrift = .false.
      overwrite_gbdrift0 = .false.
      set_bmag_const = .false.
      q_as_x = radial_variation !true by default in radial variation runs
      geo_file = 'input.geometry'

      in_file = input_unit_exist("geo_knobs", exist)
      if (exist) read (unit=in_file, nml=geo_knobs)

      ierr = error_unit()
      call get_option_value &
         (geo_option, geoopts, geo_option_switch, &
          ierr, "geo_option in geo_knobs")

      if (radial_variation .and. runtype_option_switch == runtype_multibox .and. job /= 1) then
         geo_option_switch = geo_option_multibox
      end if

      overwrite_geometry = overwrite_bmag .or. overwrite_gradpar &
                           .or. overwrite_gds2 .or. overwrite_gds21 .or. overwrite_gds22 &
                           .or. overwrite_gds23 .or. overwrite_gds24 &
                           .or. overwrite_cvdrift .or. overwrite_gbdrift .or. overwrite_gbdrift0

   end subroutine read_parameters

   subroutine broadcast_arrays

      use mp, only: broadcast

      implicit none

      call broadcast(qinp)
      call broadcast(shat)
      call broadcast(drhodpsi)
      call broadcast(drhodpsi_psi0)
      call broadcast(exb_nonlin_fac)
      call broadcast(exb_nonlin_fac_p)
      call broadcast(dIdrho)
      call broadcast(grho)
      call broadcast(grad_x)
      call broadcast(bmag)
      call broadcast(bmag_psi0)
      call broadcast(btor)
      call broadcast(rmajor)
      call broadcast(gradpar)
      call broadcast(b_dot_grad_z)
      call broadcast(gds2)
      call broadcast(gds21)
      call broadcast(gds22)
      call broadcast(gds23)
      call broadcast(gds24)
      call broadcast(gds25)
      call broadcast(gds26)
      call broadcast(dgds2dr)
      call broadcast(dgds21dr)
      call broadcast(dgds22dr)
      call broadcast(gbdrift0)
      call broadcast(gbdrift)
      call broadcast(cvdrift0)
      call broadcast(cvdrift)
      call broadcast(dgbdrift0drho)
      call broadcast(dgbdriftdrho)
      call broadcast(dcvdrift0drho)
      call broadcast(dcvdriftdrho)
      call broadcast(dBdrho)
      call broadcast(d2Bdrdth)
      call broadcast(dgradpardrho)
      call broadcast(djacdrho)

      call broadcast(geo_surf%rmaj)
      call broadcast(geo_surf%rgeo)
      call broadcast(geo_surf%kappa)
      call broadcast(geo_surf%kapprim)
      call broadcast(geo_surf%tri)
      call broadcast(geo_surf%triprim)
      call broadcast(geo_surf%rhoc)
      call broadcast(geo_surf%rhoc_psi0)
      call broadcast(geo_surf%dr)
      call broadcast(geo_surf%shift)
      call broadcast(geo_surf%qinp)
      call broadcast(geo_surf%qinp_psi0)
      call broadcast(geo_surf%shat)
      call broadcast(geo_surf%shat_psi0)
      call broadcast(geo_surf%betaprim)
      call broadcast(geo_surf%betadbprim)
      call broadcast(geo_surf%d2qdr2)
      call broadcast(geo_surf%d2psidr2)
      call broadcast(geo_surf%dpsitordrho)
      call broadcast(geo_surf%rhotor)
      call broadcast(geo_surf%psitor_lcfs)
      call broadcast(geo_surf%drhotordrho)

      call broadcast(zed_scalefac)
      call broadcast(q_as_x)
      call broadcast(set_bmag_const)
      call broadcast(alpha)
      call broadcast(zeta)
      call broadcast(dxdXcoord)
      call broadcast(dydalpha)
      call broadcast(sign_torflux)
      call broadcast(twist_and_shift_geo_fac)
      call broadcast(grad_x_grad_y_end)
      call broadcast(vmec_chosen)

      call broadcast(aref)
      call broadcast(bref)

      !!GA
      call broadcast(gradzeta_grady)
      call broadcast(gradzeta_gradx)
      call broadcast(gradpar_zeta)

   end subroutine broadcast_arrays

   subroutine communicate_geo_multibox(l_edge, r_edge)

      use millerlocal, only: communicate_parameters_multibox
      use mp, only: proc0

      implicit none

      real, intent(in) :: l_edge, r_edge

      if (proc0) then
         call communicate_parameters_multibox(geo_surf, gfac * l_edge, gfac * r_edge)
      end if

   end subroutine communicate_geo_multibox

   ! given function f(z:-pi->pi), calculate z derivative
   ! second order accurate, with equal grid spacing assumed
   ! assumes periodic in z -- may need to change this in future
   subroutine get_dzed(nz, dz, f, df)

      implicit none

      integer, intent(in) :: nz
      real, dimension(-nz:), intent(in) :: dz, f
      real, dimension(-nz:), intent(out) :: df

      df(-nz + 1:nz - 1) = (f(-nz + 2:) - f(:nz - 2)) / (dz(:nz - 2) + dz(-nz + 1:nz - 1))

      ! hack to avoid non-periodicity in full-flux-surface case
!      if (full_flux_surface .and. .not. const_alpha_geo) then
 !        df(-nz) = (f(-nz + 1) - f(-nz)) / dz(-nz)
   !      df(nz) = (f(nz) - f(nz - 1)) / dz(nz - 1)
  !    else
         ! assume periodicity in the B-field
      df(-nz) = (f(-nz + 1) - f(nz - 1)) / (dz(-nz) + dz(nz - 1))
      df(nz) = df(-nz)
   !end if

   end subroutine get_dzed

   subroutine get_gradpar_eqarc(gp, z, dz, gp_eqarc)

      use constants, only: pi
      use zgrid, only: nzgrid

      implicit none

      real, dimension(-nzgrid:), intent(in) :: gp, z, dz
      real, intent(out) :: gp_eqarc

      ! first get int dz b . grad z
      call integrate_zed(dz, 1./gp, gp_eqarc)
      ! then take (zmax-zmin)/int (dz b . gradz)
      ! to get b . grad z'
      gp_eqarc = (z(nzgrid) - z(-nzgrid)) / gp_eqarc

   end subroutine get_gradpar_eqarc

   subroutine get_zed_eqarc(gp, dz, z, gp_eqarc, z_eqarc)

      use zgrid, only: nzgrid

      implicit none

      real, dimension(-nzgrid:), intent(in) :: gp, dz, z
      real, intent(in) :: gp_eqarc
      real, dimension(-nzgrid:), intent(out) :: z_eqarc

      integer :: iz

      z_eqarc(-nzgrid) = z(-nzgrid)
      do iz = -nzgrid + 1, nzgrid
         call integrate_zed(dz(:iz), 1./gp(:iz), z_eqarc(iz))
      end do
      z_eqarc(-nzgrid + 1:) = z(-nzgrid) + z_eqarc(-nzgrid + 1:) * gp_eqarc

   end subroutine get_zed_eqarc

   ! trapezoidal rule to integrate in zed
   subroutine integrate_zed(dz, f, intf)

      use zgrid, only: nzgrid

      implicit none

      real, dimension(-nzgrid:), intent(in) :: dz
      real, dimension(-nzgrid:), intent(in) :: f
      real, intent(out) :: intf

      integer :: iz, iz_max

      iz_max = -nzgrid + size(dz) - 1

      intf = 0.
      do iz = -nzgrid + 1, iz_max
         intf = intf + dz(iz) * (f(iz - 1) + f(iz))
      end do
      intf = 0.5 * intf

   end subroutine integrate_zed

   subroutine get_x_to_rho(llim, x_in, rho_out)

      use physics_parameters, only: rhostar

      implicit none

      integer, intent(in) :: llim

      real, dimension(:), intent(in) :: x_in
      real, dimension(:), intent(out) :: rho_out

      integer :: ix, ulim

      real :: a, b, c

      ulim = size(x_in) + llim - 1

      if (q_as_x) then
         a = 0.5 * geo_surf%d2qdr2 / dqdrho
         b = 1.0
         c = -rhostar / (dqdrho * dxdXcoord)
      else
         a = 0.5 * geo_surf%d2psidr2 * drhodpsi
         b = 1.0
         c = -rhostar * drhodpsi / dxdXcoord
      end if

      do ix = llim, ulim
         if (abs(4.0 * a * c * x_in(ix)) < 1.e-6) then
            rho_out(ix) = -(c * x_in(ix)) / b - a * (c * x_in(ix))**2 / b**3
         else
            rho_out(ix) = (-b + sqrt(b**2 - 4.*a * c * x_in(ix))) / (2.*a)
         end if
      end do

   end subroutine get_x_to_rho

   subroutine write_geometric_coefficients(nalpha)

      use file_utils, only: open_output_file, close_output_file
      use zgrid, only: nzgrid, zed

      implicit none

      integer, intent(in) :: nalpha
      integer :: geometry_unit
      integer :: ia, iz

      call open_output_file(geometry_unit, '.geometry')

      write (geometry_unit, '(a1,11a14)') '#', 'rhoc', 'qinp', 'shat', 'rhotor', 'aref', 'bref', 'dxdXcoord', 'dydalpha', &
         'exb_nonlin', 'exb_nonlin_p'
      write (geometry_unit, '(a1,11e14.4)') '#', geo_surf%rhoc, geo_surf%qinp, geo_surf%shat, geo_surf%rhotor, aref, bref, &
         dxdXcoord, dydalpha, exb_nonlin_fac, exb_nonlin_fac_p * exb_nonlin_fac
      write (geometry_unit, *)

      write (geometry_unit, '(14a12)') '# alpha', 'zed', 'zeta', 'bmag', 'bdot_grad_z', 'gds2', 'gds21', 'gds22', &
         'gds23', 'gds24', 'gbdrift', 'cvdrift', 'gbdrift0', 'bmag_psi0'
      do ia = 1, nalpha
         do iz = -nzgrid, nzgrid
            !          write (geometry_unit,'(15e12.4)') alpha(ia), zed(iz), zeta(ia,iz), bmag(ia,iz), gradpar(iz), &
            write (geometry_unit, '(15e12.4)') alpha(ia), zed(iz), zeta(ia, iz), bmag(ia, iz), b_dot_grad_z(ia, iz), &
               gds2(ia, iz), gds21(ia, iz), gds22(ia, iz), gds23(ia, iz), &
               gds24(ia, iz), gbdrift(ia, iz), cvdrift(ia, iz), gbdrift0(ia, iz), &
               bmag_psi0(ia, iz)
         end do
         write (geometry_unit, *)
      end do

      call close_output_file(geometry_unit)

   end subroutine write_geometric_coefficients

   subroutine finish_init_geometry

      use mp, only: proc0
      use millerlocal, only: finish_local_geo

      implicit none

      if (proc0) then
         select case (geo_option_switch)
         case (geo_option_local)
            call finish_local_geo
         case (geo_option_multibox)
            call finish_local_geo
         case (geo_option_inputprof)
            call finish_local_geo
         case (geo_option_vmec)
         end select
      end if

   end subroutine finish_init_geometry

   subroutine finish_geometry

      implicit none

      if (allocated(zed_eqarc)) deallocate (zed_eqarc)
      if (allocated(grho)) deallocate (grho)
      if (allocated(grad_x)) deallocate (grad_x)
      if (allocated(bmag)) deallocate (bmag)
      if (allocated(bmag_psi0)) deallocate (bmag_psi0)
      if (allocated(btor)) deallocate (btor)
      if (allocated(rmajor)) deallocate (rmajor)
      if (allocated(dbdzed)) deallocate (dbdzed)
      if (allocated(jacob)) deallocate (jacob)
      if (allocated(djacdrho)) deallocate (djacdrho)
      if (allocated(gradpar)) deallocate (gradpar)
      if (allocated(b_dot_grad_z)) deallocate (b_dot_grad_z)
      if (allocated(dl_over_b)) deallocate (dl_over_b)
      if (allocated(d_dl_over_b_drho)) deallocate (d_dl_over_b_drho)
      if (allocated(gds2)) deallocate (gds2)
      if (allocated(gds21)) deallocate (gds21)
      if (allocated(gds22)) deallocate (gds22)
      if (allocated(gds23)) deallocate (gds23)
      if (allocated(gds24)) deallocate (gds24)
      if (allocated(gds25)) deallocate (gds25)
      if (allocated(gds26)) deallocate (gds26)
      if (allocated(dgds2dr)) deallocate (dgds2dr)
      if (allocated(dgds21dr)) deallocate (dgds21dr)
      if (allocated(dgds22dr)) deallocate (dgds22dr)
      if (allocated(gbdrift)) deallocate (gbdrift)
      if (allocated(gbdrift0)) deallocate (gbdrift0)
      if (allocated(cvdrift)) deallocate (cvdrift)
      if (allocated(cvdrift0)) deallocate (cvdrift0)
      if (allocated(dgbdriftdrho)) deallocate (dgbdriftdrho)
      if (allocated(dcvdriftdrho)) deallocate (dcvdriftdrho)
      if (allocated(dgbdrift0drho)) deallocate (dgbdrift0drho)
      if (allocated(dcvdrift0drho)) deallocate (dcvdrift0drho)
      if (allocated(dBdrho)) deallocate (dBdrho)
      if (allocated(d2Bdrdth)) deallocate (d2Bdrdth)
      if (allocated(dgradpardrho)) deallocate (dgradpardrho)
      if (allocated(theta_vmec)) deallocate (theta_vmec)

      if (allocated(alpha)) deallocate (alpha)
      if (allocated(zeta)) deallocate (zeta)

      if (allocated(x_displacement_fac)) deallocate (x_displacement_fac)

      !!GA
      if (allocated(gradzeta_grady)) deallocate (gradzeta_grady)
      if (allocated(gradzeta_gradx)) deallocate (gradzeta_gradx)
      if (allocated(gradpar_zeta)) deallocate (gradpar_zeta)
      geoinit = .false.

   end subroutine finish_geometry

end module stella_geometry
