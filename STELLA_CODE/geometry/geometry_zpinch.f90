!###############################################################################
!                         Define a z-pinch equilibrium                          
!###############################################################################
! 
! The parallel coordinate, z, is chosen to be the arc-length z=r0*theta,
! i.e., the physical poloidal angle times the radius of the chosen magnetic field, r0.
! the radial coordinate, x, is directed away from the middle of the circle
! and is normalised by the local magnetic field gradient scale length, L_B = (d ln B / dr)^{-1}
! the bi-normal coordinate, y, is chosen to form an orthogonal, right-handed coordinate
! system with (y,x,z) and is also normalised by L_B.
! 
! In terms of cyclindrical coordinates (r, theta, Z), we have x = r, y = Z, z = r*theta
! 
!###############################################################################
module zpinch

   implicit none

   public :: get_zpinch_geometry_coefficients

   private

   real :: betaprim
   
contains

   !****************************************************************************
   !                Calculate geometric quantities for a z-pinch                
   !****************************************************************************
   subroutine get_zpinch_geometry_coefficients(nzgrid, bmag, b_dot_gradz, grad_rho, surf, &
      grad_y_dot_grad_y, grad_x_dot_grad_y, grad_x_dot_grad_x, B_times_gradB_dot_gradx, &
      B_times_gradB_dot_grady, B_times_kappa_dot_gradx, B_times_kappa_dot_grady, btor, rmajor)

      use common_types, only: flux_surface_type
      use namelist_geometry, only: read_namelist_geometry_zpinch
     
      implicit none

      ! Arguments
      integer, intent(in) :: nzgrid
      real, dimension(-nzgrid:), intent(out) :: bmag, b_dot_gradz, grad_rho
      real, dimension(-nzgrid:), intent(out) :: grad_y_dot_grad_y, grad_x_dot_grad_y
      real, dimension(-nzgrid:), intent(out) :: grad_x_dot_grad_x, B_times_gradB_dot_gradx
      real, dimension(-nzgrid:), intent(out) :: B_times_gradB_dot_grady, B_times_kappa_dot_gradx
      real, dimension(-nzgrid:), intent(out) :: B_times_kappa_dot_grady, btor, rmajor
      type(flux_surface_type), intent(out) :: surf

      !-------------------------------------------------------------------------
      
      ! Read the zpinch namelist in the input file
      call read_namelist_geometry_zpinch (betaprim)

      ! As B is constant along radius r0, choose B_ref = B(r0), so bmag = B(r0) / B_ref = 1
      bmag = 1.0
      
      ! Define b_dot_gradz = bhat . grad z = b . (r0*grad) theta = 1
      b_dot_gradz = 1.0
      
      ! Define grad_rho = | L_B * grad (r / L_B) | = 1
      grad_rho = 1.0
      
      ! Define grad_y_dot_grad_y = 1
      grad_y_dot_grad_y = 1.0
      
      ! Define grad_x_dot_grad_y = 0 because (x,y) are orthogonal coordinates
      grad_x_dot_grad_y = 0.0
      
      ! Define grad_x_dot_grad_x = 1
      grad_x_dot_grad_x = 1.0
      
      ! Define the x-component of the grad-B drift is proportional to B_times_gradB_dot_gradx; zero in a z-pinch
      B_times_gradB_dot_gradx = 0.0
      
      ! Define the x-component of the curvature drift is proportional to B_times_kappa_dot_gradx; zero in a z-pinch
      B_times_kappa_dot_gradx = 0.0
      
      ! Define gbdrift = 2 * bhat / B_norm x (grad_norm B_norm / B_norm) . grad y
      ! = 2 * bhat / B_norm x (grad_norm x * d ln B_norm / dx) . grad y
      ! = -2 * d(ln B_norm) / dx_norm = 2 * r0 / L_B = 1 (from MHD equilibrium with beta=0)
      ! We redefined B_times_gradB_dot_grady = gbdrift / 2
      B_times_gradB_dot_grady = 1.0
      
      ! Define B_times_kappa_dot_grady = bhat / B_norm x (bhat . grad_norm bhat) . grad y = 2
      ! We redefined B_times_kappa_dot_grady = cvdrift / 2
      B_times_kappa_dot_grady = B_times_gradB_dot_grady + betaprim

      ! Note that btor and rmajor used for flow shear calculations;
      ! choosing btor = 0 means all flow shear is perpendicular rather than parallel
      btor = 0.0
      
      ! Here rmajor is set to one to avoid divide by zero and otherwise is unused
      rmajor = 1.0

      ! Note that <surf> contains artificial information about the flux surface shape
      ! that is really not needed for a z-pinch but must be given due to the use of 
      ! gds22 and gds21 as variables rather than grad_x_dot_grad_x and grad_x_dot_grad_y
      ! Setting the magnetic shear to 1, along with q and rhoc means that |grad_x| = |grad q|
      surf%shat = 1.0
      surf%rhoc = 1.0
      surf%qinp = 1.0
      surf%kappa = 1.0
      
   end subroutine get_zpinch_geometry_coefficients

end module zpinch
