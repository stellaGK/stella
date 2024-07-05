module zpinch
  
  implicit none
  
  public :: get_zpinch_geometry_coefficients

  private

contains

  ! subroutine read_zpinch_parameters (radius)

  !   use file_utils, only: input_unit_exist
    
  !   implicit none

  !   real, intent (out) :: radius
    
  !   integer :: in_file
  !   logical :: exist
    
  !   namelist /zpinch_parameters/ radius

  !   ! default value for radius
  !   radius = 0.5

  !   in_file = input_unit_exist("zpinch_parameters", exist)
  !   if (exist) read (unit=in_file, nml=zpinch_parameters)
    
  ! end subroutine read_zpinch_parameters

  ! use Z-pinch equilibrium.
  ! the parallel coordinate, z, is chosen to be the arc-length,
  ! i.e., the physical poloidal angle times the radius of the chosen magnetic field, r0.
  ! the radial coordinate, x, is directed away from the middle of the circle
  ! and is normalised by the local radius r0.
  ! the bi-normal coordinate, y, is chosen to form an orthogonal, right-handed coordinate
  ! system with (y,x,z) and is also normalised by r0.
  subroutine get_zpinch_geometry_coefficients (nzgrid, bmag, gradpar, grad_rho, &
       grad_y_dot_grad_y, grad_x_dot_grad_y, grad_x_dot_grad_x, &
       gbdrift0, gbdrift, cvdrift0, cvdrift, btor, rmajor)

    implicit none

    integer, intent(in) :: nzgrid
    real, dimension (-nzgrid:), intent(out) :: bmag, gradpar, grad_rho, &
         grad_y_dot_grad_y, grad_x_dot_grad_y, grad_x_dot_grad_x, &
         gbdrift0, gbdrift, cvdrift0, cvdrift, btor, rmajor

    ! bmag = B(r0) / B_ref
    ! as B is constant along radius r0, choose B_ref = B(r0), so bmag = 1
    bmag = 1.0
    ! gradpar = bhat . grad z = b . (r0*grad) theta = 1
    gradpar = 1.0
    ! grad_rho = | r0 * grad (r / r0) | = 1
    grad_rho = 1.0
    ! grad_y_dot_grad_y = 1
    grad_y_dot_grad_y = 1.0
    ! grad_x_dot_grad_y = 0 because (x,y) are orthogonal coordinates
    grad_x_dot_grad_y = 0.0
    ! grad_x_dot_grad_x = 1
    grad_x_dot_grad_x = 1.0
    ! the x-component of the grad-B drift is proportional to gbdrift0; zero in a z-pinch
    gbdrift0 = 0.0
    ! the x-component of the curvature drift is proportional to cvdrift0; zero in a z-pinch
    cvdrift0 = 0.0
    ! gbdrift = 2 * bhat / B_norm x (grad_norm B_norm / B_norm) . grad y
    ! = 2 * bhat / B_norm x (grad_norm x * d ln B_norm / dx) . grad y
    ! = -2 * d(ln B_norm) / dx_norm = 2 * r0 / L_B = 1 (from MHD equilibrium with beta=0)
    gbdrift = 2.0
    ! cvdrift = 2 * bhat / B_norm x (bhat . grad_norm bhat) . grad y = 2
    cvdrift = gbdrift
    
    ! btor and rmajor used for flow shear calculations;
    ! choosing btor = 0 means all flow shear is perpendicular rather than parallel
    btor = 0.0
    ! rmajor set to one to avoid divide by zero and otherwise is unused
    rmajor = 1.0
    
  end subroutine get_zpinch_geometry_coefficients
  
end module zpinch
