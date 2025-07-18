module haixuan_template

  public :: init_advection_in_zed

  private

contains

  subroutine init_advection_in_zed

    use zgrid, only: nzgrid
    use vpamu_grids, only: vpa, nvpa, mu, nmu
    use geometry, only: b_dot_grad_z, dbdzed, bmag
    
    integer :: iv, iz, imu
    integer :: ia

    ia = 1
    
    do iv = 1, nvpa
       write (*, *) 'iv: ', iv, 'vpa: ', vpa(iv)
    end do
    write (*, *)

    do iz = -nzgrid, nzgrid
       write (*, *) 'b_dot_grad_z: ', b_dot_grad_z(ia,iz), 'dbdzed: ', dbdzed(ia,iz), 'bmag: ', bmag(ia, iz)
    end do
    write (*, *)

    do imu = 1, nmu
       write (*, *) 'mu: ', mu(imu)
    end do
    
  end subroutine init_advection_in_zed

end module haixuan_template
