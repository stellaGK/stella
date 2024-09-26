module write_radial_grid

  implicit none
  
  public :: dump_radial_grid
  private
  
contains
  
    subroutine dump_radial_grid (x, rho, nx)

        use file_utils, only: run_name
        use parameters_physics, only: rhostar
        use geometry, only: q_as_x, geo_surf, dxdpsi, drhodpsip

        implicit none

        real, dimension(:), intent (in) :: x, rho
        integer, intent (in) :: nx
        
        integer :: ix
        character(300) :: filename

        filename = trim(trim(run_name)//'.radial_grid')
        open (1047, file=filename, status='unknown')
        if (q_as_x) then
        write (1047, '(1a12,1e12.4,1a12,1e12.4,1a12,1e12.4,1a12,1e12.4)') &
            '#dxdpsi = ', dxdpsi, &
            ' q    = ', geo_surf%qinp, &
            ' dqdr = ', geo_surf%shat * geo_surf%qinp / geo_surf%rhoc, &
            ' d2qdr2 = ', geo_surf%d2qdr2
        write (1047, '(3a12)') '#1.x', '2.q', '3.rho'
        do ix = 1, nx
            write (1047, '(3e12.4,i9)') &
                x(ix), &
                x(ix) / dxdpsi + geo_surf%qinp, &
                rho(ix) + geo_surf%rhoc
        end do
        else
        write (1047, '(1a12,1e12.4,1a12,1e12.4,1a12,1e12.4,1a12,1e12.4)') &
            '#dxdpsi = ', dxdpsi, &
            ' dpsidr    = ', 1.0 / drhodpsip, &
            ' d2psidr2 = ', geo_surf%d2psidr2
        write (1047, '(3a12,a9)') '#1.x', '2.psi-psi0', '3.rho'
        do ix = 1, nx
            write (1047, '(3e12.4,i9)') &
                x(ix), &
                rhostar * x(ix) / dxdpsi, &
                rho(ix) + geo_surf%rhoc
        end do
        end if

        close (1047)

    end subroutine dump_radial_grid

 end module write_radial_grid
