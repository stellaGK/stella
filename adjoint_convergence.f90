module adjoint_convergence

  implicit none

  public :: init_convergence
  public :: omega_convergence1
  public :: omega_convergence2
  public :: deallocate_convergence
  
  private
  
  real :: halfnavg
contains
  
  subroutine init_convergence
    
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky
    use stella_diagnostics, only: navg
    use stella_layouts, only: vmu_lo
    
    use adjoint_field_arrays, only: omega_g, omega
    
    implicit none
    
    real :: halfnavg

    halfnavg = navg/2

    if(.not. allocated(omega)) then 
       allocate(omega(naky,nakx))
       omega = 0.
    end if
    if (.not. allocated(omega_g)) then
       allocate(omega_g(naky,nakx))
       omega_g = 0.
    end if

  end subroutine init_convergence
  
  subroutine omega_convergence1 (istep, converged)

    use kt_grids, only: nakx, naky
    use stella_time, only: code_time
    use stella_diagnostics, only: omega_vs_time, navg

    implicit none
    
    complex, dimension(:,:), allocatable :: diff
    integer, intent (in) :: istep
    integer :: ikx, iky
    real :: max_diff
    logical, intent (out) :: converged 
   
    converged = .False.
    
    if(istep > navg) then
       allocate(diff(naky,nakx))
       do ikx = 1, nakx
          do iky = 1, naky
             if (mod(istep,navg)+1 > 1) then
                diff(iky, ikx) = omega_vs_time(mod(istep,navg)+1,iky,ikx)- omega_vs_time(mod(istep,navg),iky,ikx)
             else
                diff(iky,ikx) = omega_vs_time(1,iky,ikx)- omega_vs_time(navg,iky,ikx)
             end if
          end do
       end do
       max_diff = maxval(abs(aimag(diff)))
       if(max_diff < 5E-004) converged = .True.
      
       deallocate(diff)
    end if

  end subroutine omega_convergence1

  !!Second Convergence Test

  subroutine omega_convergence2 (istep, istep_initial,  converged)

    use stella_diagnostics, only: omega_vs_time, navg, omega_vs_time_short
    use kt_grids, only: nakx, naky
    use adjoint_field_arrays, only: omega

    implicit none
    
    integer, intent (in) :: istep, istep_initial
    
    complex, dimension (:,:), allocatable :: sum_omega, avg_omega
    complex, dimension (:,:), allocatable :: sum_omega_local, avg_omega_local
    complex, dimension(:,:), allocatable :: diff_omega
    real :: max_diff

    logical, intent (out) :: converged

    real, dimension(:,:), allocatable :: imag_omega
    real :: min_omega
    
    if (.not. allocated(sum_omega)) allocate(sum_omega(naky,nakx))
    allocate(sum_omega_local(naky,nakx))
    allocate(avg_omega(naky,nakx))
    allocate(avg_omega_local(naky,nakx))
    allocate(diff_omega(naky,nakx))

    allocate(imag_omega(naky,nakx))
    
    converged = .False.
    if ((istep-istep_initial) > navg) then
       halfnavg = navg/2
       sum_omega = sum(omega_vs_time, dim=1)
       avg_omega = sum_omega/navg
       sum_omega_local = sum(omega_vs_time_short, dim=1)
       avg_omega_local = sum_omega_local/(nint(halfnavg))
       
       imag_omega = abs(aimag(avg_omega))
       min_omega = minval(imag_omega, Mask = imag_omega .gt. 0)
       
       diff_omega = avg_omega - avg_omega_local       
       max_diff = maxval(abs(aimag(diff_omega)))/min_omega

       omega = avg_omega
       if (max_diff < 1E-005) then
          converged = .True.
       end if
    end if
       
    deallocate(avg_omega,avg_omega_local,sum_omega,sum_omega_local,diff_omega)

    deallocate(imag_omega)
  end subroutine omega_convergence2
  
  subroutine deallocate_convergence
    use adjoint_field_arrays, only: omega_g, omega

    implicit none

    deallocate(omega_g)
    deallocate(omega)
    
  end subroutine deallocate_convergence
    
end module adjoint_convergence

