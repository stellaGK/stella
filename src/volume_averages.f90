module volume_averages

  public :: init_volume_averages, finish_volume_averages
  public :: fieldline_average
  public :: volume_average
! public :: flux_surface_average
  
  public :: mode_fac
  
  private

  interface fieldline_average
     module procedure fieldline_average_real
     module procedure fieldline_average_complex
  end interface

  real, dimension (:), allocatable :: mode_fac

contains

  subroutine init_volume_averages

    use kt_grids, only: aky, naky

    implicit none

    if (.not.allocated(mode_fac)) then
       allocate (mode_fac(naky)) ; mode_fac = 2.0
       if (aky(1)<epsilon(0.)) mode_fac(1) = 1.0
    end if

  end subroutine init_volume_averages

  subroutine finish_volume_averages
    implicit none
    if (allocated(mode_fac)) deallocate (mode_fac)
  end subroutine finish_volume_averages


  !==============================================
  !============ FIELD LINE AVERAGE ==============
  !==============================================
  subroutine fieldline_average_real (unavg, avg)

    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky
    use stella_geometry, only: dl_over_b

    implicit none

    real, dimension (:,:,-nzgrid:,:), intent (in) :: unavg
    real, dimension (:,:), intent (out) :: avg

    integer :: it, ia

    ia = 1

    avg = 0.0
    do it = 1, ntubes
       avg = avg + sum(spread(spread(dl_over_b(ia,:),1,naky),2,nakx)*unavg(:,:,:,it),dim=3)
    end do
    avg = avg/real(ntubes)
    
  end subroutine fieldline_average_real

  subroutine fieldline_average_complex (unavg, avg)

    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: nakx, naky
    use stella_geometry, only: dl_over_b

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: unavg
    complex, dimension (:,:), intent (out) :: avg

    integer :: it, ia

    ia = 1

    avg = 0.0
    do it = 1, ntubes
       avg = avg + sum(spread(spread(dl_over_b(ia,:),1,naky),2,nakx)*unavg(:,:,:,it),dim=3)
    end do
    avg = avg/real(ntubes)

  end subroutine fieldline_average_complex

  !==============================================
  !============== VOLUME AVERAGE ================
  !==============================================
  subroutine volume_average (unavg, avg)

    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: naky, nakx
    use stella_geometry, only: dl_over_b

    implicit none

    complex, dimension (:,:,-nzgrid:,:), intent (in) :: unavg
    real, intent (out) :: avg

    integer :: iky, ikx, iz, it, ia

    ia = 1

    avg = 0.
    do it = 1, ntubes
       do iz = -nzgrid, nzgrid
          do ikx = 1, nakx
             do iky = 1, naky
                avg = avg + real(unavg(iky,ikx,iz,it)*conjg(unavg(iky,ikx,iz,it)))*mode_fac(iky)*dl_over_b(ia,iz)
             end do
          end do
       end do
    end do
    avg = avg/real(ntubes)

  end subroutine volume_average

end module volume_averages
