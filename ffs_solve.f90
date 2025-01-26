module ffs_solve

  implicit none 
  
!  public :: get_drifts_ffs_itteration
  public :: get_source_ffs_itteration
  private

contains

  subroutine get_source_ffs_itteration (phi, g, source) 

     use mp, only: proc0
     use stella_layouts, only: vmu_lo
     use stella_layouts, only: iv_idx, imu_idx, is_idx
     use stella_transforms, only: transform_ky2y
     use zgrid, only: nzgrid, ntubes
     use kt_grids, only: naky, naky_all, nakx, ikx_max, ny
     use kt_grids, only: swap_kxky
     use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac, maxwell_mu_avg
     use species, only: spec

     use fields, only: advance_fields, fields_updated
     use gyro_averages, only: j0_ffs, j0_const, gyro_average

     use kt_grids, only: swap_kxky_back
     use stella_transforms, only: transform_y2ky
     
     use parallel_streaming, only: center_zed, get_dgdz
     use parallel_streaming, only: stream_correction, stream, stream_store_full
     
     use species, only: has_electron_species
     use extended_zgrid, only: periodic

     use run_parameters, only: drifts_implicit
     implicit none

     complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
     complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
     complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent (out) :: source

     integer :: ivmu, iv, imu, is, ia, iz, it, ikx
     complex, dimension(:, :, :, :), allocatable :: g0
     complex, dimension(:, :, :, :), allocatable :: dgphi_dz, dphi_dz, dgdz
!     complex, dimension(:, :, :, :), allocatable :: g0y, g1y, g2y, g3y 
     complex, dimension(:, :), allocatable :: g_swap

     complex, dimension(:, :), allocatable :: g0y, g1y, g2y, g3y 
     real, dimension (:), allocatable :: coeff, coeff2

     complex, dimension(:,:), allocatable :: source_drifts
!     real, dimension (:, :), allocatable :: coeff, coeff2
     
     real :: scratch

     integer :: iky
     allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes)) ; g0 = 0.0
     
     allocate (dgphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes)) ; dgphi_dz = 0.0
     allocate (dphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes)) ; dphi_dz = 0.0
     allocate (dgdz(naky, nakx, -nzgrid:nzgrid, ntubes)) ; dgdz = 0.0

     ! allocate (dgphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes)) ; dgphi_dz = 0.0
     ! allocate (dphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes)) ; dphi_dz = 0.0
     ! allocate (dgdz(naky, nakx, -nzgrid:nzgrid, ntubes)) ; dgdz = 0.0
     
     allocate (g_swap(naky_all, ikx_max)) ; g_swap = 0.0
     ! allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g0y = 0.0
     ! allocate (g1y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g1y = 0.0
     ! allocate (g2y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g2y = 0.0 
     ! allocate (g3y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g3y = 0.0

     allocate (g0y(ny, ikx_max)) ; g0y = 0.0
     allocate (g1y(ny, ikx_max)) ; g1y = 0.0
     allocate (g2y(ny, ikx_max)) ; g2y = 0.0 
     allocate (g3y(ny, ikx_max)) ; g3y = 0.0
     
     allocate (coeff(ny)) ; coeff = 0.0 
     allocate (coeff2(ny)) ; coeff2 = 0.0

     if (drifts_implicit) allocate(source_drifts(ny, ikx_max))
     source = 0.0 
     
     do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        iv = iv_idx(vmu_lo, ivmu)
        imu = imu_idx(vmu_lo, ivmu)
        is = is_idx(vmu_lo, ivmu)
        
        !> <phi>
        call gyro_average(phi, g0, j0_ffs(:, :, :, ivmu))
        !> Full d <phi>/dz
        call get_dgdz(g0, ivmu, dgphi_dz)
        !> d phi/dz
        call get_dgdz(phi, ivmu, dphi_dz)
        !> dg/dz
        call get_dgdz(g(:, :, :, :, ivmu), ivmu, dgdz)
        
        !> get these quantities in real space 
        do it = 1, ntubes
           do iz = -nzgrid, nzgrid

              !> g0y is real space version of d<phi>/dz
              call swap_kxky(dgphi_dz(:, :, iz, it), g_swap)
              call transform_ky2y(g_swap, g0y)

              !> g1y is real space version of dphi/dz
              call swap_kxky(dphi_dz(:, :, iz, it), g_swap)
              call transform_ky2y(g_swap, g1y)
              !> g2y is real space version of dg/dz
              call swap_kxky(dgdz(:, :, iz, it), g_swap)
              call transform_ky2y(g_swap, g2y)

              scratch = maxwell_fac(is) *  maxwell_vpa(iv, is) * spec(is)%zt

              !> This is (b.grad z - avg(b.grad z)) * (dg^j/dz)
              coeff = stream_correction(:,iz,iv,is)
              g2y = spread(coeff, 2, ikx_max) * g2y

              coeff = (stream(:,iz,iv,is) + stream_correction(:,iz,iv,is))  * maxwell_mu(:, iz, imu, is) * scratch
              coeff2 = stream(:, iz, iv, is) * maxwell_mu_avg(:, iz, imu, is) * scratch

              !> This is (b.grad z * F_0 * Z/T) * d<phi>/dz - avg(b.grad z) * avg(F_0) * d(phibar)/dz
              g3y = spread(coeff, 2, ikx_max) * g0y - spread(coeff2, 2, ikx_max) * g1y

              if(drifts_implicit) then
                 call get_drifts_ffs_itteration (phi(:, :, iz, it), g(:, :, iz, it, ivmu), source_drifts, ivmu, iz, it)
                 g0y = g2y + g3y + source_drifts
              else
                 !> Add them all together and transform back to (kx, ky) space
                 g0y = g2y + g3y
              end if
              call transform_y2ky(g0y, g_swap)
              call swap_kxky_back(g_swap, source(:, :, iz, it, ivmu))
           end do
        end do
        !       !> g0y is real space version of d<phi>/dz
        !       call swap_kxky(dgphi_dz(:, :, iz, it), g_swap)
        !       call transform_ky2y(g_swap, g0y(:, :, iz, it))

        !       !> g1y is real space version of dphi/dz
        !       call swap_kxky(dphi_dz(:, :, iz, it), g_swap) 
        !       call transform_ky2y(g_swap, g1y(:, :, iz, it)) 

        !       !> g2y is real space version of dg/dz
        !       call swap_kxky(dgdz(:, :, iz, it), g_swap)
        !       call transform_ky2y(g_swap, g2y(:, :, iz, it)) 

        !    end do
        ! end do
         
        ! scratch = maxwell_fac(is) *  maxwell_vpa(iv, is) * spec(is)%zt

        ! !> This is (b.grad z - avg(b.grad z)) * (dg^j/dz + Z/T * avg(F_0) * d<phi^j> /dz)
        ! coeff = stream_correction(:,:,iv,is)
        ! coeff2 = stream_correction(:,:,iv,is) * maxwell_mu_avg(:, :, imu, is) * scratch
        ! g2y = spread(spread(coeff, 2, ikx_max), 4, ntubes) * g2y + spread(spread(coeff2, 2, ikx_max), 4, ntubes) * g1y

        ! !> This is (b.grad z) * Z/T * (F_0 - avg(F_0)) * d phi^j /dz
        ! coeff = stream_store_full (:,:,iv,is) * (maxwell_mu(:, :, imu, is) - maxwell_mu_avg(:, :, imu, is)) * scratch
        ! g3y = spread(spread(coeff, 2, ikx_max), 4, ntubes) * g1y

        ! !> This is (b.grad z) * Z/T * F_0 * d(<phi^j> - phi^j)/dz
        ! coeff = stream_store_full (:,:,iv,is) * maxwell_mu(:, :, imu, is) * scratch
        ! g0y =  spread(spread(coeff, 2, ikx_max), 4, ntubes) * (g0y - g1y)

        ! !> Add them all together and transform back to (kx, ky) space
        ! g0y = g0y + g2y + g3y

              
        ! do it = 1, ntubes
        !    do iz = -nzgrid, nzgrid
        !       call transform_y2ky(g0y(:, :, iz, it), g_swap)
        !       call swap_kxky_back(g_swap, source(:, :, iz, it, ivmu))
        !    end do
        ! end do

        source(1, 1, :, :, :) = 0.0

        do iky = 1, naky
           do ikx = 1, nakx
              call center_zed(iv, source(iky, ikx, :, 1, ivmu), -nzgrid, periodic(iky))
           end do
        end do
        
     end do
          
     deallocate(g0)
     deallocate(dgphi_dz, dphi_dz, dgdz)
     deallocate(g_swap)
     deallocate(g0y,g1y,g2y, g3y) 
     deallocate(coeff, coeff2)

     if (drifts_implicit) deallocate(source_drifts)

   end subroutine get_source_ffs_itteration

   subroutine get_drifts_ffs_itteration (phi_old, gin, sourcey, ivmu, iz, it) 

      use mp, only: proc0
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use stella_transforms, only: transform_ky2y,transform_y2ky
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, naky_all, nakx, ikx_max, ny
      use kt_grids, only: swap_kxky, swap_kxky_back
      use gyro_averages, only: j0_ffs, gyro_average

      use dist_fn_arrays, only: wdriftx_g, wdriftx_phi
      use dist_fn_arrays, only: wdrifty_g, wdrifty_phi
      use dist_fn_arrays, only: wstar

      use parallel_streaming, only: center_zed

      use extended_zgrid, only: periodic
      use fields, only: get_dchidy, get_dchidx
      use fields_arrays, only: apar, bpar 
      implicit none 

      complex, dimension(:, :), intent(in) :: phi_old
      complex, dimension(:, :), intent(in) :: gin
      complex, dimension(:, :), intent (out) :: sourcey !source_drifts

      integer, intent (in) :: ivmu, iz, it
!      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
!      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
!      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent (out) :: source

      integer :: iv, is, imu!, iky, ikx
      complex, dimension(:, :), allocatable :: g_swap
      complex, dimension(:, :), allocatable :: dgphidy, dgphidx
      complex, dimension(:, :), allocatable :: gk1, gk2, gk3, gk4
!      complex, dimension(:, :), allocatable :: sourcey
      complex, dimension(:, :), allocatable :: g1y, g2y, g3y, g4y
      
!      complex, dimension(:, :, :, :), allocatable :: gk0, gk1, gk2, gk3, gk4
!      complex, dimension(:, :, :, :), allocatable :: sourcey, g1y, g2y, g3y, g4y

      !it = 1

      !      allocate (gk0(naky, nakx, -nzgrid:nzgrid, ntubes)) ; gk0 = 0.0

      allocate (dgphidy(naky, nakx)) ; dgphidy = 0.0
      allocate (dgphidx(naky, nakx)) ; dgphidx = 0.0
!      allocate (gk0((naky, nakx)) ; gk0 = 0.0
      allocate (gk1(naky, nakx)) ; gk1 = 0.0
      allocate (gk2(naky, nakx)) ; gk2 = 0.0
      allocate (gk3(naky, nakx)) ; gk3 = 0.0
      allocate (gk4(naky, nakx)) ; gk4 = 0.0

!      allocate (sourcey(ny, ikx_max)) ; sourcey = 0.0
      allocate (g1y(ny, ikx_max)) ; g1y = 0.0
      allocate (g2y(ny, ikx_max)) ; g2y = 0.0
      allocate (g3y(ny, ikx_max)) ; g3y = 0.0
      allocate (g4y(ny, ikx_max)) ; g4y = 0.0
      
      ! allocate (gk1(naky, nakx, -nzgrid:nzgrid, ntubes)) ; gk1 = 0.0
      ! allocate (gk2(naky, nakx, -nzgrid:nzgrid, ntubes)) ; gk2 = 0.0
      ! allocate (gk3(naky, nakx, -nzgrid:nzgrid, ntubes)) ; gk3 = 0.0
      ! allocate (gk4(naky, nakx, -nzgrid:nzgrid, ntubes)) ; gk4 = 0.0

      ! allocate (sourcey(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; sourcey = 0.0
      ! allocate (g1y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g1y = 0.0
      ! allocate (g2y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g2y = 0.0
      ! allocate (g3y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g3y = 0.0
      ! allocate (g4y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g4y = 0.0

      allocate (g_swap(naky_all, ikx_max))

      apar = 0.0
      bpar = 0.0
      
!       do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
!          iv = iv_idx(vmu_lo, ivmu)
!          imu = imu_idx(vmu_lo, ivmu)
!          is = is_idx(vmu_lo, ivmu)

!          do it = 1, ntubes
!             do iz = -nzgrid, nzgrid
               
!                call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), dgphidy)
!                call get_dchidx(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), dgphidx)

!                ! call gyro_average(phi, gk0, j0_ffs(:, :, :, ivmu))
!                ! call get_dgdy(gk0, gk1)
!                ! call get_dgdx(gk0, gk2)
               
!                call gyro_average(dgphidy, gk1, j0_ffs(:, :, iz, ivmu))
!                call gyro_average(dgphidx, gk2, j0_ffs(:, :, iz, ivmu))

! !               call get_dgdy(g(:, :, iz, it, ivmu), gk3)
! !               call get_dgdx(g(:, :, iz, it, ivmu), gk4)

!                call swap_kxky(gk1, g_swap)
!                call transform_ky2y(g_swap, g1y)
! !               call swap_kxky(gk1(:, :, iz, it), g_swap)
! !               call transform_ky2y(g_swap, g1y(:, :, iz, it))

!                call swap_kxky(gk2, g_swap)
!                call transform_ky2y(g_swap, g2y)
! !               call swap_kxky(gk2(:, :, iz, it), g_swap)
! !               call transform_ky2y(g_swap, g2y(:, :, iz, it))
               
!                call swap_kxky(gk3, g_swap)
!                call transform_ky2y(g_swap, g3y)
! !               call swap_kxky(gk3(:, :, iz, it), g_swap)
! !               call transform_ky2y(g_swap, g3y(:, :, iz, it))

!                call swap_kxky(gk4, g_swap)
!                call transform_ky2y(g_swap, g4y)
! !               call swap_kxky(gk4(:, :, iz, it), g_swap)
! !               call transform_ky2y(g_swap, g4y(:, :, iz, it))      
! !            end do
! !         end do
       
!                sourcey = 0.0
!                call add_explicit_term_ffs_fields(g1y, wstar(:, iz, ivmu), sourcey)
! !               call add_explicit_term_ffs_fields(g1y, wstar(:,:,ivmu), sourcey)
!                !!
!                call add_explicit_term_ffs_fields(g4y, wdriftx_g(:, iz, ivmu), sourcey)
!                call add_explicit_term_ffs_fields(g2y, wdriftx_phi(:, iz, ivmu), sourcey)
! !         call add_explicit_term_ffs_fields(g4y, wdriftx_g(:,:,ivmu), sourcey) 
! !         call add_explicit_term_ffs_fields(g2y, wdriftx_phi(:,:,ivmu), sourcey) 
!                !!
!                call add_explicit_term_ffs_fields(g3y, wdrifty_g(:, iz, ivmu), sourcey)
!                call add_explicit_term_ffs_fields(g1y, wdrifty_phi(:, iz, ivmu), sourcey)
! !         call add_explicit_term_ffs_fields(g3y, wdrifty_g(:,:,ivmu), sourcey)
!                !         call add_explicit_term_ffs_fields(g1y, wdrifty_phi(:,:,ivmu), sourcey)

!                call transform_y2ky(sourcey, g_swap)
!                call swap_kxky_back(g_swap, source(:, :, iz, it, ivmu))
!             end do
!          end do
! !         do it = 1, ntubes
! !            do iz = -nzgrid, nzgrid
! !               call transform_y2ky(sourcey(:, :, iz, it), g_swap)
! !               call swap_kxky_back(g_swap, source(:, :, iz, it, ivmu))
! !            end do
! !         end do
!          do iky = 1, naky
!             do ikx = 1, nakx
!                call center_zed(iv, source(iky, ikx,:,1,ivmu), -nzgrid, periodic(iky) )
!             end do
!          end do
!       end do

      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      
      call get_dchidy(iz, ivmu, phi_old(:, :), apar(:, :, iz, it), bpar(:, :, iz, it), dgphidy)
      call get_dchidx(iz, ivmu, phi_old(:, :), apar(:, :, iz, it), bpar(:, :, iz, it), dgphidx)

      call gyro_average(dgphidy, gk1, j0_ffs(:, :, iz, ivmu))
      call gyro_average(dgphidx, gk2, j0_ffs(:, :, iz, ivmu))

      call get_dgdy(gin, gk3)
      call get_dgdx(gin, gk4)

      call swap_kxky(gk1, g_swap)
      call transform_ky2y(g_swap, g1y)

      call swap_kxky(gk2, g_swap)
      call transform_ky2y(g_swap, g2y)

      call swap_kxky(gk3, g_swap)
      call transform_ky2y(g_swap, g3y)
      
      call swap_kxky(gk4, g_swap)
      call transform_ky2y(g_swap, g4y)
      
      sourcey = 0.0
      call add_explicit_term_ffs_fields(g1y, wstar(:, iz, ivmu), sourcey)

      call add_explicit_term_ffs_fields(g4y, wdriftx_g(:, iz, ivmu), sourcey)
      call add_explicit_term_ffs_fields(g2y, wdriftx_phi(:, iz, ivmu), sourcey)

      call add_explicit_term_ffs_fields(g3y, wdrifty_g(:, iz, ivmu), sourcey)
      call add_explicit_term_ffs_fields(g1y, wdrifty_phi(:, iz, ivmu), sourcey)

!      call transform_y2ky(sourcey, g_swap)
!      call swap_kxky_back(g_swap, source_drifts)
            
      
!      deallocate(sourcey,
      deallocate(g_swap)
      deallocate(dgphidy, dgphidx, gk1, gk2, gk3, gk4)
      deallocate(g1y, g2y, g3y, g4y)

    end subroutine get_drifts_ffs_itteration
 
    subroutine get_dgdy(gin, dgdy)

      use constants, only: zi
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, aky

      implicit none

!      complex, dimension(:, :, -nzgrid:, :), intent(in) :: gin
!      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdy
      complex, dimension(:, :), intent(in) :: gin
      complex, dimension(:, :), intent(out) :: dgdy

      integer :: it, iz, iky

      ! do it = 1, ntubes
      !    do iz = -nzgrid, nzgrid
      !       do ikx = 1, nakx
      !          dgdy(:, ikx, iz, it) = zi * aky(:) * gin(:, ikx, iz, it)
      !       end do
      !    end do
      ! end do

      do iky = 1, naky
         dgdy(iky, :) = zi * aky(iky) * gin(iky, :)
      end do
      
   end subroutine get_dgdy

   subroutine get_dgdx(gin, dgdx)

      use constants, only: zi
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: akx, nakx

      implicit none

!      complex, dimension(:, :, -nzgrid:, :), intent(in) :: gin
!      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdx
      complex, dimension(:, :), intent(in) :: gin
      complex, dimension(:, :), intent(out) :: dgdx

      integer :: ikx, iz, it

      ! do it = 1, ntubes
      !    do iz = -nzgrid, nzgrid
      !       do ikx = 1, nakx
      !          dgdx(:, ikx, iz, it) = zi * akx(ikx) * gin(:, ikx, iz, it)
      !       end do
      !    end do
      ! end do

      do ikx = 1, nakx
         dgdx(:, ikx) = zi * akx(ikx) * gin(:, ikx)
      end do
      
   end subroutine get_dgdx

   ! subroutine add_explicit_term_ffs_fields(g, pre_factor, src)

   !    use stella_layouts, only: vmu_lo
   !    use zgrid, only: nzgrid, ntubes
   !    use kt_grids, only: ikx_max, nalpha

   !    implicit none

   !    complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
   !    real, dimension(:, -nzgrid:), intent(in) :: pre_factor
   !    complex, dimension(:, :, -nzgrid:, :), intent(in out) :: src

   !    integer :: ivmu
   !    integer :: ia, ikx, iz, it

   !    do it = 1, ntubes
   !       do iz = -nzgrid, nzgrid
   !          do ikx = 1, ikx_max
   !             do ia = 1, nalpha
   !                src(ia, ikx, iz, it) = src(ia, ikx, iz, it) + pre_factor(ia, iz) * g(ia, ikx, iz, it)
   !             end do
   !          end do
   !       end do
   !    end do

   ! end subroutine add_explicit_term_ffs_fields

   subroutine add_explicit_term_ffs_fields(g, pre_factor, src)

      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: ikx_max, nalpha

      implicit none

      complex, dimension(:, :), intent(in) :: g
      real, dimension(:), intent(in) :: pre_factor
      complex, dimension(:, :), intent(in out) :: src

      integer :: ivmu
      integer :: ia, ikx, iz, it

      do ikx = 1, ikx_max
         do ia = 1, nalpha
            src(ia, ikx) = src(ia, ikx) + pre_factor(ia) * g(ia, ikx)
         end do
      end do

   end subroutine add_explicit_term_ffs_fields


end module ffs_solve
