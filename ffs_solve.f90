module ffs_solve

  implicit none 
  
!  public :: get_drifts_ffs_itteration
  public :: get_source_ffs_itteration
  private

contains

  subroutine get_source_ffs_itteration (phi, phi_bar, g, source) 

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

     complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, phi_bar
     complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
     complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent (out) :: source

     integer :: ivmu, iv, imu, is, ia, iz, it, ikx, iky
     complex, dimension(:, :, :, :), allocatable :: g0
     complex, dimension(:, :, :, :), allocatable :: dgphi_dz, dphi_dz, dgdz
     complex, dimension(:, :), allocatable :: g_swap

     complex, dimension(:, :), allocatable :: g0y, g1y, g2y, g3y 
     real, dimension (:), allocatable :: coeff, coeff2
     real :: scratch
     
     complex, dimension(:,:), allocatable :: source_drifts


     allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes)) ; g0 = 0.0
     
     allocate (dgphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes)) ; dgphi_dz = 0.0
     allocate (dphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes)) ; dphi_dz = 0.0
     allocate (dgdz(naky, nakx, -nzgrid:nzgrid, ntubes)) ; dgdz = 0.0
     
     allocate (g_swap(naky_all, ikx_max)) ; g_swap = 0.0

     allocate (g0y(ny, ikx_max)) ; g0y = 0.0
     allocate (g1y(ny, ikx_max)) ; g1y = 0.0
     allocate (g2y(ny, ikx_max)) ; g2y = 0.0 
     allocate (g3y(ny, ikx_max)) ; g3y = 0.0
     
     allocate (coeff(ny)) ; coeff = 0.0 
     allocate (coeff2(ny)) ; coeff2 = 0.0

     if (drifts_implicit) then
        allocate(source_drifts(ny, ikx_max)); source_drifts = 0.0
     end if
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
        call get_dgdz(phi_bar, ivmu, dphi_dz)
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
                 call get_drifts_ffs_itteration (g0(:, :, iz, it), g(:, :, iz, it, ivmu), source_drifts, ivmu, iz, it)
                 g0y = g2y + g3y + source_drifts
              else
                 !> Add them all together and transform back to (kx, ky) space
                 g0y = g2y + g3y
              end if
              call transform_y2ky(g0y, g_swap)
              call swap_kxky_back(g_swap, source(:, :, iz, it, ivmu))
           end do
        end do
              
        source(1, 1, :, :, :) = 0.0

        !> Note that centering is not done in implicit_solve.f90
        
        ! do iky = 1, naky
        !    do ikx = 1, nakx
        !       call center_zed(iv, source(iky, ikx, :, 1, ivmu), -nzgrid, periodic(iky))
        !    end do
        ! end do
        
     end do
          
     deallocate(g0)
     deallocate(dgphi_dz, dphi_dz, dgdz)
     deallocate(g_swap)
     deallocate(g0y,g1y,g2y, g3y) 
     deallocate(coeff, coeff2)

     if (drifts_implicit) deallocate(source_drifts)

   end subroutine get_source_ffs_itteration

   subroutine get_drifts_ffs_itteration (gyro_phi, gin, sourcey, ivmu, iz, it) 

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

      implicit none 

      complex, dimension(:, :), intent(in) :: gyro_phi!phi_old
      complex, dimension(:, :), intent(in) :: gin
      complex, dimension(:, :), intent (out) :: sourcey 

      integer, intent (in) :: ivmu, iz, it

      integer :: iv, is, imu
      complex, dimension(:, :), allocatable :: g_swap
      complex, dimension(:, :), allocatable :: dgphidy, dgphidx
      complex, dimension(:, :), allocatable :: gk1, gk2, gk3, gk4
      complex, dimension(:, :), allocatable :: g1y, g2y, g3y, g4y
      
      allocate (dgphidy(naky, nakx)) ; dgphidy = 0.0
      allocate (dgphidx(naky, nakx)) ; dgphidx = 0.0
      
      allocate (gk1(naky, nakx)) ; gk1 = 0.0
      allocate (gk2(naky, nakx)) ; gk2 = 0.0
      allocate (gk3(naky, nakx)) ; gk3 = 0.0
      allocate (gk4(naky, nakx)) ; gk4 = 0.0

      allocate (g1y(ny, ikx_max)) ; g1y = 0.0
      allocate (g2y(ny, ikx_max)) ; g2y = 0.0
      allocate (g3y(ny, ikx_max)) ; g3y = 0.0
      allocate (g4y(ny, ikx_max)) ; g4y = 0.0
      
      allocate (g_swap(naky_all, ikx_max)) ; g_swap = 0.0

      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      
      call get_dgdy(gyro_phi, gk1)
      call get_dgdx(gyro_phi, gk2)
      
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

      call add_explicit_term_ffs_fields(g3y, wdrifty_g(:, iz, ivmu), sourcey)
      call add_explicit_term_ffs_fields(g1y, wdrifty_phi(:, iz, ivmu), sourcey)

      call add_explicit_term_ffs_fields(g4y, wdriftx_g(:, iz, ivmu), sourcey)
      call add_explicit_term_ffs_fields(g2y, wdriftx_phi(:, iz, ivmu), sourcey)

      call add_explicit_term_ffs_fields(g1y, wstar(:, iz, ivmu), sourcey)
      
      deallocate(g_swap)
      deallocate(dgphidy, dgphidx)
      deallocate(gk1, gk2, gk3, gk4)
      deallocate(g1y, g2y, g3y, g4y)

    end subroutine get_drifts_ffs_itteration
 
    subroutine get_dgdy(gin, dgdy)

      use constants, only: zi
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, aky
      use kt_grids, only: nakx
      implicit none

      complex, dimension(:, :), intent(in) :: gin
      complex, dimension(:, :), intent(out) :: dgdy

      integer :: it, iz, iky, ikx

      do ikx = 1, nakx
         dgdy(:, ikx) = zi * aky(:) * gin(:, ikx)
      end do
      
   end subroutine get_dgdy

   subroutine get_dgdx(gin, dgdx)

      use constants, only: zi
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: akx, nakx

      implicit none

      complex, dimension(:, :), intent(in) :: gin
      complex, dimension(:, :), intent(out) :: dgdx

      integer :: ikx, iz, it

      do ikx = 1, nakx
         dgdx(:, ikx) = zi * akx(ikx) * gin(:, ikx)
      end do
      
   end subroutine get_dgdx

   subroutine add_explicit_term_ffs_fields(g, pre_factor, src)

      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: ikx_max, ny

      implicit none

      complex, dimension(:, :), intent(in) :: g
      real, dimension(:), intent(in) :: pre_factor
      complex, dimension(:, :), intent(in out) :: src

      integer :: ivmu
      integer :: ia, ikx, iz, it

      do ia = 1, ny
         src(ia, :) = src(ia, :) + ( pre_factor(ia) - pre_factor(1) ) * g(ia, :)
!         src(ia, :) = src(ia, :) + pre_factor(ia) * g(ia, :)
      end do

   end subroutine add_explicit_term_ffs_fields


end module ffs_solve
