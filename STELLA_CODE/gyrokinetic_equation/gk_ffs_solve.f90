   module gk_ffs_solve

   implicit none 
   
   public :: get_drifts_ffs_itteration
   public :: get_source_ffs_itteration
   private

   contains

   !   subroutine add_correction_ffs (phiin, gin, source_out) 

   !     use grids_z, only: nzgrid, ntubes
   !     use parameters_kxky_grid, only: naky, nakx
   !     use stella_layouts, only: vmu_lo

   !     implicit none 
      
   !     complex, dimension(:, :, -nzgrid:, :), intent(in) :: phiin
   !     complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
   !     complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent (out) :: source_out
      
   !     complex, dimension(:,:,:,:,:), allocatable :: source1

   !     if(.not. allocated(source1)) allocate(source1(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

   ! !      if (drifts_implicit)
   !     source1 = 0.0
   !     source_out = 0.0 
   !     call get_drifts_ffs_itteration (phiin, gin, source1)
   !     call get_source_ffs_itteration (phiin, gin, source_out) 
      
   !     source_out = source_out + source1

   !     if(allocated(source1)) deallocate(source1) 
   !   end subroutine add_correction_ffs

   !   contains 

   subroutine get_source_ffs_itteration (phi, g, source) 

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use calculations_transforms, only: transform_ky2y
      use grids_z, only: nzgrid, ntubes
      use parameters_kxky_grid, only: naky, naky_all, nakx, ikx_max, ny
      use calculations_kxky, only: swap_kxky
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac, maxwell_mu_avg
      use grids_species, only: spec

      use fields, only: advance_fields
      use calculations_gyro_averages, only: gyro_average
      use arrays_gyro_averages, only: j0_ffs
      use calculations_kxky, only: swap_kxky_back
      use calculations_transforms, only: transform_y2ky
      
      use gk_parallel_streaming, only: center_zed, get_dgdz_centered, get_dzed
      use gk_parallel_streaming, only: stream_correction, stream_store_full
      
      use grids_species, only: has_electron_species
      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent (out) :: source

      integer :: ivmu, iv, imu, is, ia, iz, it
      complex, dimension(:, :, :, :), allocatable :: g0
      complex, dimension(:, :, :, :), allocatable :: dgphi_dz, dphi_dz, dgdz
      complex, dimension(:, :, :, :), allocatable :: g0y, g1y, g2y, g3y
      complex, dimension(:, :), allocatable :: g_swap

      real, dimension (:, :), allocatable :: coeff, coeff2
      real :: scratch 

      allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes)) ; g0 = 0.0

      allocate (dgphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes)) ; dgphi_dz = 0.0
      allocate (dphi_dz(naky, nakx, -nzgrid:nzgrid, ntubes)) ; dphi_dz = 0.0
      allocate (dgdz(naky, nakx, -nzgrid:nzgrid, ntubes)) ; dgdz = 0.0

      allocate (g_swap(naky_all, ikx_max)) ; g_swap = 0.0
      allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g0y = 0.0
      allocate (g1y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g1y = 0.0
      allocate (g2y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g2y = 0.0 
      allocate (g3y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g3y = 0.0
      
      allocate (coeff(ny, -nzgrid:nzgrid)) ; coeff = 0.0 
      allocate (coeff2(ny, -nzgrid:nzgrid)) ; coeff2 = 0.0

      source = 0.0 
      
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         
         !> <phi>
         call gyro_average(phi, g0, j0_ffs(:, :, :, ivmu))
         
         !> Full d <phi>/dz
   !!        call get_dgdz_centered(g0, ivmu, dgphi_dz)
         call get_dzed(iv, g0, dgphi_dz)
         !> d phi/dz
   !!        call get_dgdz_centered(phi, ivmu, dphi_dz)
         call get_dzed(iv, phi, dphi_dz)
         !> dg/dz
   !!        call get_dgdz(g(:, :, :, :, ivmu), ivmu, dgdz)
         call get_dzed(iv, g(:, :, :, :, ivmu), dgdz)
         !> get these quantities in real space 
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               !> g0y is real space version of d<phi>/dz
               call swap_kxky(dgphi_dz(:, :, iz, it), g_swap)
               call transform_ky2y(g_swap, g0y(:, :, iz, it))

               !> g1y is real space version of dphi/dz
               call swap_kxky(dphi_dz(:, :, iz, it), g_swap) 
               call transform_ky2y(g_swap, g1y(:, :, iz, it)) 

               !> g2y is real space version of dg/dz
               call swap_kxky(dgdz(:, :, iz, it), g_swap)
               call transform_ky2y(g_swap, g2y(:, :, iz, it)) 
               
            end do
         end do
            
         scratch = maxwell_fac(is) *  maxwell_vpa(iv, is) * spec(is)%zt

         coeff = stream_correction(:,:,iv,is)
         coeff2 = stream_correction(:,:,iv,is) * maxwell_mu_avg(:, :, imu, is) * scratch
         do ia = 1, ny
            call center_zed(iv, coeff(ia, :) ,  -nzgrid)
            call center_zed(iv, coeff2(ia, :) ,  -nzgrid)
         end do
         g2y = spread(spread(coeff, 2, ikx_max), 4, ntubes) * g2y + spread(spread(coeff2, 2, ikx_max), 4, ntubes) * g0y
   !        g2y = spread(spread(stream_correction(:,:,iv,is), 2, ikx_max), 4, ntubes) &
   !             * (g2y + g0y * spread(spread(maxwell_mu_avg(:, :, imu, is), 2, ikx_max),4, ntubes ) * scratch)

         coeff = stream_store_full (:,:,iv,is) * (maxwell_mu(:, :, imu, is) - maxwell_mu_avg(:, :, imu, is)) * scratch
         do ia = 1, ny
            call center_zed(iv, coeff(ia, :) ,  -nzgrid)
         end do
         g3y = spread(spread(coeff, 2, ikx_max), 4, ntubes) * g1y

   !        g3y = spread(spread(stream_store_full (:,:,iv,is) * (maxwell_mu(:, :, imu, is) & 
   !!             - maxwell_mu_avg(:, :, imu, is)), 2, ikx_max),4, ntubes) &
   !             * g1y * scratch 

         coeff = stream_store_full (:,:,iv,is) * maxwell_mu(:, :, imu, is)
         do ia = 1, ny 
            call center_zed(iv, coeff(ia, :) ,  -nzgrid) 
         end do
         g0y =  spread(spread(coeff, 2, ikx_max), 4, ntubes) * (g0y - g1y) * scratch 
   !        g0y = spread(spread(stream_store_full (:,:,iv,is) * maxwell_mu(:, :, imu, is) , 2, ikx_max),4, ntubes) & 
   !             * (g0y - g1y) * scratch 

   !        g0y = 0.0 
         g0y = g0y + g2y + g3y

         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               call transform_y2ky(g0y(:, :, iz, it), g_swap)
               call swap_kxky_back(g_swap, source(:, :, iz, it, ivmu))
            end do
         end do

   !!        call center_zed(iv, source(:,:,:,:,ivmu))
      end do

      deallocate(g0)
      deallocate(dgphi_dz, dphi_dz, dgdz)
      deallocate(g_swap)
      deallocate(g0y,g1y,g2y, g3y) 
      deallocate(coeff, coeff2) 
      end subroutine get_source_ffs_itteration

      subroutine get_drifts_ffs_itteration (phi, g, source) 

         use stella_layouts, only: vmu_lo
         use stella_layouts, only: iv_idx, imu_idx, is_idx
         use calculations_transforms, only: transform_ky2y,transform_y2ky
         use grids_z, only: nzgrid, ntubes
         use parameters_kxky_grid, only: naky, naky_all, nakx, ikx_max, ny
         use calculations_kxky, only: swap_kxky, swap_kxky_back
         use calculations_gyro_averages, only: gyro_average
         use arrays_gyro_averages, only: j0_ffs
         use arrays_store_useful, only: wdriftx_g, wdriftx_phi
         use arrays_store_useful, only: wdrifty_g, wdrifty_phi
         use arrays_store_useful, only: wstar

         use gk_parallel_streaming, only: center_zed

         implicit none 

         complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
         complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
         complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent (out) :: source

         integer :: iz, it, ivmu, iv, is, imu
         complex, dimension(:, :), allocatable :: g_swap
         complex, dimension(:, :, :, :), allocatable :: gk0, gk1, gk2, gk3, gk4
         complex, dimension(:, :, :, :), allocatable :: sourcey, g1y, g2y, g3y, g4y

         it = 1

         allocate (gk0(naky, nakx, -nzgrid:nzgrid, ntubes)) ; gk0 = 0.0
         allocate (gk1(naky, nakx, -nzgrid:nzgrid, ntubes)) ; gk1 = 0.0
         allocate (gk2(naky, nakx, -nzgrid:nzgrid, ntubes)) ; gk2 = 0.0
         allocate (gk3(naky, nakx, -nzgrid:nzgrid, ntubes)) ; gk3 = 0.0
         allocate (gk4(naky, nakx, -nzgrid:nzgrid, ntubes)) ; gk4 = 0.0

         allocate (sourcey(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; sourcey = 0.0
         allocate (g1y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g1y = 0.0
         allocate (g2y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g2y = 0.0
         allocate (g3y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g3y = 0.0
         allocate (g4y(ny, ikx_max, -nzgrid:nzgrid, ntubes)) ; g4y = 0.0

         allocate (g_swap(naky_all, ikx_max))
         
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)

            call gyro_average(phi, gk0, j0_ffs(:, :, :, ivmu))

            call get_dgdy(gk0, gk1)
            call get_dgdx(gk0, gk2)

            call get_dgdy(g(:,:,:,:,ivmu), gk3)
            call get_dgdx(g(:,:,:,:,ivmu), gk4)

            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call swap_kxky(gk1(:, :, iz, it), g_swap)
                  call transform_ky2y(g_swap, g1y(:, :, iz, it))

                  call swap_kxky(gk2(:, :, iz, it), g_swap)
                  call transform_ky2y(g_swap, g2y(:, :, iz, it))

                  call swap_kxky(gk3(:, :, iz, it), g_swap)
                  call transform_ky2y(g_swap, g3y(:, :, iz, it))

                  call swap_kxky(gk4(:, :, iz, it), g_swap)
                  call transform_ky2y(g_swap, g4y(:, :, iz, it))      
               end do
            end do
         
            sourcey = 0.0 
            call add_explicit_term_ffs_fields(g1y, wstar(:,:,ivmu), sourcey)
            !!
            call add_explicit_term_ffs_fields(g4y, wdriftx_g(:,:,ivmu), sourcey) 
            call add_explicit_term_ffs_fields(g2y, wdriftx_phi(:,:,ivmu), sourcey) 
            !!
            call add_explicit_term_ffs_fields(g3y, wdrifty_g(:,:,ivmu), sourcey)
            call add_explicit_term_ffs_fields(g1y, wdrifty_phi(:,:,ivmu), sourcey)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call transform_y2ky(sourcey(:, :, iz, it), g_swap)
                  call swap_kxky_back(g_swap, source(:, :, iz, it, ivmu))
               end do
            end do
            call center_zed(iv, source(:,:,:,:,ivmu))
      end do

      deallocate(sourcey, g_swap)
      deallocate(gk0,gk1,gk2,gk3,gk4)
      deallocate(g1y, g2y, g3y, g4y)

      end subroutine get_drifts_ffs_itteration
   
      subroutine get_dgdy(gin, dgdy)

         use constants, only: zi
         use grids_z, only: nzgrid, ntubes
         use parameters_kxky_grid, only: nakx
         use grids_kxky, only: aky

         implicit none

         complex, dimension(:, :, -nzgrid:, :), intent(in) :: gin
         complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdy

         integer :: it, iz, ikx

         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  dgdy(:, ikx, iz, it) = zi * aky(:) * gin(:, ikx, iz, it)
               end do
            end do
         end do

      end subroutine get_dgdy

      subroutine get_dgdx(gin, dgdx)

         use constants, only: zi
         use grids_z, only: nzgrid, ntubes
         use parameters_kxky_grid, only: nakx
         use grids_kxky, only: akx

         implicit none

         complex, dimension(:, :, -nzgrid:, :), intent(in) :: gin
         complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdx

         integer :: ikx, iz, it

         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  dgdx(:, ikx, iz, it) = zi * akx(ikx) * gin(:, ikx, iz, it)
               end do
            end do
         end do

      end subroutine get_dgdx

      subroutine add_explicit_term_ffs_fields(g, pre_factor, src)

         use grids_z, only: nzgrid, ntubes
         use parameters_kxky_grid, only: ikx_max, nalpha

         implicit none

         complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
         real, dimension(:, -nzgrid:), intent(in) :: pre_factor
         complex, dimension(:, :, -nzgrid:, :), intent(in out) :: src

         integer :: ia, ikx, iz, it

         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ikx = 1, ikx_max
                  do ia = 1, nalpha
                     src(ia, ikx, iz, it) = src(ia, ikx, iz, it) + pre_factor(ia, iz) * g(ia, ikx, iz, it)
                  end do
               end do
            end do
         end do

      end subroutine add_explicit_term_ffs_fields


   end module gk_ffs_solve
