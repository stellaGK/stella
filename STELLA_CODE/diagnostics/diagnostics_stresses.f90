


module diagnostics_stresses

  implicit none

contains

  subroutine write_stresses_to_netcdf_file(nout)

      ! Knowledge of first processor
      use mp, only: proc0
    
      ! Dimensions
      use grids_kxky, only:  nakx
      use grids_z, only: nzgrid
            
      ! Flags
      use parameters_physics, only: include_nonlinear
      use parameters_diagnostics, only: write_stresses

      !load data
      use arrays_distribution_function, only: gnew

      !Write to netCDF file
      use write_diagnostics_to_netcdf, only: write_stresses_nc

      integer, intent(in) :: nout
      !We want to write stresses(z, kx) to the netcdf file
      complex, dimension(:,:), allocatable :: phi_h_nonlin_kx, apar_h_nonlin_kx, bpar_h_nonlin_kx
      complex, dimension(:,:), allocatable :: gint

      allocate (phi_h_nonlin_kx(nakx,-nzgrid:nzgrid))
      allocate (apar_h_nonlin_kx(nakx,-nzgrid:nzgrid))
      allocate (bpar_h_nonlin_kx(nakx,-nzgrid:nzgrid))
      allocate (gint(nakx,-nzgrid:nzgrid))

      if (write_stresses .and. include_nonlinear) then
         call calculate_turbulent_stresses(gnew, phi_h_nonlin_kx,apar_h_nonlin_kx, bpar_h_nonlin_kx, gint)
      else
         phi_h_nonlin_kx = 0.
         apar_h_nonlin_kx = 0.
         bpar_h_nonlin_kx = 0.
         gint = 0.
      end if

      if (write_stresses) then
         if (proc0) call write_stresses_nc(nout, phi_h_nonlin_kx, apar_h_nonlin_kx, bpar_h_nonlin_kx, gint)
      end if

      if (allocated(phi_h_nonlin_kx)) deallocate (phi_h_nonlin_kx)
      if (allocated(apar_h_nonlin_kx)) deallocate (apar_h_nonlin_kx)
      if (allocated(bpar_h_nonlin_kx)) deallocate (bpar_h_nonlin_kx)
      if (allocated(gint)) deallocate (gint)
      
  end subroutine write_stresses_to_netcdf_file


  subroutine calculate_turbulent_stresses(g, phi_h_nonlin_kx,apar_h_nonlin_kx, bpar_h_nonlin_kx, gint)

      ! Constants
      use constants, only: pi, zi

    
      ! Parallelisation
      use mp, only: proc0, min_allreduce, sum_reduce
      use mp, only: scope, allprocs, subprocs
      use parallelisation_layouts, only: vmu_lo, imu_idx, is_idx
      use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx, kxkyz_lo,vmu_lo
    
      ! Dimensions
      use grids_kxky, only: naky, nakx
      use grids_z, only: nztot, ntubes
      use grids_species, only: nspec
      use grids_velocity, only: nvpa, nmu


      ! Flags
      use parameters_physics, only: full_flux_surface
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: include_nonlinear

      !Load data
      use arrays_distribution_function, only: gnew, gvmu
      use arrays_fields, only: phi,apar, bpar
      use arrays, only: shift_state
      use parameters_physics, only: fphi

      ! Redistribute data
      use redistribute, only: scatter
      use initialise_redistribute, only: kxkyz2vmu
      
      !Calculations
      use calculations_gyro_averages, only: gyro_average
      use calculations_tofrom_ghf, only: g_to_h
      use calculations_kxky_derivatives, only: get_dgdy, get_dgdx
      use calculations_kxky_derivatives, only: get_dchidx, get_dchidy
      use calculations_kxky, only: swap_kxky, swap_kxky_back
      use calculations_transforms, only: transform_y2ky, transform_x2kx
      use calculations_transforms, only: transform_y2ky_xfirst, transform_x2kx_xfirst

      !Grids
      use grids_kxky, only: x
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: akx, aky
      use grids_kxky, only: nakx, ikx_max, naky, naky_all, nx, ny      

      !Flow shear
      use gk_flow_shear, only: prp_shear_enabled
      use gk_flow_shear, only: hammett_flow_shear
      use gk_flow_shear, only: g_exb, g_exbfac

      !Geometry
      use geometry, only: exb_nonlin_fac, exb_nonlin_fac_p, gfac
      
      implicit none
      
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :), allocatable :: g0k, g0a
      complex, dimension(:, :), allocatable :: prefac, g0xky
      complex, dimension(:, :), allocatable :: g0kxy
      complex, dimension(:, :), allocatable :: g0k_swap
      complex, dimension(:, :), intent(in out) :: phi_h_nonlin_kx,apar_h_nonlin_kx, bpar_h_nonlin_kx, gint
      complex, dimension(:,:,:,:,:), allocatable :: phi_h_nonlin_k,apar_h_nonlin_k, bpar_h_nonlin_k
      complex, dimension(:,:,:,:), allocatable :: temp
      real, dimension(:, :), allocatable :: g0xy, g1xy
      real, dimension(:, :), allocatable :: bracket_phi_h, bracket_apar_h, bracket_bpar_h
      
      integer :: ikxkyz, iky, ikx
      integer :: ivmu, iz, it, imu, is
      complex, dimension(:, :), allocatable :: field, adjust
      logical :: yfirst
      

      ! By default, prp_shear_enabled = .false. and thus yfirst = .true., and we Fourier transform y first
      ! If perpendicular flow shear is included, it is important to Fourier transform x first
      yfirst = .not. prp_shear_enabled

      ! Allocate arrays
      allocate (g0k(naky, nakx))
      allocate (g0a(naky, nakx))
      allocate (g0xy(ny, nx))
      allocate (g1xy(ny, nx))
      allocate (bracket_phi_h(ny, nx))
      allocate (bracket_apar_h(ny, nx))
      allocate (bracket_bpar_h(ny, nx))
      allocate (prefac(naky, nx))
      allocate (temp(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (phi_h_nonlin_k(naky, nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      allocate (apar_h_nonlin_k(naky, nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      allocate (bpar_h_nonlin_k(naky, nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      ! Here naky_all = 2*naky-1
      if (yfirst) then
         allocate (g0k_swap(naky_all, ikx_max))
         allocate (g0kxy(ny, ikx_max))
      else
         allocate (g0xky(naky, nx))
      end if       

      ! Compute phase factor needed when running with equilibrium flow shear
      prefac = 1.0
      if (prp_shear_enabled .and. hammett_flow_shear) then
         prefac = exp(-zi * g_exb * g_exbfac * spread(x, 1, naky) * spread(aky * shift_state, 2, nx))
      end if      

      ! Incoming pdf is g = <f>.
      ! For EM simulations, the pdf entering the ExB nonlinearity needs to be
      ! the non-Boltzmann part of f (h = f + (Ze/T)*phi*F0)
      if (include_apar .or. include_bpar) call g_to_h(g, phi, bpar, fphi)

      if (include_nonlinear) then
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  ! Compute i*ky*g
                  call get_dgdy(g(:, :, iz, it, ivmu), g0k)
                  ! Take the FFT to get dg/dy in (y,x) space
                  call forward_transform(g0k, g0xy)                  
                  ! Compute i*kx*<phi>
                  call get_dfielddx(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k, 'phi')
                  ! If running with equilibrium flow shear, make adjustment to
                  ! The term multiplying dg/dy
                  if (prp_shear_enabled .and. hammett_flow_shear) then
                     call get_dchidy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0a)
                     g0k = g0k - g_exb * g_exbfac * spread(shift_state, 2, nakx) * g0a
                  end if
                  ! Take the FFT to get d<phi>/dx in (y,x) space
                  call forward_transform(g0k, g1xy)
                  ! Multiply by the geometric factor appearing in the Poisson bracket;
                  ! i.e., (dx/dpsi*dy/dalpha)*0.5
                  g1xy = g1xy * exb_nonlin_fac
                  ! Compute the contribution to the Poisson bracket from dg/dy*d<phi>/dx
                  bracket_phi_h = g0xy * g1xy
                  
                  !repeat the same calculation for apar and bpar terms:
                  call get_dfielddx(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k, 'apar')
                  if (prp_shear_enabled .and. hammett_flow_shear) then
                     g0k = g0k - g_exb * g_exbfac * spread(shift_state, 2, nakx) * g0a
                  end if
                  call forward_transform(g0k, g1xy)
                  g1xy = g1xy * exb_nonlin_fac
                  bracket_apar_h = g0xy * g1xy
                  
                  call get_dfielddx(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k, 'bpar')
                  if (prp_shear_enabled .and. hammett_flow_shear) then
                     g0k = g0k - g_exb * g_exbfac * spread(shift_state, 2, nakx) * g0a
                  end if
		  call forward_transform(g0k, g1xy)
                  g1xy = g1xy * exb_nonlin_fac
                  bracket_bpar_h = g0xy * g1xy

                  ! Compute dg/dx in k-space (= i*kx*g)
                  call get_dgdx(g(:, :, iz, it, ivmu), g0k)
                  ! If running with equilibrium flow shear, correct dg/dx term
                  if (prp_shear_enabled .and. hammett_flow_shear) then
                     call get_dgdy(g(:, :, iz, it, ivmu), g0a)
                     g0k = g0k - g_exb * g_exbfac * spread(shift_state, 2, nakx) * g0a
                  end if
                  ! Take the FFT to get dg/dx in (y,x) space
                  call forward_transform(g0k, g0xy)
                  ! Compute d<phi>/dy in k-space
                  call get_dfielddy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k, 'phi')
                  ! Take the FFT to get d<chi>/dy in (y,x) space
                  call forward_transform(g0k, g1xy)
                  ! Multiply by the geometric factor appearing in the Poisson bracket
                  ! i.e., (dx/dpsi*dy/dalpha)*0.5
                  g1xy = g1xy * exb_nonlin_fac
                  ! Compute the contribution to the Poisson bracket from dg/dy*d<phi>/dx
                  bracket_phi_h = bracket_phi_h - g0xy * g1xy

                  !repeat the same calculation for apar and bpar:
                  call get_dfielddy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k, 'apar')
                  call forward_transform(g0k, g1xy)
                  g1xy = g1xy * exb_nonlin_fac
                  bracket_apar_h = bracket_apar_h - g0xy * g1xy

                  call get_dfielddy(iz, ivmu, phi(:, :, iz, it), apar(:, :, iz, it), bpar(:, :, iz, it), g0k, 'bpar')
      	      	  call forward_transform(g0k, g1xy)
                  g1xy = g1xy * exb_nonlin_fac
                  bracket_bpar_h = bracket_bpar_h - g0xy * g1xy

                  !transform the brackets back to Fourier space:
                  if (yfirst) then
                     call transform_x2kx(bracket_phi_h, g0kxy)
                     call transform_y2ky(g0kxy, g0k_swap)
                     call swap_kxky_back(g0k_swap, phi_h_nonlin_k(:, :, iz, it, ivmu))

                     call transform_x2kx(bracket_apar_h, g0kxy)
                     call transform_y2ky(g0kxy, g0k_swap)
                     call swap_kxky_back(g0k_swap, apar_h_nonlin_k(:, :, iz, it, ivmu))

                     call transform_x2kx(bracket_bpar_h, g0kxy)
                     call transform_y2ky(g0kxy, g0k_swap)
                     call swap_kxky_back(g0k_swap, bpar_h_nonlin_k(:, :, iz, it, ivmu))                     
                  else
                     call transform_y2ky_xfirst(bracket_phi_h, g0xky)
                     g0xky = g0xky / prefac
                     call transform_x2kx_xfirst(g0xky, phi_h_nonlin_k(:, :, iz, it, ivmu))

                     call transform_y2ky_xfirst(bracket_apar_h, g0xky)
                     g0xky = g0xky / prefac
                     call transform_x2kx_xfirst(g0xky, apar_h_nonlin_k(:, :, iz, it, ivmu))

                     call transform_y2ky_xfirst(bracket_bpar_h, g0xky)
                     g0xky = g0xky / prefac
                     call transform_x2kx_xfirst(g0xky, bpar_h_nonlin_k(:, :, iz, it, ivmu))                     
                  end if
               end do
            end do
         end do

         !apply QN operator on the brackets: gyroaverage followed by velocity integration:
         call scatter(kxkyz2vmu, phi_h_nonlin_k, gvmu)
         call integrate(gvmu, temp)
         !take the ky = 0 component:
         phi_h_nonlin_kx(:,:) = temp(1,:,:,1)
         
         !do the same for other brackets:
         call scatter(kxkyz2vmu, apar_h_nonlin_k, gvmu)
         call integrate(gvmu, temp)
         apar_h_nonlin_kx(:,:) = temp(1,:,:,1)

         call scatter(kxkyz2vmu, bpar_h_nonlin_k, gvmu)
         call integrate(gvmu, temp)
         bpar_h_nonlin_kx(:,:) = temp(1,:,:,1)
      else
         phi_h_nonlin_kx = 0.
         apar_h_nonlin_kx = 0.
         bpar_h_nonlin_kx = 0.
      end if

      if (include_apar .or. include_bpar) call g_to_h(g, phi, bpar, -fphi)

      !calculate QN operator acting on g. This is mainly for testing the diagnostics.
      call scatter(kxkyz2vmu, g, gvmu)
      call integrate(gvmu, temp)
      gint(:,:) = temp(1,:,:,1)

      deallocate(g0k, g0a, g0xy, g1xy, bracket_phi_h, bracket_apar_h, bracket_bpar_h)
      deallocate(temp, prefac, phi_h_nonlin_k, apar_h_nonlin_k, bpar_h_nonlin_k)


     contains

  subroutine get_dfielddx(iz, ivmu, phi, apar, bpar, dfielddx, field_flag)      

     ! Constants
     use constants, only: zi
      
     ! Parallelisation
     use parallelisation_layouts, only: vmu_lo
     use parallelisation_layouts, only: is_idx, iv_idx, imu_idx
      
     ! Calculations
     use calculations_gyro_averages, only: gyro_average
     use calculations_gyro_averages, only: gyro_average_j1
      
     ! Flags
     use parameters_physics, only: fphi
      
     ! Grids
     use grids_species, only: spec
     use grids_velocity, only: vpa, mu
     use grids_kxky, only: naky, nakx, akx

     implicit none

     ! Arguments
     integer, intent(in) :: ivmu, iz
     complex, dimension(:, :), intent(in) :: phi, apar, bpar  
     complex, dimension(:, :), intent(out) :: dfielddx
     character(len=*), intent(in) :: field_flag

     ! Local variables
     integer :: iv, is, imu
     complex, dimension(:, :), allocatable :: field, gyro_tmp

     !-------------------------------------------------------------------------

     ! Allocate temporary arrays
     allocate (field(naky, nakx))
     allocate (gyro_tmp(naky, nakx))

     ! Get the (mu,vpa,s) point
     is = is_idx(vmu_lo, ivmu)
     iv = iv_idx(vmu_lo, ivmu)
     imu = imu_idx(vmu_lo, ivmu)

     ! Calculate [i kx phi(ky,kx)]
     if (field_flag == 'phi') then
        field = fphi * zi * spread(akx, 1, naky) * phi
     else if (field_flag == 'apar') then
        field = - zi * spread(akx, 1, naky) * 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
     else if (field_flag == 'bpar') then
        field = zi * spread(akx, 1, naky) * 4 * mu(imu) * (spec(is)%tz) * bpar
     end if

     if (field_flag == 'bpar') then
        call gyro_average_j1(field, iz, ivmu, dfielddx)
     else
        call gyro_average(field, iz, ivmu, dfielddx)
     end if

     ! Deallocate temporary arrays
     deallocate (field)
     deallocate (gyro_tmp)

  end subroutine get_dfielddx

  subroutine get_dfielddy(iz, ivmu, phi, apar, bpar, dfielddy, field_flag)

     ! Constants
     use constants, only: zi

     ! Parallelisation
     use parallelisation_layouts, only: vmu_lo
     use parallelisation_layouts, only: is_idx, iv_idx, imu_idx

     ! Calculations
     use calculations_gyro_averages, only: gyro_average
     use calculations_gyro_averages, only: gyro_average_j1

     ! Flags
     use parameters_physics, only: fphi

     ! Grids
     use grids_species, only: spec
     use grids_velocity, only: vpa, mu
     use grids_kxky, only: nakx, naky, aky

     implicit none

     ! Arguments
     integer, intent(in) :: ivmu, iz
     complex, dimension(:, :), intent(in) :: phi, apar, bpar
     complex, dimension(:, :), intent(out) :: dfielddy

     ! Local variables
     integer :: iv, is, imu
     complex, dimension(:, :), allocatable :: field, gyro_tmp
     character(len=*), intent(in) :: field_flag

     !-------------------------------------------------------------------------

     ! Allocate temporary arrays
     allocate (field(naky, nakx))
     allocate (gyro_tmp(naky, nakx))

     ! Get the (mu,vpa,s) point
     is = is_idx(vmu_lo, ivmu)
     iv = iv_idx(vmu_lo, ivmu)
     imu = imu_idx(vmu_lo, ivmu)

     ! Calculate [i ky phi(ky,kx)]
     if (field_flag == 'phi') then
        field = zi * spread(aky, 2, nakx) * fphi * phi
     else if (field_flag == 'apar') then
        field =  - zi * spread(aky, 2, nakx) * 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
     else if (field_flag == 'bpar') then
        field = zi * spread(aky, 2, nakx) * 4.0 * mu(imu) * (spec(is)%tz) * bpar
     end if

     if (field_flag == 'bpar') then
        call gyro_average_j1(field, iz, ivmu, dfielddy)
     else
        call gyro_average(field, iz, ivmu, dfielddy)
     end if

     ! Deallocate temporary arrays
     deallocate (field)
     deallocate (gyro_tmp)

  end subroutine get_dfielddy


  subroutine forward_transform(gk, gx)
 
     use calculations_transforms, only: transform_ky2y, transform_kx2x
     use calculations_transforms, only: transform_ky2y_xfirst, transform_kx2x_xfirst

     implicit none

     ! Arguments
     complex, dimension(:, :), intent(in) :: gk
     real, dimension(:, :), intent(out) :: gx
     
     !----------------------------------------------------------------------

     ! Fourier transform g(kx,ky) to g(kx,y) and then to g(x,y)
     if (yfirst) then
         ! We have i*ky*g(kx,ky) for ky >= 0 and all kx.
         ! We want to do 1D complex to complex transform in y,
         ! which requires i*ky*g(kx,ky) for all ky and kx >= 0 .
	 ! Use the reality condition: g(kx,-ky) = conjg(g(-kx,ky))
         ! so i*(-ky)*g(kx,-ky) = -i*ky*conjg(g(-kx,ky)) = conjg(i*ky*g(-kx,ky))
         ! and i*kx*g(kx,-ky) = i*kx*conjg(g(-kx,ky)) = conjg(i*(-kx)*g(-kx,ky))
         ! and i*(-ky)*J0(kx,-ky)*phi(kx,-ky) = conjg(i*ky*J0(-kx,ky)*phi(-kx,ky))
         ! and i*kx*J0(kx,-ky)*phi(kx,-ky) = conjg(i*(-kx)*J0(-kx,ky)*phi(-kx,ky))
         ! i.e., can calculate dg/dx, dg/dy, d<phi>/dx and d<phi>/dy
         ! on stella (kx,ky) grid, then conjugate and flip sign of (kx,ky)
         ! NB: J0(kx,ky) = J0(-kx,-ky)
         ! TODO DSO: coordinate change for shearing
         call swap_kxky(gk, g0k_swap)
         call transform_ky2y(g0k_swap, g0kxy)
         call transform_kx2x(g0kxy, gx)

     ! Fourier transform g(kx,ky) to g(x,ky) and then to g(x,y)
     else
         call transform_kx2x_xfirst(gk, g0xky)
         g0xky = g0xky * prefac
         call transform_ky2y_xfirst(g0xky, gx)
     end if

  end subroutine forward_transform


  subroutine integrate(fv, fint)


     ! Parallelisation
     use mp, only: proc0, sum_allreduce   
     use parallelisation_layouts, only: vmu_lo, imu_idx, is_idx
     use parallelisation_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx, kxkyz_lo,vmu_lo

     !Grids
     use grids_z, only: nzgrid, ntubes
     use grids_kxky, only: akx, aky
     use grids_kxky, only: nakx, ikx_max, naky, naky_all, nx, ny
     use grids_species, only: spec, nspec
     use grids_velocity, only: nvpa, nmu

     !calculations
     use constants, only: zi
     use calculations_velocity_integrals, only: integrate_vmu
     use calculations_gyro_averages, only: gyro_average

     implicit none

     complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: fv
     complex, dimension(:, :, -nzgrid:, :), intent(out) :: fint
   
     !local variables
     complex :: tmp
     real :: wgt
     complex, dimension(:, :), allocatable :: g0
     integer :: ikxkyz, iky, ikx, iz, it, is, ia

     allocate (g0(nvpa, nmu))
     fint = 0.
     do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
        iz = iz_idx(kxkyz_lo, ikxkyz)
        it = it_idx(kxkyz_lo, ikxkyz)
        ikx = ikx_idx(kxkyz_lo, ikxkyz)
        iky = iky_idx(kxkyz_lo, ikxkyz)
        is = is_idx(kxkyz_lo, ikxkyz)
        call gyro_average(fv(:, :, ikxkyz), ikxkyz, g0)
        wgt = spec(is)%z * spec(is)%dens_psi0
        call integrate_vmu(g0, iz, tmp)
        fint(iky, ikx, iz, it) = fint(iky, ikx, iz, it) + wgt * tmp
     end do
     call sum_allreduce(fint)
     deallocate (g0)
  end subroutine integrate

end subroutine calculate_turbulent_stresses
  
end module diagnostics_stresses


