# include "define.inc"

module stella_io

# ifdef NETCDF
   use netcdf, only: nf90_noerr
   use netcdf_utils, only: netcdf_error, kind_nf
# endif

   implicit none

   private

   public :: init_stella_io, finish_stella_io
   public :: write_time_nc
   public :: write_phi2_nc
   public :: write_phi_nc
   public :: write_gvmus_nc
   public :: write_gzvs_nc
   public :: write_kspectra_nc
   public :: write_omega_nc
   public :: write_moments_nc
   public :: write_radial_fluxes_nc
   public :: write_radial_moments_nc
   public :: write_fluxes_kxkyz_nc
   public :: get_nout
   public :: sync_nc

# ifdef NETCDF
   integer(kind_nf) :: ncid
   integer(kind_nf) :: naky_dim, nttot_dim, nmu_dim, nvtot_dim, nspec_dim
   integer(kind_nf) :: nakx_dim, ntubes_dim, radgridvar_dim
   integer(kind_nf) :: time_dim, char10_dim, char200_dim, ri_dim, nlines_dim, nheat_dim
   integer(kind_nf) :: nalpha_dim

   integer, dimension(7) :: moment_dim
   integer, dimension(6) :: field_dim
   integer, dimension(5) :: zvs_dim
   integer, dimension(4) :: om_dim, vmus_dim, kykxaz_dim
   integer, dimension(6) :: flx_dim
   integer, dimension(3) :: mode_dim, heat_dim, kykxz_dim, flux_x_dim
   integer, dimension(2) :: kx_dim, ky_dim, flux_dim, nin_dim, fmode_dim
   integer, dimension(2) :: flux_surface_dim, rad_grid_dim

   integer :: time_id, phi2_id
   integer :: phi_vs_t_id, phi2_vs_kxky_id
   integer :: dens_x_id, upar_x_id, temp_x_id
   integer :: pflux_x_id, vflux_x_id, qflux_x_id
   integer :: pflx_kxkyz_id, vflx_kxkyz_id, qflx_kxkyz_id
   integer :: density_id, upar_id, temperature_id, spitzer2_id
   integer :: omega_id
   integer :: gvmus_id, gzvs_id
   integer :: input_id
   integer :: bmag_id, gradpar_id, gbdrift_id, gbdrift0_id, b_dot_grad_z_id
   integer :: cvdrift_id, cvdrift0_id, gds2_id, gds21_id, gds22_id
   integer :: kperp2_id
   integer :: grho_id, jacob_id, djacdrho_id, shat_id, drhodpsi_id, q_id, jtwist_id
   integer :: d2qdr2_id, d2psidr2_id
   integer :: beta_id
   integer :: code_id
# endif

   real :: zero

!  include 'netcdf.inc'

contains

   !==============================================
   !============ INITIATE STELLA IO ==============
   !==============================================
   subroutine init_stella_io(restart, write_phi_vs_t, write_kspectra, write_gvmus, &
                             write_gzvs, write_moments, write_omega, write_radial_fluxes, write_radial_moments, write_fluxes_kxky)

      use mp, only: proc0
      use file_utils, only: run_name

# ifdef NETCDF
      use netcdf_utils, only: get_netcdf_code_precision, netcdf_real
      use neasyf, only: neasyf_open, neasyf_metadata
      use git_version, only: get_git_version
# endif

      implicit none

      logical, intent(in) :: restart
      logical, intent(in) :: write_phi_vs_t, write_kspectra, write_gvmus, write_gzvs
      logical, intent(in) :: write_moments, write_omega, write_radial_fluxes, write_radial_moments!, write_symmetry
      logical, intent(in) :: write_fluxes_kxky
# ifdef NETCDF
      character(300) :: filename

      zero = epsilon(0.0)

      if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()

      ! The netcdf file has the extension ".out.nc"
      filename = trim(trim(run_name)//'.out.nc')

      ! Only the first processor (proc0) opens the file
      if (proc0) then
         if (restart) then
            ncid = neasyf_open(trim(filename), "rw")
         else
            ncid = neasyf_open(trim(filename), "w")
         end if

         call neasyf_metadata(ncid, title="stella simulation data", software_name="stella", &
                              software_version=get_git_version(), auto_date=.true.)
         call write_grids(ncid)
         call define_vars(write_phi_vs_t, write_kspectra, write_gvmus, write_gzvs, &
                          write_moments, write_omega, write_radial_fluxes, write_radial_moments, &
                          write_fluxes_kxky)
         call nc_species(ncid)
         call nc_geo
         call save_input
      end if
# endif

   end subroutine init_stella_io

   !> Ensure the netCDF file contains all the dimensions and grids,
   !> creating them if necessary
   subroutine write_grids(file_id)
# ifdef NETCDF
      use file_utils, only: num_input_lines
      use kt_grids, only: nakx, naky, akx, aky, nalpha, theta0, phase_shift_angle, x_d, rho_d
      use zgrid, only: nzgrid, ntubes, zed
      use vpamu_grids, only: nvpa, vpa, nmu, mu
      use species, only: nspec
      use physics_flags, only: radial_variation
      use physics_parameters, only: rhostar
      use stella_geometry, only: geo_surf, dxdXcoord, q_as_x
      use mp, only: nproc
      use neasyf, only: neasyf_dim, neasyf_write
# endif

      !> NetCDF ID of the file
      integer, intent(in) :: file_id

# ifdef NETCDF
      integer :: ix
      ! Total number of mesh points
      real :: nmesh
      ! Radial grid
      real, dimension(:, :), allocatable :: rg

      ! Grids themselves
      call neasyf_dim(file_id, "ky", values=aky, dimid=naky_dim, &
                      long_name="Wavenumber perpendicular to flux surface", units="1/rho_r")
      call neasyf_dim(file_id, "kx", values=akx, dimid=nakx_dim, &
                      long_name="Wavenumber in direction of grad alpha", units="1/rho_r")
      call neasyf_dim(file_id, "tube", dim_size=ntubes, dimid=ntubes_dim)
      call neasyf_dim(file_id, "zed", values=zed, dimid=nttot_dim)
      call neasyf_dim(file_id, "alpha", dim_size=nalpha, dimid=nalpha_dim)
      call neasyf_dim(file_id, "vpa", values=vpa, dimid=nvtot_dim)
      call neasyf_dim(file_id, "mu", values=mu, dimid=nmu_dim)
      call neasyf_dim(file_id, "species", dim_size=nspec, dimid=nspec_dim)
      call neasyf_dim(file_id, "t", unlimited=.true., dimid=time_dim, &
                      long_name="Time", units="L/vt")

      ! Dimensions for various string variables
      call neasyf_dim(file_id, "char10", dim_size=10, dimid=char10_dim)
      call neasyf_dim(file_id, "char200", dim_size=200, dimid=char200_dim)
      if (num_input_lines > 0) then
         ! netCDF interprets zero-sized dimensions as unlimited, which is not what we want
         call neasyf_dim(file_id, "nlines", dim_size=num_input_lines, dimid=nlines_dim, long_name="Input file line number")
      end if

      call neasyf_dim(file_id, "ri", dim_size=2, dimid=ri_dim, &
                      long_name="Complex components", units="(real, imaginary)")

      if (radial_variation) then
         call neasyf_dim(file_id, "radgridvar", dim_size=3, dimid=radgridvar_dim, &
                         long_name="x, q/psi, rho")
      end if

      ! Grid sizes
      call neasyf_write(file_id, "nkx", nakx, long_name="Number of kx points")
      call neasyf_write(file_id, "nky", naky, long_name="Number of ky points")
      call neasyf_write(file_id, "ntubes", ntubes, long_name="Number of tubes")
      call neasyf_write(file_id, "nzed_tot", 2 * nzgrid + 1, long_name="Number of zed points")
      call neasyf_write(file_id, "theta0", theta0, dim_names=["ky", "kx"], long_name="Theta_0")
      call neasyf_write(file_id, "nspecies", nspec, long_name="Number of species")
      call neasyf_write(file_id, "nvpa_tot", nvpa, long_name="Number of vpa points")
      call neasyf_write(file_id, "nmu", nmu, long_name="Number of mu points")
      call neasyf_write(file_id, "phase_shift_angle", phase_shift_angle)

      call neasyf_write(file_id, "nproc", nproc, long_name="Number of processors")

      if (radial_variation) then
         allocate (rg(3, nakx))
         if (q_as_x) then
            do ix = 1, nakx
               rg(1, ix) = x_d(ix)
               rg(2, ix) = rhostar * x_d(ix) / dxdXcoord + geo_surf%qinp
               rg(3, ix) = rho_d(ix) + geo_surf%rhoc
            end do
         else
            do ix = 1, nakx
               rg(1, ix) = x_d(ix)
               rg(2, ix) = rhostar * x_d(ix) / dxdXcoord
               rg(3, ix) = rho_d(ix) + geo_surf%rhoc
            end do
         end if
         call neasyf_write(file_id, "rad_grid", rg, dim_names=[character(10)::"radgridvar", "kx"])
         deallocate (rg)
      end if

      nmesh = (2 * nzgrid + 1) * ntubes * nvpa * nmu * nakx * naky * nspec
      call neasyf_write(file_id, "nmesh", nmesh, long_name="Number of mesh points")
# endif
   end subroutine write_grids

   subroutine finish_stella_io
      use mp, only: proc0
# ifdef NETCDF
      use netcdf, only: nf90_close
      use netcdf_utils, only: netcdf_error

      integer :: status

      if (proc0) then
         status = nf90_close(ncid)
         if (status /= nf90_noerr) call netcdf_error(status)
      end if
# endif
   end subroutine finish_stella_io

   !> Save the input file in the NetCDF file
   subroutine save_input
# ifdef NETCDF
      use file_utils, only: num_input_lines, get_input_unit
      use netcdf, only: nf90_put_var

      character(200) line
      integer, dimension(2) :: nin_start, nin_count

      integer :: status, n, unit

      nin_start(1) = 1
      nin_start(2) = 1

      nin_count(2) = 1

      call get_input_unit(unit)
      rewind (unit=unit)
      do n = 1, num_input_lines
         read (unit=unit, fmt="(a)") line
         nin_count(1) = len(trim(line))
!       status = nf_put_vara_text (ncid, input_id, nin_start, nin_count, line)
         status = nf90_put_var(ncid, input_id, line, start=nin_start, count=nin_count)
         if (status /= nf90_noerr) call netcdf_error(status, ncid, input_id)
         nin_start(2) = nin_start(2) + 1
      end do
# endif
   end subroutine save_input

   subroutine define_vars(write_phi_vs_t, write_kspectra, write_gvmus, &
                          !       write_gzvs, write_symmetry, write_moments)
                          write_gzvs, write_moments, write_omega, write_radial_fluxes, &
                          write_radial_moments, write_fluxes_kxky)

      use run_parameters, only: fphi!, fapar, fbpar
      use physics_flags, only: radial_variation
# ifdef NETCDF
      use netcdf, only: nf90_char, nf90_int, nf90_global
      use netcdf, only: nf90_def_var, nf90_inq_varid, nf90_put_att, nf90_enddef, nf90_put_var
      use netcdf, only: nf90_inq_libvers
      use netcdf_utils, only: netcdf_real
# endif

      implicit none

      logical, intent(in) :: write_phi_vs_t, write_kspectra, write_gvmus, write_gzvs!, write_symmetry
      logical, intent(in) :: write_moments, write_omega, write_radial_fluxes, write_radial_moments
      logical, intent(in) :: write_fluxes_kxky
# ifdef NETCDF

      character(5) :: ci
      character(20) :: datestamp, timestamp, timezone

      integer :: status

      flux_surface_dim(1) = nalpha_dim
      flux_surface_dim(2) = nttot_dim

      fmode_dim(1) = naky_dim
      fmode_dim(2) = nakx_dim

      mode_dim(1) = naky_dim
      mode_dim(2) = nakx_dim
      mode_dim(3) = time_dim

      flx_dim(1) = naky_dim
      flx_dim(2) = nakx_dim
      flx_dim(3) = nttot_dim
      flx_dim(4) = ntubes_dim
      flx_dim(5) = nspec_dim
      flx_dim(6) = time_dim

      kx_dim(1) = nakx_dim
      kx_dim(2) = time_dim

      ky_dim(1) = naky_dim
      ky_dim(2) = time_dim

      om_dim(1) = ri_dim
      om_dim(2) = naky_dim
      om_dim(3) = nakx_dim
      om_dim(4) = time_dim

      nin_dim(1) = char200_dim
      nin_dim(2) = nlines_dim

      flux_dim(1) = nspec_dim
      flux_dim(2) = time_dim

      flux_x_dim(1) = nakx_dim
      flux_x_dim(2) = nspec_dim
      flux_x_dim(3) = time_dim

      rad_grid_dim(1) = radgridvar_dim
      rad_grid_dim(2) = nakx_dim

      heat_dim(1) = nspec_dim
      heat_dim(2) = nheat_dim
      heat_dim(3) = time_dim

      field_dim(1) = ri_dim
      field_dim(2) = naky_dim
      field_dim(3) = nakx_dim
      field_dim(4) = nttot_dim
      field_dim(5) = ntubes_dim
      field_dim(6) = time_dim

      moment_dim(1) = ri_dim
      moment_dim(2) = naky_dim
      moment_dim(3) = nakx_dim
      moment_dim(4) = nttot_dim
      moment_dim(5) = ntubes_dim
      moment_dim(6) = nspec_dim
      moment_dim(7) = time_dim

      vmus_dim(1) = nvtot_dim
      vmus_dim(2) = nmu_dim
      vmus_dim(3) = nspec_dim
      vmus_dim(4) = time_dim

      zvs_dim(1) = ntubes_dim
      zvs_dim(2) = nttot_dim
      zvs_dim(3) = nvtot_dim
      zvs_dim(4) = nspec_dim
      zvs_dim(5) = time_dim

      kykxz_dim(1) = naky_dim
      kykxz_dim(2) = nakx_dim
      kykxz_dim(3) = nttot_dim

      kykxaz_dim(1) = naky_dim
      kykxaz_dim(2) = nakx_dim
      kykxaz_dim(3) = nalpha_dim
      kykxaz_dim(4) = nttot_dim

      ! Write some useful general information such as the website,
      ! date and time into the NetCDF file

      datestamp(:) = ' '
      timestamp(:) = ' '
      timezone(:) = ' '
      call date_and_time(datestamp, timestamp, timezone)

      status = nf90_inq_varid(ncid, 'code_info', code_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'code_info', nf90_char, char10_dim, code_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='code_info')
      end if
      status = nf90_put_att(ncid, code_id, 'long_name', 'stella')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att='long_name')

      ci = 'c1'
      status = nf90_put_att(ncid, code_id, trim(ci), 'Date: '//trim(datestamp))
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

      ci = 'c2'
      status = nf90_put_att(ncid, code_id, trim(ci), &
                            'Time: '//trim(timestamp)//' '//trim(timezone))
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

      ci = 'c3'
      status = nf90_put_att(ncid, code_id, trim(ci), &
                            'netCDF version '//trim(nf90_inq_libvers()))
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

      ci = 'c4'
      status = nf90_put_att(ncid, code_id, trim(ci), &
                            'Units are determined with respect to reference temperature (T_ref),')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

      ci = 'c5'
      status = nf90_put_att(ncid, code_id, trim(ci), &
                            'reference charge (q_ref), reference mass (mass_ref),')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

      ci = 'c6'
      status = nf90_put_att(ncid, code_id, trim(ci), &
                            'reference field (B_ref), and reference length (a_ref)')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

      ci = 'c7'
      status = nf90_put_att(ncid, code_id, trim(ci), &
                            'from which one may construct rho_ref and vt_ref/a,')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

      ci = 'c8'
      status = nf90_put_att(ncid, code_id, trim(ci), &
                            'which are the basic units of perpendicular length and time.')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

      ci = 'c9'
      status = nf90_put_att(ncid, code_id, trim(ci), &
                            'Macroscopic lengths are normalized to the minor radius.')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

      ci = 'c10'
      status = nf90_put_att(ncid, code_id, trim(ci), &
                            'The difference between rho (normalized minor radius) and rho (gyroradius)')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

      ci = 'c11'
      status = nf90_put_att(ncid, code_id, trim(ci), &
                            'should be clear from the context in which they appear below.')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

      status = nf90_inq_varid(ncid, 'bmag', bmag_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'bmag', netcdf_real, flux_surface_dim, bmag_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='bmag')
      end if
      status = nf90_put_att(ncid, bmag_id, 'long_name', '|B|(alpha,zed)')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, bmag_id, att='long_name')
      status = nf90_put_att(ncid, bmag_id, 'units', 'B_0')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, bmag_id, att='units')

      status = nf90_inq_varid(ncid, 'b_dot_grad_z', b_dot_grad_z_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'b_dot_grad_z', netcdf_real, flux_surface_dim, b_dot_grad_z_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='b_dot_grad_z')
      end if
      status = nf90_put_att(ncid, b_dot_grad_z_id, 'long_name', 'b . grad z(alpha,zed)')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, b_dot_grad_z_id, att='long_name')
      status = nf90_put_att(ncid, b_dot_grad_z_id, 'units', 'dimensionless')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, b_dot_grad_z_id, att='units')

      status = nf90_inq_varid(ncid, 'gradpar', gradpar_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'gradpar', netcdf_real, nttot_dim, gradpar_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='gradpar')
      end if
      status = nf90_inq_varid(ncid, 'gbdrift', gbdrift_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'gbdrift', netcdf_real, flux_surface_dim, gbdrift_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='gbdrift')
      end if
      status = nf90_inq_varid(ncid, 'gbdrift0', gbdrift0_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'gbdrift0', netcdf_real, flux_surface_dim, gbdrift0_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='gbdrift0')
      end if
      status = nf90_inq_varid(ncid, 'cvdrift', cvdrift_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'cvdrift', netcdf_real, flux_surface_dim, cvdrift_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='cvdrift')
      end if
      status = nf90_inq_varid(ncid, 'cvdrift0', cvdrift0_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'cvdrift0', netcdf_real, flux_surface_dim, cvdrift0_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='cvdrift0')
      end if

      status = nf90_inq_varid(ncid, 'kperp2', kperp2_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'kperp2', netcdf_real, kykxaz_dim, kperp2_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='kperp2')
      end if
      status = nf90_inq_varid(ncid, 'gds2', gds2_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'gds2', netcdf_real, flux_surface_dim, gds2_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='gds2')
      end if
      status = nf90_inq_varid(ncid, 'gds21', gds21_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'gds21', netcdf_real, flux_surface_dim, gds21_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='gds21')
      end if
      status = nf90_inq_varid(ncid, 'gds22', gds22_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'gds22', netcdf_real, flux_surface_dim, gds22_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='gds22')
      end if
      status = nf90_inq_varid(ncid, 'grho', grho_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'grho', netcdf_real, flux_surface_dim, grho_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='grho')
      end if
      status = nf90_inq_varid(ncid, 'jacob', jacob_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'jacob', netcdf_real, flux_surface_dim, jacob_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='jacob')
      end if
      status = nf90_inq_varid(ncid, 'djacdrho', djacdrho_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'djacdrho', netcdf_real, flux_surface_dim, djacdrho_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='djacdrho')
      end if

      status = nf90_inq_varid(ncid, 'q', q_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'q', netcdf_real, q_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='q')
      end if
      status = nf90_put_att(ncid, q_id, 'long_name', 'local safety factor')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, q_id, att='long_name')
      status = nf90_inq_varid(ncid, 'beta', beta_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'beta', netcdf_real, beta_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='beta')
      end if
      status = nf90_put_att(ncid, beta_id, 'long_name', 'reference beta')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, beta_id, att='long_name')
      status = nf90_inq_varid(ncid, 'shat', shat_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'shat', netcdf_real, shat_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='shat')
      end if
      status = nf90_put_att(ncid, shat_id, 'long_name', '(rho/q) dq/drho')
      status = nf90_inq_varid(ncid, 'd2qdr2', d2qdr2_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'd2qdr2', netcdf_real, d2qdr2_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='d2qdr2')
      end if
      if (status /= nf90_noerr) call netcdf_error(status, ncid, shat_id, att='long_name')
      status = nf90_inq_varid(ncid, 'jtwist', jtwist_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'jtwist', netcdf_real, jtwist_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='jtwist')
      end if
      status = nf90_put_att(ncid, jtwist_id, 'long_name', '2*pi*shat*dky/dkx')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, jtwist_id, att='long_name')

      status = nf90_inq_varid(ncid, 'drhodpsi', drhodpsi_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'drhodpsi', netcdf_real, drhodpsi_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='drhodpsi')
      end if
      status = nf90_put_att(ncid, drhodpsi_id, 'long_name', 'drho/dPsi')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, drhodpsi_id, att='long_name')
      status = nf90_inq_varid(ncid, 'd2psidr2', d2psidr2_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'd2psidr2', netcdf_real, d2psidr2_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='d2psidr2')
      end if

      if (write_omega) then
         status = nf90_inq_varid(ncid, 'omega', omega_id)
         if (status /= nf90_noerr) then
            status = nf90_def_var(ncid, 'omega', netcdf_real, om_dim, omega_id)
            if (status /= nf90_noerr) call netcdf_error(status, var='omega')
         end if
      end if

      if (fphi > zero) then
         status = nf90_inq_varid(ncid, 'phi2', phi2_id)
         if (status /= nf90_noerr) then
            status = nf90_def_var(ncid, 'phi2', netcdf_real, time_dim, phi2_id)
            if (status /= nf90_noerr) call netcdf_error(status, var='phi2')
         end if
         status = nf90_put_att(ncid, phi2_id, 'long_name', '|Potential**2|')
         if (status /= nf90_noerr) &
            call netcdf_error(status, ncid, phi2_id, att='long_name')
         status = nf90_put_att(ncid, phi2_id, 'units', '(T/q rho/L)**2')
         if (status /= nf90_noerr) &
            call netcdf_error(status, ncid, phi2_id, att='units')

!        status = nf90_def_var &
!             (ncid, 'phi2_by_mode', netcdf_real, mode_dim, phi2_by_mode_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='phi2_by_mode')
!        if (nakx > 1) then
!           status = nf90_def_var &
!                (ncid, 'phi2_by_kx', netcdf_real, kx_dim, phi2_by_kx_id)
!           if (status /= nf90_noerr) &
!                call netcdf_error (status, var='phi2_by_kx')
!        end if

!        if (naky > 1) then
!           status = nf90_def_var &
!                (ncid, 'phi2_by_ky', netcdf_real, ky_dim, phi2_by_ky_id)
!           if (status /= nf90_noerr) &
!                call netcdf_error (status, var='phi2_by_ky')
!        end if

         if (write_phi_vs_t) then
            status = nf90_inq_varid(ncid, 'phi_vs_t', phi_vs_t_id)
            if (status /= nf90_noerr) then
               status = nf90_def_var &
                        (ncid, 'phi_vs_t', netcdf_real, field_dim, phi_vs_t_id)
               if (status /= nf90_noerr) call netcdf_error(status, var='phi_vs_t')
            end if
            status = nf90_put_att(ncid, phi_vs_t_id, 'long_name', 'Electrostatic Potential vs time')
            if (status /= nf90_noerr) call netcdf_error(status, ncid, phi_vs_t_id, att='long_name')
         end if
         if (write_radial_moments) then
            status = nf90_inq_varid(ncid, 'dens_x', dens_x_id)
            if (status /= nf90_noerr) then
               status = nf90_def_var &
                        (ncid, 'dens_x', netcdf_real, flux_x_dim, dens_x_id)
               if (status /= nf90_noerr) call netcdf_error(status, var='dens_x')
            end if
            status = nf90_inq_varid(ncid, 'upar_x', upar_x_id)
            if (status /= nf90_noerr) then
               status = nf90_def_var &
                        (ncid, 'upar_x', netcdf_real, flux_x_dim, upar_x_id)
               if (status /= nf90_noerr) call netcdf_error(status, var='upar_x')
            end if
            status = nf90_inq_varid(ncid, 'temp_x', temp_x_id)
            if (status /= nf90_noerr) then
               status = nf90_def_var &
                        (ncid, 'temp_x', netcdf_real, flux_x_dim, temp_x_id)
               if (status /= nf90_noerr) call netcdf_error(status, var='temp_x')
            end if
         end if
         if (write_radial_fluxes) then
            status = nf90_inq_varid(ncid, 'pflux_x', pflux_x_id)
            if (status /= nf90_noerr) then
               status = nf90_def_var &
                        (ncid, 'pflux_x', netcdf_real, flux_x_dim, pflux_x_id)
               if (status /= nf90_noerr) call netcdf_error(status, var='pflux_x')
            end if
            status = nf90_inq_varid(ncid, 'vflux_x', vflux_x_id)
            if (status /= nf90_noerr) then
               status = nf90_def_var &
                        (ncid, 'vflux_x', netcdf_real, flux_x_dim, vflux_x_id)
               if (status /= nf90_noerr) call netcdf_error(status, var='vflux_x')
            end if
            status = nf90_inq_varid(ncid, 'qflux_x', qflux_x_id)
            if (status /= nf90_noerr) then
               status = nf90_def_var &
                        (ncid, 'qflux_x', netcdf_real, flux_x_dim, qflux_x_id)
               if (status /= nf90_noerr) call netcdf_error(status, var='qflux_x')
            end if
         end if
         if (write_kspectra) then
            status = nf90_inq_varid(ncid, 'phi2_vs_kxky', phi2_vs_kxky_id)
            if (status /= nf90_noerr) then
               status = nf90_def_var &
                        (ncid, 'phi2_vs_kxky', netcdf_real, mode_dim, phi2_vs_kxky_id)
               if (status /= nf90_noerr) call netcdf_error(status, var='phi2_vs_kxky')
            end if
            status = nf90_put_att(ncid, phi2_vs_kxky_id, 'long_name', 'Electrostatic Potential vs (ky,kx,t)')
            if (status /= nf90_noerr) call netcdf_error(status, ncid, phi2_vs_kxky_id, att='long_name')
         end if
      end if
!
!
!
      if (write_fluxes_kxky) then
         status = nf90_inq_varid(ncid, 'pflx_kxky', pflx_kxkyz_id)
         if (status /= nf90_noerr) then
            status = nf90_def_var &
                     (ncid, 'pflx_kxky', netcdf_real, flx_dim, pflx_kxkyz_id)
            if (status /= nf90_noerr) call netcdf_error(status, var='pflx_kxky')
         end if
         status = nf90_put_att(ncid, pflx_kxkyz_id, 'long_name', 'Particle flux vs (ky,kx,spec,t)')
         if (status /= nf90_noerr) call netcdf_error(status, ncid, pflx_kxkyz_id, att='long_name')
!
         status = nf90_inq_varid(ncid, 'vflx_kxky', vflx_kxkyz_id)
         if (status /= nf90_noerr) then
            status = nf90_def_var &
                     (ncid, 'vflx_kxky', netcdf_real, flx_dim, vflx_kxkyz_id)
            if (status /= nf90_noerr) call netcdf_error(status, var='vflx_kxky')
         end if
         status = nf90_put_att(ncid, vflx_kxkyz_id, 'long_name', 'Momentum flux vs (ky,kx,spec,t)')
         if (status /= nf90_noerr) call netcdf_error(status, ncid, vflx_kxkyz_id, att='long_name')
!
         status = nf90_inq_varid(ncid, 'qflx_kxky', qflx_kxkyz_id)
         if (status /= nf90_noerr) then
            status = nf90_def_var &
                     (ncid, 'qflx_kxky', netcdf_real, flx_dim, qflx_kxkyz_id)
            if (status /= nf90_noerr) call netcdf_error(status, var='qflx_kxky')
         end if
         status = nf90_put_att(ncid, qflx_kxkyz_id, 'long_name', 'Heat flux vs (ky,kx,spec,t)')
         if (status /= nf90_noerr) call netcdf_error(status, ncid, qflx_kxkyz_id, att='long_name')
      end if
!
!
!
      if (write_moments) then
         status = nf90_inq_varid(ncid, 'density', density_id)
         if (status /= nf90_noerr) then
            status = nf90_def_var &
                     (ncid, 'density', netcdf_real, moment_dim, density_id)
            if (status /= nf90_noerr) call netcdf_error(status, var='density')
         end if
         status = nf90_put_att(ncid, density_id, 'long_name', 'perturbed density vs (ky,kx,z,t)')
         if (status /= nf90_noerr) call netcdf_error(status, ncid, density_id, att='long_name')
         status = nf90_inq_varid(ncid, 'upar', upar_id)
         if (status /= nf90_noerr) then
            status = nf90_def_var &
                     (ncid, 'upar', netcdf_real, moment_dim, upar_id)
            if (status /= nf90_noerr) call netcdf_error(status, var='upar')
         end if
         status = nf90_put_att(ncid, upar_id, 'long_name', 'perturbed parallel flow vs (ky,kx,z,t)')
         if (status /= nf90_noerr) call netcdf_error(status, ncid, upar_id, att='long_name')
         status = nf90_inq_varid(ncid, 'temperature', temperature_id)
         if (status /= nf90_noerr) then
            status = nf90_def_var &
                     (ncid, 'temperature', netcdf_real, moment_dim, temperature_id)
            if (status /= nf90_noerr) call netcdf_error(status, var='temperature')
         end if
         status = nf90_put_att(ncid, temperature_id, 'long_name', 'perturbed temperature vs (ky,kx,z,t)')
         if (status /= nf90_noerr) call netcdf_error(status, ncid, temperature_id, att='long_name')

         status = nf90_inq_varid(ncid, 'spitzer2', spitzer2_id)
         if (status /= nf90_noerr) then
            status = nf90_def_var &
                     (ncid, 'spitzer2', netcdf_real, moment_dim, spitzer2_id)
            if (status /= nf90_noerr) call netcdf_error(status, var='spitzer2')
         end if
         status = nf90_put_att(ncid, spitzer2_id, 'long_name', 'integral req. for 2. Spitzer coeff')
         if (status /= nf90_noerr) call netcdf_error(status, ncid, spitzer2_id, att='long_name')
      end if

      if (write_gvmus) then
         status = nf90_inq_varid(ncid, 'gvmus', gvmus_id)
         if (status /= nf90_noerr) then
            status = nf90_def_var &
                     (ncid, 'gvmus', netcdf_real, vmus_dim, gvmus_id)
            if (status /= nf90_noerr) call netcdf_error(status, var='gvmus')
         end if
         status = nf90_put_att(ncid, gvmus_id, 'long_name', &
                               'guiding center distribution function averaged over real space')
         if (status /= nf90_noerr) call netcdf_error(status, ncid, gvmus_id, att='long_name')
      end if

      if (write_gzvs) then
         status = nf90_inq_varid(ncid, 'gzvs', gzvs_id)
         if (status /= nf90_noerr) then
            status = nf90_def_var &
                     (ncid, 'gzvs', netcdf_real, zvs_dim, gzvs_id)
            if (status /= nf90_noerr) call netcdf_error(status, var='gzvs')
         end if
         status = nf90_put_att(ncid, gvmus_id, 'long_name', &
                               'guiding center distribution function averaged over (kx,ky,mu)')
         if (status /= nf90_noerr) call netcdf_error(status, ncid, gzvs_id, att='long_name')
      end if

!    if (write_symmetry) then
!        status = nf90_def_var &
!             (ncid, 'pflx_zvpa', netcdf_real, zvs_dim, gzvs_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='gzvs')
!        status = nf90_put_att (ncid, gvmus_id, 'long_name', &
!             'guiding center distribution function averaged over (kx,ky,mu)')
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, gzvs_id, att='long_name')
!    end if

      status = nf90_inq_varid(ncid, 'input_file', input_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'input_file', nf90_char, nin_dim, input_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='input_file')
      end if
      status = nf90_put_att(ncid, input_id, 'long_name', 'Input file')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, input_id, att='long_name')

      status = nf90_enddef(ncid)  ! out of definition mode
      if (status /= nf90_noerr) call netcdf_error(status)

# endif
   end subroutine define_vars

   subroutine write_time_nc(nout, time)

# ifdef NETCDF
      use netcdf, only: nf90_put_var, nf90_sync
# endif

      implicit none

      integer, intent(in) :: nout
      real, intent(in) :: time

# ifdef NETCDF
      integer :: status

      status = nf90_put_var(ncid, time_id, time, start=(/nout/))
      if (status /= nf90_noerr) call netcdf_error(status, ncid, time_id)

!   The two lines below are added to flush buffers to disk
      status = NF90_SYNC(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, time_id)
# endif

   end subroutine write_time_nc

   subroutine write_phi2_nc(nout, phi2)

# ifdef NETCDF
      use netcdf, only: nf90_put_var
# endif

      implicit none

      integer, intent(in) :: nout
      real, intent(in) :: phi2

# ifdef NETCDF
      integer :: status

      status = nf90_put_var(ncid, phi2_id, phi2, start=(/nout/))
      if (status /= nf90_noerr) call netcdf_error(status, ncid, phi2_id)
# endif

   end subroutine write_phi2_nc

   subroutine write_phi_nc(nout, phi)

      use convert, only: c2r
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: nakx, naky
# ifdef NETCDF
      use netcdf, only: nf90_put_var, nf90_sync
# endif

      implicit none

      integer, intent(in) :: nout
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi

# ifdef NETCDF
      integer :: status
      integer, dimension(6) :: start, count
      real, dimension(:, :, :, :, :), allocatable :: phi_ri

      start = 1
      start(6) = nout
      count(1) = 2
      count(2) = naky
      count(3) = nakx
      count(4) = 2 * nzgrid + 1
      count(5) = ntubes
      count(6) = 1

      allocate (phi_ri(2, naky, nakx, 2 * nzgrid + 1, ntubes))
      call c2r(phi, phi_ri)
      status = nf90_put_var(ncid, phi_vs_t_id, phi_ri, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, phi_vs_t_id)

!   Buffers to disk
      status = NF90_SYNC(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, phi_vs_t_id)

      deallocate (phi_ri)
# endif

   end subroutine write_phi_nc

   subroutine write_omega_nc(nout, omega)

# ifdef NETCDF
      use netcdf, only: nf90_put_var
# endif

      implicit none

      integer, intent(in) :: nout
      complex, dimension(:, :), intent(in) :: omega

# ifdef NETCDF
      integer :: status

      status = nf90_put_var(ncid, omega_id, [real(omega), aimag(omega)], start=[1, 1, 1, nout])
      if (status /= nf90_noerr) call netcdf_error(status, ncid, omega_id)
# endif

   end subroutine write_omega_nc

   subroutine write_radial_fluxes_nc(nout, pflux, vflux, qflux)

      use kt_grids, only: nakx
      use species, only: nspec
# ifdef NETCDF
      use netcdf, only: nf90_put_var
# endif

      implicit none

      integer, intent(in) :: nout
      real, dimension(:, :), intent(in) :: pflux, vflux, qflux

# ifdef NETCDF
      integer :: status
      integer, dimension(3) :: start, count

      start = 1
      start(3) = nout
      count(1) = nakx
      count(2) = nspec
      count(3) = 1

      status = nf90_put_var(ncid, pflux_x_id, pflux, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, pflux_x_id)
      status = nf90_put_var(ncid, vflux_x_id, vflux, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, vflux_x_id)
      status = nf90_put_var(ncid, qflux_x_id, qflux, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, qflux_x_id)
# endif

   end subroutine write_radial_fluxes_nc

   subroutine write_radial_moments_nc(nout, dens_x, upar_x, temp_x)

      use kt_grids, only: nakx
      use species, only: nspec
# ifdef NETCDF
      use netcdf, only: nf90_put_var
# endif

      implicit none

      integer, intent(in) :: nout
      real, dimension(:, :), intent(in) :: dens_x, upar_x, temp_x

# ifdef NETCDF
      integer :: status
      integer, dimension(3) :: start, count

      start = 1
      start(3) = nout
      count(1) = nakx
      count(2) = nspec
      count(3) = 1

      status = nf90_put_var(ncid, dens_x_id, dens_x, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, dens_x_id)
      status = nf90_put_var(ncid, upar_x_id, upar_x, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, upar_x_id)
      status = nf90_put_var(ncid, temp_x_id, temp_x, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, temp_x_id)
# endif

   end subroutine write_radial_moments_nc

   subroutine write_kspectra_nc(nout, phi2_vs_kxky)

      use kt_grids, only: nakx, naky
# ifdef NETCDF
      use netcdf, only: nf90_put_var, nf90_sync
# endif

      implicit none

      integer, intent(in) :: nout
      real, dimension(:, :), intent(in) :: phi2_vs_kxky

# ifdef NETCDF
      integer :: status
      integer, dimension(3) :: start, count

      start = 1
      start(3) = nout
      count(1) = naky
      count(2) = nakx
      count(3) = 1

      status = nf90_put_var(ncid, phi2_vs_kxky_id, phi2_vs_kxky, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, phi2_vs_kxky_id)

!   Buffers to disk
      status = NF90_SYNC(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, phi2_vs_kxky_id)
# endif

   end subroutine write_kspectra_nc
!
!
!
   subroutine write_fluxes_kxkyz_nc(nout, pflx_kxkyz, vflx_kxkyz, qflx_kxkyz)

      use kt_grids, only: nakx, naky
      use zgrid, only: nztot, ntubes
      use species, only: nspec
# ifdef NETCDF
      use netcdf, only: nf90_put_var, nf90_sync
# endif

      implicit none

      integer, intent(in) :: nout
      real, dimension(:, :, :, :, :), intent(in) :: pflx_kxkyz, vflx_kxkyz, qflx_kxkyz

# ifdef NETCDF
      integer :: status
      integer, dimension(6) :: start, count

      start = 1
      start(6) = nout
      count(1) = naky
      count(2) = nakx
      count(3) = nztot
      count(4) = ntubes
      count(5) = nspec
      count(6) = 1
!
      status = nf90_put_var(ncid, pflx_kxkyz_id, pflx_kxkyz, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, pflx_kxkyz_id)
!
      status = nf90_put_var(ncid, vflx_kxkyz_id, vflx_kxkyz, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, vflx_kxkyz_id)
!
      status = nf90_put_var(ncid, qflx_kxkyz_id, qflx_kxkyz, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, qflx_kxkyz_id)
!
!   Buffers to disk
!
      status = NF90_SYNC(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, pflx_kxkyz_id)
!
      status = NF90_SYNC(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, vflx_kxkyz_id)
!
      status = NF90_SYNC(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, qflx_kxkyz_id)
# endif

   end subroutine write_fluxes_kxkyz_nc
!
   subroutine write_moments_nc(nout, density, upar, temperature, spitzer2)

      use convert, only: c2r
      use zgrid, only: nztot, ntubes
      use kt_grids, only: nakx, naky
      use species, only: nspec
# ifdef NETCDF
      use netcdf, only: nf90_put_var, nf90_sync
# endif

      implicit none

      integer, intent(in) :: nout
      complex, dimension(:, :, :, :, :), intent(in) :: density, upar, temperature, spitzer2

# ifdef NETCDF
      integer :: status
      integer, dimension(7) :: start, count
      real, dimension(:, :, :, :, :, :), allocatable :: mom_ri

      start = 1
      start(7) = nout
      count(1) = 2
      count(2) = naky
      count(3) = nakx
      count(4) = nztot
      count(5) = ntubes
      count(6) = nspec
      count(7) = 1

      allocate (mom_ri(2, naky, nakx, nztot, ntubes, nspec))

      call c2r(density, mom_ri)
      status = nf90_put_var(ncid, density_id, mom_ri, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, density_id)

!   Buffers to disk
      status = NF90_SYNC(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, density_id)

      call c2r(upar, mom_ri)
      status = nf90_put_var(ncid, upar_id, mom_ri, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, upar_id)

!   Buffers to disk
      status = NF90_SYNC(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, upar_id)

      call c2r(temperature, mom_ri)
      status = nf90_put_var(ncid, temperature_id, mom_ri, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, temperature_id)

!   Buffers to disk
      status = NF90_SYNC(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, temperature_id)

      ! AVB: added: (move this to a separate diagnostic in the future)
      call c2r(spitzer2, mom_ri)
      status = nf90_put_var(ncid, spitzer2_id, mom_ri, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, spitzer2_id)

!   Buffers to disk
      status = NF90_SYNC(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, spitzer2_id)

      deallocate (mom_ri)

# endif

   end subroutine write_moments_nc

   subroutine write_gvmus_nc(nout, g)

      use vpamu_grids, only: nvpa, nmu
      use species, only: nspec
# ifdef NETCDF
      use netcdf, only: nf90_put_var, nf90_sync
# endif

      implicit none

      integer, intent(in) :: nout
      real, dimension(:, :, :), intent(in) :: g

# ifdef NETCDF
      integer :: status
      integer, dimension(4) :: start, count

      start(1) = 1
      start(2:3) = 1
      start(4) = nout
      count(1) = nvpa
      count(2) = nmu
      count(3) = nspec
      count(4) = 1

      status = nf90_put_var(ncid, gvmus_id, g, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, gvmus_id)

!   Buffers to disk
      status = NF90_SYNC(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, gvmus_id)

# endif

   end subroutine write_gvmus_nc

   subroutine write_gzvs_nc(nout, g)

      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: nvpa
      use species, only: nspec
# ifdef NETCDF
      use netcdf, only: nf90_put_var, nf90_sync
# endif

      implicit none

      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: g

# ifdef NETCDF
      integer :: status
      integer, dimension(5) :: start, count

      start = 1
      start(5) = nout
      count(1) = ntubes
      count(2) = 2 * nzgrid + 1
      count(3) = nvpa
      count(4) = nspec
      count(5) = 1

      status = nf90_put_var(ncid, gzvs_id, g, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, gzvs_id)

!   Buffers to disk
      status = NF90_SYNC(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, gzvs_id)

# endif

   end subroutine write_gzvs_nc

   !> Write [[species:spec]] to output netCDF file
   subroutine nc_species(file_id)
#ifdef NETCDF
      use species, only: spec, nspec
      use neasyf, only: neasyf_write
#endif
      implicit none
      !> NetCDF ID of the file to write to
      integer, intent(in) :: file_id
#ifdef NETCDF
      integer :: is
      ! FIXME: FLAG - ignoring cross-species collisions for now
      real, dimension(nspec) :: vnew
      do is = 1, nspec
         vnew(is) = spec(is)%vnew(is)
      end do

      ! Additional brackets around `(spec%z)` etc to workaround gfortran bug
      call neasyf_write(file_id, "charge", (spec%z), dim_names=["species"], &
                        long_name="Charge", units="e")
      call neasyf_write(file_id, "mass", (spec%mass), dim_names=["species"], &
                        long_name="Atomic mass", units="AMU")
      call neasyf_write(file_id, "dens", (spec%dens), dim_names=["species"], &
                        long_name="Normalised density", units="nref")
      call neasyf_write(file_id, "temp", (spec%temp), dim_names=["species"], &
                        long_name="Normalised temperature", units="Tref")
      call neasyf_write(file_id, "tprim", (spec%tprim), dim_names=["species"], &
                        long_name="Normalised temperature gradient scale length -1/rho dT/drho", units="1/aref")
      call neasyf_write(file_id, "fprim", (spec%fprim), dim_names=["species"], &
                        long_name="Normalised density gradient scale length -1/rho dn/drho", units="1/aref")
      call neasyf_write(file_id, "vnew", vnew, dim_names=["species"], &
                        long_name="Collisionality", units="vtref/aref")
      call neasyf_write(file_id, "type_of_species", (spec%type), dim_names=["species"], &
                        long_name="Species type: 1=ion, 2=electron, 3=slowing down, 4=trace")
#endif
   end subroutine nc_species

   subroutine nc_geo

      use stella_geometry, only: bmag, gradpar, gbdrift, gbdrift0, &
                                 cvdrift, cvdrift0, gds2, gds21, gds22, grho, jacob, &
                                 drhodpsi, djacdrho, b_dot_grad_z
      use stella_geometry, only: geo_surf
      use zgrid, only: nzgrid
      use physics_parameters, only: beta
      use dist_fn_arrays, only: kperp2
      use kt_grids, only: naky, nakx, nalpha, jtwist
# ifdef NETCDF
      use netcdf, only: nf90_put_var

      implicit none

      integer :: status
      integer, dimension(2) :: start, count
      integer, dimension(4) :: start2, count2

      start = 1
      count(1) = nalpha
      count(2) = 2 * nzgrid + 1

      start2 = 1
      count2(1) = naky
      count2(2) = nakx
      count2(3) = nalpha
      count2(4) = 2 * nzgrid + 1

      status = nf90_put_var(ncid, bmag_id, bmag, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, bmag_id)
      status = nf90_put_var(ncid, b_dot_grad_z_id, b_dot_grad_z, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, b_dot_grad_z_id)
      status = nf90_put_var(ncid, gradpar_id, gradpar)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, gradpar_id)
      status = nf90_put_var(ncid, gbdrift_id, gbdrift, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, gbdrift_id)
      status = nf90_put_var(ncid, gbdrift0_id, gbdrift0, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, gbdrift0_id)
      status = nf90_put_var(ncid, cvdrift_id, cvdrift, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, cvdrift_id)
      status = nf90_put_var(ncid, cvdrift0_id, cvdrift0, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, cvdrift0_id)
      status = nf90_put_var(ncid, kperp2_id, kperp2, start=start2, count=count2)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, kperp2_id)
      status = nf90_put_var(ncid, gds2_id, gds2, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, gds2_id)
      status = nf90_put_var(ncid, gds21_id, gds21, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, gds21_id)
      status = nf90_put_var(ncid, gds22_id, gds22, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, gds22_id)
      status = nf90_put_var(ncid, grho_id, grho, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, grho_id)
      status = nf90_put_var(ncid, jacob_id, jacob, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, jacob_id)
      status = nf90_put_var(ncid, djacdrho_id, djacdrho, start=start, count=count)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, djacdrho_id)

      status = nf90_put_var(ncid, beta_id, beta)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, beta_id)
      status = nf90_put_var(ncid, q_id, geo_surf%qinp)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, q_id)
      status = nf90_put_var(ncid, shat_id, geo_surf%shat)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, shat_id)
      status = nf90_put_var(ncid, d2qdr2_id, geo_surf%d2qdr2)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, d2qdr2_id)
      status = nf90_put_var(ncid, drhodpsi_id, drhodpsi)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, drhodpsi_id)
      status = nf90_put_var(ncid, d2psidr2_id, geo_surf%d2psidr2)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, d2psidr2_id)
      status = nf90_put_var(ncid, jtwist_id, jtwist)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, jtwist_id)
# endif
   end subroutine nc_geo

   !> Get the index of the time dimension in the netCDF file that corresponds to
   !> a time no larger than `tstart`
   subroutine get_nout(tstart, nout)

      use netcdf, only: nf90_inquire_dimension, nf90_get_var

      implicit none

      !> Simulation time to find
      real, intent(in) :: tstart
      !> Index of time dimension
      integer, intent(out) :: nout
      real, dimension(:), allocatable :: times
      integer :: i, length, status

      nout = 1

      status = nf90_inquire_dimension(ncid, time_dim, len=length)
      if (status /= nf90_noerr) call netcdf_error(status, ncid, dimid=time_dim)

      if (length > 0) then
         allocate (times(length))

         status = nf90_get_var(ncid, time_id, times)
         if (status /= nf90_noerr) call netcdf_error(status, ncid, dimid=time_dim)
         i = length
         do while (times(i) > tstart .and. i > 0)
            i = i - 1
         end do

         nout = i + 1

         deallocate (times)
      end if

   end subroutine get_nout

   subroutine sync_nc

# ifdef NETCDF
      use netcdf, only: nf90_sync

      implicit none
      integer :: status

      status = nf90_sync(ncid)
      if (status /= nf90_noerr) call netcdf_error(status, ncid)

# endif
   end subroutine sync_nc

end module stella_io
