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
   public :: write_apar2_nc
   public :: write_bpar2_nc
   public :: write_phi_nc
   public :: write_apar_nc
   public :: write_bpar_nc
   public :: write_gvmus_nc
   public :: write_gzvs_nc
   public :: write_kspectra_nc
   public :: write_kspectra_species_nc
   public :: write_omega_nc
   public :: write_moments_nc
   public :: write_radial_fluxes_nc
   public :: write_radial_moments_nc
   public :: write_fluxes_kxkyz_nc
   public :: write_fluxes_nc
   public :: write_stress_nc
   public :: get_nout
   public :: sync_nc

# ifdef NETCDF
   integer(kind_nf) :: ncid
   integer(kind_nf) :: char10_dim

   integer :: code_id

   !> Write a `complex` array to netcdf
   !>
   !> Converts the `complex` array to a `real` array with an extra dimension
   interface netcdf_write_complex
      module procedure write_complex_rank2, write_complex_rank4, write_complex_rank5
   end interface netcdf_write_complex
# endif

   real, parameter :: zero = epsilon(0.0)

contains

   !==============================================
   !============ INITIATE STELLA IO ==============
   !==============================================
   subroutine init_stella_io(restart)
# ifdef NETCDF
      use mp, only: proc0
      use file_utils, only: run_name
      use neasyf, only: neasyf_open, neasyf_metadata
      use git_version, only: get_git_version
# endif
      implicit none
      !> Is this run a restart?
      logical, intent(in) :: restart
# ifdef NETCDF
      character(300) :: filename

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
         call define_vars()

         call nc_species(ncid)
         call nc_geo(ncid)
         call save_input(ncid)
      end if
# endif

   end subroutine init_stella_io

   !> Ensure the netCDF file contains all the dimensions and grids,
   !> creating them if necessary
   subroutine write_grids(file_id)
# ifdef NETCDF
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
      call neasyf_dim(file_id, "ky", values=aky, long_name="Wavenumber perpendicular to flux surface", units="1/rho_r")
      call neasyf_dim(file_id, "kx", values=akx, long_name="Wavenumber in direction of grad alpha", units="1/rho_r")
      call neasyf_dim(file_id, "tube", dim_size=ntubes)
      call neasyf_dim(file_id, "zed", values=zed)
      call neasyf_dim(file_id, "alpha", dim_size=nalpha)
      call neasyf_dim(file_id, "vpa", values=vpa)
      call neasyf_dim(file_id, "mu", values=mu)
      call neasyf_dim(file_id, "species", dim_size=nspec)
      call neasyf_dim(file_id, "t", unlimited=.true., long_name="Time", units="L/vt")

      ! Dimensions for various string variables
      call neasyf_dim(file_id, "char10", dim_size=10, dimid=char10_dim)
      call neasyf_dim(file_id, "char200", dim_size=200)
      call neasyf_dim(file_id, "nlines", unlimited=.true., long_name="Input file line number")

      call neasyf_dim(file_id, "ri", dim_size=2, long_name="Complex components", units="(real, imaginary)")

      if (radial_variation) then
         call neasyf_dim(file_id, "radgridvar", dim_size=3, long_name="x, q/psi, rho")
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
      use neasyf, only: neasyf_close

      if (proc0) then
         call neasyf_close(ncid)
      end if
# endif
   end subroutine finish_stella_io

   !> Save the input file in the NetCDF file
   subroutine save_input(file_id)
#ifdef NETCDF
      use file_utils, only: num_input_lines, get_input_unit
      use neasyf, only: neasyf_write, neasyf_error
      use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension, NF90_NOERR, NF90_EBADDIM
#endif
      implicit none
      !> NetCDF ID of the file to write to
      integer, intent(in) :: file_id
#ifdef NETCDF
      integer, parameter :: line_length = 200
      character(line_length), dimension(:), allocatable ::  input_file_array
      integer :: n, unit, status, dim_id, previous_nlines

      ! Don't attempt to write zero-sized arrays
      if (num_input_lines <= 0) return

      ! If the existing input file in the output file was longer than
      ! the current one, blank out the whole thing so that we're not
      ! left with "extra" bits at the end
      status = nf90_inq_dimid(file_id, "nlines", dim_id)
      if (status == NF90_NOERR) then
         status = nf90_inquire_dimension(file_id, dim_id, len=previous_nlines)
         call neasyf_error(status, ncid=file_id, dim="nlines", dimid=dim_id)

         if (previous_nlines > num_input_lines) then
            allocate (input_file_array(previous_nlines))
            call neasyf_write(file_id, "input_file", input_file_array, &
                              long_name="Input file", dim_names=["char200", "nlines "])
            deallocate (input_file_array)
         end if
      else
         call neasyf_error(status, ncid=file_id, dim="nlines", dimid=dim_id)
      end if

      ! We need to convert the input file text into an array, one
      ! element per line
      allocate (input_file_array(num_input_lines))

      call get_input_unit(unit)
      rewind (unit=unit)
      do n = 1, num_input_lines
         read (unit=unit, fmt="(a)") input_file_array(n)
      end do

      call neasyf_write(file_id, "input_file", input_file_array, &
                        long_name="Input file", dim_names=["char200", "nlines "])
#endif
   end subroutine save_input

   subroutine define_vars()
# ifdef NETCDF
      use netcdf, only: nf90_char
      use netcdf, only: nf90_def_var, nf90_inq_varid, nf90_put_att, nf90_put_var
      use netcdf, only: nf90_inq_libvers
# endif
      implicit none
# ifdef NETCDF

      character(5) :: ci
      character(20) :: datestamp, timestamp, timezone

      integer :: status

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

# endif
   end subroutine define_vars

   !> Write simulation time to netCDF output file
   !>
   !> Also flushes to disk
   subroutine write_time_nc(nout, time)
# ifdef NETCDF
      use neasyf, only: neasyf_write, neasyf_error
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      !> Current simulation time
      real, intent(in) :: time

# ifdef NETCDF
      call neasyf_write(ncid, "t", time, dim_names=["t"], start=[nout])
# endif
   end subroutine write_time_nc

   !> Write amplitude of electrostatic potential to netCDF
   subroutine write_phi2_nc(nout, phi2)
# ifdef NETCDF
      use neasyf, only: neasyf_write
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      !> Amplitude of electrostatic potential
      real, intent(in) :: phi2

# ifdef NETCDF
      call neasyf_write(ncid, "phi2", phi2, dim_names=["t"], &
                        units="(T_ref/q rho_ref/L)**2", &
                        long_name="Amplitude of electrostatic potential", &
                        start=[nout])
# endif
   end subroutine write_phi2_nc
   

   subroutine write_stress_nc(nout, Ah_stress, phih_stress, bparh_stress, g_int)
# ifdef NETCDF
      use neasyf, only: neasyf_write
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      !> stresses:
      complex, dimension(:,:), intent(in) :: Ah_stress, phih_stress, bparh_stress
      complex, dimension(:,:), intent(in) :: g_int

# ifdef NETCDF
      call netcdf_write_complex(ncid, "phih_stress", phih_stress,  &
      	   			      		     [character(len=4)::"ri","kx", "zed", "t"], &
                                long_name="Reynolds stress as a function of kx and zed", &
                                start=[1, 1, 1, nout])
      call netcdf_write_complex(ncid, "Ah_stress", Ah_stress, &
                                [character(len=4)::"ri","kx", "zed", "t"], &
                                long_name="Maxwell stress as a function of kx and zed", &
						   	     start=[1, 1, 1, nout])
      call netcdf_write_complex(ncid, "bparh_stress", bparh_stress, &
                                [character(len=4)::"ri","kx", "zed", "t"], &
                                long_name="Bpar stress as a function of kx and zed", &
                                                             start=[1, 1, 1, nout])
      call netcdf_write_complex(ncid, "g_int", g_int, &
                                [character(len=4)::"ri","kx", "zed", "t"], &
		                long_name="integral of g as a function of kx and zed", &
                                          start=[1, 1, 1, nout])
# endif
   end subroutine write_stress_nc


   subroutine write_apar2_nc(nout, apar2)
# ifdef NETCDF
      use neasyf, only: neasyf_write
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      !> Amplitude of parallel vector potential
      real, intent(in) :: apar2

# ifdef NETCDF
      call neasyf_write(ncid, "apar2", apar2, dim_names=["t"], &
                        units="(B_ref (rho_ref)**2 / L)**2", &
                        long_name="Amplitude of parallel vector potential apar", &
                        start=[nout])
# endif
   end subroutine write_apar2_nc

   subroutine write_bpar2_nc(nout, bpar2)
# ifdef NETCDF
      use neasyf, only: neasyf_write
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      !> Amplitude of parallel vector potential
      real, intent(in) :: bpar2

# ifdef NETCDF
      call neasyf_write(ncid, "bpar2", bpar2, dim_names=["t"], &
                        units="(B_ref rho_ref / L)**2", &
                        long_name="Amplitude of parallel magnetic field fluctuation bpar", &
                        start=[nout])
# endif
   end subroutine write_bpar2_nc

   !> Write time trace of electrostatic potential to netCDF
   subroutine write_phi_nc(nout, phi)
      use zgrid, only: nzgrid
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      !> Electrostatic potential
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi

# ifdef NETCDF
      call netcdf_write_complex(ncid, "phi_vs_t", phi, &
                                [character(len=4)::"ri", "ky", "kx", "zed", "tube", "t"], &
                                long_name="Electrostatic potential", &
                                start=[1, 1, 1, 1, 1, nout])
# endif
   end subroutine write_phi_nc

   !> Write time trace of electromagnetic field A|| to netCDF
   subroutine write_apar_nc(nout, apar)
      use zgrid, only: nzgrid
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: apar

# ifdef NETCDF
      call netcdf_write_complex(ncid, "apar_vs_t", apar, &
                                [character(len=4)::"ri", "ky", "kx", "zed", "tube", "t"], &
                                long_name="Electromagnetic parallel vector potential apar", &
                                start=[1, 1, 1, 1, 1, nout])
# endif
   end subroutine write_apar_nc

   !> Write time trace of electromagnetic field B|| to netCDF
   subroutine write_bpar_nc(nout, bpar)
      use zgrid, only: nzgrid
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: bpar

# ifdef NETCDF
      call netcdf_write_complex(ncid, "bpar_vs_t", bpar, &
                                [character(len=4)::"ri", "ky", "kx", "zed", "tube", "t"], &
                                long_name="Electromagnetic field bpar", &
                                start=[1, 1, 1, 1, 1, nout])
# endif
   end subroutine write_bpar_nc

   !> Write the complex frequency to netCDF
   subroutine write_omega_nc(nout, omega, dphidt, dbpardt)
      implicit none

      integer, intent(in) :: nout
      complex, dimension(:, :), intent(in) :: omega
      complex, dimension(:, :), intent(in) :: dphidt
      complex, dimension(:, :), intent(in) :: dbpardt

# ifdef NETCDF
      call netcdf_write_complex(ncid, "omega", omega, &
                                dim_names=["ri", "ky", "kx", "t "], &
                                start=[1, 1, 1, nout], &
                                long_name="Complex frequency", units="aref/vtref")
      call netcdf_write_complex(ncid, "dphidt", dphidt, &
                                dim_names=["ri ", "kx ", "zed", "t  "], &
                                start=[1, 1, 1, nout], &
                                long_name="dphi dt", units="aref/vtref")
      call netcdf_write_complex(ncid, "dbpardt", dbpardt, &
                                dim_names=["ri ", "kx ", "zed", "t  "], &
				start=[1, 1, 1, nout], &
                                long_name="dbpar dt", units="aref/vtref")
# endif
   end subroutine write_omega_nc

   !> Write fluxes to netCDF
   subroutine write_radial_fluxes_nc(nout, pflux, vflux, qflux)
# ifdef NETCDF
      use neasyf, only: neasyf_write
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      !> Particle, velocity, heat flux
      real, dimension(:, :), intent(in) :: pflux, vflux, qflux
# ifdef NETCDF
      call neasyf_write(ncid, "pflux_x", pflux, &
                        dim_names=[character(len=7)::"kx", "species", "t"], &
                        start=[1, 1, nout])
      call neasyf_write(ncid, "vflux_x", vflux, &
                        dim_names=[character(len=7)::"kx", "species", "t"], &
                        start=[1, 1, nout])
      call neasyf_write(ncid, "qflux_x", qflux, &
                        dim_names=[character(len=7)::"kx", "species", "t"], &
                        start=[1, 1, nout])
# endif

   end subroutine write_radial_fluxes_nc

   !> Write radial moments to netCDF
   subroutine write_radial_moments_nc(nout, dens_x, upar_x, temp_x)
# ifdef NETCDF
      use neasyf, only: neasyf_write
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      !> Radial moments for density, parallel velocity, temperature
      real, dimension(:, :), intent(in) :: dens_x, upar_x, temp_x

# ifdef NETCDF
      call neasyf_write(ncid, "dens_x", dens_x, &
                        dim_names=[character(len=7)::"kx", "species", "t"], &
                        start=[1, 1, nout])
      call neasyf_write(ncid, "upar_x", upar_x, &
                        dim_names=[character(len=7)::"kx", "species", "t"], &
                        start=[1, 1, nout])
      call neasyf_write(ncid, "temp_x", temp_x, &
                        dim_names=[character(len=7)::"kx", "species", "t"], &
                        start=[1, 1, nout])
# endif
   end subroutine write_radial_moments_nc

   subroutine write_kspectra_nc(nout, field_vs_kxky, keyname, longname)
# ifdef NETCDF
      use neasyf, only: neasyf_write
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      real, dimension(:, :), intent(in) :: field_vs_kxky
      character(len=*), intent(in) :: keyname, longname
# ifdef NETCDF
      call neasyf_write(ncid, keyname, field_vs_kxky, &
                        dim_names=["ky", "kx", "t "], &
                        start=[1, 1, nout], &
                        long_name=longname)
# endif
   end subroutine write_kspectra_nc

   subroutine write_kspectra_species_nc(nout, field_vs_kxkys, keyname, longname)
# ifdef NETCDF
      use neasyf, only: neasyf_write
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      real, dimension(:, :, :), intent(in) :: field_vs_kxkys
      character(len=*), intent(in) :: keyname, longname
# ifdef NETCDF
      call neasyf_write(ncid, keyname, field_vs_kxkys, &
                        dim_names=[character(len=7)::"ky", "kx", "species", "t "], &
                        start=[1, 1, 1, nout], &
                        long_name=longname)
# endif
   end subroutine write_kspectra_species_nc

   subroutine write_fluxes_nc(nout, pflx, vflx, qflx, qflx_phi, qflx_apar, qflx_bpar)
# ifdef NETCDF
      use neasyf, only: neasyf_write
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      real, dimension(:), intent(in) :: pflx, vflx, qflx, qflx_phi, qflx_apar, qflx_bpar

# ifdef NETCDF
      call neasyf_write(ncid, "pflx", pflx, &
                        dim_names=[character(len=7)::"species", "t"], &
                        start=[1, nout], &
                        units="TBD", &
                        long_name="Particle flux")
      call neasyf_write(ncid, "vflx", vflx, &
                        dim_names=[character(len=7)::"species", "t"], &
                        start=[1, nout], &
                        units="TBD", &
                        long_name="Momentum flux")
      call neasyf_write(ncid, "qflx", qflx, &
                        dim_names=[character(len=7)::"species", "t"], &
                        start=[1, nout], &
                        units="TBD", &
                        long_name="Total heat flux")
      call neasyf_write(ncid, "qflx_phi", qflx_phi, &
                        dim_names=[character(len=7)::"species", "t"], &
			start=[1, nout], &
			units="TBD", &
			long_name="Electrostatic heat flux")
      call neasyf_write(ncid, "qflx_apar", qflx_apar, &
                        dim_names=[character(len=7)::"species", "t"], &
                        start=[1, nout], &
                        units="TBD", &
                        long_name="Apar contribution to heat flux")
      call neasyf_write(ncid, "qflx_bpar", qflx_bpar, &
                        dim_names=[character(len=7)::"species", "t"], &
                        start=[1, nout], &
                        units="TBD", &
                        long_name="Bpar contribution to heat flux")
# endif
   end subroutine write_fluxes_nc

   subroutine write_fluxes_kxkyz_nc(nout, pflx_kxkyz, vflx_kxkyz, qflx_kxkyz, qflx_kxkyz_phi, qflx_kxkyz_apar, qflx_kxkyz_bpar)
# ifdef NETCDF
      use neasyf, only: neasyf_write
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      real, dimension(:, :, :, :, :), intent(in) :: pflx_kxkyz, vflx_kxkyz, qflx_kxkyz
      real, dimension(:, :, :, :, :), intent(in) :: qflx_kxkyz_phi, qflx_kxkyz_apar, qflx_kxkyz_bpar
# ifdef NETCDF
      call neasyf_write(ncid, "pflx_kxky", pflx_kxkyz, &
                        dim_names=[character(len=7)::"ky", "kx", "zed", "tube", "species", "t"], &
                        start=[1, 1, 1, 1, 1, nout], &
                        long_name="Particle flux")
      call neasyf_write(ncid, "vflx_kxky", vflx_kxkyz, &
                        dim_names=[character(len=7)::"ky", "kx", "zed", "tube", "species", "t"], &
                        start=[1, 1, 1, 1, 1, nout], &
                        long_name="Momentum flux")
      call neasyf_write(ncid, "qflx_kxky", qflx_kxkyz, &
                        dim_names=[character(len=7)::"ky", "kx", "zed", "tube", "species", "t"], &
			start=[1, 1, 1, 1, 1, nout], &
			long_name="Total heat flux")
      call neasyf_write(ncid, "qflx_kxky_phi", qflx_kxkyz_phi, &
                        dim_names=[character(len=7)::"ky", "kx", "zed", "tube", "species", "t"], &
                        start=[1, 1, 1, 1, 1, nout], &
                        long_name="Electrostatic heat flux")
      call neasyf_write(ncid, "qflx_kxky_apar", qflx_kxkyz_apar, &
                        dim_names=[character(len=7)::"ky", "kx", "zed", "tube", "species", "t"], &
                        start=[1, 1, 1, 1, 1, nout], &
                        long_name="Apar contribution to heat flux")
      call neasyf_write(ncid, "qflx_kxky_bpar", qflx_kxkyz_bpar, &
                        dim_names=[character(len=7)::"ky", "kx", "zed", "tube", "species", "t"], &
                        start=[1, 1, 1, 1, 1, nout], &
                        long_name="Bpar contribution to heat flux")
# endif
   end subroutine write_fluxes_kxkyz_nc

   subroutine write_moments_nc(nout, density, upar, temperature, spitzer2)
      implicit none

      integer, intent(in) :: nout
      complex, dimension(:, :, :, :, :), intent(in) :: density, upar, temperature, spitzer2

# ifdef NETCDF
      character(*), dimension(*), parameter :: dims = [character(7)::"ri", "ky", "kx", "zed", "tube", "species", "t"]
      integer, dimension(7) :: start

      start = [1, 1, 1, 1, 1, 1, nout]

      call netcdf_write_complex(ncid, "density", density, &
                                dim_names=dims, start=start, &
                                long_name="Perturbed density")

      call netcdf_write_complex(ncid, "upar", upar, &
                                dim_names=dims, start=start, &
                                long_name="Perturbed upar")

      call netcdf_write_complex(ncid, "temperature", temperature, &
                                dim_names=dims, start=start, &
                                long_name="Perturbed temperature")

      ! AVB: added: (move this to a separate diagnostic in the future)
      call netcdf_write_complex(ncid, "spitzer2", spitzer2, &
                                dim_names=dims, start=start, &
                                long_name="Integral required for second Spitzer coefficient")
# endif
   end subroutine write_moments_nc

   !> Write guiding center distribution function averaged over real space
   subroutine write_gvmus_nc(nout, gvmus)
# ifdef NETCDF
      use neasyf, only: neasyf_write
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      !> Guiding centre distribution function
      real, dimension(:, :, :), intent(in) :: gvmus

# ifdef NETCDF
      call neasyf_write(ncid, "gvmus", gvmus, &
                        dim_names=[character(len=7)::"vpa", "mu", "species", "t"], &
                        start=[1, 1, 1, nout], &
                        long_name="Guiding center distribution function averaged over real space")
# endif
   end subroutine write_gvmus_nc

   !> Write guiding center distribution function averaged over (kx, ky, mu)
   subroutine write_gzvs_nc(nout, gzvs)
# ifdef NETCDF
      use neasyf, only: neasyf_write
# endif
      implicit none
      !> Current timestep
      integer, intent(in) :: nout
      !> Guiding centre distribution function
      real, dimension(:, :, :, :), intent(in) :: gzvs

# ifdef NETCDF
      call neasyf_write(ncid, "gzvs", gzvs, &
                        dim_names=[character(len=7)::"tube", "zed", "vpa", "species", "t"], &
                        start=[1, 1, 1, 1, nout], &
                        long_name="Guiding center distribution function averaged over (kx, ky, mu)")
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

   !> Write various geometric quantities to output netCDF file
   subroutine nc_geo(file_id)
# ifdef NETCDF
      use neasyf, only: neasyf_write
      use stella_geometry, only: bmag, gradpar, gbdrift, gbdrift0, &
                                 cvdrift, cvdrift0, gds2, gds21, gds22, grho, jacob, &
                                 drhodpsi, djacdrho, b_dot_grad_z
      use stella_geometry, only: geo_surf
      use physics_parameters, only: beta
      use dist_fn_arrays, only: kperp2
      use kt_grids, only: naky, nakx, jtwist
#endif
      implicit none
      !> NetCDF ID of the file to write to
      integer, intent(in) :: file_id

# ifdef NETCDF
      character(len=*), dimension(*), parameter :: flux_surface_dim = ["alpha", "zed  "]

      call neasyf_write(file_id, "bmag", bmag, dim_names=flux_surface_dim, &
                        long_name="Magnitude of magnetic field", units="B_0")
      call neasyf_write(file_id, "b_dot_grad_z", b_dot_grad_z, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "gradpar", gradpar, dim_names=["zed"], &
                        long_name="Parallel derivative multiplier")
      call neasyf_write(file_id, "gbdrift", gbdrift, dim_names=flux_surface_dim, &
                        long_name="Magnetic gradient drift")
      call neasyf_write(file_id, "gbdrift0", gbdrift0, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "cvdrift", cvdrift, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "cvdrift0", cvdrift0, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "kperp2", kperp2, dim_names=[character(len=5)::"ky", "kx", "alpha", "zed"])
      call neasyf_write(file_id, "gds2", gds2, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "gds21", gds21, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "gds22", gds22, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "grho", grho, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "jacob", jacob, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "djacdrho", djacdrho, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "beta", beta, &
                        long_name="Reference beta", units="2.mu0.nref.Tref/B_a**2")
      call neasyf_write(file_id, "q", geo_surf%qinp, &
                        long_name="Local safety factor")
      call neasyf_write(file_id, "shat", geo_surf%shat, &
                        long_name="(rho/q) dq/drho")
      call neasyf_write(file_id, "d2qdr2", geo_surf%d2qdr2)
      call neasyf_write(file_id, "drhodpsi", drhodpsi, &
                        long_name="drho/dPsi")
      call neasyf_write(file_id, "d2psidr2", geo_surf%d2psidr2)
      call neasyf_write(file_id, "jtwist", jtwist, &
                        long_name="2*pi*shat*dky/dkx")
# endif
   end subroutine nc_geo

   !> Get the index of the time dimension in the netCDF file that corresponds to
   !> a time no larger than `tstart`
   subroutine get_nout(tstart, nout)

      use netcdf, only: nf90_inquire_dimension, nf90_inq_dimid
      use neasyf, only: neasyf_read, neasyf_error

      implicit none

      !> Simulation time to find
      real, intent(in) :: tstart
      !> Index of time dimension
      integer, intent(out) :: nout
      real, dimension(:), allocatable :: times
      integer :: i, length, time_dim

      nout = 1

      call neasyf_error(nf90_inq_dimid(ncid, "t", time_dim), ncid)
      call neasyf_error(nf90_inquire_dimension(ncid, time_dim, len=length), ncid)

      if (length > 0) then
         allocate (times(length))
         call neasyf_read(ncid, "t", times)

         i = length
         do while (times(i) > tstart .and. i > 0)
            i = i - 1
         end do

         nout = i + 1

         deallocate (times)
      end if

   end subroutine get_nout

   !> Flush netCDF file to disk
   subroutine sync_nc
#ifdef NETCDF
      use netcdf, only: nf90_sync
      use neasyf, only: neasyf_error

      call neasyf_error(nf90_sync(ncid), ncid=ncid, message="Couldn't flush to disk")
#endif
   end subroutine sync_nc

   subroutine write_complex_rank2(parent_id, name, values, dim_names, units, long_name, start)
      use neasyf, only: neasyf_write
      use convert, only: c2r
      !> Name of the variable
      character(len=*), intent(in) :: name
      !> NetCDF ID of the parent group/file
      integer, intent(in) :: parent_id
      !> Array to be written
      complex, dimension(:, :), intent(in) :: values
      !> Array of dimension names
      character(len=*), dimension(:), intent(in) :: dim_names
      !> Units of coordinate
      character(len=*), optional, intent(in) :: units
      !> Long descriptive name
      character(len=*), optional, intent(in) :: long_name
      integer, dimension(:), optional, intent(in) :: start

      real, dimension(2, &
                      size(values, 1), &
                      size(values, 2) &
                      ) :: real_values

      call c2r(values, real_values)
      call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
   end subroutine write_complex_rank2

   subroutine write_complex_rank4(parent_id, name, values, dim_names, units, long_name, start)
      use neasyf, only: neasyf_write
      use convert, only: c2r
      !> Name of the variable
      character(len=*), intent(in) :: name
      !> NetCDF ID of the parent group/file
      integer, intent(in) :: parent_id
      !> Array to be written
      complex, dimension(:, :, :, :), intent(in) :: values
      !> Array of dimension names
      character(len=*), dimension(:), intent(in) :: dim_names
      !> Units of coordinate
      character(len=*), optional, intent(in) :: units
      !> Long descriptive name
      character(len=*), optional, intent(in) :: long_name
      integer, dimension(:), optional, intent(in) :: start

      real, dimension(2, &
                      size(values, 1), &
                      size(values, 2), &
                      size(values, 3), &
                      size(values, 4) &
                      ) :: real_values

      call c2r(values, real_values)
      call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
   end subroutine write_complex_rank4

   subroutine write_complex_rank5(parent_id, name, values, dim_names, units, long_name, start)
      use neasyf, only: neasyf_write
      use convert, only: c2r
      !> Name of the variable
      character(len=*), intent(in) :: name
      !> NetCDF ID of the parent group/file
      integer, intent(in) :: parent_id
      !> Array to be written
      complex, dimension(:, :, :, :, :), intent(in) :: values
      !> Array of dimension names
      character(len=*), dimension(:), intent(in) :: dim_names
      !> Units of coordinate
      character(len=*), optional, intent(in) :: units
      !> Long descriptive name
      character(len=*), optional, intent(in) :: long_name
      integer, dimension(:), optional, intent(in) :: start

      real, dimension(2, &
                      size(values, 1), &
                      size(values, 2), &
                      size(values, 3), &
                      size(values, 4), &
                      size(values, 5) &
                      ) :: real_values

      call c2r(values, real_values)
      call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
   end subroutine write_complex_rank5
end module stella_io