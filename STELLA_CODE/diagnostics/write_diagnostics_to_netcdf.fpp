# include "define.inc"
!###############################################################################
!######################### WRITE OUTPUT TO NETCDF FILE #########################
!###############################################################################
! 
! Various diagnostics are written to the NetCDF file, depending on which 
! diagnostics flags have been turned on in the input file.
! 
! The variable <nout> keeps track of the time step pointer in the NetCDF file.
! 
! Some stella information is written to the NetCDF file:
!   - GitHub commit number
!   - GitHub commit date
!   - Simulation date
!   - Simulation time
!   - NetCDF version
! 
! The complete input file is written to the NetCDF file.
!     TODO - Perhaps we should write the full input file, including default variables.
! 
! The following variables, related to the dimensions, are written to the NetCDF file:
!   - Dimensions (kx, ky, z, mu, vpa, s, t)
!   - Other dimensions like (tube, alpha, ri)
!   - Dimension variables like (nkx, nky, ntubes, nzed, nspecies, nvpa, nmu)
!   - Other variables like (phase_shift_angle, theta0)
!   - Number of processors (nproc)
!   - The number of mesh points, i.e., nmesh = (2*nzgrid + 1) * ntubes * nvpa * nmu * nakx * naky * nspec
! 
! The species characteristics are written to the NetCDF file:
!   - charge, mass, dens, temp, vnew, tprim, fprim, type
! 
! The following geometric quantities are written to the NetCDF file:
!   - beta, q, shat, drhodpsi, jtwist, d2psidr2, d2qrd2,
!   - B_times_gradB_dot_grady, bmag, B_times_gradB_dot_gradx, B_times_kappa_dot_grady, 
!   - B_times_kappa_dot_gradx, gds2, gds21, gradx_dot_gradx, 
!   - grho, jacob, djacdrho, b_dot_grad_z
!   - gradpar, kperp2
! 
! The fields are written to the NetCDF file:
!   - |phi|^2(t)                    -->      phi2
!   - |apar|^2(t)                   -->      apar2
!   - |bpar|^2(t)                   -->      bpar2
!   - |phi|^2(t, kx, ky)            -->      phi2_vs_kxky
!   - |apar|^2(t, kx, ky)           -->      apar2_vs_kxky
!   - |bpar|^2(t, kx, ky)           -->      bpar2_vs_kxky
!   - phi(t, kx, ky, z, ri)         -->      phi_vs_t
!   - apar(t, kx, ky, z, ri)        -->      apar_vs_t
!   - bpar(t, kx, ky, z, ri)        -->      bpar_vs_t
! 
! The complex frequency (Omega = <omega>+i<gamma>) is written to the NetCDF file:
!   - Omega(t, kx, ky, ri)          -->      omega
! 
! The fluxes are written to the NetCDF file:
!   - pflux(t, s)                   -->      pflux_vs_s
!   - qflux(t, s)                   -->      qflux_vs_s
!   - vflux(t, s)                   -->      vflux_vs_s
!   - pflux(t, kx, ky, s)           -->      pflux_vs_kxkys
!   - qflux(t, kx, ky, s)           -->      qflux_vs_kxkys
!   - vflux(t, kx, ky, s)           -->      vflux_vs_kxkys
!   - pflux(t, kx, ky, z, s)        -->      pflux_vs_kxkyzs
!   - qflux(t, kx, ky, z, s)        -->      qflux_vs_kxkyzs
!   - vflux(t, kx, ky, z, s)        -->      vflux_vs_kxkyzs
!   - pflux(t, kx, s)               -->      pflux_x
!   - qflux(t, kx, s)               -->      qflux_x
!   - vflux(t, kx, s)               -->      vflux_x
! 
! The moments are written to the NetCDF file:
!   - dens(t, kx, ky, z, s, ri)     -->      density
!   - temp(t, kx, ky, z, s, ri)     -->      temperature
!   - upar(t, kx, ky, z, s, ri)     -->      upar
!   - dens(t, kx, s, ri)            -->      dens_x
!   - temp(t, kx, s, ri)            -->      upar_x
!   - upar(t, kx, s, ri)            -->      temp_x
! 
! The distribution functions are written to the NetCDF file:
!   - |g|^2(t, mu, vpa, s)             -->   g2_vs_vpamus
!   - |g|^2(t, z, vpa, s)              -->   g2_vs_zvpas
!   - |g|^2(t, z, mu, s)               -->   g2_vs_zmus
!   - |g|^2(t, kx, ky, z, s)           -->   g2_vs_zkykxs
!   - |g|^2(t, z, vpa, mu, s)          -->   g2_vs_zvpamus
!   - |g_nozonal|^2(t, mu, vpa, s)     -->   g2nozonal_vs_vpamus
!   - |g_nozonal|^2(t, z, vpa, s)      -->   g2nozonal_vs_zvpas
!   - |g_nozonal|^2(t, z, mu, s)       -->   g2nozonal_vs_zmus
!   - |g_nozonal|^2(t, kx, ky, z, s)   -->   g2nozonal_vs_zkykxs
!   - |g_nozonal|^2(t, z, vpa, mu, s)  -->   g2nozonal_vs_zvpamus
!   - |h|^2(t, mu, vpa, s)             -->   h2_vs_vpamus
!   - |h|^2(t, z, vpa, s)              -->   h2_vs_zvpas
!   - |h|^2(t, z, mu, s)               -->   h2_vs_zmus
!   - |h|^2(t, kx, ky, z, s)           -->   h2_vs_zkykxs
!   - |h|^2(t, z, vpa, mu, s)          -->   h2_vs_zvpamus
!   - |h_nozonal|^2(t, mu, vpa, s)     -->   h2nozonal_vs_vpamus
!   - |h_nozonal|^2(t, z, vpa, s)      -->   h2nozonal_vs_zvpas
!   - |h_nozonal|^2(t, z, mu, s)       -->   h2nozonal_vs_zmus
!   - |h_nozonal|^2(t, kx, ky, z, s)   -->   h2nozonal_vs_zkykxs
!   - |h_nozonal|^2(t, z, vpa, mu, s)  -->   h2nozonal_vs_zvpamus
!   - |f|^2(t, mu, vpa, s)             -->   f2_vs_vpamus
!   - |f|^2(t, z, vpa, s)              -->   f2_vs_zvpas
!   - |f|^2(t, z, mu, s)               -->   f2_vs_zmus
!   - |f|^2(t, kx, ky, z, s)           -->   f2_vs_zkykxs
!   - |f|^2(t, z, vpa, mu, s)          -->   f2_vs_zvpamus
!   - |f_nozonal|^2(t, mu, vpa, s)     -->   f2nozonal_vs_vpamus
!   - |f_nozonal|^2(t, z, vpa, s)      -->   f2nozonal_vs_zvpas
!   - |f_nozonal|^2(t, z, mu, s)       -->   f2nozonal_vs_zmus
!   - |f_nozonal|^2(t, kx, ky, z, s)   -->   f2nozonal_vs_zkykxs
!   - |f_nozonal|^2(t, z, vpa, mu, s)  -->   f2nozonal_vs_zvpamus
! 
! TODO-GA: Print input parameters to netcdf file
! 
!###############################################################################
module write_diagnostics_to_netcdf

   ! Load the NetCDF variables which are used in every routine
#ifdef NETCDF
   use netcdf, only: nf90_noerr
   use netcdf_utils, only: netcdf_error, kind_nf
   use neasyf, only: neasyf_write, neasyf_error
#endif

   implicit none

   private

   ! Make the routines available to other modules
   public :: init_write_diagnostics_to_netcdf
   public :: finish_write_diagnostics_to_netcdf
   
   ! Write time traces
   public :: write_time_nc
   public :: write_phi2_nc
   public :: write_apar2_nc
   public :: write_bpar2_nc
   public :: write_fluxes_vs_time_nc
   
   ! Write full (kx,ky,z,t) data
   public :: write_phi_nc
   public :: write_apar_nc
   public :: write_bpar_nc
   public :: write_g2_vs_vpamus_nc
   public :: write_g2_vs_zvpas_nc
   public :: write_g2_vs_zmus_nc
   public :: write_g2_vs_zkykxs_nc
   public :: write_g2_vs_zvpamus_nc
   public :: write_g2nozonal_vs_vpamus_nc
   public :: write_g2nozonal_vs_zvpas_nc
   public :: write_g2nozonal_vs_zmus_nc
   public :: write_g2nozonal_vs_zvpamus_nc
   public :: write_h2_vs_vpamus_nc
   public :: write_h2_vs_zvpas_nc
   public :: write_h2_vs_zmus_nc
   public :: write_h2_vs_zkykxs_nc
   public :: write_h2_vs_zvpamus_nc
   public :: write_h2nozonal_vs_vpamus_nc
   public :: write_h2nozonal_vs_zvpas_nc
   public :: write_h2nozonal_vs_zmus_nc
   public :: write_h2nozonal_vs_zvpamus_nc
   public :: write_f2_vs_vpamus_nc
   public :: write_f2_vs_zvpas_nc
   public :: write_f2_vs_zmus_nc
   public :: write_f2_vs_zkykxs_nc
   public :: write_f2_vs_zvpamus_nc
   public :: write_f2nozonal_vs_vpamus_nc
   public :: write_f2nozonal_vs_zvpas_nc
   public :: write_f2nozonal_vs_zmus_nc
   public :: write_f2nozonal_vs_zvpamus_nc
   public :: write_kspectra_nc
   public :: write_omega_nc
   public :: write_moments_nc
   public :: write_radial_fluxes_nc
   public :: write_radial_moments_nc
   public :: write_fluxes_kxkyzs_nc
   public :: write_fluxes_kxkys_nc
   public :: get_nout
   public :: sync_nc

   ! Variables for the netcdf files
#ifdef NETCDF
   integer(kind_nf) :: ncid
   integer(kind_nf) :: char10_dim
   integer :: code_id
#endif

   ! Converts a complex array to a real array with an extra dimension
#ifdef NETCDF
   interface netcdf_write_complex
      module procedure write_complex_rank2
      module procedure write_complex_rank4
      module procedure write_complex_rank5
   end interface netcdf_write_complex
#endif

   ! Define a zero
   real, parameter :: zero = epsilon(0.0)

contains

!###############################################################################
!############################ INITIALISE NETCDF FILE ###########################
!###############################################################################

   !****************************************************************************
   !                         Initialise the netcdf file                         
   !****************************************************************************
   subroutine init_write_diagnostics_to_netcdf(restart, git_commit, git_date)

      implicit none
   
      ! Arguments
      logical, intent(in) :: restart
      character(len=40), intent(in) :: git_commit
      character(len=10), intent(in) :: git_date

      !-------------------------------------------------------------------------
      
      ! Only continue if the netcdf library is available
#ifdef NETCDF
      call init_write_diagnostics_to_netcdf_file(restart, git_commit, git_date)
#else
      write (*,*) 'Warning: The NetCDF library is not available. Diagnostics can not be written.'
      write (*,*) 'Diagnostics can not be written to a NetCDF file.'
#endif
   
   end subroutine init_write_diagnostics_to_netcdf
 
#ifdef NETCDF
   !****************************************************************************
   !                         Initialise the netcdf file                         
   !****************************************************************************
   subroutine init_write_diagnostics_to_netcdf_file(restart, git_commit, git_date)

      ! Parallelisation
      use mp, only: proc0
      
      ! Read and write files
      use file_utils, only: run_name
      use neasyf, only: neasyf_open
      use neasyf, only: neasyf_metadata
      
      ! GitHub version
      use git_version, only: get_git_version

      implicit none
      
      ! If the simulation has been restarted, we continue writing on the existing netcdf file
      logical, intent(in) :: restart

      ! Print git information to netcdf file
      character(len=40), intent(in) :: git_commit
      character(len=10), intent(in) :: git_date

      ! Local variables
      character(300) :: filename

      !-------------------------------------------------------------------------

      ! The netcdf file has the extension ".out.nc"
      filename = trim(trim(run_name)//'.out.nc')

      ! Only the first processor (proc0) opens the file
      if (proc0) then
      
         ! If the simulation is being restarted, continue writing on the existing file
         ! Otherwise open a new file from scratch.
         if (restart) then
            ncid = neasyf_open(trim(filename), "rw")
         else
            ncid = neasyf_open(trim(filename), "w")
         end if

         ! Add meta data
         call neasyf_metadata(ncid, title="stella simulation data", software_name="stella", &
            software_version=get_git_version(), auto_date=.true.)

         ! Write the vectors corresponding to the dimensions
         call write_grids(ncid)

         ! Write the code info
         call write_text_to_netcdf_file(git_commit, git_date)
         
         ! Write the input file
         call write_input_file_to_netcdf_file(ncid)

         ! Write constants to the netcdf file
         call nc_species(ncid)
         call nc_geo(ncid)
         
      end if

   end subroutine init_write_diagnostics_to_netcdf_file
#endif

#ifdef NETCDF
   !****************************************************************************
   !                         Initialise the dimensions                          
   !****************************************************************************
   subroutine write_grids(file_id)
   
      ! Parallelisation
      use mp, only: nproc

      ! Grids
      use grids_kxky, only: nakx, naky, nalpha, phase_shift_angle
      use grids_kxky, only: x_d, rho_d, akx, aky, theta0
      use grids_z, only: nzgrid, ntubes, zed
      use grids_velocity, only: nvpa, vpa, nmu, mu
      use grids_species, only: nspec
      
      ! Geometry
      use parameters_physics, only: rhostar
      use geometry, only: geo_surf, dxdpsi, q_as_x
      use neasyf, only: neasyf_dim
      
      implicit none

      ! NetCDF ID of the file
      integer, intent(in) :: file_id

      ! Local variables
      integer :: ix
      real :: nmesh  ! Total number of mesh points


      !-------------------------------------------------------------------------
      !                                DIMENSIONS                               
      !-------------------------------------------------------------------------

      ! Dimensions (kx,ky,z,mu,vpa,s,t) as well as (alpha, tube) and (real, imaginary)
      call neasyf_dim(file_id, "kx", values=akx, long_name="Wavenumber in direction of grad alpha", units="1/rho_ref")
      call neasyf_dim(file_id, "ky", values=aky, long_name="Wavenumber perpendicular to flux surface", units="1/rho_ref")
      call neasyf_dim(file_id, "zed", values=zed)
      call neasyf_dim(file_id, "mu", values=mu)
      call neasyf_dim(file_id, "vpa", values=vpa)
      call neasyf_dim(file_id, "tube", dim_size=ntubes)
      call neasyf_dim(file_id, "alpha", dim_size=nalpha)
      call neasyf_dim(file_id, "species", dim_size=nspec)
      call neasyf_dim(file_id, "t", unlimited=.true., long_name="Time", units="a_ref/v_ref")
      call neasyf_dim(file_id, "ri", dim_size=2, long_name="Complex components", units="(real, imaginary)")

      ! Dimensions for various string variables
      call neasyf_dim(file_id, "char10", dim_size=10, dimid=char10_dim)
      call neasyf_dim(file_id, "char200", dim_size=200)
      call neasyf_dim(file_id, "nlines", unlimited=.true., long_name="Input file line number")

      !-------------------------------------------------------------------------
      !                                VARIABLES                                
      !-------------------------------------------------------------------------

      ! Integers related to the dimensions
      call neasyf_write(file_id, "nkx", nakx, long_name="Number of kx points")
      call neasyf_write(file_id, "nky", naky, long_name="Number of ky points")
      call neasyf_write(file_id, "ntubes", ntubes, long_name="Number of tubes")
      call neasyf_write(file_id, "nzed_tot", 2 * nzgrid + 1, long_name="Number of zed points")
      call neasyf_write(file_id, "nspecies", nspec, long_name="Number of species")
      call neasyf_write(file_id, "nvpa_tot", nvpa, long_name="Number of vpa points")
      call neasyf_write(file_id, "nmu", nmu, long_name="Number of mu points")
      
      ! Other variables related to the dimensions
      call neasyf_write(file_id, "phase_shift_angle", phase_shift_angle)
      call neasyf_write(file_id, "theta0", theta0, dim_names=["ky", "kx"], long_name="Theta_0")

      ! Number of processors
      call neasyf_write(file_id, "nproc", nproc, long_name="Number of processors")

      ! Radial variation
      call write_radial_grid

      ! Number of mesh points
      nmesh = (2 * nzgrid + 1) * ntubes * nvpa * nmu * nakx * naky * nspec
      call neasyf_write(file_id, "nmesh", nmesh, long_name="Number of mesh points")
      
   contains
   
      !-------------------------------------------------------------------------
      !                             RADIAL VARIATION                            
      !-------------------------------------------------------------------------
   
      subroutine write_radial_grid
      
         ! Flags
         use parameters_physics, only: radial_variation
         
         implicit none
         
         ! Local variables
         real, dimension(:, :), allocatable :: rg
         
         !-------------------------------------------------------------------------
         
         ! Only continue for radial variation runs
         if (.not. radial_variation) return
         
         ! Create the dimension of the radial grid
         call neasyf_dim(file_id, "radgridvar", dim_size=3, long_name="x, q/psi, rho")
         
         ! Add the values of the radial grid
         if (radial_variation) then
            allocate (rg(3, nakx))
            if (q_as_x) then
               do ix = 1, nakx
                  rg(1, ix) = x_d(ix)
                  rg(2, ix) = rhostar * x_d(ix) / dxdpsi + geo_surf%qinp
                  rg(3, ix) = rho_d(ix) + geo_surf%rhoc
               end do
            else
               do ix = 1, nakx
                  rg(1, ix) = x_d(ix)
                  rg(2, ix) = rhostar * x_d(ix) / dxdpsi
                  rg(3, ix) = rho_d(ix) + geo_surf%rhoc
               end do
            end if
            call neasyf_write(file_id, "rad_grid", rg, dim_names=[character(10)::"radgridvar", "kx"])
            deallocate (rg)
         end if
      
      end subroutine write_radial_grid

   end subroutine write_grids
#endif

#ifdef NETCDF
   !****************************************************************************
   !                  Write some text data to the NetCDF file
   !****************************************************************************
   subroutine write_text_to_netcdf_file(git_commit, git_date)

      ! NetCDF variables
      use netcdf, only: nf90_char
      use netcdf, only: nf90_def_var
      use netcdf, only: nf90_inq_varid
      use netcdf, only: nf90_put_att
      use netcdf, only: nf90_put_var
      use netcdf, only: nf90_inq_libvers

      implicit none

      ! Print git information to netcdf file
      character(len=40), intent(in) :: git_commit
      character(len=10), intent(in) :: git_date

      ! Local variables
      character(5) :: ci
      character(20) :: datestamp, timestamp, timezone
      integer :: status

      !-------------------------------------------------------------------------

      ! Write some useful general information such as the website, date and time into the NetCDF file
      datestamp(:) = ' '; timestamp(:) = ' '; timezone(:) = ' '
      call date_and_time(datestamp, timestamp, timezone)

      ! Create a viarble to hold the code info
      status = nf90_inq_varid(ncid, 'code_info', code_id)
      if (status /= nf90_noerr) then
         status = nf90_def_var(ncid, 'code_info', nf90_char, char10_dim, code_id)
         if (status /= nf90_noerr) call netcdf_error(status, var='code_info')
      end if
      status = nf90_put_att(ncid, code_id, 'long_name', 'stella')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att='long_name')

      ! Add strings with data as attributes
      ci = 'c1'; status = nf90_put_att(ncid, code_id, trim(ci), 'stella git commit: '//trim(git_commit))
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)
      ci = 'c2'; status = nf90_put_att(ncid, code_id, trim(ci), 'stella git date: '//trim(git_date))
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)
      ci = 'c3'; status = nf90_put_att(ncid, code_id, trim(ci), 'Simulation date: '//trim(datestamp))
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)
      ci = 'c4'; status = nf90_put_att(ncid, code_id, trim(ci), 'Simulation time: '//trim(timestamp)//' '//trim(timezone))
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)
      ci = 'c5'; status = nf90_put_att(ncid, code_id, trim(ci), 'netCDF version: '//trim(nf90_inq_libvers()))
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)
      ci = 'c6'; status = nf90_put_att(ncid, code_id, trim(ci), ' ')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

      ! Write a little explanation about the units, print each line to the netcdf file
      ci = 'c7'; status = nf90_put_att(ncid, code_id, trim(ci), 'Units are determined with respect to reference temperature (T_ref),')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)
      ci = 'c8'; status = nf90_put_att(ncid, code_id, trim(ci), 'reference charge (q_ref), reference mass (mass_ref),')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)
      ci = 'c9'; status = nf90_put_att(ncid, code_id, trim(ci), 'reference field (B_ref), and reference length (a_ref)')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)
      ci = 'c10'; status = nf90_put_att(ncid, code_id, trim(ci), 'from which one may construct rho_ref and vt_ref/a,')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)
      ci = 'c11'; status = nf90_put_att(ncid, code_id, trim(ci), 'which are the basic units of perpendicular length and time.')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)
      ci = 'c12'; status = nf90_put_att(ncid, code_id, trim(ci), 'Macroscopic lengths are normalized to the minor radius.')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)
      ci = 'c13'; status = nf90_put_att(ncid, code_id, trim(ci), 'The difference between rho (normalized minor radius) and rho (gyroradius)')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)
      ci = 'c14'; status = nf90_put_att(ncid, code_id, trim(ci), 'should be clear from the context in which they appear below.')
      if (status /= nf90_noerr) call netcdf_error(status, ncid, code_id, att=ci)

   end subroutine write_text_to_netcdf_file
#endif

#ifdef NETCDF
   !****************************************************************************
   !                   Save the input file to the NetCDF file                   
   !****************************************************************************
   subroutine write_input_file_to_netcdf_file(file_id)

      ! Read and write files
      use file_utils, only: num_input_lines
      use file_utils, only: get_input_unit
      
      ! NetCDF variables
      use netcdf, only: nf90_inq_dimid
      use netcdf, only: nf90_inquire_dimension
      use netcdf, only: NF90_NOERR, NF90_EBADDIM

      implicit none

      ! NetCDF ID of the file to write to
      integer, intent(in) :: file_id
      
      ! Local variables
      integer, parameter :: line_length = 200
      character(line_length), dimension(:), allocatable ::  namelist_array
      integer :: n, unit, status, dim_id, previous_nlines
      
      !-------------------------------------------------------------------------

      ! Dont attempt to write zero-sized text files
      if (num_input_lines <= 0) return

      ! If the existing input file in the output file was longer than
      ! the current one, blank out the whole thing so that we are not
      ! left with "extra" bits at the end
      status = nf90_inq_dimid(file_id, "nlines", dim_id)
      if (status == NF90_NOERR) then
      
         ! Get the current dimension of the input file
         status = nf90_inquire_dimension(file_id, dim_id, len=previous_nlines)
         call neasyf_error(status, ncid=file_id, dim="nlines", dimid=dim_id)

         ! Blank out the input file
         if (previous_nlines > num_input_lines) then
            allocate (namelist_array(previous_nlines))
            call neasyf_write(file_id, "input_file", namelist_array, long_name="Input file", dim_names=["char200", "nlines "])
            deallocate (namelist_array)
         end if
      else
         call neasyf_error(status, ncid=file_id, dim="nlines", dimid=dim_id)
      end if

      ! Convert the input file text into an array
      allocate (namelist_array(num_input_lines))
      call get_input_unit(unit)
      rewind (unit=unit)
      do n = 1, num_input_lines
         read (unit=unit, fmt="(a)") namelist_array(n)
      end do

      ! Write the input file to the NetCDF file
      call neasyf_write(file_id, "input_file", namelist_array, long_name="Input file", dim_names=["char200", "nlines "])

   end subroutine write_input_file_to_netcdf_file
#endif

!###############################################################################
!############################# WRITE CONSTANT DATA #############################
!###############################################################################
! Variables that do not change are written once at the beginning. 
!###############################################################################

#ifdef NETCDF
   !****************************************************************************
   !                                  SPECIES                                   
   !****************************************************************************
   subroutine nc_species(file_id)

      ! Grids
      use grids_species, only: spec, nspec

      implicit none
      
      ! NetCDF ID of the file to write to
      integer, intent(in) :: file_id

      ! Local variables
      integer :: is
      real, dimension(nspec) :: vnew
      
      !-------------------------------------------------------------------------

      ! TODO: FLAG - Ignoring cross-species collisions for now
      do is = 1, nspec
         vnew(is) = spec(is)%vnew(is)
      end do

      ! Save the characteristics of the species (charge, mass, temp, ...) to the NetCDF file.
      ! Note that we use additional brackets around e.g. '(spec%z)' to work around a gfortran bug
      call neasyf_write(file_id, "charge", (spec%z), dim_names=["species"], long_name="Charge", units="e")
      call neasyf_write(file_id, "mass", (spec%mass), dim_names=["species"], long_name="Atomic mass", units="AMU")
      call neasyf_write(file_id, "dens", (spec%dens), dim_names=["species"], long_name="Normalised density", units="nref")
      call neasyf_write(file_id, "temp", (spec%temp), dim_names=["species"], long_name="Normalised temperature", units="Tref")
      call neasyf_write(file_id, "vnew", vnew, dim_names=["species"], long_name="Collisionality", units="vtref/aref")
      call neasyf_write(file_id, "tprim", (spec%tprim), dim_names=["species"], long_name="Normalised temperature gradient scale length -1/rho dT/drho", units="1/aref")
      call neasyf_write(file_id, "fprim", (spec%fprim), dim_names=["species"], long_name="Normalised density gradient scale length -1/rho dn/drho", units="1/aref")
      call neasyf_write(file_id, "type_of_species", (spec%type), dim_names=["species"], long_name="Species type: 1=ion, 2=electron, 3=slowing down, 4=trace")

   end subroutine nc_species
#endif

#ifdef NETCDF
   !****************************************************************************
   !                                 GEOMETRY                                   
   !****************************************************************************
   subroutine nc_geo(file_id)
      
      ! Geometric quantities
      use geometry, only: bmag, gradpar, B_times_gradB_dot_grady, B_times_gradB_dot_gradx
      use geometry, only: B_times_kappa_dot_grady, B_times_kappa_dot_gradx
      use geometry, only: gds2, gds21, gradx_dot_gradx
      use geometry, only: grho, jacob
      use geometry, only: drhodpsi, djacdrho, b_dot_grad_z, geo_surf
      use parameters_physics, only: beta
      use arrays, only: kperp2
      use grids_kxky, only: jtwist

      implicit none

      ! NetCDF ID of the file to write to
      integer, intent(in) :: file_id

      ! Dimensions of the geometry arrays
      character(len=*), dimension(*), parameter :: flux_surface_dim = [character(5)::"alpha", "zed"]
      
      !-------------------------------------------------------------------------

      ! Variables
      call neasyf_write(file_id, "beta", beta, long_name="Reference beta", units="2.mu0.nref.Tref/B_a**2")
      call neasyf_write(file_id, "q", geo_surf%qinp, long_name="Local safety factor")
      call neasyf_write(file_id, "shat", geo_surf%shat, long_name="(rho/q) dq/drho")
      call neasyf_write(file_id, "drhodpsi", drhodpsi, long_name="drho/dPsi")
      call neasyf_write(file_id, "jtwist", jtwist, long_name="2*pi*shat*dky/dkx")
      call neasyf_write(file_id, "d2psidr2", geo_surf%d2psidr2)
      call neasyf_write(file_id, "d2qdr2", geo_surf%d2qdr2)

      ! Vectors along the field line
      call neasyf_write(file_id, "gradpar", gradpar, dim_names=["zed"], long_name="Parallel derivative multiplier")

      ! Vectors on the flux surface
      call neasyf_write(file_id, "B_times_gradB_dot_grady", B_times_gradB_dot_grady, dim_names=flux_surface_dim, long_name="Magnetic gradient drift")
      call neasyf_write(file_id, "bmag", bmag, dim_names=flux_surface_dim, long_name="Magnitude of magnetic field", units="B_0")
      call neasyf_write(file_id, "B_times_gradB_dot_gradx", B_times_gradB_dot_gradx, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "B_times_kappa_dot_grady", B_times_kappa_dot_grady, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "B_times_kappa_dot_gradx", B_times_kappa_dot_gradx, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "gds2", gds2, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "gds21", gds21, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "gradx_dot_gradx", gradx_dot_gradx, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "grho", grho, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "jacob", jacob, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "djacdrho", djacdrho, dim_names=flux_surface_dim)
      call neasyf_write(file_id, "b_dot_grad_z", b_dot_grad_z, dim_names=flux_surface_dim)

      ! Perpendicular wavevector depends on (kx,ky) and the flux surface (alpha,z)
      call neasyf_write(file_id, "kperp2", kperp2, dim_names=[character(len=5)::"ky", "kx", "alpha", "zed"])

   end subroutine nc_geo
#endif


!###############################################################################
!######################## WRITE DATA AT EVERY TIME STEP ########################
!###############################################################################
! For each variable, create a routine to write to the netcdf file
! 
! Arguments
! ---------
!    nout : current timestep
! 
!###############################################################################

   !----------------------------------- time -----------------------------------
   subroutine write_time_nc(nout, time)
      implicit none
      integer, intent(in) :: nout
      real, intent(in) :: time
#ifdef NETCDF
      call neasyf_write(ncid, "t", time, dim_names=["t"], start=[nout])
#endif
   end subroutine write_time_nc


   !============================================================================
   !================================ POTENTIAL =================================
   !============================================================================

   !--------------------------------- phi2(t) ----------------------------------
   subroutine write_phi2_nc(nout, phi2)
      implicit none
      integer, intent(in) :: nout
      real, intent(in) :: phi2
#ifdef NETCDF
      call neasyf_write(ncid, "phi2", phi2, dim_names=["t"], start=[nout], &
         units="(T_ref/q rho_ref/a_ref)**2", &
         long_name="Amplitude of electrostatic potential")
#endif
   end subroutine write_phi2_nc

!--------------------------------- apar2(t) ----------------------------------
   subroutine write_apar2_nc(nout, apar2)
      implicit none
      integer, intent(in) :: nout
      real, intent(in) :: apar2
#ifdef NETCDF
      call neasyf_write(ncid, "apar2", apar2, dim_names=["t"], start=[nout], &
         units="(B_ref (rho_ref)**2 / a_ref)**2", &
         long_name="Amplitude of parallel vector potential apar")
#endif
   end subroutine write_apar2_nc

!--------------------------------- bpar2(t) ----------------------------------
   subroutine write_bpar2_nc(nout, bpar2)
      implicit none
      integer, intent(in) :: nout 
      real, intent(in) :: bpar2
#ifdef NETCDF
      call neasyf_write(ncid, "bpar2", bpar2, dim_names=["t"], start=[nout], &
         units="(B_ref rho_ref / a_ref)**2", &
         long_name="Amplitude of parallel magnetic field fluctuation bpar")
#endif
   end subroutine write_bpar2_nc

   !--------------- phi2(ky,kx,t), apar2(ky,kx,t), bpar2(ky,kx,t) --------------
   subroutine write_kspectra_nc(nout, phi2_vs_kxky, keyname, longname)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :), intent(in) :: phi2_vs_kxky
      character(len=*), intent(in) :: keyname, longname
#ifdef NETCDF
      call neasyf_write(ncid, keyname, phi2_vs_kxky, &
         dim_names=[character(len=2)::"ky", "kx", "t "], start=[1, 1, nout], &
         units=" ", &
         long_name=longname)
#endif
   end subroutine write_kspectra_nc

   !-------------------------- phi(ri,ky,kx,z,tube,t) --------------------------
   subroutine write_phi_nc(nout, phi)
      use grids_z, only: nzgrid
      implicit none
      integer, intent(in) :: nout
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi
#ifdef NETCDF
      call netcdf_write_complex(ncid, "phi_vs_t", phi, &
         dim_names=[character(len=4)::"ri", "ky", "kx", "zed", "tube", "t"], start=[1, 1, 1, 1, 1, nout], &
         units=" ", &
         long_name="Electrostatic potential")
#endif
   end subroutine write_phi_nc

   !-------------------------- apar(ri,ky,kx,z,tube,t) --------------------------
   subroutine write_apar_nc(nout, apar)
      use grids_z, only: nzgrid
      implicit none
      integer, intent(in) :: nout
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: apar
#ifdef NETCDF
      call netcdf_write_complex(ncid, "apar_vs_t", apar, &
         dim_names=[character(len=4)::"ri", "ky", "kx", "zed", "tube", "t"], start=[1, 1, 1, 1, 1, nout], &
         units=" ", &
         long_name="Electromagnetic parallel vector potential apar")
#endif
   end subroutine write_apar_nc

   !-------------------------- bpar(ri,ky,kx,z,tube,t) --------------------------
   subroutine write_bpar_nc(nout, bpar)
      use grids_z, only: nzgrid
      implicit none
      integer, intent(in) :: nout
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: bpar
#ifdef NETCDF
      call netcdf_write_complex(ncid, "bpar_vs_t", bpar, &
         dim_names=[character(len=4)::"ri", "ky", "kx", "zed", "tube", "t"], start=[1, 1, 1, 1, 1, nout], &
         units=" ", &
         long_name="Electromagnetic field bpar")
# endif
   end subroutine write_bpar_nc

   !---------------------------- omega(ri,ky,kx,t) -----------------------------
   subroutine write_omega_nc(nout, omega)
      implicit none
      integer, intent(in) :: nout
      complex, dimension(:, :), intent(in) :: omega
#ifdef NETCDF
      call netcdf_write_complex(ncid, "omega", omega, &
         dim_names=[character(len=2)::"ri", "ky", "kx", "t"], start=[1, 1, 1, nout], &
         units="a_ref/v_ref", &
         long_name="Complex frequency Omega = omega+i*gamma")
#endif
   end subroutine write_omega_nc


   !============================================================================
   !================================== FLUXES ==================================
   !============================================================================
   
   !--------------------------------- flux(s,t) --------------------------------
   subroutine write_fluxes_vs_time_nc(nout, pflux_vs_s, vflux_vs_s, qflux_vs_s)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:), intent(in) :: pflux_vs_s, vflux_vs_s, qflux_vs_s
#ifdef NETCDF
      call neasyf_write(ncid, "pflux_vs_s", pflux_vs_s, &
         dim_names=[character(len=7)::"species", "t"], start=[1, nout], &
         units="n_ref * v_ref * (rho_ref/a_ref)^2 (with v_ref = sqrt(2 T_ref/m_ref))", &
         long_name="Particle flux(s,t)")
      call neasyf_write(ncid, "vflux_vs_s", vflux_vs_s, &
         dim_names=[character(len=7)::"species", "t"], start=[1, nout], &
         units="m_ref * n_ref * (v_ref)^2 * (rho_ref/a_ref)^2", &
         long_name="Momentum flux(s,t)")
      call neasyf_write(ncid, "qflux_vs_s", qflux_vs_s, &
         dim_names=[character(len=7)::"species", "t"], start=[1, nout], &
         units="n_ref * T_ref * v_ref * (rho_ref/a_ref)^2", &
         long_name="Heat flux(s,t)")
#endif
   end subroutine write_fluxes_vs_time_nc

   !-------------------------- flux(ky,kx,s,t) --------------------------
   subroutine write_fluxes_kxkys_nc(nout, pflux_kxkys, vflux_kxkys, qflux_kxkys)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :), intent(in) :: pflux_kxkys, vflux_kxkys, qflux_kxkys
#ifdef NETCDF
      call neasyf_write(ncid, "pflux_vs_kxkys", pflux_kxkys, &
         dim_names=[character(len=7)::"ky", "kx", "species", "t"], start=[1, 1, 1, nout], &
         units="n_ref * v_ref * (rho_ref/a_ref)^2 (with v_ref = sqrt(2 T_ref/m_ref))", &
         long_name="Particle flux(ky,kx,s,t)")
      call neasyf_write(ncid, "vflux_vs_kxkys", vflux_kxkys, &
         dim_names=[character(len=7)::"ky", "kx", "species", "t"], start=[1, 1, 1, nout], &
         units="m_ref * n_ref * (v_ref)^2 * (rho_ref/a_ref)^2", &
         long_name="Momentum flux(ky,kx,s,t)")
      call neasyf_write(ncid, "qflux_vs_kxkys", qflux_kxkys, &
         dim_names=[character(len=7)::"ky", "kx", "species", "t"], start=[1, 1, 1, nout], &
         units="n_ref * T_ref * v_ref * (rho_ref/a_ref)^2", &
         long_name="Heat flux(ky,kx,s,t)")
#endif
   end subroutine write_fluxes_kxkys_nc

   !-------------------------- flux(ky,kx,z,tube,s,t) --------------------------
   subroutine write_fluxes_kxkyzs_nc(nout, pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :, :), intent(in) :: pflux_kxkyzts, vflux_kxkyzts, qflux_kxkyzts
#ifdef NETCDF
      call neasyf_write(ncid, "pflux_vs_kxkyzs", pflux_kxkyzts, &
         dim_names=[character(len=7)::"ky", "kx", "zed", "tube", "species", "t"], start=[1, 1, 1, 1, 1, nout], &
         units="n_ref * v_ref * (rho_ref/a_ref)^2 (with v_ref = sqrt(2 T_ref/m_ref))", &
         long_name="Particle flux(ky,kx,z,tube,s,t)")
      call neasyf_write(ncid, "vflux_vs_kxkyzs", vflux_kxkyzts, &
         dim_names=[character(len=7)::"ky", "kx", "zed", "tube", "species", "t"], start=[1, 1, 1, 1, 1, nout], &
         units="m_ref*n_ref*(v_ref)^2*(rho_ref/a_ref)^2", &
         long_name="Momentum flux(ky,kx,z,tube,s,t)")
      call neasyf_write(ncid, "qflux_vs_kxkyzs", qflux_kxkyzts, &
         dim_names=[character(len=7)::"ky", "kx", "zed", "tube", "species", "t"], start=[1, 1, 1, 1, 1, nout], &
         units="n_ref * T_ref * v_ref * (rho_ref/a_ref)^2", &
         long_name="Heat flux(ky,kx,z,tube,s,t)")
#endif

   end subroutine write_fluxes_kxkyzs_nc

   !------------------------------ flux_x(kx,s,t) ------------------------------
   subroutine write_radial_fluxes_nc(nout, pflux, vflux, qflux)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :), intent(in) :: pflux, vflux, qflux
#ifdef NETCDF

      ! Define the dimensions and starting pointer
      character(*), dimension(*), parameter :: dims = [character(len=7)::"kx", "species", "t"]
      integer, dimension(3) :: start
      start = [1, 1, nout]

      ! Write the radial flux(kx,s,t)
      call neasyf_write(ncid, "pflux_x", pflux, dim_names=dims, start=start)
      call neasyf_write(ncid, "vflux_x", vflux, dim_names=dims, start=start)
      call neasyf_write(ncid, "qflux_x", qflux, dim_names=dims, start=start)
#endif
   end subroutine write_radial_fluxes_nc

   !============================================================================
   !================================= MOMENTS ==================================
   !============================================================================

   !----------------------- density(kx,ky,z,tube,s,t,ri) -----------------------
   subroutine write_moments_nc(nout, density, upar, temperature, spitzer2)
      implicit none
      integer, intent(in) :: nout
      complex, dimension(:, :, :, :, :), intent(in) :: density, upar, temperature, spitzer2
#ifdef NETCDF

      ! Define the dimensions and starting pointer
      character(*), dimension(*), parameter :: dims = [character(7)::"ri", "ky", "kx", "zed", "tube", "species", "t"]
      integer, dimension(7) :: start
      start = [1, 1, 1, 1, 1, 1, nout]

      ! Write the moments(kx,ky,z,tube,s,t,ri)
      call netcdf_write_complex(ncid, "density", density, dim_names=dims, start=start, long_name="Perturbed density")
      call netcdf_write_complex(ncid, "upar", upar, dim_names=dims, start=start, long_name="Perturbed upar")
      call netcdf_write_complex(ncid, "temperature", temperature, dim_names=dims, start=start, long_name="Perturbed temperature") 

      ! AVB: added: (move this to a separate diagnostic in the future)
      call netcdf_write_complex(ncid, "spitzer2", spitzer2, dim_names=dims, start=start, &
         long_name="Integral required for second Spitzer coefficient")
#endif
   end subroutine write_moments_nc

   !------------------------------ dens_x(kx,s,t) ------------------------------
   subroutine write_radial_moments_nc(nout, dens_x, upar_x, temp_x)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :), intent(in) :: dens_x, upar_x, temp_x
#ifdef NETCDF

      ! Define the dimensions and starting pointer
      character(*), dimension(*), parameter :: dims = [character(len=7)::"kx", "species", "t"]
      integer, dimension(3) :: start
      start = [1, 1, nout]

      ! Write the radial moments for density, parallel velocity, temperature
      call neasyf_write(ncid, "dens_x", dens_x, dim_names=dims, start=start)
      call neasyf_write(ncid, "upar_x", upar_x, dim_names=dims, start=start)
      call neasyf_write(ncid, "temp_x", temp_x, dim_names=dims, start=start)

#endif
   end subroutine write_radial_moments_nc

   !============================================================================
   !=============================== DISTRIBUTION ===============================
   !============================================================================

   !---------------------------- g2_vs_vpamus(vpa,mu,s,t) -----------------------------
   subroutine write_g2_vs_vpamus_nc(nout, g2_vs_vpamus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :), intent(in) :: g2_vs_vpamus
#ifdef NETCDF
      call neasyf_write(ncid, "g2_vs_vpamus", g2_vs_vpamus, &
         dim_names=[character(len=7)::"vpa", "mu", "species", "t"], start=[1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over real space")
#endif
   end subroutine write_g2_vs_vpamus_nc

   !--------------------------- g2_vs_zvpas(tube,z,vpa,s,t) ---------------------------
   subroutine write_g2_vs_zvpas_nc(nout, g2_vs_zvpas)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: g2_vs_zvpas
#ifdef NETCDF
      call neasyf_write(ncid, "g2_vs_zvpas", g2_vs_zvpas, &
         dim_names=[character(len=7)::"tube", "zed", "vpa", "species", "t"], start=[1, 1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over (kx, ky, mu)")
#endif
   end subroutine write_g2_vs_zvpas_nc

   !------------------------- g2_vs_zmus(tube,z,mu,s,t) -------------------------
   subroutine write_g2_vs_zmus_nc(nout, g2_vs_zmus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: g2_vs_zmus
#ifdef NETCDF
      call neasyf_write(ncid, "g2_vs_zmus", g2_vs_zmus, &
         dim_names=[character(len=7)::"tube", "zed", "mu", "species", "t"], start=[1, 1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over (kx, ky, vpa)")
#endif
   end subroutine write_g2_vs_zmus_nc

   !------------------------- g2_vs_zkykxs(tube,z,kx,ky,s,t) -------------------------
   subroutine write_g2_vs_zkykxs_nc(nout, g2_vs_zkykxs)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :, :), intent(in) :: g2_vs_zkykxs
#ifdef NETCDF
      call neasyf_write(ncid, "g2_vs_zkykxs", g2_vs_zkykxs, &
         dim_names=[character(len=7)::"tube", "zed", "ky", "kx", "species", "t"], start=[1, 1, 1, 1, 1, nout])
#endif
   end subroutine write_g2_vs_zkykxs_nc

   !------------------------- g2_vs_zvpamus(tube,z,vpa,mu,s,t) -------------------------
   subroutine write_g2_vs_zvpamus_nc(nout, g2_vs_zvpamus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :, :), intent(in) :: g2_vs_zvpamus
#ifdef NETCDF
      call neasyf_write(ncid, "g2_vs_zvpamus", g2_vs_zvpamus, &
         dim_names=[character(len=7)::"zed", "tube", "vpa", "mu", "species", "t"], start=[1, 1, 1, 1, 1, nout])
#endif
   end subroutine write_g2_vs_zvpamus_nc

   !---------------------------- g2nozonal_vs_vpamus(vpa,mu,s,t) -----------------------------
   subroutine write_g2nozonal_vs_vpamus_nc(nout, g2nozonal_vs_vpamus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :), intent(in) :: g2nozonal_vs_vpamus
#ifdef NETCDF
      call neasyf_write(ncid, "g2nozonal_vs_vpamus", g2nozonal_vs_vpamus, &
         dim_names=[character(len=7)::"vpa", "mu", "species", "t"], start=[1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over real space")
#endif
   end subroutine write_g2nozonal_vs_vpamus_nc

   !--------------------------- g2nozonal_vs_zvpas(tube,z,vpa,s,t) ---------------------------
   subroutine write_g2nozonal_vs_zvpas_nc(nout, g2nozonal_vs_zvpas)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: g2nozonal_vs_zvpas
#ifdef NETCDF
      call neasyf_write(ncid, "g2nozonal_vs_zvpas", g2nozonal_vs_zvpas, &
         dim_names=[character(len=7)::"tube", "zed", "vpa", "species", "t"], start=[1, 1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over (kx, ky, mu)")
#endif
   end subroutine write_g2nozonal_vs_zvpas_nc

   !------------------------- g2nozonal_vs_zmus(tube,z,mu,s,t) -------------------------
   subroutine write_g2nozonal_vs_zmus_nc(nout, g2nozonal_vs_zmus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: g2nozonal_vs_zmus
#ifdef NETCDF
      call neasyf_write(ncid, "g2nozonal_vs_zmus", g2nozonal_vs_zmus, &
         dim_names=[character(len=7)::"tube", "zed", "mu", "species", "t"], start=[1, 1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over (kx, ky, vpa)")
#endif
   end subroutine write_g2nozonal_vs_zmus_nc 

   !------------------------- g2nozonal_vs_zvpamus(tube,z,vpa,mu,s,t) -------------------------
   subroutine write_g2nozonal_vs_zvpamus_nc(nout, g2nozonal_vs_zvpamus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :, :), intent(in) :: g2nozonal_vs_zvpamus
#ifdef NETCDF
      call neasyf_write(ncid, "g2nozonal_vs_zvpamus", g2nozonal_vs_zvpamus, &
         dim_names=[character(len=7)::"zed", "tube", "vpa", "mu", "species", "t"], start=[1, 1, 1, 1, 1, nout])
#endif
   end subroutine write_g2nozonal_vs_zvpamus_nc

   !---------------------------- h2_vs_vpamus(vpa,mu,s,t) -----------------------------
   subroutine write_h2_vs_vpamus_nc(nout, h2_vs_vpamus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :), intent(in) :: h2_vs_vpamus
#ifdef NETCDF
      call neasyf_write(ncid, "h2_vs_vpamus", h2_vs_vpamus, &
         dim_names=[character(len=7)::"vpa", "mu", "species", "t"], start=[1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over real space")
#endif
   end subroutine write_h2_vs_vpamus_nc

   !--------------------------- h2_vs_zvpas(tube,z,vpa,s,t) ---------------------------
   subroutine write_h2_vs_zvpas_nc(nout, h2_vs_zvpas)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: h2_vs_zvpas
#ifdef NETCDF
      call neasyf_write(ncid, "h2_vs_zvpas", h2_vs_zvpas, &
         dim_names=[character(len=7)::"tube", "zed", "vpa", "species", "t"], start=[1, 1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over (kx, ky, mu)")
#endif
   end subroutine write_h2_vs_zvpas_nc

   !------------------------- h2_vs_zmus(tube,z,mu,s,t) -------------------------
   subroutine write_h2_vs_zmus_nc(nout, h2_vs_zmus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: h2_vs_zmus
#ifdef NETCDF
      call neasyf_write(ncid, "h2_vs_zmus", h2_vs_zmus, &
         dim_names=[character(len=7)::"tube", "zed", "mu", "species", "t"], start=[1, 1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over (kx, ky, vpa)")
#endif
   end subroutine write_h2_vs_zmus_nc

   !------------------------- h2_vs_zkykxs(tube,z,kx,ky,s,t) -------------------------
   subroutine write_h2_vs_zkykxs_nc(nout, h2_vs_zkykxs)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :, :), intent(in) :: h2_vs_zkykxs
#ifdef NETCDF
      call neasyf_write(ncid, "h2_vs_zkykxs", h2_vs_zkykxs, &
         dim_names=[character(len=7)::"tube", "zed", "ky", "kx", "species", "t"], start=[1, 1, 1, 1, 1, nout])
#endif
   end subroutine write_h2_vs_zkykxs_nc

   !------------------------- h2_vs_zvpamus(tube,z,vpa,mu,s,t) -------------------------
   subroutine write_h2_vs_zvpamus_nc(nout, h2_vs_zvpamus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :, :), intent(in) :: h2_vs_zvpamus
#ifdef NETCDF
      call neasyf_write(ncid, "h2_vs_zvpamus", h2_vs_zvpamus, &
         dim_names=[character(len=7)::"zed", "tube", "vpa", "mu", "species", "t"], start=[1, 1, 1, 1, 1, nout])
#endif
   end subroutine write_h2_vs_zvpamus_nc

   !---------------------------- h2nozonal_vs_vpamus(vpa,mu,s,t) -----------------------------
   subroutine write_h2nozonal_vs_vpamus_nc(nout, h2nozonal_vs_vpamus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :), intent(in) :: h2nozonal_vs_vpamus
#ifdef NETCDF
      call neasyf_write(ncid, "h2nozonal_vs_vpamus", h2nozonal_vs_vpamus, &
         dim_names=[character(len=7)::"vpa", "mu", "species", "t"], start=[1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over real space")
#endif
   end subroutine write_h2nozonal_vs_vpamus_nc

   !--------------------------- h2nozonal_vs_zvpas(tube,z,vpa,s,t) ---------------------------
   subroutine write_h2nozonal_vs_zvpas_nc(nout, h2nozonal_vs_zvpas)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: h2nozonal_vs_zvpas
#ifdef NETCDF
      call neasyf_write(ncid, "h2nozonal_vs_zvpas", h2nozonal_vs_zvpas, &
         dim_names=[character(len=7)::"tube", "zed", "vpa", "species", "t"], start=[1, 1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over (kx, ky, mu)")
#endif
   end subroutine write_h2nozonal_vs_zvpas_nc

   !------------------------- h2nozonal_vs_zmus(tube,z,mu,s,t) -------------------------
   subroutine write_h2nozonal_vs_zmus_nc(nout, h2nozonal_vs_zmus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: h2nozonal_vs_zmus
#ifdef NETCDF
      call neasyf_write(ncid, "h2nozonal_vs_zmus", h2nozonal_vs_zmus, &
         dim_names=[character(len=7)::"tube", "zed", "mu", "species", "t"], start=[1, 1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over (kx, ky, vpa)")
#endif
   end subroutine write_h2nozonal_vs_zmus_nc 

   !------------------------- h2nozonal_vs_zvpamus(tube,z,vpa,mu,s,t) -------------------------
   subroutine write_h2nozonal_vs_zvpamus_nc(nout, h2nozonal_vs_zvpamus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :, :), intent(in) :: h2nozonal_vs_zvpamus
#ifdef NETCDF
      call neasyf_write(ncid, "h2nozonal_vs_zvpamus", h2nozonal_vs_zvpamus, &
         dim_names=[character(len=7)::"zed", "tube", "vpa", "mu", "species", "t"], start=[1, 1, 1, 1, 1, nout])
#endif
   end subroutine write_h2nozonal_vs_zvpamus_nc

   !---------------------------- f2_vs_vpamus(vpa,mu,s,t) -----------------------------
   subroutine write_f2_vs_vpamus_nc(nout, f2_vs_vpamus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :), intent(in) :: f2_vs_vpamus
#ifdef NETCDF
      call neasyf_write(ncid, "f2_vs_vpamus", f2_vs_vpamus, &
         dim_names=[character(len=7)::"vpa", "mu", "species", "t"], start=[1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over real space")
#endif
   end subroutine write_f2_vs_vpamus_nc

   !--------------------------- f2_vs_zvpas(tube,z,vpa,s,t) ---------------------------
   subroutine write_f2_vs_zvpas_nc(nout, f2_vs_zvpas)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: f2_vs_zvpas
#ifdef NETCDF
      call neasyf_write(ncid, "f2_vs_zvpas", f2_vs_zvpas, &
         dim_names=[character(len=7)::"tube", "zed", "vpa", "species", "t"], start=[1, 1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over (kx, ky, mu)")
#endif
   end subroutine write_f2_vs_zvpas_nc

   !------------------------- f2_vs_zmus(tube,z,mu,s,t) -------------------------
   subroutine write_f2_vs_zmus_nc(nout, f2_vs_zmus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: f2_vs_zmus
#ifdef NETCDF
      call neasyf_write(ncid, "f2_vs_zmus", f2_vs_zmus, &
         dim_names=[character(len=7)::"tube", "zed", "mu", "species", "t"], start=[1, 1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over (kx, ky, vpa)")
#endif
   end subroutine write_f2_vs_zmus_nc

   !------------------------- f2_vs_zkykxs(tube,z,kx,ky,s,t) -------------------------
   subroutine write_f2_vs_zkykxs_nc(nout, f2_vs_zkykxs)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :, :), intent(in) :: f2_vs_zkykxs
#ifdef NETCDF
      call neasyf_write(ncid, "f2_vs_zkykxs", f2_vs_zkykxs, &
         dim_names=[character(len=7)::"tube", "zed", "ky", "kx", "species", "t"], start=[1, 1, 1, 1, 1, nout])
#endif
   end subroutine write_f2_vs_zkykxs_nc

   !------------------------- f2_vs_zvpamus(tube,z,vpa,mu,s,t) -------------------------
   subroutine write_f2_vs_zvpamus_nc(nout, f2_vs_zvpamus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :, :), intent(in) :: f2_vs_zvpamus
#ifdef NETCDF
      call neasyf_write(ncid, "f2_vs_zvpamus", f2_vs_zvpamus, &
         dim_names=[character(len=7)::"zed", "tube", "vpa", "mu", "species", "t"], start=[1, 1, 1, 1, 1, nout])
#endif
   end subroutine write_f2_vs_zvpamus_nc

   !---------------------------- f2nozonal_vs_vpamus(vpa,mu,s,t) -----------------------------
   subroutine write_f2nozonal_vs_vpamus_nc(nout, f2nozonal_vs_vpamus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :), intent(in) :: f2nozonal_vs_vpamus
#ifdef NETCDF
      call neasyf_write(ncid, "f2nozonal_vs_vpamus", f2nozonal_vs_vpamus, &
         dim_names=[character(len=7)::"vpa", "mu", "species", "t"], start=[1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over real space")
#endif
   end subroutine write_f2nozonal_vs_vpamus_nc

   !--------------------------- f2nozonal_vs_zvpas(tube,z,vpa,s,t) ---------------------------
   subroutine write_f2nozonal_vs_zvpas_nc(nout, f2nozonal_vs_zvpas)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: f2nozonal_vs_zvpas
#ifdef NETCDF
      call neasyf_write(ncid, "f2nozonal_vs_zvpas", f2nozonal_vs_zvpas, &
         dim_names=[character(len=7)::"tube", "zed", "vpa", "species", "t"], start=[1, 1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over (kx, ky, mu)")
#endif
   end subroutine write_f2nozonal_vs_zvpas_nc

   !------------------------- f2nozonal_vs_zmus(tube,z,mu,s,t) -------------------------
   subroutine write_f2nozonal_vs_zmus_nc(nout, f2nozonal_vs_zmus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :), intent(in) :: f2nozonal_vs_zmus
#ifdef NETCDF
      call neasyf_write(ncid, "f2nozonal_vs_zmus", f2nozonal_vs_zmus, &
         dim_names=[character(len=7)::"tube", "zed", "mu", "species", "t"], start=[1, 1, 1, 1, nout], &
         long_name="Guiding center distribution function averaged over (kx, ky, vpa)")
#endif
   end subroutine write_f2nozonal_vs_zmus_nc 

   !------------------------- f2nozonal_vs_zvpamus(tube,z,vpa,mu,s,t) -------------------------
   subroutine write_f2nozonal_vs_zvpamus_nc(nout, f2nozonal_vs_zvpamus)
      implicit none
      integer, intent(in) :: nout
      real, dimension(:, :, :, :, :), intent(in) :: f2nozonal_vs_zvpamus
#ifdef NETCDF
      call neasyf_write(ncid, "f2nozonal_vs_zvpamus", f2nozonal_vs_zvpamus, &
         dim_names=[character(len=7)::"zed", "tube", "vpa", "mu", "species", "t"], start=[1, 1, 1, 1, 1, nout])
#endif
   end subroutine write_f2nozonal_vs_zvpamus_nc

!###############################################################################
!################################## ROUTINES ###################################
!###############################################################################

   !****************************************************************************
   !                                 Get <nout>                                 
   !****************************************************************************
   ! Get the index of the time dimension in the netCDF file that corresponds to
   ! a time no larger than 'tstart'. Important when restarting simulations.
   !****************************************************************************
   subroutine get_nout(tstart, nout)

      ! NetCDF variables
      use netcdf, only: nf90_inquire_dimension, nf90_inq_dimid
      use neasyf, only: neasyf_read, neasyf_error

      implicit none

      ! Arguments
      real, intent(in) :: tstart ! Simulation time to find
      integer, intent(out) :: nout ! Index of time dimension
      
      ! Local variables
      real, dimension(:), allocatable :: times
      integer :: i, length, time_dim

      !-------------------------------------------------------------------------

      ! Initialise the time counter
      nout = 1

      ! Read the time dimension in the NetCDF file
      call neasyf_error(nf90_inq_dimid(ncid, "t", time_dim), ncid)
      call neasyf_error(nf90_inquire_dimension(ncid, time_dim, len=length), ncid)

      ! If we have restarted the simulation, the time vector will have a non-zero length.
      if (length > 0) then
      
         ! Allocate a temporary array for the time vector
         allocate (times(length))
         
         ! Read the time dimension in the NetCDF file
         call neasyf_read(ncid, "t", times)

         ! Find the time index where time > tstart
         i = length
         do while (times(i) > tstart .and. i > 0)
            i = i - 1
         end do

         ! Set the time counter
         nout = i + 1

         ! Deallocate the temporary array
         deallocate (times)
         
      end if

   end subroutine get_nout

   !****************************************************************************
   !                         Flush NetCDF file to disk                          
   !****************************************************************************
   subroutine sync_nc
#ifdef NETCDF
      use netcdf, only: nf90_sync
      use neasyf, only: neasyf_error
      call neasyf_error(nf90_sync(ncid), ncid=ncid, message="Couldn't flush to disk")
#endif
   end subroutine sync_nc

   !****************************************************************************
   !                             WRITE COMPLEX DATA                             
   !****************************************************************************
   subroutine write_complex_rank2(parent_id, name, values, dim_names, units, long_name, start)
   
      use neasyf, only: neasyf_write
      use convert, only: c2r
      
      implicit none
      
      ! Arguments
      character(len=*), intent(in) :: name ! Name of the variable
      integer, intent(in) :: parent_id ! NetCDF ID of the parent group/file
      complex, dimension(:, :), intent(in) :: values ! Array to be written
      character(len=*), dimension(:), intent(in) :: dim_names ! Array of dimension names
      character(len=*), optional, intent(in) :: units ! Units of coordinate
      character(len=*), optional, intent(in) :: long_name ! Long descriptive name
      integer, dimension(:), optional, intent(in) :: start

      ! Local variables
      real, dimension(2, size(values, 1), size(values, 2)) :: real_values

      !-------------------------------------------------------------------------

      ! Convert from a complex array to a real array with two dimensions,
      ! seperating the complex values in (real, imaginary) values
      call c2r(values, real_values)
      
      ! Write the real array to the NetCDF file
      call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
      
   end subroutine write_complex_rank2

   !----------------------------------------------------------------------------
   subroutine write_complex_rank4(parent_id, name, values, dim_names, units, long_name, start)
   
      use neasyf, only: neasyf_write
      use convert, only: c2r
      
      implicit none
       
      ! Arguments
      character(len=*), intent(in) :: name ! Name of the variable
      integer, intent(in) :: parent_id ! NetCDF ID of the parent group/file
      complex, dimension(:, :, :, :), intent(in) :: values ! Array to be written
      character(len=*), dimension(:), intent(in) :: dim_names ! Array of dimension names
      character(len=*), optional, intent(in) :: units ! Units of coordinate
      character(len=*), optional, intent(in) :: long_name ! Long descriptive name
      integer, dimension(:), optional, intent(in) :: start

      ! Local variables
      real, dimension(2, size(values, 1), size(values, 2), size(values, 3), size(values, 4)) :: real_values

      !-------------------------------------------------------------------------

      ! Convert from a complex array to a real array with two dimensions,
      ! seperating the complex values in (real, imaginary) values
      call c2r(values, real_values)
      
      ! Write the real array to the NetCDF file
      call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
      
   end subroutine write_complex_rank4

   !----------------------------------------------------------------------------
   subroutine write_complex_rank5(parent_id, name, values, dim_names, units, long_name, start)
   
      use neasyf, only: neasyf_write
      use convert, only: c2r
      
      implicit none
       
      ! Arguments
      character(len=*), intent(in) :: name ! Name of the variable
      integer, intent(in) :: parent_id ! NetCDF ID of the parent group/file
      complex, dimension(:, :, :, :, :), intent(in) :: values ! Array to be written
      character(len=*), dimension(:), intent(in) :: dim_names ! Array of dimension names
      character(len=*), optional, intent(in) :: units ! Units of coordinate
      character(len=*), optional, intent(in) :: long_name ! Long descriptive name
      integer, dimension(:), optional, intent(in) :: start

      ! Local variables
      real, dimension(2, size(values, 1), size(values, 2), size(values, 3), size(values, 4), size(values, 5)) :: real_values

      !-------------------------------------------------------------------------

      ! Convert from a complex array to a real array with two dimensions,
      ! seperating the complex values in (real, imaginary) values
      call c2r(values, real_values)
      
      ! Write the real array to the NetCDF file
      call neasyf_write(parent_id, name, real_values, dim_names=dim_names, units=units, long_name=long_name, start=start)
      
   end subroutine write_complex_rank5
   
!###############################################################################
!############################## FINISH NETCDF FILE #############################
!###############################################################################

   !****************************************************************************
   !                           FINISH THE NETCDF FILE                           
   !****************************************************************************
   subroutine finish_write_diagnostics_to_netcdf
#ifdef NETCDF
      use mp, only: proc0
      use neasyf, only: neasyf_close

      if (proc0) then
         call neasyf_close(ncid)
      end if
#endif
   end subroutine finish_write_diagnostics_to_netcdf

end module write_diagnostics_to_netcdf
