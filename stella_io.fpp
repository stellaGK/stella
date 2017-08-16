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
  public :: write_phi_nc
  public :: write_gvmus_nc
  public :: write_gzvs_nc

# ifdef NETCDF
  integer (kind_nf) :: ncid

  integer (kind_nf) :: naky_dim, nakx_dim, nttot_dim, nmu_dim, nvtot_dim, nspec_dim, ncoord_dim, ncoordt_dim
  integer (kind_nf) :: time_dim, char10_dim, char200_dim, ri_dim, nlines_dim, nheat_dim
  integer (kind_nf) :: nttotext_dim, time_big_dim

  integer, dimension (6) :: mom_t_dim
  integer, dimension (5) :: field_dim, final_mom_dim, heatk_dim
  integer, dimension (4) :: vmus_dim
  integer, dimension (4) :: zvs_dim
  integer, dimension (5) :: phi_corr_dim
  integer, dimension (4) :: omega_dim, fluxk_dim, final_field_dim, loop_mom_dim
  integer, dimension (4) :: phi_corr_2pi_dim
  integer, dimension (3) :: fluxx_dim
  integer, dimension (3) :: mode_dim, phase_dim, loop_phi_dim, heat_dim
  integer, dimension (2) :: kx_dim, ky_dim, om_dim, flux_dim, nin_dim, fmode_dim

  integer :: nakx_id, naky_id, nttot_id, akx_id, aky_id, theta_id, nspec_id
  integer :: nmu_id, nvtot_id, mu_id, vpa_id
  integer :: time_id, phi2_id, apar2_id, bpar2_id, theta0_id, nproc_id, nmesh_id
  integer :: phi2_by_mode_id, apar2_by_mode_id, bpar2_by_mode_id
  integer :: phtot_id, dmix_id, kperpnorm_id
  integer :: phi2_by_kx_id, apar2_by_kx_id, bpar2_by_kx_id
  integer :: phi2_by_ky_id, apar2_by_ky_id, bpar2_by_ky_id
  integer :: phi0_id, apar0_id, bpar0_id
  integer :: omega_id, omegaavg_id, phase_id
  integer :: es_heat_flux_id, es_mom_flux_id, es_part_flux_id, es_energy_exchange_id
  integer :: es_heat_par_id, es_heat_perp_id
  integer :: apar_heat_flux_id, apar_mom_flux_id, apar_part_flux_id
  integer :: apar_heat_par_id, apar_heat_perp_id
  integer :: bpar_heat_flux_id, bpar_mom_flux_id, bpar_part_flux_id
  integer :: bpar_heat_par_id, bpar_heat_perp_id
  integer :: hflux_tot_id, zflux_tot_id, vflux_tot_id
  integer :: es_heat_by_k_id, es_mom_by_k_id, es_part_by_k_id
  integer :: es_parmom_by_k_id, es_perpmom_by_k_id, es_mom0_by_k_id, es_mom1_by_k_id
  integer :: phi_corr_id, phi_corr_2pi_id
  integer :: apar_heat_by_k_id, apar_mom_by_k_id, apar_part_by_k_id
  integer :: apar_heat_by_x_id
  integer :: bpar_heat_by_k_id, bpar_mom_by_k_id, bpar_part_by_k_id
  integer :: phi_vs_t_id
  integer :: gvmus_id, gzvs_id
  integer :: apar_t_id, bpar_t_id
  integer :: ntot_t_id
  integer :: phi_norm_id, apar_norm_id, bpar_norm_id
  integer :: phi_id, apar_id, bpar_id, epar_id
  integer :: antot_id, antota_id, antotp_id
  integer :: ntot_id, density_id, upar_id, tpar_id, tperp_id
  integer :: qparflux_id, pperpj1_id, qpperpj1_id
  integer :: ntot2_id, ntot2_by_mode_id, ntot20_id, ntot20_by_mode_id
  integer :: tpar2_by_mode_id, tperp2_by_mode_id
  integer :: phi00_id, ntot00_id, density00_id, upar00_id, tpar00_id, tperp00_id
  integer :: input_id
  integer :: charge_id, mass_id, dens_id, temp_id, tprim_id, fprim_id
  integer :: vnewk_id, spec_type_id
  integer :: bmag_id, gradpar_id, gbdrift_id, gbdrift0_id
  integer :: cdrift_id, cdrift0_id
  integer :: cvdrift_id, cvdrift0_id, gds2_id, gds21_id, gds22_id
  integer :: grho_id, jacob_id, shat_id, eps_id, drhodpsi_id, q_id, surfarea_id
  integer :: beta_id
  integer :: code_id, datestamp_id, timestamp_id, timezone_id
  integer, dimension (5) :: mom_dim
  integer :: ntot0_id, density0_id, upar0_id, tpar0_id, tperp0_id
  logical :: write_apar_t, write_bpar_t ! Should the fields be written out every nwrite?

# endif
  real :: zero
  
!  include 'netcdf.inc'
  
contains

  subroutine init_stella_io (write_phi_vs_t, write_gvmus, write_gzvs)

    use mp, only: proc0
    use file_utils, only: run_name
# ifdef NETCDF
    use netcdf, only: nf90_clobber, nf90_create
    use netcdf_utils, only: get_netcdf_code_precision, netcdf_real
# endif

    implicit none

    logical, intent(in) :: write_phi_vs_t, write_gvmus, write_gzvs
# ifdef NETCDF
    character (300) :: filename
    integer :: status

    zero = epsilon(0.0)

    if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()
    status = nf90_noerr

    filename = trim(trim(run_name)//'.out.nc')
    ! only proc0 opens the file:
    if (proc0) then
       status = nf90_create (trim(filename), nf90_clobber, ncid) 
       if (status /= nf90_noerr) call netcdf_error (status, file=filename)

       call define_dims
       call define_vars (write_phi_vs_t, write_gvmus, write_gzvs)
       call nc_grids
       call nc_species
       call nc_geo
    end if
# endif

  end subroutine init_stella_io

  subroutine define_dims

    use file_utils, only: num_input_lines
    use kt_grids, only: naky, nakx
    use theta_grid, only: ntgrid
    use vpamu_grids, only: nvgrid, nmu
    use species, only: nspec
# ifdef NETCDF
    use netcdf, only: nf90_unlimited
    use netcdf, only: nf90_def_dim
# endif

# ifdef NETCDF
    integer :: status

    ! Associate the grid variables, e.g. ky, kx, with their size, e.g. naky, ntheta0 (= nakx),
    ! and a variable which is later used to store these sizes in the NetCDF file, e.g. naky_dim, nakx_dim
    status = nf90_def_dim (ncid, 'ky', naky, naky_dim)
    if (status /= nf90_noerr) call netcdf_error (status, dim='ky')
    status = nf90_def_dim (ncid, 'kx', nakx, nakx_dim)
    if (status /= nf90_noerr) call netcdf_error (status, dim='kx')
    status = nf90_def_dim (ncid, 'theta', 2*ntgrid+1, nttot_dim)
    if (status /= nf90_noerr) call netcdf_error (status, dim='theta')
    status = nf90_def_dim (ncid, 'vpa', 2*nvgrid+1, nvtot_dim)
    if (status /= nf90_noerr) call netcdf_error (status, dim='vpa')
    status = nf90_def_dim (ncid, 'mu', nmu, nmu_dim)
    if (status /= nf90_noerr) call netcdf_error (status, dim='mu')
    status = nf90_def_dim (ncid, 'species', nspec, nspec_dim)
    if (status /= nf90_noerr) call netcdf_error (status, dim='species')
    status = nf90_def_dim (ncid, 't', nf90_unlimited, time_dim)
    if (status /= nf90_noerr) call netcdf_error (status, dim='t')
    status = nf90_def_dim (ncid, 'char10', 10, char10_dim)
    if (status /= nf90_noerr) call netcdf_error (status, dim='char10')
    status = nf90_def_dim (ncid, 'char200', 200, char200_dim)
    if (status /= nf90_noerr) call netcdf_error (status, dim='char200')
    status = nf90_def_dim (ncid, 'nlines', num_input_lines, nlines_dim)
    if (status /= nf90_noerr) call netcdf_error (status, dim='nlines')
    status = nf90_def_dim (ncid, 'ri', 2, ri_dim)
    if (status /= nf90_noerr) call netcdf_error (status, dim='ri')
# endif
  end subroutine define_dims

  subroutine nc_grids

    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, nakx, theta0, akx, aky
    use species, only: nspec
    use vpamu_grids, only: nvgrid, nmu, vpa, mu
!    use nonlinear_terms, only: nonlin
# ifdef NETCDF
    use netcdf, only: nf90_put_var
    use constants, only: pi
    
    integer :: status
    real :: nmesh

    ! Store the size of the grid dimensions (as defined in def_dims), in the NetCDF file
    status = nf90_put_var (ncid, nttot_id, 2*ntgrid+1)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, nttot_id)
    status = nf90_put_var (ncid, naky_id, naky)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, naky_id)
    status = nf90_put_var (ncid, nakx_id, nakx)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, nakx_id)
    status = nf90_put_var (ncid, nspec_id, nspec)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, nspec_id)
    status = nf90_put_var (ncid, nmu_id, nmu)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, nmu_id)
    status = nf90_put_var (ncid, nvtot_id, 2*nvgrid+1)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, nvtot_id)

    status = nf90_put_var (ncid, akx_id, akx)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, akx_id)
    status = nf90_put_var (ncid, aky_id, aky)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, aky_id)
    status = nf90_put_var (ncid, theta_id, theta)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, theta_id)
    status = nf90_put_var (ncid, theta0_id, theta0)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, theta0_id)
    status = nf90_put_var (ncid, mu_id, mu)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, mu_id)
    status = nf90_put_var (ncid, vpa_id, vpa)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, vpa_id)

!    if (nonlin) then
!       nmesh = (2*ntgrid+1)*(2*nvgrid+1)*nmu*nx*ny*nspec
!    else
       nmesh = (2*ntgrid+1)*(2*nvgrid+1)*nmu*nakx*naky*nspec
!    end if

    status = nf90_put_var (ncid, nmesh_id, nmesh)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, nmesh_id)

# endif
  end subroutine nc_grids

  subroutine finish_stella_io
    use mp, only: proc0
# ifdef NETCDF
    use netcdf, only: nf90_close
    use netcdf_utils, only: netcdf_error

    integer :: status

    if (proc0) then
       call save_input
       status = nf90_close (ncid)
       if (status /= nf90_noerr) call netcdf_error (status)
    end if
# endif
  end subroutine finish_stella_io

  subroutine save_input
    !<doc> Save the input file in the NetCDF file </doc>
# ifdef NETCDF    
    use file_utils, only: num_input_lines, get_input_unit
    use netcdf, only: nf90_put_var

    character(200) line
    integer, dimension (2) :: nin_start, nin_count

    integer :: status, n, unit

    nin_start(1) = 1
    nin_start(2) = 1

    nin_count(2) = 1

    call get_input_unit (unit)
    rewind (unit=unit)
    do n = 1, num_input_lines
       read (unit=unit, fmt="(a)") line
       nin_count(1) = len(trim(line))
!       status = nf_put_vara_text (ncid, input_id, nin_start, nin_count, line)
       status = nf90_put_var (ncid, input_id, line, start=nin_start, count=nin_count)
       if (status /= nf90_noerr) call netcdf_error (status, ncid, input_id)
       nin_start(2) = nin_start(2) + 1
    end do
# endif
  end subroutine save_input

  subroutine define_vars (write_phi_vs_t, write_gvmus, write_gzvs)

    use mp, only: nproc
    use species, only: nspec
    use kt_grids, only: naky, nakx
    use run_parameters, only: fphi, fapar, fbpar
# ifdef NETCDF
    use netcdf, only: nf90_char, nf90_int, nf90_global
    use netcdf, only: nf90_def_var, nf90_put_att, nf90_enddef, nf90_put_var
    use netcdf, only: nf90_inq_libvers
    use netcdf_utils, only: netcdf_real
# endif

    implicit none

    logical, intent(in) :: write_phi_vs_t, write_gvmus, write_gzvs
# ifdef NETCDF
    character (5) :: ci
    character (20) :: datestamp, timestamp, timezone
    
    integer :: status

    fmode_dim(1) = nakx_dim
    fmode_dim(2) = naky_dim

    mode_dim (1) = nakx_dim
    mode_dim (2) = naky_dim
    mode_dim (3) = time_dim

    kx_dim (1) = nakx_dim
    kx_dim (2) = time_dim
    
    ky_dim (1) = naky_dim
    ky_dim (2) = time_dim
    
    om_dim (1) = ri_dim
    om_dim (2) = time_dim

    omega_dim (1) = ri_dim
    omega_dim (2) = nakx_dim
    omega_dim (3) = naky_dim
    omega_dim (4) = time_dim

    phase_dim (1) = ri_dim
    phase_dim (2) = nakx_dim
    phase_dim (3) = naky_dim
    
    nin_dim(1) = char200_dim
    nin_dim(2) = nlines_dim
    
    flux_dim (1) = nspec_dim
    flux_dim (2) = time_dim

    fluxk_dim (1) = nakx_dim
    fluxk_dim (2) = naky_dim
    fluxk_dim (3) = nspec_dim
    fluxk_dim (4) = time_dim

    fluxx_dim (1) = nakx_dim
    fluxx_dim (2) = nspec_dim
    fluxx_dim (3) = time_dim

    heat_dim (1) = nspec_dim
    heat_dim (2) = nheat_dim
    heat_dim (3) = time_dim

    heatk_dim (1) = nakx_dim
    heatk_dim (2) = naky_dim
    heatk_dim (3) = nspec_dim
    heatk_dim (4) = nheat_dim
    heatk_dim (5) = time_dim

    field_dim (1) = ri_dim
    field_dim (2) = naky_dim
    field_dim (3) = nakx_dim
    field_dim (4) = nttot_dim
    field_dim (5) = time_dim

    vmus_dim (1) = nvtot_dim
    vmus_dim (2) = nmu_dim
    vmus_dim (3) = nspec_dim
    vmus_dim (4) = time_dim

    zvs_dim (1) = nttot_dim
    zvs_dim (2) = nvtot_dim
    zvs_dim (3) = nspec_dim
    zvs_dim (4) = time_dim
    
    mom_t_dim (1) = ri_dim
    mom_t_dim (2) = nttot_dim
    mom_t_dim (3) = nakx_dim
    mom_t_dim (4) = naky_dim
    mom_t_dim (5) = nspec_dim
    mom_t_dim (6) = time_dim
    
    final_field_dim (1) = ri_dim
    final_field_dim (2) = nttot_dim
    final_field_dim (3) = nakx_dim
    final_field_dim (4) = naky_dim

    final_mom_dim (1) = ri_dim
    final_mom_dim (2) = nttot_dim
    final_mom_dim (3) = nakx_dim
    final_mom_dim (4) = naky_dim
    final_mom_dim (5) = nspec_dim

    loop_mom_dim (1) = ri_dim
    loop_mom_dim (2) = nakx_dim
    loop_mom_dim (3) = nspec_dim
    loop_mom_dim (4) = time_dim

    loop_phi_dim (1) = ri_dim
    loop_phi_dim (2) = nakx_dim
    loop_phi_dim (3) = time_dim

    mom_dim(1) = ri_dim
    mom_dim(2) = nakx_dim
    mom_dim(3) = naky_dim
    mom_dim(4) = nspec_dim
    mom_dim(5) = time_dim

    phi_corr_2pi_dim(1) = ri_dim
    phi_corr_2pi_dim(2) = nttot_dim
    phi_corr_2pi_dim(3) = naky_dim
    phi_corr_2pi_dim(4) = time_dim

    ! Write some useful general information such as the website,
    ! date and time into the NetCDF file
    status = nf90_put_att (ncid, nf90_global, 'title', 'GS2 Simulation Data')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, nf90_global, att='title')
    status = nf90_put_att (ncid, nf90_global, 'Conventions', &
         'http://gs2.sourceforge.net')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, nf90_global, att='Conventions')

    datestamp(:) = ' '
    timestamp(:) = ' '
    timezone(:) = ' '
    call date_and_time (datestamp, timestamp, timezone)
    
    status = nf90_def_var (ncid, 'code_info', nf90_char, char10_dim, code_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='code_info')
    status = nf90_put_att (ncid, code_id, 'long_name', 'GS2')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, code_id, att='long_name')

    ci = 'c1'
    status = nf90_put_att (ncid, code_id, trim(ci), 'Date: '//trim(datestamp))
    if (status /= nf90_noerr) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c2'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'Time: '//trim(timestamp)//' '//trim(timezone))
    if (status /= nf90_noerr) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c3'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'netCDF version '//trim(nf90_inq_libvers()))
    if (status /= nf90_noerr) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c4'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'Units are determined with respect to reference temperature (T_ref),')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c5'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'reference charge (q_ref), reference mass (mass_ref),')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c6'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'reference field (B_ref), and reference length (a_ref)')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c7'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'from which one may construct rho_ref and vt_ref/a,')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c8'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'which are the basic units of perpendicular length and time.')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c9'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'Macroscopic lengths are normalized to the minor radius.')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c10'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'The difference between rho (normalized minor radius) and rho (gyroradius)')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, code_id, att=ci)

    ci = 'c11'
    status = nf90_put_att (ncid, code_id, trim(ci), &
         'should be clear from the context in which they appear below.')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, code_id, att=ci)

    ! Write lots of input variables (e.g. nproc, nkx, nky)
    ! into the NetCDF file
    status = nf90_def_var (ncid, 'nproc', nf90_int, nproc_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='nproc')
    status = nf90_put_att (ncid, nproc_id, 'long_name', 'Number of processors')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, nproc_id, att='long_name')

    status = nf90_def_var (ncid, 'nmesh', netcdf_real, nmesh_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='nmesh')
    status = nf90_put_att (ncid, nmesh_id, 'long_name', 'Number of meshpoints')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, nmesh_id, att='long_name')

    status = nf90_def_var (ncid, 'nkx', nf90_int, nakx_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='nkx')
    status = nf90_def_var (ncid, 'nky', nf90_int, naky_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='nky')
    status = nf90_def_var (ncid, 'ntheta_tot', nf90_int, nttot_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='ntheta_tot')
    status = nf90_def_var (ncid, 'nspecies', nf90_int, nspec_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='nspecies')
    status = nf90_def_var (ncid, 'nmu', nf90_int, nmu_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='nmu')
    status = nf90_def_var (ncid, 'nvpa_tot', nf90_int, nvtot_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='nvpa_tot')

    status = nf90_def_var (ncid, 't', netcdf_real, time_dim, time_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='t')
    status = nf90_put_att (ncid, time_id, 'long_name', 'Time')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, time_id, att='long_name')
    status = nf90_put_att (ncid, time_id, 'units', 'L/vt')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, time_id, att='units')

    status = nf90_def_var (ncid, 'charge', nf90_int, nspec_dim, charge_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='charge')
    status = nf90_put_att (ncid, charge_id, 'long_name', 'Charge')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, charge_id, att='long_name')
    status = nf90_put_att (ncid, charge_id, 'units', 'q')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, charge_id, att='units')

    status = nf90_def_var (ncid, 'mass', netcdf_real, nspec_dim, mass_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='mass')
    status = nf90_put_att (ncid, mass_id, 'long_name', 'Atomic mass')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, mass_id, att='long_name')
    status = nf90_put_att (ncid, mass_id, 'units', 'm')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, mass_id, att='units')

    status = nf90_def_var (ncid, 'dens', netcdf_real, nspec_dim, dens_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='dens')
    status = nf90_put_att (ncid, dens_id, 'long_name', 'Density')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, dens_id, att='long_name')
    status = nf90_put_att (ncid, dens_id, 'units', 'n_e')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, dens_id, att='units')

    status = nf90_def_var (ncid, 'temp', netcdf_real, nspec_dim, temp_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='temp')
    status = nf90_put_att (ncid, temp_id, 'long_name', 'Temperature')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, temp_id, att='long_name')
    status = nf90_put_att (ncid, temp_id, 'units', 'T')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, temp_id, att='units')

    status = nf90_def_var (ncid, 'tprim', netcdf_real, nspec_dim, tprim_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='tprim')
    status = nf90_put_att (ncid, tprim_id, 'long_name', '-1/rho dT/drho')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, tprim_id, att='long_name')

    status = nf90_def_var (ncid, 'fprim', netcdf_real, nspec_dim, fprim_id) 
    if (status /= nf90_noerr) call netcdf_error (status, var='fprim')
    status = nf90_put_att (ncid, fprim_id, 'long_name', '-1/rho dn/drho')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, fprim_id, att='long_name')

    status = nf90_def_var (ncid, 'vnewk', netcdf_real, nspec_dim, vnewk_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='vnewk')
    status = nf90_put_att (ncid, vnewk_id, 'long_name', 'Collisionality')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, vnewk_id, att='long_name')
    status = nf90_put_att (ncid, vnewk_id, 'units', 'v_t/L')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, vnewk_id, att='units')
    
    status = nf90_def_var (ncid, 'type_of_species', nf90_int, nspec_dim, spec_type_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='type_of_species')

    status = nf90_def_var (ncid, 'theta0', netcdf_real, fmode_dim, theta0_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='theta0')
    status = nf90_put_att (ncid, theta0_id, 'long_name', 'Theta_0')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, theta0_id, att='long_name')

    status = nf90_def_var (ncid, 'kx', netcdf_real, nakx_dim, akx_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='kx')
    status = nf90_put_att (ncid, akx_id, 'long_name', 'kx rho')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, akx_id, att='long_name')

    status = nf90_def_var (ncid, 'ky', netcdf_real, naky_dim, aky_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='ky')
    status = nf90_put_att (ncid, aky_id, 'long_name', 'ky rho')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, aky_id, att='long_name')

    status = nf90_def_var (ncid, 'mu', netcdf_real, nmu_dim, mu_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='mu')
    status = nf90_def_var (ncid, 'vpa', netcdf_real, nvtot_dim, vpa_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='vpa')

    status = nf90_def_var (ncid, 'theta', netcdf_real, nttot_dim, theta_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='theta')

    status = nf90_def_var (ncid, 'bmag', netcdf_real, nttot_dim, bmag_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='bmag')
    status = nf90_put_att (ncid, bmag_id, 'long_name', '|B|(theta)')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, bmag_id, att='long_name')
    status = nf90_put_att (ncid, bmag_id, 'units', 'B_0')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, bmag_id, att='units')

    status = nf90_def_var (ncid, 'gradpar', netcdf_real, nttot_dim, gradpar_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='gradpar')
    status = nf90_def_var (ncid, 'gbdrift', netcdf_real, nttot_dim, gbdrift_id) 
    if (status /= nf90_noerr) call netcdf_error (status, var='gbdrift')
    status = nf90_def_var (ncid, 'gbdrift0', netcdf_real, nttot_dim, gbdrift0_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='gbdrift0')
    status = nf90_def_var (ncid, 'cvdrift', netcdf_real, nttot_dim, cvdrift_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='cvdrift')
    status = nf90_def_var (ncid, 'cvdrift0', netcdf_real, nttot_dim, cvdrift0_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='cvdrift0')
    status = nf90_def_var (ncid, 'cdrift', netcdf_real, nttot_dim, cdrift_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='cdrift')
    status = nf90_def_var (ncid, 'cdrift0', netcdf_real, nttot_dim, cdrift0_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='cdrift0')

    status = nf90_def_var (ncid, 'gds2', netcdf_real, nttot_dim, gds2_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='gds2')
    status = nf90_def_var (ncid, 'gds21', netcdf_real, nttot_dim, gds21_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='gds21')
    status = nf90_def_var (ncid, 'gds22', netcdf_real, nttot_dim, gds22_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='gds22')
    status = nf90_def_var (ncid, 'grho', netcdf_real, nttot_dim, grho_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='grho')
    status = nf90_def_var (ncid, 'jacob', netcdf_real, nttot_dim, jacob_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='jacob')

    status = nf90_def_var (ncid, 'q', netcdf_real, q_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='q')
    status = nf90_put_att (ncid, q_id, 'long_name', 'local safety factor')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, q_id, att='long_name')
    status = nf90_def_var (ncid, 'eps', netcdf_real, eps_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='eps')
    status = nf90_def_var (ncid, 'beta', netcdf_real, beta_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='beta')
    status = nf90_put_att (ncid, beta_id, 'long_name', 'reference beta')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, beta_id, att='long_name')
    status = nf90_def_var (ncid, 'shat', netcdf_real, shat_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='shat')
    status = nf90_put_att (ncid, shat_id, 'long_name', '(rho/q) dq/drho')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, shat_id, att='long_name')
    
    status = nf90_def_var (ncid, 'drhodpsi', netcdf_real, drhodpsi_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='drhodpsi')
    status = nf90_put_att (ncid, drhodpsi_id, 'long_name', 'drho/dPsi')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, drhodpsi_id, att='long_name')
    
    if (fphi > zero) then
       status = nf90_def_var (ncid, 'phi2', netcdf_real, time_dim, phi2_id)
       if (status /= nf90_noerr) call netcdf_error (status, var='phi2')
       status = nf90_put_att (ncid, phi2_id, 'long_name', '|Potential**2|')
       if (status /= nf90_noerr) &
            call netcdf_error (status, ncid, phi2_id, att='long_name')
       status = nf90_put_att (ncid, phi2_id, 'units', '(T/q rho/L)**2')
       if (status /= nf90_noerr) &
            call netcdf_error (status, ncid, phi2_id, att='units')
       
       status = nf90_def_var &
            (ncid, 'phi2_by_mode', netcdf_real, mode_dim, phi2_by_mode_id)
       if (status /= nf90_noerr) call netcdf_error (status, var='phi2_by_mode')
       if (nakx > 1) then
          status = nf90_def_var &
               (ncid, 'phi2_by_kx', netcdf_real, kx_dim, phi2_by_kx_id)
          if (status /= nf90_noerr) &
               call netcdf_error (status, var='phi2_by_kx')
       end if
       
       if (naky > 1) then
          status = nf90_def_var &
               (ncid, 'phi2_by_ky', netcdf_real, ky_dim, phi2_by_ky_id)
          if (status /= nf90_noerr) &
               call netcdf_error (status, var='phi2_by_ky')
       end if
       
       status = nf90_def_var (ncid, 'phi0', netcdf_real, omega_dim, phi0_id)
       if (status /= nf90_noerr) call netcdf_error (status, var='phi0')

       if (write_phi_vs_t) then
          status = nf90_def_var &
               (ncid, 'phi_vs_t', netcdf_real, field_dim, phi_vs_t_id)
          if (status /= nf90_noerr) call netcdf_error (status, var='phi_vs_t')
          status = nf90_put_att (ncid, phi_vs_t_id, 'long_name', 'Electrostatic Potential vs time')
          if (status /= nf90_noerr) call netcdf_error (status, ncid, phi_vs_t_id, att='long_name')
       end if
    end if

    if (write_gvmus) then
       status = nf90_def_var &
            (ncid, 'gvmus', netcdf_real, vmus_dim, gvmus_id)
       if (status /= nf90_noerr) call netcdf_error (status, var='gvmus')
       status = nf90_put_att (ncid, gvmus_id, 'long_name', &
            'guiding center distribution function averaged over real space')
       if (status /= nf90_noerr) call netcdf_error (status, ncid, gvmus_id, att='long_name')
    end if
    
    if (write_gzvs) then
       status = nf90_def_var &
            (ncid, 'gzvs', netcdf_real, zvs_dim, gzvs_id)
       if (status /= nf90_noerr) call netcdf_error (status, var='gzvs')
       status = nf90_put_att (ncid, gvmus_id, 'long_name', &
            'guiding center distribution function averaged over (kx,ky,mu)')
       if (status /= nf90_noerr) call netcdf_error (status, ncid, gzvs_id, att='long_name')
    end if


!        if (write_nl_flux) then
!           status = nf90_def_var (ncid, 'es_heat_par',  netcdf_real, flux_dim, es_heat_par_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_heat_par')
!           status = nf90_def_var (ncid, 'es_heat_perp', netcdf_real, flux_dim, es_heat_perp_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_heat_perp')
!           status = nf90_def_var (ncid, 'es_heat_flux', netcdf_real, flux_dim, es_heat_flux_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_heat_flux')
!           status = nf90_def_var (ncid, 'es_mom_flux',  netcdf_real, flux_dim, es_mom_flux_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_mom_flux')
!           status = nf90_def_var (ncid, 'es_part_flux', netcdf_real, flux_dim, es_part_flux_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_part_flux')
!           status = nf90_def_var (ncid, 'es_energy_exchange', netcdf_real, flux_dim, es_energy_exchange_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_energy_exchange')
!           status = nf90_def_var (ncid, 'es_heat_by_k', netcdf_real, fluxk_dim, es_heat_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_heat_by_k')
!           status = nf90_def_var (ncid, 'es_mom_by_k',  netcdf_real, fluxk_dim, es_mom_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_mom_by_k')
!           status = nf90_def_var (ncid, 'es_part_by_k', netcdf_real, fluxk_dim, es_part_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_part_by_k')
!           status = nf90_def_var (ncid, 'es_parmom_by_k', netcdf_real, fluxk_dim, es_parmom_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_parmom_by_k')
!           status = nf90_def_var (ncid, 'es_perpmom_by_k', netcdf_real, fluxk_dim, es_perpmom_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_perpmom_by_k')
!           status = nf90_def_var (ncid, 'es_mom0_by_k', netcdf_real, fluxk_dim, es_mom0_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_mom0_by_k')
!           status = nf90_def_var (ncid, 'es_mom1_by_k', netcdf_real, fluxk_dim, es_mom1_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='es_mom1_by_k')
!        end if
!        if (write_correlation) then
!           status = nf90_def_var (ncid, 'phi_corr_2pi',  netcdf_real, phi_corr_2pi_dim, phi_corr_2pi_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='phi_corr_2pi')
!        end if
!        if (write_eigenfunc) then
!           status = nf90_def_var (ncid, 'phi_norm',  netcdf_real, final_field_dim, phi_norm_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='phi_norm')
!        endif
!        status = nf90_def_var (ncid, 'phi', netcdf_real, final_field_dim, phi_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='phi')
!        status = nf90_put_att &
!             (ncid, phi_id, 'long_name', 'Electrostatic Potential')
!        if (status /= nf90_noerr) &
!             call netcdf_error (status, ncid, phi_id, att='long_name')
!        status = nf90_put_att (ncid, phi_id, 'idl_name', '!7U!6')
!        if (status /= nf90_noerr) &
!             call netcdf_error (status, ncid, phi_id, att='idl_name')
!        status = nf90_put_att (ncid, phi_id, 'units', 'T/q rho/L')
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, phi_id, att='units')
!        if (write_final_antot) then
!           status = nf90_def_var (ncid, 'antot', netcdf_real, final_field_dim, antot_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='antot')
!        endif
! !CMR
!        if (write_moments) then
!           status = nf90_def_var &
!                (ncid, 'ntot_t', netcdf_real, mom_t_dim, ntot_t_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='ntot_t')
!           status = nf90_put_att (ncid, ntot_t_id, 'long_name', 'Total perturbed density over time')
!           if (status /= nf90_noerr) call netcdf_error (status, ncid, ntot_t_id, att='long_name')
!        end if
! !CMRend
!     end if

!     if (fapar > zero) then
!        status = nf90_def_var (ncid, 'apar2', netcdf_real, time_dim, apar2_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='apar2')
!        status = nf90_def_var &
!             (ncid, 'apar2_by_mode', netcdf_real, mode_dim, apar2_by_mode_id)
!        if (status /= nf90_noerr) &
!             call netcdf_error (status, var='apar2_by_mode')
!        if (ntheta0 > 1) then
!           status = nf90_def_var &
!                (ncid, 'apar2_by_kx', netcdf_real, kx_dim, apar2_by_kx_id)
!           if (status /= nf90_noerr) &
!                call netcdf_error (status, var='apar2_by_kx')
!        end if
!        if (naky > 1) then
!           status = nf90_def_var &
!                (ncid, 'apar2_by_ky', netcdf_real, ky_dim, apar2_by_ky_id)
!           if (status /= nf90_noerr) &
!                call netcdf_error (status, var='apar2_by_ky')
!        end if

!        status = nf90_def_var (ncid, 'apar0', netcdf_real, omega_dim, apar0_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='apar0')
!        if (write_nl_flux) then
!           status = nf90_def_var &
!                (ncid,'apar_heat_flux',netcdf_real, flux_dim, apar_heat_flux_id)
!           if (status /= nf90_noerr) &
!                call netcdf_error (status, var='apar_heat_flux')
!           status = nf90_def_var &
!                (ncid, 'apar_heat_par', netcdf_real, flux_dim, apar_heat_par_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='apar_heat_par')
!           status = nf90_def_var (ncid, 'apar_heat_perp', netcdf_real, flux_dim, apar_heat_perp_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='apar_heat_perp')
!           status = nf90_def_var (ncid, 'apar_mom_flux',  netcdf_real, flux_dim, apar_mom_flux_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='apar_mom_flux')
!           status = nf90_def_var (ncid, 'apar_part_flux', netcdf_real, flux_dim, apar_part_flux_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='apar_part_flux')
!           status = nf90_def_var (ncid, 'apar_heat_by_k', netcdf_real, fluxk_dim, apar_heat_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='apar_heat_by_k')
!           status = nf90_def_var (ncid, 'apar_heat_by_x', netcdf_real, fluxx_dim, apar_heat_by_x_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='apar_heat_by_x')
!           status = nf90_def_var (ncid, 'apar_mom_by_k',  netcdf_real, fluxk_dim, apar_mom_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='apar_mom_by_k')
!           status = nf90_def_var (ncid, 'apar_part_by_k', netcdf_real, fluxk_dim, apar_part_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='apar_part_by_k')
!        end if
!        if (write_eigenfunc) then
!           status = nf90_def_var (ncid, 'apar_norm', netcdf_real, final_field_dim, apar_norm_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='apar_norm')
!        endif
!        status = nf90_def_var (ncid, 'apar', netcdf_real, final_field_dim, apar_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='apar')
!        if (write_final_antot) then
!           status = nf90_def_var (ncid, 'antota', netcdf_real, final_field_dim, antota_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='antota')
!        endif
!        if (write_apar_t) then
!           status = nf90_def_var (ncid, 'apar_t', netcdf_real, field_dim, apar_t_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='apar_t')
!           status = nf90_put_att (ncid, apar_t_id, 'long_name', 'Parallel Magnetic Potential over Time')
!           if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_t_id, att='long_name')
!        end if
!        status = nf90_put_att (ncid, apar2_by_mode_id, 'long_name', 'Apar squared')
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar2_by_mode_id, att='long_name')
!        status = nf90_put_att (ncid, apar_id, 'long_name', 'Parallel Magnetic Potential')
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_id, att='long_name')
!        status = nf90_put_att (ncid, apar_id, 'idl_name', '!6A!9!D#!N!6')
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_id, att='idl_name')
!        status = nf90_put_att (ncid, apar2_id, 'long_name', 'Total A_par squared')
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar2_id, att='long_name')
!     end if

!     if (fbpar > zero) then
!        status = nf90_def_var (ncid, 'bpar2', netcdf_real, time_dim, bpar2_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='bpar2')
!        status = nf90_def_var (ncid, 'bpar2_by_mode', netcdf_real, mode_dim, bpar2_by_mode_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='bpar2_by_mode')
!        if (ntheta0 > 1) then
!           status = nf90_def_var (ncid, 'bpar2_by_kx', netcdf_real, kx_dim, bpar2_by_kx_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='bpar2_by_kx')
!        end if
!        if (naky > 1) then
!           status = nf90_def_var (ncid, 'bpar2_by_ky', netcdf_real, ky_dim, bpar2_by_ky_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='bpar2_by_ky')
!        end if
!        status = nf90_def_var (ncid, 'bpar0', netcdf_real, omega_dim, bpar0_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='bpar0')
!        if (write_nl_flux) then
!           status = nf90_def_var (ncid, 'bpar_heat_flux', netcdf_real, flux_dim, bpar_heat_flux_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='bpar_heat_flux')
!           status = nf90_def_var (ncid, 'bpar_heat_par', netcdf_real, flux_dim, bpar_heat_par_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='bpar_heat_par')
!           status = nf90_def_var (ncid, 'bpar_heat_perp', netcdf_real, flux_dim, bpar_heat_perp_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='bpar_heat_perp')
!           status = nf90_def_var (ncid, 'bpar_mom_flux', netcdf_real, flux_dim, bpar_mom_flux_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='bpar_mom_flux')
!           status = nf90_def_var (ncid, 'bpar_part_flux', netcdf_real, flux_dim, bpar_part_flux_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='bpar_part_flux')
!           status = nf90_def_var (ncid, 'bpar_heat_by_k', netcdf_real, fluxk_dim, bpar_heat_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='bpar_heat_by_k')
!           status = nf90_def_var (ncid, 'bpar_mom_by_k', netcdf_real, fluxk_dim, bpar_mom_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='bpar_mom_by_k')
!           status = nf90_def_var (ncid, 'bpar_part_by_k', netcdf_real, fluxk_dim, bpar_part_by_k_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='bpar_part_by_k')
!        end if
!        if (write_eigenfunc) then
!           status = nf90_def_var (ncid, 'bpar_norm', netcdf_real, final_field_dim, bpar_norm_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='bpar_norm')
!        endif
!        status = nf90_def_var (ncid, 'bpar', netcdf_real, final_field_dim, bpar_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='bpar')
!        if (write_final_antot) then
!           status = nf90_def_var (ncid, 'antotp', netcdf_real, final_field_dim, antotp_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='antotp')
!        endif
!        if (write_bpar_t) then
!           status = nf90_def_var (ncid, 'bpar_t', netcdf_real, field_dim, bpar_t_id)
!           if (status /= nf90_noerr) call netcdf_error (status, var='bpar_t')
!           status = nf90_put_att (ncid, bpar_t_id, 'long_name', 'delta B Parallel over time')
!           if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_t_id, att='long_name')
!        end if

!        status = nf90_put_att (ncid, bpar2_by_mode_id, 'long_name', 'A_perp squared')
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar2_by_mode_id, att='long_name')
!        status = nf90_put_att (ncid, bpar_id, 'long_name', 'delta B Parallel')
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_id, att='long_name')
!        status = nf90_put_att (ncid, bpar_id, 'idl_name', '!6B!9!D#!N!6')
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_id, att='idl_name')
!        status = nf90_put_att (ncid, bpar2_id, 'long_name', 'Total A_perp squared')
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar2_id, att='long_name')
!     end if

!     status = nf90_def_var (ncid, 'phase', netcdf_real, phase_dim, phase_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='phase')
!     status = nf90_put_att (ncid, phase_id, 'long_name', 'Normalizing phase')
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, phase_id, att='long_name')

! !    status = nf90_def_var (ncid, 'phtot', netcdf_real, mode_dim, phtot_id)
! !    if (status /= nf90_noerr) call netcdf_error (status, var='phtot')
! !    status = nf90_def_var (ncid, 'dmix', netcdf_real, mode_dim, dmix_id)
! !    if (status /= nf90_noerr) call netcdf_error (status, var='dmix')
! !    status = nf90_def_var (ncid, 'kperpnorm', netcdf_real, mode_dim, kperpnorm_id)
! !    if (status /= nf90_noerr) call netcdf_error (status, var='kperpnorm')

!     if (write_omega) then
!        status = nf90_def_var (ncid, 'omega', netcdf_real, omega_dim, omega_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='omega')
!        status = nf90_def_var (ncid, 'omegaavg', netcdf_real, omega_dim, omegaavg_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='omegaavg')
!     end if

!     if (write_nl_flux) then
!        status = nf90_def_var (ncid, 'hflux_tot', netcdf_real, time_dim, hflux_tot_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='kflux_tot')
!        status = nf90_def_var (ncid, 'vflux_tot', netcdf_real, time_dim, vflux_tot_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='vflux_tot')
!        status = nf90_def_var (ncid, 'zflux_tot', netcdf_real, time_dim, zflux_tot_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='zflux_tot')
!     end if

!     status = nf90_def_var (ncid, 'epar', netcdf_real, final_field_dim, epar_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='epar')

! !    <phi>

!     status = nf90_def_var (ncid, 'ntot2', netcdf_real, flux_dim,  ntot2_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='ntot2')
!     status = nf90_def_var (ncid, 'ntot2_by_mode', netcdf_real, fluxk_dim, ntot2_by_mode_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='ntot2_by_mode')
!     status = nf90_def_var (ncid, 'tpar2_by_mode', netcdf_real, fluxk_dim, tpar2_by_mode_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='tpar2_by_mode')
!     status = nf90_def_var (ncid, 'tperp2_by_mode', netcdf_real, fluxk_dim, tperp2_by_mode_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='tperp2_by_mode')

!     status = nf90_def_var (ncid, 'ntot20', netcdf_real, flux_dim,  ntot20_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='ntot20')
!     status = nf90_put_att (ncid, ntot20_id, 'long_name', 'Density**2 at theta=0')
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, ntot20_id, att='long_name')
!     status = nf90_def_var (ncid, 'ntot20_by_mode', netcdf_real, fluxk_dim, ntot20_by_mode_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='ntot20_by_mode')

!     status = nf90_def_var (ncid, 'ntot', netcdf_real, final_mom_dim, ntot_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='ntot')
!     status = nf90_def_var (ncid, 'density', netcdf_real, final_mom_dim, density_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='density')
!     status = nf90_def_var (ncid, 'upar', netcdf_real, final_mom_dim, upar_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='upar')
!     status = nf90_def_var (ncid, 'tpar', netcdf_real, final_mom_dim, tpar_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='tpar')
!     status = nf90_def_var (ncid, 'tperp', netcdf_real, final_mom_dim, tperp_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='tperp')
!     status = nf90_def_var (ncid, 'qparflux', netcdf_real, final_mom_dim, qparflux_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='qparflux')
!     status = nf90_def_var (ncid, 'pperpj1', netcdf_real, final_mom_dim, pperpj1_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='pperpj1')
!     status = nf90_def_var (ncid, 'qpperpj1', netcdf_real, final_mom_dim, qpperpj1_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='qpperpj1')

!     status = nf90_def_var (ncid, 'phi00', netcdf_real, loop_phi_dim, phi00_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='phi00')
!     status = nf90_def_var (ncid, 'ntot00', netcdf_real, loop_mom_dim, ntot00_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='ntot00')
!     status = nf90_def_var (ncid, 'density00', netcdf_real, loop_mom_dim, density00_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='density00')
!     status = nf90_def_var (ncid, 'upar00', netcdf_real, loop_mom_dim, upar00_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='upar00')
!     status = nf90_def_var (ncid, 'tpar00', netcdf_real, loop_mom_dim, tpar00_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='tpar00')
!     status = nf90_def_var (ncid, 'tperp00', netcdf_real, loop_mom_dim, tperp00_id)
!     if (status /= nf90_noerr) call netcdf_error (status, var='tperp00')
    
!     ! RN> not guiding center moments
!     if (write_full_moments_notgc) then
!        status = nf90_def_var (ncid, 'ntot0', netcdf_real, mom_dim, ntot0_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='ntot0')
!        status = nf90_def_var (ncid, 'density0', netcdf_real, mom_dim, density0_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='density0')
!        status = nf90_def_var (ncid, 'upar0', netcdf_real, mom_dim, upar0_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='upar0')
!        status = nf90_def_var (ncid, 'tpar0', netcdf_real, mom_dim, tpar0_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='tpar0')
!        status = nf90_def_var (ncid, 'tperp0', netcdf_real, mom_dim, tperp0_id)
!        if (status /= nf90_noerr) call netcdf_error (status, var='tperp0')
!     endif

    status = nf90_def_var (ncid, 'input_file', nf90_char, nin_dim, input_id)
    if (status /= nf90_noerr) call netcdf_error (status, var='input_file')
    status = nf90_put_att (ncid, input_id, 'long_name', 'Input file')
    if (status /= nf90_noerr) call netcdf_error (status, ncid, input_id, att='long_name')

    status = nf90_enddef (ncid)  ! out of definition mode
    if (status /= nf90_noerr) call netcdf_error (status)

    status = nf90_put_var (ncid, nproc_id, nproc)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, nproc_id)

# endif
  end subroutine define_vars

!   subroutine nc_eigenfunc (phase)

!     use convert, only: c2r
!     use run_parameters, only: fphi, fapar, fbpar
!     use fields_arrays, only: phiold!, aparold, bparold
!     use theta_grid, only: ntgrid
!     use kt_grids, only: naky, ntheta0
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var
! # endif
!     complex, dimension(:,:), intent (in) :: phase
! # ifdef NETCDF
!     complex, dimension(-ntgrid:ntgrid, ntheta0, naky) :: tmp
!     real, dimension(2, ntheta0, naky) :: ri2
!     real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
!     integer :: status, ig

!     call c2r (phase, ri2)
!     status = nf90_put_var (ncid, phase_id, ri2)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, phase_id)

!     if (fphi > zero) then
!        do ig = -ntgrid, ntgrid
!           tmp(ig,:,:) = phiold(ig,:,:)/phase(:,:)
!        end do
!        call c2r (tmp, ri3)
!        status = nf90_put_var(ncid, phi_norm_id, ri3)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, phi_norm_id)
!     end if

! !     if (fapar > zero) then
! !        do ig = -ntgrid, ntgrid
! !           tmp(ig,:,:) = aparold(ig,:,:)/phase(:,:)
! !        end do
! !        call c2r (tmp, ri3)
! !        status = nf90_put_var(ncid, apar_norm_id, ri3)
! !        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_norm_id)
! !     end if

! !     if (fbpar > zero) then
! !        do ig = -ntgrid, ntgrid
! !           tmp(ig,:,:) = bparold(ig,:,:)/phase(:,:)
! !        end do
! !        call c2r (tmp, ri3)
! !        status = nf90_put_var(ncid, bpar_norm_id, ri3)
! !        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_norm_id)
! !     end if
! # endif
!   end subroutine nc_eigenfunc

!   subroutine nc_write_fields (nout, phinew, aparnew, bparnew)
!     use convert, only: c2r
!     use run_parameters, only: fphi, fapar, fbpar
!     use theta_grid, only: ntgrid
!     use kt_grids, only: naky, ntheta0
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var
! # endif
!     complex, dimension (:,:,:), intent (in) :: phinew, aparnew, bparnew
!     integer, intent (in) :: nout
! !    real, dimension (2, 2*ntgrid+1, ntheta0, naky, 1) :: ri4
!     real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
!     integer, dimension (5) :: start5, count5
!     integer :: status
! # ifdef NETCDF
!     start5(1) = 1
!     start5(2) = 1
!     start5(3) = 1
!     start5(4) = 1
!     start5(5) = nout
    
!     count5(1) = 2
!     count5(2) = 2*ntgrid+1
!     count5(3) = ntheta0
!     count5(4) = naky
!     count5(5) = 1

!     if (fapar > zero) then
!        call c2r (aparnew, ri3)
!        status = nf90_put_var(ncid, apar_t_id, ri3, start=start5, count=count5)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_t_id)
!     end if

!     if (fbpar > zero) then
!        call c2r (bparnew, ri3)
!        status = nf90_put_var(ncid, bpar_t_id, ri3, start=start5, count=count5)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_t_id)
!     end if
! # endif
!   end subroutine nc_write_fields

!   subroutine nc_write_moments (nout, ntot)
!     use convert, only: c2r
!     use theta_grid, only: ntgrid
!     use kt_grids, only: naky, ntheta0
!     use species, only: nspec
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var
! # endif
!     complex, dimension (:,:,:,:), intent (in) :: ntot
!     integer, intent (in) :: nout
!     real, dimension (2, 2*ntgrid+1, ntheta0, naky, nspec) :: ri4
!     integer, dimension (6) :: start6, count6
!     integer :: status
! # ifdef NETCDF
!     start6(1) = 1
!     start6(2) = 1
!     start6(3) = 1
!     start6(4) = 1
!     start6(5) = 1
!     start6(6) = nout
    
!     count6(1) = 2
!     count6(2) = 2*ntgrid+1
!     count6(3) = ntheta0
!     count6(4) = naky
!     count6(5) = nspec
!     count6(6) = 1

!     call c2r (ntot, ri4)
!     status = nf90_put_var(ncid, ntot_t_id, ri4, start=start6, count=count6)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, ntot_t_id)

! # endif
!   end subroutine nc_write_moments

!   subroutine nc_final_fields

!     use convert, only: c2r
!     use run_parameters, only: fphi, fapar, fbpar
!     use fields_arrays, only: phiold!, aparold, bparold
!     use theta_grid, only: ntgrid
!     use kt_grids, only: naky, ntheta0
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var

!     real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
!     integer :: status

!     if (fphi > zero) then
!        call c2r (phiold, ri3)
!        status = nf90_put_var (ncid, phi_id, ri3)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, phi_id)
!     end if

! !     if (fapar > zero) then
! !        call c2r (aparold, ri3)
! !        status = nf90_put_var (ncid, apar_id, ri3)
! !        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_id)
! !     end if

! !     if (fbpar > zero) then
! !        call c2r (bparold, ri3)
! !        status = nf90_put_var (ncid, bpar_id, ri3)
! !        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_id)
! !     end if
! # endif
!   end subroutine nc_final_fields

!   subroutine nc_final_epar (epar)

!     use convert, only: c2r
!     use run_parameters, only: fphi, fapar, fbpar
!     use theta_grid, only: ntgrid
!     use kt_grids, only: naky, ntheta0
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var
! # endif
!     complex, dimension (:,:,:), intent (in) :: epar
! # ifdef NETCDF
!     real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
!     integer :: status

!     call c2r (epar, ri3)
!     status = nf90_put_var (ncid, epar_id, ri3)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, epar_id)
! # endif
!   end subroutine nc_final_epar

!   subroutine nc_final_moments (ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)

!     use convert, only: c2r
!     use theta_grid, only: ntgrid
!     use kt_grids, only: naky, ntheta0
!     use species, only: nspec
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var
! # endif
!     complex, dimension (:,:,:,:), intent (in) :: ntot, density, upar, tpar, tperp
!     complex, dimension (:,:,:,:), intent (in) :: qparflux, pperpj1, qpperpj1
! # ifdef NETCDF
!     real, dimension (2, 2*ntgrid+1, ntheta0, naky, nspec) :: ri4
!     integer :: status

!     call c2r (ntot, ri4)
!     status = nf90_put_var (ncid, ntot_id, ri4)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, ntot_id)

!     call c2r (density, ri4)
!     status = nf90_put_var (ncid, density_id, ri4)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, density_id)

!     call c2r (upar, ri4)
!     status = nf90_put_var (ncid, upar_id, ri4)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, upar_id)

!     call c2r (tpar, ri4)
!     status = nf90_put_var (ncid, tpar_id, ri4)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, tpar_id)

!     call c2r (tperp, ri4)
!     status = nf90_put_var (ncid, tperp_id, ri4)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, tperp_id)

!     call c2r (qparflux, ri4)
!     status = nf90_put_var (ncid, qparflux_id, ri4)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, qparflux_id)

!     call c2r (pperpj1, ri4)
!     status = nf90_put_var (ncid, pperpj1_id, ri4)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, pperpj1_id)

!     call c2r (qpperpj1, ri4)
!     status = nf90_put_var (ncid, qpperpj1_id, ri4)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, qpperpj1_id)

! # endif
!   end subroutine nc_final_moments

!   subroutine nc_loop_moments (nout, ntot2, ntot2_by_mode, ntot20, ntot20_by_mode, &
!        phi00, ntot00, density00, upar00, tpar00, tperp00, tpar2_by_mode, tperp2_by_mode)

!     use convert, only: c2r
!     use theta_grid, only: ntgrid
!     use kt_grids, only: naky, ntheta0
!     use species, only: nspec
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var
! # endif
!     integer, intent (in) :: nout
!     real, dimension (:), intent (in) :: ntot2, ntot20
!     real, dimension (:,:,:), intent (in) :: ntot2_by_mode, ntot20_by_mode
!     real, dimension (:,:,:), intent (in) :: tpar2_by_mode, tperp2_by_mode
!     complex, dimension (:), intent (in) :: phi00
!     complex, dimension (:,:), intent (in) :: ntot00, density00, upar00, tpar00, tperp00
! # ifdef NETCDF
!     real, dimension (2, ntheta0, nspec) :: ri2
!     real, dimension (2, ntheta0) :: ri1
!     integer, dimension (2) :: start, count
!     integer, dimension (3) :: start3, count3
!     integer, dimension (4) :: start00, count00, start4, count4
!     integer :: status

!     start00(1) = 1
!     start00(2) = 1
!     start00(3) = 1
!     start00(4) = nout
    
!     count00(1) = 2
!     count00(2) = ntheta0
!     count00(3) = nspec
!     count00(4) = 1

!     start3(1) = 1
!     start3(2) = 1
!     start3(3) = nout
    
!     count3(1) = 2
!     count3(2) = ntheta0
!     count3(3) = 1

!     start(1) = 1
!     start(2) = nout
    
!     count(1) = nspec
!     count(2) = 1

!     start4(1) = 1
!     start4(2) = 1
!     start4(3) = 1
!     start4(4) = nout
    
!     count4(1) = ntheta0
!     count4(2) = naky
!     count4(3) = nspec
!     count4(4) = 1

!     status = nf90_put_var (ncid, ntot2_id, ntot2, start=start, count=count)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, ntot2_id)
!     status = nf90_put_var (ncid, ntot2_by_mode_id, ntot2_by_mode, start=start4, count=count4)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, ntot2_by_mode_id)
!     status = nf90_put_var (ncid, tpar2_by_mode_id, tpar2_by_mode, start=start4, count=count4)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, tpar2_by_mode_id)
!     status = nf90_put_var (ncid, tperp2_by_mode_id, tperp2_by_mode, start=start4, count=count4)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, tperp2_by_mode_id)

!     status = nf90_put_var (ncid, ntot20_id, ntot20, start=start, count=count)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, ntot20_id)
!     status = nf90_put_var (ncid, ntot20_by_mode_id, ntot20_by_mode, start=start4, count=count4)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, ntot20_by_mode_id)

!     call c2r (phi00, ri1)
!     status = nf90_put_var (ncid, phi00_id, ri1, start=start3, count=count3)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, phi00_id)
    
!     call c2r (ntot00, ri2)
!     status = nf90_put_var (ncid, ntot00_id, ri2, start=start00, count=count00)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, ntot00_id)
    
!     call c2r (density00, ri2)
!     status = nf90_put_var (ncid, density00_id, ri2, start=start00, count=count00)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, density00_id)
    
!     call c2r (upar00, ri2)
!     status = nf90_put_var (ncid, upar00_id, ri2, start=start00, count=count00)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, upar00_id)
    
!     call c2r (tpar00, ri2)
!     status = nf90_put_var (ncid, tpar00_id, ri2, start=start00, count=count00)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, tpar00_id)

!     call c2r (tperp00, ri2)
!     status = nf90_put_var (ncid, tperp00_id, ri2, start=start00, count=count00)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, tperp00_id)
! # endif
!   end subroutine nc_loop_moments

!   subroutine nc_final_an (antot, antota, antotp)

!     use convert, only: c2r
!     use run_parameters, only: fphi, fapar, fbpar
!     use theta_grid, only: ntgrid
!     use kt_grids, only: naky, ntheta0
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var
! # endif
!     complex, dimension (:,:,:) :: antot, antota, antotp ! intent?
! # ifdef NETCDF
!     real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
!     integer :: status

!     if (fphi > zero) then
!        call c2r (antot, ri3)
!        status = nf90_put_var(ncid, antot_id, ri3)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, antot_id)
!     end if

!     if (fapar > zero) then
!        call c2r (antota, ri3)
!        status = nf90_put_var(ncid, antota_id, ri3)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, antota_id)
!     end if

!     if (fbpar > zero) then
!        call c2r (antotp, ri3)
!        status = nf90_put_var(ncid, antotp_id, ri3)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, antotp_id)
!     end if
! # endif
!   end subroutine nc_final_an

!   subroutine nc_qflux (nout, qheat, qmheat, qbheat, &
!        heat_par,  mheat_par,  bheat_par, &
!        heat_perp, mheat_perp, bheat_perp, &
!        heat_fluxes, mheat_fluxes, bheat_fluxes, x_qmflux, hflux_tot, &
!        energy_exchange)

!     use species, only: nspec
!     use kt_grids, only: naky, ntheta0
!     use run_parameters, only: fphi, fapar, fbpar
! !    use nf90_mod, only: nf90_put_vara, nf90_put_var1
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var
! # endif
!     integer, intent (in) :: nout
!     real, dimension (:,:,:), intent (in) :: qheat, qmheat, qbheat
!     real, dimension (:), intent (in) :: heat_par, mheat_par, bheat_par
!     real, dimension (:), intent (in) :: heat_perp, mheat_perp, bheat_perp
!     real, dimension (:), intent (in) :: heat_fluxes, mheat_fluxes, bheat_fluxes
!     real, dimension (:), intent (in) :: energy_exchange
!     real, dimension (:,:), intent (in) :: x_qmflux
!     real, intent (in) :: hflux_tot
! # ifdef NETCDF
!     integer, dimension (2) :: start, count
!     integer, dimension (3) :: start3, count3
!     integer, dimension (4) :: start4, count4
!     integer :: status

!     start3(1) = 1
!     start3(2) = 1
!     start3(3) = nout

!     count3(1) = ntheta0
!     count3(2) = nspec
!     count3(3) = 1    

!     start4(1) = 1
!     start4(2) = 1
!     start4(3) = 1
!     start4(4) = nout
    
!     count4(1) = ntheta0
!     count4(2) = naky
!     count4(3) = nspec
!     count4(4) = 1

!     start(1) = 1
!     start(2) = nout
    
!     count(1) = nspec
!     count(2) = 1

!     if (fphi > zero) then
!        status = nf90_put_var (ncid, es_heat_flux_id, heat_fluxes, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_heat_flux_id)
!        status = nf90_put_var (ncid, es_heat_par_id, heat_par, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_heat_par_id)
!        status = nf90_put_var (ncid, es_heat_perp_id, heat_perp, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_heat_perp_id)
!        status = nf90_put_var (ncid, es_heat_by_k_id, qheat, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_heat_by_k_id)
!        status = nf90_put_var (ncid, es_energy_exchange_id, energy_exchange, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_energy_exchange_id)
!     end if

!     if (fapar > zero) then
!        status = nf90_put_var (ncid, apar_heat_flux_id, mheat_fluxes, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_heat_flux_id)
!        status = nf90_put_var (ncid, apar_heat_par_id, mheat_par, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_heat_par_id)
!        status = nf90_put_var (ncid, apar_heat_perp_id, mheat_perp, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_heat_perp_id)
!        status = nf90_put_var (ncid, apar_heat_by_k_id, qmheat, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_heat_by_k_id)
!        status = nf90_put_var (ncid, apar_heat_by_x_id, qmheat, start=start3, count=count3)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_heat_by_x_id)
!     end if

!     if (fbpar > zero) then
!        status = nf90_put_var (ncid, bpar_heat_flux_id, bheat_fluxes, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_heat_flux_id)
!        status = nf90_put_var (ncid, bpar_heat_par_id, bheat_par, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_heat_par_id)
!        status = nf90_put_var (ncid, bpar_heat_perp_id, bheat_perp, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_heat_perp_id)
!        status = nf90_put_var (ncid, bpar_heat_by_k_id, qbheat, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_heat_by_k_id)
!     end if
    
!     status = nf90_put_var (ncid, hflux_tot_id, hflux_tot, start=(/ nout /))
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, hflux_tot_id)
! # endif
!   end subroutine nc_qflux

!   subroutine nc_pflux (nout, pflux, pmflux, pbflux, &
!        part_fluxes, mpart_fluxes, bpart_fluxes, zflux_tot)

!     use species, only: nspec
!     use kt_grids, only: naky, ntheta0
!     use run_parameters, only: fphi, fapar, fbpar
! !    use nf90_mod, only: nf90_put_vara, nf90_put_var1
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var
! # endif
!     integer, intent (in) :: nout
!     real, dimension (:,:,:), intent (in) :: pflux, pmflux, pbflux
!     real, dimension(:), intent (in) :: part_fluxes, mpart_fluxes, bpart_fluxes
!     real, intent (in) :: zflux_tot
! # ifdef NETCDF
!     integer, dimension (2) :: start, count
!     integer, dimension (4) :: start4, count4
!     integer :: status

!     start4(1) = 1
!     start4(2) = 1
!     start4(3) = 1
!     start4(4) = nout
    
!     count4(1) = ntheta0
!     count4(2) = naky
!     count4(3) = nspec
!     count4(4) = 1

!     start(1) = 1
!     start(2) = nout
    
!     count(1) = nspec
!     count(2) = 1

!     if (fphi > zero) then
!        status = nf90_put_var (ncid, es_part_flux_id, part_fluxes, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_part_flux_id)
!        status = nf90_put_var (ncid, es_part_by_k_id, pflux, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_part_by_k_id)
!     end if

!     if (fapar > zero) then
!        status = nf90_put_var (ncid, apar_part_flux_id, mpart_fluxes, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_part_flux_id)
!        status = nf90_put_var (ncid, apar_part_by_k_id, pmflux, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_part_by_k_id)
!     end if

!     if (fbpar > zero) then
!        status = nf90_put_var (ncid, bpar_part_flux_id, bpart_fluxes, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_part_flux_id)
!        status = nf90_put_var (ncid, bpar_part_by_k_id, pbflux, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_part_by_k_id)
!     end if
    
!     status = nf90_put_var (ncid, zflux_tot_id, zflux_tot, start=(/ nout /))
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, zflux_tot_id)
! # endif
!   end subroutine nc_pflux

!   subroutine nc_vflux (nout, vflux, vmflux, vbflux, &
!        mom_fluxes, mmom_fluxes, bmom_fluxes, vflux_tot, &
!        vflux_par, vflux_perp, vflux0, vflux1)

!     use species, only: nspec
!     use kt_grids, only: naky, ntheta0
!     use run_parameters, only: fphi, fapar, fbpar
! !    use nf90_mod, only: nf90_put_vara, nf90_put_var1
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var
! # endif
!     integer, intent (in) :: nout
!     real, dimension (:,:,:), intent (in) :: vflux, vmflux, vbflux
!     real, dimension (:,:,:), intent (in) :: vflux_par, vflux_perp
!     real, dimension (:,:,:), intent (in) :: vflux0, vflux1
!     real, dimension(:), intent (in) :: mom_fluxes, mmom_fluxes, bmom_fluxes
!     real, intent (in) :: vflux_tot
! # ifdef NETCDF
!     integer, dimension (2) :: start, count
!     integer, dimension (4) :: start4, count4
!     integer :: status

!     start4(1) = 1
!     start4(2) = 1
!     start4(3) = 1
!     start4(4) = nout
    
!     count4(1) = ntheta0
!     count4(2) = naky
!     count4(3) = nspec
!     count4(4) = 1

!     start(1) = 1
!     start(2) = nout
    
!     count(1) = nspec
!     count(2) = 1

!     if (fphi > zero) then
!        status = nf90_put_var (ncid, es_mom_flux_id, mom_fluxes, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_mom_flux_id)
!        status = nf90_put_var (ncid, es_mom_by_k_id, vflux, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_mom_by_k_id)
!        status = nf90_put_var (ncid, es_parmom_by_k_id, vflux_par, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_parmom_by_k_id)
!        status = nf90_put_var (ncid, es_perpmom_by_k_id, vflux_perp, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_perpmom_by_k_id)
!        status = nf90_put_var (ncid, es_mom0_by_k_id, vflux0, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_mom0_by_k_id)
!        status = nf90_put_var (ncid, es_mom1_by_k_id, vflux1, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, es_mom1_by_k_id)
!     end if

!     if (fapar > zero) then
!        status = nf90_put_var (ncid, apar_mom_flux_id, mmom_fluxes, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_mom_flux_id)
!        status = nf90_put_var (ncid, apar_mom_by_k_id, vmflux, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_mom_by_k_id)
!     end if

!     if (fbpar > zero) then
!        status = nf90_put_var (ncid, bpar_mom_flux_id, bmom_fluxes, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_mom_flux_id)
!        status = nf90_put_var (ncid, bpar_mom_by_k_id, vbflux, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_mom_by_k_id)
!     end if
    
!     status = nf90_put_var (ncid, vflux_tot_id, vflux_tot, start=(/ nout /))
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, vflux_tot_id)
! # endif
!   end subroutine nc_vflux

!   subroutine nc_loop_corr (nout,phi_2pi_corr)

!     use convert, only: c2r
!     use run_parameters, only: fphi
!     use kt_grids, only: ntheta0, naky, jtwist_out
!     use theta_grid, only: ntgrid
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var
! # endif
!     integer, intent (in) :: nout
!     complex, dimension (-ntgrid:,:), intent (in) :: phi_2pi_corr
! # ifdef NETCDF
!     integer, dimension (4) :: start4, count4
!     integer :: status

!     real, dimension (2, size(phi_2pi_corr,1), size(phi_2pi_corr,2)) :: ri2

!     start4(1) = 1
!     start4(2) = 1
!     start4(3) = 1
!     start4(4) = nout

!     count4(1) = 2
!     count4(2) = 2*ntgrid+1
!     count4(3) = naky
!     count4(4) = 1

!     if (fphi > zero) then
!        call c2r (phi_2pi_corr, ri2)
!        status = nf90_put_var (ncid, phi_corr_2pi_id, ri2, start=start4, count=count4)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, phi_corr_2pi_id)
!     end if

! # endif
!   end subroutine nc_loop_corr

  subroutine write_time_nc (nout, time)

# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif

    implicit none

    integer, intent (in) :: nout
    real, intent (in) :: time

# ifdef NETCDF
    integer :: status

    status = nf90_put_var (ncid, time_id, time, start=(/ nout /))
    if (status /= nf90_noerr) call netcdf_error (status, ncid, time_id)
# endif

  end subroutine write_time_nc

  subroutine write_phi_nc (nout, phi)

    use convert, only: c2r
    use theta_grid, only: ntgrid
    use kt_grids, only: nakx, naky
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif

    implicit none

    integer, intent (in) :: nout
    complex, dimension (:,:,-ntgrid:), intent (in) :: phi

# ifdef NETCDF
    integer :: status
    integer, dimension (5) :: start, count
    real, dimension (:,:,:,:), allocatable :: phi_ri

    start = 1
    start(5) = nout
    count(1) = 2
    count(2) = naky
    count(3) = nakx
    count(4) = 2*ntgrid+1
    count(5) = 1

    allocate (phi_ri(2, naky, nakx, 2*ntgrid+1))
    call c2r (phi, phi_ri)
    status = nf90_put_var (ncid, phi_vs_t_id, phi_ri, start=start, count=count)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, phi_vs_t_id)
    deallocate (phi_ri)
# endif

  end subroutine write_phi_nc

  subroutine write_gvmus_nc (nout, g)

    use vpamu_grids, only: nvgrid, nmu
    use species, only: nspec
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif

    implicit none

    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: g

# ifdef NETCDF
    integer :: status
    integer, dimension (4) :: start, count

    start(1) = 1
    start(2:3) = 1
    start(4) = nout
    count(1) = 2*nvgrid+1
    count(2) = nmu
    count(3) = nspec
    count(4) = 1

    status = nf90_put_var (ncid, gvmus_id, g, start=start, count=count)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, gvmus_id)
# endif

  end subroutine write_gvmus_nc

  subroutine write_gzvs_nc (nout, g)

    use theta_grid, only: ntgrid
    use vpamu_grids, only: nvgrid
    use species, only: nspec
# ifdef NETCDF
    use netcdf, only: nf90_put_var
# endif

    implicit none

    integer, intent (in) :: nout
    real, dimension (:,:,:), intent (in) :: g

# ifdef NETCDF
    integer :: status
    integer, dimension (4) :: start, count

    start(1:3) = 1
    start(4) = nout
    count(1) = 2*ntgrid+1
    count(2) = 2*nvgrid+1
    count(3) = nspec
    count(4) = 1

    status = nf90_put_var (ncid, gzvs_id, g, start=start, count=count)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, gzvs_id)
# endif

  end subroutine write_gzvs_nc

!   subroutine nc_loop (nout, time, fluxfac, &
!        phi0,   phi2,   phi2_by_mode, &
!        apar0,  apar2,  apar2_by_mode, &
!        bpar0, bpar2, bpar2_by_mode, &
!        omega, omegaavg, woutunits, phitot, write_omega)

!     use run_parameters, only: fphi, fapar, fbpar
!     use kt_grids, only: naky, ntheta0
!     use theta_grid, only: ntgrid
!     use species, only: nspec
!     use convert, only: c2r
!     use fields_arrays, only: phiold!, aparold, bparold

! !    use nf90_mod, only: nf90_put_var1, nf90_put_vara
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var, nf90_sync
! # endif

!     integer, intent (in) :: nout
!     real, intent (in) :: time, phi2, apar2, bpar2
!     real, dimension (:), intent (in) :: fluxfac, woutunits
!     complex, dimension(:,:), intent (in) :: phi0, apar0, bpar0, omega, omegaavg
!     real, dimension(:,:), intent (in) :: phi2_by_mode, apar2_by_mode, bpar2_by_mode, phitot
!     logical :: write_omega
! # ifdef NETCDF
!     real, dimension (ntheta0) :: field2_by_kx
!     real, dimension (naky) :: field2_by_ky
!     real, dimension (2, ntheta0, naky) :: ri2
!     real, dimension (ntheta0, naky, nspec) :: tmps
!     complex, dimension (ntheta0, naky) :: tmp
!     integer, dimension (4) :: start0, count0, start4, count4
!     integer, dimension (3) :: start, count, starth, counth
!     integer, dimension (2) :: startx, countx, starty, county, starts, counts
!     integer :: status, it, ik

!     real, dimension (2, 2*ntgrid+1, ntheta0, naky) :: ri3
!     integer, dimension (5) :: start5, count5

!     start5(1) = 1
!     start5(2) = 1
!     start5(3) = 1
!     start5(4) = 1
!     start5(5) = nout
  
!     count5(1) = 2
!     count5(2) = 2*ntgrid+1
!     count5(3) = ntheta0
!     count5(4) = naky
!     count5(5) = 1
  
!     !End EGH



!     status = nf90_put_var (ncid, time_id, time, start=(/ nout /))

!     start(1) = 1
!     start(2) = 1
!     start(3) = nout

!     count(1) = ntheta0
!     count(2) = naky
!     count(3) = 1

!     start0(1) = 1
!     start0(2) = 1
!     start0(3) = 1
!     start0(4) = nout

!     count0(1) = 2
!     count0(2) = ntheta0
!     count0(3) = naky
!     count0(4) = 1

!     start4(1) = 1
!     start4(2) = 1
!     start4(3) = 1
!     start4(4) = nout

!     count4(1) = ntheta0
!     count4(2) = naky
!     count4(3) = nspec
!     count4(4) = 1

!     starty(1) = 1
!     starty(2) = nout

!     county(1) = naky
!     county(2) = 1

!     startx(1) = 1
!     startx(2) = nout

!     countx(1) = ntheta0
!     countx(2) = 1

!     starts(1) = 1
!     starts(2) = nout

!     counts(1) = nspec
!     counts(2) = 1

!     starth(1) = 1 
!     starth(2) = 1
!     starth(3) = nout

!     counth(1) = nspec
!     counth(2) = 7
!     counth(3) = 1
    
!     if (fphi > zero) then

! !        if(write_apar_t) then
! !           call c2r (aparold, ri3)
! !           !ri_apar_t(:,:,:,:,1) = ri3(:,:,:,:)
! !           status = nf90_put_var (ncid, apar_t_id, ri3, start=start5, count=count5)
! !           if (status /= nf90_noerr) call netcdf_error (status, ncid, apar_id)
! !        end if
       
! !        if(write_bpar_t) then
! !           call c2r (bparold, ri3)
! !           status = nf90_put_var (ncid, bpar_t_id, ri3, start=start5, count=count5)
! !           if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar_id)
! !        end if

!        if (ntheta0 > 1) then
!           do it = 1, ntheta0
!              field2_by_kx(it) = sum(phi2_by_mode(it,:)*fluxfac(:))
!           end do
!           status = nf90_put_var (ncid, phi2_by_kx_id, field2_by_kx, start=startx, count=countx)
!           if (status /= nf90_noerr) call netcdf_error (status, ncid, phi2_by_kx_id)
!        end if
       
!        if (naky > 1) then
!           do ik = 1, naky
!              field2_by_ky(ik) = sum(phi2_by_mode(:,ik)*fluxfac(ik))
!           end do
!           status = nf90_put_var (ncid, phi2_by_ky_id, field2_by_ky, start=starty, count=county)
!           if (status /= nf90_noerr) call netcdf_error (status, ncid, phi2_by_ky_id)
!        end if
       
!        call c2r (phi0, ri2) 
!        status = nf90_put_var (ncid, phi0_id, ri2, start=start0, count=count0)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, phi0_id)
!        status = nf90_put_var (ncid, phi2_by_mode_id, phi2_by_mode, start=start, count=count)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, phi2_by_mode_id)
!        status = nf90_put_var (ncid, phi2_id, phi2, start=(/nout/))
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, phi2_id)
!     end if
    
! !     if (fapar > zero) then
       
! !        if (ntheta0 > 1) then
! !           do it = 1, ntheta0
! !              field2_by_kx(it) = sum(apar2_by_mode(it,:)*fluxfac(:))
! !           end do
! !           status = nf90_put_var (ncid, apar2_by_kx_id, field2_by_kx, start=startx, count=countx)
! !           if (status /= nf90_noerr) call netcdf_error (status, ncid, apar2_by_kx_id)
! !        end if
       
! !        if (naky > 1) then
! !           do ik = 1, naky
! !              field2_by_ky(ik) = sum(apar2_by_mode(:,ik)*fluxfac(ik))
! !           end do
! !           status = nf90_put_var (ncid, apar2_by_ky_id, field2_by_ky, start=starty, count=county)
! !           if (status /= nf90_noerr) call netcdf_error (status, ncid, apar2_by_ky_id)
! !        end if

! !        call c2r (apar0, ri2) 
! !        status = nf90_put_var (ncid, apar0_id, ri2, start=start0, count=count0)
! !        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar0_id)
! !        status = nf90_put_var (ncid, apar2_by_mode_id, apar2_by_mode, start=start, count=count)
! !        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar2_by_mode_id)
! !        status = nf90_put_var (ncid, apar2_id, apar2, start=(/nout/))
! !        if (status /= nf90_noerr) call netcdf_error (status, ncid, apar2_id)
! !     end if

! !     if (fbpar > zero) then

! !        if (ntheta0 > 1) then
! !           do it = 1, ntheta0
! !              field2_by_kx(it) = sum(bpar2_by_mode(it,:)*fluxfac(:))
! !           end do
! !           status = nf90_put_var (ncid, bpar2_by_kx_id, field2_by_kx, start=startx, count=countx)
! !           if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar2_by_kx_id)
! !        end if

! !        if (naky > 1) then
! !           do ik = 1, naky
! !              field2_by_ky(ik) = sum(bpar2_by_mode(:,ik)*fluxfac(ik))
! !           end do
! !           status = nf90_put_var (ncid, bpar2_by_ky_id, field2_by_ky, start=starty, count=county)
! !           if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar2_by_ky_id)
! !        end if

! !        call c2r (bpar0, ri2) 
! !        status = nf90_put_var (ncid, bpar0_id, ri2, start=start0, count=count0)
! !        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar0_id)
! !        status = nf90_put_var (ncid, bpar2_by_mode_id, bpar2_by_mode, start=start, count=count)
! !        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar2_by_mode_id)
! !        status = nf90_put_var (ncid, bpar2_id, bpar2, start=(/nout/))
! !        if (status /= nf90_noerr) call netcdf_error (status, ncid, bpar2_id)
! !     end if

!     if (write_omega) then
!        do it = 1, ntheta0
!           tmp(it, :) = omega(it, :) * woutunits
!        end do
       
!        call c2r (tmp, ri2)
!        status = nf90_put_var (ncid, omega_id, ri2, start=start0, count=count0)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, omega_id)
       
!        do it = 1, ntheta0
!           tmp(it, :) = omegaavg(it, :) * woutunits
!        end do

!        call c2r (tmp, ri2)
!        status = nf90_put_var (ncid, omegaavg_id, ri2, start=start0, count=count0)
!        if (status /= nf90_noerr) call netcdf_error (status, ncid, omegaavg_id)
!     end if

!     if (mod(nout, 10) == 0) then
!        status = nf90_sync (ncid)
!        if (status /= nf90_noerr) call netcdf_error (status)
!     end if
! # endif
!   end subroutine nc_loop

!   ! RN> output full not guiding center moments 
!   subroutine nc_loop_fullmom (nout, time, &
!        & ntot0, density0, upar0, tpar0, tperp0)
!     use kt_grids, only: naky, nakx => ntheta0
!     use species, only: nspec
!     use convert, only: c2r
! # ifdef NETCDF
!     use netcdf, only: nf90_put_var, nf90_sync
!     use netcdf_utils, only: netcdf_error
! # endif
!     implicit none

!     integer, intent (in) :: nout
!     real, intent (in) :: time
!     complex, intent(in) :: ntot0(:,:,:), density0(:,:,:), upar0(:,:,:)
!     complex, intent(in) :: tpar0(:,:,:), tperp0(:,:,:)
! # ifdef NETCDF
!     real, dimension (2, nakx, naky, nspec) :: ri3
!     integer, dimension (5) :: startmom, countmom
!     integer :: status, it, ik, is

!     status = nf90_put_var (ncid, time_id, time, start=(/nout/))
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, time_id)

!     startmom(1) = 1;    countmom(1)=2
!     startmom(2) = 1;    countmom(2)=nakx
!     startmom(3) = 1;    countmom(3)=naky
!     startmom(4) = 1;    countmom(4)=nspec
!     startmom(5) = nout; countmom(5)=1

!     call c2r (ntot0, ri3)
!     status = nf90_put_var (ncid, ntot0_id, ri3, start=startmom, count=countmom)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, ntot0_id)
!     call c2r (density0, ri3)
!     status = nf90_put_var (ncid, density0_id, ri3, start=startmom, count=countmom)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, density0_id)
!     call c2r (upar0, ri3)
!     status = nf90_put_var (ncid, upar0_id, ri3, start=startmom, count=countmom)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, upar0_id)
!     call c2r (tpar0, ri3)
!     status = nf90_put_var (ncid, tpar0_id, ri3, start=startmom, count=countmom)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, tpar0_id)
!     call c2r (tperp0, ri3)
!     status = nf90_put_var (ncid, tperp0_id, ri3, start=startmom, count=countmom)
!     if (status /= nf90_noerr) call netcdf_error (status, ncid, tperp0_id)

!     if (mod(nout, 10) == 0) then
!        status = nf90_sync (ncid)
!        if (status /= nf90_noerr) call netcdf_error (status)
!     end if
! # endif
!   end subroutine nc_loop_fullmom

  subroutine nc_species

    use run_parameters, only: beta
    use species, only: spec
# ifdef NETCDF
    use netcdf, only: nf90_put_var

    integer :: status

    status = nf90_put_var (ncid, charge_id, spec%z)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, charge_id)
    status = nf90_put_var (ncid, mass_id, spec%mass)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, mass_id)
    status = nf90_put_var (ncid, dens_id, spec%dens)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, dens_id)
    status = nf90_put_var (ncid, temp_id, spec%temp)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, temp_id)
    status = nf90_put_var (ncid, tprim_id, spec%tprim)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, tprim_id)
    status = nf90_put_var (ncid, fprim_id, spec%fprim)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, fprim_id)
    status = nf90_put_var (ncid, vnewk_id, spec%vnewk)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, vnewk_id)
    status = nf90_put_var (ncid, spec_type_id, spec%type)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, spec_type_id)

# endif
  end subroutine nc_species

  subroutine nc_geo

    use theta_grid, only: bmag, gradpar, gbdrift, gbdrift0, &
         cvdrift, cvdrift0, gds2, gds21, gds22, grho, jacob, &
         shat, drhodpsi, eps, cdrift, cdrift0, qval
    use run_parameters, only: beta
# ifdef NETCDF
    use netcdf, only: nf90_put_var

    implicit none

    integer :: status

    status = nf90_put_var (ncid, bmag_id, bmag)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, bmag_id)
    status = nf90_put_var (ncid, gradpar_id, gradpar)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, gradpar_id)
    status = nf90_put_var (ncid, gbdrift_id, gbdrift)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, gbdrift_id)
    status = nf90_put_var (ncid, gbdrift0_id, gbdrift0)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, gbdrift0_id)
    status = nf90_put_var (ncid, cvdrift_id, cvdrift)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, cvdrift_id)
    status = nf90_put_var (ncid, cvdrift0_id, cvdrift0)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, cvdrift0_id)
    status = nf90_put_var (ncid, cdrift_id, cdrift)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, cdrift_id)
    status = nf90_put_var (ncid, cdrift0_id, cdrift0)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, cdrift0_id)
    status = nf90_put_var (ncid, gds2_id, gds2)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, gds2_id)
    status = nf90_put_var (ncid, gds21_id, gds21)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, gds21_id)
    status = nf90_put_var (ncid, gds22_id, gds22)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, gds22_id)
    status = nf90_put_var (ncid, grho_id, grho)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, grho_id)
    status = nf90_put_var (ncid, jacob_id, jacob)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, jacob_id)

    status = nf90_put_var (ncid, beta_id, beta)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, beta_id)
    status = nf90_put_var (ncid, q_id, qval)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, q_id)
    status = nf90_put_var (ncid, shat_id, shat)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, shat_id)
    status = nf90_put_var (ncid, eps_id, eps)
    if (status /= nf90_noerr) call netcdf_error (status, ncid, eps_id)
    status = nf90_put_var (ncid, drhodpsi_id, drhodpsi)   
    if (status /= nf90_noerr) call netcdf_error (status, ncid, drhodpsi_id)
# endif
  end subroutine nc_geo

end module stella_io
