# Diagnostics

Contains the diagnostic scripts which analyse the distribution function
and electrostatic potential. Moreover, it contains all writing scripts.
Therefore, the script to write the stella data to a netcdf file to allow
the user to restart the simulation is also in this directory.


## Input file

The diagnostics are writen to text files after every <nwrite> time steps. Moreover, 
they are written to the NetCDF file at every <nwrite>*<nc_mult> time steps. At
every <nsave> time steps, the distribution function is saved to a NetCDF file, 
which allows the user to restart a simulation in the future, if <save_for_restart>
if toggled to True. The frequency and growth rates are calculated for linear 
simulations at every time step, and are also averaged over <navg> time steps.

&diagnostics
    nwrite = 50.0
    navg = 50.0
    nsave = -1.0
    nc_mult = 1.0
    save_for_restart = .false.
    write_all = .false.
    write_all_time_traces = .true.
    write_all_spectra_kxkyz = .false.
    write_all_spectra_kxky = .false.
    write_all_velocity_space = .false.
    write_all_potential = .false.
    write_all_omega = .false.
    write_all_distribution = .false.
    write_all_fluxes = .false.
    write_all_moments = .false.
/
&diagnostics_potential
    write_all_potential_time_traces = .false.
    write_all_potential_spectra = .false.
    write_phi2_vs_time = .true.
    write_apar2_vs_time = .true.
    write_bpar2_vs_time = .true.
    write_phi_vs_kxkyz = .false.
    write_apar_vs_kxkyz = .false.
    write_bpar_vs_kxkyz = .false.
    write_phi2_vs_kxky = .false.
    write_apar2_vs_kxky = .false.
    write_bpar2_vs_kxky = .false.
/
&diagnostics_omega
    write_omega_vs_kxky = .not. include_nonlinear
    write_omega_avg_vs_kxky = .not. include_nonlinear
/
&diagnostics_distribution
    write_g2_vs_vpamus = .false.
    write_g2_vs_zvpas = .false.
    write_g2_vs_zmus = .false.
    write_g2_vs_kxkyzs = .false.
    write_g2_vs_zvpamus = .false.
    write_distribution_g = .true.
    write_distribution_h = .false.
    write_distribution_f = .false.
/
&diagnostics_fluxes
    flux_norm = .true.
    write_fluxes_vs_time = .true.
    write_radial_fluxes = radial_variation
    write_fluxes_kxkyz = .false.
    write_fluxes_kxky = .false.
/
&diagnostics_moments
    write_moments = .false.
    write_radial_moments = .false.
/

## Time traces

By default, the following time traces are always calculated:
  - |phi|^2(t)
  - |apar|^2(t)
  - |bpar|^2(t)

Moreover, for linear simulations the frequency and growth rate are calculated by default:
  - omega(t)
  - gamma(t)
  
## Potential

The following fields are evolved within stella:
  - The perturbed electrostatic potential (phi)
  - The perturbed parallel magnetic potential (apar)
  - The perturbed parallel magnetic field (bpar)

If <write_all_potential_time_traces> = .true., or <write_phi2_vs_time>, <write_apar2_vs_time>
and <write_bpar2_vs_time> are True, then the following time traces are calculated:
  - |phi|^2(t)   -->  phi2 in the NetCDF file, and written to *.out
  - |apar|^2(t)  -->  apar2 in the NetCDF file, and written to *.out
  - |bpar|^2(t)  -->  bpar2 in the NetCDF file, and written to *.out
  
If <write_phi2_vs_kxky>, <write_apar2_vs_kxky> and <write_bpar2_vs_kxky> are
True, then the following quantities are calculated:
  - |phi|^2(t, kx, ky)   -->  phi2_vs_kxky in the NetCDF file
  - |apar|^2(t, kx, ky)  -->  apar2_vs_kxky in the NetCDF file
  - |bpar|^2(t, kx, ky)  -->  bpar2_vs_kxky in the NetCDF file
  
If <write_phi_vs_kxkyz>, <write_apar_vs_kxkyz> and <write_bpar_vs_kxkyz> are
True, then the following complex quantities are calculated:
  - phi(t, kx, ky, z, ri)   -->  phi_vs_t in the NetCDF file
  - apar(t, kx, ky, z, ri)  -->  apar_vs_t in the NetCDF file
  - bpar(t, kx, ky, z, ri)  -->  bpar_vs_t in the NetCDF file

## Distribution

The following distribution functions are evolved within stella:
  - The perturbed distribution function (f)
  - The perturbed gyroaveraged distribution function (g)
  - The non-adiabatic part of the perturbed distribution function (h)

If <write_g2_vs_vpamus> = .true. the following quantities are calculated:
  - |g|^2(t, mu, vpa, s)           -->  g2_vs_vpamus in the NetCDF file
  - |g_nozonal|^2(t, mu, vpa, s)   -->  g2nozonal_vs_vpamus in the NetCDF file

If <write_g2_vs_zvpas> = .true. the following quantities are calculated:
  - |g|^2(t, z, vpa, s)            -->  g2_vs_zvpas in the NetCDF file
  - |g_nozonal|^2(t, z, vpa, s)    -->  g2nozonal_vs_zvpas in the NetCDF file

If <write_g2_vs_zmus> = .true. the following quantities are calculated:
  - |g|^2(t, z, mu, s)             -->  g2_vs_zmus in the NetCDF file
  - |g_nozonal|^2(t, z, mu, s)     -->  g2nozonal_vs_zmus in the NetCDF file

If <write_g2_vs_kxkyzs> = .true. the following quantities are calculated:
  - |g|^2(t, kx, ky, z, s)          -->  g2_vs_zkykxs in the NetCDF file
  - |g_nozonal|^2(t, kx, ky, z, s)  -->  g2nozonal_vs_zkykxs in the NetCDF file

If <write_g2_vs_zvpamus> = .true. the following quantities are calculated:
  - |g|^2(t, z, vpa, mu, s)          -->  g2_vs_zvpamus in the NetCDF file
  - |g_nozonal|^2(t, z, vpa, mu, s)  -->  g2nozonal_vs_zvpamus in the NetCDF file

If <write_distribution_h> = .true., and depending on the previous flags, we calculate:
  - |h|^2(t, mu, vpa, s)             -->  h2_vs_vpamus in the NetCDF file
  - |h|^2(t, z, vpa, s)              -->  h2_vs_zvpas in the NetCDF file
  - |h|^2(t, z, mu, s)               -->  h2_vs_zmus in the NetCDF file
  - |h|^2(t, kx, ky, z, s)           -->  h2_vs_zkykxs in the NetCDF file
  - |h|^2(t, z, vpa, mu, s)          -->  h2_vs_zvpamus in the NetCDF file
  - |h_nozonal|^2(t, mu, vpa, s)     -->  h2nozonal_vs_vpamus in the NetCDF file
  - |h_nozonal|^2(t, z, vpa, s)      -->  h2nozonal_vs_zvpas in the NetCDF file
  - |h_nozonal|^2(t, z, mu, s)       -->  h2nozonal_vs_zmus in the NetCDF file
  - |h_nozonal|^2(t, kx, ky, z, s)   -->  h2nozonal_vs_zkykxs in the NetCDF file
  - |h_nozonal|^2(t, z, vpa, mu, s)  -->  h2nozonal_vs_zvpamus in the NetCDF file

If <write_distribution_f> = .true., and depending on the previous flags, we calculate:
  - |f|^2(t, mu, vpa, s)             -->  f2_vs_vpamus in the NetCDF file
  - |f|^2(t, z, vpa, s)              -->  f2_vs_zvpas in the NetCDF file
  - |f|^2(t, z, mu, s)               -->  f2_vs_zmus in the NetCDF file
  - |f|^2(t, kx, ky, z, s)           -->  f2_vs_zkykxs in the NetCDF file
  - |f|^2(t, z, vpa, mu, s)          -->  f2_vs_zvpamus in the NetCDF file
  - |f_nozonal|^2(t, mu, vpa, s)     -->  f2nozonal_vs_vpamus in the NetCDF file
  - |f_nozonal|^2(t, z, vpa, s)      -->  f2nozonal_vs_zvpas in the NetCDF file
  - |f_nozonal|^2(t, z, mu, s)       -->  f2nozonal_vs_zmus in the NetCDF file
  - |f_nozonal|^2(t, kx, ky, z, s)   -->  f2nozonal_vs_zkykxs in the NetCDF file
  - |f_nozonal|^2(t, z, vpa, mu, s)  -->  f2nozonal_vs_zvpamus in the NetCDF file

## Omega

The complex frequency Omega is calculated within stella for linear simulations:
  - <phi> = exp(-i*<0mega>*t) = exp(-i*(<omega>+i<gamma>)*t)
  - <Omega> = log(dphi)/(-i*dt) = i*log(dphi)/dt
  - <omega> = real(<Omega>)
  - <gamma> = imag(<Omega>)
  
The frequency and growth rate are calculated automatically for linear simulations:
  - omega(t, kx, ky)      -->  Written to the *.omega text file
  - gamma(t, kx, ky)      -->  Written to the *.omega text file
  - Omega(t, kx, ky, ri)  -->  Written as "omega" to the NetCDF file

## Fluxes

The following fluxes can be calculated within stella:
  - Particle flux (pflux)
  - Heat flux (qflux)
  - Momentum flux (vflux)

If <write_fluxes_vs_time> = .true. the following time traces are calculated:
  - pflux(t, s)  -->  pflux_vs_s in the NetCDF file, and written to *.fluxes
  - qflux(t, s)  -->  qflux_vs_s in the NetCDF file, and written to *.fluxes
  - vflux(t, s)  -->  vflux_vs_s in the NetCDF file, and written to *.fluxes

If <write_fluxes_kxky> = .true. the flux spectra are calculated, note that these
do not represent the Fourier components of the fluxes, but rather the contribution 
of each (kx,ky) mode to the flux. To obtain the total flux, simply sum over all 
contributions, since the reality condition has already been taking into account.
  - pflux(t, kx, ky, s)  -->  pflux_vs_kxkys in the NetCDF file
  - qflux(t, kx, ky, s)  -->  qflux_vs_kxkys in the NetCDF file
  - vflux(t, kx, ky, s)  -->  vflux_vs_kxkys in the NetCDF file
  
If <write_fluxes_kxkyz> = .true. the following quantities are calculated:
  - pflux(t, kx, ky, z, s)  -->  pflux_vs_kxkyzs in the NetCDF file
  - qflux(t, kx, ky, z, s)  -->  qflux_vs_kxkyzs in the NetCDF file
  - vflux(t, kx, ky, z, s)  -->  vflux_vs_kxkyzs in the NetCDF file
  
If <write_radial_fluxes> = .true. the following quantities are calculated:
  - pflux(t, kx, s)  -->  pflux_x in the NetCDF file
  - qflux(t, kx, s)  -->  qflux_x in the NetCDF file
  - vflux(t, kx, s)  -->  vflux_x in the NetCDF file

## Moments

The following moments can be calculated within stella:
  - Density (dens)
  - Temperature (temp)
  - Parallel velocity (upar)

If <write_moments> = .true. the following complex quantities are calculated:
  - dens(t, kx, ky, z, s, ri)  -->  density in the NetCDF file
  - temp(t, kx, ky, z, s, ri)  -->  temperature in the NetCDF file
  - upar(t, kx, ky, z, s, ri)  -->  upar in the NetCDF file

If <write_radial_moments> = .true. the following complex quantities are calculated:
  - dens(t, kx, s, ri)  -->  dens_x in the NetCDF file
  - temp(t, kx, s, ri)  -->  upar_x in the NetCDF file
  - upar(t, kx, s, ri)  -->  temp_x in the NetCDF file


