""" ============================================================================================================================================================================= """
""" -------------------------------------------------------------- Plot mode frequency and growth rate against modenumber. -------------------------------------------------------"""
""" =============================================================================================================================================================================="""

import numpy as np
import matplotlib.pyplot as plt
import os
from ncdf2dict import ncdf2dict

# ------------------------------------------------------------------------------------------------------------------------------------ #

# aky values to scan.
aky_vals = np.logspace(-1, 1, 20)[2:]


# ------------------------------------------------------------------------------------------------------------------------------------ #
# Function to read in stella results. 

def get_final_omega_gamma(output_dir):
    """Read input.omega in the specified directory and return the final omega and gamma."""
    omega_file = os.path.join(output_dir, "CBC.omega")
    if not os.path.exists(omega_file):
        print(f"{omega_file} not found!")
        return None, None

    data = np.loadtxt(omega_file, comments='#')
    omega = data[:, 3]  # Re[om].
    gamma = data[:, 4]  # Im[om].

    final_omega = omega[-1]
    final_gamma = gamma[-1]

    return final_omega, final_gamma


# ------------------------------------------------------------------------------------------------------------------------------------ #
# Function to read gs2 results.

def read_gs2_file(nc_file):
    """Read a GS2 input.out.nc file and return final omega and gamma."""
    data = ncdf2dict(nc_file)
    if data is None or 'omega' not in data:
        print(f"⚠️ {nc_file} not found or missing 'omega' key!")
        return None, None

    omega = np.vstack(data['omega']).reshape(-1)
    final_freq = omega.real[-1]
    final_growth = omega.imag[-1]
    return final_freq, final_growth

def read_gs2_dir(base_dir):
    """Loop over subdirectories, read GS2 results, and normalize."""
    final_freqs, final_growths, found_aky_vals = [], [], []

    for subdir in sorted(os.listdir(base_dir)):
        full_path = os.path.join(base_dir, subdir)
        if not os.path.isdir(full_path):
            continue

        nc_files = [f for f in os.listdir(full_path) if f.endswith('.out.nc')]
        if not nc_files:
            print(f"⚠️ No .out.nc file found in {subdir}")
            continue

        nc_path = os.path.join(full_path, nc_files[0])
        final_freq, final_growth = read_gs2_file(nc_path)
        if final_freq is not None:
            final_freqs.append(final_freq)
            final_growths.append(final_growth)
            # Extract aky value from folder name (assumes 'aky_X.X')
            try:
                aky_val = float(subdir.split("_")[-1])
                found_aky_vals.append(aky_val)
            except ValueError:
                print(f"⚠️ Could not parse aky value from {subdir}")
        else:
            print(f"⚠️ Could not read GS2 results in {subdir}")

    # Convert to arrays and normalize.
    final_freqs = np.array(final_freqs) / 1.0       # 2.77778
    final_growths = np.array(final_growths) / 1.0   # 2.77778
    found_aky_vals = np.array(found_aky_vals)

    return found_aky_vals, final_freqs, final_growths

 
# ------------------------------------------------------------------------------------------------------------------------------------ #
# Read in data.

# final_omegas_exp = []
# final_gammas_exp = []

# for aky in aky_vals:
    # output_dir = f"/users/rjs659/NEO_stella/CBC_electrostatic_aky_scan_rhostar_0_01_no_neoclassics_adiabatic_e_explicit/aky_{aky}"
    # omega_f, gamma_f = get_final_omega_gamma(output_dir)
    # if omega_f is not None:
        # final_omegas_exp.append(omega_f)
        # final_gammas_exp.append(gamma_f)
    # else:
        # final_omegas_exp.append(np.nan)
        # final_gammas_exp.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------ #

final_omegas_imp = []
final_gammas_imp = []

for aky in aky_vals:
    output_dir = f"/users/rjs659/NEO_stella/CBC_akx_0.1_electrostatic_rhostar_0_01_aky_scans/NEO_stella_CBC_electrostatic_aky_scan_akx_0.1_rhostar_0_01_no_neoclassics_kinetic_e/aky_{aky}"
    omega_f, gamma_f = get_final_omega_gamma(output_dir)
    if omega_f is not None:
        final_omegas_imp.append(omega_f)
        final_gammas_imp.append(gamma_f)
    else:
        final_omegas_imp.append(np.nan)
        final_gammas_imp.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------ #

# final_omegas_exp_neo = []
# final_gammas_exp_neo = []

# for aky in aky_vals:
    # output_dir = f"/users/rjs659/NEO_stella/CBC_electrostatic_aky_scan_rhostar_0_01_neoclassics_adiabatic_e_explicit/aky_{aky}"
    # omega_f, gamma_f = get_final_omega_gamma(output_dir)
    # if omega_f is not None:
        # final_omegas_exp_neo.append(omega_f)
        # final_gammas_exp_neo.append(gamma_f)
    # else:
        # final_omegas_exp_neo.append(np.nan)
        # final_gammas_exp_neo.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------ #

final_omegas_imp_neo = []
final_gammas_imp_neo = []

for aky in aky_vals:
    output_dir = f"/users/rjs659/NEO_stella/CBC_akx_0.1_electrostatic_rhostar_0_01_aky_scans/NEO_stella_CBC_electrostatic_aky_scan_akx_0.1_rhostar_0_01_neoclassics_kinetic_e/aky_{aky}"
    omega_f, gamma_f = get_final_omega_gamma(output_dir)
    if omega_f is not None:
        final_omegas_imp_neo.append(omega_f)
        final_gammas_imp_neo.append(gamma_f)
    else:
        final_omegas_imp_neo.append(np.nan)
        final_gammas_imp_neo.append(np.nan)


# ------------------------------------------------------------------------------------------------------------------------------------ #
# Read in gs2 data.

aky_gs2, freq_gs2, growth_gs2         = read_gs2_dir("/users/rjs659/NEO_stella/CBC_akx_0.1_electrostatic_rhostar_0_01_aky_scans/NEO_gs2_CBC_aky_scan_akx_0.1_electrostatic_rhostar_0_01_no_neoclassics_kinetic_e")
aky_gs2, freq_neo_gs2, growth_neo_gs2 = read_gs2_dir("/users/rjs659/NEO_stella/CBC_akx_0.1_electrostatic_rhostar_0_01_aky_scans/NEO_gs2_CBC_aky_scan_akx_0.1_electrostatic_rhostar_0_01_neoclassics_kinetic_e")

aky_gs2 = aky_gs2[2:]
freq_gs2 = freq_gs2[2:]
growth_gs2 = growth_gs2[2:]
freq_neo_gs2 = freq_neo_gs2[2:]
growth_neo_gs2 = growth_neo_gs2[2:]

# ------------------------------------------------------------------------------------------------------------------------------------ #
# Plot omega and gamma vs aky.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

# Define the stride (how often to place a marker). 
stride = 4 

label_fontsize = 14
tick_fontsize = 12

# Stella results.
# ax1.plot(aky_vals, final_omegas_exp, 'o-', lw=2, markerfacecolor='dodgerblue', color='dodgerblue', label=r"stella, explicit", markevery=(1, stride))
ax1.plot(aky_vals, final_omegas_imp, 'o-', lw=1.2, markerfacecolor='dodgerblue', color='dodgerblue', label=r"stella, implicit", markevery=(1, stride))
# ax1.plot(aky_vals, final_omegas_exp_neo, 'o-', lw=1, markerfacecolor='forestgreen', color='forestgreen', label=r"NEO-stella, explicit", markevery=(3, stride))
ax1.plot(aky_vals, final_omegas_imp_neo, 'o-', lw=1, markerfacecolor='forestgreen', color='forestgreen', label=r"NEO-stella, implicit", markevery=(3, stride))

# gs2 results.
ax1.plot(aky_vals, freq_gs2, 'o-', lw=1.5, markerfacecolor='crimson', color='crimson', label=r"gs2", markevery=(2, stride))
ax1.plot(aky_vals, freq_neo_gs2, 'o-', lw=0.5, markerfacecolor='gold', color='gold', label=r"NEO-gs2", markevery=(4, stride))

# Formatting. 
ax1.set_xscale('log')
ax1.set_xlabel(r'$k_y\rho_i$', fontsize=label_fontsize)
ax1.set_ylabel(r'$\omega \; , \; a/v_{th,i}$', fontsize=label_fontsize)
ax1.grid(True, which='major', linestyle='--', alpha=0.7)
ax1.grid(False, which='minor')
ax1.legend(fontsize=tick_fontsize) 

# Stella results. 
# ax2.plot(aky_vals, final_gammas_exp, 'o-', lw=2, markerfacecolor='dodgerblue', color='dodgerblue', label=r"stella, explicit", markevery=(1, stride))
ax2.plot(aky_vals, final_gammas_imp, 'o-', lw=2, markerfacecolor='dodgerblue', color='dodgerblue', label=r"stella, implicit", markevery=(1, stride)) 
# ax2.plot(aky_vals, final_gammas_exp_neo, 'o-', lw=1, markerfacecolor='forestgreen', color='forestgreen', label=r"NEO-stella, explicit", markevery=(3, stride))
ax2.plot(aky_vals, final_gammas_imp_neo, 'o-', lw=1, markerfacecolor='forestgreen', color='forestgreen', label=r"NEO-stella, implicit", markevery=(3, stride))

# gs2 results. 
ax2.plot(aky_vals, growth_gs2, 'o-', lw=1.5, markerfacecolor='crimson', color='crimson', label=r"gs2", markevery=(2, stride))
ax2.plot(aky_vals, growth_neo_gs2, 'o-', lw=0.5, markerfacecolor='gold', color='gold', label=r"NEO-gs2", markevery=(4, stride))

# Formatting.
ax2.set_xscale('log')
ax2.set_xlabel(r'$k_y\rho_i$', fontsize=label_fontsize)
ax2.set_ylabel(r'$\gamma \; , \; a/v_{th,i}$', fontsize=label_fontsize)
ax2.grid(True, which='major', linestyle='--', alpha=0.7)
ax2.grid(False, which='minor')
ax2.legend(fontsize=tick_fontsize) 

fig.suptitle(r'CBC, $\beta = 0.0\%$, $k_x = 0.1$, $\rho_\star = 1e^{-2}$, Kinetic Electrons', fontsize=label_fontsize + 2)

output_filename = "/users/rjs659/NEO_stella/CBC_akx_0.1_electrostatic_rhostar_0_01_aky_scans/NEO_gs2_CBC_aky_scan_akx_0.1_electrostatic_rhostar_0_01_neoclassics_kinetic_e/omega_gamma_vs_aky.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"✅ Figure saved as '{output_filename}'")

