""" ============================================================================================================================================================================= """
""" -------------------------------------------------- Plot mode frequency and growth rate against the gyrokinetic parameter rhostar. --------------------------------------------"""
""" ------------- We should find that the complex frequency does not change for conventional gyrokinetic runs, while it does change in the presence of neoclassics. --------------""" 
""" =============================================================================================================================================================================="""

import numpy as np
import matplotlib.pyplot as plt
import os


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

# Scan over rhostar values.

rhostar_vals = np.logspace(-7, -1, 20)[:-2]

final_omegas = []
final_gammas = []

for rhostar in rhostar_vals:
    output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_rhostar_scan_no_neoclassics_kinetic_e/rhostar_{rhostar}"
    omega_f, gamma_f = get_final_omega_gamma(output_dir)
    if omega_f is not None:
        final_omegas.append(omega_f)
        final_gammas.append(gamma_f)
    else:
        final_omegas.append(np.nan)
        final_gammas.append(np.nan)

# Repeat this for the NEO_stella results.
# ------------------------------------------------------------------------------------------------------------------------------------------------------- #

final_omegas_neo = []
final_gammas_neo = []

for rhostar in rhostar_vals:
    output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_rhostar_scan_neoclassics_kinetic_e/rhostar_{rhostar}"
    omega_f, gamma_f = get_final_omega_gamma(output_dir)
    if omega_f is not None:
        final_omegas_neo.append(omega_f)
        final_gammas_neo.append(gamma_f)
    else:
        final_omegas_neo.append(np.nan)
        final_gammas_neo.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------------------------- #

# final_omegas_neo_dchidz = []
# final_gammas_neo_dchidz = []

# for rhostar in rhostar_vals:
    # output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_rhostar_scan_neoclassics_adiabatic_e_dchidz_terms_only/rhostar_{rhostar}"
    # omega_f, gamma_f = get_final_omega_gamma(output_dir)
    # if omega_f is not None:
        # final_omegas_neo_dchidz.append(omega_f)
        # final_gammas_neo_dchidz.append(gamma_f)
    # else:
        # final_omegas_neo_dchidz.append(np.nan)
        # final_gammas_neo_dchidz.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------------------------- #

# final_omegas_neo_curv = []
# final_gammas_neo_curv = []

# for rhostar in rhostar_vals:
    # output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_rhostar_scan_neoclassics_adiabatic_e_curv_drift_terms_only/rhostar_{rhostar}"
    # omega_f, gamma_f = get_final_omega_gamma(output_dir)
    # if omega_f is not None:
        # final_omegas_neo_curv.append(omega_f)
        # final_gammas_neo_curv.append(gamma_f)
    # else:
        # final_omegas_neo_curv.append(np.nan)
        # final_gammas_neo_curv.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------------------------- #

# final_omegas_neo_mag = []
# final_gammas_neo_mag = []

# for rhostar in rhostar_vals:
    # output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_rhostar_scan_neoclassics_adiabatic_e_mag_drift_terms_only/rhostar_{rhostar}"
    # omega_f, gamma_f = get_final_omega_gamma(output_dir)
    # if omega_f is not None:
        # final_omegas_neo_mag.append(omega_f)
        # final_gammas_neo_mag.append(gamma_f)
    # else:
        # final_omegas_neo_mag.append(np.nan)
        # final_gammas_neo_mag.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------------------------- #

# final_omegas_neo_drive = []
# final_gammas_neo_drive = []

# for rhostar in rhostar_vals:
    # output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_rhostar_scan_neoclassics_adiabatic_e_drive_terms_only/rhostar_{rhostar}"
    # omega_f, gamma_f = get_final_omega_gamma(output_dir)
    # if omega_f is not None:
       #  final_omegas_neo_drive.append(omega_f)
       #  final_gammas_neo_drive.append(gamma_f)
    # else:
        # final_omegas_neo_drive.append(np.nan)
        # final_gammas_neo_drive.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------------------------- #

# final_omegas_neo_all = []
# final_gammas_neo_all = []

# for rhostar in rhostar_vals:
    # output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_rhostar_scan_neoclassics_adiabatic_e/rhostar_{rhostar}"
    # omega_f, gamma_f = get_final_omega_gamma(output_dir)
    # if omega_f is not None:
        # final_omegas_neo_all.append(omega_f)
        # final_gammas_neo_all.append(gamma_f)
    # else:
        # final_omegas_neo_all.append(np.nan)
        # final_gammas_neo_all.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------------------------- #
# --- Plot omega and gamma vs aky (side-by-side with log x-axis) ---.

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

label_fontsize = 14
tick_fontsize = 12

ax1.plot(rhostar_vals, final_omegas, 'o-', color='dodgerblue', label="conventional")
ax1.plot(rhostar_vals, final_omegas_neo, 'o-', color='crimson', label=r"neoclassical")
# ax1.plot(rhostar_vals, final_omegas_neo_dchidz, 'o-', color='forestgreen', label=r"neoclassical, $d\tilde{\chi}'/dz$ terms only")
# ax1.plot(rhostar_vals, final_omegas_neo_curv, 'o-', color='orange', label=r"neoclassical, $\tilde{\omega}_\kappa$ terms only")
# ax1.plot(rhostar_vals, final_omegas_neo_mag, 'o-', color='mediumslateblue', label=r"neoclassical, $\tilde{\omega}_d$ terms only")
# ax1.plot(rhostar_vals, final_omegas_neo_drive, 'o-', color='dimgrey', label=r"neoclassical, $\tilde{\omega}_{*,1}$ only")
# ax1.plot(rhostar_vals, final_omegas_neo_all, 'x-', color='deeppink', label="neoclassical, all terms")
ax1.set_xscale('log')
ax1.set_xlabel(r'$\rho_\star$', fontsize=label_fontsize)
ax1.set_ylabel(r'$\omega \; , \; a/v_{th,i}$', fontsize=label_fontsize)
ax1.grid(True, which='major', linestyle='--', alpha=0.7)
ax1.grid(False, which='minor')

ax2.plot(rhostar_vals, final_gammas, 'o-', color='dodgerblue', label="conventional")
ax2.plot(rhostar_vals, final_gammas_neo, 'o-', color='crimson', label=r"neoclassical")
# ax2.plot(rhostar_vals, final_gammas_neo_dchidz, 'o-', color='forestgreen', label=r"neoclassical, $d\tilde{\chi}'/dz$ terms only")
# ax2.plot(rhostar_vals, final_gammas_neo_curv, 'o-', color='orange', label=r"neoclassical, $\tilde{\omega}_\kappa$ terms only")
# ax2.plot(rhostar_vals, final_gammas_neo_mag, 'o-', color='mediumslateblue', label=r"neoclassical, $\tilde{\omega}_d$ terms only")
# ax2.plot(rhostar_vals, final_gammas_neo_drive, 'o-', color='dimgrey', label=r"neoclassical, $\tilde{\omega}_{*,1}$ only")
# ax2.plot(rhostar_vals, final_gammas_neo_all, 'x-', color='deeppink', label="neoclassical, all terms")
ax2.set_xscale('log')
ax2.set_xlabel(r'$\rho_\star$', fontsize=label_fontsize)
ax2.set_ylabel(r'$\gamma \; , \; a/v_{th,i}$', fontsize=label_fontsize)
ax2.grid(True, which='major', linestyle='--', alpha=0.7)
ax2.grid(False, which='minor')
ax2.legend(loc='lower center', bbox_to_anchor=(-0.1, 1.05), ncol=2, fontsize=11, frameon=True)
plt.subplots_adjust(top=0.85)

fig.suptitle(r'CBC, $\beta = 0.0\%$, $k_y = 0.3$, Kinetic Electrons', fontsize=label_fontsize + 2)

output_filename = "/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_rhostar_scan_neoclassics_kinetic_e/omega_gamma_vs_rhostar.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"✅ Figure saved as '{output_filename}'")

