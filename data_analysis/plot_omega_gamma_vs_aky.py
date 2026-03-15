""" ============================================================================================================================================================================= """
""" -------------------------------------------------------------- Plot mode frequency and growth rate against modenumber. -------------------------------------------------------"""
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

aky_vals = np.linspace(0.1, 1.0, 50)

final_omegas = []
final_gammas = []

for aky in aky_vals:
    output_dir = f"/users/rjs659/NEO_stella/CBC_electromagnetic_beta_0_02_aky_scan_rhostar_0_01_neoclassics_adiabatic_e_explicit/aky_{aky}"
    omega_f, gamma_f = get_final_omega_gamma(output_dir)
    if omega_f is not None:
        final_omegas.append(omega_f)
        final_gammas.append(gamma_f)
    else:
        final_omegas.append(np.nan)
        final_gammas.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------ #

# final_omegas_neo_rhostar_1e_3 = []
# final_gammas_neo_rhostar_1e_3 = []

# for aky in aky_vals:
#     output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_scan_electrostatic_rhostar_1e-3_neoclassics_adiabatic_e/aky_{aky}"
#     omega_f, gamma_f = get_final_omega_gamma(output_dir)
#     if omega_f is not None:
#         final_omegas_neo_rhostar_1e_3.append(omega_f)
#         final_gammas_neo_rhostar_1e_3.append(gamma_f)
#     else:
#         final_omegas_neo_rhostar_1e_3.append(np.nan)
#         final_gammas_neo_rhostar_1e_3.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------ #

final_omegas_neo_rhostar_1e_2 = []
final_gammas_neo_rhostar_1e_2 = []

for aky in aky_vals:
    output_dir = f"/users/rjs659/NEO_stella/CBC_electromagnetic_beta_0_02_aky_scan_rhostar_0_01_no_neoclassics_adiabatic_e_explicit/aky_{aky}"
    omega_f, gamma_f = get_final_omega_gamma(output_dir)
    if omega_f is not None:
        final_omegas_neo_rhostar_1e_2.append(omega_f)
        final_gammas_neo_rhostar_1e_2.append(gamma_f)
    else:
        final_omegas_neo_rhostar_1e_2.append(np.nan)
        final_gammas_neo_rhostar_1e_2.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------ #
# --- Plot omega and gamma vs aky (side-by-side with log x-axis) ---.

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

label_fontsize = 14
tick_fontsize = 12

ax1.plot(aky_vals, final_omegas, 'o-', markerfacecolor='none', color='dodgerblue', label=r"neoclassics, explicit")
# ax1.plot(aky_vals, final_omegas_neo_rhostar_1e_3, 'o-', markerfacecolor='none', color='forestgreen', label=r"neoclassical, $\rho_\star = 10e^{-3}$")
ax1.plot(aky_vals, final_omegas_neo_rhostar_1e_2, 'o-', markerfacecolor='none', color='orange', label=r"no neoclassics, explicit$")
ax1.set_xscale('log')
# ax1.set_yscale('log')
ax1.set_xlabel(r'$k_y\rho_i$', fontsize=label_fontsize)
ax1.set_ylabel(r'$\omega \; , \; a/v_{th,i}$', fontsize=label_fontsize)
ax1.grid(True, which='major', linestyle='--', alpha=0.7)
ax1.grid(False, which='minor')
ax1.legend(fontsize=tick_fontsize) 

ax2.plot(aky_vals, final_gammas, 'o-', markerfacecolor='none', color='dodgerblue', label=r"neoclassics, explicit")
# ax2.plot(aky_vals, final_gammas_neo_rhostar_1e_3, 'o-', markerfacecolor='none', color='forestgreen', label=r"neoclassical, $\rho_\star = 10e^{-3}$")
ax2.plot(aky_vals, final_gammas_neo_rhostar_1e_2, 'o-', markerfacecolor='none', color='orange', label=r"no neoclassics, explicit")
ax2.set_xscale('log')
# ax2.set_yscale('log')
ax2.set_xlabel(r'$k_y\rho_i$', fontsize=label_fontsize)
ax2.set_ylabel(r'$\gamma \; , \; a/v_{th,i}$', fontsize=label_fontsize)
ax2.grid(True, which='major', linestyle='--', alpha=0.7)
ax2.grid(False, which='minor')
ax2.legend(fontsize=tick_fontsize) 

fig.suptitle(r'CBC, $\beta = 2.0\%$, $\rho_\star = 10e^{-2}$, Adiabatic Electrons', fontsize=label_fontsize + 2)

output_filename = "/users/rjs659/NEO_stella/CBC_electromagnetic_beta_0_02_aky_scan_rhostar_0_01_neoclassics_adiabatic_e_explicit/omega_gamma_vs_aky.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"✅ Figure saved as '{output_filename}'")

