""" ============================================================================================================================================================================= """
""" -------------------------------------------------- Plot mode frequency and growth rate against the gyrokinetic parameter rhostar. --------------------------------------------"""
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

fprim_vals = [2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50]

final_omegas = []
final_gammas = []

for fprim in fprim_vals:
    output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_fprim_scan_no_neoclassics_adiabatic_e_rhostar_1e-2/fprim_{fprim}"
    omega_f, gamma_f = get_final_omega_gamma(output_dir)
    if omega_f is not None:
        final_omegas.append(omega_f)
        final_gammas.append(gamma_f)
    else:
        final_omegas.append(np.nan)
        final_gammas.append(np.nan)

# Repeat this for each directory of NEO_stella results.

# final_omegas_neo_rhostar_1e_4 = []
# final_gammas_neo_rhostar_1e_4 = []

# for fprim in fprim_vals:
#     output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_fprim_scan_neoclassics_adiabatic_e_rhostar_1e-4/fprim_{fprim}"
#     omega_f, gamma_f = get_final_omega_gamma(output_dir)
#     if omega_f is not None:
#         final_omegas_neo_rhostar_1e_4.append(omega_f)
#         final_gammas_neo_rhostar_1e_4.append(gamma_f)
#     else:
#         final_omegas_neo_rhostar_1e_4.append(np.nan)
#         final_gammas_neo_rhostar_1e_4.append(np.nan)

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# final_omegas_neo_rhostar_1e_3 = []
# final_gammas_neo_rhostar_1e_3 = []

# for fprim in fprim_vals:
#     output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_fprim_scan_neoclassics_adiabatic_e_rhostar_1e-3/fprim_{fprim}"
#     omega_f, gamma_f = get_final_omega_gamma(output_dir)
#     if omega_f is not None:
#         final_omegas_neo_rhostar_1e_3.append(omega_f)
#         final_gammas_neo_rhostar_1e_3.append(gamma_f)
#     else:
#         final_omegas_neo_rhostar_1e_3.append(np.nan)
#         final_gammas_neo_rhostar_1e_3.append(np.nan)

# -------------------------------------------------------------------------------------------------------------------------------------------- #

# final_omegas_neo_rhostar_1e_2 = []
# final_gammas_neo_rhostar_1e_2 = []

# for fprim in fprim_vals:
#     output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_fprim_scan_neoclassics_adiabatic_e_rhostar_1e-2/fprim_{fprim}"
#     omega_f, gamma_f = get_final_omega_gamma(output_dir)
#     if omega_f is not None:
#         final_omegas_neo_rhostar_1e_2.append(omega_f)
#         final_gammas_neo_rhostar_1e_2.append(gamma_f)
#     else:
#         final_omegas_neo_rhostar_1e_2.append(np.nan)
#         final_gammas_neo_rhostar_1e_2.append(np.nan)

# ---------------------------------------- Plot omega and gamma vs aky (side-by-side with log x-axis) ---------------------------------------- #

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

label_fontsize = 14
tick_fontsize = 12

ax1.plot(fprim_vals, final_omegas, 'o-', color='dodgerblue', label="")
# ax1.plot(fprim_vals, final_omegas_neo_rhostar_1e_4, 'o-', color='crimson', label="")
# ax1.plot(fprim_vals, final_omegas_neo_rhostar_1e_3, 'o-', color='forestgreen', label="")
# ax1.plot(fprim_vals, final_omegas_neo_rhostar_1e_2, 'o-', color='gold', label="")
ax1.set_xlabel(r'$a/|L_n|$', fontsize=label_fontsize)
ax1.set_ylabel(r'$\omega \; , \; a/v_{th,i}$', fontsize=label_fontsize)
ax1.grid(True, which='major', linestyle='--', alpha=0.7)
ax1.grid(False, which='minor')
ax2.legend(fontsize=tick_fontsize) 

ax2.plot(fprim_vals, final_gammas, 'o-', color='dodgerblue', label="")
# ax2.plot(fprim_vals, final_gammas_neo_rhostar_1e_4, 'o-', color='crimson', label="")
# ax1.plot(fprim_vals, final_gammas_neo_rhostar_1e_3, 'o-', color='forestgreen', label="")
# ax1.plot(fprim_vals, final_gammas_neo_rhostar_1e_2, 'o-', color='gold', label="")
ax2.set_xlabel(r'$a/|L_n|$', fontsize=label_fontsize)
ax2.set_ylabel(r'$\gamma \; , \; a/v_{th,i}$', fontsize=label_fontsize)
ax2.grid(True, which='major', linestyle='--', alpha=0.7)
ax2.grid(False, which='minor')
ax2.legend(fontsize=tick_fontsize) 

fig.suptitle(r'CBC, $\beta = 0.0\%$, $k_y = 0.3$, Adiabatic Electrons', fontsize=label_fontsize + 2)

output_filename = "/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_electrostatic_fprim_scan_no_neoclassics_adiabatic_e_rhostar_1e-2/omega_gamma_vs_rhostar.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"✅ Figure saved as '{output_filename}'")

