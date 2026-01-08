""" ============================================================================================================================================================================= """
""" ----------------- Plot mode frequency and growth rate against the parameter nperiod (controls number of poloidal turns in the field aligned coordinate, z). ----------------- """
""" ------------------------------------------------- We should find that the complex frequency converges for increasing nperiod. ----------------------------------------------- """ 
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

nperiod_vals = [1, 2, 3, 4, 5]

final_omegas = []
final_gammas = []

for nperiod in nperiod_vals:
    output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_beta_0.005_nperiod_scan_no_neoclassics/nperiod_{nperiod}"
    omega_f, gamma_f = get_final_omega_gamma(output_dir)
    if omega_f is not None:
        final_omegas.append(omega_f)
        final_gammas.append(gamma_f)
    else:
        final_omegas.append(np.nan)
        final_gammas.append(np.nan)

# Repeat this for the NEO_stella results.

final_omegas_neo = []
final_gammas_neo = []

for nperiod in nperiod_vals:
    output_dir = f"/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_beta_0.005_nperiod_scan_neoclassics/nperiod_{nperiod}"
    omega_f, gamma_f = get_final_omega_gamma(output_dir)
    if omega_f is not None:
        final_omegas_neo.append(omega_f)
        final_gammas_neo.append(gamma_f)
    else:
        final_omegas_neo.append(np.nan)
        final_gammas_neo.append(np.nan)


# Plot omega and gamma vs nperiod.

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

label_fontsize = 14
tick_fontsize = 12

ax1.plot(nperiod_vals, final_omegas, 'o-', color='dodgerblue', label="conventional")
ax1.plot(nperiod_vals, final_omegas_neo, 'o-', color='crimson', label="neoclassical")
ax1.set_xlabel('nperiod', fontsize=label_fontsize)
ax1.set_ylabel(r'$\omega \; , \; a/v_{th,i}$', fontsize=label_fontsize)
ax1.grid(True, which='major', linestyle='--', alpha=0.7)
ax1.grid(False, which='minor')
ax2.legend(fontsize=tick_fontsize) 

ax2.plot(nperiod_vals, final_gammas, 'o-', color='dodgerblue', label="conventional")
ax2.plot(nperiod_vals, final_gammas_neo, 'o-', color='crimson', label="neoclassical")
ax2.set_xlabel(r'nperiod', fontsize=label_fontsize)
ax2.set_ylabel(r'$\gamma \; , \; a/v_{th,i}$', fontsize=label_fontsize)
ax2.grid(True, which='major', linestyle='--', alpha=0.7)
ax2.grid(False, which='minor')
ax2.legend(fontsize=tick_fontsize) 

fig.suptitle(r'CBC, $\beta = 0.5\%$, $a k_y = 0.3$, $\rho_\star$ = 5e-3, Large Grids', fontsize=label_fontsize + 2)

output_filename = "/users/rjs659/NEO_stella/NEO_stella_CBC_aky_0.3_beta_0.005_nperiod_scan_neoclassics/omega_gamma_vs_nperiod.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"✅ Figure saved as '{output_filename}'")

