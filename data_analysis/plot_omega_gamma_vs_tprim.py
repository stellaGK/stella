""" ============================================================================================================================================================================= """
""" -------------------------------------------------- Plot mode frequency and growth rate against the gyrokinetic parameter rhostar. --------------------------------------------"""
""" =============================================================================================================================================================================="""

import numpy as np
import matplotlib.pyplot as plt
import os
from ncdf2dict import ncdf2dict

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

tprims = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]

final_omegas = []
final_gammas = []

# Read in stella/NEO-stella data. 
# ------------------------------------------------------------------------------------------------------------------------------------------- #

for tprim in tprims:
    output_dir = f"/users/rjs659/NEO_stella/ky_0_3_tprim_scans/NEO_stella_CBC_aky_0.3_electrostatic_tprim_scan_no_neoclassics_kinetic_e_rhostar_1e-2/tprim_{tprim}"
    omega_f, gamma_f = get_final_omega_gamma(output_dir)
    if omega_f is not None:
        final_omegas.append(omega_f)
        final_gammas.append(gamma_f)
    else:
        final_omegas.append(np.nan)
        final_gammas.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------------- #

final_omegas_neo_rhostar_1e_2 = []
final_gammas_neo_rhostar_1e_2 = []

for tprim in tprims:
    output_dir = f"/users/rjs659/NEO_stella/ky_0_3_tprim_scans/NEO_stella_CBC_aky_0.3_electrostatic_tprim_scan_neoclassics_kinetic_e_rhostar_1e-2/tprim_{tprim}"
    omega_f, gamma_f = get_final_omega_gamma(output_dir)
    if omega_f is not None:
        final_omegas_neo_rhostar_1e_2.append(omega_f)
        final_gammas_neo_rhostar_1e_2.append(gamma_f)
    else:
        final_omegas_neo_rhostar_1e_2.append(np.nan)
        final_gammas_neo_rhostar_1e_2.append(np.nan)

# ------------------------------------------------------------------------------------------------------------------------------------------- #

final_omegas_neo_rhostar_1e_3 = []
final_gammas_neo_rhostar_1e_3 = []

for tprim in tprims:
    output_dir = f"/users/rjs659/NEO_stella/ky_0_3_tprim_scans/NEO_stella_CBC_aky_0.3_electrostatic_tprim_scan_neoclassics_kinetic_e_rhostar_1e-3/tprim_{tprim}"
    omega_f, gamma_f = get_final_omega_gamma(output_dir)
    if omega_f is not None:
        final_omegas_neo_rhostar_1e_3.append(omega_f)
        final_gammas_neo_rhostar_1e_3.append(gamma_f)
    else:
        final_omegas_neo_rhostar_1e_3.append(np.nan)
        final_gammas_neo_rhostar_1e_3.append(np.nan)


# ------------------------------------------------------------------------------------------------------------------------------------------- #
# Read in gs2/NEO_gs2 data. 

freq = []
freq_neo_gs2 = []

growth = []
growth_neo_gs2 = []

# Define a function for reading in a file and extracting angular frequency and growth rate of the dominant mode. 
def read_file(x, freq, growth):
    print(" ... read GS2 output files ... ")
    data_main = ncdf2dict(x)
    
    # Read in relevant parameters.
    omega   = np.vstack(data_main['omega'])
    omega   = omega.reshape(-1)            

    freq.append(omega.real[len(omega) -1])
    growth.append(omega.imag[len(omega) -1])
    return freq, growth

for i in range(len(tprims)): 
    # Read in the  files.   
    freq, growth                 = read_file(f'/users/rjs659/NEO_stella/ky_0_3_tprim_scans/NEO_gs2_CBC_aky_0.3_electrostatic_tprim_scan_no_neoclassics_kinetic_e_rhostar_1e-2/tprim_{tprims[i]}/CBC.out.nc', freq, growth)
    freq_neo_gs2, growth_neo_gs2 = read_file(f'/users/rjs659/NEO_stella/ky_0_3_tprim_scans/NEO_gs2_CBC_aky_0.3_electrostatic_tprim_scan_neoclassics_kinetic_e_rhostar_1e-2/tprim_{tprims[i]}/CBC.out.nc', freq_neo_gs2, growth_neo_gs2)

# ---------------------------------------- Plot omega and gamma vs aky (side-by-side with log x-axis) ---------------------------------------- #

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
stride = 4

label_fontsize = 14
tick_fontsize = 12

ax1.plot(tprims, final_omegas, 'o-', color='dodgerblue', label="stella", markevery=(1, stride))
ax1.plot(tprims, final_omegas_neo_rhostar_1e_2, 'o-', color='gold', label=r"NEO-stella $\rho_\star = 1e-2$", markevery=(2, stride))
ax1.plot(tprims, final_omegas_neo_rhostar_1e_3, 'o-', color='purple', label=r"NEO-stella $\rho_\star = 1e-3$", markevery=(2, stride))
ax1.plot(tprims, freq, 'o-', color='crimson', label="gs2", markevery=(3, stride))
# ax1.plot(tprims, freq_neo_gs2, 'o-', color='forestgreen', label=r"NEO-gs2 $\rho_\star = 1e-2$", markevery=(4, stride))
ax1.set_xlabel(r'$a/|L_T|$', fontsize=label_fontsize)
ax1.set_ylabel(r'$\omega \; , \; a/v_{th,i}$', fontsize=label_fontsize)
ax1.grid(True, which='major', linestyle='--', alpha=0.7)
ax1.grid(False, which='minor')
# ax1.set_yscale('log')

ax2.plot(tprims, final_gammas, 'o-', color='dodgerblue', label="stella", markevery=(1, stride))
ax2.plot(tprims, final_gammas_neo_rhostar_1e_2, 'o-', color='gold', label=r"NEO-stella $\rho_\star = 1e-2$", markevery=(2, stride))
ax2.plot(tprims, final_gammas_neo_rhostar_1e_3, 'o-', color='purple', label=r"NEO-stella $\rho_\star = 1e-3$", markevery=(2, stride))
ax2.plot(tprims, growth, 'o-', color='crimson', label="gs2", markevery=(3, stride))
# ax2.plot(tprims, growth_neo_gs2, 'o-', color='forestgreen', label=r"NEO-gs2 $\rho_\star = 1e-2$", markevery=(4, stride))
ax2.set_xlabel(r'$a/|L_T|$', fontsize=label_fontsize)
ax2.set_ylabel(r'$\gamma \; , \; a/v_{th,i}$', fontsize=label_fontsize)
ax2.grid(True, which='major', linestyle='--', alpha=0.7)
ax2.grid(False, which='minor')
ax2.legend(fontsize=tick_fontsize) 
# ax2.set_yscale('log')

fig.suptitle(r'CBC, $\beta = 0.0\%$, $k_y = 0.3$, Kinetic Electrons', fontsize=label_fontsize + 2)

output_filename = "/users/rjs659/NEO_stella/ky_0_3_tprim_scans/NEO_stella_CBC_aky_0.3_electrostatic_tprim_scan_neoclassics_kinetic_e_rhostar_1e-2/omega_gamma_vs_tprim_no_neo_gs2.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"✅ Figure saved as '{output_filename}'")

