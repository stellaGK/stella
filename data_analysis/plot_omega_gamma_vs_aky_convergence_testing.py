import numpy as np
import matplotlib.pyplot as plt
import os

def get_final_omega_gamma(output_dir):
    """
    Read input.omega in the specified directory and return the final omega and gamma.
    """
    omega_file = os.path.join(output_dir, "input.omega")
    if not os.path.exists(omega_file):
        print(f"{omega_file} not found!")
        return None, None

    try:
        data = np.loadtxt(omega_file, comments='#')
        omega = data[:, 3]  # Re[om]
        gamma = data[:, 4]  # Im[om]
        return omega[-1], gamma[-1]
    except Exception as e:
        print(f"Error reading {omega_file}: {e}")
        return None, None


# --- User Configuration ---

# Same aky values for all directories
aky_values = [0.5, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]

# Four separate base directories
base_dirs = ["nzed_local_16", "nzed_local_32", "nzed_local_64", "nzed_local_128"]

# --- Plot setup ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
label_fontsize = 14
tick_fontsize = 12
colors = plt.cm.tab10(np.linspace(0, 1, len(base_dirs)))

# --- Process each directory ---
for base_dir, color in zip(base_dirs, colors):
    final_omegas = []
    final_gammas = []

    for aky in aky_values:
        output_dir = os.path.join(base_dir, f"aky_{aky}")
        omega_f, gamma_f = get_final_omega_gamma(output_dir)
        if omega_f is not None:
            final_omegas.append(omega_f)
            final_gammas.append(gamma_f)
        else:
            final_omegas.append(np.nan)
            final_gammas.append(np.nan)

    label = os.path.basename(base_dir)
    ax1.plot(aky_values, final_omegas, 'o-', label=label, color=color)
    ax2.plot(aky_values, final_gammas, 'o-', label=label, color=color)

# --- Formatting ---
for ax, ylabel in zip([ax1, ax2], [r'$\omega \; , \; a/v_{th,i}$', r'$\gamma \; , \; a/v_{th,i}$']):
    # ax.set_xscale('log')
    ax.set_xlabel(r'$k_y\rho_i$', fontsize=label_fontsize)
    ax.set_ylabel(ylabel, fontsize=label_fontsize)
    ax.grid(True, which='both', linestyle='--', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)
    ax.tick_params(axis='x', which='minor', bottom=True)
    ax.legend(fontsize=10)

plt.tight_layout()

# --- Save the figure ---#
output_filename = "convergence_test_nzed_local.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"✅ Figure saved as '{output_filename}'")

