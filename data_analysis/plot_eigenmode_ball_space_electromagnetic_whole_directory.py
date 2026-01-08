import numpy as np
import matplotlib.pyplot as plt
import os

base_dir="/users/rjs659/NEO_stella/NEO_stella_CBC_aky_2.0_beta_0.015_rhostar_scan_neoclassics"

def plot_phi_eigenmode(sim_dir):
    input_file = os.path.join(sim_dir, "CBC.final_fields")
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        return

    # Load data
    data = np.loadtxt(input_file, comments='#')

    # Extract relevant data.
    z_zed0 = data[:, 1]      
    phi_real = data[:, 4]    
    phi_imag = data[:, 5]    
    apar_real = data[:, 6]    
    apar_imag = data[:, 7]    
    bpar_real = data[:, 8]    
    bpar_imag = data[:, 9]    

    # Normalize all datasets to the peak of Phi' (either real or imaginary parts, whichever is largest).
    phi_max = max(np.max(np.abs(phi_real)), np.max(np.abs(phi_imag)))
    phi_real_norm = phi_real / phi_max
    phi_imag_norm = phi_imag / phi_max
    apar_real_norm = apar_real / phi_max
    apar_imag_norm = apar_imag / phi_max
    bpar_real_norm = bpar_real / phi_max
    bpar_imag_norm = bpar_imag / phi_max

    # Create side-by-side subplots.
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True)
    fields = [
        (phi_real_norm, phi_imag_norm, r"$\phi'$", "tab:blue", "tab:orange"),
        (apar_real_norm, apar_imag_norm, r"$A_\parallel$", "tab:blue", "tab:orange"),
        (bpar_real_norm, bpar_imag_norm, r"$B_\parallel$", "tab:blue", "tab:orange"),
    ]

    # Plot each field in its own subplot.
    for ax, (real_part, imag_part, label, color_r, color_i) in zip(axes, fields):
        ax.plot(z_zed0, real_part, 'o-', mfc='none', color=color_r, label=f"Re({label})")
        ax.plot(z_zed0, imag_part, 'o-', mfc='none', color=color_i, label=f"Im({label})")
        ax.set_xlabel(r'$z - z_0$', fontsize=13)
        ax.set_ylabel(f"{label} / φ'_max", fontsize=13)
        ax.legend(fontsize=10)
        ax.grid(True)

    # Overall title and layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Save figure
    output_file = os.path.join(sim_dir, "electromagnetic_eigenmodes.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_file}")
    plt.close()


# Loop through all subdirectories in base_dir.
for run_dir in os.listdir(base_dir):
    sim_dir = os.path.join(base_dir, run_dir)
    if os.path.isdir(sim_dir):
        plot_phi_eigenmode(sim_dir)
