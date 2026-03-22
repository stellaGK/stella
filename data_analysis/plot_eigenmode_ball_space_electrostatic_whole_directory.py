import numpy as np
import matplotlib.pyplot as plt
import os

base_dir="/users/rjs659/NEO_stella/CBC_electrostatic_aky_scan_rhostar_0_01_no_neoclassics_adiabatic_e_implicit"

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
 
    # Calculate the absolute value (magnitude)
    phi_abs = np.sqrt(phi_real**2 + phi_imag**2)

    # Normalize the absolute value.
    phi_abs_norm = phi_abs / np.max(phi_abs)

    # Create side-by-side subplots.
    fig, ax = plt.subplots(figsize=(8, 5))

    # Plot.
    ax.plot(z_zed0, phi_abs_norm, color='black', lw=2.5, label=r"$|\phi'|/|\phi'_{\text{max}}|$")
    ax.set_xlabel(r'$z - z_0$ (Field Line Coordinate)', fontsize=13)
    ax.set_ylabel(r"Normalized Potential $|\phi'|$", fontsize=13)
    ax.set_ylim(0, 1.1) 
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=11, loc='upper right')

    plt.tight_layout()

    # Save figure
    output_file = os.path.join(sim_dir, "electrostatic_eigenmode.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_file}")
    plt.close()


# Loop through all subdirectories in base_dir.
for run_dir in os.listdir(base_dir):
    sim_dir = os.path.join(base_dir, run_dir)
    if os.path.isdir(sim_dir):
        plot_phi_eigenmode(sim_dir)
