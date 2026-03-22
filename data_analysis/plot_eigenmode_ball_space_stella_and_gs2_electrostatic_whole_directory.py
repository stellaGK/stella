import numpy as np
import matplotlib.pyplot as plt
import os
from ncdf2dict import ncdf2dict

stella_base_dir = "/users/rjs659/NEO_stella/fprim_scans/fprim_scans/ky_0.3/kinetic_electrons/electrostatic/NEO_stella_no_neoclassics_data"
gs2_base_dir    = "/users/rjs659/NEO_stella/fprim_scans/fprim_scans/ky_0.3/kinetic_electrons/electrostatic/NEO_gs2_no_neoclassics_data"

def plot_phi_eigenmode(stella_sim_dir, gs2_sim_dir):
    # ---------------------------------------------------------------------------------------------------------------------- #
    # Get stella result.

    input_file = os.path.join(stella_sim_dir, "CBC.final_fields")
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        return

    # Load data.
    data = np.loadtxt(input_file, comments='#')

    # Extract relevant data.
    z_zed0 = data[:, 1]      
    phi_real = data[:, 4]    
    phi_imag = data[:, 5]    
 
    # Calculate the absolute value (magnitude).
    phi_abs = np.sqrt(phi_real**2 + phi_imag**2)

    # Normalize the absolute value.
    phi_abs_norm_stella = phi_abs / np.max(phi_abs)

    # ---------------------------------------------------------------------------------------------------------------------- #
    # Get gs2 result.

    ncfile = os.path.join(gs2_sim_dir, "CBC.out.nc")
    if not os.path.exists(ncfile):
        print(f"Missing file: {ncfile}")
        return

    print(f" ... processing {gs2_sim_dir} ... ")
    data_main = ncdf2dict(ncfile)

    theta_theta0 = data_main['theta']
    omega = np.vstack(data_main['omega']).reshape(-1)

    phi = np.vstack(data_main['phi']).reshape(-1)
    phi_real = phi.real
    phi_imag = phi.imag   

    # Calculate the absolute value and normalise to the peak. 
    phi_abs = np.sqrt(phi_real**2 + phi_imag**2)
    phi_abs_max =np.max(phi_abs)
    phi_abs_norm_gs2  = phi_abs  / phi_abs_max

    # ---------------------------------------------------------------------------------------------------------------------- #
    # Plot results.
    fig, ax = plt.subplots(figsize=(8, 5))

    # Plot.
    ax.plot(z_zed0, phi_abs_norm_stella, color='dodgerblue', lw=2.5, label=r"stella")
    ax.plot(theta_theta0, phi_abs_norm_gs2, color='crimson', lw=2.5, label=r"gs2")
    ax.set_xlabel(r'$z - z_0$, $\theta - \theta_0$, (Field Line Coordinate)', fontsize=13)
    ax.set_ylabel(r"$|\phi'|/|\phi'_{\text{max}}|$", fontsize=13)
    ax.set_ylim(0, 1.1) 
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(fontsize=11, loc='upper right')

    plt.tight_layout()

    # Save figure
    output_file = os.path.join(stella_sim_dir, "electrostatic_eigenmode.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_file}")
    plt.close()

# ---------------------------------------------------------------------------------------- #
# Loop through all subdirectories in stella and gs2 base_dir.
for run_dir in os.listdir(stella_base_dir):
    # Construct the full path for the stella subdirectory.
    sim_dir_stella = os.path.join(stella_base_dir, run_dir)
    
    # Construct the matching path for the gs2 subdirectory.
    sim_dir_gs2 = os.path.join(gs2_base_dir, run_dir)

    # Only run the function if it's actually a directory and exists in both places.
    if os.path.isdir(sim_dir_stella) and os.path.isdir(sim_dir_gs2):
        plot_phi_eigenmode(sim_dir_stella, sim_dir_gs2)
    else:
        print(f"Skipping {run_dir}: matching directory not found in both locations.")
