import numpy as np
import matplotlib.pyplot as plt
import os

base_dir="/users/rjs659/NEO_stella/NEO_stella_CBC_aky_2.0_beta_0.015_rhostar_scan_no_neoclassics"

def plot_omega_gamma_vs_time(sim_dir):
    input_file = os.path.join(sim_dir, "CBC.omega")
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        return

    # Load data
    data = np.loadtxt(input_file, comments='#')

    # Extract columns
    time = data[:, 0]
    omega = data[:, 3]       # Re[om]
    gamma = data[:, 4]       # Im[om]
    omega_avg = data[:, 5]   # Re[omavg]
    gamma_avg = data[:, 6]   # Im[omavg]

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 5))

    # Plot ω and γ
    ax.plot(time, omega, marker='o', markersize=4, mfc='none', color='forestgreen', label=r'$\omega$')
    ax.plot(time, omega_avg, linestyle='--', color='forestgreen', label=r'$\omega_{\mathrm{avg}}$')

    ax.plot(time, gamma, marker='o', markersize=4, mfc='none', color='crimson', label=r'$\gamma$')
    ax.plot(time, gamma_avg, linestyle='--', color='crimson', label=r'$\gamma_{\mathrm{avg}}$')

    # Axis labels and title
    ax.set_xlabel('Time', fontsize=14)
    ax.set_ylabel(r'$\omega,\, \gamma \; , \; a / v_{th,i}$', fontsize=14)
    ax.set_xscale('log')
    ax.legend()
    ax.grid(True)

    # Add box with final ω, γ values
    final_omega = omega[-1]
    final_gamma = gamma[-1]
    textstr = f'ω = {final_omega:.3f}\nγ = {final_gamma:.3f}'

    # Text box style
    props = dict(boxstyle='round,pad=0.5', facecolor='#FFFDD0', edgecolor='tab:blue', linewidth=1.5, alpha=0.9)
    ax.text(0.85, 0.05, textstr, transform=ax.transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='right', bbox=props)

    plt.tight_layout()

    output_file = os.path.join(sim_dir, "omega_gamma_vs_time.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_file}")
    plt.close()


# Loop through all subdirectories in base_dir.
for run_dir in os.listdir(base_dir):
    sim_dir = os.path.join(base_dir, run_dir)
    if os.path.isdir(sim_dir):
        plot_omega_gamma_vs_time(sim_dir)

