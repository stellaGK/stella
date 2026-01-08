import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt('/users/rjs659/stella-v1.0/input.omega', comments='#')

# Extract columns
time = data[:, 0]
omega = data[:, 3]       # Re[om]
gamma = data[:, 4]       # Im[om]
omega_avg = data[:, 5]   # Re[omavg]
gamma_avg = data[:, 6]   # Im[omavg]

# Create figure
fig, ax = plt.subplots(figsize=(8, 5))

# --- Plot ω and γ ---
ax.plot(time, omega, marker='o', markersize=5, mfc='none', color='forestgreen', label=r'$\omega$')
ax.plot(time, omega_avg, linestyle='--', color='forestgreen', label=r'$\omega_{\mathrm{avg}}$')

ax.plot(time, gamma, marker='o', markersize=5, mfc='none', color='crimson', label=r'$\gamma$')
ax.plot(time, gamma_avg, linestyle='--', color='crimson', label=r'$\gamma_{\mathrm{avg}}$')

# Axis labels and title
ax.set_xlabel('Time')
ax.set_ylabel(r'$\omega,\, \gamma \; , \; a / v_{th,i}$')
ax.set_title(r'Electrostatic CBC, $k_y = 0.1$')
ax.legend()
ax.grid(True)

# --- Add box with final values ---
final_omega = omega[-1]
final_gamma = gamma[-1]

textstr = f'ω = {final_omega:.3f}\nγ = {final_gamma:.3f}'

# Text box properties.
props = dict(boxstyle='round,pad=0.5', facecolor='#FFFDD0', edgecolor='tab:blue', linewidth=1.5, alpha=0.9)

# Place text box in bottom-right (axes coordinates)
ax.text(0.85, 0.05, textstr, transform=ax.transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='right', bbox=props)

# Layout and save
plt.tight_layout()
output_filename = 'omega_gamma_vs_time_aky_0.3.png'
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
plt.show()

print(f"Plot saved as {output_filename}")

