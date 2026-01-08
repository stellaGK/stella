import numpy as np
import matplotlib.pyplot as plt

# Load data (skip header lines starting with '#').
data = np.loadtxt('input.final_fields', comments='#')

# Extract columns
z_zed0 = data[:, 1]      # z - zed0.
phi_real = data[:, 4]    # real(phi).
phi_imag = data[:, 5]    # imag(phi).

# Find the normalization factor (max absolute value in real or imag). 
phi_max = max(np.max(np.abs(phi_real)), np.max(np.abs(phi_imag)))

# Normalize the data.
phi_real_norm = phi_real / phi_max
phi_imag_norm = phi_imag / phi_max

# Plotting.
plt.figure(figsize=(8,5))
plt.plot(z_zed0, phi_real_norm, 'o-', mfc='none', color='tab:blue', label=r"Re($\phi'$)")
plt.plot(z_zed0, phi_imag_norm, 'o-', mfc='none', color='tab:orange', label=r"Im($\phi'$)")

# Labels and title.
plt.xlabel(r'$z - z_0$')
plt.ylabel(r"$\phi' / \phi'_{\mathrm{max}}$")
plt.title(r'Electrostatic CBC, $k_y = 0.1$')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save figure.
output_filename = 'electrostatic_CBC_aky_0.3_eigenmode.png'
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved as {output_filename}")
