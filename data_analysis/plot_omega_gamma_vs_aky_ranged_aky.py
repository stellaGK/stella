import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# === User settings ===
filename = "input.omega"
save_formats = ["png", "pdf"]   # can include "png", "pdf", etc.

# === Read the file ===
data = pd.read_csv(
    filename,
    delim_whitespace=True,
    comment="#",
    names=["time", "ky", "kx", "Re_om", "Im_om", "Re_omavg", "Im_omavg"],
)

# === Identify final time step ===
final_time = data["time"].max()
df_final = data[data["time"] == final_time]

print(f"Plotting and saving results for final time step: t = {final_time:.5e}")

# === Plot ω and γ vs ky ===
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)

# ω (real frequency)
ax1.plot(df_final["ky"], df_final["Re_om"], "o-", color="C0")
ax1.set_ylabel(r"$\omega$ (Re[$\omega$])")
ax1.grid(True)

# γ (growth rate)
ax2.plot(df_final["ky"], df_final["Im_om"], "o-", color="C1")
ax2.set_ylabel(r"$\gamma$ (Im[$\omega$])")
ax2.set_xlabel(r"$k_y \rho_s$")
ax2.grid(True)

ax1.set_title(f"stella output: ω and γ vs ky at t = {final_time:.3e}")

plt.tight_layout()

# === Save figure ===
out_base = Path(filename).stem + "_omega_vs_ky_final"
for fmt in save_formats:
    out_path = f"{out_base}.{fmt}"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {out_path}")

plt.show()

