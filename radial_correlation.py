import os
os.chdir("/Users/hojin/Library/CloudStorage/OneDrive-SharedLibraries-GeorgiaInstituteofTechnology/Harris, Tequila A - Hojin/Literature and Data/_hydro_chitosan_kMC")
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate

# === Load hydration map ===
hydration_map = np.loadtxt("/Users/hojin/Library/CloudStorage/OneDrive-SharedLibraries-GeorgiaInstituteofTechnology/Harris, Tequila A - Hojin/Literature and Data/_hydro_chitosan_kMC/Radial Correlation/n20_E0_5_hydration_map.csv", delimiter=",")
grid_shape = (48, 28, 10)  # ← UPDATE to your actual voxel dimensions
hydration_3d = hydration_map.reshape(grid_shape)

# === Normalize ===
hydration_3d -= np.mean(hydration_3d)
hydration_3d /= np.std(hydration_3d)

# === Compute 3D autocorrelation ===
corr = correlate(hydration_3d, hydration_3d, mode='full')
center = tuple(s // 2 for s in corr.shape)
corr_y = corr[center[0], :, center[2]]  # example: slice along Y

# === Normalize correlation ===
corr_y /= np.max(corr_y)

# === Plot ===
y_range = np.arange(-grid_shape[1] + 1, grid_shape[1])
plt.figure(figsize=(8, 5))
plt.plot(y_range, corr_y, label="Y-axis radial correlation")
plt.axhline(0, color='gray', linestyle='--')
plt.title("Radial Spatial Correlation (Y Direction)")
plt.xlabel("Voxel Distance")
plt.ylabel("Normalized Correlation")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()