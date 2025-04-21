import numpy as np
import matplotlib.pyplot as plt

# === Configuration ===
hydration_file = " " # hydration map file
grid_shape = (48, 28, 10)  # (X, Y, Z)

# === Load and reshape ===
hydration_map = np.loadtxt(hydration_file, delimiter=",")
hydration_3d = hydration_map.reshape(grid_shape)
hydration_3d /= np.max(hydration_3d)  # Normalize

# === Find max location ===
max_idx = np.unravel_index(np.argmax(hydration_3d), hydration_3d.shape)
x0, y0, z0 = max_idx

# === Extract decay profiles ===
x_profile = hydration_3d[:, y0, z0]
y_profile = hydration_3d[x0, :, z0]
z_profile = hydration_3d[x0, y0, :]

# === Mirror around center for symmetry ===
x_range = np.arange(-x0, grid_shape[0] - x0)
y_range = np.arange(-y0, grid_shape[1] - y0)
z_range = np.arange(-z0, grid_shape[2] - z0)

# === Plot ===
plt.figure(figsize=(8, 5))
plt.plot(x_range, x_profile, label='X-axis', linewidth=2)
plt.plot(y_range, y_profile, label='Y-axis', linewidth=2)
plt.plot(z_range, z_profile, label='Z-axis', linewidth=2)
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)

plt.xlabel("Voxel Distance from Max Hydration")
plt.ylabel("Normalized Hydration")
plt.title("Directional Hydration Decay from Max Site (Anisotropy)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
