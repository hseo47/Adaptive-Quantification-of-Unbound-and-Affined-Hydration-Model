import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# === Load hydration map ===
base_path = "" # load base path
files = [
    (base_path + "n20_E0_5_hydration_map.csv", 0.5),
    (base_path + "n20_E0_75_hydration_map.csv", 0.75),
    (base_path + "n20_E1_0_hydration_map.csv", 1.0),
    (base_path + "n20_E1_25_hydration_map.csv", 1.25),
    (base_path + "n20_E1_5_hydration_map.csv", 1.5),
    (base_path + "n20_E1_75_hydration_map.csv", 1.75),
    (base_path + "n20_E2_0_hydration_map.csv", 2.0),
]

grid_shape = (48, 28, 10)  # Ensure this matches your simulation

# === Compute global min and max ===
global_min = float('inf')
global_max = float('-inf')
for fname, _ in files:
    data = np.loadtxt(fname, delimiter=",")
    global_min = min(global_min, np.min(data))
    global_max = max(global_max, np.max(data))

from matplotlib.colors import LinearSegmentedColormap

hydration_blues = LinearSegmentedColormap.from_list("hydration_blues", ["#a6c8ff", "#003f5c"])
cmap = hydration_blues
norm = plt.Normalize(vmin=0.5, vmax=2.0)

plt.figure(figsize=(8, 5))

for fname, E in files:
    hydration_map = np.loadtxt(fname, delimiter=",")
    hydration_3d = hydration_map.reshape(grid_shape)

    hydration_3d = (hydration_3d - global_min) / (global_max - global_min)

    max_idx = np.unravel_index(np.argmax(hydration_3d), hydration_3d.shape)
    max_coords = np.array(max_idx)

    zz, yy, xx = np.indices(hydration_3d.shape)
    coords = np.stack((xx, yy, zz), axis=-1)
    radii = np.linalg.norm(coords - max_coords, axis=-1)

    radii_flat = radii.flatten()
    hydration_flat = hydration_3d.flatten()
    max_radius = int(np.max(radii_flat)) + 1

    bin_means = np.zeros(max_radius)
    bin_counts = np.zeros(max_radius)

    for r, h in zip(radii_flat, hydration_flat):
        idx = int(np.floor(r))
        bin_means[idx] += h
        bin_counts[idx] += 1

    bin_counts[bin_counts == 0] = 1
    bin_means /= bin_counts

    linewidth = 2
    linestyle = '-'

    color = cmap(norm(E))
    plt.plot(np.arange(max_radius), bin_means, label=f"E = {E} MPa", color=color, linewidth=linewidth, linestyle=linestyle)

plt.xlabel("Radius (voxels)")
plt.ylabel("Normalized Hydration")
plt.title("3D Radial Hydration Decay (All Moduli)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
