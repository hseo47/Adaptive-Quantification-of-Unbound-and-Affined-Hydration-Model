import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# === Load voxel hydration data ===
base_dir = os.path.dirname(os.path.abspath(__file__))
n_units = 20
E_PDMS = 0.5
filename_suffix = f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_')

hydration_path = os.path.join(base_dir, f"{filename_suffix}_hydration_map.csv")
class_map_path = os.path.join(base_dir, f"{filename_suffix}_hydration_class_map.csv")

hydration_flat = np.loadtxt(hydration_path, delimiter=",")
class_flat = np.loadtxt(class_map_path, delimiter=",", dtype=int)

# === Reshape based on saved voxel dimensions (assume stored in distrib txt) ===
distrib_path = os.path.join(base_dir, f"{filename_suffix}_hyd_distrib.txt")
reference_distrib_path = os.path.join(base_dir, "n20_E2_0_hyd_distrib.txt")
global_max = None
with open(reference_distrib_path) as f:
    for line in f:
        if line.startswith("max_hydration"):
            global_max = float(line.strip().split(":")[1])
            break

x_voxels = y_voxels = z_voxels = None
with open(distrib_path) as f:
    for line in f:
        if line.startswith("x_voxels"):
            x_voxels = int(line.strip().split(":")[1])
        elif line.startswith("y_voxels"):
            y_voxels = int(line.strip().split(":")[1])
        elif line.startswith("z_voxels"):
            z_voxels = int(line.strip().split(":")[1])

shape = (x_voxels, y_voxels, z_voxels)
hydration_grid = hydration_flat.reshape(shape)
class_grid = class_flat.reshape(shape)

# === Create a central slice (XY) ===
central_z = z_voxels // 2
hydration_slice = hydration_grid[:, :, central_z]
class_slice = class_grid[:, :, central_z]

# === Find center of maximum hydration in slice ===
max_idx = np.unravel_index(np.argmax(hydration_slice), hydration_slice.shape)
center_x, center_y = max_idx

# === Generate radial profile from center ===
def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int32)
    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / np.maximum(nr, 1)
    return radialprofile, np.arange(len(radialprofile))

hydration_radial, r_vals = radial_profile(hydration_slice, (center_x, center_y))

# === Fit Gaussian to radial profile ===
def gaussian(r, A, sigma, offset):
    return A * np.exp(-(r**2) / (2 * sigma**2)) + offset

p0 = [hydration_radial.max(), 5, hydration_radial.min()]
popt, _ = curve_fit(gaussian, r_vals, hydration_radial, p0=p0)

# === Plot results ===
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
vlim = global_max if global_max else None
plt.imshow(hydration_slice.T, cmap='Blues', origin='lower', vmin=0, vmax=vlim)
plt.title("Hydration Slice (Z = middle)")
plt.colorbar(label="Hydration Potential")
plt.scatter(center_x, center_y, color='red', label='Max Hydration')
plt.xlabel("X (Å)")
plt.ylabel("Y (Å)")
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(r_vals, hydration_radial, 'b-', label='Radial Profile')
plt.plot(r_vals, gaussian(r_vals, *popt), 'r--', label=f'Fit: σ = {popt[1]:.2f} Å')
plt.xlabel("Radius (voxels)")
plt.ylabel("Hydration")
plt.title("Radial Hydration Decay")
plt.legend()
plt.tight_layout()

# Save and show
output_path = os.path.join(base_dir, f"{filename_suffix}_shielding_radius_fit.png")
plt.savefig(output_path, dpi=300)
print(f"✅ Saved shielding radius visualization to: {output_path}")
plt.show()