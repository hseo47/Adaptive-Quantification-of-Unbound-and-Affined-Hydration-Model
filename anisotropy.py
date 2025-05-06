import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# === Configuration ===
n_units = 10
E_PDMS = 2.0

base_dir = os.path.dirname(os.path.abspath(__file__))
output_folder = os.path.join(base_dir, f"output_n{n_units}_E{E_PDMS:.1f}".replace(".", "_"))
hydration_file = os.path.join(output_folder, f"n{n_units}_E{E_PDMS:.1f}".replace(".", "_") + "_hydration_map.csv")
distrib_file = os.path.join(output_folder, f"n{n_units}_E{E_PDMS:.1f}".replace(".", "_") + "_hyd_distrib.txt")
try:
    with open(distrib_file) as f:
        lines = f.readlines()
        x_voxels = int([l for l in lines if "x_voxels" in l][0].split(":")[1])
        y_voxels = int([l for l in lines if "y_voxels" in l][0].split(":")[1])
        z_voxels = int([l for l in lines if "z_voxels" in l][0].split(":")[1])
        grid_shape = (x_voxels, y_voxels, z_voxels)
except Exception as e:
    raise RuntimeError(f"Failed to parse grid shape from distribution file: {e}")

# === Font Setup ===
font_dir = os.path.join(base_dir, "fonts")
possible_fonts = ["Helvetica.ttf", "Helvetica-Regular.ttf", "Helvetica.ttc"]
font_path = None
for fname in possible_fonts:
    test_path = os.path.join(font_dir, fname)
    if os.path.exists(test_path):
        font_path = test_path
        break

if font_path:
    try:
        helvetica_prop = fm.FontProperties(fname=font_path)
        plt.rcParams['font.family'] = helvetica_prop.get_name()
    except Exception:
        plt.rcParams['font.family'] = 'sans-serif'
else:
    plt.rcParams['font.family'] = 'sans-serif'

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
# Limit ranges and profiles to values between -5 and 5
x_mask = (x_range >= -5) & (x_range <= 5)
y_mask = (y_range >= -5) & (y_range <= 5)
z_mask = (z_range >= -5) & (z_range <= 5)
x_range, x_profile = x_range[x_mask], x_profile[x_mask]
y_range, y_profile = y_range[y_mask], y_profile[y_mask]
z_range, z_profile = z_range[z_mask], z_profile[z_mask]

# === Plot ===
plt.figure(figsize=(8, 5))
muted_colors = ['#4A6FA5', '#A55A5A', '#D4AF37']
plt.plot(x_range, x_profile, label='X-axis', linewidth=2, color=muted_colors[0])
plt.plot(y_range, y_profile, label='Y-axis', linewidth=2, color=muted_colors[1])
plt.plot(z_range, z_profile, label='Z-axis', linewidth=2, color=muted_colors[2])
plt.axhline(0, color='gray', linestyle='--', linewidth=2, label='Zero Baseline')

plt.xlabel("Voxel Distance from Max Hydration")
plt.ylabel("Normalized Hydration")
plt.title("Directional Hydration Decay from Max Site (Anisotropy)")
plt.legend(frameon=False)
plt.tight_layout()
plt.show()