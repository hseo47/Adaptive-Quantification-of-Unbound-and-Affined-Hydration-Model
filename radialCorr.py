import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate
import matplotlib.font_manager as fm

# === Parameters ===
n_units = 10
E_PDMS = 2.0

# === Directory Setup ===
base_dir = os.path.dirname(os.path.abspath(__file__))
output_folder = f"output_n{n_units}_E{E_PDMS:.1f}".replace('.', '_')
output_path = os.path.join(base_dir, output_folder)

distrib_filename = os.path.join(output_path, f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_') + "_hyd_distrib.txt")
try:
    with open(distrib_filename) as f:
        lines = f.readlines()
        x_voxels = int([l for l in lines if "x_voxels" in l][0].split(":")[1])
        y_voxels = int([l for l in lines if "y_voxels" in l][0].split(":")[1])
        z_voxels = int([l for l in lines if "z_voxels" in l][0].split(":")[1])
    grid_shape = (x_voxels, y_voxels, z_voxels)
except Exception as e:
    raise ValueError(f"⚠️ Could not determine grid shape from {distrib_filename}: {e}")

hydration_filename = os.path.join(output_path, f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_') + "_hydration_map.csv")

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
    except Exception as e:
        print(f"⚠️ Failed to load Helvetica: {e}")
        plt.rcParams['font.family'] = 'sans-serif'
else:
    print("⚠️ Helvetica font not found. Falling back to sans-serif.")
    plt.rcParams['font.family'] = 'sans-serif'

# === Load hydration map ===
hydration_map = np.loadtxt(hydration_filename, delimiter=",")
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
plt.plot(y_range, corr_y, label="Y-axis radial correlation", color='black')
plt.axhline(0, color='gray', linestyle='--')
plt.title("Radial Spatial Correlation (Y Direction)")
plt.xlabel("Voxel Distance")
plt.ylabel("Normalized Correlation")
plt.legend()
plt.tight_layout()

save_path = os.path.join(output_path, f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_') + "_radialCorr.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"✅ Saved radial correlation plot to: {save_path}")
plt.show()