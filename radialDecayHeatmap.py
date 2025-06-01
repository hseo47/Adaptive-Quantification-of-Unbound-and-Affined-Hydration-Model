import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.cm as cm
from matplotlib import ticker

import os
base_dir = os.path.dirname(os.path.abspath(__file__))
font_path = os.path.join(base_dir, "fonts", "Computer Modern.ttf")
cm_font = fm.FontProperties(fname=font_path)
if os.path.exists(font_path):
    plt.rcParams['font.family'] = cm_font.get_name()
plt.rcParams.update({
    'axes.titlesize': 28,
    'axes.labelsize': 24,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 20,
    'figure.titlesize': 28,
    'axes.unicode_minus': False
})

import glob

base_path = os.path.dirname(__file__)
files = []

target_Es = {0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0}
import re
for fpath in sorted(glob.glob(os.path.join(base_path, "output_n1_E*","n1_E*_hydration_map.csv"))):
    try:
        match = re.search(r"_E([0-9_]+)", fpath)
        if not match:
            raise ValueError("Could not extract modulus from filename")
        E_str = match.group(1).replace("_", ".")
        E_value = float(E_str)
        if E_value not in target_Es:
            continue
        files.append((fpath, E_value))
    except Exception as e:
        print(f"Skipping {fpath}: {e}")

def load_voxel_shape(hydration_file):
    base_name = os.path.basename(hydration_file).replace("_hydration_map.csv", "_hyd_distrib.txt")
    distrib_path = os.path.join(os.path.dirname(hydration_file), base_name)
    with open(distrib_path) as f:
        lines = f.readlines()
        x_voxels = int([l for l in lines if "x_voxels" in l][0].split(":")[1])
        y_voxels = int([l for l in lines if "y_voxels" in l][0].split(":")[1])
        z_voxels = int([l for l in lines if "z_voxels" in l][0].split(":")[1])
    return (x_voxels, y_voxels, z_voxels)

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
cmap = plt.cm.inferno  

import pandas as pd

modulus_list = sorted([E for _, E in files])
max_radii = []

for fname, _ in files:
    grid_shape = load_voxel_shape(fname)
    hydration_map = np.loadtxt(fname, delimiter=",")
    hydration_3d = hydration_map.reshape(grid_shape)
    hydration_3d = (hydration_3d - global_min) / (global_max - global_min)

    max_idx = np.unravel_index(np.argmax(hydration_3d), hydration_3d.shape)
    zz, yy, xx = np.indices(hydration_3d.shape)
    coords = np.stack((xx, yy, zz), axis=-1)
    radii = np.linalg.norm(coords - np.array(max_idx), axis=-1)

    radii_flat = radii.flatten()
    hydration_flat = hydration_3d.flatten()
    max_radius = int(np.max(radii_flat)) + 1
    max_radii.append(max_radius)

heatmap_data = np.zeros((len(modulus_list), max(max_radii)))

for idx, (fname, E) in enumerate(sorted(files, key=lambda x: x[1])):
    grid_shape = load_voxel_shape(fname)
    hydration_map = np.loadtxt(fname, delimiter=",")
    hydration_3d = hydration_map.reshape(grid_shape)
    hydration_3d = (hydration_3d - global_min) / (global_max - global_min)

    max_idx = np.unravel_index(np.argmax(hydration_3d), hydration_3d.shape)
    zz, yy, xx = np.indices(hydration_3d.shape)
    coords = np.stack((xx, yy, zz), axis=-1)
    radii = np.linalg.norm(coords - np.array(max_idx), axis=-1)

    radii_flat = radii.flatten()
    hydration_flat = hydration_3d.flatten()
    max_radius = int(np.max(radii_flat)) + 1

    bin_means = np.zeros(max_radius)
    bin_counts = np.zeros(max_radius)

    for r, h in zip(radii_flat, hydration_flat):
        idx_r = int(np.floor(r))
        bin_means[idx_r] += h
        bin_counts[idx_r] += 1

    bin_counts[bin_counts == 0] = 1
    bin_means /= bin_counts
    heatmap_data[idx, :len(bin_means)] = bin_means

plt.figure(figsize=(8, 8))
im = plt.imshow(
    heatmap_data,
    aspect='auto',
    cmap="Blues",
    origin='lower',
    extent=[0, heatmap_data.shape[1], modulus_list[0], modulus_list[-1]],
    interpolation='nearest'
)

cbar = plt.colorbar(im)
cbar.set_label("Normalized Hydration", fontproperties=cm_font, fontsize=24)
cbar.ax.tick_params(labelsize=20)
for label in cbar.ax.get_yticklabels():
    label.set_fontproperties(cm_font)
    label.set_fontsize(20)
plt.xlabel("Radius (Voxels)", fontproperties=cm_font, fontsize=24)
plt.xticks(fontproperties=cm_font, fontsize=20)
plt.ylabel("Elastic Modulus (MPa)", fontproperties=cm_font, fontsize=24)
plt.yticks(fontproperties=cm_font, fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.title("Hydration Field Heatmap", fontproperties=cm_font, fontsize=28)
plt.tight_layout()
plt.savefig("hydration_heatmap.png", dpi=1200, bbox_inches='tight', pad_inches=0.05, transparent=True)
plt.show()
