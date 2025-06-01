import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy.stats import entropy
import matplotlib
matplotlib.use('Agg')

base_dir = os.path.dirname(os.path.abspath(__file__))
font_path = os.path.join(base_dir, "fonts", "Computer Modern.ttf")

if os.path.exists(font_path):
    cm_font = fm.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = cm_font.get_name()
else:
    cm_font = None

plt.rcParams.update({
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'figure.titlesize': 18,
    'axes.unicode_minus': False
})

modulus_dirs = sorted([d for d in os.listdir(base_dir) if d.startswith("output_n1_E")])
moduli = [float(d.split("_E")[1].replace("_", ".")) for d in modulus_dirs]

hydration_data = {}

for mod, d in zip(moduli, modulus_dirs):
    map_file = os.path.join(base_dir, d, f"n1_E{str(mod).replace('.', '_')}_hydration_map.csv")
    if os.path.exists(map_file):
        values = np.loadtxt(map_file, delimiter=",")
        hydration_data[mod] = values

shannon_entropies = {}
all_values = np.concatenate(list(hydration_data.values()))
bins = np.histogram_bin_edges(all_values, bins=50)

for mod, values in hydration_data.items():
    hist, _ = np.histogram(values, bins=bins, density=True)
    hist = hist[hist > 0]  # Avoid log(0)
    shannon_entropies[mod] = entropy(hist)

reference_mod = 2.0
ref_values = hydration_data[reference_mod]
E_cutoff = np.percentile(ref_values, 95)

high_energy_prevalence = {}
for mod, values in hydration_data.items():
    high_voxel_count = np.sum(values > E_cutoff)
    high_energy_prevalence[mod] = high_voxel_count / len(values)

plt.figure(figsize=(5, 5))
plt.plot(
    list(shannon_entropies.keys()),
    list(shannon_entropies.values()),
    marker='o',
    color='black',
    label='Shannon Entropy'
)
plt.xlabel("Elastic Modulus (MPa)", fontproperties=cm_font)
plt.ylabel("Entropy", fontproperties=cm_font)
plt.title("Global Shannon Entropy", fontproperties=cm_font)
plt.grid(False)
plt.xticks(fontproperties=cm_font)
plt.yticks(fontproperties=cm_font)
plt.tight_layout()
plt.savefig(os.path.join(base_dir, "shannon_entropy_vs_modulus.png"), dpi=900)

plt.figure(figsize=(5, 5))
plt.plot(
    list(high_energy_prevalence.keys()),
    list(high_energy_prevalence.values()),
    marker='o',
    color='black',
    label='High-Energy Prevalence'
)
plt.xlabel("Elastic Modulus (MPa)", fontproperties=cm_font)
plt.ylabel("High-Energy Prevalence", fontproperties=cm_font)
plt.title("Tail Behavior", fontproperties=cm_font)
plt.grid(False)
plt.xticks(fontproperties=cm_font)
plt.yticks(fontproperties=cm_font)
plt.tight_layout()
plt.savefig(os.path.join(base_dir, "high_energy_prevalence_vs_modulus.png"), dpi=900)
