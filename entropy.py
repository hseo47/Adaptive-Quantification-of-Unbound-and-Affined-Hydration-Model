import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import entropy

# === Parameters ===
n_units = 10
E_PDMS = 2.0
foldername = f"output_n{n_units}_E{E_PDMS}".replace(".", "_")

# === Base path ===
base_path = os.path.join(os.path.dirname(__file__), foldername)

# === Files ===
files = [
    (f"n{n_units}_E{E_PDMS}_hydration_map.csv".replace(".", "_").replace("_csv", ".csv"), E_PDMS),
]

E_vals = []
entropy_vals = []

for fname, E in files:
    path = os.path.join(base_path, fname)
    try:
        data = np.loadtxt(path, delimiter=",")
        grid_size = data.size
        inferred_shape = (grid_size, )
        if grid_size == np.prod((48, 28, 10)):
            data = data.reshape((48, 28, 10))
        else:
            print(f"⚠️ Warning: Expected {np.prod((48, 28, 10))} elements but got {grid_size}. Using flat shape.")
        flattened = data.flatten()
        hist, _ = np.histogram(flattened, bins=100, density=True)
        hist += 1e-12  # prevent log(0)
        shannon = entropy(hist)

        E_vals.append(E)
        entropy_vals.append(shannon)
        print(f"✅ {fname} -> Entropy = {shannon:.4f}")

    except FileNotFoundError:
        print(f"❌ File not found: {fname}")

# === Helvetica Font ===
from matplotlib import font_manager
font_path = os.path.join(os.path.dirname(__file__), "fonts", "Helvetica.ttc")
helvetica_prop = font_manager.FontProperties(fname=font_path)
plt.rcParams['font.family'] = helvetica_prop.get_name()

# === Plot ===
plt.figure(figsize=(7, 4))
plt.plot(E_vals, entropy_vals, marker='o', linestyle='-', color='#005F9E')
plt.xlabel("Substrate Modulus (MPa)")
plt.ylabel("Shannon Entropy")
plt.title("Hydration Field Entropy vs. Modulus")
plt.grid(True)
plt.tight_layout()
os.makedirs(base_path, exist_ok=True)
save_path = os.path.join(base_path, "entropy_vs_modulus.png")
plt.savefig(save_path, dpi=300)
print(f"✅ Saved plot to: {save_path}")