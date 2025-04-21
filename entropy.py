import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import entropy

# === Base path ===
base_path = "/Users/hojin/Library/CloudStorage/OneDrive-SharedLibraries-GeorgiaInstituteofTechnology/Harris, Tequila A - Hojin/Literature and Data/_hydro_chitosan_kMC/Entropy/"

# === Files ===
files = [
    ("n20_E0_5_hydration_map.csv", 0.5),
    ("n20_E0_75_hydration_map.csv", 0.75),
    ("n20_E1_0_hydration_map.csv", 1.0),
    ("n20_E1_25_hydration_map.csv", 1.25),
    ("n20_E1_5_hydration_map.csv", 1.5),
    ("n20_E1_75_hydration_map.csv", 1.75),
    ("n20_E2_0_hydration_map.csv", 2.0),
]

grid_shape = (48, 28, 10)

E_vals = []
entropy_vals = []

for fname, E in files:
    path = os.path.join(base_path, fname)
    try:
        data = np.loadtxt(path, delimiter=",").reshape(grid_shape)
        flattened = data.flatten()
        hist, _ = np.histogram(flattened, bins=100, density=True)
        hist += 1e-12  # prevent log(0)
        shannon = entropy(hist)

        E_vals.append(E)
        entropy_vals.append(shannon)
        print(f"✅ {fname} -> Entropy = {shannon:.4f}")

    except FileNotFoundError:
        print(f"❌ File not found: {fname}")

# === Plot ===
plt.figure(figsize=(7, 4))
plt.plot(E_vals, entropy_vals, marker='o', linestyle='-', color='#005F9E')
plt.xlabel("Substrate Modulus (MPa)")
plt.ylabel("Shannon Entropy")
plt.title("Hydration Field Entropy vs. Modulus")
plt.grid(True)
plt.tight_layout()
plt.show()