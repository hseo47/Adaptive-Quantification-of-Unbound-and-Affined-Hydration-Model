import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import os

# Define your folder and files
n_units = 10
E_PDMS = 2.0
output_folder = f"output_n{n_units}_E{E_PDMS:.1f}".replace('.', '_')
prefix = f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_') + "_hydration_"
base_dir = os.path.join(os.path.dirname(__file__), output_folder)
states = ["free", "weak", "tight"]

# Plotting
plt.figure(figsize=(10, 6))

font_dir = os.path.join(os.path.dirname(__file__), "fonts")
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
        print(f"⚠️ Failed to set Helvetica from {font_path}: {e}")
        plt.rcParams['font.family'] = 'sans-serif'
else:
    print("⚠️ Helvetica font file not found. Using fallback.")
    plt.rcParams['font.family'] = 'sans-serif'

for state in states:
    filepath = os.path.join(base_dir, f"{prefix}{state}.csv")
    if not os.path.exists(filepath):
        print(f"❌ File missing: {filepath}")
        continue

    data = np.loadtxt(filepath, delimiter=",")
    entropy = -np.sum(data * np.log2(data + 1e-12))
    mean = np.mean(data)
    std = np.std(data)

    print(f"{state.capitalize()} - Entropy: {entropy:.4f}, Mean: {mean:.4f}, Std: {std:.4f}")
    plt.hist(data, bins=50, alpha=0.5, label=f"{state.capitalize()} (S={entropy:.2f})")

plt.title("Hydration Distribution per State")
plt.xlabel("Hydration Value")
plt.ylabel("Voxel Count")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(base_dir, "hydration_distribution_diagnostics.png"))
print(f"🖼 Saved plot to: {os.path.join(base_dir, 'hydration_distribution_diagnostics.png')}")
plt.show()