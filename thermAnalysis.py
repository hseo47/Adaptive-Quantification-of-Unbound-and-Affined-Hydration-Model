import numpy as np
import os
from scipy.interpolate import make_interp_spline

# === Configuration parameters (manually set here) ===
n_units = 10  # Number of chitosan units
E_PDMS = 2.0  # Elastic modulus of PDMS

# === Configuration ===
base_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(base_dir, f"output_n{n_units}_E{E_PDMS:.1f}".replace('.', '_'))

# === Load state-separated hydration data ===
prefix = f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_')
free_file = os.path.join(output_dir, f"{prefix}_hydration_free.csv")
weak_file = os.path.join(output_dir, f"{prefix}_hydration_weak.csv")
tight_file = os.path.join(output_dir, f"{prefix}_hydration_tight.csv")

print("Looking for files:")
print("Free:", free_file)
print("Weak:", weak_file)
print("Tight:", tight_file)

hydration_file = os.path.join(output_dir, f"{prefix}_hydration_map.csv")
entropy_output_path = os.path.join(output_dir, "shannon_entropy.txt")

# === Load hydration data ===
hydration_values = np.loadtxt(hydration_file, delimiter=",", skiprows=1)
hydration_values = hydration_values[hydration_values > 0]

# === Normalization for probability density ===
p = hydration_values / np.sum(hydration_values)
epsilon = 1e-12  # to prevent log(0)
shannon_entropy = -np.sum(p * np.log(p + epsilon))

# === Save entropy value ===
with open(entropy_output_path, "w") as f:
    f.write(f"Continuous Shannon entropy (normalized): {shannon_entropy:.6f}\n")
    f.write(f"Total voxels (nonzero): {len(hydration_values)}\n")

print(f"✅ Saved Shannon entropy to: {entropy_output_path}")

# === Compute entropy and smearing for each state ===
def compute_entropy_and_smearing(filepath):
    try:
        values = np.loadtxt(filepath, delimiter=",", skiprows=1)
        values = values[values > 0]
        p_vals = values / np.sum(values)
        entropy = -np.sum(p_vals * np.log(p_vals + epsilon))
        smearing = np.var(values)
        return entropy, smearing
    except Exception as e:
        print(f"❌ Error reading {filepath}: {e}")
        return None, None

for label, path in zip(["free", "weak", "tight"], [free_file, weak_file, tight_file]):
    if not os.path.exists(path):
        print(f"❌ File missing: {path}")
    else:
        print(f"✅ File found: {path}")

# === Compute pressure for each state (mean hydration value) ===
state_pressures = {}
for label, path in zip(["free", "weak", "tight"], [free_file, weak_file, tight_file]):
    if os.path.exists(path):
        try:
            values = np.loadtxt(path, delimiter=",", skiprows=1)
            values = values[values > 0]
            pressure = np.mean(values)  # average energy density
            state_pressures[label] = pressure
        except Exception as e:
            print(f"❌ Error reading {path} for pressure: {e}")

# Normalize pressures
max_pressure = max(state_pressures.values()) if state_pressures else 1.0
for k in state_pressures:
    state_pressures[k] /= max_pressure

state_entropies = {}
state_temperatures = {}
for label, path in zip(["free", "weak", "tight"], [free_file, weak_file, tight_file]):
    if os.path.exists(path):
        ent, smear = compute_entropy_and_smearing(path)
        if ent is not None:
            state_entropies[label] = ent
            state_temperatures[label] = smear

# Normalize temperatures
max_temp = max(state_temperatures.values()) if state_temperatures else 1.0
for k in state_temperatures:
    state_temperatures[k] /= max_temp

# Save individual state entropies
with open(entropy_output_path, "a") as f:
    for state, ent in state_entropies.items():
        f.write(f"{state.title()} state Shannon entropy: {ent:.6f}\n")

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
font_dir = os.path.join(base_dir, "fonts")
possible_fonts = ["Helvetica.ttf", "Helvetica-Regular.ttf", "Helvetica.ttc"]
font_path = None
for fname in possible_fonts:
    test_path = os.path.join(font_dir, fname)
    if os.path.exists(test_path):
        font_path = test_path
        break
# already set above
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

from scipy.interpolate import interp1d

cycle_order = ["free", "weak", "tight"]
S_values = [state_entropies[s] for s in cycle_order]
T_values = [state_temperatures[s] for s in cycle_order]
labels = [s.title() for s in cycle_order]

color_map = {'Free': '#000000', 'Weak': '#4C6EA9', 'Tight': '#6EA8DC'}  # black to calm blue

plt.figure(figsize=(6, 4))
markers = {'Free': 'o', 'Weak': 's', 'Tight': 'D'}  # circle, square, diamond
for s, T, label in zip(S_values, T_values, labels):
    plt.scatter(s, T, label=label, marker=markers[label], color=color_map[label])
    plt.annotate(label, (s, T), textcoords="offset points", xytext=(5,5), ha='left', fontsize=9, color=color_map[label])

# Connect points with lines using segmented interpolation
for i in range(len(S_values) - 1):
    s_pair = [S_values[i], S_values[i+1]]
    T_pair = [T_values[i], T_values[i+1]]
    if s_pair[0] < s_pair[1]:
        log_s_interp = np.linspace(np.log10(s_pair[0]), np.log10(s_pair[1]), 100, endpoint=False)
    else:
        log_s_interp = np.linspace(np.log10(s_pair[1]), np.log10(s_pair[0]), 100, endpoint=False)[::-1]
    s_interp = 10 ** log_s_interp
    linear_interp = interp1d(s_pair, T_pair, kind='linear')
    T_interp = linear_interp(s_interp)
    plt.plot(s_interp, T_interp, linestyle='--', color='black', alpha=0.6)

# Draw arrows between points
for i in range(len(S_values) - 1):
    plt.annotate("", xy=(S_values[i+1], T_values[i+1]), xytext=(S_values[i], T_values[i]),
                 arrowprops=dict(arrowstyle="->", color='black', lw=1))

plt.xlabel("ln(Shannon Entropy)")
plt.ylabel("ln(Dimensionless Temperature)")
plt.title("Hydro-softening T-dS Relationship")
plt.legend()

# Proper log axis formatting with ticks
from matplotlib.ticker import LogLocator, NullFormatter

plt.xscale("log")
plt.yscale("log")
ax = plt.gca()
ax.invert_yaxis()
ax.xaxis.set_major_locator(LogLocator(base=10.0))
ax.yaxis.set_major_locator(LogLocator(base=10.0))
ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=10))
ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=10))
ax.xaxis.set_minor_formatter(NullFormatter())
ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_xticklabels([])
ax.set_yticklabels([])

plot_path = os.path.join(output_dir, "T_vs_S_plot.png")
plt.savefig(plot_path)
plt.show()
plt.close()

print(f"✅ Saved T vs S plot to: {plot_path}")

# === Plot Pressure vs Volume (P-v diagram) ===
plt.figure(figsize=(6, 4))
volume_values = [
    len(np.loadtxt(path, delimiter=",", skiprows=1)[np.loadtxt(path, delimiter=",", skiprows=1) > 0])
    if os.path.exists(path) else 1 for path in [free_file, weak_file, tight_file]
]
V_values = volume_values
P_values = [state_pressures[state] for state in cycle_order]

markers = {'Free': 'o', 'Weak': 's', 'Tight': 'D'}  # circle, square, diamond
for v, p, label in zip(V_values, P_values, labels):
    plt.scatter(v, p, label=label, marker=markers[label], color=color_map[label])
    plt.annotate(label, (v, p), textcoords="offset points", xytext=(5,5), ha='left', fontsize=9, color=color_map[label])

for i in range(len(V_values) - 1):
    plt.annotate("", xy=(V_values[i+1], P_values[i+1]), xytext=(V_values[i], P_values[i]),
                 arrowprops=dict(arrowstyle="->", color='black', lw=1))

plt.xlabel("Volume (number of hydrated voxels)")
plt.ylabel("Normalized Pressure (mean hydration value)")
plt.title("Hydro-softening P-v Relationship")
plt.legend()
pv_plot_path = os.path.join(output_dir, "P_vs_V_plot.png")
plt.savefig(pv_plot_path)
plt.show()
plt.close()

print(f"✅ Saved P vs V plot to: {pv_plot_path}")
print(f"✅ Saved T vs S plot to: {plot_path} (n={n_units}, E={E_PDMS})")
print(f"✅ Saved P vs V plot to: {pv_plot_path} (n={n_units}, E={E_PDMS})")