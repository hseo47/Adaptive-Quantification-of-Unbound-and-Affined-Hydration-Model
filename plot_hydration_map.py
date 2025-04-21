import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from rdkit import Chem
from rdkit.Chem import AllChem
from matplotlib.lines import Line2D
import os

base_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(base_dir)  # Make sure we are operating in the script's directory

# === Load molecule and generate coordinates ===
repeat_unit = "NC(CO)C(CO)"
n_units = 
E_PDMS =   # MPa — must match the value used in the generator script
smiles = repeat_unit * n_units
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
conf = mol.GetConformer()

# === Identify primary hydroxyls and bonded subset ===
primary_oh_smarts = Chem.MolFromSmarts("[CH2][OH]")
matches = mol.GetSubstructMatches(primary_oh_smarts)
primary_oh_indices = [match[1] for match in matches]

np.random.seed(42)
bonded_oh_indices = set(np.random.choice(primary_oh_indices, size=int(0.4 * len(primary_oh_indices)), replace=False))

atom_coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
atom_symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(mol.GetNumAtoms())]

min_coord = atom_coords.min(axis=0) - 2.0
max_coord = atom_coords.max(axis=0) + 2.0
grid_size = 1.0

x_range = np.arange(min_coord[0], max_coord[0], grid_size)
y_range = np.arange(min_coord[1], max_coord[1], grid_size)
z_range = np.arange(min_coord[2], max_coord[2], grid_size)

# === Load hydration field ===
hydration_filename = f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_') + "_hydration_map.csv"
hydration_flat = np.loadtxt(hydration_filename, delimiter=",")
# Load expected voxel shape from corresponding distribution file
distrib_filename = os.path.join(base_dir, f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_') + "_hyd_distrib.txt")
print(f"🧭 Looking for hydration distribution file at: {distrib_filename}")
print(f"📂 Current working directory: {os.getcwd()}")
print(f"📁 Files in script directory: {os.listdir(base_dir)}")
try:
    with open(distrib_filename) as f:
        lines = f.readlines()
        x_voxels = int([l for l in lines if "x_voxels" in l][0].split(":")[1])
        y_voxels = int([l for l in lines if "y_voxels" in l][0].split(":")[1])
        z_voxels = int([l for l in lines if "z_voxels" in l][0].split(":")[1])
except Exception as e:
    print(f"⚠️ Failed to parse voxel shape from distribution file. Trying to infer from CSV...")
    x_voxels = len(x_range)
    y_voxels = len(y_range)
    z_voxels = len(z_range)
expected_size = x_voxels * y_voxels * z_voxels
if hydration_flat.shape[0] != expected_size:
    raise ValueError(f"Expected {expected_size} values but got {hydration_flat.shape[0]}")
print(f"📈 Local hydration grid max: {hydration_flat.max()}")
hydration_grid = hydration_flat.reshape((x_voxels, y_voxels, z_voxels))

# === Real-space voxel positions ===
x, y, z = np.meshgrid(x_range, y_range, z_range, indexing="ij")
mask = hydration_grid > 0.1

# === Plot ===
# === Reference maximum hydration from E = 2.0 simulation ===
ref_filename = os.path.join(base_dir, "n20_E2_0_hyd_distrib.txt")
try:
    with open(ref_filename) as f:
        lines = f.readlines()
        ref_max_line = [l.strip() for l in lines if "max_hydration" in l]
        if not ref_max_line:
            raise ValueError("max_hydration not found in reference file.")
        ref_max_hydration = float(ref_max_line[0].split(":")[1])
        print(f"📊 Using global max hydration for color scale: {ref_max_hydration}")
except Exception as e:
    print(f"⚠️ Failed to read global max hydration. Using local max. Reason: {e}")
    ref_max_hydration = hydration_grid.max()
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=25, azim=125)  # Set consistent 3D view angle
sc = ax.scatter(x[mask], y[mask], z[mask], c=hydration_grid[mask], cmap='Blues', alpha=0.3, s=20, vmin=0, vmax=ref_max_hydration)

fig.colorbar(sc, ax=ax, shrink=0.6, label="Hydration Value")
ax.set_title(f"n = {n_units}; E = {E_PDMS} MPa")
ax.set_xlabel("X (Å)")
ax.set_ylabel("Y (Å)")
ax.set_zlabel("Z (Å)")

# PDMS substrate z-level
pdms_plane_z = atom_coords[:, 2].min() - 3.0
xlim = ax.get_xlim()
ylim = ax.get_ylim()
xx, yy = np.meshgrid(np.linspace(*xlim, num=2), np.linspace(*ylim, num=2))
zz = np.full_like(xx, pdms_plane_z)
ax.plot_surface(xx, yy, zz, alpha=0.2, color='gray', linewidth=0, antialiased=False)

# === Plot molecule atoms ===
mol_center = np.mean(atom_coords, axis=0)
grid_center = (max_coord + min_coord) / 2
shift = grid_center - mol_center
atom_coords += shift

for i, (x_val, y_val, z_val) in enumerate(atom_coords):
    if i in bonded_oh_indices:
        color = 'blue'
    elif atom_symbols[i] in ['O', 'N']:
        color = 'red'
    else:
        color = 'black'
    ax.scatter(x_val, y_val, z_val, color=color, s=40, edgecolors='white', linewidths=0.5)

legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Polar Regimes', markerfacecolor='red', markersize=6),
    Line2D([0], [0], marker='o', color='w', label='Non-Polar Regimes', markerfacecolor='black', markersize=6),
    Line2D([0], [0], marker='o', color='w', label='Bonded to Substrate', markerfacecolor='blue', markersize=6)
]
ax.legend(handles=legend_elements, loc='upper right')

for i in bonded_oh_indices:
    x_val, y_val, z_val = atom_coords[i]
    ax.plot([x_val, x_val], [y_val, y_val], [z_val, pdms_plane_z],
            color='blue', linestyle='--', linewidth=1, alpha=0.6)
plt.rcParams['font.family'] = 'Helvetica'
plt.tight_layout()
plot_filename = f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_') + "_hydration_plot.png"
print("🔥 Reached the savefig line.")
print(f"Saving to → {os.path.abspath(plot_filename)}")
print(f"Figure has {len(ax.collections)} collections and {len(ax.texts)} text annotations.")
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
print(f"🖼 Saved plot to: {plot_filename}")
plt.show()
# plt.close()
