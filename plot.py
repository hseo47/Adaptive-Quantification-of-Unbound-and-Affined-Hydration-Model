import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from rdkit import Chem
from rdkit.Chem import AllChem
from matplotlib.lines import Line2D
import os
import matplotlib.font_manager as fm

base_dir = os.path.dirname(os.path.abspath(__file__))

font_path = os.path.join(base_dir, "fonts", "Computer Modern.ttf")
cm_font = fm.FontProperties(fname=font_path)
plt.rcParams['font.family'] = cm_font.get_name()
plt.rcParams.update({
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 16
})
import matplotlib.ticker as mticker
plt.rcParams['axes.unicode_minus'] = False

os.chdir(base_dir)  

repeat_unit = "YOUR SMILES STRING GOES HERE"
n_units = 1
E_PDMS = 0.5  
smiles = repeat_unit * n_units
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
conf = mol.GetConformer()

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

output_folder = f"output_n{n_units}_E{E_PDMS:.1f}".replace('.', '_')
hydration_filename = os.path.join(output_folder, f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_') + "_hydration_map.csv")
hydration_flat = np.loadtxt(hydration_filename, delimiter=",")

distrib_filename = os.path.join(base_dir, output_folder, f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_') + "_hyd_distrib.txt")
print(f"Looking for hydration distribution file at: {distrib_filename}")
print(f"Current working directory: {os.getcwd()}")
print(f"Files in script directory: {os.listdir(base_dir)}")
try:
    with open(distrib_filename) as f:
        lines = f.readlines()
        x_voxels = int([l for l in lines if "x_voxels" in l][0].split(":")[1])
        y_voxels = int([l for l in lines if "y_voxels" in l][0].split(":")[1])
        z_voxels = int([l for l in lines if "z_voxels" in l][0].split(":")[1])
except Exception as e:
    print(f"Failed to parse voxel shape from distribution file. Trying to infer from CSV...")
    x_voxels = len(x_range)
    y_voxels = len(y_range)
    z_voxels = len(z_range)
expected_size = x_voxels * y_voxels * z_voxels
if hydration_flat.shape[0] != expected_size:
    raise ValueError(f"Expected {expected_size} values but got {hydration_flat.shape[0]}")
print(f"Local hydration grid max: {hydration_flat.max()}")
hydration_grid = hydration_flat.reshape((x_voxels, y_voxels, z_voxels))

x, y, z = np.meshgrid(x_range, y_range, z_range, indexing="ij")
mask = hydration_grid > 0.1

ref_output_folder = f"output_n{n_units}_E2_0".replace('.', '_')
ref_filename = os.path.join(base_dir, ref_output_folder, f"n{n_units}_E2_0_hyd_distrib.txt")
try:
    with open(ref_filename) as f:
        lines = f.readlines()
        ref_max_line = [l.strip() for l in lines if "max_hydration" in l]
        if not ref_max_line:
            raise ValueError("max_hydration not found in reference file.")
        ref_max_hydration = float(ref_max_line[0].split(":")[1])
        print(f"Using global max hydration for color scale: {ref_max_hydration}")
except Exception as e:
    print(f"Failed to read global max hydration. Using local max. Reason: {e}")
    ref_max_hydration = hydration_grid.max()
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=25, azim=125)  
sc = ax.scatter(x[mask], y[mask], z[mask], c=hydration_grid[mask], cmap='Blues', alpha=0.3, s=20, vmin=0, vmax=ref_max_hydration)

cbar = fig.colorbar(sc, ax=ax, shrink=0.6)
cbar.set_label("Hydration Value", fontproperties=cm_font, fontsize=14)
cbar.ax.tick_params(labelsize=14)
for tick in cbar.ax.get_yticklabels():
    tick.set_fontproperties(cm_font)
    tick.set_fontsize(14)

ax.set_xlabel("X (Å)", fontproperties=cm_font, fontsize=16)
ax.set_ylabel("Y (Å)", fontproperties=cm_font, fontsize=16)
ax.set_zlabel("Z (Å)", fontproperties=cm_font, fontsize=16)

for label in (ax.get_xticklabels() + ax.get_yticklabels() + ax.get_zticklabels()):
    label.set_fontproperties(cm_font)
    label.set_fontsize(14)

ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))
ax.zaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))

pdms_plane_z = atom_coords[:, 2].min() - 3.0
xlim = ax.get_xlim()
ylim = ax.get_ylim()
xx, yy = np.meshgrid(np.linspace(*xlim, num=2), np.linspace(*ylim, num=2))
zz = np.full_like(xx, pdms_plane_z)
ax.plot_surface(xx, yy, zz, alpha=0.2, color='gray', linewidth=0, antialiased=False)

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
leg = ax.legend(handles=legend_elements, loc='upper right')
for text in leg.get_texts():
    text.set_fontproperties(cm_font)
    text.set_fontsize(14)

leg.get_frame().set_linewidth(0)
leg.get_frame().set_facecolor('none')

for i in bonded_oh_indices:
    x_val, y_val, z_val = atom_coords[i]
    ax.plot([x_val, x_val], [y_val, y_val], [z_val, pdms_plane_z],
            color='blue', linestyle='--', linewidth=1, alpha=0.6)

plt.tight_layout()
plot_filename = os.path.join(output_folder, f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_') + "_hydration_plot.png")
print("Reached the savefig line.")
print(f"Saving to → {os.path.abspath(plot_filename)}")
print(f"Figure has {len(ax.collections)} collections and {len(ax.texts)} text annotations.")
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
print(f"Saved plot to: {plot_filename}")
plt.show()
# plt.close()
