import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D

n_units    = 1   # REPEATING UNIT FOR YOUR POLYMER OF INTEREST
E_PDMS     = 0.5 # SUBSTRATE, REPLACE AS NECESSARY  
REF_E_PDMS = 2.0 # SUBSTRATE, REPLACE AS NECESSARY (FOR GLOBAL MAXIMA)  

def _parse_line(text, key):
    for line in text.splitlines():
        if key in line:
            return line.split(":")[1].strip()
    raise KeyError(f"{key} not found in distrib file")

def load_hydration_map(n_units, e_value, base_dir):
    out_dir = f"output_n{n_units}_E{e_value:.1f}".replace('.', '_')
    csv_file = os.path.join(base_dir, out_dir,
                            f"n{n_units}_E{e_value:.1f}".replace('.', '_')
                            + "_hydration_map.csv")
    flat = np.loadtxt(csv_file, delimiter=",")
    distrib_txt = os.path.join(base_dir, out_dir,
                              f"n{n_units}_E{e_value:.1f}".replace('.', '_')
                              + "_hyd_distrib.txt")
    text = open(distrib_txt).read()
    xv = int(_parse_line(text, "x_voxels"))
    yv = int(_parse_line(text, "y_voxels"))
    zv = int(_parse_line(text, "z_voxels"))
    if flat.size != xv*yv*zv:
        raise ValueError(f"Expected {xv*yv*zv} values but got {flat.size}")
    return flat.reshape((xv, yv, zv))

def compute_global_thresholds(base_dir):
    ref_grid = load_hydration_map(n_units, REF_E_PDMS, base_dir)
    gmax = ref_grid.max()

    return {
        0: 0.2 * gmax,   
        1: 0.4 * gmax    
    }

def reclassify_grid(grid, thresholds):
    def cls(h):
        if   h < thresholds[0]: return 0
        elif h < thresholds[1]: return 1
        else:                   return 2
    return np.vectorize(cls)(grid)

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
    'figure.titlesize': 16,
    'axes.unicode_minus': False
})

filename_suffix = f"n1_E{E_PDMS:.1f}".replace(".", "_")
output_folder = f"output_n1_E{E_PDMS:.1f}".replace(".", "_")

thresholds = compute_global_thresholds(base_dir)

local_grid = load_hydration_map(n_units, E_PDMS, base_dir)
x_voxels, y_voxels, z_voxels = local_grid.shape

try:
    distrib_txt = os.path.join(base_dir, 
        f"output_n{n_units}_E{E_PDMS:.1f}".replace('.', '_'),
        f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_') + "_hyd_distrib.txt"
    )
    text = open(distrib_txt).read()
    x_min = float(_parse_line(text, "x_min"))
    x_max = float(_parse_line(text, "x_max"))
    y_min = float(_parse_line(text, "y_min"))
    y_max = float(_parse_line(text, "y_max"))
    z_min = float(_parse_line(text, "z_min"))
    z_max = float(_parse_line(text, "z_max"))
    x_range = np.linspace(x_min, x_max, x_voxels)
    y_range = np.linspace(y_min, y_max, y_voxels)
    z_range = np.linspace(z_min, z_max, z_voxels)
except (KeyError, FileNotFoundError):
    print("Real-space bounds not found; using grid indices.")
    x_range = np.arange(x_voxels)
    y_range = np.arange(y_voxels)
    z_range = np.arange(z_voxels)

class_grid = reclassify_grid(local_grid, thresholds)

# === Create voxel grid ===
x, y, z = np.meshgrid(x_range, y_range, z_range, indexing="ij")
x = x.flatten()
y = y.flatten()
z = z.flatten()
c = class_grid.flatten()

# === Individual class figures ===
 # softer, lighter base colors for classes
base_colors = {
    0: (0.6, 0.8, 0.6),  # lighter muted green
    1: (0.9, 0.9, 0.7),  # lighter muted yellow
    2: (0.9, 0.7, 0.7)   # lighter muted red
}
labels = {0: "Unbound", 1: "Weakly Bound", 2: "Tightly Bound"}

class_alphas = {
    0: 0.1,   
    1: 0.2,   
    2: 0.3    
}

from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

chitosan_1mer_smiles = "YOUR DESIRED SMILES STRING"
mol = Chem.MolFromSmiles(chitosan_1mer_smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
conf = mol.GetConformer()
atom_coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
atom_symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(mol.GetNumAtoms())]

mol_center = atom_coords.mean(axis=0)
grid_center = np.array([x_range.mean(), y_range.mean(), z_range.mean()])
atom_coords += (grid_center - mol_center)

bonded_indices = []
for bond in mol.GetBonds():
    a1 = bond.GetBeginAtomIdx()
    a2 = bond.GetEndAtomIdx()
    if atom_symbols[a1] == "O" and atom_symbols[a2] == "H":
        bonded_indices.append(a1)
    elif atom_symbols[a2] == "O" and atom_symbols[a1] == "H":
        bonded_indices.append(a2)

substrate_z = 0

fig = plt.figure(figsize=(18, 10))
for row, E in enumerate([0.5, 2.0]):
    filename_suffix = f"n1_E{E:.1f}".replace(".", "_")
    output_folder = f"output_n1_E{E:.1f}".replace(".", "_")
    class_map_path = os.path.join(base_dir, output_folder, f"{filename_suffix}_hydration_class_map.csv")
    distrib_path = os.path.join(base_dir, output_folder, f"{filename_suffix}_hyd_distrib.txt")

    try:
        with open(distrib_path) as f:
            lines = f.readlines()
    except FileNotFoundError:
        alt_suffix = f"n{n_units}_E{E}".replace(".", "_")
        alt_path = os.path.join(base_dir, f"output_n{n_units}_E{E}".replace(".", "_"), f"{alt_suffix}_hyd_distrib.txt")
        with open(alt_path) as f:
            lines = f.readlines()
    x_voxels = int([l for l in lines if "x_voxels" in l][0].split(":")[1])
    y_voxels = int([l for l in lines if "y_voxels" in l][0].split(":")[1])
    z_voxels = int([l for l in lines if "z_voxels" in l][0].split(":")[1])

    hydration_filename = os.path.join(base_dir, output_folder, f"{filename_suffix}_hydration_map.csv")
    hydration_flat = np.loadtxt(hydration_filename, delimiter=",")
    hydration_grid = hydration_flat.reshape((x_voxels, y_voxels, z_voxels))
    kernel_mask = hydration_grid > 0.1
    kernel_mask_flat = kernel_mask.flatten()
    kernel_alpha = np.clip(hydration_grid[kernel_mask] / hydration_grid.max(), 0.1, 1.0)
    kernel_colors = cm.Blues(kernel_alpha)

    kernel_x, kernel_y, kernel_z = np.where(kernel_mask)

    thresholds = compute_global_thresholds(base_dir)
    class_grid = reclassify_grid(hydration_grid, thresholds)
    c = class_grid.flatten()

    try:
        distrib_txt = os.path.join(base_dir, 
            f"output_n{n_units}_E{E:.1f}".replace('.', '_'),
            f"n{n_units}_E{E:.1f}".replace('.', '_') + "_hyd_distrib.txt"
        )
        text = open(distrib_txt).read()
        x_min = float(_parse_line(text, "x_min"))
        x_max = float(_parse_line(text, "x_max"))
        y_min = float(_parse_line(text, "y_min"))
        y_max = float(_parse_line(text, "y_max"))
        z_min = float(_parse_line(text, "z_min"))
        z_max = float(_parse_line(text, "z_max"))
        x_local = np.linspace(x_min, x_max, x_voxels)
        y_local = np.linspace(y_min, y_max, y_voxels)
        z_local = np.linspace(z_min, z_max, z_voxels)
    except (KeyError, FileNotFoundError):
        print("Real-space bounds not found; using grid indices for local axes.")
        x_local = np.arange(x_voxels)
        y_local = np.arange(y_voxels)
        z_local = np.arange(z_voxels)
    x_local, y_local, z_local = np.meshgrid(x_local, y_local, z_local, indexing="ij")
    x_local = x_local.flatten()
    y_local = y_local.flatten()
    z_local = z_local.flatten()

    for col, value in enumerate([0, 1, 2]):
        ax = fig.add_subplot(2, 3, row * 3 + col + 1, projection='3d')
        ax.set_box_aspect([1,1,1])
        ax.set_facecolor('white')
        try:
            ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
            ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
            ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        except Exception:
            pass
        ax.view_init(elev=25, azim=120)
        mask = c == value

        grid_mask_flat = class_grid.flatten() == value
        ax.scatter(
            x_local[grid_mask_flat],
            y_local[grid_mask_flat],
            z_local[grid_mask_flat],
            color=base_colors[value],
            alpha=class_alphas[value],
            s=2
        )

        ax.scatter(
            x_local[mask],
            y_local[mask],
            z_local[mask],
            color=base_colors[value],
            s=12,
            alpha=class_alphas[value]
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.grid(False)
        ax.set_axis_off()

        for i, (x_val, y_val, z_val) in enumerate(atom_coords):
            color = 'red' if atom_symbols[i] in ['O', 'N'] else 'black'
            ax.scatter(x_val, y_val, z_val, color=color, s=20, edgecolors='white', linewidths=0.2)

        for i in bonded_indices:
            x_val, y_val, z_val = atom_coords[i]
            ax.plot([x_val, x_val], [y_val, y_val], [z_val, substrate_z],
                    color='blue', linestyle='--', linewidth=0.6, alpha=0.5)

        if col == 0:
            ax.text2D(-0.1, 0.5, f"E = {E:.2f} MPa", transform=ax.transAxes,
                      fontproperties=cm_font, fontsize=16, rotation=90, va='center')


from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Unbound', markerfacecolor=base_colors[0], markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Weakly Bound', markerfacecolor=base_colors[1], markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Tightly Bound', markerfacecolor=base_colors[2], markersize=8),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=3, 
           frameon=False, fontsize=14, bbox_to_anchor=(0.5, 0.02))

from matplotlib.patches import Rectangle

cb_ax = fig.add_axes([0.75, 0.05, 0.2, 0.05])  

class_labels = [
    f"< {thresholds[0]:.1f}",
    f"< {thresholds[1]:.1f}",
    f"≥ {thresholds[1]:.1f}"
]
for i, cls in enumerate([0, 1, 2]):
    rect = Rectangle((i * (1/3), 0), 1/3, 1, color=base_colors[cls], alpha=class_alphas[cls], transform=cb_ax.transAxes)
    cb_ax.add_patch(rect)

cb_ax.set_xticks([1/6, 3/6, 5/6])
cb_ax.set_xticklabels(class_labels, fontproperties=cm_font, fontsize=12)
cb_ax.set_yticks([])
cb_ax.set_xlim(0, 1)
cb_ax.set_ylim(0, 1)
cb_ax.set_xlabel("Hydration potential", fontproperties=cm_font, fontsize=12)
cb_ax.set_frame_on(False)

plt.subplots_adjust(wspace=0.3, hspace=0.3)
save_path = os.path.join(base_dir, "hydration_class_3D_grid.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Grid saved: {save_path}")
plt.close()
