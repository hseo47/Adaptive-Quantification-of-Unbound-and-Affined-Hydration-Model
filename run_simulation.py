import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.special import expi
import os

n_units = 10  # Adjust this to control chain length
E_PDMS = 2  # Adjust this to control substrate modulus
base_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(base_dir, f"output_n{n_units}_E{E_PDMS:.1f}".replace('.', '_'))
os.makedirs(output_dir, exist_ok=True)

# === Molecule setup ===
repeat_unit = "NC(CO)C(CO)"
smiles = repeat_unit * n_units
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
conf = mol.GetConformer()

# === Functional group hydration weights ===
hydration_weights = {
    "O": 1.0,     # e.g. hydroxyl
    "N": 1.4,     # protonated amines (NH3+)
    "C": 0.3,     # aliphatic C
    "H": 0.1      # small contribution
}

# === Identify primary hydroxyls (C-O-H where C is primary carbon) ===
primary_oh_smarts = Chem.MolFromSmarts("[CH2][OH]")
matches = mol.GetSubstructMatches(primary_oh_smarts)
primary_oh_indices = [match[1] for match in matches]

# Randomly select 40% of these to simulate PDMS bonding
np.random.seed(42)
bonded_oh_indices = set(np.random.choice(primary_oh_indices, size=int(0.4 * len(primary_oh_indices)), replace=False))

# === Bounding box calculation ===
atom_coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
min_coord = atom_coords.min(axis=0) - 2.0
max_coord = atom_coords.max(axis=0) + 2.0

# === Voxel grid setup ===
grid_size = 1.0  # Å per voxel
x_range = np.arange(min_coord[0], max_coord[0], grid_size)
y_range = np.arange(min_coord[1], max_coord[1], grid_size)
z_range = np.arange(min_coord[2], max_coord[2], grid_size)

voxel_grid = np.zeros((len(x_range), len(y_range), len(z_range)))

# Add per-class contributions dictionary
class_contributions = {
    "N": np.zeros_like(voxel_grid),
    "O": np.zeros_like(voxel_grid),
    "C": np.zeros_like(voxel_grid),
    "H": np.zeros_like(voxel_grid)
}

# === Hydration model (Gaussian kernel) ===
sigma = 1.5  # Å, width of Gaussian kernel
sigma_theta = 25  # degrees, angular width of directional hydration
theta_pref = {
    "O": 110,   # OH hydrogen bonding angle
    "N": 109    # NH3+ hydrogen bonding lobe (idealized tetrahedral)
}

# Assign preferred direction vectors per atom
preferred_directions = {}
for atom in mol.GetAtoms():
    idx = atom.GetIdx()
    symbol = atom.GetSymbol()
    if symbol in ["O", "N"]:
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
        if neighbors:
            vecs = [atom_coords[idx] - atom_coords[nbr_idx] for nbr_idx in neighbors]
            avg_vec = np.mean(vecs, axis=0)
            norm_vec = avg_vec / np.linalg.norm(avg_vec)
            preferred_directions[idx] = norm_vec
 
# === PDMS modulus-driven shielding factor ===
hydration_boost = max(expi(E_PDMS), 0)
for i, x in enumerate(x_range):
    for j, y in enumerate(y_range):
        for k, z in enumerate(z_range):
            point = np.array([x, y, z])
            total = 0
            for atom_idx, atom_pos in enumerate(atom_coords):
                if atom_idx in bonded_oh_indices:
                    continue  # Skip bonded OHs
                dist = np.linalg.norm(point - atom_pos)
                if dist < 0.5:
                    dist = 0.5
                atom_symbol = mol.GetAtomWithIdx(atom_idx).GetSymbol()
                w = hydration_weights.get(atom_symbol, 0.2)
                # Steric shielding attenuation from bonded OHs
                steric_attenuation = 1.0
                for shield_idx in bonded_oh_indices:
                    shield_pos = atom_coords[shield_idx]
                    shield_dist = np.linalg.norm(point - shield_pos)
                    steric_attenuation *= 1 + hydration_boost * np.exp(-shield_dist**2 / (2 * (2.5**2)))  # hydration-enhancing effect
                direction_factor = 1.0
                if atom_idx in preferred_directions:
                    atom_to_point = point - atom_pos
                    atom_to_point /= np.linalg.norm(atom_to_point)
                    angle_cos = np.clip(np.dot(atom_to_point, preferred_directions[atom_idx]), -1.0, 1.0)
                    angle_deg = np.degrees(np.arccos(angle_cos))
                    theta0 = theta_pref.get(atom_symbol, 90)
                    direction_factor = np.exp(-(angle_deg - theta0)**2 / (2 * sigma_theta**2))
                contrib = steric_attenuation * w * np.exp(-dist**2 / (2 * sigma**2)) * direction_factor
                if atom_symbol in ["N", "O", "C", "H"]:
                    class_contributions[atom_symbol][i, j, k] += contrib
                total += contrib
                voxel_grid[i, j, k] = total
                        
# === Flatten and save ===
flat = voxel_grid.flatten()
filename_suffix = f"n{n_units}_E{E_PDMS:.1f}".replace('.', '_')
output_path = os.path.join(output_dir, f"{filename_suffix}_hydration_map.csv")
# Save hydration distribution with unique filename
hydration_mean = voxel_grid.mean()
hydration_std = voxel_grid.std()
distrib_path = os.path.join(output_dir, f"{filename_suffix}_hyd_distrib.txt")
with open(distrib_path, "w") as f:
    f.write(f"mean: {hydration_mean:.4f}\n")
    f.write(f"std: {hydration_std:.4f}\n")
    f.write(f"n_units: {n_units}\n")
    f.write(f"E_PDMS: {E_PDMS:.2f} MPa\n")
    f.write(f"total_voxels: {voxel_grid.size}\n")
    np.savetxt(output_path, flat, delimiter=",")
    
for atom_type, grid in class_contributions.items():
    atom_output_path = os.path.join(output_dir, f"{filename_suffix}_hydration_{atom_type}.csv")
    np.savetxt(atom_output_path, grid.flatten(), delimiter=",")
    print(f"🧪 Saved per-atom hydration: {atom_output_path}")
    
# === Voxel classification based on hydration strength ===
hydration_classes = np.zeros_like(voxel_grid, dtype=int)
v90 = np.percentile(voxel_grid, 90)
v40 = np.percentile(voxel_grid, 40)

hydration_classes[voxel_grid >= v90] = 2  # Tightly bound
hydration_classes[(voxel_grid >= v40) & (voxel_grid < v90)] = 1  # Weakly bound
hydration_classes[voxel_grid < v40] = 0  # Free/unbound

# === Save hydration class map for voxel-by-voxel analysis ===
class_map_output = os.path.join(output_dir, f"{filename_suffix}_hydration_class_map.csv")
np.savetxt(class_map_output, hydration_classes.flatten(), delimiter=",", fmt='%d')
print(f"📊 Hydration class map saved to: {class_map_output}")

total_voxels = hydration_classes.size
counts = [(hydration_classes == i).sum() for i in range(3)]
labels = ["free", "weakly bound", "tightly bound"]
for i, count in enumerate(counts):
    percent = 100 * count / total_voxels
    print(f"{labels[i].capitalize()} voxels: {percent:.2f}%")

with open(distrib_path, "a") as f:
    f.write(f"x_voxels: {len(x_range)}\n")
    f.write(f"y_voxels: {len(y_range)}\n")
    f.write(f"z_voxels: {len(z_range)}\n")
    f.write("voxel_values:\n")
    np.savetxt(f, flat, delimiter=",", fmt="%.6f")
    f.write(f"max_hydration: {voxel_grid.max():.6f}\n")
print(f"✅ Hydration map saved to: {output_path}")
print(f"📊 Hydration class map saved to: {class_map_output}")
print(f"📁 All files saved in: {output_dir}")

# === Export per-class voxel values for streamlined entropy/TdS analysis ===
free_voxels = voxel_grid[hydration_classes == 0].flatten()
weak_voxels = voxel_grid[hydration_classes == 1].flatten()
tight_voxels = voxel_grid[hydration_classes == 2].flatten()

# Ensure all files include the correct filename_suffix and are saved to output_dir
free_path = os.path.join(output_dir, f"{filename_suffix}_hydration_free.csv")
weak_path = os.path.join(output_dir, f"{filename_suffix}_hydration_weak.csv")
tight_path = os.path.join(output_dir, f"{filename_suffix}_hydration_tight.csv")

np.savetxt(free_path, free_voxels, delimiter=",")
np.savetxt(weak_path, weak_voxels, delimiter=",")
np.savetxt(tight_path, tight_voxels, delimiter=",")

print(f"✅ Saved free voxel hydration to: {free_path}")
print(f"✅ Saved weak voxel hydration to: {weak_path}")
print(f"✅ Saved tight voxel hydration to: {tight_path}")
