import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import matplotlib.font_manager as fm

script_dir = os.path.dirname(os.path.abspath(__file__))
font_path = os.path.join(script_dir, "fonts", "Computer Modern.ttf")
font_prop = fm.FontProperties(fname=font_path)

dcd_path = os.path.join(script_dir, "trajectory.dcd")
pdb_path = os.path.join(script_dir, "merged_system.pdb")
traj = md.load(dcd_path, top=pdb_path)

print(f"Loaded trajectory with {traj.n_frames} frames, {traj.n_atoms} atoms.")
print(f"Simulation time: {traj.n_frames * traj.timestep:.1f} ps")

top = traj.topology
water_oxygen = top.select("water and name O")
chitosan_heavy = top.select("not water and not element H")

print("Finding water-chitosan pairs within 1.2 nm...")
pairs = []
for frame in traj:
    neighbors = md.compute_neighbors(frame, 1.2, chitosan_heavy)
    for atom_i, atom_list in zip(chitosan_heavy, neighbors):
        for atom_j in atom_list:
            if atom_j in water_oxygen:
                pairs.append([atom_i, atom_j])
pairs = np.unique(pairs, axis=0)

print(f"Found {len(pairs)} unique pairs for RDF.")
traj.unitcell_lengths = np.ones((traj.n_frames, 3)) * 10.0  # 10 nm box
traj.unitcell_angles = np.ones((traj.n_frames, 3)) * 90.0
rdf, r = md.compute_rdf(traj, pairs=pairs, r_range=(0.0, 1.2))
rdf = rdf.flatten()

plt.figure(figsize=(4, 4))
fm.fontManager.addfont(font_path)
plt.rcParams['font.family'] = font_prop.get_name()

plt.plot(r, rdf, color="black")
plt.grid(False)
plt.xlabel("Distance (nm)", fontsize=10)
plt.ylabel("g(r)", fontsize=10)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.title("Water-Chitosan RDF", fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(script_dir, "rdf_water_chitosan.png"), dpi=900)
print("RDF plot saved to rdf_water_chitosan.png")

print("Calculating RMSD...")
chitosan_atoms = top.select("not water")
rmsd = md.rmsd(traj, traj, frame=0, atom_indices=chitosan_atoms)

plt.figure(figsize=(4, 4))
fm.fontManager.addfont(font_path)
plt.rcParams['font.family'] = font_prop.get_name()
plt.plot(rmsd, color="black")
plt.grid(False)
plt.xlabel("Frame")
plt.ylabel("RMSD (nm)")
plt.title("Chitosan RMSD Over Time")
plt.tight_layout()
plt.savefig(os.path.join(script_dir, "rmsd_chitosan.png"), dpi=900)
print("RMSD plot saved to rmsd_chitosan.png")
