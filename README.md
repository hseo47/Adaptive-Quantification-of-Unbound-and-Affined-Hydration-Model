# hydratedPolymer: A Voxel-Based Model of Steric-Confined Water Binding in Polymeric Networks

> _Simulating steric shielding and directional hydration in chitosan films bonded to soft or stiff PDMS substrates (Nat. Commun. In Review)._

## 🧠 Project Overview
HydratedPolymer is a  simulation framework that models hydration potential in flexible hydrophilic polymer chains (e.g., chitosan) using a voxelized Gaussian kernel influenced by interfacial steric effects. The model maps 3D hydration fields with respect to:

- Local atomic polarity (NH3+, OH, CH)
- Substrate stiffness (via PDMS modulus in MPa)
- Spatial steric shielding from bonded sites
- Directional confinement effects

The model is grounded in voxel-based distance weighting, Gaussian interaction kernels, and physics-informed scaling (via exponential integrals).

---

## 🔬 Motivation
Hydro-softening (Nat. Commun., In Review) significantly changes the mechanical and biological responses (Matter, In Review) of polymers. However, few models capture how hydration evolves under spatially resolved steric effects from substrate interfacial stiffness (Driving mechanism of hydro-softening).

HydratedPolymer introduces a voxel-based strategy to:
- Simulate and visualize hydration potential in real-space
- Connect substrate stiffness to steric confinement
- Predict water-binding distributions in bonded polymer chains

---

## 📦 Features
- **SMILES-to-3D** structure generation via RDKit  
- **Hydration kernels** influenced by PDMS-induced shielding  
- **Voxel grid mapping** of hydration potentials  
- **Substrate Elastic Modulus scaling** with exponential integral modifiers  
- **Exportable CSV & 3D plots** for hydration distribution and molecular structure  

---

## 🚀 Quickstart

1. Clone this repo (private for now):
```bash
git clone https://github.com/your-org/HydrationSpace.git
```

2. Install dependencies
```bash
pip install -r requirements.txt
```
4. Run generator script (starting from max. elastic modulus)
```bash
python smiles_to_hydration_map.py
```

5. Visualize results
```bash
python plot_hydration_map.py
```

6. Post-processing generated data
```bash
avg_hydration.m
```
```bash
hydration_by_class.m
```
```bash
visualize_shielding_radius.py
```
```bash
radial_decay.py
```
```bash
radial_correlation.py
```
```bash
anisotropy.py
```
```bash
entropy.py
```

## 📁 File Structure

hydratedPolymer
smiles_to_hydration_map.py   
plot_hydration_map.py
avg_hydration.m
hydration_by_class.m
visualize_shielding_radius.py
radial_decay.py
radial_correlation.py
anisotropy.py
entropy.py
*.csv / *.txt / *.png     
README.md # you are here                   
