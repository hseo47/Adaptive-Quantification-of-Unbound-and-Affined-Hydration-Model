# hydratedPolymer: A Voxel-Based Model of Steric-Confined Water Binding in Polymeric Networks

> _Simulating steric shielding and directional hydration in chitosan films bonded to soft or stiff PDMS substrates._

## 🧠 Project Overview
HydrationSpace is a Python-based simulation framework that models hydration potential in flexible hydrophilic polymer chains (e.g., chitosan) using a voxelized Gaussian kernel influenced by interfacial steric effects. The model maps 3D hydration fields with respect to:

- Local atomic polarity (NH3+, OH, CH)
- Substrate stiffness (via PDMS modulus in MPa)
- Spatial steric shielding from bonded sites
- Directional confinement effects

The model is grounded in voxel-based distance weighting, Gaussian interaction kernels, and physics-informed scaling (via exponential integrals). Visualization tools render hydration fields alongside molecular structure using RDKit and matplotlib.

---

## 🔬 Motivation
The hydration state of biopolymer films like chitosan significantly impacts their mechanical, biological, and drug-release performance. However, few models capture how hydration evolves under spatially resolved steric effects from substrate bonding or interfacial stiffness.

HydrationSpace introduces a voxel-based strategy to:
- Simulate and visualize hydration potential in real-space
- Connect substrate stiffness to steric confinement
- Predict water-binding distributions in bonded polymer chains

---

## 📦 Features
- **SMILES-to-3D** structure generation via RDKit  
- **Hydration kernels** influenced by PDMS-induced shielding  
- **Voxel grid mapping** of hydration potentials  
- **E_PDMS scaling** with exponential integral modifiers  
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

## 📁 File Structure

hydratedPolymer
smiles_to_hydration_map.py   
plot_hydration_map.py    
*.csv / *.txt / *.png     
README.md # you are here                   
