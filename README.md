# 💧 AQUA Hydration Model 💧

## Project Motivation and Overview

***Mechanical properties of polymeric materials are determined by the disproportionately small, tightly-bonded water particles*** 💧 ***near the polymer backbone.*** 

**AQUA** (Adaptive Quantification of Unbound and Affined) is a voxel-based computational framework designed to interpret hydration behavior in polymer chains under the influence of steric confinement from bonded interfaces. It models local hydration potential based on:

- Atomic polarity from chemical structures  
- Substrate-induced steric shielding effects 
- Spatial anisotropy and directional decay  
- Voxel-based hydration fields using Gaussian kernels  

We have previously demonstrated softening of rigid materials with nanoscale confinement of water.

For more information on **Hydro-Softening** and Water-Mediated Material Property Changes:
- Skin-Inspired Hydro-Softening Enables Flexible Chitosan Films (*Nat. Commun., In Review*). 🧪
- Enhancing Bacterial Adhesion with Hydro-Softened Chitosan Films (*In Review*). 🦠
- Hydration Confinement as a Thermodynamic Driver of Hydro-Softening 
  via Adaptive-Quantification-of-Unbound-and-Affined (AQUA) Hydration Model
  **(Manuscript Linked to this Work)** 🖥️

---

## Requirements

To run the Hydro-Softening-Simulation, you’ll need the following:

### Python Environment
- `Python ≥ 3.8 (recommended: 3.11+)`

### Required Python Packages
- `numpy`
- `matplotlib`
- `scipy`
- `rdkit`
- `fonttools`

### OPTIONALLY: To Run MD as a Comparison to AQUA
- `openMM`
- `openFF-toolkit`
- `mdtraj`

---

## File Structure

| Script | Description |
|--------|-------------|
| `run_simulation.py` | Core script that takes a SMILES string and generates a voxelized 3D hydration map under steric constraints. |
| `plot.py` | Visualizes the output hydration map with atomic overlays and substrate planes. |
| `thermAnalysis.py` | Generates the T–dS diagram and Shannon entropy values based on different hydration states. |
| `thermDiag.py` | Plots histogram distributions of hydration values per state for diagnostic analysis. |
| `entropy.py` | Correlates substrate modulus and average hydration to entropy per voxel. |
| `radialDecay.py` | Extracts radial decay of hydration centered at max voxel and plots directionally. |
| `radialCorr.py` | Computes radial correlation function for hydration fluctuations around max voxel. |
| `shieldRad.py` | Calculates shielding radius effects from bonded hydroxyls to substrate. |
| `anisotropy.py` | Extracts anisotropic decay along X, Y, and Z axes relative to max hydration site. |
| `avgHydration.m` | MATLAB script for computing average hydration per class (if class labels exist). |
| `hydrationClass.m` | MATLAB script for class-based hydration assignment. |
| `fonts/` | Contains Helvetica typeface for consistent visuals. |
| `output_n##_E#_#` | Output folder for a given number of repeat units (n) and PDMS modulus (E). |

---

## Operation

### 1. Clone the Repository
```bash
git clone https://github.com/your-org/Hydro-Softening-Simulation.git
cd Hydro-Softening-Simulation
```

### 2. Install Dependencies
```bash
pip install -r requirements.txt
```

### 3. Run the Hydration Simulation
Generate the 3D hydration map for a polymer with `n_units` repeat units and PDMS stiffness `E_PDMS`:
```bash
python run_simulation.py
```

### 4. Visualize Hydration Field
Render a 3D hydration plot overlaid on the polymer structure:
```bash
python plot.py
```

### 5. Thermodynamic Analysis
Generate T–dS plots and entropy diagnostics for different hydration states:
```bash
python thermAnalysis.py
python thermDiag.py
```

### 6. Post-Processing Scripts
Analyze hydration distribution using the following scripts:
```bash
python entropy.py         # Shannon entropy vs modulus
python anisotropy.py      # Directional hydration decay
python radialDecay.py     # Radial decay from max hydration site
python radialCorr.py      # Radial correlation function
python shieldRad.py       # Steric shielding radius visualization
```
