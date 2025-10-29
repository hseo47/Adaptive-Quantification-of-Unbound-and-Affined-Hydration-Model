# 💧 AQUA Hydration Model 💧

## Project Motivation and Overview

***Mechanical properties of polymeric materials are determined by the disproportionately small, tightly-bonded water particles*** ***near the polymer backbone.*** 

**AQUA** (Adaptive Quantification of Unbound and Affined) Hydration Model is a voxel-based computational framework designed to interpret hydration behavior in polymer chains under the influence of steric confinement from bonded interfaces. It models local hydration potential based on:

- Atomic polarity from chemical structures  
- Substrate-induced steric shielding effects 
- Spatial anisotropy and directional decay  
- Voxel-based hydration fields using Gaussian kernels  

For more information on **Hydro-Softening** and Water-Mediated Material Property Changes:
- Skin-Inspired Hydro-Softening Enables Flexible Chitosan Films (**Manuscript Linked to this Work**). 🧪
- Enhancing Bacterial Adhesion with Hydro-Softened Chitosan Films (*Published, ACS Macro Letters 2025*). 🦠

The data supporting the findings of this study will be made available with the paper and its supplementary materials. In compliance with institutional data management policies, the full codebase is maintained on Georgia Tech’s GitHub Enterprise instance. This is a stable public version is available at and is archived at Zenodo. This public version includes usage instructions and basic dependencies for replication.

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

---

## File Structure

| Script | Description |
|--------|-------------|
| `runAQUA.py` | Core script. |
| `plot.py` | Hydration map visualization. |
| `plotClass.py` | Hydration microstates visualization. |
