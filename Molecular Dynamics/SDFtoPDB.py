from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles

sdf_file = "File_Name_Here"

def convert_sdf_to_pdb(sdf_file, pdb_file):
    supplier = Chem.SDMolSupplier(sdf_file)
    mol = supplier[0]


    mol = Chem.AddHs(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=1, params=params)

    AllChem.UFFOptimizeMolecule(mol, confId=conf_ids[0])

    rdmolfiles.MolToPDBFile(mol, pdb_file)
    print(f"Conversion successful. PDB file saved as '{pdb_file}'.")

sdf_file_path = sdf_file
pdb_output_path = "Output_Name_Here"
convert_sdf_to_pdb(sdf_file_path, pdb_output_path)
