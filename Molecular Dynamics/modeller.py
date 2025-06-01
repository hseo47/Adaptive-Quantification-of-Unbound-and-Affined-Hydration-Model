from openmm.app import Topology
from openmm.app import PDBFile, Modeller, Simulation, ForceField, PME, HBonds, DCDReporter
from openmm import LangevinIntegrator, Platform, System
from openmm import Vec3
from openmm.unit import nanometer, nanometers, picosecond, femtoseconds, kelvin, molar
from openff.toolkit.topology import Molecule
from openff.toolkit.topology import Topology as OFFTopology
from openff.toolkit.typing.engines.smirnoff import ForceField as OFFForceField
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
pdb_path = os.path.join(script_dir, "Chitosan.pdb")

print("Loading Chitosan PDB...")
pdb = PDBFile(pdb_path)
modeller = Modeller(pdb.topology, pdb.positions)

print("Creating pure TIP3P water box...")
empty_top = Topology()
empty_modeller = Modeller(empty_top, [])
tip3p_ff = ForceField("amber14/tip3p.xml")
empty_modeller.addSolvent(
    tip3p_ff,
    model='tip3p',
    boxSize=Vec3(4.0, 4.0, 4.0) * nanometer,
    neutralize=False,
    ionicStrength=0.0 * molar
)

water_positions = empty_modeller.positions
water_topology = empty_modeller.topology

solvated_pdb_path = os.path.join(script_dir, "solvated.pdb")
with open(solvated_pdb_path, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"Saved solvated PDB to {solvated_pdb_path}")

print("Assigning OpenFF parameters to Chitosan...")
sdf_path = os.path.join(script_dir, "Chitosan.sdf")
mol = Molecule.from_file(sdf_path, file_format="sdf", allow_undefined_stereo=True)
mol.assign_partial_charges("gasteiger")
openff_top = mol.to_topology()
openff_ff = OFFForceField("openff_unconstrained-2.0.0.offxml")
chitosan_system = openff_ff.create_openmm_system(openff_top, charge_from_molecules=[mol])

print("Generating water parameters...")
water_ff = ForceField("amber14/tip3p.xml")
water_system = water_ff.createSystem(water_topology, nonbondedMethod=PME, constraints=HBonds)

print("Merging Chitosan and Water topologies...")
modeller = Modeller(openff_top.to_openmm(), modeller.positions)
modeller.add(water_topology, water_positions)

chitosan_mol = mol
water_mol = Molecule.from_smiles("O")

print("Building full system with merged force fields...")
off_top = OFFTopology.from_openmm(modeller.topology, unique_molecules=[chitosan_mol, water_mol])
full_system = openff_ff.create_openmm_system(off_top, charge_from_molecules=[mol])
water_ff = ForceField("amber14/tip3p.xml")
water_system = water_ff.createSystem(water_topology, nonbondedMethod=PME, constraints=HBonds)

merged_pdb_path = os.path.join(script_dir, "merged_system.pdb")
with open(merged_pdb_path, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("Minimizing energy...")
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(modeller.topology, full_system, integrator)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()

print("Running short MD (10 ps)...")
dcd_path = os.path.join(script_dir, "trajectory.dcd")
simulation.reporters.append(DCDReporter(dcd_path, 100))
simulation.step(5000)

print("Simulation complete. Output saved to trajectory.dcd")
