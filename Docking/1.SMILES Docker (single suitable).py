import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import vina
import meeko

def read_pdb(pdb_file):
    # pdb_file = open(pdb_file, "r").read()
    # prot = Chem.MolFromMolBlock(pdb_file, sanitize=False, removeHs=False)
    prot = Chem.MolFromPDBFile(
        pdb_file,
        removeHs=False,
        sanitize=False,
        proximityBonding=True)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SETAROMATICITY)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SETCONJUGATION)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SETHYBRIDIZATION)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SYMMRINGS)
    # Chem.SanitizeMol(prot, Chem.SANITIZE_PROPERTIES)
    Chem.SanitizeMol(prot, Chem.SANITIZE_CLEANUP)
    # prot = Chem.AddHs(prot, addCoords=True)
    return prot
  
def parse_poses(poses):
    pmol = meeko.PDBQTMolecule(poses)
    return meeko.RDKitMolCreate.from_pdbqt_mol(pmol)[0]

protein_file_pdbqt = "./data/Receptor files/7vhy-prepped.pdbqt"
protein_file_pdb = "./data/Receptor files/7vhy-prepped.pdb" # just for complex depiction
ligs = {
  "6QI": "COc1ccc(cc1)C2(CCCC2)C(=O)N3CCOC[C@@H](C3)Cc4cccc5c4cn[nH]5"
#  "QRY": "C1=CC=C(C=C1)[C@H](C(=O)C2=CNC3=CC=CC=C32)NCCC4=CC=C(C=C4)S(=O)(=O)N"
#  "YJY": "Cc1noc(C)c1c1cc2nc(C3CCCC(=O)N3c3cc(cc(c3)C(F)(F)F)C(F)(F)F)n(C3CCC(OC)CC3)c2cc1"
#  "99E": "CCC1=C(NC(=C1C(=O)C)C)C(=O)NCC2=C(C=CC3=C2C(=CC=C3)O)O"
}

for lig_name, lig in ligs.items():
  lig = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(lig), canonical=True, isomericSmiles=True))
  lig = Chem.AddHs(lig)
  Chem.AllChem.EmbedMolecule(lig, randomSeed=42)
  meeko_prep = meeko.MoleculePreparation()
  meeko_prep.prepare(lig)
  lig_pdbqt = meeko_prep.write_pdbqt_string()
  v = vina.Vina(sf_name='vina', seed=42)
  v.set_receptor(protein_file_pdbqt)
  v.set_ligand_from_string(lig_pdbqt)
  # previous_ligand = next(Chem.SDMolSupplier('5t1a_C_VT5.sdf')) # coords from 5t1a
  # centroid = Chem.rdMolTransforms.ComputeCentroid(previous_ligand.GetConformer())
  # v.compute_vina_maps(center=[centroid.x, centroid.y,
  #                             centroid.z], box_size=[30, 30, 30])
  v.compute_vina_maps(
    # this box is used in the drugex project
    center=[-12.2, -16.1, 13.3], 
    box_size=[16.5, 20.6, 16.7]
  )
  v.dock(exhaustiveness=64)
  output_pdbqt = v.poses()
  with open(f"{lig_name}_poses.pdbqt", "w") as pout:
      pout.write(output_pdbqt)
  poses = parse_poses(output_pdbqt)
  f = Chem.SDWriter(f'{lig_name}_poses.sdf')
  for pose_id in range(poses.GetNumConformers()):
      # conf = poses.GetConformer(pose_id)
      # conf = Chem.Mol(conf, 0)
      f.write(poses, confId=pose_id)
  f.close()
  # combine into a complex
  for idx, pose in enumerate(Chem.SDMolSupplier(f'{lig_name}_poses.sdf')):
    protein = read_pdb(protein_file_pdb)
    complex = Chem.CombineMols(protein, pose)
    Chem.MolToPDBFile(complex, lig_name + f"_complex_{idx}.pdb", flavor=4)
