from rdkit import Chem
from rdkit.Chem import AllChem

# SMILES
smiles = "CN1CCN(CCCOc2cc3c(Nc4ccc(O)c(Cl)c4)ncnc3cc2CO)CC1"  
mol = Chem.MolFromSmiles(smiles)

# 1. Hydrogens
mol_h = Chem.AddHs(mol)

# 2. 3D conformation
AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())

# 3. Minimization
AllChem.MMFFOptimizeMolecule(mol_h)

# 4. Save to file .xyz
Chem.MolToXYZFile(mol_h, "3D_conformation.xyz")
print("¡XYZ file generated and ready for ORCA!")
