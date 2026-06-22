import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen
from rdkit.Chem import rdMolDescriptors

# =========================================================================
# 1. POST-DOCKING SCAFFOLDS (RENAMED WITH EXPLICIT SUBSTITUENTS)
# =========================================================================
SCAFFOLDS = {
    # Piperazine Derivatives (16)
    "Pip_OH_F": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(O)=C(F)C=C4)=C1",
    "Pip_F_F": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(F)=C(F)C=C4)=C1",
    "Pip_Cl_F": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(Cl)=C(F)C=C4)=C1",
    "Pip_Cl_OH": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(Cl)=C(O)C=C4)=C1",
    "Pip_F_Cl": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(F)=C(Cl)C=C4)=C1",
    "Pip_Cl_H": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(Cl)=CC=C4)=C1",
    "Pip_Br_F": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(Br)=C(F)C=C4)=C1",
    "Pip_OH_Cl": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(O)=C(Cl)C=C4)=C1",
    "Pip_F_OH": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(F)=C(O)C=C4)=C1",
    "Pip_OH_Br": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(O)=C(Br)C=C4)=C1",
    "Pip_OH_H": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(O)=CC=C4)=C1",
    "Pip_F_Br": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(F)=C(Br)C=C4)=C1",
    "Pip_F_H": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(F)=CC=C4)=C1",
    "Pip_OH_OH": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(O)=C(O)C=C4)=C1",
    "Pip_Cl_Cl": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC(Cl)=C(Cl)C=C4)=C1",
    "Pip_H_Cl": "[*]C1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC=C(Cl)C=C4)=C1",
    
    # Morpholine Derivatives (6)
    "Mor_Cl_F": "[*]C1=C(OCCCN2CCOCC2)C=C3C(N=CN=C3NC4=CC(Cl)=C(F)C=C4)=C1",
    "Mor_H_F": "[*]C1=C(OCCCN2CCOCC2)C=C3C(N=CN=C3NC4=CC=C(F)C=C4)=C1",
    "Mor_OH_F": "[*]C1=C(OCCCN2CCOCC2)C=C3C(N=CN=C3NC4=CC(O)=C(F)C=C4)=C1",
    "Mor_Cl_H": "[*]C1=C(OCCCN2CCOCC2)C=C3C(N=CN=C3NC4=CC(Cl)=CC=C4)=C1",
    "Mor_F_F": "[*]C1=C(OCCCN2CCOCC2)C=C3C(N=CN=C3NC4=CC(F)=C(F)C=C4)=C1",
    "Mor_F_Br": "[*]C1=C(OCCCN2CCOCC2)C=C3C(N=CN=C3NC4=CC(F)=C(Br)C=C4)=C1",

    # N-Methylpiperazine Derivatives (8)
    "NMP_OH_F": "[*]C1=C(OCCCN2CCN(C)CC2)C=C3C(N=CN=C3NC4=CC(O)=C(F)C=C4)=C1",
    "NMP_Cl_F": "[*]C1=C(OCCCN2CCN(C)CC2)C=C3C(N=CN=C3NC4=CC(Cl)=C(F)C=C4)=C1",
    "NMP_F_F": "[*]C1=C(OCCCN2CCN(C)CC2)C=C3C(N=CN=C3NC4=CC(F)=C(F)C=C4)=C1",
    "NMP_F_Cl": "[*]C1=C(OCCCN2CCN(C)CC2)C=C3C(N=CN=C3NC4=CC(F)=C(Cl)C=C4)=C1",
    "NMP_F_OH": "[*]C1=C(OCCCN2CCN(C)CC2)C=C3C(N=CN=C3NC4=CC(F)=C(O)C=C4)=C1",
    "NMP_OH_H": "[*]C1=C(OCCCN2CCN(C)CC2)C=C3C(N=CN=C3NC4=CC(O)=CC=C4)=C1",
    "NMP_OH_Cl": "[*]C1=C(OCCCN2CCN(C)CC2)C=C3C(N=CN=C3NC4=CC(O)=C(Cl)C=C4)=C1",
    "NMP_H_F": "[*]C1=C(OCCCN2CCN(C)CC2)C=C3C(N=CN=C3NC4=CC=C(F)C=C4)=C1"
}

# =========================================================================
# 2. TARGET SUBSTITUENTS
# =========================================================================
SUBSTITUENTS = {
    "OCH3": "OC",
    "OH": "O",
    "NH2": "N",
    "NO2": "[N+](=O)[O-]"
}

# =========================================================================
# 3. FRAGMENT COUPLING ENGINE (SINGLE DUMMY ATOM)
# =========================================================================
def generate_molecule(base_smiles, sub_smiles):
    """
    Couples a substituent to the scaffold by finding the dummy atom [*].
    Uses the Chem.CombineMols and Chem.RWMol template.
    """
    mol_base = Chem.MolFromSmiles(base_smiles)
    mol_sub = Chem.MolFromSmiles(sub_smiles)

    if mol_base is None or mol_sub is None:
        return None

    # Combine molecules into an editable object
    combo = Chem.CombineMols(mol_base, mol_sub)
    rw = Chem.RWMol(combo)

    dummy_idx = None
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() == 0:
            dummy_idx = atom.GetIdx()
            break

    if dummy_idx is None:
        return None

    # Find the scaffold atom connected to the dummy atom
    dummy_neighbors = [n.GetIdx() for n in rw.GetAtomWithIdx(dummy_idx).GetNeighbors()]

    if len(dummy_neighbors) != 1:
        return None

    scaffold_anchor_idx = dummy_neighbors[0]
    
    # The substituent starts precisely at the index equal to the base molecule's atom count
    sub_start_idx = mol_base.GetNumAtoms()

    # Form the bond and drop the dummy atom
    rw.AddBond(scaffold_anchor_idx, sub_start_idx, Chem.BondType.SINGLE)
    rw.RemoveAtom(dummy_idx)

    final_mol = rw.GetMol()

    try:
        Chem.SanitizeMol(final_mol)
        return final_mol
    except:
        return None

# =========================================================================
# 4. PHYSICOCHEMICAL & ADMET FILTERS
# =========================================================================
def passes_admet_filters(mol):
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot = Lipinski.NumRotatableBonds(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)

    # Lipinski Rule of 5 Constraints
    lipinski_ok = (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)

    # ADMET Constraints
    admet_ok = (tpsa <= 140 and rot <= 10)

    return lipinski_ok and admet_ok

# =========================================================================
# 5. PIPELINE EXECUTION & EXPORT
# =========================================================================
if __name__ == "__main__":
    print("Initializing Library Generation for Post-Docking Scaffolds...")

    sdf_writer = Chem.SDWriter("gefitinib_aromatic.sdf")
    master_data = []

    for scaffold_name, scaffold_smiles in SCAFFOLDS.items():
        for sub_name, sub_smiles in SUBSTITUENTS.items():
            
            mol = generate_molecule(scaffold_smiles, sub_smiles)

            if mol is None:
                continue

            # Check properties before heavy 3D calculations
            if not passes_admet_filters(mol):
                continue

            # Prepare for 3D structure generation by adding explicit Hydrogens
            mol = Chem.AddHs(mol)

            # Embed 3D coordinates using ETKDG
            embed_status = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            if embed_status != 0:
                continue

            # Energy optimization via MMFF Force Field
            try:
                AllChem.MMFFOptimizeMolecule(mol)
            except:
                continue

            # Assign structured nomenclature
            compound_id = f"GF_{scaffold_name}_{sub_name}"
            mol.SetProp("_Name", compound_id)

            # Save structural data
            sdf_writer.write(mol)

            # Extract chemical properties for the final report
            master_data.append({
                "Compound_ID": compound_id,
                "Scaffold": scaffold_name,
                "Substituent": sub_name,
                "SMILES": Chem.MolToSmiles(Chem.RemoveHs(mol)),
                "MW": round(Descriptors.MolWt(mol), 2),
                "LogP": round(Crippen.MolLogP(mol), 2),
                "HBD": Lipinski.NumHDonors(mol),
                "HBA": Lipinski.NumHAcceptors(mol),
                "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
                "RotB": Lipinski.NumRotatableBonds(mol)
            })
            
            print(f"[+] Successfully optimized and stored: {compound_id}")

    sdf_writer.close()

    # =========================================================================
    # 6. MASTER REPORT EXPORTATION (CSV)
    # =========================================================================
    if master_data:
        df_master = pd.DataFrame(master_data)
        output_csv = "gefitinib_aromatic.csv"
        df_master.to_csv(output_csv, index=False)
        
        print("\n" + "="*60)
        print(f"Pipeline executed successfully!")
        print(f"Total passed compounds generated: {len(df_master)}")
        print(f"Structural SDF exported to: 'gefitinib_aromatic.sdf'")
        print(f"Descriptors CSV exported to: '{output_csv}'")
        print("="*60)
    else:
        print("\n[-] Error: No compounds passed the pipeline constraints.")
