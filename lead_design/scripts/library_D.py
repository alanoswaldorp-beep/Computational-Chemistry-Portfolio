import os
import pandas as pd
import itertools
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen
from rdkit.Chem import rdMolDescriptors

# =========================================================================
# 1. CORE SCAFFOLDS (CORRECTED WITH AROMATIC SMILES TO PREVENT 3D DISTORTION)
# =========================================================================
SCAFFOLDS = {
    # Clean aromatic 6-6 fused rings
    "Quinazoline": "[4*]c1cc2ncnc(Nc3cc([3*])c([2*])cc3)c2cc1OCCC[1*]",
    "Pyridopyrimidine": "[4*]c1nc2ncnc(Nc3cc([3*])c([2*])cc3)c2cc1OCCC[1*]",
    "Cinnoline": "[4*]c1cc2nncc(Nc3cc([3*])c([2*])cc3)c2cc1OCCC[1*]"
}

# =========================================================================
# 2. TARGET SUBSTITUENTS (MAPPED TO ATTACHMENT POINTS)
# =========================================================================
SUBSTITUENTS_R1 = {
    "Pip": "N1CCNCC1",
    "Mor": "N1CCOCC1",
    "NMP": "N1CCN(C)CC1"
}

SUBSTITUENTS_R2 = {
    "F": "F",
    "OH": "O",
    "Cl": "Cl"
}

SUBSTITUENTS_R3 = {
    "OH": "O",
    "F": "F",
    "Cl": "Cl"
}

SUBSTITUENTS_R4 = {
    "OCH3": "OC",
    "NO2": "[N+](=O)[O-]"
}

# =========================================================================
# 3. FRAGMENT COUPLING ENGINE (STRING TO RDKIT OBJECT)
# =========================================================================
def generate_molecule(scaffold_smiles, r1, r2, r3, r4):
    """
    Couples multiple substituents to the scaffold via targeted string replacement
    before creating the RDKit molecule object.
    """
    # Replace attachment points with corresponding SMILES fragments
    assembled_smiles = scaffold_smiles.replace("[1*]", r1)
    assembled_smiles = assembled_smiles.replace("[2*]", r2)
    assembled_smiles = assembled_smiles.replace("[3*]", r3)
    assembled_smiles = assembled_smiles.replace("[4*]", r4)

    # Convert the assembled string into an RDKit Mol object
    final_mol = Chem.MolFromSmiles(assembled_smiles)
    
    if final_mol is None:
        return None

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
    print("Initializing Combinatorial Library Generation...")

    output_sdf = "gefitinib_derivatives_3D.sdf"
    sdf_writer = Chem.SDWriter(output_sdf)
    master_data = []

    # Create all possible combinations of Scaffolds and Substituents
    combinations = itertools.product(
        SCAFFOLDS.items(),
        SUBSTITUENTS_R1.items(),
        SUBSTITUENTS_R2.items(),
        SUBSTITUENTS_R3.items(),
        SUBSTITUENTS_R4.items()
    )

    total_generated = 0

    for scaffold, r1, r2, r3, r4 in combinations:
        scaffold_name, scaffold_smiles = scaffold
        r1_name, r1_smiles = r1
        r2_name, r2_smiles = r2
        r3_name, r3_smiles = r3
        r4_name, r4_smiles = r4

        # Assemble the molecule
        mol = generate_molecule(scaffold_smiles, r1_smiles, r2_smiles, r3_smiles, r4_smiles)

        if mol is None:
            continue

        # 1. Filter Check (Fast calculation before 3D embedding)
        if not passes_admet_filters(mol):
            continue

        # 2. Prepare 3D structure
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

        # 3. Create ID and save data
        compound_id = f"GF_{scaffold_name}_{r1_name}_{r2_name}_{r3_name}_{r4_name}"
        mol.SetProp("_Name", compound_id)

        sdf_writer.write(mol)

        master_data.append({
            "Compound_ID": compound_id,
            "Scaffold": scaffold_name,
            "R1_Tail": r1_name,
            "R2_Sub": r2_name,
            "R3_Sub": r3_name,
            "R4_Head": r4_name,
            "SMILES": Chem.MolToSmiles(Chem.RemoveHs(mol)),
            "MW": round(Descriptors.MolWt(mol), 2),
            "LogP": round(Crippen.MolLogP(mol), 2),
            "HBD": Lipinski.NumHDonors(mol),
            "HBA": Lipinski.NumHAcceptors(mol),
            "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
            "RotB": Lipinski.NumRotatableBonds(mol)
        })
        
        total_generated += 1
        print(f"[+] Successfully optimized and stored: {compound_id}")

    sdf_writer.close()

    # =========================================================================
    # 6. MASTER REPORT EXPORTATION (CSV)
    # =========================================================================
    if master_data:
        df_master = pd.DataFrame(master_data)
        output_csv = "gefitinib_derivatives_report.csv"
        df_master.to_csv(output_csv, index=False)
        
        print("\n" + "="*60)
        print(f"Pipeline executed successfully!")
        print(f"Total theoretical combinations: 162")
        print(f"Total passed constraints and exported: {total_generated}")
        print(f"Structural SDF exported to: '{output_sdf}'")
        print(f"Descriptors CSV exported to: '{output_csv}'")
        print("="*60)
    else:
        print("\n[-] Error: No compounds passed the pipeline constraints.")
