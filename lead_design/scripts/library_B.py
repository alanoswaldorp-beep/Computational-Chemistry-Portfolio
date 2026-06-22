import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen
from rdkit.Chem import rdMolDescriptors

# =========================================================================
# 1. SCAFFOLDS WITH TWO DUMMY ATOMS ([1*] and [2*])
# =========================================================================
SCAFFOLDS = {
    "GF_Piperazine": "COC1=C(OCCCN2CCNCC2)C=C3C(N=CN=C3NC4=CC([2*])=C([1*])C=C4)=C1",
    "GF_N-Methylpiperazine": "COC1=C(OCCCN2CCN(C)CC2)C=C3C(N=CN=C3NC4=CC([2*])=C([1*])C=C4)=C1",
    "GF_Morpholine": "COC1=C(OCCCN2CCOCC2)C=C3C(N=CN=C3NC4=CC([2*])=C([1*])C=C4)=C1"
}

# =========================================================================
# 2. SUBSTITUENT LISTS FOR EACH POSITION
# =========================================================================
# Target substituents for dummy atom [1*]
SUBSTITUENTS_1 = {
    "F": "F",
    "H": "[H]",
    "Cl": "Cl",
    "Br": "Br",
    "OH": "O"
}

# Target substituents for dummy atom [2*]
SUBSTITUENTS_2 = {
    "F": "F",
    "H": "[H]",
    "Cl": "Cl",
    "Br": "Br",
    "OH": "O"
}

# =========================================================================
# 3. FRAGMENT COUPLING ENGINE (DUAL SUBSTITUTION)
# =========================================================================
def generate_disubstituted_molecule(base_smiles, sub1_smiles, sub2_smiles):
    """
    Sequentially attaches two substituents to a scaffold based on isotopic 
    dummy atoms ([1*] and [2*]).
    """
    mol_base = Chem.MolFromSmiles(base_smiles)
    mol1 = Chem.MolFromSmiles(sub1_smiles)
    mol2 = Chem.MolFromSmiles(sub2_smiles)

    if None in (mol_base, mol1, mol2):
        return None

    # --- STEP 1: Attach substituent 1 to [1*] ---
    combo1 = Chem.CombineMols(mol_base, mol1)
    rw1 = Chem.RWMol(combo1)
    
    dummy1_idx = None
    for atom in rw1.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == 1:
            dummy1_idx = atom.GetIdx()
            break
            
    if dummy1_idx is not None:
        neighbors = rw1.GetAtomWithIdx(dummy1_idx).GetNeighbors()
        if neighbors:
            scaffold_anchor_1 = neighbors[0].GetIdx()
            sub1_start_idx = mol_base.GetNumAtoms()
            
            rw1.AddBond(scaffold_anchor_1, sub1_start_idx, Chem.BondType.SINGLE)
            rw1.RemoveAtom(dummy1_idx)
            
    intermediate_mol = rw1.GetMol()
    
    # --- STEP 2: Attach substituent 2 to [2*] ---
    combo2 = Chem.CombineMols(intermediate_mol, mol2)
    rw2 = Chem.RWMol(combo2)
    
    dummy2_idx = None
    for atom in rw2.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == 2:
            dummy2_idx = atom.GetIdx()
            break
            
    if dummy2_idx is not None:
        neighbors2 = rw2.GetAtomWithIdx(dummy2_idx).GetNeighbors()
        if neighbors2:
            scaffold_anchor_2 = neighbors2[0].GetIdx()
            sub2_start_idx = intermediate_mol.GetNumAtoms()
            
            rw2.AddBond(scaffold_anchor_2, sub2_start_idx, Chem.BondType.SINGLE)
            rw2.RemoveAtom(dummy2_idx)

    final_mol = rw2.GetMol()
    
    try:
        Chem.SanitizeMol(final_mol)
        # Normalize molecule to remove explicit [H] atoms introduced by the dictionary
        final_mol = Chem.RemoveHs(final_mol)
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

    # Lipinski Rule of 5 parameters
    lipinski_ok = (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)

    # ADMET constraints
    admet_ok = (tpsa <= 140 and rot <= 10)

    return lipinski_ok and admet_ok

# =========================================================================
# 5. PIPELINE EXECUTION & EXPORT
# =========================================================================
if __name__ == "__main__":
    print("Initializing Combinatorial Library Generation...")

    sdf_writer = Chem.SDWriter("gefitinib_combinatorial_library.sdf")
    master_data = []

    # Nested loop to cross every scaffold with every possible substitution pattern
    for scaffold_name, scaffold_smiles in SCAFFOLDS.items():
        for name_1, smi_1 in SUBSTITUENTS_1.items():
            for name_2, smi_2 in SUBSTITUENTS_2.items():
                
                # Generate molecule
                mol = generate_disubstituted_molecule(scaffold_smiles, smi_1, smi_2)

                if mol is None:
                    continue

                # Apply chemical filters
                if not passes_admet_filters(mol):
                    continue

                # Prepare for 3D coordinate generation
                mol = Chem.AddHs(mol)

                # Embed 3D coordinates
                embed_status = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                if embed_status != 0:
                    continue

                # Energy minimization via MMFF
                try:
                    AllChem.MMFFOptimizeMolecule(mol)
                except:
                    continue

                # Automated Nomenclature Setup (e.g., GF_Morpholine_1-F_2-Cl)
                compound_id = f"{scaffold_name}_1-{name_1}_2-{name_2}"
                mol.SetProp("_Name", compound_id)

                sdf_writer.write(mol)

                # Store physicochemical data for Pandas
                master_data.append({
                    "Compound_ID": compound_id,
                    "Scaffold": scaffold_name,
                    "Sub_1": name_1,
                    "Sub_2": name_2,
                    "SMILES": Chem.MolToSmiles(Chem.RemoveHs(mol)),
                    "MW": round(Descriptors.MolWt(mol), 2),
                    "LogP": round(Crippen.MolLogP(mol), 2),
                    "HBD": Lipinski.NumHDonors(mol),
                    "HBA": Lipinski.NumHAcceptors(mol),
                    "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
                    "RotB": Lipinski.NumRotatableBonds(mol)
                })
                
                print(f"[+] Optimized & Stored: {compound_id}")

    sdf_writer.close()

    # =========================================================================
    # 6. CSV REPORT EXPORTATION
    # =========================================================================
    if master_data:
        df_master = pd.DataFrame(master_data)
        output_csv = "gefitinib_combinatorial_properties.csv"
        df_master.to_csv(output_csv, index=False)
        
        print("\n" + "="*60)
        print(f"Combinatorial Pipeline executed successfully!")
        print(f"Total novel derivatives generated: {len(df_master)}")
        print(f"3D structures exported to: 'gefitinib_combinatorial_library.sdf'")
        print(f"Descriptors exported to: '{output_csv}'")
        print("="*60)
    else:
        print("\n[-] Error: No compounds passed the pipeline constraints.")
