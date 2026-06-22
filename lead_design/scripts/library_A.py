import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen
from rdkit.Chem import rdMolDescriptors

# =========================================================================
# 1. INPUT DATA: GEFITINIB SCAFFOLD & TARGET SUBSTITUENTS
# =========================================================================

# Core structure of Gefitinib with a manual dummy atom [*] as the attachment point
GEFITINIB_SCAFFOLD = "COC1=C(OCCC[*])C=C2C(N=CN=C2NC3=CC(Cl)=C(F)C=C3)=C1"

# Target amine substituents (Starting with the connecting Nitrogen atom)
SUBSTITUENTS = {
    "Morpholine": "N1CCOCC1",
    "Piperidine": "N1CCCCC1",
    "Piperazine": "N1CCNCC1",
    "N-Methylpiperazine": "N1CCN(C)CC1",
    "Pyrrolidine": "N1CCCC1",
    "Azepane": "N1CCCCCC1",
    "Dimethylamine": "N(C)C",         
    "Diethylamine": "N(CC)CC"        
}

# =========================================================================
# 2. FRAGMENT COUPLING ENGINE (COMBINE & EDIT BONS)
# =========================================================================

def generate_molecule(base_scaffold, substituent):
    """
    Couples the scaffold and the substituent by finding the dummy atom,
    forming a single bond with the substituent's first atom, and removing the dummy.
    """
    mol_base = Chem.MolFromSmiles(base_scaffold)
    mol_sub = Chem.MolFromSmiles(substituent)

    if mol_base is None or mol_sub is None:
        return None

    # Combine both structures into a single editable molecule container
    combo = Chem.CombineMols(mol_base, mol_sub)
    rw = Chem.RWMol(combo)

    # Find the dummy atom (Atomic Number = 0)
    dummy_idx = None
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() == 0:
            dummy_idx = atom.GetIdx()
            break

    if dummy_idx is None:
        return None

    # Identify the scaffold neighbor atom connected to the dummy atom
    dummy_neighbors = [n.GetIdx() for n in rw.GetAtomWithIdx(dummy_idx).GetNeighbors()]
    if len(dummy_neighbors) != 1:
        return None

    scaffold_attach_idx = dummy_neighbors[0]
    
    # The first atom of the substituent starts exactly where the base ends
    sub_start_idx = mol_base.GetNumAtoms()

    # Create the new single bond and discard the dummy atom
    rw.AddBond(scaffold_attach_idx, sub_start_idx, Chem.BondType.SINGLE)
    rw.RemoveAtom(dummy_idx)

    final_mol = rw.GetMol()

    try:
        Chem.SanitizeMol(final_mol)
        return final_mol
    except:
        return None

# =========================================================================
# 3. PHYSICOCHEMICAL & ADMET FILTERS
# =========================================================================

def passes_admet_filters(mol):
    """
    Evaluates strict Lipinski Rule of 5 and basic ADMET criteria.
    Returns True only if all compound thresholds are met.
    """
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot = Lipinski.NumRotatableBonds(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)

    # Lipinski Rule of 5 parameters
    lipinski_ok = (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)

    # Drug-likeness ADMET constraints (Veber / Ghose filters)
    admet_ok = (tpsa <= 140 and rot <= 10)

    return lipinski_ok and admet_ok

# =========================================================================
# 4. PIPELINE EXECUTION & 3D OPTIMIZATION
# =========================================================================

if __name__ == "__main__":
    print("Initializing chemical library generation pipeline...")

    # Set up the SDWriter to save the 3D molecular conformations
    sdf_writer = Chem.SDWriter("gefitinib_library.sdf")
    master_data = []
    counter = 1

    # Iterate through the dictionary of substituents to derive compounds
    for name, sub_smiles in SUBSTITUENTS.items():
        mol = generate_molecule(GEFITINIB_SCAFFOLD, sub_smiles)

        if mol is None:
            print(f"[-] Failed to couple substituent: {name}")
            continue

        # Apply molecular filters before heavy 3D calculations
        if not passes_admet_filters(mol):
            print(f"[-] Molecule with {name} failed ADMET/Lipinski thresholds.")
            continue

        # Prepare molecule for 3D structure generation by adding Hydrogens
        mol = Chem.AddHs(mol)

        # Embed 3D coordinates using the experimental ETKDG algorithm
        embed_status = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if embed_status != 0:
            print(f"[-] 3D Embedding failed for derivative: {name}")
            continue

        # Optimize the 3D structure using the Molecular Mechanics Force Field (MMFF)
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except:
            print(f"[-] MMFF Optimization failed for derivative: {name}")
            continue

        # Naming convention
        compound_id = f"GF_{counter}_{name}"
        mol.SetProp("_Name", compound_id)

        # Write 3D coordinates to the structural file
        sdf_writer.write(mol)

        # Extract molecular descriptors for tracking and logging
        master_data.append({
            "Compound_ID": compound_id,
            "SMILES": Chem.MolToSmiles(Chem.RemoveHs(mol)),
            "MW": round(Descriptors.MolWt(mol), 2),
            "LogP": round(Crippen.MolLogP(mol), 2),
            "HBD": Lipinski.NumHDonors(mol),
            "HBA": Lipinski.NumHAcceptors(mol),
            "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
            "RotB": Lipinski.NumRotatableBonds(mol)
        })
        
        print(f"[+] Successfully optimized and stored: {compound_id}")
        counter += 1

    sdf_writer.close()

    # =========================================================================
    # 5. MASTER REPORT EXPORTATION (PANDAS DATA)
    # =========================================================================
    if master_data:
        df_master = pd.DataFrame(master_data)
        output_csv = "gefitinib_derivatives_properties.csv"
        df_master.to_csv(output_csv, index=False)
        
        print("\n" + "="*60)
        print(f"Pipeline executed successfully! Total compounds generated: {len(df_master)}")
        print(f"Structural outputs exported to: 'gefitinib_library.sdf'")
        print(f"Tabular descriptors exported to: '{output_csv}'")
        print("="*60)
    else:
        print("\n[-] Error: No compounds passed the pipeline constraints.")
