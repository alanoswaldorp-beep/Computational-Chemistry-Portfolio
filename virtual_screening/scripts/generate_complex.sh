#!/bin/bash

echo "=== PHASE 1: GENERATING COMPLEXES (WITHOUT END TAG) ==="

# Check if the target receptor protein exists in the working directory
if [ ! -f "4WKQ_clean.pdb" ]; then
    echo "ERROR: '4WKQ_clean.pdb' not found."
    exit 1
fi

mkdir -p complexes_pdb

# --- THE TRUCK IS HERE ---
# Create a copy of the receptor protein but strip the "END" record from the file
grep -v "^END" 4WKQ_clean.pdb > protein_no_end.pdb

for folder in Ligand_*/; do
    ligand=$(basename "$folder")
    out_file="${folder}${ligand}_out.pdbqt"
    
    if [ -f "$out_file" ]; then
        echo "Preparing complex for: $ligand"
        
        # Convert the ligand from PDBQT to PDB (suppressing OpenBabel warnings to avoid clutter)
        obabel -ipdbqt "$out_file" -opdb -O "${folder}${ligand}_temp.pdb" 2>/dev/null
        
        # Concatenate the protein (without END tag) and the ligand into a single complex file
        cat protein_no_end.pdb "${folder}${ligand}_temp.pdb" > "complexes_pdb/complex_${ligand}.pdb"
        
        # Remove the temporary ligand PDB file from the subfolder
        rm "${folder}${ligand}_temp.pdb"
    fi
done

# Clean up the temporary intermediate protein file
rm protein_no_end.pdb

echo "=== PHASE 1 COMPLETE: Ready for PLIP ==="
