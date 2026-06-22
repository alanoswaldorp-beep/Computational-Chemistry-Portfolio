#!/bin/bash

# Create the CSV file and write the column headers
echo "Ligand,Affinity_kcal_mol" > screening_results.csv

for f in GF_*.pdbqt; do
    # Validate if matching files exist to avoid loop errors
    [ -e "$f" ] || continue

    b=$(basename "$f" .pdbqt)
    echo "Processing ligand $b"
    
    mkdir -p "$b"
    
    # Run Vina and redirect the entire output to a temporary log file using '>'
    vina --config config.txt --ligand "$f" --out "${b}/${b}_out.pdbqt" > "${b}/temp_log.txt"
    
    # --- EXTRACT DATA FOR THE CSV HERE ---
    # Search for the line starting with '1' in the energy table and extract the 2nd column
    affinity=$(grep -E '^[[:space:]]*1[[:space:]]+' "${b}/temp_log.txt" | awk '{print $2}')
    
    # If it failed or wasn't calculated for some reason, assign "ERROR" so it doesn't stay empty
    if [ -z "$affinity" ]; then
        affinity="ERROR"
    fi
    
    # Append the ligand name and its affinity to the CSV file (comma-separated)
    echo "$b,$affinity" >> screening_results.csv
    
    # Optional: Uncomment the following line to delete the temporary log and save space
    # rm "${b}/temp_log.txt"
    
done

echo "Screening complete! Results have been saved to screening_results.csv"
