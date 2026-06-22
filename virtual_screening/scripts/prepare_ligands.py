import os
import subprocess
from rdkit import Chem

# =========================================================
# Archive SDF
# =========================================================

sdf_file = "gefitinib_library_A.sdf"

# =========================================================
# Output
# =========================================================

output_dir = "pdbqt_ligandos"

os.makedirs(output_dir, exist_ok=True)

# =========================================================
# Molecules
# =========================================================

supplier = Chem.SDMolSupplier(sdf_file)

for mol in supplier:

    if mol is None:
        continue

    # =====================================================
    # Name
    # =====================================================

    nombre = mol.GetProp("_Name")

    sdf_temp = os.path.join(output_dir, f"{nombre}.sdf")
    pdbqt_out = os.path.join(output_dir, f"{nombre}.pdbqt")

    # =====================================================
    # SDF
    # =====================================================

    writer = Chem.SDWriter(sdf_temp)
    writer.write(mol)
    writer.close()

    # =====================================================
    # OpenBabel
    # =====================================================

    comando = [
        "obabel",
        sdf_temp,
        "-O",
        pdbqt_out,
        "--gen3d",
        "-p",
        "7.4"
    ]

    try:

        subprocess.run(
            comando,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        print(f"Convertido: {nombre}")

    except subprocess.CalledProcessError as e:

        print(f"Error con {nombre}")
        print(e.stderr.decode())

    # =====================================================
    # Remove SDF 
    # =====================================================

    os.remove(sdf_temp)

print("\nConversión terminada.")
