import os
import pandas as pd


def merge_affinities():
    # Automatically get the directory where this script (bind.py) is located
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 1. Define file names with absolute paths
    file_library_d = os.path.join(script_dir, "screening_results_library_D.xlsx")
    file_gefitinib = os.path.join(script_dir, "gefitinib_derivates_report.xlsx")
    output_file = os.path.join(script_dir, "final_affinities.xlsx")

    # Verify that both Excel files exist
    if not os.path.exists(file_library_d) or not os.path.exists(file_gefitinib):
        print(
            f"Error: Make sure '{os.path.basename(file_library_d)}' and "
            f"'{os.path.basename(file_gefitinib)}' are in: {script_dir}"
        )
        return

    print("Reading Excel files...")
    df_library_d = pd.read_excel(file_library_d)
    df_gefitinib = pd.read_excel(file_gefitinib)

    print("Merging datasets...")
    # Merge using 'Ligand' for the left df and 'Compound_ID' for the right df
    df_merged = pd.merge(
        df_library_d,
        df_gefitinib,
        left_on="Ligand",
        right_on="Compound_ID",
        how="outer",
        suffixes=("_libD", "_gefitinib"),
    )

    # Clean up redundant ID columns from the outer join
    df_merged["Ligand"] = df_merged["Ligand"].fillna(df_merged["Compound_ID"])
    df_merged.drop(columns=["Compound_ID"], inplace=True)

    # 3. Save the results to a new Excel file
    df_merged.to_excel(output_file, index=False)
    print(f"Success! File saved as: '{output_file}'")
    print(f"Total ligands processed: {len(df_merged)}")


if __name__ == "__main__":
    merge_affinities()
