import zipfile
import csv
from pathlib import Path
from Bio.PDB import PDBParser

def unzip_results(zip_path: Path, extract_to: Path):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)
    print(f" ZIP entpackt nach: {extract_to}")

def extract_structure_features(pdb_dir: Path, output_csv: Path):
    parser = PDBParser(QUIET=True)
    features = []

    for pdb_file in pdb_dir.glob("*.pdb"):
        structure_id = pdb_file.stem
        structure = parser.get_structure(structure_id, pdb_file)

        atoms = list(structure.get_atoms())
        chains = list(structure.get_chains())
        residues = list(structure.get_residues())

        features.append({
            "structure": structure_id,
            "num_atoms": len(atoms),
            "num_chains": len(chains),
            "num_residues": len(residues),
        })

    if not features:
        print(" Keine .pdb-Dateien gefunden. Feature-Liste leer.")
        return

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=features[0].keys())
        writer.writeheader()
        writer.writerows(features)

    print(f" Struktur-Features extrahiert für {len(features)} Dateien → {output_csv}")

def main():
    zip_path = Path("data/colabfold_results.zip")
    pdb_dir = Path("output/colabfold_results/content/results")
    feature_output = Path("output/structure_features.csv")

    unzip_results(zip_path, pdb_dir.parent.parent)
    extract_structure_features(pdb_dir, feature_output)

if __name__ == "__main__":
    main()