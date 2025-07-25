""" 
Extrahiert strukturelle Merkmale aus Proteinmodellen (Pdb Dateien), 
die mithilfe von AlphaFold bzw ColabFold generiert wurden.

 Extrahierte Struktur-Features:
1. `num_atoms`: Anzahl aller Atome im Modell (räumliche Komplexität).
2. `num_chains`: Anzahl getrennter Ketten - relevant bei Multimeren.
3. `num_residues`: Anzahl Aminosäurereste - korreliert mit Sequenzlänge.

 Sekundärstruktur-Merkmale (via DSSP):
4. `percent_helix`: Anteil an Helices (H, G, I) - z.B. Alpha-Helices.
5. `percent_sheet`: Anteil an β-Faltblatt-Strukturen (E, B).
6. `percent_coil`: Anteil ungeordneter/loopartiger Bereiche - strukturelle Flexibilität.
7. `ASA_mean`: Mittlere lösliche Oberfläche (Accessible Surface Area) in Å² 
   Rückschluss auf Oberflächenexposition und potenzielle Bindungsregionen.


"""

import zipfile
import csv
from pathlib import Path
from Bio.PDB import PDBParser, DSSP
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

# DSSP benötigt eine gültige Header Zeile in .pdb Dateien
# AlphaFold generierte Dateien haben diese nicht
# Die Funktion prüft das und fügt bei Bedarf einen Header ein
def fix_pdb_header(pdb_path: Path):
    content = pdb_path.read_text()
    if not content.startswith("HEADER"):
        fixed = "HEADER    DUMMY HEADER\n" + content
        pdb_path.write_text(fixed) 


# Entpackt die ZIP Datei mit den ColabFold Ergebnissen
def unzip_results(zip_path: Path, extract_to: Path):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)
    print(f" ZIP entpackt nach: {extract_to}")


#Extrahiert einfache strukturelle Merkmale aus PDB Dateien:
# Anzahl Atome, Ketten, Aminosäuren
#Speichert die Ergebnisse als CSV
def extract_structure_features(pdb_dir: Path, output_csv: Path):
    parser = PDBParser(QUIET=True)
    features = []

    # Oteration über alle .pdb Dateien 
    for pdb_file in pdb_dir.glob("*.pdb"):
        structure_id = pdb_file.stem
        structure = parser.get_structure(structure_id, pdb_file)

        atoms = list(structure.get_atoms())
        chains = list(structure.get_chains())
        residues = list(structure.get_residues())

        #Initialwerte für Sekundärstruktur und ASA Features
        percent_helix = percent_sheet = percent_coil = asa_mean = 0.0


        # DSSP: Tool zur Analyse von Sekundärstruktur & Oberfläche
        try:
            model = structure[0]
            fix_pdb_header(pdb_file)  
            dssp = DSSP(model, pdb_file)

            ss_counts = {"H": 0, "E": 0, "C": 0}
            asa_total = 0.0
            asa_residues = 0 

            for key in dssp.keys():

                ss  = dssp[key][1]          
                asa = dssp[key][2] 

                #Klassifikation der Sekundärstruktur
                if ss in ("H", "G", "I"): #Alpha-Helix, 3₁₀-Helix, π-Helix
                    ss_counts["H"] += 1 #Helix
                elif ss in ("E", "B"): #β-Strang, Brücke
                    ss_counts["E"] += 1 #Beta-Faltblatt
                else:
                    ss_counts["C"] += 1 #Kurven, Turns, ungeordnet

                 # ASA addieren – nur wenn numerisch
                try:
                    asa_val = float(asa)
                    asa_total += asa_val
                    asa_residues += 1
                except ValueError:
                    pass                    # '-' wird übersprungen

            # Prozentuale Verteilung und mittlere ASA berechnen
            total = sum(ss_counts.values())
            if total == 0:
                print("DSSP WARN:", pdb_file, "enthält 0 Residues - wurde DSSP korrekt ausgeführt?")

            if total > 0:
                percent_helix = ss_counts["H"] / total 
                percent_sheet = ss_counts["E"] / total 
                percent_coil = ss_counts["C"] / total 
                asa_mean = asa_total / asa_residues if asa_residues else 0.0

               
        except Exception as e:
            print(f"DSSP fehlgeschlagen für {pdb_file.name}: {e}")

        #Zeile für Features zur Liste hinzufügen
        features.append({
            "GeneID": structure_id,
            "num_atoms": len(atoms),
            "num_chains": len(chains),
            "num_residues": len(residues),
            "percent_helix": round(percent_helix, 3),
            "percent_sheet": round(percent_sheet, 3),
            "percent_coil": round(percent_coil, 3),
            #"ASA_mean": round(asa_mean, 3),
        })
    
    #Falls keine .pdb Dateien gefunden werden
    if not features:
        print(" Keine .pdb-Dateien gefunden.")
        return

    #Speichern der extrahierten Features in einer CSV Datei
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=features[0].keys())
        writer.writeheader()
        writer.writerows(features)

    print(f" Struktur-Features extrahiert für {len(features)} Dateien → {output_csv}")


def main():
    #Pfade zur ZIP Datei und zum entpacken PDB Datei
    zip_path = Path("data/colabfold_results.zip")
    pdb_dir = Path("output/colabfold_results/content/results")
    feature_output = Path("output/structure_features.csv")

    #Entpacke ZIP   
    unzip_results(zip_path, pdb_dir.parent.parent)
    
    #Extrahiere alle Struktur-Features
    extract_structure_features(pdb_dir, feature_output)

if __name__ == "__main__":
    main()