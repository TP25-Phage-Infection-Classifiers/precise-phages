import files
from Bio.PDB import PDBParser, DSSP
import subprocess
import os
import pandas as pd
import warnings #*
from pathlib import Path 
import shutil

""" # 1) Binary & DB finden
PSIPRED_BIN = os.getenv("PSIPRED_BIN") or shutil.which("run_psipred.pl")
PSIPRED_DB  = os.getenv("PSIPRED_DB")

if not PSIPRED_BIN or not PSIPRED_DB:
    raise RuntimeError("PSIPRED nicht gefunden - bitte Env-Variablen setzen oder bootstrap-Script ausführen.")

# 2) Pfade
FASTA_DIR = files.output_dir / "database_Protein"
OUT_DIR   = Path("output/structure_features")
OUT_DIR.mkdir(parents=True, exist_ok=True) """




# PSIPRED - Aufruf 
#Starte PSIPRED und liefere die erzeugte *.ss2 Datei zurück.
def run_psipred(fasta_path: Path, out_prefix: Path) -> Path:
    ss2_file = out_prefix.with_suffix(".ss2")
    if ss2_file.exists():   # Überprüfe zunächst, ob das Ziel-.ss2-File bereits existiert
        return ss2_file
                             # Falls nicht, führe PSIPRED aus und erzeuge die .ss2-Datei.
    cmd = [
        PSIPRED_BIN,  #Pfad zu run_psipred.pl
        "-d", PSIPRED_DB, #Datenbank - BLAST-Datenbank im FASTA-Format
        str(fasta_path)
    ]
    print(f"Running PSIPRED: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return ss2_file

#Parsing der SS2-Datei für H E C Verteilung
# H = Helix
# E = Sheet
# S = Coil
# Ziel: Zähle, wie viele Reste als Helix (H), Strand (E) oder Coil (C) vorhergesagt wurden.
#  Am Ende werden relative Häufigkeiten (Anteile) für jede Sekundärstruktur-Klasse berechnet.
# Bei leerer oder fehlerhafter .ss2-Datei (total==0) liefern wir None-Werte zurück.

def parse_psipred(ss2_file: Path) -> dict:
    helix = sheet = coil = 0

    with open(ss2_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            _, _, sec, *_ = line.split()
            if sec == "H":
                helix += 1
            elif sec == "E":
                sheet += 1
            else:
                coil += 1
        
        total = helix + coil + sheet

        if total == 0:
            return {"sec_helix_frac": None, "sec_sheet_frac": None, "sec_coil_frac": None}
        return {
            "sec_helix_frac" : round(helix / total, 4),
            "sec_sheet_frac" : round(sheet / total, 4),
            "sec_coil_frac" : round(coil / total, 4),


        }
    

# Tertiärstruktur
# Erstellt mit ColabFold ein PDB-Modell 
# ColabFold erzeugt beim Lauf einen neuen Ordner und schreibt dort die Modelle hinein.
# Wir suchen nach dem ersten *_model_1.pdb

def run_colabfold(fasta_path, out_prefix):
    
    pdb_file = out_prefix.with_suffix(".pdb")
    if pdb_file.exists():
        return pdb_file

    #ColabFold Aufruf
    cmd = [
        COLABFOLD_BIN, 
        fasta_path,         #Eingabe
        out_prefix.parent,  #Output-Ordner
        "--model-type", "alphafold2_multimer" #Modus (Mono/Multi)
        ]
    subprocess.run(cmd, check=True)

    # ColabFold legt Datei unter out_prefix.parent/<seq>_model_1.pdb ab
    generated = next(out_prefix.parent.glob("*model_1.pdb"))
    generated.rename(pdb_file)
    return pdb_file



# Extraktion von 3D-Features aus dem PDB-Modell:
# - pLDDT: Mittlere Vorhersagequalität, gespeichert als B-Factor (Spalte in PDB)
# - rel_acc_mean: Durchschnittliche relative Zugänglichkeit (ASA) der Aminosäuren,
#   berechnet mit DSSP (Spalte 3)
# - num_domains: Anzahl der Ketten als grobe Schätzung für die Anzahl von Domänen

def extract_tertiary_features(pdb_file):
  
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file)

    # pLDDT Mittelwert (AlphaFold speichert ihn in B-Factor-Spalte)
    b_factors = [atom.bfactor for atom in structure.get_atoms()]
    plddt_mean = sum(b_factors) / len(b_factors) if b_factors else 0

    # DSSP für Oberfläche & Sek-Struktur
    model = structure[0]
    dssp = DSSP(model, pdb_file, dssp='mkdssp')
    acc = [entry[3] for entry in dssp]  # relative ASA
    rel_acc_mean = sum(acc) / len(acc) if acc else 0

    # einfache Domänen-Abschätzung: zähle zusammenhängende Ketten
    chains = {chain.id for chain in model}
    num_domains = len(chains)

    return {
        "plddt_mean": round(plddt_mean, 2),
        "rel_acc_mean": round(rel_acc_mean, 3),
        "num_domains": num_domains,
    }


# Berechne PSIPRED Features für alle FASTA-Dateien und schreibe eine CSV

# Hauptfunktion:
# Durchläuft alle FASTA-Dateien im Datenverzeichnis.
# Führt für jedes Protein sowohl Sekundärstruktur- (PSIPRED) als auch
# Tertiärstruktur-Analysen (ColabFold + DSSP) durch.
# Ergebnisse werden zu einer Tabelle zusammengeführt und als CSV gespeichert.
# Fehler bei einzelnen Proteinen werden abgefangen, das restliche Processing läuft weiter.

def build_structure_feature_table():
    rows = []

    fasta_dir = files.output_dir / "database_Protein"
    for fasta_path in fasta_dir.glob("*.fasta"):
        gene_id = fasta_path.stem.replace("_genes_protein", "")
        out_prefix = OUT_DIR / gene_id

        # --- Sekundärstruktur ---
        try:
            ss2 = run_psipred(str(fasta_path), str(out_prefix))
            sec_feats = parse_psipred(ss2)
        except Exception as e:
            warnings.warn(f"PSIPRED failed for {gene_id}: {e}")
            sec_feats = {"sec_helix_frac": None,
                         "sec_sheet_frac": None,
                         "sec_coil_frac": None}

        # --- Tertiärstruktur ---
        try:
            pdb = run_colabfold(str(fasta_path), out_prefix)
            tert_feats = extract_tertiary_features(pdb)
        except Exception as e:
            warnings.warn(f"ColabFold/DSSP failed for {gene_id}: {e}")
            tert_feats = {"plddt_mean": None,
                          "rel_acc_mean": None,
                          "num_domains": None}

        row = {"GeneID": gene_id, **sec_feats, **tert_feats}
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "structure_features.csv", index=False)
    print(" Struktur-Features gespeichert:", OUT_DIR / "structure_features.csv")

