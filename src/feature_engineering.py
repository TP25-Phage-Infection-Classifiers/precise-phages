"""

Extrahiert Features  aus FASTA Einzelndateien oder Batch.

Features
---------
* Länge des Gens
* GC-Gehalt (%)
* Motiv-Vorkommen (TATA-Box, Sigma-70-Motife)
* K-Mer-Frequenzen (variable k-Liste)

Wie soll man das ausführen:
----
* Einzelfile: `python feature_engineering.py path/to/file.fasta [-o out.csv] [--transpose]`
* Batch:      `python feature_engineering.py --batch [--transpose]`
  → nutzt `input_dir` und `my_genome_files` aus `src/files.py` und legt alle CSVs unter
    `output/feature-engineering/` ab.
"""

import argparse
import importlib.util
import re
from collections import Counter
from itertools import product
from pathlib import Path

import pandas as pd
from Bio import SeqIO

#Versuche normalen Import
try:
    from src.files import input_dir, my_genome_files

#Wenn das Paket nicht erkannt wird:
except ImportError:

#src/files.py manuell laden
    spec = importlib.util.spec_from_file_location("files", Path(__file__).parent / "files.py")
    files = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(files)  
    #Zugriff auf die Variablen
    input_dir = files.input_dir
    my_genome_files = files.my_genome_files

# -----------------------------------------------------------------------------

#Reguläre Ausdrücke zur Erkennung bakterielle Promotor Motive
#Sigma 70 abhängige Transkriptionsstartstellen

MOTIFS = {
    "TATA": re.compile(r"TATA[AT]A[AT]"), #TATA-Box (mögliche Matches:TATAAA, TATAAT, TATATA, TATATT)
    "sigma35": re.compile(r"TTGACA"),     # -35 Region bakterielle Promotoren
    "sigma10": re.compile(r"TATAAT"),     # Pribnow Box oder -10 Region
}


#GC Verhätlnis
def compute_gc(seq: str) -> float:
    #Sequenz in Großbuchstaben umwandeln, damit man es richtig zählen kann
    seq = seq.upper()

    #Anzahl der G und C Basen
    gc = seq.count("G") + seq.count("C")

    #GC Gehalt in Prozent berechnen, auf 2 Nachkommastellen gerundet
    return round(100 * gc / len(seq), 2) if seq else 0.0


#Erzeugt alle möglichen Kombinationen von A, C, G, T der Länge k
# Beispiel: k=2 →'AA', 'AC', 'AG' 
def all_kmers(k: int):
    return ["".join(p) for p in product("ACGT", repeat=k)]


def kmer_freqs(seq: str, k: int):
    # In Großbuchtaben für Zählung
    seq = seq.upper()

    #Gesamtzahl möglicher K-Mere
    total = max(len(seq) - k + 1, 1)

    #Zähle alle K-Mere in der Sequenz
    counts = Counter(seq[i : i + k] for i in range(len(seq) - k + 1))

    #Normalisiere K-Mer Häufigkeiten auf (0,1), auf 4 Nachkommastellen runden
    return {
        f"k{k}_{mer}": round(counts.get(mer, 0) / total, 4) 
        for mer in all_kmers(k)
        }


def motif_counts(seq: str):
    #DNA Sequenz in Großbuchtaben für Mustererkennung
    seq = seq.upper()

    #Durchlaufe alle Motive und zähle,
    #wie oft jedes Motiv vorkommt
    #Als Ergebnis hat man Dictionary mit Namen und Vorkommenshäufigkeit
    return {
        f"motif_{name}": len(p.findall(seq))
          for name, p in MOTIFS.items()
          }


def extract_features(records, ks):
    #Initialisiere leeres Dictionary für alle Sequenzen
    data = {}

    #Durchlaufe alle FASTA Records (Liste von Bio.SeqRecord-Objekten aus FASTA)
    for r in records:

        # Länge und GC Gehalt
        row = {
            "length": len(r.seq),
            "GC_content": compute_gc(str(r.seq)),
        }

        #Zähle Promotor Motive
        row.update(motif_counts(str(r.seq)))

        #Berechne K-Mer Frequenzen für alle k-Werte
        for k in ks:
            row.update(kmer_freqs(str(r.seq), k))

        #Weise die Feature Zeile der  Sequenz-ID zu
        data[r.id] = row

  #Erstelle ein DataFrame: Zeilen = Gene, Spalten = Features
    return pd.DataFrame.from_dict(data, orient="index")


# -----------------------------------------------------------------------------
# Hilfsfunktion: eine Datei verarbeiten
# -----------------------------------------------------------------------------

def process_one(fasta_path: Path, out_path: Path, ks, transpose):
    #Ausgabe: Dateiname der FASTA Datei
    print(f"[•] {fasta_path} …")

    #Lade alle Sequenzen aus der FASTA Datei
    records = list(SeqIO.parse(fasta_path, "fasta"))

    #Wenn keine Sequenzen gefunden wurden
    if not records:
        print(f" Keine Sequenzen gefunden in {fasta_path}")
        return

    #Features extrahieren
    df = extract_features(records, ks)
    #Transponiere die Matrix Features = Zeilen, Gene = Spalten
    if transpose:
        df = df.T

    #Sicherstellen, dass der Zielordner existiert
    out_path.parent.mkdir(parents=True, exist_ok=True)
    #Speichere die Feature Tabelle als CSV Datei 
    df.to_csv(out_path)

    #Damit der Pfad in der Ausgabe besser aussieht
    try:
        rel = out_path.relative_to(Path.cwd())
    except ValueError:
        rel = out_path
    print(f"   gespeichert  {rel}")


# CLI

def main():
    ap = argparse.ArgumentParser()

    #Optionale Argumnte
    ap.add_argument("fasta", nargs="?", help="Einzelne FASTA‑Datei (optional, wenn --batch)")
    ap.add_argument("--output", "-o", help="Ausgabe‑CSV (nur Einzelmodus)")
    ap.add_argument("--k", nargs="+", type=int, default=[3, 4], help="K‑Mer‑Größen")
    ap.add_argument("--transpose", action="store_true", help="Matrix transponieren (Features=Zeilen)")
    ap.add_argument("--batch", action="store_true", help="Alle FASTA‑Dateien verarbeiten")
    args = ap.parse_args()


    #Batch: Alle Dateien aus files.py durchlaufen
    if args.batch:
        if not input_dir or not my_genome_files:
            raise RuntimeError("input_dir oder my_genome_files fehlen in files.py")
        for rel in my_genome_files:
            fasta_path = input_dir / rel
            out_path = Path("output/feature-engineering") / f"{Path(rel).stem}_features.csv"
            process_one(fasta_path, out_path, args.k, args.transpose)

    #Nur eine FASTA Datei verarbeiten
    else:
        if not args.fasta:
            ap.error("FASTA-Datei fehlt (oder --batch verwenden)")
        fasta_path = Path(args.fasta).expanduser().resolve()
        if not fasta_path.exists():
            ap.error(f"FASTA-Datei nicht gefunden: {fasta_path}")
        out = Path(args.output) if args.output else Path("output/feature-engineering/features.csv")
        process_one(fasta_path, out, args.k, args.transpose)


if __name__ == "__main__":
    main()
