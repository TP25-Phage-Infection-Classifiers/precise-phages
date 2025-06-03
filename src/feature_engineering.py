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
    for i, r in enumerate(records):

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
        uid = f"{r.id}_{i}"  # verhindert Überschreiben
        data[uid] = row

  #Erstelle ein DataFrame: Zeilen = Gene, Spalten = Features
    return pd.DataFrame.from_dict(data, orient="index")


# -----------------------------------------------------------------------------
# Hilfsfunktion: eine Datei verarbeiten
# -----------------------------------------------------------------------------

def extract_all_features():
    input_dir = Path("output/database_DNA")
    output_dir = Path("output/feature_engineering")
    output_dir.mkdir(parents=True, exist_ok=True)

    ks = [3, 4]
    transpose = False

    for fasta_path in input_dir.glob("*.fasta"):
        out_path = output_dir / f"{fasta_path.stem}_features.csv"

        # Wenn bereits vorhanden, überspringen
        if out_path.exists():
            continue

        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            continue

        df = extract_features(records, ks)
        if transpose:
            df = df.T

        df.to_csv(out_path)